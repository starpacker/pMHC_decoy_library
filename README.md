# Decoy Library — TCR 脱靶毒性 pMHC 负样本库

## 问题

工程化 TCR-T 疗法的核心安全风险是**交叉反应性**：TCR 错误识别正常组织细胞表面的 pMHC，攻击健康组织，导致严重不良事件甚至死亡。

> MAGE-A3 靶向 TCR 识别肿瘤肽 `EVDPIGHLY`，但交叉识别了心肌蛋白 Titin 的肽段 `ESDPIVAQY`（仅 1 个氨基酸差异），导致患者心脏毒性死亡。—— Cameron et al., Sci Transl Med, 2013

**Decoy Library 的目标**：构建一个全面的负样本库——收集所有可能被 TCR 错误识别的 pMHC，用于临床前安全性筛查和 IND 申报。

## 三层互补体系

Decoy Library 通过三条互补路线构建负样本，覆盖不同机制的交叉反应风险：

```
                        ┌─────────────────────────────┐
                        │       Decoy Library          │
                        │   TCR 脱靶毒性 pMHC 负样本库   │
                        └──────┬──────┬──────┬─────────┘
                               │      │      │
               ┌───────────────┘      │      └───────────────┐
               ▼                      ▼                      ▼
     ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
     │    Decoy A       │   │    Decoy B       │   │    Decoy C       │
     │   序列同源扫描    │   │  结构相似 + 逆向  │   │   文献挖掘        │
     │                  │   │   设计            │   │                  │
     │  "长得像的"       │   │  "形状像的"       │   │ "文献报道过的"     │
     │  Hamming ≤ 2     │   │  3D 表面相似      │   │  临床/实验证据     │
     │  50-200 条/靶标  │   │  100-500 条/靶标  │   │  105 条（已整理）  │
     │                  │   │                  │   │                  │
     │  确定性: 高       │   │  确定性: 中       │   │  确定性: 最高      │
     │  覆盖面: 窄       │   │  覆盖面: 最广     │   │  覆盖面: 受限      │
     └─────────────────┘   └─────────────────┘   └─────────────────┘
             │                      │                      │
             └──────────────────────┼──────────────────────┘
                                    ▼
                            统一风险评分与排名
                            → 湿实验验证 / IND 申报
```

---

## Decoy C — 文献与临床证据驱动（已完成）

### 原理

从 PubMed、bioRxiv、ClinicalTrials.gov 等 6 个数据源中挖掘已有实验证据或临床报告的脱靶毒性 pMHC。这些不是计算预测，而是真实发生过的交叉反应事件。

### 当前状态

- **105 条高质量条目**（有明确证据等级分类）
- 250 条总条目（含 LLM 自动提取的待审条目）
- 288 条广谱交叉反应参考（`cross_reactivity_reference.json`）

### 证据等级

| 等级 | 含义 | 数量 | 典型案例 |
|------|------|------|---------|
| Level 1: Clinical Fatal | 临床致死/FDA 暂停 | 5 | MAGE-A3 TCR → Titin → 心脏死亡 |
| Level 2: In Vitro Confirmed | 体外细胞毒性验证 | 61 | DMF5 TCR → MLANA → 黑色素细胞破坏 |
| Level 3: High-Throughput Screened | 高通量筛查命中 | 27 | X-scan / 酵母展示发现的交叉肽 |
| Level 4: In Silico High-Risk | 仅计算预测 | 12 | ARDitox / EpiTox 预测的高风险肽 |

### 实现架构

```
PubMed / bioRxiv / arXiv / ClinicalTrials.gov / Semantic Scholar / Europe PMC
   │         │         │           │                  │               │
   └─────────┴─────────┴───────────┴──────────────────┴───────────────┘
                                   │
                        Fetcher Agent (363+ 搜索策略)
                        ├── Strategy A: 149 关键词布尔查询
                        ├── Strategy B: 126 高风险基因专项查询
                        ├── Strategy C: 20 MeSH 受控词查询
                        ├── Strategy D: 48 年份分层查询
                        ├── Strategy E: LLM 动态生成查询
                        └── Strategy F: 77 多源特化查询
                                   │
                        Extractor Agent (LLM gpt-4o, T=0.1)
                        ├── 结构化 JSON 输出
                        └── thought_process 防幻觉机制
                                   │
                        Validator Agent
                        ├── UniProt: gene_symbol ↔ uniprot_id 验证
                        └── IEDB: 免疫表位数据库交叉验证
                                   │
                        去重 + 质量过滤 → decoy_library.json
```

### 核心代码

| 文件 | 职责 |
|------|------|
| `fetcher.py` | PubMed/PMC 文献获取 |
| `multi_source_fetcher.py` | arXiv/bioRxiv/medRxiv/ClinicalTrials.gov/Semantic Scholar/Europe PMC |
| `extractor.py` | LLM 驱动的结构化信息提取 + 规则回退 |
| `validator.py` | UniProt + IEDB 交叉验证 |
| `orchestrator.py` | 管线编排 + 持久化 |
| `scale_up.py` | 自动化扩量（6 种查询策略、断点续传、无限模式） |

### 数据 Schema

每条 Decoy C 记录包含 6 个核心字段组：

| 字段组 | 内容 |
|--------|------|
| `peptide_info` | decoy_sequence, hla_allele, source_protein, gene_symbol, uniprot_id |
| `discovery_context` | original_target_sequence, original_target_protein, tcr_name_or_id |
| `risk_profile` | evidence_level (4-tier), critical_organs_affected, expression_pattern |
| `experimental_evidence` | mass_spec_confirmed, assays_performed, cross_reactivity_affinity |
| `provenance` | pmid, clinical_trial_id, evidence_summary |
| `source` | title, authors, journal, year, doi, citation |

---

## Decoy A — 序列同源扫描（已完成，全管线已部署 + Presentation Score）

### 原理

从人类蛋白组中找出与靶标肽段在序列上高度相似（Hamming distance ≤ 2）且能被目标 HLA 呈递的内源性肽段。序列越相似，TCR 交叉识别的概率越大。

### 已部署数据资产

| 资产 | 规模 | 文件 |
|------|------|------|
| 人类 K-mer 字典 | **41,684,988** 条唯一 8-11mer | `human_kmer_db.parquet` (380 MB) |
| 组织表达谱 | **20,151** 基因 × 50 组织 | `gene_expression.parquet` (1.9 MB) |
| HLA-A\*02:01 binders (affinity) | **1,280,221** 条 (SB: 292,479 + WB: 987,742) | `hla_filtered_HLA-A0201.parquet` (27 MB) |
| HLA-A\*02:01 binders (presentation) | **971,925** 条 (SB: 337,584 + WB: 634,341) | `hla_filtered_HLA-A0201_presentation.parquet` (48 MB) |

### Presentation Score Model 升级

原始过滤使用 `Class1AffinityPredictor`（仅结合亲和力）。已升级为 `Class1PresentationPredictor`，整合抗原处理通路：

```
Affinity Model:   肽段 + HLA → 结合力 (IC50)           → "能否结合 MHC"
Presentation Model: 肽段 + HLA → 结合力 × 加工效率      → "能否真正呈递到细胞表面"
                                  ↑
                    蛋白酶体切割 + TAP 转运效率 (processing_score)
```

**Presentation Score 重新评分结果 (HLA-A\*02:01)**:

| 变化 | 数量 | 说明 |
|------|------|------|
| Weak → Strong (升级) | 91,586 | 高 processing score 补偿了中等亲和力 |
| Strong → Weak (降级) | 46,481 | affinity 虽强但 processing 差 |
| Weak → Non-Binder (移出) | **308,296** | processing score 极低，无法被呈递 |
| **有效 binder 池** | **971,925** | 从 128 万精确缩减至 97 万 (-24.1%) |

### 实现

```
Step 1: Swiss-Prot 20,431 蛋白 → 滑动窗口 → 41.7M 唯一 8-11mer
        + HPA RNA 组织表达谱（20,151 基因 × 50 组织，标注核心脏器 TPM）
        构建方式: SQLite 磁盘方案（适配 16GB RAM）

Step 2: mhcflurry 2.2.0 Class1PresentationPredictor
        → presentation_percentile ≤ 2.0 → 971,925 HLA-A*02:01 可呈递肽
        输出: affinity_nM + processing_score + presentation_score + presentation_percentile
        后端优先级: NetMHCpan 4.1 > mhcflurry > IEDB API
        结果缓存为 Parquet，每个 HLA 型别只需运行一次

Step 3A: NumPy 向量化 Hamming distance 扫描（~6 秒 / 75 万条 9-mer）
         → 筛选 distance ≤ 2 的候选
         → 标注每个错配: 锚定位(p2/p9) vs TCR接触面(p4-p8)
         → 关联组织表达谱（心脏/大脑/肺 TPM）
```

### 输出示例

**靶标 GILGFVFTL（流感 M1）/ HLA-A\*02:01** — 从 128 万可呈递肽中搜索：

```
#  Sequence      HD  EL%Rank  Genes     Mismatches
1  GLLGFVGTL      2    0.17   TAP2      p2:I>L^  p7:F>G*
2  GILLFLFTL      2    0.48   PKD1L1    p4:G>L*  p6:V>L*

*= TCR 接触面   ^= HLA 锚定位
```

TAP2（抗原加工相关转运体）和 PKD1L1（多囊肾蛋白）均为正常组织表达蛋白，且都是 HLA-A\*02:01 强结合肽（%Rank < 0.5）。

**靶标 EVDPIGHLY（MAGE-A3）/ 全蛋白组扫描**（无 HLA 过滤，10.4M 9-mer）：

```
#  Sequence      HD  Genes         Vital TPM  Category
1  EVDPIGHVY      1  MAGEA6              0.1  restricted
2  EVGPIFHLY      2  FGD5               20.0  high_risk     ← 肺/心肌高表达
3  EVVRIGHLY      2  MAGEA12             0.3  restricted
4  EVDPAGHSY      2  MAGEA8,MAGEA9       0.0  restricted
5  EVDPIRHYY      2  MAGEB18             0.0  silent
6  EVVPISHLY      2  MAGEA2              0.0  restricted
```

注：EVDPIGHLY 是 HLA-A\*01:01 限制性肽段，在 HLA-A\*02:01 binder 池中无命中是预期结果。对其他 HLA 型别扫描需先运行对应的 HLA 过滤。

### 核心代码

| 文件 | 职责 |
|------|------|
| `decoy_a/kmer_builder.py` | Swiss-Prot 下载 → k-mer 生成（SQLite 磁盘方案）→ HPA 表达谱 |
| `decoy_a/hla_filter.py` | HLA 呈递门控（三级回退: NetMHCpan > mhcflurry > IEDB） |
| `decoy_a/scanner.py` | Hamming 扫描 + 错配标注 + 表达谱关联 |
| `decoy_a/run.py` | 独立运行器（demo / 自定义肽段池 / 全管线） |
| `decoy_a/tools/mhcflurry.py` | mhcflurry 2.2.0 封装（本地 Python HLA 预测） |
| `decoy_a/tools/netmhcpan.py` | NetMHCpan 4.1 封装（可选，需学术许可） |
| `scripts/rescore_presentation.py` | Presentation Score 重新评分脚本 |

---

## Decoy B — 结构相似扫描 + 逆向设计（代码完成，待部署 GPU 工具）

### 原理

捕获序列差异大（Hamming > 2）但 TCR 接触面的 3D 构象和静电分布高度相似的肽段——Titin 致死案例就属于此类，是最危险也最难预测的一类脱靶。

### 实现：四阶段漏斗 + MPNN 逆向设计

```
Stage 1 — Atchley 理化初筛（分钟级，CPU）
  TCR 接触面残基 → Atchley Factors 25 维向量 → 余弦相似度 → Top 5000

Stage 2 — tFold 批量结构预测（天级，GPU）
  Top 5000 → tFold pMHC 3D 结构预测 → PDB

Stage 3 — AF3 精细预测（可选，小时级，GPU）
  Top 200 → AlphaFold3 高精度结构 → mmCIF

Stage 4 — 结构比对与评分
  肽段骨架 RMSD + 表面相似度 → 综合评分

MPNN 逆向设计分支（并行）
  靶标 pMHC 结构 → ProteinMPNN 固定锚定位、设计 TCR 接触面
  → 生成 ~1000 候选 → 回溯人类蛋白组 → NetMHCpan 过滤
```

### 核心代码

| 文件 | 职责 |
|------|------|
| `decoy_b/scanner.py` | 四阶段主管线 + MPNN 分支 |
| `decoy_b/tools/tfold.py` | tFold pMHC 结构预测封装 |
| `decoy_b/tools/alphafold3.py` | AlphaFold3 封装（Docker + 本地） |
| `decoy_b/tools/proteinmpnn.py` | ProteinMPNN 逆向序列设计封装 |
| `decoy_a/tools/netmhcpan.py` | NetMHCpan 4.1 HLA 结合预测封装 |

---

## 统一风险评分

无论来自 Decoy A、B 还是 C，所有候选通过统一公式排名：

```
Total Risk = Similarity × (1 / EL_Rank) × Vital_Organ_TPM_Weight
```

| 因子 | Decoy A | Decoy B | Decoy C |
|------|---------|---------|---------|
| Similarity | 1 - HD/len | 0.4×cosine + 0.6×surface | evidence_level 映射 |
| EL_Rank | NetMHCpan %Rank | 同左 | 文献值 |
| TPM_Weight | 心脑 TPM>10 → ×10 | 同左 | 同左 |

---

## 快速开始

### 环境

```bash
pip install pandas pyarrow numpy pydantic requests biopython mhcflurry
PYTHONUTF8=1 mhcflurry-downloads fetch models_class1_pan    # HLA 预测模型（~157 MB）
```

### Decoy C — 文献挖掘

```bash
# 查看已有的 105 条高质量负样本
python -m decoy_c show
python -m decoy_c stats

# 从特定论文提取
python -m decoy_c fetch --pmids 23926201 23863783

# 自动扩量到 N 条
cd decoy_c && python scale_up.py --required_num 500 --resume
```

### Decoy A — 序列同源扫描（全管线已部署）

```bash
# Demo 模式（秒级，无需任何下载）
python -m decoy_a demo

# 在 128 万 HLA-A*02:01 binders 中扫描（~6 秒）
PYTHONUTF8=1 python -m decoy_a scan --target GILGFVFTL --hla HLA-A*02:01
PYTHONUTF8=1 python -m decoy_a scan --target RMFPNAPYL --hla HLA-A*02:01

# 构建 k-mer 数据库（一次性，已完成，无需重复）
python -m decoy_a build-kmer

# 对其他 HLA 型别运行过滤（一次性，~6-7 小时/型别）
PYTHONUTF8=1 python -m decoy_a hla-filter --hla HLA-A*01:01

# 自定义肽段池
python -m decoy_a scan --target EVDPIGHLY --pool my_peptides.csv
```

### Decoy A + B 全管线

```bash
python -m decoy_b run \
  --target EVDPIGHLY \
  --hla HLA-A*02:01 \
  --protein MAGE-A3 \
  --top-n 100
```

### 检查外部工具状态

```python
from decoy_a.tools.mhcflurry import check_available as mhcf_ok
from decoy_a.tools.netmhcpan import check_available as nmhc_ok
from decoy_b.tools.tfold import check_available as tf_ok
from decoy_b.tools.alphafold3 import check_available as af3_ok
from decoy_b.tools.proteinmpnn import check_available as mpnn_ok

print(f"mhcflurry:    {mhcf_ok()}")   # Decoy A HLA 过滤（已部署）
print(f"NetMHCpan:    {nmhc_ok()}")   # Decoy A HLA 过滤（可选升级）
print(f"tFold:        {tf_ok()}")     # Decoy B 结构预测
print(f"AlphaFold3:   {af3_ok()}")    # Decoy B 结构精筛
print(f"ProteinMPNN:  {mpnn_ok()}")   # Decoy B 逆向设计
```

---

## TODO: 外部工具部署

Decoy A 已完全部署（mhcflurry 作为 HLA 预测后端）。以下为 Decoy B 所需的外部工具：

### 已完成: mhcflurry 2.2.0 — HLA 呈递过滤

**状态**: 已安装并验证。HLA-A\*02:01 过滤已完成（1,280,221 binders）。

### 可选升级: NetMHCpan 4.1 — 更精确的 HLA 预测

**用于**: 替代 mhcflurry 获得更精确的 EL %Rank（Eluted Ligand 模型）

| 项目 | 内容 |
|------|------|
| 申请 | https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/ (学术免费) |
| 下载 | `netMHCpan-4.1.Linux.tar.gz` + `data.tar.gz` |
| 当前状态 | mhcflurry 已作为替代方案部署，功能等价 |

```bash
export NETMHCPAN_DIR=~/tools/netMHCpan-4.1
```

### P0: tFold — pMHC 结构预测

**用于**: Decoy B Stage 2（核心）

| 项目 | 内容 |
|------|------|
| 代码 | `pip install tfold` 或 `git clone https://github.com/TencentAI4S/tfold.git` |
| 权重 | Zenodo: https://zenodo.org/records/12602915 |
| 无该工具时 | Decoy B 只能用 Atchley 理化特征，无法做 3D 结构比对 |

```
~/tools/tfold/
└── weights/
    ├── esm_ppi_650m_tcr.pt
    ├── tfold_tcr_trunk.pt
    └── tfold_tcr_pmhc_trunk.pt
```

```bash
export TFOLD_DIR=~/tools/tfold
export TFOLD_WEIGHTS_DIR=~/tools/tfold/weights
```

### P1: ProteinMPNN — 逆向序列设计

**用于**: Decoy B MPNN 分支

| 项目 | 内容 |
|------|------|
| 代码 | `git clone https://github.com/dauparas/ProteinMPNN.git` |
| 权重 | **已随仓库提供**（`vanilla_model_weights/v_48_020.pt`） |
| 依赖 | `pip install torch` |
| 无该工具时 | 缺少"主动生成"能力，仅靠搜索 |

```
~/tools/ProteinMPNN/
├── protein_mpnn_run.py
├── helper_scripts/
└── vanilla_model_weights/    ← 权重已包含
    └── v_48_020.pt
```

```bash
export PROTEINMPNN_DIR=~/tools/ProteinMPNN
```

### P2: AlphaFold3 — 结构精筛（可选）

**用于**: Decoy B Stage 3（Top 200 精细验证）

| 项目 | 内容 |
|------|------|
| 代码 | `git clone https://github.com/google-deepmind/alphafold3.git` |
| 权重 | 申请: https://forms.gle/svvpY4u2jsHEwWYS6 (约 2-3 天) |
| 数据库 | `bash fetch_databases.sh` (~630 GB) |
| 无该工具时 | tFold 结果已足够，AF3 仅提升 Top 200 精度 |

```
~/tools/alphafold3/
├── models/        ← 权重
├── databases/     ← 遗传数据库 (630 GB)
└── run_alphafold.py
```

```bash
export AF3_DIR=~/tools/alphafold3
export AF3_MODEL_DIR=~/tools/alphafold3/models
export AF3_DB_DIR=~/tools/alphafold3/databases
```

### 环境变量汇总（~/.bashrc）

```bash
export NETMHCPAN_DIR=~/tools/netMHCpan-4.1
export TFOLD_DIR=~/tools/tfold
export TFOLD_WEIGHTS_DIR=~/tools/tfold/weights
export PROTEINMPNN_DIR=~/tools/ProteinMPNN
export AF3_DIR=~/tools/alphafold3
export AF3_MODEL_DIR=~/tools/alphafold3/models
export AF3_DB_DIR=~/tools/alphafold3/databases
```

---

## 可视化

每条管线均配有可视化脚本，运行即可生成全套图表：

```bash
python scripts/visualize_decoy_a.py       # Decoy A: 4 张 matplotlib + 1 个 Plotly HTML
python scripts/visualize_decoy_c.py       # Decoy C: 4 张 matplotlib + 1 个 Plotly HTML
python scripts/visualize_presentation.py  # Presentation Score 对比: 5 面板 matplotlib
```

| 文件 | 内容 |
|------|------|
| `figures/decoy_a_overview.png` | EL%Rank 分布、Binder 比例、肽段长度、各长度 violin |
| `figures/decoy_a_funnel.png` | Pipeline 生成 + 过滤漏斗 (含 pass rate) |
| `figures/decoy_a_scatter.png` | EL%Rank vs Hamming distance 散点图 |
| `figures/decoy_a_mismatch_positions.png` | 突变位点分布 (TCR/Anchor/Other) |
| `figures/decoy_a_interactive.html` | Plotly: EPS8L2 案例 + 50 组织表达热图 |
| `figures/decoy_a_presentation_rescore.png` | Affinity vs Presentation 重新分类对比 (5 面板) |
| `figures/decoy_c_pipeline.png` | 4-stage 管线流程图 |
| `figures/decoy_c_overview.png` | 证据等级、验证状态、HLA 分布、受累器官 |
| `figures/decoy_c_filtering.png` | 筛选漏斗 (~200 papers → 5 fatal) |
| `figures/decoy_c_assay_year.png` | 实验方法频次 + 发表年份时间线 |
| `figures/decoy_c_interactive.html` | Plotly: DC-0001 Case Study + MAGE-A3 网络 + 旭日图 |

---

## 项目结构

```
decoy_library/
│
├── README.md
├── .env.example
│
├── decoy_a/                         # Decoy A: 序列同源扫描 + 共享基座
│   ├── PROGRESS_REPORT.md           #   技术进展报告
│   ├── config.py                    #   路径、阈值、Atchley 因子
│   ├── models.py                    #   14 个 Pydantic v2 数据模型（A/B 共用）
│   ├── kmer_builder.py              #   Step 1: Swiss-Prot → 41.7M k-mer + 表达谱
│   ├── hla_filter.py                #   Step 2: HLA 呈递门控（mhcflurry/NetMHCpan/IEDB）
│   ├── scanner.py                   #   Step 3A: Hamming 序列同源扫描
│   ├── run.py                       #   独立运行器（demo / 自定义 / 全管线）
│   └── tools/
│       ├── mhcflurry.py             #   mhcflurry 2.2.0 封装（已部署）
│       └── netmhcpan.py             #   NetMHCpan 4.1 封装（可选）
│
├── decoy_b/                         # Decoy B: 结构相似 + 逆向设计
│   ├── PROGRESS_REPORT.md           #   技术进展报告
│   ├── DEPLOYMENT.md                #   部署指南
│   ├── scanner.py                   #   Step 3B: Atchley → tFold → AF3 → 比对
│   ├── risk_scorer.py               #   Step 4: 统一风险评分（合并 A+B）
│   ├── orchestrator.py              #   全管线编排（Step 1-4 + 断点续传）
│   ├── main.py                      #   CLI（python -m decoy_b）
│   └── tools/
│       ├── tfold.py                 #   tFold pMHC 结构预测封装
│       ├── alphafold3.py            #   AlphaFold3 封装
│       └── proteinmpnn.py           #   ProteinMPNN 逆向设计封装
│
├── decoy_c/                         # Decoy C: 文献挖掘
│   ├── PROGRESS_REPORT.md           #   技术进展报告
│   ├── config.py                    #   API 配置（NCBI/UniProt/IEDB/LLM）
│   ├── models.py                    #   Pydantic v2 schema (DecoyEntry)
│   ├── fetcher.py                   #   PubMed/PMC 文献获取
│   ├── multi_source_fetcher.py      #   arXiv/bioRxiv/ClinicalTrials.gov 等 6 源
│   ├── extractor.py                 #   LLM 结构化信息提取
│   ├── validator.py                 #   UniProt + IEDB 交叉验证
│   ├── orchestrator.py              #   管线编排
│   ├── seed_data.py                 #   51 条手工标注种子
│   ├── scale_up.py                  #   自动扩量（6 策略 / 无限模式）
│   └── main.py                      #   CLI（python -m decoy_c）
│
├── scripts/                         # 工具脚本
│   ├── visualize_decoy_a.py         #   Decoy A 可视化 (matplotlib + plotly)
│   ├── visualize_decoy_c.py         #   Decoy C 可视化 (matplotlib + plotly)
│   ├── visualize_presentation.py    #   Presentation Score 对比可视化
│   ├── rescore_presentation.py      #   Presentation Score 重新评分脚本
│   ├── setup_decoy_b.py             #   Decoy B 外部工具自动安装
│   └── verify_decoy_b.py            #   Decoy B 工具可用性验证
│
├── figures/                         # 生成的可视化图表
│   ├── decoy_a_*.png / .html        #   Decoy A 可视化 (6 文件)
│   ├── decoy_c_*.png / .html        #   Decoy C 可视化 (5 文件)
│   └── decoy_a_presentation_rescore.png  # Presentation Score 对比
│
└── data/
    ├── decoy_a/                     # Decoy A 数据（共享基座 + A 输出）
    │   ├── human_kmer_db.parquet    #   41.7M 唯一人类短肽（380 MB, git-ignored）
    │   ├── gene_expression.parquet  #   20,151 基因 × 50 组织（1.9 MB, git-ignored）
    │   ├── hla_filtered_HLA-A0201.parquet           # 1.28M binders, affinity（27 MB）
    │   ├── hla_filtered_HLA-A0201_presentation.parquet  # 同上 + presentation score（48 MB）
    │   └── decoy_a_results.json     #   最近一次扫描结果
    ├── decoy_b/                     # Decoy B 数据（结构模型 + 排名结果，待生成）
    └── decoy_c/                     # Decoy C 数据（文献库）
        ├── decoy_library.json       #   主库（250 条，105 条有证据等级）
        ├── seed_entries.json        #   51 条手工标注种子
        ├── cross_reactivity_reference.json  #   288 条广谱交叉反应参考
        └── removed_entries.json     #   46 条已排除条目
```

---

## 当前进度总览

| 管线 | 状态 | 说明 |
|------|------|------|
| **Decoy A** | ✅ 已完成 | 全管线部署 + Presentation Score 升级 + 可视化 |
| **Decoy B** | ⏳ 代码完成 | Atchley 物化筛选可用；tFold/AF3 待部署 GPU 工具 |
| **Decoy C** | ✅ 已完成 | 250 条入库 (105 有证据等级) + 可视化 |
| **可视化** | ✅ 已完成 | A + C + Presentation Score 全套图表 |

## 计算资源

| 步骤 | 工具 | GPU | 耗时 | 状态 |
|------|------|-----|------|------|
| K-mer 构建 | SQLite + Python | 否 | ~45 分钟 | ✅ 已完成 |
| 表达谱构建 | HPA 数据 | 否 | ~2 分钟 | ✅ 已完成 |
| HLA 过滤 (affinity) | mhcflurry AffinityPredictor | 否 | ~6 小时 | ✅ HLA-A\*02:01 已完成 |
| HLA 重评 (presentation) | mhcflurry PresentationPredictor | 否 | ~16 分钟 | ✅ HLA-A\*02:01 已完成 |
| **Decoy A 扫描** | **NumPy** | **否** | **~6 秒** | **✅ 可用** |
| Decoy B Atchley 初筛 | NumPy | 否 | 1-5 分钟 | ✅ 可用 |
| Decoy B tFold 批量 | A100/V100 | 是 | 1-2 天 (5K 条) | ⏳ 待部署 tFold |
| Decoy B AF3 精筛 | A100 80GB | 是 | 7-17 小时 (200 条) | ⏳ 待部署 AF3 |
| Decoy B MPNN 设计 | 任意 GPU / CPU | 可选 | ~10 分钟 | ⏳ 待部署 MPNN |
| Decoy C 文献挖掘 | LLM API | 否 | 持续运行 | ✅ 可用 |

## 参考文献

- Cameron et al. (2013) Sci Transl Med. PMID: 23926201 — MAGE-A3/Titin 心脏毒性致死
- Linette et al. (2013) Blood. PMID: 23863783 — MAGE-A3 TCR 脑毒性
- Morgan et al. (2013) J Immunother. PMID: 24475783 — CEA TCR 结肠毒性
- Raman et al. (2016) J Biol Chem. PMID: 26457759 — TCR 交叉反应结构基础
- Parkhurst et al. (2011) Mol Ther. PMID: 21282551 — DMF5 TCR / MART-1 黑色素细胞毒性
