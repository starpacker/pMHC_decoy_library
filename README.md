# Decoy Library — TCR 脱靶毒性 pMHC 负样本库

## 问题

工程化 TCR-T 疗法的核心安全风险是**交叉反应性**：TCR 错误识别正常组织细胞表面的 pMHC，攻击健康组织，导致严重不良事件甚至死亡。

**Decoy Library 的目标**：构建一个全面的负样本库——收集所有可能被 TCR 错误识别的 pMHC，用于临床前安全性筛查和 IND 申报。所有的计算资源、数据、以及外部依赖库均被严格限制部署在 `/share/liuyutian` 目录下。

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
     │   序列同源扫描    │   │  结构相似筛选     │   │   文献挖掘        │
     │                  │   │                  │   │                  │
     │  "长得像的"       │   │  "形状像的"       │   │ "文献报道过的"     │
     │  Hamming ≤ 2     │   │  肽段 RMSD 相似   │   │  临床/实验证据     │
     │  50-200 条/靶标  │   │  100-500 条/靶标  │   │  105 条（已整理）  │
     │                  │   │                  │   │                  │
     │  确定性: 高       │   │  确定性: 中       │   │  确定性: 最高      │
     │  覆盖面: 窄       │   │  覆盖面: 最广     │   │  覆盖面: 受限      │
     └─────────────────┘   └─────────────────┘   └─────────────────┘
             │                      │                      │
             └──────────────────────┼──────────────────────┘
                                    ▼
                            统一风险评分与排名
```

---

## 核心管线与使用方法

所有模块均通过标准的 Python 模块运行 (`python -m`)。

### Decoy A — 序列同源扫描

**功能**: 从人类 K-mer 数据库中搜索与靶标 Hamming 距离 $\le 2$ 且满足 HLA 呈递过滤的同源肽段。

**子命令** (`python -m decoy_a <command>`):
- `demo`: 快速演示（无需完整数据库，使用内置样本）
- `build-kmer`: Step 1: 构建 K-mer + 表达谱数据库
- `hla-filter`: Step 2: 运行 HLA 呈递过滤门控
- `scan`: 运行靶标同源扫描

**使用示例**:
```bash
# 快速演示模式
python -m decoy_a demo --target EVDPIGHLY

# 正式运行扫描 (默认靶标 EVDPIGHLY, HLA-A*02:01)
python -m decoy_a scan --target EVDPIGHLY --hla "HLA-A*02:01"
```

### Decoy B — 结构相似筛选

**功能**: 从 HLA 可呈递肽段库中筛选与靶标 TCR 接触面结构高度相似的肽段。**支持 8-11mer 全长度候选**（不仅限于与靶标同长度），采用四阶段管线：Atchley 理化因子初筛 → tFold 批量结构预测 → 肽段 RMSD 精筛 → 综合评分排名。

**四阶段管线**:
1. **Atchley 理化因子初筛**: 从 ~128 万 8-11mer 候选中，排除 Decoy A 领地 (同长度 HD≤2)，按 TCR 接触核心 (中央 5 残基) 的 Atchley 向量余弦相似度筛选 Top 5000
2. **tFold 批量结构预测**: 对 Top 候选批量生成 pMHC 三维结构 (PDB)，支持 8-11mer
3. **肽段 RMSD 精筛**: 以 MHC 重链+β2m (Chain M+N, 301 CA) 为对齐基准，计算肽段 RMSD：
   - **同长度**: 全肽段 (Chain P) CA 原子 RMSD
   - **跨长度**: TCR 接触核心 (中央 5 个 CA 原子) RMSD
4. **综合评分**: 融合 Atchley 余弦相似度 (40%) 与结构相似度 (60%)，输出最终排名

> **RMSD 计算方法**: 先在 MHC 骨架 (Chain M + Chain N, 共 301 个 CA 原子) 上做 Superimpose 对齐，然后计算肽段 RMSD。同长度时使用全肽段 CA 原子；**跨长度时仅比较 TCR 接触核心** (中央 5 个 CA 原子)，因为 TCR 识别主要依赖这些残基的空间构象。典型的肽段 RMSD 范围为 0.2–2.5 Å（同长度），2.0–3.0 Å（跨长度）。

**子命令** (`python -m decoy_b <command>`):
- `run`: 运行完整的 Decoy B 管线
- `build-kmer`: 构建 K-mer 库
- `hla-filter`: HLA 呈递过滤
- `scan-a`: 单独运行 Decoy A 扫描
- `scan-b`: 单独运行 Decoy B 结构扫描
- `score`: 对已有结果进行综合风险打分
- `show`: 显示排序结果 (支持 table / json)
- `stats`: 显示管线运行统计
- `clear`: 清除中间检查点，重新开始

**使用示例**:
```bash
# 运行完整的 Decoy B 管线
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# 仅运行 Decoy B 并在缺少结构预测环境时跳过 3D 建模
python -m decoy_b scan-b --target GILGFVFTL --skip-structural

# 显示排名前列的风险结果
python -m decoy_b show --format table
```

### Decoy C — 文献与临床证据挖掘

**功能**: 从 PubMed 等文献库自动化搜集已有报导的脱靶毒性 pMHC 数据。

**子命令** (`python -m decoy_c <command>`):
- `cold-start`: 初始化数据并从 PubMed 挖掘
- `search`: 搜索 PubMed 提取 Decoy
- `fetch`: 根据具体 PMID 提取
- `validate`: 验证库中记录的基因和蛋白信息
- `show`: 显示文献挖掘库内容
- `stats`: 显示文献库统计数据

**使用示例**:
```bash
# 查看已有数据统计
python -m decoy_c stats

# 显示记录内容
python -m decoy_c show
```

---

## 环境配置与外部工具部署

本项目已高度 **self-contained（自包含）**。核心的预测和设计代码（如 tFold 和 ProteinMPNN）均已复制并重构至 `decoy_b/external` 目录。
代码库外**只允许存在**大模型权重文件与个别独立容器（AlphaFold 3）。

所有的依赖和外部工具路径均统一在 `decoy_a/config.py` 中管理，使用者只需修改这一个配置文件即可指向自己的模型权重。

### 1. tFold — 批量结构预测 (Decoy B 必选)
tFold 的推理代码已经完全内置到本项目中 (`decoy_b/external/tfold`)。
- **模型检查点**: `esm_ppi_650m_tcr.pth` + `tfold_pmhc_trunk.pth`
- **权重位置**: 默认需要放置在 `~/.cache/torch/hub/checkpoints/`，或者通过 `decoy_a/config.py` 中的设置进行自定义配置。
- **输出 PDB 链约定**: Chain M = MHC 重链 (201 CA), Chain N = β2m (100 CA), Chain P = 肽段 (9 CA for 9-mer)
- **性能**: ~2s/结构 (GPU), 批量 5000 条约 3 小时

### 2. AlphaFold3 — 精细结构验证 (Decoy B 可选)
AF3 由于其庞大的模型权重与特殊的运行环境，**不在内置代码内**。
- **配置方式**: 通过 `decoy_a/config.py` 中的 `AF3_DIR` 指向您的本地 AF3 安装目录。默认指向 `/share/liuyutian/alphafold3`。
- **模型**: `af3.bin.zst`

### 3. ProteinMPNN — 逆向设计 (Decoy B 扩展)
ProteinMPNN 的推理代码已经内置到本项目中 (`decoy_b/external/proteinmpnn`)。
- **模型检查点**: 内置在 `vanilla_model_weights` 目录中。

### 验证工具可用性
你可以使用自带脚本检查外部工具的状态：
```bash
python scripts/verify_decoy_b.py
```

---

## 风险评分模型

### Decoy A 评分
$$ Risk_A = \left(1 - \frac{\text{Hamming}}{\text{Length}}\right) \times \frac{1}{EL\_Rank} \times TPM\_Weight $$

### Decoy B 评分
Decoy B 的综合评分融合理化相似度与结构相似度：

$$ Combined = 0.4 \times CosSim + 0.6 \times StructSim $$

其中结构相似度基于**肽段 RMSD**：

$$ StructSim = \max\left(0,\ 1 - \frac{PeptideRMSD}{3.0}\right) $$

- **CosSim**: TCR 接触核心 (p4-p8) 的 Atchley 向量余弦相似度
- **PeptideRMSD**: 在 MHC 骨架 (Chain M+N) 对齐后计算的肽段 RMSD (Å)。同长度→全肽段 CA；跨长度→TCR 接触核心 (中央 5 CA)
- **3.0 Å**: 归一化常数，RMSD ≥ 3.0 Å 视为完全不相似

### 通用权重因子
- **EL_Rank**: mhcflurry presentation_percentile，越小亲和力越强
- **TPM_Weight**: 表达谱权重。心脏、大脑等核心致命器官 ($TPM > 10$) 会有 $10\times$ 的加权惩罚

## 数据来源

### HLA 可呈递肽段库

Decoy A 和 Decoy B 共享同一份 HLA 可呈递肽段库，由 **mhcflurry** (`Class1PresentationPredictor`) 预计算生成。

- **数据文件**: `data/decoy_a/hla_filtered_HLA-A0201_presentation.parquet`
- **HLA-A\*02:01 总量**: 1,280,221 条可呈递肽段
- **长度分布**: 8-mer 19,592 (1.5%) / 9-mer 757,053 (59.1%) / 10-mer 364,839 (28.5%) / 11-mer 138,737 (10.8%)
- **筛选标准**: `presentation_percentile ≤ 2.0`（mhcflurry presentation score 前 2%）
- **包含字段**: `peptide`, `allele`, `presentation_score`, `presentation_percentile`, `processing_score`, `best_allele`, `best_transcript`, `gene_name`

> 加载优先级：系统会自动优先加载 `_presentation.parquet` 文件（mhcflurry presentation 模型），若不存在则回退到普通 `hla_filtered_*.parquet`（affinity 模型）。

## 输出与可视化

所有的运行结果与缓存均存放在 `data/` 目录下：
- `data/decoy_a/`: k-mer 数据库、表达谱、HLA 过滤缓存、**HLA 可呈递肽段 parquet 文件**。
- `data/decoy_b/`: tFold 生成的 PDB 结构、结构比对结果、最终排名 `final_ranked_decoys.json`。

可视化图表脚本存放于 `scripts/`：
```bash
python scripts/visualize_decoy_a.py
python scripts/visualize_decoy_c.py
python scripts/visualize_presentation.py
```
生成的图表保存在 `figures/` 目录下，包含特征散点图、漏斗图及 Plotly 交互式表达热图。

---

## Decoy B 逐步演示

`scripts/run_decoy_b_stepwise.py` 提供了一个完整的 Decoy B 四阶段管线演示，以 GILGFVFTL / HLA-A\*02:01 为例：

```bash
python scripts/run_decoy_b_stepwise.py
```

**输出目录**: `data/decoy_b/stepwise_demo/`

| 阶段 | 输出文件 | 说明 |
|------|---------|------|
| Stage 1 | `stage1_atchley_candidates.csv` | Atchley 理化因子初筛 Top 5000 |
| Stage 2 | `tfold_pdbs/*.pdb` | tFold 生成的 pMHC 三维结构 |
| Stage 3 | `stage3_structure_comparison.json` | 肽段 RMSD 比对结果 |
| Stage 4 | `final_decoy_b_results.json` | 最终综合评分排名 |

**示例结果** (GILGFVFTL / HLA-A\*02:01, Top 5):

| Rank | Sequence | Len | HD | CosSim | PepRMSD (Å) | Score | Gene |
|------|----------|-----|----|--------|-------------|-------|------|
| 1 | LVLGFVFML | 9 | 3 | 1.000 | 0.179 | 0.989 | SLC39A9 |
| 2 | AYLGFVFYL | 9 | 3 | 1.000 | 0.270 | 0.974 | SLC11A2 |
| 3 | ASLGFVFSA | 9 | 4 | 1.000 | 0.500 | 0.900 | SLC16A11 |
| 4 | FLLGFVIMP | 9 | 5 | 0.987 | 0.664 | 0.862 | TAAR3P |
| 5 | LVLGFVFMLL | 10 | ×10 | 1.000 | 2.148* | 0.698 | SLC39A9 |

\* 跨长度比较使用 TCR 接触核心 RMSD (中央 5 CA)

---

## MPNN 逆向设计分支 — 两层资格认证管线

### 原理

Decoy B 的主管线（Atchley → tFold → RMSD）是"从人类蛋白组中搜索"。MPNN 分支则反过来——**直接在靶标 pMHC 的 3D 骨架上设计全新序列**，然后验证这些序列是否真的危险。

ProteinMPNN 固定 HLA 锚定位（p2/p9），重新设计 TCR 接触面（p4-p8），生成 ~1000 条能完美折叠成靶标 3D 形状的候选序列。但"能折叠"不等于"能被呈递到细胞表面"——只有通过 HLA 呈递验证的序列才是真正的威胁。

### 两层资格认证

```
ProteinMPNN 生成 ~1000 条设计序列
        │
        ├── Tier 1: 精确匹配人类蛋白组 ──→ 最高置信度
        │   (设计序列恰好存在于 HLA 可呈递肽段库中)
        │
        └── Tier 2: mhcflurry HLA 呈递过滤 ──→ 理论风险
            (不在蛋白组中，但 mhcflurry 预测
             presentation_percentile ≤ 2.0)
            可能来源：体细胞突变 / 非经典 ORF /
            Swiss-Prot 未收录蛋白
```

| 层级 | 含义 | 置信度 | 典型数量 |
|------|------|--------|---------|
| **Tier 1** | 设计序列在人类蛋白组中精确匹配 | 最高 | 0-10 条 |
| **Tier 2** | 无蛋白组匹配，但 mhcflurry 预测可被 HLA 呈递 | 中等 | 50-300 条 |
| 淘汰 | mhcflurry 预测不可呈递 (percentile > 2.0) | — | ~600-900 条 |

### 数据模型

每个 MPNN 来源的 `DecoyBHit` 额外携带两个字段：

| 字段 | 类型 | 说明 |
|------|------|------|
| `mpnn_source` | `str` | `"proteome_matched"` (Tier 1) / `"hla_qualified_synthetic"` (Tier 2) / `""` (非 MPNN) |
| `mpnn_score` | `float` | ProteinMPNN 负对数概率分数（越低 = 折叠置信度越高） |

### 使用方法

MPNN 分支默认关闭，通过 `--mpnn` 标志启用：

```bash
# 完整管线 + MPNN 逆向设计
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --mpnn

# 仅 Decoy B 扫描 + MPNN（跳过 3D 建模）
python -m decoy_b scan-b --target GILGFVFTL --mpnn --skip-structural

# 完整管线 + MPNN + 3D 建模
python -m decoy_b run --target EVDPIGHLY --hla "HLA-A*02:01" --protein MAGE-A3 --mpnn
```

### 实现代码

| 文件 | 职责 |
|------|------|
| `decoy_b/tools/proteinmpnn.py` | ProteinMPNN 封装：PDB 解析 → 锚定位固定 → 序列设计 → FASTA 解析 |
| `decoy_b/scanner.py :: run_mpnn_design()` | 两层资格认证主逻辑：Tier 1 蛋白组匹配 + Tier 2 mhcflurry 过滤 |
| `decoy_b/scanner.py :: scan_decoy_b()` | 主管线末尾合并 MPNN 分支结果（去重 + 排序） |
| `decoy_a/tools/mhcflurry.py` | mhcflurry 封装：Tier 2 的 HLA 呈递预测后端 |
| `decoy_a/models.py :: DecoyBHit` | 数据模型：`mpnn_source` + `mpnn_score` 字段 |

### 前置条件

MPNN 分支需要：
1. **ProteinMPNN** 已安装（`/share/liuyutian/S3AI/rebuttal/ProteinMPNN`）
2. **靶标 PDB 结构**（由 tFold 或 AF3 生成，或手动提供）
3. **mhcflurry** 已安装（Tier 2 HLA 过滤所需，已部署）

---

## 项目结构

```
pMHC_decoy_library/
├── README.md
├── model_path.md                    # 外部模型路径说明
├── .env.example
│
├── decoy_a/                         # Decoy A: 序列同源扫描 + 共享基座
│   ├── config.py                    #   全局路径、阈值、Atchley 因子
│   ├── models.py                    #   Pydantic v2 数据模型（A/B 共用）
│   ├── kmer_builder.py              #   Step 1: Swiss-Prot → 41.7M k-mer + 表达谱
│   ├── hla_filter.py                #   Step 2: HLA 呈递门控（mhcflurry/NetMHCpan）
│   ├── scanner.py                   #   Step 3A: Hamming 序列同源扫描
│   ├── run.py                       #   独立运行器（demo / 自定义 / 全管线）
│   └── tools/
│       ├── mhcflurry.py             #   mhcflurry 2.2.0 封装（已部署）
│       └── netmhcpan.py             #   NetMHCpan 4.1 封装（可选）
│
├── decoy_b/                         # Decoy B: 结构相似 + 逆向设计
│   ├── scanner.py                   #   四阶段管线 + MPNN 两层资格认证分支
│   ├── risk_scorer.py               #   统一风险评分（合并 A+B）
│   ├── orchestrator.py              #   全管线编排（Step 1-4 + 断点续传）
│   ├── main.py                      #   CLI（python -m decoy_b）
│   └── tools/
│       ├── tfold.py                 #   tFold pMHC 结构预测封装
│       ├── alphafold3.py            #   AlphaFold3 封装（可选精筛）
│       └── proteinmpnn.py           #   ProteinMPNN 逆向设计 + HLA 过滤
│
├── decoy_c/                         # Decoy C: 文献挖掘
│   ├── fetcher.py                   #   PubMed/PMC 文献获取
│   ├── multi_source_fetcher.py      #   arXiv/bioRxiv/ClinicalTrials.gov 等
│   ├── extractor.py                 #   LLM 结构化信息提取
│   ├── validator.py                 #   UniProt + IEDB 交叉验证
│   ├── orchestrator.py              #   管线编排
│   ├── scale_up.py                  #   自动扩量（6 策略 / 无限模式）
│   └── main.py                      #   CLI（python -m decoy_c）
│
├── scripts/                         # 工具脚本
│   ├── run_decoy_b_stepwise.py      #   Decoy B 四阶段逐步演示
│   ├── build_demo_data.py           #   构建 demo 数据集
│   ├── setup_decoy_b.py             #   Decoy B 外部工具自动安装
│   ├── verify_decoy_b.py            #   Decoy B 工具可用性验证
│   ├── rescore_presentation.py      #   Presentation Score 重新评分
│   ├── generate_decoy_b_viewer.py   #   3D 结构对比查看器
│   ├── visualize_decoy_a.py         #   Decoy A 可视化
│   ├── visualize_decoy_c.py         #   Decoy C 可视化
│   └── visualize_presentation.py    #   Presentation Score 对比可视化
│
├── data/
│   ├── decoy_a/                     #   K-mer 库 + HLA 过滤 + 扫描结果
│   ├── decoy_b/                     #   PDB 结构 + RMSD 比对 + 排名结果
│   │   └── stepwise_demo/           #   逐步演示输出
│   └── decoy_c/                     #   文献库 (250 条)
│
└── figures/                         #   生成的可视化图表
    ├── decoy_a_*.png / .html        #   Decoy A 可视化 (6 文件)
    ├── decoy_b_viewer.html          #   Decoy B 3D 结构对比
    └── decoy_c_*.png / .html        #   Decoy C 可视化 (5 文件)
```

---

## 当前进度总览

| 管线 | 状态 | 说明 |
|------|------|------|
| **Decoy A** | ✅ 已完成 | 全管线部署 + Presentation Score 升级 + 可视化 |
| **Decoy B 主管线** | ✅ 已完成 | Atchley 初筛 + tFold 结构预测 + 肽段 RMSD 精筛 + 综合评分 |
| **Decoy B MPNN** | ✅ 已完成 | 两层资格认证（Tier 1 蛋白组匹配 + Tier 2 mhcflurry HLA 过滤） |
| **Decoy C** | ✅ 已完成 | 250 条入库 (105 有证据等级) + 可视化 |
| **可视化** | ✅ 已完成 | A + B + C 全套图表 |
