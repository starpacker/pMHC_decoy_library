# Decoy Library — TCR 脱靶毒性 pMHC 负样本库

## 问题

工程化 TCR-T 疗法的核心安全风险是**交叉反应性**：TCR 错误识别正常组织细胞表面的 pMHC，攻击健康组织，导致严重不良事件甚至死亡。

**Decoy Library 的目标**：构建一个全面的负样本库——收集所有可能被 TCR 错误识别的 pMHC，用于临床前安全性筛查和 IND 申报。所有的计算资源、数据、以及外部依赖库均被严格限制部署在 `/share/liuyutian` 目录下。

## 四层互补体系

```
                        ┌─────────────────────────────┐
                        │       Decoy Library          │
                        │   TCR 脱靶毒性 pMHC 负样本库   │
                        └──────┬──────┬──────┬──────┬──┘
                               │      │      │      │
               ┌───────────────┘      │      │      └───────────────┐
               ▼                      ▼      │                      ▼
     ┌─────────────────┐   ┌─────────────────┤   ┌─────────────────┐
     │    Decoy A       │   │    Decoy B       │   │    Decoy C       │
     │   序列同源扫描    │   │  结构相似筛选     │   │   文献挖掘        │
     │                  │   │                  │   │                  │
     │  "长得像的"       │   │  "形状像的"       │   │ "文献报道过的"     │
     │  Hamming ≤ 4     │   │  双叠合 RMSD      │   │  临床/实验证据     │
     │  50-200 条/靶标  │   │  + 5维界面描述符   │   │  105 条（已整理）  │
     │                  │   │  100-500 条/靶标  │   │                  │
     │  确定性: 高       │   │  确定性: 中       │   │  确定性: 最高      │
     │  覆盖面: 窄       │   │  覆盖面: 广       │   │  覆盖面: 受限      │
     └─────────────────┘   └─────────────────┘   └─────────────────┘
                                     │
                             ┌───────▼─────────┐
                             │    Decoy D       │
                             │  MPNN 逆向设计    │
                             │                  │
                             │  "AI 生成的"      │
                             │  固定锚定位设计    │
                             │  mhcflurry 过滤   │
                             │                  │
                             │  确定性: 理论风险   │
                             │  覆盖面: 预测未来   │
                             └─────────────────┘
```

---

## 使用方法

```bash
# Decoy A — 序列同源扫描
python -m decoy_a scan --target GILGFVFTL --hla "HLA-A*02:01"

# Decoy B — 结构相似筛选（完整管线）
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# Decoy B + MPNN 逆向设计
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --mpnn

# Decoy C — 文献挖掘
python -m decoy_c stats
python -m decoy_c show
```

---

## 风险评分模型

### Decoy A
$$ Risk_A = \left(1 - \frac{\text{Hamming}}{\text{Length}}\right) \times \frac{1}{EL\_Rank} \times TPM\_Weight $$

### Decoy B
$$ Combined = 0.4 \times CosSim + 0.6 \times StructSim + CV\_Boost $$

其中 `CV_Boost = 0.10 × cross_validation_agreement`（Boltz-2 交叉验证一致性加成，最高 +10%）。

**StructSim** 由双叠合 RMSD 和 5 维界面描述符综合得出：

$$StructSim = 0.50 \times RMSD\_Geo + 0.50 \times Interface\_Combined$$

| 组件 | 计算方式 |
|------|---------|
| **RMSD_Geo** | avg(Method A: MHC叠合→肽段RMSD, Method B: 肽段叠合→groove RMSD) |
| **Interface_Combined** | 0.25×PLIP + 0.10×BSA + 0.20×PRODIGY + 0.25×ESP + 0.20×PeSTo |

5 个界面描述符：

| 描述符 | 度量方式 | 权重 | 捕捉维度 |
|--------|---------|------|---------|
| PLIP | 非共价相互作用 Tanimoto | 0.25 | 侧链 H-bond/疏水/盐桥 |
| APBS/ESP | 静电势 Pearson 相关 | 0.25 | 电荷互补性 (molecular mimicry 核心) |
| PRODIGY | 结合亲和力 ΔG 差值 | 0.20 | 热力学稳定性 |
| PeSTo | 界面嵌入 Cosine 相似度 | 0.20 | learned 界面兼容性 |
| FreeSASA/BSA | 埋藏面积差值 | 0.10 | 界面几何 |

### 通用权重因子
- **EL_Rank**: mhcflurry presentation_percentile，越小亲和力越强
- **TPM_Weight**: 心脏/大脑等致命器官 ($TPM > 10$) 有 $10\times$ 加权惩罚

---

## 环境配置

### Conda 环境

本项目**仅使用一个 conda 环境**运行所有管线：

| 环境名 | 路径 | 用途 |
|--------|------|------|
| **`base`** | `/home/liuyutian/server/miniconda3` | **唯一运行环境** |

### 关键 Python 包

| 包名 | 用途 |
|------|------|
| `mhcflurry` | HLA 呈递预测 |
| `torch` (PyTorch) | tFold/ProteinMPNN/PeSTo 推理 |
| `biopython` | PDB 解析、结构叠合 |
| `numpy`, `scipy` | 数值计算 |
| `pandas`, `pyarrow` | 数据处理 |
| `pydantic` | 数据模型验证 |
| `plip` | PLIP 非共价相互作用分析 |
| `freesasa` | 溶剂可及表面积 |
| `prodigy-prot` | 接触基结合亲和力预测 |

### 外部工具

所有外部依赖部署在 `/share/liuyutian`，路径统一在 `decoy_a/config.py` 管理：

| 工具 | 用途 | 加载方式 |
|------|------|---------|
| **tFold** | pMHC 结构预测 | vendored: `decoy_b/external/tfold/` |
| **ProteinMPNN** | 逆向序列设计 | vendored: `decoy_d/external/proteinmpnn/` |
| **AlphaFold 3** | 精细结构验证 | subprocess 调用 |
| **Boltz-2** | 结构交叉验证 | CLI / Python API |
| **pdb2pqr + APBS** | 静电势计算 | subprocess (有 BioPython fallback) |
| **PeSTo** | 界面嵌入 | PyTorch (有 contact proxy fallback) |

### 环境变量

```bash
export MHCFLURRY_DATA_DIR="/share/liuyutian/mhcflurry_data/4"
export PYTHONUTF8=1
export TFOLD_DIR="<project>/decoy_b/external/tfold"
```

---

## 项目结构

```
pMHC_decoy_library/
├── decoy_a/                         # Decoy A: 序列同源扫描 + 共享基座
│   ├── config.py                    #   全局配置
│   ├── models.py                    #   Pydantic v2 数据模型（A/B/D 共用）
│   ├── kmer_builder.py              #   Swiss-Prot → 41.7M k-mer
│   ├── hla_filter.py                #   HLA 呈递门控
│   ├── scanner.py                   #   Hamming 序列同源扫描
│   └── tools/                       #   mhcflurry + NetMHCpan
│
├── decoy_b/                         # Decoy B: 结构相似筛选
│   ├── scanner.py                   #   五阶段管线 + 双叠合 + 界面描述符
│   ├── risk_scorer.py               #   统一风险评分（合并 A+B）
│   ├── orchestrator.py              #   全管线编排 + 断点续传
│   ├── tools/
│   │   ├── interface_descriptors.py #   5-descriptor 界面描述符
│   │   ├── tfold.py                 #   tFold 封装
│   │   ├── alphafold3.py            #   AF3 封装
│   │   └── proteinmpnn.py           #   ProteinMPNN 封装
│   └── DEPLOYMENT.md                #   外部工具部署指南
│
├── decoy_c/                         # Decoy C: 文献挖掘
│   ├── extractor.py                 #   LLM 结构化提取
│   ├── validator.py                 #   UniProt + IEDB 验证
│   └── ...
│
├── decoy_d/                         # Decoy D: MPNN 逆向设计
│   └── scanner.py                   #   MPNN + mhcflurry + tFold
│
├── scripts/                         # 工具脚本
│   ├── test_descriptors_on_real_data.py  # 5-descriptor 真实数据验证
│   ├── benchmark_interface_descriptors.py
│   ├── visualize_*.py
│   └── ...
│
├── tests/                           # 测试
│   └── test_interface_descriptors.py  # 46 tests (36 pass + 10 integration skip)
│
├── data/                            # 运行数据与结果
├── figures/                         # 可视化图表
├── docs/                            # 补充文档
│   └── deep-research_decoy_b.md     # 文献综述 (界面描述符选型依据)
├── progress_and_report.md           # 详细技术报告与进展
└── README.md
```

---

## 当前进度

| 管线 | 状态 | 说明 |
|------|------|------|
| **Decoy A** | ✅ 完成 | 全管线 + Presentation Score 升级 |
| **Decoy B** | ✅ 完成 | 5-descriptor 界面评分 + Boltz 交叉验证 |
| **Decoy C** | ✅ 完成 | 250 条入库 (105 有证据等级) |
| **Decoy D** | ✅ 完成 | MPNN 逆向设计 + mhcflurry 过滤 |

详细技术报告见 → [`progress_and_report.md`](progress_and_report.md)
