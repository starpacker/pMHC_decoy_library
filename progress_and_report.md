# Decoy Library — 技术进展与详细报告

**最后更新**: 2026-04-03
**项目**: TCR 脱靶毒性 pMHC 负样本库 (Decoy Library)

---

## 目录

1. [架构重构与 Self-Contained 改造](#1-架构重构与-self-contained-改造)
2. [核心算法修复](#2-核心算法修复)
3. [模型部署状态](#3-模型部署状态)
4. [靶标测试结果 (GILGFVFTL)](#4-靶标测试结果-gilgfvftl)
5. [双叠合结构比较方法](#5-双叠合结构比较方法)
6. [界面描述符升级 (5-Descriptor)](#6-界面描述符升级-5-descriptor)
7. [界面描述符实验结果](#7-界面描述符实验结果)
8. [全管线 Bug 修复记录](#8-全管线-bug-修复记录)
9. [部署教程 (Linux Server)](#9-部署教程-linux-server)

---

## 1. 架构重构与 Self-Contained 改造

### 1.1 代码内置化

- 将 `tfold` 源码完整迁移至 `decoy_b/external/tfold`
- 将 `ProteinMPNN` 源码迁移至 `decoy_b/external/proteinmpnn`
- 通过 `sys.path` 注入 vendored 代码，**彻底移除对外部环境的依赖**

### 1.2 配置解耦

所有权重路径与工具路径收敛到 `decoy_a/config.py`。使用者只需修改一个配置文件即可完成整个管线的迁移。

---

## 2. 核心算法修复

### 2.1 候选池扩充：解除同长度限制

**问题**: 原 pipeline 存在 `len(candidate) == len(target)` 硬限制，丢失 40.9% 的 8/10/11-mer 候选。

**修复**: Atchley TCR 接触核心（中央 5 残基，25 维）天然支持跨长度余弦相似度。候选池从 75.7 万 9-mer 扩大到 **128 万条 8-11mer**。

### 2.2 RMSD 算法修正

**问题**: 原 RMSD 在 310 个 CA 原子（MHC+β2m+肽段）上计算，MHC 保守性稀释了肽段差异，所有打分集中在 ~0.15Å。

**修复**:
1. 全局对齐：MHC 骨架 (Chain M+N, 301 CA) 做 Superimpose
2. 局部计算：同长度肽段算全长 CA RMSD；跨长度取中央 5 CA RMSD

修复后 RMSD 分布 0.18–3.0+ Å，获得有效区分度。

---

## 3. 模型部署状态

### 3.1 数据源

| 文件 | 总量 |
|------|------|
| `hla_filtered_HLA-A0201_presentation.parquet` | 1,280,221 条 (EL%Rank ≤ 2.0) |
| 长度分布 | 8-mer 1.5% / 9-mer 59.1% / 10-mer 28.5% / 11-mer 10.8% |

### 3.2 预测模型

| 工具 | 状态 |
|------|------|
| tFold | 部署成功，~2s/结构 |
| ProteinMPNN | 权重挂载成功 |
| AlphaFold 3 | 权重存在，留作精筛 |
| Boltz-2 | 待部署 |

---

## 4. 靶标测试结果 (GILGFVFTL)

使用流感 M1 经典靶标 **GILGFVFTL / HLA-A\*02:01** 进行全链路测试。

### 4.1 Decoy A

Hamming ≤ 2 在同长度 75.7 万候选中命中 **2 条** (TAP2: GLLGFVGTL, PKD1L1: GILLFLFTL)。

### 4.2 Decoy B Top 5

| Rank | Sequence | HD | CosSim | RMSD (Å) | Score | Gene |
|------|----------|----|--------|-----------|-------|------|
| 1 | LVLGFVFML | 3 | 1.000 | 0.179 | 0.989 | SLC39A9 |
| 2 | AYLGFVFYL | 3 | 1.000 | 0.270 | 0.974 | SLC11A2 |
| 3 | ASLGFVFSA | 4 | 1.000 | 0.500 | 0.900 | SLC16A11 |
| 4 | FLLGFVIMP | 5 | 0.987 | 0.664 | 0.862 | TAAR3P |
| 5 | LVLGFVFMLL | ×10 | 1.000 | 2.148\* | 0.698 | SLC39A9 |

\*跨长度 10-mer 与 9-mer 中央 5 CA RMSD。

**结论**: LVLGFVFML (Hamming=3) 与靶标 TCR 接触核心 `LGFVF` 极度相似，肽段 RMSD 仅 0.179Å，属传统安全盲区中的高危候选。

---

## 5. 双叠合结构比较方法

### 5.1 方法

| | Method A | Method B |
|--|---------|---------|
| **叠合基准** | MHC 骨架 (Chain M+N, 301 CA) | 肽段骨架 (Chain P CA) |
| **测量对象** | 肽段 RMSD | MHC groove 螺旋 RMSD (α1+α2) |
| **物理含义** | 在同一 HLA 凹槽中，肽段构象偏差 | 在同一肽段构象下，groove 包裹偏差 |

最终 StructSim = 两种方法 similarity 的平均值。

### 5.2 验证结果 (50 个 decoy)

- Pearson r = **0.9514**, Spearman ρ = **0.8955**
- Rank 31-50 **完全一致**

### 5.3 有意义的排名分歧

**A好B差** (肽段近似但 groove 偏差大): YLLGFLFNY (+12), LLLGILFLI (+11)
**B好A差** (groove 保守但肽段位移大): RALGFLIGL (-12), KTLGIVIGV (-11)

**结论**: 取平均值可同时捕获两类风险。

---

## 6. 界面描述符升级 (5-Descriptor)

### 6.1 Motivation

Cα RMSD 的三个盲区：(1) 侧链相互作用不可见 (2) 埋藏面积未捕获 (3) 无结合亲和力代理。

文献确认序列相似性对识别结构 mimics **几乎没有预测价值** (MimicryDB-Auto)。cross-reactivity 更好地由界面相互作用 profile 和热力学签名预测。

### 6.2 五个描述符

| # | 描述符 | 方法 | 相似度度量 | 权重 |
|---|--------|------|-----------|------|
| 1 | **PLIP** | 非共价相互作用指纹 (H-bond, 疏水, 盐桥, π-stacking, cation-π) | Generalized Tanimoto $T(a,b) = \frac{a \cdot b}{\|a\|^2 + \|b\|^2 - a \cdot b}$ | 0.25 |
| 2 | **FreeSASA** | BSA = SASA(pep) + SASA(MHC) - SASA(complex) | $1 - \|BSA_{tgt} - BSA_{cand}\| / 800$ | 0.10 |
| 3 | **PRODIGY** | 接触极性分类 → ΔG 线性模型 | $1 - \|\Delta G_{tgt} - \Delta G_{cand}\| / 5.0$ | 0.20 |
| 4 | **APBS/ESP** | 静电表面电位 (Poisson-Boltzmann / Coulomb proxy) | Pearson 相关映射 [-1,1]→[0,1] | 0.25 |
| 5 | **PeSTo** | 几何深度学习界面嵌入 / contact-propensity proxy | Cosine 相似度映射 [-1,1]→[0,1] | 0.20 |

### 6.3 Adaptive Weighting

如果某个描述符计算失败（工具不可用或 PDB 解析错误），其权重**自动重分配**给成功的描述符：

$$interface\_combined = \frac{\sum_{i \in active} w_i \cdot s_i}{\sum_{i \in active} w_i}$$

### 6.4 评分公式

```
surface_correlation = 0.50 × rmsd_geo + 0.50 × interface_combined

rmsd_geo = max(avg(sim_A, sim_B), bf_corr)            # unchanged
interface_combined = adaptive_weighted(
    0.25×PLIP + 0.10×BSA + 0.20×PRODIGY + 0.25×ESP + 0.20×PeSTo
)

similarity_score = 0.4 × cos_sim + 0.6 × surface_correlation
risk = similarity_score × (1/EL_Rank) × TPM_Weight
```

### 6.5 实现架构

```
decoy_b/tools/interface_descriptors.py (~900 lines)
├── InteractionFingerprint / InterfaceDescriptors / InterfaceSimilarity
├── Section 1: PLIP (compute_plip_fingerprint, compute_plip_tanimoto)
├── Section 2: FreeSASA (compute_bsa, compute_bsa_similarity)
├── Section 3: PRODIGY (compute_prodigy_affinity, compute_prodigy_similarity)
├── Section 4: APBS (compute_esp_vector → full APBS / Coulomb proxy)
├── Section 5: PeSTo (compute_pesto_embedding → full PeSTo / contact proxy)
└── Composite (compute_all_descriptors, compute_interface_similarity)

decoy_a/models.py
├── StructuralScore: +plip_tanimoto, bsa_*, prodigy_*, esp_similarity, pesto_similarity, interface_combined
└── DecoyABEntry: +plip_tanimoto, bsa_similarity, prodigy_similarity, esp_similarity, pesto_similarity, interface_combined

decoy_b/scanner.py: compute_structure_similarity() 集成全部 5 个描述符
decoy_b/risk_scorer.py: 传递新字段到最终输出
```

### 6.6 Proxy 实现

在没有 C 扩展的环境（如 Windows），每个描述符有 BioPython-only 的 fallback：

| 描述符 | Full | Proxy |
|--------|------|-------|
| PLIP | `plip` Python 包 | N/O < 3.5Å (H-bond), C-C < 4.5Å (hydrophobic) |
| FreeSASA | `freesasa` C 库 | 界面原子计数 (< 4.0Å) |
| PRODIGY | `prodigy-prot` | CA-CA 8.5Å 接触极性分类 + 线性模型 |
| APBS | pdb2pqr → APBS → DX grid | Coulomb: $\sum q_j/r_{ij}$ (AMBER 偏电荷) |
| PeSTo | PyTorch PeSTo 模型 | 接触密度 + 残基类型 + proximity → sigmoid interface prob |

### 6.7 测试状态

- **46 tests**: 36 passed + 10 skipped (integration tests 需要 PDB 文件)
- Unit tests 覆盖: Tanimoto (7), BSA similarity (6), PRODIGY similarity (5), ESP similarity (7), PeSTo similarity (7), InterfaceSimilarity (1), backward compatibility (3)

---

## 7. 界面描述符实验结果

### 7.1 数据集

27 个 9-mer GILGFVFTL decoy pMHC 结构 (tFold 预测)，基于 BioPython proxy 实现测试。

### 7.2 靶标 GILGFVFTL 描述符

| 描述符 | 值 |
|--------|-----|
| Contact fingerprint | HB=9, Hydro=137, SaltBr=4, Arom=8 |
| BSA proxy | 106 buried atoms |
| PRODIGY proxy | dG = -9.25 kcal/mol, 37 CA contacts |

### 7.3 关键发现

**RMSD 与界面描述符呈负相关**（r = -0.34 to -0.41），验证了多描述符方法的正交性：

| 相关对 | r |
|--------|---|
| RMSD vs PLIP Tanimoto | -0.34 |
| RMSD vs PRODIGY similarity | -0.41 |
| RMSD vs BSA similarity | -0.18 |

**最大排名变化**:

| 肽段 | 旧排名 | 新排名 | 变化 | 原因 |
|------|--------|--------|------|------|
| YLLGFLFNY | 3 | 27 | -24 | RMSD 好但 PRODIGY similarity 仅 0.178 |
| RALGFLIGL | 22 | 5 | +17 | RMSD 差但界面描述符全面匹配 |
| LLVGFVFVV | 1 | 10 | -9 | 几何完美但电荷/亲和力不匹配 |

**结论**: YLLGFLFNY 是几何上贴合但热力学上不同的"假阳性"；RALGFLIGL 是几何偏移但功能等价的"被遗漏的 mimic"。多描述符方法显著提高了筛选质量。

---

## 8. 全管线 Bug 修复记录

### 2026-04-02

#### Decoy A
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | Non_Binder 过滤增加 `binding` 列检查 | `scanner.py` |
| **严重** | 基因注释 LEFT JOIN → INNER JOIN + drop_duplicates | `hla_filter.py` |
| 中等 | mhcflurry 优先用 `Class1PresentationPredictor` | `tools/mhcflurry.py` |

#### Decoy B
| 严重性 | 修复 | 文件 |
|--------|------|------|
| 高 | BOTH 来源 risk 改为 0.5×seq + 0.5×physchem | `risk_scorer.py` |
| 高 | `get_gene_expression()` None 防护 | `scanner.py` |
| 中等 | EL rank 倒数增加 clamp (100.0) | `risk_scorer.py` |

#### Decoy C
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | PMID 字段格式化（非 PubMed 源返回 `source:id`） | `multi_source_fetcher.py` |
| **严重** | 去重改为 (sequence, hla_allele) 双键 | `models.py` |

#### Decoy D
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | anchor 位点改用 HLA-allele 特异查找表 | `scanner.py` |
| **严重** | `_parse_pdb()` 死代码删除 + 输出验证 | `tools/proteinmpnn.py` |

### 2026-04-02 晚 — 新增优化

1. **Decoy A 扩展**: Hamming ≤ 4 + Non_Binder 过滤，GILGFVFTL 命中 590+ 条
2. **Decoy D 独立模块**: MPNN 逆向设计从 Decoy B 分支独立为 Decoy D
3. **3D 可视化**: 新增 `visualize_3d.py` (3Dmol.js)
4. **全样本归档**: `data/GILGFVFTL_summary/` 统一归档 A/B/D 输出

---

## 9. 部署教程 (Linux Server)

### 9.1 基础环境

```bash
# 服务器路径
cd /share/liuyutian/pMHC_decoy_library

# Conda 环境 (使用 base)
conda activate base
```

### 9.2 界面描述符依赖安装

```bash
# PLIP (需要 OpenBabel)
conda install -c conda-forge openbabel -y
pip install plip

# FreeSASA
pip install freesasa

# PRODIGY
pip install prodigy-prot

# APBS + pdb2pqr
conda install -c conda-forge apbs pdb2pqr -y

# PeSTo (需要 PyTorch)
pip install pesto-protein  # 或从源码安装
```

### 9.3 验证安装

```bash
python -c "
from decoy_b.tools.interface_descriptors import (
    compute_plip_fingerprint, compute_bsa,
    compute_prodigy_affinity, compute_esp_vector,
    compute_pesto_embedding
)
print('All descriptors available')
"
```

### 9.4 运行完整管线

```bash
PYTHONUTF8=1 python -m decoy_b run \
    --target GILGFVFTL \
    --hla "HLA-A*02:01"
```

### 9.5 环境变量

```bash
export MHCFLURRY_DATA_DIR="/share/liuyutian/mhcflurry_data/4"
export PYTHONUTF8=1
export TFOLD_DIR="/share/liuyutian/pMHC_decoy_library/decoy_b/external/tfold"
```

详细的外部工具部署（tFold/AF3/Boltz/ProteinMPNN）见 → [`decoy_b/DEPLOYMENT.md`](decoy_b/DEPLOYMENT.md)

---

## 附录 A: 文件清单

### 根目录

| 文件 | 内容 |
|------|------|
| `README.md` | 项目总览 + 评分公式 + 环境配置 |
| `progress_and_report.md` | 本文件：详细技术报告 |

### 每个模块

| 模块 | README | 详细部署 |
|------|--------|---------|
| `decoy_a/` | `README.md` | — |
| `decoy_b/` | `README.md` | `DEPLOYMENT.md` |
| `decoy_c/` | `README.md` | — |
| `decoy_d/` | `README.md` | — |

### 参考文献

| 文件 | 内容 |
|------|------|
| `docs/deep-research_decoy_b.md` | 界面描述符文献综述（选型依据） |
| `docs/TUTORIAL_AF3_LINUX_INSTALL.md` | AlphaFold 3 安装教程 |
