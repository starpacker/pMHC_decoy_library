# Decoy Library — 技术进展与详细报告

**最后更新**: 2026-04-07
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
8. [TCR-Facing 表面描述符重设计 (v2)](#8-tcr-facing-表面描述符重设计-v2)
9. [全管线 Bug 修复记录](#9-全管线-bug-修复记录)
10. [部署教程 (Linux Server)](#10-部署教程-linux-server)
11. [3D 结构可视化最佳实践](#11-3d-结构可视化最佳实践)
12. [PLIP 部署与集成](#12-plip-部署与集成)
13. [统一 CLI 入口与批量运行](#13-统一-cli-入口与批量运行)
14. [Decoy C 三重验证 + IEDB 采掘升级](#14-decoy-c-三重验证--iedb-采掘升级)

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

## 8. TCR-Facing 表面描述符重设计 (v2)

### 8.1 问题发现

通过对 GILGFVFTL 的 50 个候选 PDB 进行 **solo rank analysis**（每个描述符独立排名 → Spearman 相关矩阵），发现 v1 描述符存在严重问题：

**v1 Spearman 相关矩阵**:

| 对比 | ρ | 问题 |
|------|---|------|
| Atchley ↔ sim_A | -0.01 | **完全无关**，不应混入最终评分 |
| ESP Coulomb ↔ rmsd_geo | 0.43 | 噪声过大 (proxy 精度不足) |
| bf_corr ↔ rmsd_geo | 0.96 | **完全冗余** |
| ESP ↔ PeSTo | 0.13 | 两者都不可靠 |

**更根本的问题**: 所有 v1 描述符（PLIP、BSA、PRODIGY、ESP、PeSTo）度量的是 **peptide-MHC groove 结合界面**（peptide 埋在 groove 里的面），而 TCR 交叉反应取决于 **pMHC 暴露给 TCR 的上表面**，这是两个不同的界面。

### 8.2 v2 TCR-Facing 描述符

新描述符通过 **per-atom SASA > 0 过滤**，自动排除 groove 埋藏的锚定残基，只度量 TCR 实际接触的溶剂暴露面：

| 描述符 | 计算方法 | 度量 | 权重 |
|--------|---------|------|------|
| **TCR-facing ESP** | 暴露肽原子处的 Coulomb 电位 (来自周围所有原子) | Hodgkin SI | 0.30 |
| **Exposed Shape** | MHC 叠合后暴露侧链质心坐标 (Kabsch RMSD) | 1-RMSD/3 | 0.30 |
| **rSASA Profile** | per-position 相对溶剂可及面积 (complex/free) | Cosine sim | 0.20 |
| **Exposed Hydrophobicity** | rSASA 加权 Kyte-Doolittle 疏水性 | Cosine sim | 0.20 |

**关键技术**:
- **ShrakeRupley SASA** (BioPython): 区分溶剂暴露原子 vs groove 埋藏原子
- **Hodgkin SI**: `2·dot(A,B) / (dot(A,A)+dot(B,B))`，PIPSA/MatchTope 领域标准
- **Kabsch alignment**: SVD 最优旋转对齐暴露侧链质心
- **rSASA**: `SASA_complex(i) / SASA_free(i)`，验证: P4-P8 rSASA=0.36-0.76 (暴露), P2/P9 rSASA≈0 (锚定)

### 8.3 v2 验证结果

**v2 Spearman 相关矩阵** (50 candidates, GILGFVFTL):

| | tcr_esp | tcr_sasa | tcr_hydro | tcr_shape | tcr_combined | sim_A |
|---|---|---|---|---|---|---|
| **tcr_esp** | 1.00 | 0.65 | 0.33 | 0.44 | 0.67 | 0.43 |
| **tcr_sasa** | 0.65 | 1.00 | 0.40 | 0.31 | 0.59 | 0.72 |
| **tcr_hydro** | 0.33 | 0.40 | 1.00 | 0.02 | 0.43 | 0.15 |
| **tcr_shape** | 0.44 | 0.31 | 0.02 | 1.00 | 0.85 | 0.22 |
| **sim_A** | 0.43 | 0.72 | 0.15 | 0.22 | 0.38 | 1.00 |

**改进**:
- 4 个描述符之间无冗余（最高 ρ=0.65，v1 bf_corr↔rmsd 达 0.96）
- 每个描述符与 RMSD 的相关性适中 (0.15-0.72)，提供独立信息
- `tcr_shape` (暴露侧链形状) 是最有区分力的新维度，v1 完全缺失

### 8.4 v2 评分公式

```
similarity_score = surface_correlation    # 纯结构 (Atchley 移除)

surface_correlation = 0.50 * rmsd_geo + 0.50 * tcr_combined
rmsd_geo = avg(sim_A, sim_B)             # B-factor floor 移除
tcr_combined = 0.30*ESP + 0.30*Shape + 0.20*SASA + 0.20*Hydro
```

### 8.5 代码变更

| 文件 | 变更 |
|------|------|
| `decoy_b/tools/tcr_surface_descriptors.py` | **新增** — 4 个 TCR-facing 描述符 |
| `decoy_b/scanner.py` | 评分公式重写: 移除 Atchley/bf_corr，接入 TCR 描述符 |
| `decoy_b/tools/interface_descriptors.py` | 旧 groove 描述符，标记弃用 |
| `docs/deep-research_tcr-facing-descriptors.md` | 文献调研: MatchTope/PepSim/MaSIF/IMPRINT/dMaSIF |

### 8.6 参考文献

| 方法 | 论文 |
|------|------|
| MatchTope (PIPSA/MEP) | Mendes et al. Front Immunol 13:930590, 2022 |
| PepSim (表面+序列) | Hall-Swan et al. Front Immunol 14:1108303, 2023 |
| MaSIF (surface fingerprint) | Gainza et al. Nature Methods 17:184-192, 2020 |
| IMPRINT (immunological fingerprint) | Shang et al. Brief Bioinform 27(2):bbag048, 2026 |
| Titin/MAGE-A3 案例 | Cameron et al. Sci Transl Med 5:197ra103, 2013 |

---

## 9. 全管线 Bug 修复记录

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

## 10. 部署教程 (Linux Server)

### 9.1 基础环境

```bash
# 服务器路径
cd /share/liuyutian/pMHC_decoy_library

# Conda 环境 (使用 base)
conda activate base
```

### 10.2 TCR-Facing 描述符依赖

v2 TCR-facing 描述符仅依赖 BioPython (ShrakeRupley SASA)，**无额外安装**：

```bash
pip install biopython  # 已在 base 环境中
```

旧 v1 groove 描述符依赖（PLIP/FreeSASA/PRODIGY/APBS/PeSTo）已弃用，无需安装。

### 10.3 验证安装

```bash
python -c "
from decoy_b.tools.tcr_surface_descriptors import compute_tcr_facing_similarity
print('TCR-facing descriptors available')
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

## 11. 3D 结构可视化最佳实践

### 10.1 技术选型

| 方案 | 优点 | 缺点 | **结论** |
|------|------|------|---------|
| PyMOL 静态图 | 出版级质量 | 不可交互、需安装、不可分享 | 仅用于论文 figure |
| py3Dmol / Jupyter | 快速原型 | 绑定 notebook 环境 | 开发调试用 |
| **3Dmol.js 自包含 HTML** | 零依赖、可交互、浏览器直接打开、可嵌入报告 | 文件较大 (~3MB) | **Pipeline 标准输出** |
| NGL Viewer | 功能丰富 | API 复杂、加载慢 | 不推荐 |
| Mol\* (RCSB) | 专业级 | 过于重量级 | 不推荐 |

**决策**: 所有管线的结构可视化统一使用 **3Dmol.js 自包含 HTML**，PDB 数据内嵌 JSON string，无需后端服务器。

### 10.2 Decoy B/D 3D Viewer (`decoy_3d_comparison.html`)

#### 10.2.1 四种展示模式

| 模式 | 快捷键 | MHC 显示 | 背景色 | 肽段渲染 | 用途 |
|------|--------|---------|--------|---------|------|
| **Side-by-Side** | `1` | 半透明 cartoon | `#0d1117` (暗黑) | 绿/橙 Stick + Tube + Surface (0.28 opacity) | 默认浏览，比较整体构象 |
| **Superpose** | `2` | 仅 target 的 | `#0d1117` | 两条肽段叠合，mismatch 红色高亮 | **核心分析模式**，判断结构相似性 |
| **pLDDT Confidence** | `3` | **隐藏** | `#0d1117` | ROYGB 渐变 (红→黄→绿→蓝) + 半透明 Surface | 判断预测可信度 |
| **Surface Compare** | `4` | **隐藏** | `#2b2b2b` (深灰) | 按氨基酸电荷性质着色 + 不透明 Surface (0.88) | 比较电荷分布 |

**关键设计决策**:
- pLDDT 和 Surface 模式下 **MHC 完全隐藏**。MHC cartoon 颜色 (`#30363d`) 与暗色背景撞色，保留会产生"幽灵蛋白"干扰。
- Surface 模式背景改为 `#2b2b2b` 深灰（非纯黑），让白色中性残基可见。切回其他模式时自动恢复暗黑背景。

#### 10.2.2 颜色体系

**基础配色 (Side-by-Side / Superpose)**:

| 元素 | 颜色 | Hex | 渲染方式 |
|------|------|-----|---------|
| Target 肽段 | 绿色 | `#7ee787` | Stick (r=0.16) + Cartoon tube + Surface (0.28) |
| Decoy 肽段 | 橙色 | `#ffa657` | Stick (r=0.16) + Cartoon tube + Surface (0.28) |
| 突变残基 | 红色 | `#ff7b72` | 加粗 Stick (r=0.28)，覆盖原色 |
| MHC groove | 深灰 | `#30363d` | Cartoon, opacity 0.20 |
| β2m | 深灰 | `#30363d` | Cartoon, opacity 0.10 |

**pLDDT 配色** (ROYGB 渐变，映射 B-factor 0.5→1.0):

| pLDDT 值 | 颜色 | 含义 |
|-----------|------|------|
| > 0.90 | 蓝色 | 高置信度 |
| 0.70-0.90 | 绿色 | 中等置信度 |
| 0.50-0.70 | 黄色 | 低置信度 |
| < 0.50 | 红色 | 极低置信度 |

**Surface Compare 配色** (按氨基酸残基电荷，`colorfunc` 逐原子回调):

| 残基类型 | 颜色 | Hex |
|----------|------|-----|
| D, E (酸性) | 红色 | `#ff4444` |
| K, R (碱性) | 蓝色 | `#4466ff` |
| H (弱碱性) | 浅蓝 | `#7799ff` |
| G, A, P (小/中性) | 白色 | `#ffffff` / `#f5f5f5` |
| V, I, L, M, F (疏水) | 暖色 | `#ffddaa` / `#ffcc88` |
| S, T, N, Q (极性) | 浅灰 | `#e8e8e8` |

**注意**: 不能使用 PDB 的 B-factor 做电荷着色。tFold 输出的 B-factor 是 pLDDT（全部 ~0.93-0.98），用 RWB 渐变会映射到同一颜色。必须根据残基类型重新赋色。

#### 10.2.3 色标图例 (Color Scale)

pLDDT 和 Surface 模式在 **viewer 右下角**叠加一个半透明色标面板:

```
┌─────────────────────────────┐
│ pLDDT Confidence            │
│ ████████████████████████    │
│ 0.5 (low)       1.0 (high) │
└─────────────────────────────┘
```

- 使用 `position:absolute; bottom:14px; right:14px` 定位
- 背景 `rgba(0,0,0,0.82)` + `backdrop-filter:blur(8px)`
- 切换到 Side-by-Side / Superpose 时自动隐藏

#### 10.2.4 残基级热力图 (Heatmap Bar)

Viewer 下方展示 **per-residue comparison panel**:

```
         p1    p2    p3    p4    p5    p6    p7    p8    p9
Target:  [G]   [I]   [L]   [G]   [F]   [V]   [F]   [T]   [L]    (绿底)
Decoy:   [L]   [L]   [V]   [G]   [F]   [V]   [F]   [V]   [V]    (匹配绿/错配红底)
pLDDT:   ███   ████  ████  ████  ███   ████  ███   ████  ████   (按置信度着色)
```

- 匹配位: 绿色底 `rgba(126,231,135,0.15)`
- 错配位: 红色底 `rgba(255,123,114,0.15)`
- pLDDT 条: 宽度 = pLDDT%，颜色 = 蓝(>0.9)/青(>0.7)/黄(>0.5)/红(<0.5)
- Hover 显示精确数值

#### 10.2.5 交互功能

| 操作 | 效果 |
|------|------|
| 左键拖拽 | 旋转 |
| 滚轮 | 缩放 |
| 右键拖拽 | 平移 |
| `↑` / `↓` 方向键 | 切换上/下一个 decoy |
| `1` `2` `3` `4` 键 | 切换四种显示模式 |
| 侧边栏 Filter pills | 筛选 Decoy B / Decoy D / All |
| 点击 decoy 卡片 | 加载对应结构，自动滚动到当前卡片 |

#### 10.2.6 信息面板

底部 info-bar 包含:

| 字段 | 来源 |
|------|------|
| Source | Decoy B / Decoy D (MPNN) |
| Hamming Distance | 序列比对 |
| Physicochemical Similarity | Atchley cosine |
| Structural Similarity | Dual-superposition score |
| Risk Score | Composite risk |
| Gene | 基因名 |
| ipTM badge | tFold REMARK 提取，颜色编码 (绿>0.85, 黄>0.7, 红<0.7) |

#### 10.2.7 Chain 约定

Pipeline 输出的 PDB 文件统一使用以下 chain ID:

| Chain | 内容 |
|-------|------|
| `M` | MHC heavy chain (α1+α2+α3) |
| `N` | β2-microglobulin |
| `P` | Peptide |

### 10.3 Decoy A Overview (`decoy_a_overview.html`)

Decoy A 无 3D 结构，使用 **SVG 气泡图 + 序列对齐面板** 的纯 2D 可视化。

#### 10.3.1 气泡图

| 轴 / 属性 | 映射 |
|-----------|------|
| X 轴 | Hamming distance (0-5) |
| Y 轴 | EL%Rank (lower = stronger HLA binder) |
| 气泡大小 | Risk score |
| 气泡颜色 | 关键器官表达: 红 (TPM>10), 黄 (1-10), 蓝 (低/未知) |
| 气泡边框 | 黄 = TCR contact mismatch, 紫 = anchor mismatch |

- Hover 显示 tooltip (序列对齐 + 基因名)
- Click 打开右侧 detail panel
- 重叠点自动 jitter 避免遮挡

#### 10.3.2 Detail Panel (右侧)

- **Stats grid**: HD, EL%Rank, TCR contact mismatches, Anchor mismatches
- **Position-level alignment**: 每个位置一个方块
  - 绿底 = 匹配，红底 = 错配
  - 黄色边框 = TCR 接触位，紫色边框 = anchor 位
- **Expression bar**: 基因名 + TPM 数值 + 颜色条

### 10.4 构建流程

```
输入:
  ├── target PDB (chain M/N/P)
  ├── decoy PDB 列表 (同 chain 规范)
  ├── Decoy B ranked results JSON
  ├── Decoy D results CSV
  └── Decoy A results JSON

构建脚本: build_viz.py
  1. 读取 target PDB + 提取 REMARK (lDDT/pTM/ipTM) + per-residue B-factor
  2. 遍历 Decoy B/D PDB，同样提取 remarks + bfactors
  3. 组装 entry 列表，内嵌 PDB text
  4. 生成 decoy_3d_comparison.html (3Dmol.js)
  5. 读取 Decoy A JSON，生成 decoy_a_overview.html (SVG + CSS)

输出:
  ├── decoy_3d_comparison.html  (~3MB, Decoy B/D 3D 交互)
  └── decoy_a_overview.html     (~300KB, Decoy A 气泡图)
```

### 10.5 使用示例

```bash
# 构建指定靶标的全部可视化
PYTHONUTF8=1 python build_viz.py --target GILGFVFTL

# 指定自定义目录
PYTHONUTF8=1 python build_viz.py --target NLVPMVATV --dir data/NLVPMVATV_summary

# 输出文件（双击即可在浏览器中打开）
data/GILGFVFTL_summary/decoy_3d_comparison.html
data/GILGFVFTL_summary/decoy_a_overview.html
```

### 10.6 性能注意事项

| 限制 | 建议 |
|------|------|
| 3D viewer 的 decoy 数 | ≤ 20 个 (每个 PDB ~150KB, 20 个 ≈ 3MB HTML) |
| Decoy A 气泡图 | ≤ 600 个 hit (SVG 渲染上限，超出后浏览器卡顿) |
| 侧边栏选择制 | 每次只加载 1 个 decoy vs target，不同时渲染多个 |
| 大文件加载 | 3Dmol.js 解析 PDB ~200ms, 渲染 ~100ms |
| 浏览器兼容 | Chrome / Edge / Firefox 均支持 WebGL; Safari 可能较慢 |
| Surface 模式背景 | 必须用 `#2b2b2b` 深灰，不能用纯黑 (白色残基不可见) |

### 10.7 踩坑记录

| 问题 | 原因 | 解决方案 |
|------|------|---------|
| Surface Compare 全部显示紫色 | tFold PDB 的 B-factor 是 pLDDT (0.93-0.98)，RWB 渐变映射到同一色段 | 不用 B-factor，改用 `colorfunc` 按残基类型赋色 |
| pLDDT/Surface 模式有"幽灵蛋白"环绕 | MHC cartoon `#30363d` 与暗色背景撞色 | 这两个模式完全隐藏 MHC (`setStyle({}, {})`) |
| Surface 白色残基在纯黑背景上不可见 | 背景 `#0d1117` 与 `#ffffff` surface 对比太强 | Surface 模式动态切换背景为 `#2b2b2b` |
| `selectedAtoms()` 修改 B-factor 无效 | 3Dmol.js `selectedAtoms` 返回副本，不是引用 | 改用 `colorfunc` 回调直接按 `atom.resn` 返回颜色 |

---

## 12. PLIP 部署与集成

### 12.1 背景

PLIP (Protein-Ligand Interaction Profiler) 是一个自动化的非共价相互作用分析工具，可检测蛋白质-配体之间的氢键、疏水接触、盐桥、π-stacking 等相互作用。

在 Decoy B 中，PLIP 用于分析 peptide-MHC **binding groove 界面**的非共价相互作用模式。这与 TCR-facing 描述符互补：
- **TCR-facing 描述符 (v2)**: 度量 pMHC 暴露给 TCR 的**上表面**特征
- **PLIP**: 度量 peptide 锚定在 MHC groove 中的**结合模式**特征

两者都是判断交叉反应性的重要维度——相似的 groove 结合模式意味着相似的肽段呈递方式，结合 TCR-facing 表面相似性，可以更全面地评估交叉反应风险。

### 12.2 技术挑战：ATOM → HETATM 转换

**问题**: PLIP 设计用于 protein-ligand（蛋白质-小分子）分析，只检测 HETATM 记录（非蛋白分子）与蛋白质的相互作用。在 pMHC 结构中，peptide 是标准蛋白链（ATOM 记录），PLIP 的 `analyze()` 返回 0 个 binding site。

**解决方案**: 新增 `_prepare_pdb_for_plip()` 函数，在分析前将 peptide chain P 的 ATOM 记录转换为 HETATM 记录：

```python
# interface_descriptors.py
def _prepare_pdb_for_plip(pdb_path, peptide_chain="P"):
    """将 peptide chain 的 ATOM → HETATM，使 PLIP 识别为 'ligand'"""
    # 逐行读取 PDB，将 chain P 的 ATOM 替换为 HETATM
    # 写入临时文件，分析后自动清理
```

转换后 PLIP 成功检测到**每个肽段残基**的 binding site，包括：
- H-bond（氢键）: 主要出现在 N 端 (G1) 和 C 端 (T8, L9) 的锚定位
- Hydrophobic contacts（疏水接触）: 集中在中间残基 (I2, L3, F5, F7)
- Salt bridges（盐桥）和 pi-stacking: 较少出现

### 12.3 部署环境

PLIP 依赖 OpenBabel C++ 库（SWIG 绑定），在 Windows 上无法通过 pip 直接安装。采用独立 conda 环境：

```bash
conda create -n plip_env python=3.12 -y
conda install -n plip_env -c conda-forge openbabel=3.1.1 -y
conda run -n plip_env pip install plip --no-deps
conda run -n plip_env pip install numpy lxml
```

**验证结果** (GILGFVFTL / HLA-A*02:01):

```
Binding sites found: 9
  GLY:P:1: H-bonds=3, Hydrophobic=0
  ILE:P:2: H-bonds=1, Hydrophobic=6
  LEU:P:3: H-bonds=1, Hydrophobic=6
  GLY:P:4: H-bonds=0, Hydrophobic=0
  PHE:P:5: H-bonds=0, Hydrophobic=7
  VAL:P:6: H-bonds=0, Hydrophobic=2
  PHE:P:7: H-bonds=0, Hydrophobic=4
  THR:P:8: H-bonds=3, Hydrophobic=1
  LEU:P:9: H-bonds=3, Hydrophobic=5
```

### 12.4 Pipeline 集成

PLIP 集成到 `scanner.py` 的 Stage 5 结构比较中，位于 TCR-facing 描述符计算之后：

```
Stage 5 计算流程:
  1. 双叠合 RMSD (RMSD_Geo)
  2. TCR-facing 描述符 (ESP, Shape, SASA, Hydro) → TCR_Combined
  3. StructSim = 0.50 × RMSD_Geo + 0.50 × TCR_Combined
  4. PLIP Tanimoto similarity (补充信息，记录在 StructuralScore.plip_tanimoto)
```

PLIP Tanimoto 使用广义 Tanimoto 系数比较两个 pMHC 的相互作用指纹：

$$T(A, B) = \frac{A \cdot B}{|A|^2 + |B|^2 - A \cdot B}$$

其中 $A$, $B$ 是 5 维计数向量 `[hbond, hydrophobic, salt_bridge, pi_stack, pi_cation]`。

当前 PLIP Tanimoto 作为**补充描述符**记录在输出中，不参与 StructSim 权重计算。当 PLIP 不可用时，Pipeline 自动跳过该步骤。

### 12.5 测试验证

```bash
# 部署测试
conda run -n plip_env python scripts/test_plip_deploy.py

# 单元测试 (36 passed)
python -m pytest tests/test_interface_descriptors.py -x -q
```

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

### 可视化输出

| 文件 | 内容 |
|------|------|
| `build_viz.py` | 可视化构建脚本 (3D viewer + Decoy A overview) |
| `data/GILGFVFTL_summary/decoy_3d_comparison.html` | Decoy B/D 交互式 3D 比较 (4 模式: Side-by-Side / Superpose / pLDDT / Surface) |
| `data/GILGFVFTL_summary/decoy_a_overview.html` | Decoy A 气泡图 + 序列对齐 detail panel |

### 参考文献

| 文件 | 内容 |
|------|------|
| `docs/deep-research_decoy_b.md` | 界面描述符文献综述（选型依据） |
| `docs/TUTORIAL_AF3_LINUX_INSTALL.md` | AlphaFold 3 安装教程 |

---

## 13. 统一 CLI 入口与批量运行

### 13.1 背景

此前各 Decoy 策略的入口分散：Decoy A 用 `python -m decoy_a`，Decoy B 用 `python -m decoy_b`，Decoy D 需要通过 `python -m decoy_b run --mpnn` 间接触发。缺乏统一的命令行入口和批量运行机制。

### 13.2 新增文件

| 文件 | 用途 |
|------|------|
| `run_decoy.py` | 项目根目录统一 CLI 入口，支持 `python run_decoy.py <SEQUENCE> <strategies...>` |
| `decoy_d/__main__.py` | 使 `python -m decoy_d` 可用 |
| `scripts/batch_decoy_d.sh` | 批量对 `candidate_targets.json` 中所有靶标运行 Decoy D |
| `scripts/visualize_decoy_d_detail.py` | 通用 Decoy D 3-panel 可视化脚本（非 GILGFVFTL 硬编码） |

### 13.3 统一 CLI 设计

`run_decoy.py` 接受位置参数 `<target> <strategies...>`：

```bash
python run_decoy.py GILGFVFTL d              # 单策略
python run_decoy.py GILGFVFTL a b d          # 多策略组合
python run_decoy.py EVDPIGHLY d --hla HLA-A*01:01  # 指定 HLA allele
python run_decoy.py GILGFVFTL d --designs 2000 --top-k 20
```

策略路由：
- `a` → `decoy_b.scanner.scan_decoy_a()`
- `b` → `decoy_b.scanner.scan_decoy_b()`
- `d` → `decoy_d.scanner.run_decoy_d()`
- `c` → 待集成（当前跳过并提示）

输出目录：`data/decoy_{a,b,d}/{TARGET}/`

### 13.4 批量运行脚本

`scripts/batch_decoy_d.sh` 流程：
1. 从 `data/candidate_targets.json` 读取所有靶标（含 existing + proposed，共 12 个）
2. 逐个运行 `python run_decoy.py <SEQ> d --hla <HLA>`
3. 运行完成后调用 `scripts/visualize_decoy_d_detail.py --all` 生成所有 Pareto 图
4. 输出：12 张 `figures/{TARGET}_decoy_d_detail.png`

额外参数可透传：`bash scripts/batch_decoy_d.sh --designs 2000`

### 13.5 可视化脚本

`scripts/visualize_decoy_d_detail.py` 是通用版的 Decoy D 3-panel 图表生成器（从原 `visualize_GILGFVFTL_summary.py` 的 Figure 3 提取并参数化）：

- Panel a: Top 15 by MPNN Score（附 EL%Rank 标注）
- Panel b: Top 15 by EL%Rank（附 MPNN Score 标注）
- Panel c: Pareto front — MPNN Score vs EL%Rank 散点图，红色菱形标出 Pareto 最优解

支持单靶标和批量模式：
```bash
python scripts/visualize_decoy_d_detail.py GILGFVFTL
python scripts/visualize_decoy_d_detail.py --all
```

---

## 14. Decoy C 三重验证 + IEDB 采掘升级

### 14.1 背景与问题

对 Decoy C 管线的系统审计发现了两类关键问题：

**召回问题（Sequence 遗漏）**：
- IEDB 有约 100 万条 T 细胞 epitope 记录，但仅被用作验证工具，从未作为主数据源
- 论文 48k 字符截断导致大型 immunopeptidome 论文的 peptide table 被砍掉
- 每篇论文只做一次 LLM 提取，无重试机制

**精度问题（虚假 Sequence 入库）**：
- LLM 提取的序列从不检查是否真正出现在原论文中（幻觉风险）
- UniProt 验证仅检查 gene→accession 映射，不验证 peptide ⊂ protein
- 验证状态过于宽松：UniProt NOT_FOUND + IEDB NOT_FOUND = PARTIAL 仍入库
- Quality filter 条件太松：`gene_ok OR hla_ok OR summary_ok`

### 14.2 精度修复：三重验证 + 硬拒绝

#### 14.2.1 源文本验证 (`extractor.py: verify_peptide_in_source`)

提取后自动检查每个 peptide 序列是否出现在原论文全文中：
- 直接字符串匹配 + 去除空格/连字符后匹配
- 设置 `validation_flags["source_text_match"] = FOUND | NOT_FOUND`
- NOT_FOUND 意味着 LLM 可能幻觉出了不存在的序列

#### 14.2.2 蛋白子序列验证 (`validator.py: verify_peptide_in_protein`)

从 UniProt REST API 拉取全长蛋白 FASTA 序列，验证 `peptide in protein_seq`：
- 使用 `https://rest.uniprot.org/uniprotkb/{ID}.fasta` 获取全长蛋白
- 设置 `validation_flags["protein_containment"] = CONFIRMED | NOT_FOUND | SKIPPED | ERROR`
- NOT_FOUND 是强幻觉信号：声称来自某蛋白但序列不在蛋白中

#### 14.2.3 硬拒绝门控 (`validator.py: hard_reject`)

Entry 满足以下任一条件即被拒绝入库：

1. **protein_containment = NOT_FOUND** → 直接拒绝（序列不在蛋白中 = 幻觉）
2. **三重失败**：source_text_match NOT_FOUND + iedb_match NOT_FOUND + protein_containment NOT_FOUND/SKIPPED → 拒绝（零外部证据）
3. **overall_status = REJECTED**（3+ 项验证失败）

例外：seed/manual/iedb_miner 来源的条目不受硬拒绝影响。

#### 14.2.4 双次共识提取 (`extractor.py: _extract_consensus`)

同一论文运行两次 LLM 提取（temperature=0.1 和 0.4），仅保留两次都提取到的序列：
- 大幅降低幻觉率（两次独立调用同时幻觉出相同序列的概率极低）
- 代价：2× API 调用成本
- 通过 `--consensus` CLI flag 启用

### 14.3 召回修复：IEDB 系统性采掘

新增 `iedb_miner.py` 模块，将 IEDB 从验证工具升级为主数据源：

| 策略 | 方法 | 预期产出 |
|------|------|---------|
| Protein | 对 60+ 危险蛋白（TTN, MAGEA3, cardiac, neural）查询所有 T-cell assay | 高置信度脱靶肽 |
| Self-antigen | 人类自身抗原 × 11 个常见 HLA allele 交叉查询 | 潜在脱靶池 |
| Mass-spec | 仅保留有 MHC-elution 质谱确认的天然呈递肽 | 最高置信度 |
| Cross-reactive | 查询自身免疫/交叉反应/分子拟态标注的记录 | 已知交叉反应 |

### 14.4 质量过滤升级 (`scale_up.py: _quality_filter`)

旧标准（任一即可）：`gene_ok OR hla_ok OR summary_ok`

新标准（全部满足 + 证据要求）：
- **必须全部满足**：非 UNKNOWN gene + 有效 HLA allele + evidence_summary > 20 chars
- **至少一项证据**：trusted_method / VALIDATED status / source_text_match / protein_containment / iedb_match

### 14.5 文件变更清单

| 文件 | 变更 |
|------|------|
| `extractor.py` | +`verify_peptide_in_source()`, +`_extract_consensus()`, +`_call_llm_t2()`, `extract_from_paper()` 增加 consensus 参数 |
| `validator.py` | +`_fetch_protein_sequence()`, +`verify_peptide_in_protein()`, +`hard_reject()`, `validate_entry()` 增加蛋白验证步骤 |
| `orchestrator.py` | `process_papers()` 集成硬拒绝门控和 consensus 参数 |
| `iedb_miner.py` | **新文件** — IEDB 系统性采掘模块（4 策略） |
| `scale_up.py` | 集成 IEDB miner, 传递 consensus flag, 增强质量过滤 |

### 14.6 CLI 用法

```bash
# 标准扩库（含 IEDB 采掘 + 三重验证）
python scale_up.py -n 10000

# 高精度模式（双次共识 + 三重验证）
python scale_up.py -n 10000 --consensus

# 仅 IEDB 采掘（跳过 PubMed）
python scale_up.py -n 10000 --no-pubmed
```

### 14.7 IEDB 大规模采掘结果 (2026-04-07)

#### 挖掘过程
1. **策略 A (人类自身抗原 T 细胞检测)**：11 个常见 HLA allele × Homo sapiens 阳性记录 → **757 个候选条目**
2. **策略 B (自身免疫疾病)**：0 个新增（被策略 A 覆盖）
3. **三重验证**：757 候选 → 745 存活初轮验证
4. **UniProt ID 修复**：IEDB 记录 `parent_source_antigen_iri` 字段含 UniProt accession（格式 `UNIPROT:P60709`），原代码未提取。修复后 675 条获得正确 uniprot_id
5. **蛋白子序列重验证**：675 条用正确 UniProt ID 重跑 `verify_peptide_in_protein()`
6. **清退**：33 条蛋白子序列验证失败（peptide 不在声称的蛋白中）被移除

#### 关键发现：gene_symbol 错误
IEDB 记录中 `gene_symbol` 字段通常为空，代码回退使用 `source_antigen_name`（完整蛋白名如 "Glypican-3"），导致 UniProt 基因查询全部失败。

**修复**：
- `iedb_miner.py`：从 `parent_source_antigen_iri` 提取 UniProt accession
- `iedb_miner.py`：新增 `_protein_to_gene()` 映射（~60 个常见蛋白名→基因符号）
- `fix_and_revalidate.py`：从 `iedb_source_antigen` validation flag 中提取 UniProt ID 并重跑蛋白验证

#### 最终库统计

| 指标 | 数值 |
|------|------|
| 总条目数 | 817 |
| 唯一序列 | 817 |
| VALIDATED | 662 (81.0%) |
| PARTIAL | 111 (13.6%) |
| NEEDS_REVIEW | 44 (5.4%) |
| 蛋白子序列确认 | 667/817 (81.6%) |
| IEDB 匹配 | 768/817 (94.0%) |
| 双重验证（金标准） | 662/817 (81.0%) |
| 质谱确认 | 476/817 (58.3%) |
| 唯一 HLA allele | 46 |
| 唯一基因 | 310 |

| HLA Allele | 数量 |
|------------|------|
| HLA-A*02:01 | 223 |
| HLA-B*07:02 | 133 |
| HLA-A*24:02 | 131 |
| HLA-A*03:01 | 101 |
| HLA-A*11:01 | 61 |
| HLA-A*01:01 | 41 |

#### 数据来源
- `iedb_miner` (IEDB 直接采掘): 712 条
- `llm` (PubMed 文献 + GPT-4o 提取): 60 条
- `llm_batch_migrated` (旧批次迁移): 45 条
