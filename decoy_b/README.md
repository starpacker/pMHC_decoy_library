# Decoy B — 结构相似筛选

**状态**: 代码框架完成，5-descriptor 界面评分已实现

## 核心思想

放宽序列相似性要求，寻找 **物化性质相似** 和 **三维结构相似** 的脱靶候选。经典案例：MAGE-A3 靶标 EVDPIGHLY 与 Titin 的 ESDPIVAQY (Hamming=5) 在 HLA-A\*01:01 上 TCR 接触面高度相似，导致临床致死。

## 五阶段管线

```
Stage 1: 物化筛选 (CPU, 秒级)
  Atchley factor → 25维向量 → cosine_sim ≥ 0.70
  排除 Hamming ≤ 2 (已由 Decoy A 覆盖)
  → 2,000-5,000 candidates

Stage 2: tFold 批量结构预测 (GPU, 分钟级)
  pMHC 3D 结构建模 → PDB + 置信度

Stage 3: AF3 精细化 (GPU, 可选)
  Top-200 → AlphaFold3 高精度重建

Stage 4: Boltz 交叉验证 (GPU, 可选)
  独立结构预测 → tFold/AF3 一致性评估

Stage 5: 结构比较 & 综合评分
  双叠合 RMSD + 5-descriptor 界面评分
```

## 评分公式

```
surface_correlation = 0.50 * rmsd_geo + 0.50 * interface_combined

rmsd_geo = max(avg(sim_A, sim_B), bf_corr)
  sim_A = 1 - peptide_RMSD / 3.0    (MHC 叠合后)
  sim_B = 1 - groove_RMSD / 3.0     (肽段叠合后)

interface_combined = adaptive_weighted_avg(
    0.25 * plip_tanimoto,       # PLIP 非共价相互作用
    0.10 * bsa_similarity,      # FreeSASA 埋藏面积
    0.20 * prodigy_similarity,  # PRODIGY 结合亲和力
    0.25 * esp_similarity,      # APBS 静电势
    0.20 * pesto_similarity,    # PeSTo 界面嵌入
)

similarity_score = 0.4 * cosine_sim + 0.6 * surface_correlation
risk = similarity_score * (1/EL_Rank) * TPM_Weight
```

## 5 个界面描述符

| 描述符 | 方法 | 度量 | 权重 |
|--------|------|------|------|
| **PLIP** | 非共价相互作用指纹 (H-bond, 疏水, 盐桥, pi-stacking) | Tanimoto 系数 | 0.25 |
| **FreeSASA** | 埋藏表面积 BSA = SASA(pep) + SASA(MHC) - SASA(complex) | 差值归一化 | 0.10 |
| **PRODIGY** | 基于接触的结合亲和力 ΔG 预测 | ΔG 差值归一化 | 0.20 |
| **APBS** | 静电表面电位 (Poisson-Boltzmann / Coulomb proxy) | Pearson 相关 | 0.25 |
| **PeSTo** | 几何深度学习界面嵌入 / contact-propensity proxy | Cosine 相似度 | 0.20 |

所有描述符支持 **adaptive weighting**：如果某个描述符计算失败，权重自动重分配给成功的。

## 部署依赖

| 工具 | 用途 | 必需? |
|------|------|-------|
| **tFold** | pMHC 结构预测 | 推荐 (Stage 2 核心) |
| **AlphaFold3** | 高精度精修 | 可选 |
| **Boltz-2** | 结构交叉验证 | 可选 |
| **ProteinMPNN** | 逆向序列设计 | 可选 |
| **PLIP** | 相互作用指纹 | 推荐 (界面描述符) |
| **FreeSASA** | 埋藏面积 | 推荐 |
| **PRODIGY** | 结合亲和力 | 推荐 |
| **pdb2pqr + APBS** | 静电势 | 推荐 (有 BioPython Coulomb fallback) |
| **PeSTo + PyTorch** | 界面嵌入 | 推荐 (有 BioPython contact proxy fallback) |

Pipeline 会自动检测可用工具并降级。即使只有 BioPython，也能运行全部 5 个描述符的 proxy 版本。

## 使用方法

```bash
# 仅物化筛选 (纯 CPU, 无外部依赖)
python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01 --skip-structural

# 完整管线
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# 含 MPNN 逆向设计
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --mpnn
```

## 文件清单

| 文件 | 作用 |
|------|------|
| `scanner.py` | 五阶段管线 + 双叠合 RMSD + 界面描述符集成 |
| `risk_scorer.py` | A+B 统一风险评分 |
| `orchestrator.py` | 全管线编排 + 断点续传 |
| `tools/interface_descriptors.py` | 5-descriptor 界面描述符模块 |
| `tools/tfold.py` | tFold pMHC 预测封装 |
| `tools/alphafold3.py` | AF3 高精度精修封装 |
| `tools/proteinmpnn.py` | ProteinMPNN 逆向设计封装 |
| `DEPLOYMENT.md` | 外部工具部署指南 |
