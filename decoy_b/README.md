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

rmsd_geo = avg(sim_A, sim_B)
  sim_A = 1 - peptide_RMSD / 3.0    (MHC 叠合后)
  sim_B = 1 - groove_RMSD / 3.0     (肽段叠合后)

interface_combined = adaptive_weighted_avg(
    0.30 * tcr_esp_similarity,       # TCR-facing 静电势 (Hodgkin SI)
    0.30 * tcr_shape_similarity,     # 暴露侧链形状 (Kabsch RMSD)
    0.20 * tcr_sasa_similarity,      # rSASA 暴露模式 (cosine)
    0.20 * tcr_hydro_similarity,     # 暴露疏水性 (cosine)
)

similarity_score = surface_correlation   # 纯结构评分
risk = similarity_score * (1/EL_Rank) * TPM_Weight
```

## 4 个 TCR-facing 表面描述符

> **关键设计变更**: 所有描述符现在度量 pMHC **暴露给 TCR 的上表面**（溶剂可及面），
> 而非 peptide-MHC 结合 groove 界面。通过 per-atom SASA > 0 过滤，自动排除
> groove 内埋藏的锚定残基，只保留 TCR 实际接触的原子。

| 描述符 | 方法 | 度量 | 权重 |
|--------|------|------|------|
| **TCR-facing ESP** | 暴露肽原子处的 Coulomb 电位 (来自周围所有原子) | Hodgkin SI | 0.30 |
| **Exposed Shape** | MHC 叠合后暴露侧链质心的 3D 坐标 (Kabsch RMSD) | 1-RMSD/3.0 | 0.30 |
| **rSASA Profile** | 每个肽位点的相对溶剂可及面积 (complex/free) | Cosine 相似度 | 0.20 |
| **Exposed Hydrophobicity** | rSASA 加权的 Kyte-Doolittle 疏水性 | Cosine 相似度 | 0.20 |

**已替换的旧描述符** (度量方向错误 — peptide-MHC groove 界面):
- ~~PLIP~~, ~~FreeSASA/BSA~~, ~~PRODIGY~~, ~~PeSTo~~, ~~APBS/ESP~~

**Spearman ρ 验证** (GILGFVFTL 50 candidates):
- tcr_esp ↔ sim_A: 0.43 (独立信息)
- tcr_shape ↔ sim_A: 0.22 (高度独立)
- tcr_sasa ↔ sim_A: 0.72 (中等相关)
- tcr_hydro ↔ sim_A: 0.15 (高度独立)

所有描述符支持 **adaptive weighting**：如果某个描述符计算失败，权重自动重分配给成功的。

## 部署依赖

| 工具 | 用途 | 必需? |
|------|------|-------|
| **tFold** | pMHC 结构预测 | 推荐 (Stage 2 核心) |
| **AlphaFold3** | 高精度精修 | 可选 |
| **Boltz-2** | 结构交叉验证 | 可选 |
| **ProteinMPNN** | 逆向序列设计 | 可选 |
| **BioPython ShrakeRupley** | SASA 计算 (TCR-facing 描述符核心) | 必需 |

Pipeline 仅依赖 BioPython 即可运行全部 4 个 TCR-facing 描述符，无外部工具依赖。

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
| `tools/tcr_surface_descriptors.py` | 4-descriptor TCR-facing 表面描述符模块 |
| `tools/interface_descriptors.py` | (旧) peptide-MHC groove 描述符 (已弃用) |
| `tools/tfold.py` | tFold pMHC 预测封装 |
| `tools/alphafold3.py` | AF3 高精度精修封装 |
| `tools/proteinmpnn.py` | ProteinMPNN 逆向设计封装 |
| `DEPLOYMENT.md` | 外部工具部署指南 |

## 参考文献与工具来源

### Stage 1: 序列级描述符

| 工具/方法 | 论文 | 链接 |
|-----------|------|------|
| **Atchley Factors** | Atchley WR et al. "Solving the protein sequence metric problem." *PNAS* 102(18):6395–6400, 2005 | [DOI](https://doi.org/10.1073/pnas.0408677102) |

### Stage 2–4: 结构预测工具

| 工具 | 论文 | 代码仓库 |
|------|------|---------|
| **tFold** | Wu L et al. "tFold-Ab/Ag: Fast and Accurate Antibody and Antigen Structure Prediction." 2024 | [GitHub](https://github.com/TencentAI4S/tfold) |
| **AlphaFold3** | Abramson J et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493–500, 2024 | [GitHub](https://github.com/google-deepmind/alphafold3) |
| **Boltz-2** | Wohlwend J et al. "Boltz-1: Democratizing Biomolecular Interaction Modeling." 2024; Boltz-2 (2025) | [GitHub](https://github.com/jwohlwend/boltz) · [HuggingFace](https://huggingface.co/boltz-community/boltz-2) |
| **ProteinMPNN** | Dauparas J et al. "Robust deep learning–based protein sequence design using ProteinMPNN." *Science* 378(6615):49–56, 2022 | [GitHub](https://github.com/dauparas/ProteinMPNN) |

### Stage 5: 界面描述符 (5-descriptor)

| 描述符 | 论文 | 代码/工具 |
|--------|------|----------|
| **PLIP** | Adasme MF et al. "PLIP 2021: expanding the scope of the protein–ligand interaction profiler to DNA and RNA." *Nucleic Acids Res.* 49(W1):W530–W534, 2021 | [GitHub](https://github.com/pharmai/plip) · `pip install plip` |
| **FreeSASA** | Mitternacht S. "FreeSASA: An open source C library for solvent accessible surface area calculations." *F1000Research* 5:189, 2016 | [GitHub](https://github.com/mitternacht/freesasa) · `pip install freesasa` |
| **PRODIGY** | Xue LC et al. "PRODIGY: a web server for predicting the binding affinity of protein–protein complexes." *Bioinformatics* 32(23):3676–3678, 2016 | [GitHub](https://github.com/haddocking/prodigy) · `pip install prodigy-prot` |
| **APBS** | Jurrus E et al. "Improvements to the APBS biomolecular solvation software suite." *Protein Sci.* 27(1):112–128, 2018 | [GitHub](https://github.com/Electrostatics/apbs) · `conda install -c conda-forge apbs` |
| **pdb2pqr** | Dolinsky TJ et al. "PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations." *Nucleic Acids Res.* 35(suppl 2):W522–W525, 2007 | [GitHub](https://github.com/Electrostatics/pdb2pqr) · `conda install -c conda-forge pdb2pqr` |
| **PeSTo** | Krapp LF et al. "PeSTo: parameter-free geometric deep learning for accurate prediction of protein binding interfaces." *Nature Commun.* 14:2175, 2023 | [GitHub](https://github.com/LBM-EPFL/PeSTo) · `pip install pesto-protein` |

### Stage 5: 结构叠合 RMSD

| 工具 | 论文 | 代码 |
|------|------|------|
| **BioPython Superimposer** | Cock PJA et al. "Biopython: freely available Python tools for computational molecular biology and bioinformatics." *Bioinformatics* 25(11):1422–1423, 2009 | [GitHub](https://github.com/biopython/biopython) · `pip install biopython` |

### 其他关键参考

| 主题 | 论文 |
|------|------|
| **Titin 致死案例 (Decoy B 核心动机)** | Cameron BJ et al. "Identification of a Titin-derived HLA-A1–presented peptide as a cross-reactive target for engineered MAGE A3–directed T cells." *Sci. Transl. Med.* 5(197):197ra103, 2013. [DOI](https://doi.org/10.1126/scitranslmed.3006034) |
| **TCR 交叉反应性综述** | Sewell AK. "Why must T cells be cross-reactive?" *Nature Rev. Immunol.* 12:669–677, 2012. [DOI](https://doi.org/10.1038/nri3279) |
| **MatchTope (备选，未采用)** | Mendes MFA et al. "MatchTope: A tool to predict the cross reactivity of peptides complexed with MHC I." *Front. Immunol.* 13:930590, 2022. [DOI](https://doi.org/10.3389/fimmu.2022.930590) · [GitHub](https://github.com/Marcus-Mendes/MatchTope) — 与 APBS/ESP 描述符功能重叠，未集成 |
