# Deep Research: TCR-Facing Surface Descriptors for Decoy B

## 问题

当前 Decoy B 管线的界面描述符（PLIP、BSA、PRODIGY、ESP、PeSTo）全部在度量 **peptide-MHC 结合界面**（peptide 埋在 groove 里的面），而 TCR 交叉反应取决于 **pMHC 暴露给 TCR 的上表面**。这是两个不同的界面：

```
        TCR
    ─────┬─────
         │  ← TCR-pMHC 界面 (TCR 真正看到的面)  ← 应该比较的
    ═══peptide═══
         │  ← peptide-MHC 界面 (peptide 埋在groove里)  ← 之前在比较的
    ───MHC groove───
```

需要重新设计描述符，使其度量 TCR-facing surface 的相似度。

---

## 文献综述

### 1. 专用 pMHC 表面交叉反应性工具

#### 1.1 MatchTope (PIPSA-based)
- **论文**: Mendes et al. (2022) Front Immunol 13:930590
- **GitHub**: https://github.com/Marcus-Mendes/MatchTope
- **方法**: 在 pMHC 上表面放置圆柱体 (R=40Å, L=33Å)，仅在此区域内计算 MEP (UHBD → PIPSA → Hodgkin SI)
- **TCR-facing 定义**: 圆柱体自动排除 groove 内部，只保留上表面
- **输出**: Hodgkin Similarity Index, 范围 [-1, 1]
- **局限**: 仅支持 HLA-A*02:01 等少数等位基因；只有静电势

#### 1.2 PepSim
- **论文**: Hall-Swan et al. (2023) Front Immunol 14:1108303
- **Web**: https://pepsim.kavrakilab.org (Kavraki Lab, Rice)
- **方法**: MSMS 三角化分子表面 → 以 peptide 质心为中心提取 patch (16 edges) → 三通道注释 (ESP + 疏水性 + 氢键势)
- **TCR-facing 定义**: 以 peptide 质心为中心的圆形 patch，天然捕获 peptide + 侧翼 MHC 螺旋的上表面
- **输出**: 综合 PepSim score (结构 + 序列)

#### 1.3 CrossDome
- **论文**: Borden et al. (2023) Front Immunol
- **GitHub**: https://github.com/AntunesLab/crossdome
- **方法**: 12 维生化性质空间，按 TCR 接触位点加权 (P4,P5,P6,P8 在 HLA-A*02:01)
- **局限**: 纯序列工具，非结构

### 2. 几何深度学习表面方法

#### 2.1 MaSIF / MaSIF-search
- **论文**: Gainza et al. (2020) Nature Methods 17:184-192
- **GitHub**: https://github.com/LPDI-EPFL/masif
- **方法**: MSMS 分子表面 → 5 通道特征 (shape index, curvature, ESP, hydropathy, H-bond potential) → geodesic CNN → **80 维 fingerprint per patch**
- **MaSIF-search**: Siamese 网络训练，使互作 patch 的 descriptor 距离小
- **适用性**: 高 — 可提取 TCR-facing patch 并比较 fingerprint
- **依赖**: MSMS, APBS, PyMesh, TensorFlow 1.x → Docker 推荐

#### 2.2 dMaSIF
- **论文**: Sverrisson et al. (2021) CVPR/NeurIPS
- **GitHub**: https://github.com/FreyrS/dMaSIF
- **方法**: 无需 MSMS/APBS 预计算，直接从原子坐标生成点云表面，端到端学习
- **速度**: 比 MaSIF 快 600x
- **依赖**: PyTorch, PyTorch Geometric
- **适用性**: 最实用的 MaSIF 替代，PyTorch 原生

#### 2.3 IMPRINT (最直接相关!)
- **论文**: Shang, Chan & Zhou (2026) Briefings in Bioinformatics 27(2):bbag048
- **GitHub**: https://github.com/Xiyougailv/IMPRINT
- **方法**: 基于 MaSIF 架构，专门在 pMHC 的 **TCR-facing 表面** 上计算 "immunological fingerprint"
  - 界面定义: peptide 原子 4Å 内的表面点
  - Patch 大小: 12Å geodesic radius, 32 patches/structure
  - 特征: 5 通道 MaSIF 特征
  - 输出: 80 维 descriptor → 80×80 fingerprint matrix
- **验证**: 在 40 个 HLA-A*02-peptide-TCR 晶体结构上训练，成功预测 TCR 结合偏好
- **适用性**: 极高 — 这正是我们需要的

#### 2.4 PeSTo
- **论文**: Krapp et al. (2023) Nature Communications 14:2175
- **GitHub**: https://github.com/LBM-EPFL/PeSTo
- **输出**: per-atom interface probability，不是固定长度 fingerprint
- **适用性**: 中等 — 需要自定义聚合方案

#### 2.5 PEP-Patch
- **论文**: J Chem Inf Model 2023
- **GitHub**: https://github.com/liedllab/surface_analyses (MIT)
- **方法**: APBS → 表面静电势 patch 分析（正/负 patch 面积和强度）
- **适用性**: 轻量级 ESP patch 替代方案

### 3. 经典方法

#### 3.1 PIPSA / webPIPSA
- **论文**: Wade et al. (2008) Nucleic Acids Research
- **Web**: https://pipsa.h-its.org/pipsa/
- **方法**: 分子表面 "skin" 上的 ESP 两两比较 → Hodgkin SI / Carbo SI
- **ROI 支持**: 可定义特定区域 (如 MatchTope 的圆柱体)

#### 3.2 3D Zernike Descriptors
- **论文**: Kihara et al. (2011) Curr Protein Pept Sci; Sael & Kihara (2009)
- **方法**: 旋转不变的表面形状 + 静电势 → 121 维 compact vector
- **优势**: 超快比较（向量距离），无需预对齐
- **适用性**: 适合大规模筛选

---

## TCR-Facing Surface 的定义

### 9-mer in HLA-A*02:01 的位点分配

| 位置 | 角色 | rSASA 范围 | 说明 |
|------|------|-----------|------|
| P1 | 部分埋藏 | 0.15-0.30 | A pocket |
| **P2** | **主锚定** | 0.00-0.05 | B pocket (Leu/Met) |
| P3 | 部分暴露 | 0.10-0.25 | |
| **P4** | **TCR 接触** | 0.40-0.70 | CDR3 接触 |
| **P5** | **最暴露** | 0.50-0.80 | CDR3 核心接触 |
| **P6** | **TCR 接触** | 0.40-0.70 | CDR3 接触 |
| P7 | 可变 | 0.20-0.50 | |
| **P8** | **TCR 接触** | 0.30-0.60 | CDR3 接触 |
| **P9** | **主锚定** | 0.00-0.10 | F pocket (Val/Leu) |

**注意**: 位点分配随 HLA 等位基因和肽段长度变化！对于 10-11mer，中央区域 (P4-P6) 会形成 bulge 更加暴露。

### 计算提取方法

1. **SASA 过滤法**: 计算 pMHC 复合物中 peptide 每个原子的 SASA，SASA > 0 的原子为溶剂暴露（TCR 可及）
2. **MatchTope 圆柱法**: 在 groove 上方放置 R=40Å 圆柱体，只取柱内的分子表面
3. **PepSim Patch 法**: MSMS 三角化表面 → 以 peptide 质心为中心提取 geodesic patch
4. **法向量法**: 选取表面法向量 z 分量 > 0 的点（指向上方/溶剂）

### TCR footprint 参数

- 典型 BSA: 1600-2000 Å²
- Peptide 贡献: 25-40% BSA
- Shape complementarity (Sc): 0.55-0.70
- MHC 接触残基: α1 helix R65,K66,A69,Q72,T73,R75; α2 helix A150,H151,A152,E154,Q155,W167

---

## Titin/MAGE-A3 案例的表面分析

Cameron et al. (2013) Sci Transl Med 5:197ra103:

- **EVDPIGHLY** (MAGE-A3) vs **ESDPIVAQY** (Titin), Hamming=5
- 共享: P1(E), P3-P4(DP), P5(I), P9(Y) — 锚定 + TCR 核心接触位保守
- 差异: P2(V/S), P6(G/V), P7(H/A), P8(L/Q) — 部分是 groove 埋藏位
- **关键**: Cα RMSD = 0.285Å，backbone 构象几乎完全相同
- P4 的 Pro 产生特征性 kink，使 P5 在两个复合物中位置完全一致
- **TCR-facing 表面的形状和电荷分布高度相似，尽管序列差异大**

**教训**: Hamming 距离、BLOSUM 等序列指标完全无法预测此交叉反应。只有 TCR-facing surface 的结构比较才能捕获。

---

## 描述符重设计方案

### 方案比较

| 方案 | 描述符 | 优势 | 劣势 | 复杂度 |
|------|--------|------|------|--------|
| **A: SASA-filtered ESP + PeSTo** | (1) rSASA per-position (2) TCR-facing ESP (Hodgkin SI) (3) PeSTo on exposed atoms | 实现简单，仅修改现有代码 | 精度有限 | 低 |
| **B: dMaSIF surface fingerprint** | dMaSIF 在 TCR-facing patch 上生成学习的 surface embedding | 端到端学习，PyTorch 原生，无需 MSMS | 需要适配 patch 提取 | 中 |
| **C: IMPRINT** | MaSIF 架构的 immunological fingerprint | 最直接，文献验证 | 依赖 MSMS/PyMesh/TF1.x，Docker | 高 |
| **D: 混合** | rSASA + TCR-facing ESP + dMaSIF/PeSTo | 覆盖面广 | 维护多套代码 | 中-高 |

### 推荐: 方案 D (分层实施)

**Layer 1 — 快速修复 (修改现有代码):**
- 将 ESP 采样点从 "interface atoms" 改为 "SASA > 0 的 peptide 原子"（TCR-facing atoms）
- 使用 Hodgkin SI 替代 cosine similarity
- PeSTo proxy 改为只考虑暴露残基的 contact density

**Layer 2 — 新增 TCR-facing 描述符:**
- per-position rSASA profile (9 维向量)
- TCR-facing shape descriptor (3D Zernike 或 dMaSIF embedding)

**Layer 3 — 高精度 (可选):**
- dMaSIF 或 IMPRINT 的 80 维 surface fingerprint

### 新评分公式 (提案)

```
surface_correlation = 0.50 * rmsd_geo + 0.50 * interface_combined

rmsd_geo = avg(sim_A, sim_B)    # 不变

interface_combined = adaptive_weighted_avg(
    0.30 * tcr_facing_esp,      # APBS/Hodgkin SI on SASA>0 atoms (替代旧 ESP)
    0.25 * tcr_facing_shape,    # dMaSIF embedding 或 rSASA profile cosine
    0.25 * plip_tanimoto,       # PLIP (保留，但需确认是否有 TCR-relevant interactions)
    0.20 * bsa_exposed,         # TCR-exposed SASA (不是 buried SA)
)
```

---

## 关键参考文献

| 论文 | 主题 | DOI/Link |
|------|------|----------|
| Cameron et al. 2013, Sci Transl Med | Titin/MAGE-A3 致死案例 | 10.1126/scitranslmed.3006034 |
| Mendes et al. 2022, Front Immunol | MatchTope (PIPSA-based MEP) | 10.3389/fimmu.2022.930590 |
| Hall-Swan et al. 2023, Front Immunol | PepSim (表面 + 序列) | 10.3389/fimmu.2023.1108303 |
| Gainza et al. 2020, Nature Methods | MaSIF (surface fingerprint) | 10.1038/s41592-019-0666-6 |
| Shang et al. 2026, Brief Bioinform | IMPRINT (immunological fingerprint) | 10.1093/bib/bbag048 |
| Sverrisson et al. 2021, NeurIPS | dMaSIF (fast surface learning) | N/A (CVPR 2021) |
| Krapp et al. 2023, Nat Commun | PeSTo (interface prediction) | 10.1038/s41467-023-37701-8 |
| Rossjohn et al. 2015, Annu Rev Immunol | TCR-pMHC 结构生物学综述 | 10.1146/annurev-immunol-032414-112334 |
| Borden et al. 2023, Front Immunol | CrossDome (TCR-centered weighting) | 10.3389/fimmu.2023.1142835 |
| Wade et al. 2008, NAR | PIPSA/webPIPSA | 10.1093/nar/gkn181 |
