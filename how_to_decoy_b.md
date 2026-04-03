  两个 PDB 文件：target pMHC 和 candidate pMHC（都由 tFold 或 AF3 预测生成）

  第一步：分链提取 CA 原子

  tFold 输出的 PDB 按链组织：
  - Chain M = MHC heavy chain（重链）
  - Chain N = β2-microglobulin（β2m）
  - Chain P = Peptide（肽段）

  代码提取每条链的 Cα 原子（蛋白质骨架的代表原子）。

  第二步：基于 MHC 的结构叠合（Superimposition）

  scanner.py:333-339
  用 BioPython 的 Superimposer 对 MHC 链（M+N）的 Cα 原子 进行最小二乘叠合：
  1. 以 target 的 MHC 为参考，旋转/平移 candidate 的整个结构
  2. 这样做的目的是：对齐 HLA 凹槽，然后在同一坐标系下比较肽段构象

  关键点：不是直接叠合肽段，而是先叠合 MHC，再测量肽段的偏差。这模拟的是 TCR 从上方"看"MHC-肽段复合物的视角。

  第三步：Peptide RMSD 计算

  scanner.py:342-369
  MHC 叠合之后，不再做额外拟合，直接测量两个肽段 Cα 原子坐标之间的偏差：

  同长度肽段（如两个 9-mer）：
  diff = target_pep_coords - cand_pep_coords
  peptide_rmsd = sqrt(mean(sum(diff², axis=1)))
  逐原子坐标差 → 平方和 → 均值 → 开根号 = 标准 RMSD

  不同长度肽段（如 9-mer vs 10-mer）：
  取两者的 TCR-contact core（中央 5 个残基的 Cα），再算 RMSD：
  n_contact = min(5, n - 2)
  start = (n - n_contact) // 2
  core_cas = cas[start : start + n_contact]
  这样 9-mer 取 p3-p7，10-mer 取 p3-p7 或 p4-p8，保证比较的是 TCR 接触面。

  第四步：Surface Correlation（B-factor 相关性）

  scanner.py:372-380
  用 B-factor 相关系数作为表面柔性相似度的代理指标：
  corr = np.corrcoef(target_pep_bfactors, cand_pep_bfactors)[0, 1]
  B-factor 反映原子的热运动幅度 → 间接反映表面暴露程度和动态性。如果两个肽段在相同位置上有相似的 B-factor 分布，说明它们的 表面柔性特征 相似。

  第五步：RMSD → Similarity 转换 + 取最大值

  scanner.py:385-393
  如果 peptide RMSD < 3.0 Å，把它转换为 0-1 的相似度分数：
  rmsd_similarity = max(0.0, 1.0 - rmsd / 3.0)
  # RMSD=0 → similarity=1.0
  # RMSD=1.5 → similarity=0.5
  # RMSD=3.0 → similarity=0.0
  然后取 max(B-factor correlation, RMSD similarity) 作为最终的 surface_correlation。

  第六步：Combined Score（最终相似度）

  scanner.py:839-842（修复后）
  在构建 DecoyBHit 时，如果有结构数据：
  combined = 0.4 * cosine_similarity + 0.6 * surface_correlation
  - 40% 权重：Atchley factor 理化相似度（序列层面，TCR 接触区氨基酸的理化性质）
  - 60% 权重：3D 结构相似度（surface_correlation）

  ---
  总结

  tFold/AF3 预测 pMHC 3D 结构
           ↓
      提取 MHC 链 Cα → Superimpose (叠合 HLA 凹槽)
           ↓
      在叠合后坐标系下测量 Peptide Cα RMSD
           ↓
      RMSD → 0-1 similarity (线性映射, 3Å cutoff)
           ↓
      B-factor 相关系数 (表面柔性)
           ↓
      surface_correlation = max(rmsd_similarity, bfactor_corr)
           ↓
      combined = 0.4 × Atchley cosine + 0.6 × surface_correlation
