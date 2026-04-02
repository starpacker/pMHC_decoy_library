# Decoy Library 部署与优化工作进展报告
**汇报日期**: 2026年4月2日
**项目**: TCR 脱靶毒性 pMHC 负样本库 (Decoy Library)

本报告详细总结了我们在部署、重构及优化 Decoy Library（特别是结构建模主导的 Decoy B 管线）方面的最新进展。当前代码库已经达到高度自包含（Self-Contained）状态，核心逻辑已打通，并成功运行了基于 128 万候选多肽的全链路测试。

---

## 1. 架构重构与 Self-Contained 改造

为了降低后期的维护成本并解决依赖冲突，我们对 Decoy B 的外部依赖架构进行了深度重构：

*   **代码内置化 (Internalization)**: 
    *   将原本作为外部依赖的 `tfold` (Tencent AI4S) 源码完整迁移至 `decoy_b/external/tfold`。
    *   将逆向设计模块 `ProteinMPNN` 的源码迁移至 `decoy_b/external/proteinmpnn`。
    *   通过动态修改 `sys.path` 的方式，使得工程可以直接在内部执行 `import tfold`，**彻底移除了对外部 Python 环境/包的侵入性修改要求**。
*   **配置解耦**: 
    *   外部环境中唯一保留的是极其庞大的**模型权重文件**与个别不可拆分的工具（如 AlphaFold3）。
    *   所有的权重路径与工具路径被统一收敛到 `decoy_a/config.py` 中进行管理（例如 `TFOLD_DIR`, `AF3_DIR`）。使用者只需修改一个配置即可完成整个管线的迁移。

---

## 2. 核心算法修复与多长度多肽 (Multi-length) 支持

在之前的版本中，Decoy B 管线存在严重的逻辑缺陷。我们针对性地进行了两大修复，大幅提高了筛选的科学性与准确度。

### 2.1 候选池扩充：解除同长度限制
**问题**: 之前的初筛管线存在 `len(candidate) == len(target)` 的硬性限制。对于靶标 9-mer，直接丢失了高达 40.9%（约 52 万条）的 8-mer、10-mer 和 11-mer 候选。
**修复**: 我们去除了同长度过滤限制。因为 Atchley 理化因子提取的是 TCR 接触核心（中央 5 个残基，25 维向量），天然支持跨长度的余弦相似度计算。
**结果**: 候选池从 75.7 万条 9-mer 扩大到完整的 **128 万条 8-11mer**（HLA-I 结合范围），防止漏报长短肽导致的交叉反应。

### 2.2 RMSD 算法修正：TCR 核心构象对比
**问题**: 之前的结构相似度直接计算了整个 pMHC 复合物（310个 CA 原子）的 RMSD。由于 MHC 重链（201 CA）和 β2m（100 CA）极为庞大且保守，严重稀释了多肽（9 CA）的构象差异。导致所有多肽的结构打分都集中在 `~0.15 Å` 的极低区间，失去了区分度。
**修复**: 我们重写了 RMSD 对比逻辑：
1.  **全局对齐**: 首先在 MHC 骨架 (Chain M + Chain N, 301 CA) 上做 Superimpose 空间对齐。
2.  **局部计算**: **同长度**多肽只计算全长多肽 (Chain P) 的 CA 原子 RMSD；**跨长度**多肽（例如 9-mer vs 10-mer）仅提取并对比**TCR 接触核心**（中央 5 个 CA 原子）的 RMSD。
**结果**: 修复后，肽段的 RMSD 呈现出 `0.18 Å ~ 3.0+ Å` 的真实分布，具有了极强的筛选鉴别力。

---

## 3. 模型部署状态与数据源确认

### 3.1 真实数据源切入
将之前代码中误用的 1万条 demo 数据，正式切换至基于 mhcflurry 生成的完整全人类 HLA 呈递肽段库：
*   **文件**: `data/decoy_a/hla_filtered_HLA-A0201_presentation.parquet`
*   **总量**: 1,280,221 条 (经过 `presentation_percentile ≤ 2.0` 门控)
*   **长度分布**: 
    *   8-mer: 19,592 (1.5%)
    *   9-mer: 757,053 (59.1%)
    *   10-mer: 364,839 (28.5%)
    *   11-mer: 138,737 (10.8%)

### 3.2 预测模型就绪
*   **tFold (主轴)**: 使用 `esm_ppi_650m_tcr.pth` 配合 `tfold_pmhc_trunk.pth`。部署成功，单张 GPU 生成 1 个 pMHC 结构耗时约 2 秒。
*   **ProteinMPNN (逆向)**: Vanilla weight 挂载成功。
*   **AlphaFold 3 (精筛)**: 100GB 级别的模型权重已存在于 `/share/liuyutian/alphafold3`。留作未来高精度验证。

---

## 4. 靶标测试结果 (Stepwise Demo)

我们使用流感 M1 经典靶标 **GILGFVFTL / HLA-A\*02:01** 运行了 `scripts/run_decoy_b_stepwise.py` 完整测试。

### 4.1 Decoy A 筛查情况
*   **机制**: Hamming 距离 $\le 2$ （"长得像的"）。
*   **结果**: 在同长度 75.7万 条候选中，仅命中 **2** 条多肽归属 Decoy A 领地。

### 4.2 Decoy B 筛查情况
剩余的 128 万条序列进入 Decoy B，寻找序列差异大但结构/电荷相似的“替身”。
*   **Stage 1 (理化初筛)**: 按 Atchley 因子计算 TCR 接触核心余弦相似度，过滤得到 Top 5000。
*   **Stage 2 (结构预测)**: tFold 在 40 秒内完成了 11 条代表性多肽（含跨长度序列）的三维结构生成。
*   **Stage 3 & 4 (综合评分)**: 以靶标 `pmhc_GILGFVFTL_HLA-A0201.pdb` 为基准。

**最终 Top 5 排序样本 (摘录)**：

| Rank | Sequence | Len | HD | CosSim | Core/Pep RMSD (Å) | Combined Score | Gene |
|------|----------|-----|----|--------|-------------|-------|------|
| 1 | LVLGFVFML | 9 | 3 | 1.000 | 0.179 | 0.989 | SLC39A9 |
| 2 | AYLGFVFYL | 9 | 3 | 1.000 | 0.270 | 0.974 | SLC11A2 |
| 3 | ASLGFVFSA | 9 | 4 | 1.000 | 0.500 | 0.900 | SLC16A11 |
| 4 | FLLGFVIMP | 9 | 5 | 0.987 | 0.664 | 0.862 | TAAR3P |
| 5 | LVLGFVFMLL | 10 | ×10 | 1.000 | 2.148* | 0.698 | SLC39A9 |

*(注：带有 `*` 为跨长度 10-mer 与靶标 9-mer 进行中央 5 CA 原子的 RMSD 对比，显示其核心构象偏离达到 2.148Å)*

**完整输出文件路径**:
本次 Stepwise Demo 产生的所有中间过程与最终结果文件均已保存在以下目录中，供 PI 详细查阅完整列表和三维结构：
*   **输出主目录**: `data/decoy_b/stepwise_demo/`
*   **Stage 1 初筛候选 (Top 5000)**: `data/decoy_b/stepwise_demo/stage1_atchley_candidates.csv`
*   **Stage 2 pMHC 预测结构**: `data/decoy_b/stepwise_demo/tfold_pdbs/` (包含 11 个 `.pdb` 结构文件)
*   **Stage 3 结构比对详情**: `data/decoy_b/stepwise_demo/stage3_structure_comparison.json`
*   **Stage 4 最终风险排名**: `data/decoy_b/stepwise_demo/final_decoy_b_results.json` (包含完整的 JSON 格式报告)

**分析**:
可以看到，筛选出的多肽（如 LVLGFVFML）与靶标的序列相差 3 个残基（Hamming=3），属于传统的安全盲区。但其 TCR 接触核心 `LGFVF` 与靶标极度相似，且肽段三维构象 RMSD 仅为 0.179Å。此类样本具备极高的结构致敏风险，证明了 Decoy B 管线运行极其有效。

---

## 5. 下一步建议 (Next Steps)

提交 PI 审核，请指示以下行动方向：

1.  **结构筛选规模扩大**: 是否对前 5000 个 Atchley Hit 全部运行 tFold 批量结构生成？（预期耗时：约 3 小时）。
2.  **AlphaFold 3 校验**: 针对 tFold 选出的 Top 100 高危结构，是否引入 AlphaFold 3 进行复核，以抵消轻量级模型可能的折叠偏差？
3.  **计算资源调度**: 如果进入万级别高通量筛选（如批处理多个靶标），当前的单卡能否满足需求，是否需要配置多进程推理队列？

## 6. 新增优化记录 (2026年4月2日晚)

1. **Decoy A 筛选增强**：将 Hamming 距离条件扩展至 `≤ 4`，并在搜索时直接过滤了 `mhcflurry` 中被判定为 `Non_Binder` 的无效多肽，大幅度提升了高维扫描效率，并成功跑通靶标 `GILGFVFTL`（590多条同源匹配）。
2. **Decoy D 模块化独立**：将原本作为 Decoy B 附加分支的 MPNN 逆向设计管线完全独立为 **Decoy D**。支持独立调用 ProteinMPNN 固定锚定位进行靶向设计，结合 mhcflurry 亲和力门控过滤，并利用 tFold 为 Top-K AI 生成序列预测真实 3D 结构，形成自洽且完整的纯合成风险预测分支。
3. **结构可视化比对工具 (`visualize_3d.py`)**：新增轻量级结构比对查看工具，采用 3Dmol.js，可在一键生成的 HTML 文件中同时渲染 Target 骨架与批量预测 Candidate 肽段，直观验证 TCR 结合核心的结构异同。
4. **完整样本聚合输出**：已针对 `GILGFVFTL` (HLA-A*02:01) 真实执行完整评估流，并将 Decoy A、B、D 的核心输出文件（含表格排名及 `.pdb` 结构模型和 HTML 可视化报告）统一归档至 `data/GILGFVFTL_summary/` 目录下。

---

## 7. 全管线 Bug 修复 (2026年4月2日)

对四条管线进行了系统性代码审查，修复了以下关键缺陷：

### Decoy A
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | Non_Binder 过滤只检查 `presentation_binding` 列（不存在），增加 `binding` 列检查 | `scanner.py:213` |
| **严重** | 基因注释 LEFT JOIN 改为 INNER JOIN + drop_duplicates，防止数据丢失和行膨胀 | `hla_filter.py:377` |
| 中等 | mhcflurry 优先使用 `Class1PresentationPredictor` 获取真实 presentation_score | `tools/mhcflurry.py` |
| 中等 | 排序逻辑区分"缺失表达数据"(unknown) 和 TPM=0 (confirmed low risk) | `scanner.py:313` |

### Decoy B
| 严重性 | 修复 | 文件 |
|--------|------|------|
| 高 | BOTH 来源 risk 计算改为加权组合 (0.5×seq + 0.5×physchem)，不再 max() 混用不同度量 | `risk_scorer.py:195` |
| 高 | `get_gene_expression()` 返回 None 时的 NoneType crash | `scanner.py:596` |
| 中等 | EL rank 倒数增加上限 clamp (100.0) 防止极端 risk score | `risk_scorer.py:121` |
| 中等 | 结构数据只要存在就使用 (`structural is not None`)，不再要求 correlation > 0 | `scanner.py:839` |

### Decoy C
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | `PaperRecord.pmid` 对非 PubMed 来源返回 `source:id` 格式，防止 PMID 字段污染 | `multi_source_fetcher.py:57` |
| **严重** | 去重改为 `(sequence, hla_allele)` 双键，同序列不同 HLA 不再被误删 | `models.py:221` |
| 高 | LLM 失败时回退到 `extract_rule_based()`，不再 crash | `extractor.py:364` |

### Decoy D
| 严重性 | 修复 | 文件 |
|--------|------|------|
| **严重** | anchor 位点改用 HLA-allele 特异查找表，而非 `[2, len(seq)]` | `scanner.py:44` |
| **严重** | `_parse_pdb()` return 后的死代码删除，增加输出验证 | `tools/proteinmpnn.py:82` |
| 高 | tFold 失败的预测 (pdb_path=None) 不再写入结果 CSV | `scanner.py:91` |
| 中等 | 替换 stale orchestrator.py（原文件有不存在模块的 import） | `orchestrator.py` |

---

## 8. 双叠合结构比较方法 (Dual Superposition)

### 背景

原有的 Decoy B 结构比较只使用一种策略：在 MHC 上叠合后测量肽段 RMSD (Method A)。我们新增了一种互补策略：在肽段上叠合后测量 MHC groove 螺旋的 RMSD (Method B)，用于验证两者是否给出一致结论。

### 方法

| | Method A（原有） | Method B（新增） |
|--|----------------|----------------|
| **叠合基准** | MHC 骨架 (Chain M+N, 301 CA) | 肽段骨架 (Chain P CA) |
| **测量对象** | 肽段 RMSD (Chain P) | MHC groove 螺旋 RMSD (α1: res 50-85, α2: res 138-175) |
| **物理含义** | 在同一个 HLA 凹槽中，肽段的构象偏差有多大 | 在同一个肽段构象下，MHC groove 的包裹方式偏差有多大 |
| **归一化** | $\max(0, 1 - RMSD/3.0)$ | $\max(0, 1 - RMSD/3.0)$ |

最终 StructSim = 两种方法 similarity 的平均值。

### 验证结果 (GILGFVFTL, 50 个 decoy)

**整体一致性极高**：
- Pearson r = **0.9514** (p = 3.4e-26)
- Spearman ρ = **0.8955** (p = 1.7e-18)
- Kendall τ = **0.7318** (p = 1.3e-12)

Rank 31-50 **完全一致**（差的就是差的，没有分歧）。

### 有意义的排名分歧

**A 好 B 差**（肽段近似但 groove 偏差大）：

| Sequence | Rank A | Rank B | ΔRank | 解释 |
|----------|:------:|:------:|:-----:|------|
| YLLGFLFNY | 3 | 15 | +12 | 肽段骨架位置相近，但 groove 螺旋有扭转 |
| LLLGILFLI | 8 | 19 | +11 | 同上 |
| ILLGFLFCA | 10 | 21 | +11 | 同上 |
| FLLGFLFWN | 1 | 9 | +8 | 同上 |
| LLVGFVFVV | 7 | 14 | +7 | 同上 |

**B 好 A 差**（groove 保守但肽段位移大）：

| Sequence | Rank A | Rank B | ΔRank | 解释 |
|----------|:------:|:------:|:-----:|------|
| RALGFLIGL | 13 | 1 | -12 | groove 包裹几乎完全一致，但肽段有整体平移 |
| KTLGIVIGV | 21 | 10 | -11 | 同上 |
| KTLGIVMGV | 23 | 12 | -11 | 同上 |
| SALGFVIPA | 24 | 13 | -11 | 同上 |
| YGLGIVFPI | 9 | 2 | -7 | 同上 |

### 结论

单独使用 Method A 会漏掉 groove 上下文保守的候选（肽段平移但整体 TCR 呈递表面保守），单独使用 Method B 会漏掉肽段构象保守但 groove 有偏转的候选。**取平均值可同时捕获两类风险**。

验证脚本：`scripts/compare_superposition_methods.py`
验证图表：`figures/superposition_method_comparison.png`
