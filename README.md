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
     │   序列同源扫描    │   │  结构相似筛选     │   │   文献+IEDB采掘    │
     │                  │   │                  │   │                  │
     │  "长得像的"       │   │  "形状像的"       │   │ "文献报道过的"     │
     │  Hamming ≤ 4     │   │  双叠合 RMSD      │   │  三重验证+硬拒绝     │
     │  50-200 条/靶标  │   │  + 5维界面描述符   │   │  817 条入库   │
     │                  │   │  100-500 条/靶标  │   │                  │
     │  确定性: 高       │   │  确定性: 中       │   │  确定性: 最高      │
     │  覆盖面: 窄       │   │  覆盖面: 广       │   │  覆盖面: IEDB扩展      │
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

### 统一入口（推荐）

```bash
# 单策略运行
python run_decoy.py GILGFVFTL d                          # Decoy D only
python run_decoy.py GILGFVFTL a                          # Decoy A only

# 多策略组合
python run_decoy.py GILGFVFTL a b d                      # A + B + D

# 指定 HLA / MPNN 参数
python run_decoy.py EVDPIGHLY d --hla HLA-A*01:01
python run_decoy.py GILGFVFTL d --designs 2000 --top-k 20

# 批量运行全部候选靶标 Decoy D + 生成图表
bash scripts/batch_decoy_d.sh
```

### 模块入口

```bash
# Decoy A — 序列同源扫描
python -m decoy_a scan --target GILGFVFTL --hla "HLA-A*02:01"

# Decoy B — 结构相似筛选（完整管线）
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# Decoy D — MPNN 逆向设计
python -m decoy_d --target GILGFVFTL --hla "HLA-A*02:01"

# Decoy C — 文献+IEDB采掘（三重验证）
python -m decoy_c stats
python -m decoy_c show
python scale_up.py -n 10000 --consensus    # 扩库（含IEDB采掘+双次共识提取）
```

---

## 风险评分模型

### Decoy A
$$ Risk_A = \left(1 - \frac{\text{Hamming}}{\text{Length}}\right) \times \frac{1}{EL\_Rank} \times TPM\_Weight $$

### Decoy B
$$ Combined = StructSim + CV\_Boost $$

其中 `CV_Boost = 0.10 × cross_validation_agreement`（Boltz-2 交叉验证一致性加成，最高 +10%）。

> **v2 变更**: Atchley cosine 从最终评分移除（solo rank 分析 ρ≈0 与结构排名无关），仅保留在 Stage 1 物化筛选。

**StructSim** 由双叠合 RMSD 和 **TCR-facing 表面描述符** 综合得出：

$$StructSim = 0.50 \times RMSD\_Geo + 0.50 \times TCR\_Combined$$

| 组件 | 计算方式 |
|------|---------|
| **RMSD_Geo** | avg(Method A: MHC叠合→肽段RMSD, Method B: 肽段叠合→groove RMSD) |
| **TCR_Combined** | 0.30×ESP + 0.30×Shape + 0.20×SASA + 0.20×Hydrophobicity |

4 个 **TCR-facing** 表面描述符（度量 pMHC 暴露给 TCR 的上表面）：

| 描述符 | 度量方式 | 权重 | 捕捉维度 |
|--------|---------|------|---------|
| TCR-facing ESP | 暴露原子 Coulomb 势 → Hodgkin SI | 0.30 | TCR 看到的电荷分布 |
| Exposed Shape | 暴露侧链质心 Kabsch RMSD | 0.30 | TCR 接触的 3D 形状 |
| rSASA Profile | per-position 相对 SASA → Cosine | 0.20 | 哪些位点面向 TCR |
| Exposed Hydrophobicity | rSASA 加权 Kyte-Doolittle → Cosine | 0.20 | TCR 面的疏水性 |

补充描述符 **PLIP**（Protein-Ligand Interaction Profiler）——度量 peptide-MHC **结合槽界面**的非共价相互作用：

| 描述符 | 度量方式 | 用途 |
|--------|---------|------|
| PLIP Tanimoto | H-bond/疏水/盐桥/π-stacking 计数向量 → 广义 Tanimoto 系数 | 相似结合模式 → 相似呈递 → 更高交叉反应风险 |

> PLIP 与 TCR-facing 描述符互补：TCR-facing 度量 TCR 看到的上表面，PLIP 度量 peptide 锚定在 MHC groove 中的结合模式。

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
| `plip` | PLIP 非共价相互作用分析（需 conda 环境 `plip_env` + OpenBabel） |
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
│   │   ├── tcr_surface_descriptors.py #  4 TCR-facing 表面描述符 (v2)
│   │   ├── interface_descriptors.py #   PLIP + groove 界面描述符
│   │   ├── tfold.py                 #   tFold 封装
│   │   ├── alphafold3.py            #   AF3 封装
│   │   ├── boltz.py                 #   Boltz-2 封装
│   │   └── proteinmpnn.py           #   ProteinMPNN 封装
│   └── DEPLOYMENT.md                #   外部工具部署指南（含 PLIP）
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
| **Decoy C** | ✅ 完成 | 817 条入库 (662 VALIDATED, 81% 蛋白确认, 94% IEDB 匹配) |
| **Decoy D** | ✅ 完成 | MPNN 逆向设计 + mhcflurry 过滤 |

详细技术报告见 → [`progress_and_report.md`](progress_and_report.md)

---

## 参考文献与工具来源 (References & Tools)

本项目集成了多种前沿的计算生物学工具和数据库，以下是核心依赖的文献与代码仓库链接：

### 1. 结构预测与逆向设计 (Structure Prediction & Design)
| 工具 | 论文 | 链接 |
|------|------|------|
| **tFold** | Wu L et al. "tFold-Ab/Ag: Fast and Accurate Antibody and Antigen Structure Prediction." 2024 | [GitHub](https://github.com/TencentAI4S/tfold) |
| **AlphaFold 3** | Abramson J et al. "Accurate structure prediction of biomolecular interactions with AlphaFold 3." *Nature* 630:493–500, 2024 | [GitHub](https://github.com/google-deepmind/alphafold3) |
| **Boltz-2** | Wohlwend J et al. "Boltz-1: Democratizing Biomolecular Interaction Modeling." 2024 | [GitHub](https://github.com/jwohlwend/boltz) |
| **ProteinMPNN** | Dauparas J et al. "Robust deep learning–based protein sequence design using ProteinMPNN." *Science* 378(6615):49–56, 2022 | [GitHub](https://github.com/dauparas/ProteinMPNN) |

### 2. 免疫信息学与序列分析 (Immunoinformatics & Sequence Analysis)
| 工具/方法 | 论文 | 链接 |
|-----------|------|------|
| **mhcflurry** | O'Donnell TJ et al. "MHCflurry 2.0: Improved Pan-Allele Prediction of MHC Class I-Presented Peptides..." *Cell Systems* 11(1):42-48.e7, 2020 | [GitHub](https://github.com/openvax/mhcflurry) |
| **NetMHCpan** | Reynisson B et al. "NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation..." *Nucleic Acids Res.* 48(W1):W449-W454, 2020 | [Web](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) |
| **Atchley Factors** | Atchley WR et al. "Solving the protein sequence metric problem." *PNAS* 102(18):6395–6400, 2005 | [DOI](https://doi.org/10.1073/pnas.0408677102) |
| **BioPython** | Cock PJA et al. "Biopython: freely available Python tools for computational molecular biology..." *Bioinformatics* 25(11):1422–1423, 2009 | [GitHub](https://github.com/biopython/biopython) |

### 3. 数据库 (Databases)
| 数据库 | 论文 | 链接 |
|--------|------|------|
| **IEDB** | Vita R et al. "The Immune Epitope Database (IEDB): 2018 update." *Nucleic Acids Res.* 47(D1):D339-D343, 2019 | [Web](https://www.iedb.org/) |
| **UniProt** | The UniProt Consortium. "UniProt: the Universal Protein Knowledgebase in 2023." *Nucleic Acids Res.* 51(D1):D523-D531, 2023 | [Web](https://www.uniprot.org/) |

### 4. 关键背景文献 (Key Literature)
| 主题 | 论文 | 链接 |
|------|------|------|
| **Titin 致死案例** | Cameron BJ et al. "Identification of a Titin-derived HLA-A1–presented peptide as a cross-reactive target..." *Sci. Transl. Med.* 5(197):197ra103, 2013 | [DOI](https://doi.org/10.1126/scitranslmed.3006034) |
| **TCR 交叉反应性综述** | Sewell AK. "Why must T cells be cross-reactive?" *Nature Rev. Immunol.* 12:669–677, 2012 | [DOI](https://doi.org/10.1038/nri3279) |

---

## 常见问题 (Q&A)

**Q: 为什么需要四个不同的 Decoy 管线？**
A: TCR 交叉反应性是一个复杂的多维度问题。Decoy A 解决序列高度相似的"近亲"脱靶；Decoy B 解决序列不同但三维结构和理化性质相似的"远亲"脱靶（如著名的 MAGE-A3/Titin 致死案例）；Decoy C 提供真实世界已发生的临床和实验证据作为基准；Decoy D 则通过 AI 逆向设计探索未知的理论风险空间。四者互补，构成严密的防御网。

**Q: 风险评分模型中的 TPM_Weight 是如何工作的？**
A: 即使某个多肽与靶标极其相似，如果它在人体内根本不表达，或者只在非重要组织低水平表达，其引发致命副作用的风险也较低。TPM_Weight 引入了人类蛋白质图谱（HPA）的表达数据，对在心脏、大脑等致命器官中高表达（TPM > 10）的基因给予 10 倍的风险加权惩罚，从而将计算资源和注意力集中在真正具有临床危险的脱靶候选上。
