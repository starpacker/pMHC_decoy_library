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
     │  50-200 条/靶标  │   │  100-500 条/靶标  │   │  105 条（已整理）  │
     │                  │   │                  │   │                  │
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

所有模块均通过 `python -m` 运行。

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
$$ Combined = 0.4 \times CosSim + 0.6 \times StructSim $$

其中 **StructSim** 由双叠合方法综合得出：

$$StructSim = \frac{1}{2}\left[\max\left(0,\ 1 - \frac{PeptideRMSD_A}{3.0}\right) + \max\left(0,\ 1 - \frac{GrooveRMSD_B}{3.0}\right)\right]$$

| 方法 | 叠合基准 | 测量对象 | 回答的问题 |
|------|---------|---------|-----------|
| **Method A** | MHC 骨架 (Chain M+N, 301 CA) | 肽段 RMSD | "在同一个 HLA 凹槽中，肽段构象偏差多大？" |
| **Method B** | 肽段骨架 (Chain P CA) | MHC groove 螺旋 RMSD (α1+α2) | "在同一个肽段构象下，groove 包裹偏差多大？" |

两种方法在 50 个 GILGFVFTL decoy 上验证：**Pearson r = 0.95, Spearman ρ = 0.90**，高度一致但存在有意义的互补差异。详细分析见 `progress_and_report.md`。

### 通用权重因子
- **EL_Rank**: mhcflurry presentation_percentile，越小亲和力越强
- **TPM_Weight**: 心脏/大脑等致命器官 ($TPM > 10$) 有 $10\times$ 加权惩罚

---

## 环境配置

所有依赖部署在 `/share/liuyutian`，路径统一在 `decoy_a/config.py` 管理：

| 工具 | 用途 | 路径 |
|------|------|------|
| tFold | pMHC 结构预测（Decoy B 必选） | `/share/liuyutian/tfold` |
| AlphaFold 3 | 精细结构验证（可选） | `/share/liuyutian/alphafold3` |
| ProteinMPNN | 逆向序列设计（Decoy D） | `/share/liuyutian/S3AI/rebuttal/ProteinMPNN` |
| mhcflurry | HLA 呈递预测 | 已 pip 安装 |

```bash
# 验证工具可用性
python scripts/verify_decoy_b.py
```

---

## 数据来源

### HLA 可呈递肽段库
- **文件**: `data/decoy_a/hla_filtered_HLA-A0201_presentation.parquet`
- **总量**: 1,280,221 条（HLA-A*02:01, presentation_percentile ≤ 2.0）
- **长度分布**: 8-mer 1.5% / 9-mer 59.1% / 10-mer 28.5% / 11-mer 10.8%

---

## 项目结构

```
pMHC_decoy_library/
├── decoy_a/                         # Decoy A: 序列同源扫描 + 共享基座
│   ├── config.py                    #   全局配置（路径、阈值、Atchley 因子）
│   ├── models.py                    #   Pydantic v2 数据模型（A/B/D 共用）
│   ├── kmer_builder.py              #   Swiss-Prot → 41.7M k-mer + 表达谱
│   ├── hla_filter.py                #   HLA 呈递门控（mhcflurry/NetMHCpan）
│   ├── scanner.py                   #   Hamming 序列同源扫描
│   └── tools/                       #   mhcflurry + NetMHCpan 封装
│
├── decoy_b/                         # Decoy B: 结构相似筛选
│   ├── scanner.py                   #   四阶段管线 + 双叠合结构比较
│   ├── risk_scorer.py               #   统一风险评分（合并 A+B）
│   ├── orchestrator.py              #   全管线编排 + 断点续传
│   └── tools/                       #   tFold + AF3 + ProteinMPNN 封装
│
├── decoy_c/                         # Decoy C: 文献挖掘
│   ├── extractor.py                 #   LLM 结构化提取 + 规则回退
│   ├── validator.py                 #   UniProt + IEDB 交叉验证
│   └── ...
│
├── decoy_d/                         # Decoy D: MPNN 逆向设计
│   ├── scanner.py                   #   MPNN + mhcflurry + tFold 管线
│   └── tools/                       #   ProteinMPNN + tFold 封装
│
├── scripts/                         # 工具脚本
│   ├── run_decoy_b_stepwise.py      #   Decoy B 四阶段逐步演示
│   ├── compare_superposition_methods.py  # 双叠合方法对比验证
│   └── visualize_*.py               #   各类可视化
│
├── data/                            # 运行数据与结果
│   ├── decoy_a/                     #   K-mer 库 + HLA 过滤缓存
│   ├── decoy_b/                     #   PDB 结构 + 排名结果
│   ├── decoy_c/                     #   文献库
│   └── GILGFVFTL_summary/           #   GILGFVFTL 完整评估归档
│
├── figures/                         # 可视化图表
├── progress_and_report.md           # 详细技术报告与进展
└── README.md
```

---

## 当前进度

| 管线 | 状态 | 说明 |
|------|------|------|
| **Decoy A** | ✅ 完成 | 全管线 + Presentation Score 升级 |
| **Decoy B** | ✅ 完成 | Atchley 初筛 + tFold + 双叠合 RMSD + 综合评分 |
| **Decoy C** | ✅ 完成 | 250 条入库 (105 有证据等级) |
| **Decoy D** | ✅ 完成 | MPNN 逆向设计 + mhcflurry 过滤 |

详细技术报告、算法修复记录、测试结果见 → [`progress_and_report.md`](progress_and_report.md)
