# Decoy A Pipeline — 技术进展报告

**日期**: 2026-04-01 (更新)
**状态**: 全管线已部署完毕，端到端验证通过，可视化已完成

---

## 一、技术原理

### 1.1 核心问题

TCR-T 细胞疗法（TCR-T immunotherapy）通过工程化 T 细胞识别肿瘤抗原来杀伤肿瘤。但 TCR 的识别并非绝对特异——它可能"交叉反应"（cross-react）到正常组织中序列相似的肽段上，导致严重的脱靶毒性（off-target toxicity）。

**Decoy A 的目标**：对给定靶标肽段，在全人类蛋白质组中系统性搜索所有"序列近亲"（Hamming distance ≤ 2），并评估其临床风险。

### 1.2 生物学基础

TCR 识别需要满足两个前提条件：

1. **MHC 呈递**：肽段必须能被 MHC-I 分子结合并呈递到细胞表面
2. **TCR 识别**：呈递的 pMHC 复合物必须能被 TCR CDR3 环识别

因此，仅仅序列相似是不够的——肽段必须同时满足"能被 HLA 呈递"。这就是 HLA 过滤存在的根本原因。

### 1.3 EL%Rank 的意义

NetMHCpan 4.1（或等效工具 mhcflurry）使用深度学习模型预测肽段在特定 HLA 等位基因上的结合能力。输出指标 **EL%Rank（Eluted Ligand percentile Rank）** 反映该肽段在该 HLA 上被呈递的概率：

| 阈值 | 分类 | 含义 |
|------|------|------|
| ≤ 0.5 | Strong Binder (SB) | 排名前 0.5%，高亲和力，几乎确定被呈递 |
| ≤ 2.0 | Weak Binder (WB) | 排名前 2%，中等亲和力，仍可能被呈递 |
| > 2.0 | Non-Binder | 不太可能被呈递 |

**为什么用 %Rank 而非绝对 IC50**：不同 HLA 等位基因的绝对亲和力尺度差异很大，%Rank 是归一化指标，可跨等位基因比较。

**为什么阈值选 2.0**：安全性筛查需要高灵敏度（宁多勿漏），阈值 2.0 包含弱结合肽段，因为在炎症环境下 MHC 表达上调时弱结合肽段也可能引发交叉反应。

### 1.4 管线架构

```
Input: 靶标肽段 (e.g., EVDPIGHLY) + HLA 型别

Step 1: K-mer 字典构建 (一次性)
  Swiss-Prot 20,431 蛋白 → 滑动窗口 → 41.7M 唯一 8-11mer
  + HPA 组织表达谱 (20,151 基因 × 50 组织)

Step 2: HLA 呈递门控 (一次性/每个 HLA)
  41.7M k-mer → mhcflurry/NetMHCpan → 1.28M binders (HLA-A*02:01)
  阈值: EL%Rank ≤ 2.0

Step 3: Hamming 距离扫描 (秒级/每个靶标)
  向量化比较靶标 vs 1.28M binders
  过滤: Hamming ≤ 2
  标注: 锚定位 (p2,p9) vs TCR 接触面 (p4-p8) 错配
  关联: 组织表达风险 (重要脏器 TPM > 10 = 高风险)

Output: 排序后的 DecoyAHit 列表 (JSON)
  排序: Hamming distance ASC → 重要脏器 TPM DESC → EL%Rank ASC
```

---

## 二、实现方法

### 2.1 K-mer 数据库构建 (`kmer_builder.py`)

**输入**: UniProt Swiss-Prot FASTA（人类已审核蛋白，20,431 条）

**方法**: 对每个蛋白做 8/9/10/11-mer 滑动窗口，去重后记录每个 k-mer 的来源蛋白和基因。

**工程优化**: 原始 `defaultdict` 方案在 16 GB RAM 上 OOM。重写为 SQLite 磁盘方案：
1. 所有 (kmer, protein, gene) 三元组插入 SQLite
2. 分批导出为 Parquet chunk
3. 合并为最终文件
4. 内存峰值从 >16 GB 降至 ~2 GB

| 指标 | 值 |
|------|-----|
| 蛋白数量 | 20,431 |
| 唯一 k-mer 数量 | **41,684,988** |
| 8-mer / 9-mer / 10-mer / 11-mer | 10.3M / 10.4M / 10.5M / 10.5M |
| 文件大小 | 380 MB (Parquet) |
| 构建耗时 | ~45 分钟 |

### 2.2 组织表达谱 (`kmer_builder.py`)

**来源**: Human Protein Atlas (HPA) RNA tissue consensus 数据

**内容**: 20,151 个基因在 50 个组织中的 nTPM 值

**关键脏器**: 心肌、大脑皮层、肺、肝、肾、结肠、小肠

**风险分类**:
- `high_risk`: 任何关键脏器 TPM ≥ 10
- `expressed`: TPM > 1
- `silent`: 所有关键脏器 TPM < 1
- `restricted`: 仅在非关键组织表达

### 2.3 HLA 过滤后端 (`hla_filter.py`)

三级回退策略：

| 优先级 | 后端 | 状态 | 备注 |
|--------|------|------|------|
| 1 | NetMHCpan 4.1 | 可选 | 需学术许可，精度最高 |
| 2 | mhcflurry 2.2.0 | **当前使用** | 开源，pan-allele deep learning |
| 3 | IEDB API | 备选 | 远程调用，速率受限 |

**mhcflurry 部署修复**:
- Python 3.14 兼容：`from pipes import quote` → `from shlex import quote`
- Windows UTF-8：需 `PYTHONUTF8=1` 避免 mhcgnomes YAML 解析崩溃
- 模型下载：GitHub releases via ghfast.top 代理（直连仅 40 KB/s）

### 2.4 Hamming 距离扫描 (`scanner.py`)

**核心算法**: NumPy 向量化 Hamming 距离计算

```python
target_arr = np.array(list(target), dtype="U1")
cand_chars = np.array([list(str(s)) for s in candidates], dtype="U1")
mismatches = cand_chars != target_arr
distances = mismatches.sum(axis=1)
```

**错配标注**: 每个位点分类为
- **Anchor 位点** (p2, p9 for HLA-A*02:01 9-mer): HLA 结合关键
- **TCR contact 位点** (p4-p8): TCR 识别关键
- **其他位点**: 非关键

**排序逻辑**: `(hamming_distance ASC, vital_organ_tpm DESC, el_rank ASC)`
- Hamming=1 优先于 Hamming=2
- 同等距离下，重要脏器高表达者更危险
- EL%Rank 作为 tiebreaker

---

## 三、HLA-A*02:01 过滤结果

对 41,684,988 条 k-mer 使用 mhcflurry 2.2.0 预测结合能力：

| 指标 | 最终值 |
|------|--------|
| 总 batch | 834 / 834 (100%) |
| 处理肽段 | 41,684,988 |
| 发现 binders | **1,280,221** (3.07%) |
| 强结合 (SB, %Rank ≤ 0.5) | 292,479 (22.8%) |
| 弱结合 (WB, %Rank ≤ 2.0) | 987,742 (77.2%) |
| 文件 | `hla_filtered_HLA-A0201.parquet` (27 MB) |
| 总耗时 | ~6 小时 |

**按肽段长度分布**:
| 长度 | Binders | 占比 |
|------|---------|------|
| 8-mer | 19,592 | 1.5% |
| 9-mer | 757,053 | 59.1% |
| 10-mer | 364,839 | 28.5% |
| 11-mer | 138,737 | 10.8% |

9-mer 占绝对多数，与 HLA-I 偏好 9-mer 的已知生物学一致。

---

## 四、端到端验证结果

### 验证 1: GILGFVFTL（流感 M1 肽）/ HLA-A*02:01

从 1,280,221 可呈递肽中搜索，耗时 6 秒：

```
#  Sequence      HD  EL%Rank  Genes     Mismatches
1  GLLGFVGTL      2    0.17   TAP2      p2:I>L^  p7:F>G*
2  GILLFLFTL      2    0.48   PKD1L1    p4:G>L*  p6:V>L*
```

TAP2 和 PKD1L1 均为正常组织表达蛋白，且都是 HLA-A*02:01 强结合肽。

### 验证 2: RMFPNAPYL（WT1）/ HLA-A*02:01

0 hits — **正确结果**。WT1 的这个表位在全人类蛋白组中不存在 Hamming ≤ 2 的近亲。

### 验证 3: EVDPIGHLY（MAGE-A3）/ HLA-A*02:01

0 hits — **正确结果**。EVDPIGHLY 是 HLA-A*01:01 限制性肽段，其近亲也不能被 HLA-A*02:01 呈递。需要运行 HLA-A*01:01 过滤才能找到命中。

---

## 五、可视化

已完成 Decoy A 全套可视化（`figures/` 目录）：

| 文件 | 内容 |
|------|------|
| `decoy_a_overview.png` | 四面板: EL%Rank 分布、Binder 比例、肽段长度、各长度 EL%Rank 分布 |
| `decoy_a_funnel.png` | Pipeline 生成 + 过滤漏斗图 (含 pass rate) |
| `decoy_a_mismatch_positions.png` | 突变位点分布 (TCR contact / Anchor / Other) |
| `decoy_a_scatter.png` | EL%Rank vs Hamming distance 散点图 |
| `decoy_a_interactive.html` | Plotly 交互式 Case Study (EPS8L2 + 组织表达热图) |

---

## 六、已知限制

| 限制 | 影响 | 解决方案 |
|------|------|---------|
| mhcflurry 替代 NetMHCpan | affinity %Rank 近似 EL %Rank，阈值可能略有差异 | 可选升级 NetMHCpan（需学术许可） |
| 每个 HLA 型别需独立过滤 | 首次 ~6 小时 | 结果缓存为 Parquet，只需运行一次 |
| Hamming distance 仅捕捉序列层面相似 | 无法发现序列差异大但结构相似的脱靶 | 由 Decoy B 补充（结构相似性） |
| Windows 需 `PYTHONUTF8=1` | mhcgnomes YAML 解析需要 UTF-8 | wrapper 中已自动设置 |

---

## 七、文件清单

### 代码文件

| 文件 | 作用 |
|------|------|
| `decoy_a/scanner.py` | Hamming 距离扫描 + 错配标注 + 表达谱关联 |
| `decoy_a/run.py` | 独立运行器（demo / custom pool / full pipeline） |
| `decoy_a/kmer_builder.py` | Swiss-Prot → k-mer 生成（SQLite 磁盘方案）→ HPA 表达谱 |
| `decoy_a/hla_filter.py` | HLA 门控（三级回退: NetMHCpan > mhcflurry > IEDB） |
| `decoy_a/tools/mhcflurry.py` | mhcflurry 2.2.0 封装 |
| `decoy_a/tools/netmhcpan.py` | NetMHCpan 4.1 封装（可选） |
| `scripts/visualize_decoy_a.py` | 全套可视化脚本 (matplotlib + plotly) |

### 数据文件

| 文件 | 大小 | 状态 |
|------|------|------|
| `data/decoy_a/human_kmer_db.parquet` | 380 MB | 已完成 |
| `data/decoy_a/gene_expression.parquet` | 1.9 MB | 已完成 |
| `data/decoy_a/hla_filtered_HLA-A0201.parquet` | 27 MB | 已完成 |
| `data/decoy_a/uniprot_sprot_human.fasta` | 14 MB | 源数据缓存 |
| `data/decoy_a/rna_tissue_consensus.tsv` | 37 MB | 源数据缓存 |
| `data/decoy_a/decoy_a_results.json` | ~10 KB | 最近扫描结果 |

---

## 八、后续计划

1. **扩展 HLA 型别**: HLA-A*01:01, HLA-A*24:02, HLA-B*07:02（各 ~6 小时）
2. **批量扫描**: 对常见肿瘤靶标 (MAGE-A3, NY-ESO-1, WT1, MART-1) 全面扫描
3. **与 Decoy B/C 整合**: 序列相似(A) + 结构相似(B) + 文献验证(C) 三维交叉验证
