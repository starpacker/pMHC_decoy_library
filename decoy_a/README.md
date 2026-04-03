# Decoy A — 序列同源扫描

**状态**: 全管线已部署完毕，端到端验证通过

## 核心思想

在全人类蛋白质组中搜索与靶标肽段序列相似的"近亲"（Hamming distance ≤ 2），评估其被 HLA 呈递后引发 TCR 交叉反应的风险。

## 管线架构

```
Step 1: K-mer 字典构建 (一次性)
  Swiss-Prot 20,431 蛋白 → 滑动窗口 → 41.7M 唯一 8-11mer
  + HPA 组织表达谱 (20,151 基因 × 50 组织)

Step 2: HLA 呈递门控 (一次性/每个 HLA)
  41.7M k-mer → mhcflurry/NetMHCpan → 1.28M binders (HLA-A*02:01)
  阈值: EL%Rank ≤ 2.0

Step 3: Hamming 距离扫描 (秒级/每个靶标)
  向量化比较靶标 vs 1.28M binders → 过滤 Hamming ≤ 2
  标注: 锚定位 (p2,p9) vs TCR 接触面 (p4-p8) 错配
  关联: 组织表达风险 (重要脏器 TPM > 10 = 高风险)
```

## 风险评分

$$Risk_A = \left(1 - \frac{Hamming}{Length}\right) \times \frac{1}{EL\_Rank} \times TPM\_Weight$$

## 使用方法

```bash
python -m decoy_a scan --target GILGFVFTL --hla "HLA-A*02:01"
```

## 文件清单

| 文件 | 作用 |
|------|------|
| `scanner.py` | Hamming 距离扫描 + 错配标注 + 表达谱关联 |
| `kmer_builder.py` | Swiss-Prot → 41.7M k-mer + HPA 表达谱 |
| `hla_filter.py` | HLA 门控（三级回退: NetMHCpan > mhcflurry > IEDB） |
| `models.py` | Pydantic v2 数据模型（A/B/D 共用） |
| `config.py` | 全局配置（路径、阈值、Atchley 因子） |
| `tools/` | mhcflurry + NetMHCpan 封装 |

## 数据文件

| 文件 | 大小 |
|------|------|
| `data/decoy_a/human_kmer_db.parquet` | 380 MB |
| `data/decoy_a/gene_expression.parquet` | 1.9 MB |
| `data/decoy_a/hla_filtered_HLA-A0201.parquet` | 27 MB |

## 验证结果

- **GILGFVFTL** (流感 M1): 发现 TAP2 (GLLGFVGTL, HD=2), PKD1L1 (GILLFLFTL, HD=2)
- **RMFPNAPYL** (WT1): 0 hits — 正确，无 Hamming ≤ 2 近亲
