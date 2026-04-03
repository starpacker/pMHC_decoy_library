# Decoy C — 文献挖掘

**状态**: 全管线已部署运行，250 条记录入库

## 核心思想

从 PubMed 已发表论文中系统性提取所有已知的 TCR 交叉反应案例，构建有实验证据支撑的参考数据库——回答"哪些脱靶已经在实验或临床中被证实？"

## 管线架构

```
Stage 1: PubMed 文献发现 (fetcher.py)
  8 组冷启动关键词 → NCBI E-utilities → ~200 papers

Stage 2: LLM 结构化提取 (extractor.py)
  Title + Abstract → GPT-4o (T=0.1) → JSON array of DecoyEntry

Stage 3: 外部数据库验证 (validator.py)
  UniProt 基因验证 + IEDB 表位数据库交叉验证

Stage 4: 库管理 (orchestrator.py)
  去重 → 分配 DC-NNNN ID → 持久化 JSON
```

## 证据分级 (4-Tier)

| 等级 | 名称 | 示例 |
|------|------|------|
| Level 1 | Clinical Fatal | MAGE-A3/Titin 心脏致死 |
| Level 2 | In Vitro Confirmed | Cr-release, IFN-gamma assay |
| Level 3 | HT Screened | X-scan, yeast display |
| Level 4 | In Silico | ARDitox, EpiTox 预测 |

## 库统计

| 指标 | 值 |
|------|-----|
| 总条目 | 250 |
| 唯一基因 | 79 |
| 唯一 HLA alleles | 15 |
| Mass-spec 确认 | 62 (24.8%) |
| Level 1 (Clinical Fatal) | 5 |

## 使用方法

```bash
python -m decoy_c stats    # 库统计
python -m decoy_c show     # 查看条目
```

## 标志性案例: MAGE-A3 / Titin

EVDPIGHLY (MAGE-A3) → ESDPIVAQY (Titin), HLA-A*01:01. 亲和力增强型 TCR 导致 2 名患者心脏骤停死亡 (2013, Adaptimmune)。Hamming distance = 5，远超 Decoy A 阈值，但 TCR 接触面高度相似。

## 文件清单

| 文件 | 作用 |
|------|------|
| `fetcher.py` | PubMed 文献获取 |
| `extractor.py` | LLM 结构化提取 |
| `validator.py` | UniProt + IEDB 双重验证 |
| `orchestrator.py` | 管线编排 + 去重 |
| `models.py` | Pydantic 数据模型 |
| `scale_up.py` | 自动化扩库 |
