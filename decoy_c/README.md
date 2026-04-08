# Decoy C — 文献挖掘 + IEDB 系统性采掘

**状态**: 全管线运行完成，**817 条记录入库**（662 条 VALIDATED），三重验证 + 蛋白子序列确认 + 硬拒绝机制已上线

## 核心思想

从 PubMed 已发表论文 + IEDB 数据库中系统性提取所有已知的 TCR 交叉反应案例，构建有实验证据支撑的参考数据库——回答"哪些脱靶已经在实验或临床中被证实？"

## 管线架构 (v2 — 三重验证)

```
Stage 1: 多源文献发现 (fetcher.py + multi_source_fetcher.py)
  376+ 关键词 × PubMed/arXiv/bioRxiv/ClinicalTrials → papers

Stage 2: LLM 结构化提取 (extractor.py)
  Title + Abstract → GPT-4o → JSON array of DecoyEntry
  NEW: 双次提取共识模式 (--consensus) — 2 次 LLM 调用，仅保留两次一致的序列
  NEW: 源文本验证 — 检查提取的序列是否真正出现在原论文中

Stage 3: 三重外部验证 (validator.py)
  ① UniProt 基因验证 (gene_symbol → accession)
  ② Peptide-in-Protein 子序列验证 (从 UniProt 拉全长蛋白，验证 peptide ⊂ protein)
  ③ IEDB 表位数据库交叉验证

Stage 4: 硬拒绝门控 (orchestrator.py)
  - protein_containment = NOT_FOUND → 硬拒绝（序列不在蛋白中 = 幻觉）
  - 源文本 + IEDB + 蛋白三重失败 → 硬拒绝（零外部证据）

Stage 5: IEDB 直接采掘 (iedb_miner.py) — NEW
  4 大策略系统性扫描 IEDB：
    ① 危险蛋白查询 (TTN, MAGEA3, cardiac, neural etc.)
    ② 人类自身抗原 × 常见 HLA 交叉查询
    ③ 质谱确认的天然呈递肽
    ④ 自身免疫/交叉反应标注记录

Stage 6: 库管理 (orchestrator.py)
  去重 → 分配 DC-NNNN ID → 质量过滤 → 持久化 JSON
```

## 验证层级

| 验证 | 方法 | 防御目标 |
|------|------|---------|
| 源文本匹配 | 检查序列是否出现在原文中 | LLM 幻觉序列 |
| 蛋白子序列 | 从 UniProt 拉全长蛋白做 `peptide in protein_seq` | LLM 编造不存在的肽段 |
| IEDB 交叉验证 | IEDB epitope search API | 确认实验证据 |
| 硬拒绝门控 | 三重验证全失败 → 拒绝入库 | 无任何证据支撑的条目 |

## 证据分级 (4-Tier)

| 等级 | 名称 | 示例 |
|------|------|------|
| Level 1 | Clinical Fatal | MAGE-A3/Titin 心脏致死 |
| Level 2 | In Vitro Confirmed | Cr-release, IFN-gamma assay |
| Level 3 | HT Screened | X-scan, yeast display |
| Level 4 | In Silico | ARDitox, EpiTox 预测 |

## 使用方法

```bash
python -m decoy_c stats    # 库统计
python -m decoy_c show     # 查看条目

# 扩库 (含三重验证 + IEDB 采掘)
python scale_up.py -n 10000 --strategy all --consensus

# 仅 IEDB 采掘
python -c "from decoy_c.iedb_miner import mine_iedb; from decoy_c.orchestrator import load_library; mine_iedb(load_library())"
```

## 文件清单

| 文件 | 作用 |
|------|------|
| `fetcher.py` | PubMed 文献获取 |
| `extractor.py` | LLM 结构化提取 + 源文本验证 + 双次共识 |
| `validator.py` | UniProt + 蛋白子序列 + IEDB 三重验证 + 硬拒绝 |
| `orchestrator.py` | 管线编排 + 硬拒绝门控 + 去重 |
| `iedb_miner.py` | IEDB 系统性采掘 (4 策略) + UniProt ID 提取 |
| `revalidate_and_mine.py` | 旧 schema 迁移 + 全库重验证 + IEDB 采掘 |
| `fix_and_revalidate.py` | UniProt ID 修复 + 蛋白子序列重验证 + 清退 |
| `models.py` | Pydantic 数据模型 |
| `scale_up.py` | 自动化扩库 + IEDB 集成 |
| `split_library.py` | 库分类 (strict / broad / removed) |
| `multi_source_fetcher.py` | arXiv/bioRxiv/ClinicalTrials 多源获取 |

## 标志性案例: MAGE-A3 / Titin

EVDPIGHLY (MAGE-A3) → ESDPIVAQY (Titin), HLA-A*01:01. 亲和力增强型 TCR 导致 2 名患者心脏骤停死亡 (2013, Adaptimmune)。Hamming distance = 5，远超 Decoy A 阈值，但 TCR 接触面高度相似。

## 当前库统计 (2026-04-07)

| 指标 | 数值 |
|------|------|
| 总条目 | 817 |
| VALIDATED (双重验证通过) | 662 (81.0%) |
| 蛋白子序列确认 | 667 (81.6%) |
| IEDB 匹配 | 768 (94.0%) |
| 质谱确认 | 476 (58.3%) |
| HLA allele 覆盖 | 46 个 |
| 基因覆盖 | 310 个 |

## 常见问题 (Q&A)

**Q: 三重验证的通过率大约是多少？**
A: 实测 IEDB 采掘条目 92.5% 通过蛋白子序列验证（667/720），94% 有 IEDB 匹配。LLM 提取条目约 60-70% 通过。

**Q: IEDB 采掘能带来多少新条目？**
A: 实测单次采掘（11 个 HLA allele × Homo sapiens 阳性记录）获得 757 候选，最终 712 条入库。主要来源是 T 细胞检测阳性的人类自身抗原表位。

**Q: 为什么需要 peptide-in-protein 验证？**
A: GPT-4o 在提取时可能"幻觉"出不存在的肽段序列——看起来合理但实际不在蛋白质中。从 UniProt 拉取全长蛋白序列做子串匹配，是防止虚假条目入库的最可靠手段。IEDB 采掘中也有约 4.6% (33/712) 的条目被此验证淘汰。
