# Decoy C Library — 技术进展报告

**日期**: 2026-04-01
**状态**: 全管线已部署运行，250 条记录入库，可视化已完成

---

## 一、技术原理

### 1.1 为什么需要 Decoy C

Decoy A（序列同源）和 Decoy B（结构相似）都是**计算预测**——它们回答"理论上可能存在哪些脱靶"。但临床安全性评估需要一个更关键的问题："**哪些脱靶已经在实验或临床中被证实？**"

Decoy C 是一个 **文献挖掘系统**（literature mining system），从 PubMed 已发表论文中系统性提取所有已知的 TCR 交叉反应案例，构建一个有实验证据支撑的参考数据库。

### 1.2 核心设计思想

```
计算预测 (Decoy A/B)  ←→  实验验证 (Decoy C)
      ↓                          ↓
  "可能的脱靶"              "已证实的脱靶"
      ↓                          ↓
         → 交叉验证 → 综合安全性评估
```

**三者的互补关系**:
- Decoy A 发现序列近亲（Hamming ≤ 2）
- Decoy B 发现结构近亲（cosine ≥ 0.70）
- Decoy C 提供实验"金标准"——如果 A/B 的预测命中了 C 库中的已知案例，置信度大幅提升

### 1.3 信息提取的挑战

生物医学文献中的交叉反应数据以自然语言散落在论文的不同段落中：
- 肽段序列可能出现在正文、图表、补充材料中
- HLA 型别格式不统一（"HLA-A2", "HLA-A*02:01", "A*0201"）
- 证据等级需要从实验描述中推断
- 器官毒性信息可能在临床报告的叙述段落中

传统的正则表达式匹配无法可靠处理这种非结构化数据。Decoy C 使用 **大语言模型（LLM）** 作为信息提取引擎，将非结构化文本转换为结构化 JSON。

---

## 二、管线架构

### 2.1 四阶段管线

```
Stage 1: PubMed 文献发现 (fetcher.py)
  ├─ 8 组冷启动关键词查询
  ├─ NCBI E-utilities API (ESearch + EFetch)
  ├─ 获取 Title + Abstract + PMC 全文（如可用）
  └─ Rate limit: 0.34s/request (3 req/sec)
  → ~200 papers

Stage 2: LLM 结构化提取 (extractor.py)
  ├─ 输入: Title + Abstract (截断至 48K 字符)
  ├─ LLM: GPT-4o (T=0.1, strict mode)
  ├─ 输出: JSON array of DecoyEntry
  ├─ Pydantic 严格 enum 验证
  └─ 规则提取后备 (regex, 仅 demo)
  → ~250 raw entries

Stage 3: 外部数据库验证 (validator.py)
  ├─ UniProt: gene symbol → accession (human, reviewed)
  │   状态: OK / ENRICHED / MISMATCH / NOT_FOUND
  ├─ IEDB: peptide sequence → epitope DB
  │   确认: mass-spec elution / assay types / HLA alleles
  └─ 综合: VALIDATED / PARTIAL / NEEDS_REVIEW
  → Validated entries

Stage 4: 库管理 (orchestrator.py)
  ├─ 51 条种子数据 (手工策展)
  ├─ 按 peptide sequence 去重
  ├─ 分配 DC-NNNN 顺序 ID
  └─ 持久化为 JSON
  → decoy_library.json (250 entries)
```

### 2.2 冷启动查询策略

```python
COLD_START_QUERIES = [
    '"TCR" AND "cross-reactivity" AND ("fatal" OR "toxicity" OR "cardiac")',
    '"affinity-enhanced TCR" AND "off-target"',
    '"MAGE-A3" AND "TCR" AND ("Titin" OR "cross-reactive")',
    '"adoptive cell therapy" AND "off-target" AND "peptide"',
    '"engineered T cell" AND "cross-reactivity" AND "safety"',
    '"ARDitox" OR "EpiTox" AND "TCR" AND "safety"',
    '"X-scan" AND "TCR" AND "peptide"',
    '"alanine scanning" AND "TCR" AND "cross-reactivity"',
]
```

### 2.3 证据分级体系 (4-Tier)

| 等级 | 名称 | 含义 | 示例 |
|------|------|------|------|
| Level 1 | Clinical Fatal | 患者死亡或器官衰竭 | MAGE-A3/Titin 心脏致死 |
| Level 2 | In Vitro Confirmed | 细胞毒性实验验证 | Cr-release, IFN-gamma assay |
| Level 3 | High-Throughput Screened | 高通量筛选发现 | X-scan, yeast display, FACS |
| Level 4 | In Silico High Risk | 纯计算预测 | ARDitox, EpiTox 预测 |

---

## 三、实现方法

### 3.1 文献获取 (`fetcher.py`)

**NCBI E-utilities 调用链**:
1. `ESearch`: keyword → PMID list
2. `EFetch`: PMID → XML metadata (title, abstract, authors, MeSH)
3. `PMC OA`: 尝试获取全文（PMC Open Access 子集）

**Rate limiting**: 0.34 秒/请求，无 API key 下 NCBI 允许 ~3 req/sec

### 3.2 LLM 结构化提取 (`extractor.py`)

**模型选择** (按性能排序):

| 模型 | 平均延迟 | 成功率 | 推荐度 |
|------|---------|--------|--------|
| GPT-4o | 1.48s | 100% | ★★★★★ |
| Gemini 2.5 Flash | 2.60s | 100% | ★★★★ |
| Claude Sonnet 4.5 | 3.19s | 100% | ★★★★ |
| GPT-5.2 | 3.99s | 100% | ★★★ |

**System Prompt 设计要点**:
- 要求提取所有 distinct decoy peptides
- 严格 enum 验证 (evidence_level, expression_pattern)
- 反幻觉机制: `thought_process` 字段要求引用原文
- 输出: JSON array，Pydantic 即时验证

**后处理**:
- HLA allele 格式标准化: "HLA-A2" → "HLA-A*02:01"
- 器官名称大写化
- decoy_id 重新分配

### 3.3 外部数据库验证 (`validator.py`)

**UniProt 验证**:
- 查询: `gene:{GENE_SYMBOL} AND organism_id:9606`
- 功能: 验证基因是否存在于人类蛋白组（reviewed=true）
- 输出: 填充 UniProt accession, 标记 mismatch

**IEDB 验证**:
- 查询: `linear_sequence eq.{PEPTIDE_SEQ}`
- 功能:
  - 确认肽段是否在 IEDB 中有记录
  - 检查 mass-spec (elution assay) 确认
  - 提取关联的 assay 类型、HLA alleles、TCR names
- 输出: FOUND/NOT_FOUND + 富集信息

**验证状态**:
| 状态 | 含义 |
|------|------|
| `VALIDATED` | UniProt OK + IEDB 找到 |
| `PARTIAL` | UniProt OK, IEDB 未找到 |
| `NEEDS_REVIEW` | UniProt mismatch, 需人工检查 |

### 3.4 规模化扩展 (`scale_up.py`)

支持自动化扩库至目标数量：

```bash
python scale_up.py --required_num 10000 --strategy all --resume
```

**策略扩展** (A-H):
- A: 临床毒性关键词
- B: 蛋白家族名 (TTN, MAGE, HLA)
- C: 实验方法 (X-scan, Cr-release, mass-spec)
- D: 计算工具 (ARDitox, EpiTox)
- E: TCR 克隆名 (1G4, DMF5, a3a)
- F-H: 更广泛的免疫学/肿瘤学查询

---

## 四、当前进展

### 4.1 库统计

| 指标 | 值 |
|------|-----|
| 总条目 | **250** |
| 种子数据 (手工策展) | 51 |
| PubMed 自动发现 | ~199 |
| 唯一基因 | 79 |
| 唯一 HLA alleles | 15 |
| 受累器官 | 33 |
| Mass-spec 确认 | 62 (24.8%) |

### 4.2 证据等级分布

| 等级 | 数量 | 占比 |
|------|------|------|
| Level 1: Clinical Fatal | 5 | 2.0% |
| Level 2: In Vitro Confirmed | 61 | 24.4% |
| Level 3: HT Screened | 27 | 10.8% |
| Level 4: In Silico | 12 | 4.8% |
| Unknown (待完善) | 145 | 58.0% |

### 4.3 验证状态分布

| 状态 | 数量 | 占比 |
|------|------|------|
| VALIDATED | 63 | 25.2% |
| PARTIAL | 30 | 12.0% |
| NEEDS_REVIEW | 12 | 4.8% |
| Unknown (未验证) | 145 | 58.0% |

### 4.4 Top HLA alleles

| HLA | 数量 |
|-----|------|
| HLA-A*02:01 | 70 |
| HLA-A*01:01 | 5 |
| HLA-A*24:02 | 4 |
| HLA-DQ*03:02 | 4 |
| HLA-B*27:05 | 3 |

### 4.5 Top 受累器官

| 器官 | 数量 | 临床严重性 |
|------|------|-----------|
| Skin | 12 | 中等 |
| Eye | 11 | 高 |
| Liver | 10 | 高 (重要脏器) |
| Pancreas | 8 | 高 |
| Brain | 5 | 极高 (重要脏器) |
| Joints | 5 | 中等 |
| Lung | 4 | 高 (重要脏器) |
| Heart | 3 | 极高 (重要脏器) |

### 4.6 Top 实验方法

| 方法 | 出现次数 |
|------|---------|
| qualitative binding | 187 |
| IFNg release | 158 |
| ligand presentation (mass-spec) | 87 |
| cytotoxicity | 79 |
| TNFa release | 58 |
| dissociation constant KD (SPR) | 54 |
| IL-2 release | 54 |

---

## 五、标志性案例: DC-0001 (MAGE-A3 / Titin)

**这是 TCR-T 安全性领域最重要的案例**：

| 字段 | 值 |
|------|-----|
| Decoy ID | DC-0001 |
| 靶标肽段 | EVDPIGHLY (MAGE-A3) |
| 脱靶肽段 | ESDPIVAQY (Titin) |
| HLA | HLA-A*01:01 |
| TCR | MAGE-A3 a3a TCR (affinity-enhanced) |
| 证据等级 | Level 1: Clinical Fatal |
| 受累器官 | Heart, Skeletal Muscle |
| Mass-spec | ✅ 确认 (HLA elution) |
| 表达模式 | Essential Organ Specific (心肌高表达) |

**事件经过**: 2013 年，Adaptimmune 的亲和力增强型 MAGE-A3 TCR 在临床试验中导致 2 名患者心脏骤停死亡。后续研究（Cameron et al., Sci Transl Med, 2013）发现该 TCR 交叉反应到心肌蛋白 Titin 的 ESDPIVAQY 肽段上。

**文献溯源**: 10 篇 PMID (23926201, 23863783, 26758806, ...) + 2 个临床试验 (NCT01273181, NCT01350401)

---

## 六、可视化

已完成 Decoy C 全套可视化（`figures/` 目录）：

| 文件 | 内容 |
|------|------|
| `decoy_c_pipeline.png` | 4-stage 管线流程图 (Discovery → Extraction → Validation → Curation) |
| `decoy_c_overview.png` | 四面板: 证据等级、验证状态、HLA 分布、受累器官 |
| `decoy_c_filtering.png` | 筛选漏斗图 (~200 papers → 250 entries → 63 validated → 5 fatal) |
| `decoy_c_assay_year.png` | 实验方法频次 + 发表年份时间线 |
| `decoy_c_interactive.html` | Plotly 交互式: DC-0001 Case Study + MAGE-A3 网络 + 旭日图 |

---

## 七、数据模型

### DecoyEntry 核心结构

```
DecoyEntry
├── decoy_id: "DC-NNNN"
├── peptide_info
│   ├── decoy_sequence: "ESDPIVAQY" (8-15 AA)
│   ├── hla_allele: "HLA-A*01:01"
│   ├── source_protein: "Titin"
│   ├── gene_symbol: "TTN"
│   └── uniprot_id: "Q8WZ42"
├── discovery_context
│   ├── original_target_sequence: "EVDPIGHLY"
│   ├── original_target_protein: "MAGE-A3"
│   └── tcr_name_or_id: "MAGE-A3 a3a TCR"
├── risk_profile
│   ├── evidence_level: Level_1_Clinical_Fatal
│   ├── critical_organs_affected: ["Heart", "Skeletal Muscle"]
│   └── expression_pattern: Essential_Organ_Specific
├── experimental_evidence
│   ├── mass_spec_confirmed: true
│   ├── assays_performed: ["Clinical Trial", "Cytotoxicity", ...]
│   └── cross_reactivity_affinity: "High"
├── provenance
│   ├── pmid: ["23926201", ...]
│   ├── clinical_trial_id: ["NCT01273181", ...]
│   └── evidence_summary: "Affinity-enhanced MAGE-A3 TCR caused fatal cardiac shock..."
├── source (论文元数据)
├── thought_process (LLM 推理过程, 反幻觉)
└── validation_flags
    ├── uniprot_match: "OK"
    ├── iedb_match: "FOUND"
    └── overall_status: "VALIDATED"
```

---

## 八、已知限制

| 限制 | 影响 | 解决方案 |
|------|------|---------|
| 145 条 Unknown 证据等级 | 大量条目缺乏证据分级 | 重新运行 LLM 提取或人工标注 |
| LLM 可能产生幻觉 | 少数提取结果不准确 | UniProt/IEDB 双重验证 + thought_process 审计 |
| IEDB 覆盖率有限 | 部分已知肽段不在 IEDB 中 | 标记为 PARTIAL, 不丢弃 |
| Rate limiting | PubMed/UniProt/IEDB 调用慢 | 批量处理 + 结果缓存 |
| HLA 格式不统一 | 文献中 HLA 命名混乱 | 标准化后处理 (→ HLA-X*NN:NN) |

---

## 九、文件清单

### 代码文件

| 文件 | 作用 |
|------|------|
| `decoy_c/fetcher.py` | PubMed 文献获取 (NCBI E-utilities) |
| `decoy_c/extractor.py` | LLM 结构化提取 (GPT-4o) |
| `decoy_c/validator.py` | UniProt + IEDB 双重验证 |
| `decoy_c/orchestrator.py` | 管线编排 + 库管理 + 去重 |
| `decoy_c/models.py` | Pydantic 数据模型 (DecoyEntry, DecoyLibrary) |
| `decoy_c/config.py` | API endpoints, 延迟配置, 路径 |
| `decoy_c/scale_up.py` | 自动化扩库 (多策略查询 + 断点续传) |
| `decoy_c/seed_data.py` | 种子数据加载器 |
| `scripts/visualize_decoy_c.py` | 全套可视化脚本 (matplotlib + plotly) |

### 数据文件

| 文件 | 大小 | 状态 |
|------|------|------|
| `data/decoy_c/decoy_library.json` | ~459 KB | 250 条 (当前库) |
| `data/decoy_c/seed_entries.json` | ~72 KB | 51 条种子数据 |
| `data/decoy_c/decoy_library_full_backup.json` | ~763 KB | 备份 |
| `data/decoy_c/removed_entries.json` | ~42 KB | 验证移除的条目 |
| `data/decoy_c/cross_reactivity_reference.json` | ~264 KB | 交叉反应参考 |
| `data/decoy_c/scale_up_checkpoint.json` | ~83 KB | 扩库断点 |

---

## 十、后续计划

1. **完善 Unknown 条目**: 对 145 条缺乏证据等级的条目重新提取或人工标注
2. **扩库至 1,000+**: 使用 scale_up.py 多策略查询扩展
3. **与 Decoy A/B 交叉验证**: 检查 A/B 的计算预测是否命中 C 库中的已知案例
4. **IEDB 批量注入** (`iedb_inject.py`): 直接从 IEDB 数据库批量导入已知表位
5. **定期更新**: 建立 PubMed alert 机制，自动追踪新发表的交叉反应论文
