# Decoy B Pipeline — 技术进展报告

**日期**: 2026-04-01
**状态**: 代码框架完成，物化筛选可用；结构建模依赖外部工具权重（tFold/AF3 待部署）

---

## 一、技术原理

### 1.1 为什么需要 Decoy B

Decoy A 通过 Hamming distance ≤ 2 找到序列近亲，但 TCR 的识别并非纯序列匹配——它识别的是 pMHC 复合物表面的 **三维构象和电荷分布**。这意味着：

- 两个序列差异很大的肽段（Hamming > 2），可能因为氨基酸理化性质相近，折叠出相似的 TCR 接触面
- 经典案例：MAGE-A3 靶标 EVDPIGHLY 与 Titin 的 ESDPIVAQY，Hamming distance = 5（远超 Decoy A 的阈值 2），但两者在 HLA-A*01:01 上呈递后 TCR 接触面高度相似，导致了临床致死事件

**Decoy B 的核心思想**：放宽序列相似性的要求，转而寻找 **物化性质相似** 和 **三维结构相似** 的脱靶候选。

### 1.2 物化性质相似性 — Atchley Factors

Atchley factors 是通过对氨基酸 494 项理化属性做因子分析（factor analysis）提取出的 5 个正交维度：

| 因子 | 含义 | 生物学意义 |
|------|------|-----------|
| Factor I | 极性 (polarity) | 影响溶剂可及性和氢键网络 |
| Factor II | 二级结构倾向 | α-helix / β-sheet / coil 偏好 |
| Factor III | 分子体积 | 侧链大小，影响空间位阻 |
| Factor IV | 密码子多样性 | 进化保守性的代理指标 |
| Factor V | 静电荷 | 影响 TCR CDR3 环的电荷互补 |

每个氨基酸被表示为一个 5 维向量。一条肽段的 TCR 接触面（9-mer 的 p4-p8，共 5 个残基）被展开为 25 维向量。

**核心度量**: 余弦相似度 (cosine similarity)

```
cos_sim(A, B) = (A · B) / (||A|| × ||B||)
```

阈值 ≥ 0.70：表示两条肽段在 TCR 接触面上理化性质高度相似。

### 1.3 三维结构相似性

物化筛选后，对候选集进行 pMHC 三维结构建模和比较：

| 步骤 | 工具 | 作用 |
|------|------|------|
| 快速建模 | tFold | TCR-pMHC 结构预测，GPU 上 ~10-30 秒/肽 |
| 精细化 | AlphaFold3 | 高精度 refinement，~2-5 分钟/肽 |
| 逆向设计 | ProteinMPNN | 固定 MHC+锚定位，设计保持 TCR 接触面的新序列 |

**结构比较指标**:
- **Backbone RMSD**: 肽段主链叠合后的均方根偏差（< 2.0 Å = 高度相似）
- **B-factor 相关**: 表面柔性分布的皮尔逊相关（反映动态构象相似性）

### 1.4 管线架构

```
Input: 靶标肽段 + HLA 型别 + HLA 过滤后的候选池 (来自 Decoy A Step 2)

Stage 1: 物化筛选 (CPU, 秒级)
  ├─ 提取 TCR 接触残基 (central 5 residues)
  ├─ 转换为 Atchley 向量 (25维)
  ├─ 批量余弦相似度 (NumPy 向量化)
  ├─ 排除 Hamming ≤ 2 (已由 Decoy A 覆盖)
  └─ 阈值: cosine_sim ≥ 0.70
  → 2,000-5,000 candidates

Stage 2: tFold 批量结构预测 (GPU, 分钟级)
  ├─ 对靶标 + Top-K 候选建模 pMHC 结构
  ├─ 输出: PDB + 置信度分数
  └─ 吞吐: ~25x faster than AF3
  → PDB structures

Stage 3: AF3 精细化 (GPU, 小时级, 可选)
  ├─ Top-200 候选 → AlphaFold3 精细建模
  ├─ 输出: pLDDT, pAE, ranking score
  └─ 适合需要高精度的最终验证
  → Refined structures

Stage 4: 结构比较 & 评分
  ├─ 肽段主链 RMSD (BioPython superposition)
  ├─ B-factor 相关 (表面柔性)
  └─ 组合评分: 40% 物化相似 + 60% 结构相似
  → DecoyBHit 列表

Optional: MPNN 逆向设计分支
  ├─ ProteinMPNN 固定 HLA 锚定位 (p2, p9)
  ├─ 设计 ~1,000 保持 TCR 面的新序列
  ├─ 映射回人类蛋白组
  └─ 补充发现非同源但结构相似的候选

Final: 与 Decoy A 结果合并 → 综合风险评分
```

---

## 二、实现方法

### 2.1 物化筛选引擎 (`scanner.py`)

**批量余弦相似度计算**:

```python
# O(n×25) 向量化, 10万条候选 <1秒
target_norm = target_vec / (np.linalg.norm(target_vec) + 1e-12)
norms = np.linalg.norm(candidate_vecs, axis=1, keepdims=True) + 1e-12
cand_normed = candidate_vecs / norms
similarities = cand_normed @ target_norm
```

**TCR 接触残基提取**:
- 9-mer: 取中间 5 个残基 (p4-p8, 0-indexed 3-7)
- 其他长度: 居中对齐取 min(5, n-2) 个残基

### 2.2 tFold 封装 (`tools/tfold.py`)

三种执行模式（按优先级回退）：

| 模式 | 方式 | 适用场景 |
|------|------|---------|
| Python API | `import tfold` 直接调用 | 安装了 tfold Python 包 |
| CLI | `python tfold/predict.py` | 只 clone 了仓库 |
| HTTP Server | `POST /predict` | 远程 GPU 服务 |

**输入**: 肽段序列 + HLA allele → **输出**: PDB 结构 + 置信度

### 2.3 AlphaFold3 封装 (`tools/alphafold3.py`)

两种执行模式：

| 模式 | 方式 | 适用场景 |
|------|------|---------|
| Docker | `docker run --gpus=all alphafold3` | 推荐（隔离环境） |
| Local Python | `python run_alphafold.py` | 本地 GPU + 手动安装依赖 |

**输入**: JSON 格式的 pMHC 输入规范 → **输出**: CIF/PDB + 质量指标

### 2.4 ProteinMPNN 封装 (`tools/proteinmpnn.py`)

**流程**:
1. 解析靶标 pMHC PDB 结构
2. 标记固定位点（HLA 链 + 锚定位 p2, p9）
3. ProteinMPNN 在不同温度下生成 ~1,000 条序列
4. 过滤：保留 HLA 结合 + 人类蛋白组匹配的序列

### 2.5 综合风险评分 (`risk_scorer.py`)

**公式**:
```
Risk = Similarity × (1 / EL_Rank) × TPM_Weight

其中:
  Similarity = 0.4 × cosine_sim + 0.6 × surface_correlation  (有结构时)
             = cosine_sim                                      (仅物化时)

  TPM_Weight = 10.0  (重要脏器 TPM > 10)
             = 1.0   (普通组织)
             = 0.1   (testis/placenta)
```

---

## 三、当前进展

### 3.1 已完成

| 模块 | 状态 | 说明 |
|------|------|------|
| Atchley 物化筛选引擎 | ✅ 完成 | 纯 NumPy, 无外部依赖 |
| BioPython 结构比较 | ✅ 完成 | RMSD + B-factor 相关 |
| Pydantic 数据模型 | ✅ 完成 | DecoyBHit, PhysicochemFeatures, StructuralScore |
| 综合风险评分 | ✅ 完成 | A+B 合并排名 |
| tFold wrapper | ✅ 完成 | 三模式回退 |
| AF3 wrapper | ✅ 完成 | Docker + Local 双模式 |
| ProteinMPNN wrapper | ✅ 完成 | 逆向设计全流程 |
| CLI 接口 | ✅ 完成 | `python -m decoy_b scan-b / run / show` |
| 断点续传 | ✅ 完成 | PipelineState checkpoint |
| 部署脚本 | ✅ 完成 | `scripts/setup_decoy_b.py` + `scripts/verify_decoy_b.py` |
| 部署文档 | ✅ 完成 | `decoy_b/DEPLOYMENT.md` |

### 3.2 待部署

| 依赖 | 状态 | 获取方式 |
|------|------|---------|
| tFold 模型权重 | ⏳ 待下载 | `github.com/TencentAI4S/tfold` |
| AF3 模型权重 | ⏳ 待申请 | Google DeepMind 学术申请 |
| ProteinMPNN 权重 | ✅ 已含 | 仓库自带 `vanilla_model_weights/` |
| GPU (NVIDIA) | ⏳ 需配置 | tFold/AF3 需要 GPU 加速 |

### 3.3 可立即使用的功能

即使没有外部工具，Decoy B 的 **Stage 1 物化筛选** 可以独立运行：

```bash
# 仅物化筛选 (纯 CPU, 无外部依赖)
python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01 --skip-structural
```

这会返回所有 cosine_sim ≥ 0.70 且 Hamming > 2 的候选，按物化相似度排序。

---

## 四、关键配置参数

| 参数 | 值 | 作用 |
|------|-----|------|
| `PHYSICOCHEMICAL_COSINE_THRESHOLD` | 0.70 | Atchley 余弦相似度阈值 |
| `PHYSICOCHEMICAL_TOP_K` | 5,000 | 进入结构建模的候选数 |
| `SURFACE_SIMILARITY_THRESHOLD` | 0.75 | 结构表面相似性阈值 |
| `ORGAN_TPM_WEIGHT_HIGH` | 10.0 | 重要脏器 TPM 权重 |
| `ORGAN_TPM_WEIGHT_LOW` | 0.1 | Testis/Placenta 权重 |
| `FINAL_TOP_N` | 100 | 最终输出的 Top-N |

---

## 五、与 Decoy A 的互补关系

| 维度 | Decoy A | Decoy B |
|------|---------|---------|
| 核心方法 | Hamming distance (序列匹配) | Atchley cosine + 结构叠合 |
| 搜索范围 | Hamming ≤ 2 (近亲) | Hamming > 2 (远亲) |
| 计算开销 | 秒级 (NumPy 向量化) | 分钟-小时级 (结构建模) |
| 发现范围 | 1-2 个氨基酸差异的近亲 | 序列远但结构近的"隐性风险" |
| 已知案例 | GILGFVFTL → TAP2, PKD1L1 | EVDPIGHLY → ESDPIVAQY (Titin) |
| 互补性 | 高灵敏度序列筛查 | 深层结构风险发现 |

**整合策略**: Decoy A + B 的结果合并后由 `risk_scorer.py` 统一排名，输出 Top-100。

---

## 六、文件清单

### 代码文件

| 文件 | 作用 |
|------|------|
| `decoy_b/scanner.py` | 物化筛选 + 结构比较 + MPNN 分支 |
| `decoy_b/orchestrator.py` | 全管线编排 (A+B) + 断点续传 |
| `decoy_b/risk_scorer.py` | 综合风险评分 (A+B 合并) |
| `decoy_b/main.py` | CLI 入口 |
| `decoy_b/tools/tfold.py` | tFold pMHC 预测封装 |
| `decoy_b/tools/alphafold3.py` | AF3 高精度 refinement 封装 |
| `decoy_b/tools/proteinmpnn.py` | ProteinMPNN 逆向设计封装 |
| `scripts/setup_decoy_b.py` | 自动化部署脚本 |
| `scripts/verify_decoy_b.py` | 工具可用性验证 |

### 数据文件

| 文件 | 位置 | 状态 |
|------|------|------|
| HLA 过滤候选池 | `data/decoy_a/hla_filtered_*.parquet` | 共享自 Decoy A |
| 组织表达谱 | `data/decoy_a/gene_expression.parquet` | 共享自 Decoy A |
| Decoy B 结果 | `data/decoy_b/decoy_b_results.json` | 待生成 |
| A+B 合并排名 | `data/decoy_b/final_ranked_decoys.json` | 待生成 |
| pMHC 结构缓存 | `data/decoy_b/pmhc_models/` | 待生成 |

---

## 七、后续计划

1. **部署 tFold**: 下载模型权重，配置 GPU 环境，运行验证
2. **申请 AF3 权重**: 提交 Google DeepMind 学术申请表
3. **端到端测试**: 对 MAGE-A3 / Titin case 运行完整 4-stage 管线
4. **可视化**: 物化相似度分布图、结构叠合可视化、A+B 合并排名散点图
