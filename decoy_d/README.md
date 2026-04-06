# Decoy D — MPNN 逆向设计

**状态**: 代码框架完成，依赖 ProteinMPNN 权重部署

## 核心思想

使用 ProteinMPNN 从靶标 pMHC 结构出发，固定 HLA 锚定位 (p2, p9)，逆向设计保持 TCR 接触面的新肽段序列。发现非同源但结构功能等价的候选——覆盖计算同源搜索无法触及的"暗物质"空间。

## 管线架构

```
Input: 靶标 pMHC 结构 (tFold/AF3 预测)

1. 固定位点标记
   HLA 链 (M, N) 全部固定 + 锚定位 (p2, p9) 固定
   仅 TCR 接触位 (p4-p8) 允许设计

2. ProteinMPNN 序列生成
   多温度采样 (T=0.1-0.5) → ~1,000 设计序列

3. 双层过滤
   Tier 1: 人类蛋白组精确匹配 → 最高可信度
   Tier 2: mhcflurry HLA 呈递预测 → 理论风险候选

4. 评分与排名
   与 Decoy B 主流程合并 → risk_scorer 统一排名
```

## 使用方法

```bash
# 统一入口（推荐）
python run_decoy.py GILGFVFTL d
python run_decoy.py GILGFVFTL d --hla HLA-A*02:01 --designs 2000

# 模块入口
python -m decoy_d --target GILGFVFTL --hla "HLA-A*02:01"

# 集成到 Decoy B 管线
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --mpnn

# 批量运行全部候选靶标 + 生成 Pareto 图表
bash scripts/batch_decoy_d.sh
```

## 文件清单

| 文件 | 作用 |
|------|------|
| `scanner.py` | MPNN 逆向设计管线主逻辑 |
| `__main__.py` | `python -m decoy_d` 入口 |
| `main.py` | CLI argparse 定义 |
| `orchestrator.py` | 与 Decoy B 主流程集成 |
| `tools/` | ProteinMPNN / tFold 封装 |
| `../run_decoy.py` | 项目根目录统一 CLI |
| `../scripts/batch_decoy_d.sh` | 批量运行脚本 |
| `../scripts/visualize_decoy_d_detail.py` | Decoy D 3-panel 可视化 |

## 常见问题 (Q&A)

**Q: Decoy D 的 MPNN 逆向设计与 Decoy A/B 有什么本质区别？**
A: Decoy A 和 B 都是在**已知**的人类蛋白质组中进行搜索（自上而下）。而 Decoy D 是自下而上的生成式方法：它固定了 HLA 结合所需的锚定位，让 AI 自由生成能形成相似 TCR 接触面的新序列。这不仅能发现人类基因组中罕见的变异或未被充分注释的蛋白，还能预测未来可能出现的病毒突变带来的交叉反应风险。

**Q: 为什么在 MPNN 设计后还需要 mhcflurry 过滤？**
A: ProteinMPNN 主要关注三维结构的合理性和序列的折叠能力，但它并不直接优化多肽在细胞内的抗原加工和 HLA 呈递效率。生成的序列即使在结构上能完美契合 TCR，如果不能被细胞有效呈递到表面，也不会引发实际的毒性。因此，必须使用 mhcflurry 进行二次过滤，确保候选多肽具有真实的生物学呈递潜力。
