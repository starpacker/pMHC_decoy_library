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
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --mpnn
```

## 文件清单

| 文件 | 作用 |
|------|------|
| `scanner.py` | MPNN 逆向设计管线主逻辑 |
| `orchestrator.py` | 与 Decoy B 主流程集成 |
| `tools/` | ProteinMPNN 封装 |
