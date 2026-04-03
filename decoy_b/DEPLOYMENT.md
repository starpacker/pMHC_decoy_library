# Decoy B 部署指南

**状态**: 代码已就绪，等待模型权重部署
**最后更新**: 2026-04-03

---

## 概览

Decoy B 依赖四个外部工具进行结构预测、交叉验证和序列设计：

| 工具 | 用途 | 阶段 | 必需? |
|------|------|------|-------|
| **tFold** | pMHC 结构预测（批量筛选） | Stage 2 | 推荐（Stage 2 核心） |
| **AlphaFold3** | pMHC 高精度结构预测（精修） | Stage 3 | 可选（提升精度） |
| **Boltz-2** | pMHC 结构交叉验证 | Stage 4 | 推荐（提升可靠性） |
| **ProteinMPNN** | 逆向序列设计 | MPNN 分支 | 可选（发现新候选） |

Pipeline 会自动检测可用工具并跳过不可用的阶段。即使只有 tFold，也能运行完整的 Decoy B 筛选。

---

## 一键安装

```bash
# 安装全部工具（克隆 repo + 创建目录 + 安装依赖）
python scripts/setup_decoy_b.py --all

# 如果在中国大陆，使用 GitHub 代理加速
python scripts/setup_decoy_b.py --all --proxy https://ghfast.top/

# 只安装某个工具
python scripts/setup_decoy_b.py --tfold
python scripts/setup_decoy_b.py --af3
python scripts/setup_decoy_b.py --mpnn

# 指定安装目录
python scripts/setup_decoy_b.py --all --tools-dir D:/bio_tools
```

---

## 1. tFold 部署

### 1.1 克隆代码

```bash
git clone https://github.com/TencentAI4S/tfold.git ~/tools/tfold
cd ~/tools/tfold
pip install -e .
```

### 1.2 放置权重

将 tFold 模型权重文件放到以下路径：

```
~/tools/tfold/weights/
├── esm_ppi_650m/          # ESM-PPI 预训练模型
│   └── *.pt               # PyTorch 权重文件
└── tfold_tcr_pmhc/        # tFold TCR-pMHC 主干模型
    └── *.pt               # PyTorch 权重文件
```

**权重来源**:
- DP Tech 内部模型仓库
- Zenodo (公开版本)
- 微云 / Google Drive (根据 tFold README)

### 1.3 环境变量

```bash
export TFOLD_DIR=~/tools/tfold
export TFOLD_WEIGHTS_DIR=~/tools/tfold/weights
```

### 1.4 验证

```bash
python -c "from decoy_b.tools.tfold import check_available; print(check_available())"
```

### 1.5 可选：推理服务模式

如果已部署 tFold 推理服务（HTTP API），可以设置：

```bash
export TFOLD_SERVER_URL=http://localhost:8501
```

Wrapper 会自动尝试 Python API -> CLI -> Server 三种模式。

---

## 2. AlphaFold3 部署

### 2.1 克隆代码

```bash
git clone https://github.com/google-deepmind/alphafold3.git ~/tools/alphafold3
```

### 2.2 放置权重

**申请**: https://forms.gle/svvpY4u2jsHEwWYS6

将权重放到：

```
~/tools/alphafold3/models/
├── af3.bin.zst            # 或其他格式的模型权重
└── ...
```

### 2.3 选择运行模式

**模式 A: Docker（推荐用于生产环境）**

```bash
cd ~/tools/alphafold3
docker build -t alphafold3 -f docker/Dockerfile .

# 下载序列数据库
bash fetch_databases.sh ~/tools/alphafold3/databases
```

环境变量：
```bash
export AF3_USE_DOCKER=true
export AF3_DOCKER_IMAGE=alphafold3
export AF3_MODEL_DIR=~/tools/alphafold3/models
export AF3_DB_DIR=~/tools/alphafold3/databases
```

**模式 B: 本地 Python**

```bash
cd ~/tools/alphafold3
pip install -e .
```

环境变量：
```bash
export AF3_USE_DOCKER=false
export AF3_DIR=~/tools/alphafold3
export AF3_MODEL_DIR=~/tools/alphafold3/models
```

### 2.4 验证

```bash
python -c "from decoy_b.tools.alphafold3 import check_available; print(check_available())"
```

### 2.5 注意事项

- AF3 每个预测需要 2-5 分钟（GPU），用于 Top 200 候选的精修阶段
- 如果没有 AF3，Pipeline 会跳过 Stage 3 refinement，只使用 tFold 结构
- Docker 模式需要 NVIDIA GPU + nvidia-docker

---

## 3. ProteinMPNN 部署

### 3.1 克隆代码

```bash
git clone https://github.com/dauparas/ProteinMPNN.git ~/tools/ProteinMPNN
```

**注意**: ProteinMPNN 的权重已包含在 repo 中（`vanilla_model_weights/` 目录），无需单独下载。

### 3.2 依赖

```bash
pip install torch  # ProteinMPNN 需要 PyTorch
```

### 3.3 目录结构

确认以下文件存在：

```
~/tools/ProteinMPNN/
├── protein_mpnn_run.py              # 主运行脚本
├── helper_scripts/
│   ├── parse_multiple_chains.py     # PDB 解析
│   ├── make_fixed_positions_dict.py # 固定位置生成
│   └── assign_fixed_chains.py       # 链指定
└── vanilla_model_weights/
    ├── v_48_002.pt                  # 模型权重文件
    ├── v_48_010.pt
    └── v_48_020.pt                  # 默认使用
```

### 3.4 环境变量

```bash
export PROTEINMPNN_DIR=~/tools/ProteinMPNN
```

### 3.5 验证

```bash
python -c "from decoy_b.tools.proteinmpnn import check_available; print(check_available())"
```

---

## 4. Boltz-2 部署

### 4.1 克隆代码

```bash
git clone <boltz-repo> ~/tools/boltz
cd ~/tools/boltz
pip install -e .
```

**注意**: Boltz-2 模型权重会在首次运行时自动下载到 `~/.boltz/`（约 2-3 GB）。

### 4.2 验证安装

```bash
# 检查 CLI
boltz --help

# 检查 Python API
python -c "from decoy_b.tools.boltz import check_available; print(check_available())"

# Quick smoke test
boltz predict ~/tools/boltz/examples/prot.yaml --out_dir /tmp/boltz_test --devices 1
```

### 4.3 环境变量

```bash
export BOLTZ_DIR=~/tools/boltz          # Boltz 代码目录
export BOLTZ_CACHE=~/.boltz             # 权重和缓存目录
export BOLTZ_MODEL=boltz2               # 模型版本 (boltz1 或 boltz2)
export BOLTZ_DEVICE=gpu                 # 设备 (gpu 或 cpu)
```

### 4.4 注意事项

- Boltz-2 首次运行需下载 CCD 字典 + 模型权重（~2-3 GB）
- 预测速度与 AF3 相当（每个 pMHC 复合体 2-5 分钟，GPU）
- 用于 Top 200 候选的交叉验证阶段
- 如果没有 Boltz，Pipeline 会跳过 Stage 4 交叉验证，不影响其他阶段
- MIT 开源许可，可用于商业用途

---

## 5. 完整验证

运行验证脚本检查所有工具：

```bash
python scripts/verify_decoy_b.py

# 详细输出
python scripts/verify_decoy_b.py --verbose

# 包含 smoke test（需要已安装工具）
python scripts/verify_decoy_b.py --smoke-test
```

预期输出示例：

```
== tFold ==
  [+] TFOLD_DIR exists: OK
  [+] Weights directory: OK
  [+] tFold available (any mode): OK

== AlphaFold3 ==
  [+] AF3_DIR exists: OK
  [+] Model weights present: OK
  [+] AF3 available: OK

== ProteinMPNN ==
  [+] PROTEINMPNN_DIR exists: OK
  [+] Run script: OK
  [+] Weight files (.pt): 3 found: OK
  [+] ProteinMPNN available: OK

== Boltz-2 ==
  [+] BOLTZ_DIR exists: OK
  [+] Boltz CLI available: OK
  [+] Boltz available (any mode): OK

Summary:
  tFold                : READY
  AlphaFold3           : READY
  Boltz-2              : READY
  ProteinMPNN          : READY
  HLA backend          : READY
```

---

## 6. 运行 Decoy B

### 6.1 物理化学筛选（不需要外部工具）

Stage 1 使用 Atchley Factor 余弦相似度，纯 NumPy 计算，不需要任何外部工具：

```bash
PYTHONUTF8=1 python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01 --skip-structural
```

### 6.2 完整 Pipeline（含结构预测）

```bash
PYTHONUTF8=1 python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01
```

### 6.3 完整 A+B Pipeline

```bash
PYTHONUTF8=1 python -m decoy_b run --target GILGFVFTL --hla HLA-A*02:01
```

### 6.4 Python API

```python
import os; os.environ["PYTHONUTF8"] = "1"

from decoy_b.scanner import scan_decoy_b

# 仅物理化学筛选
hits = scan_decoy_b("GILGFVFTL", "HLA-A*02:01", run_structural=False)

# 完整 5 阶段 (含 Boltz 交叉验证)
hits = scan_decoy_b("GILGFVFTL", "HLA-A*02:01", run_structural=True)

# 不含 Boltz 交叉验证
hits = scan_decoy_b("GILGFVFTL", "HLA-A*02:01", run_boltz_crossval=False)

# 含 MPNN 逆向设计分支
hits = scan_decoy_b("GILGFVFTL", "HLA-A*02:01", run_mpnn=True)

for h in hits[:5]:
    print(f"{h.sequence}  cos={h.physicochemical.cosine_similarity:.3f}  "
          f"HD={h.hamming_distance}  genes={h.gene_symbols}")
```

---

## 7. 环境变量汇总

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `TFOLD_DIR` | `~/tools/tfold` | tFold 代码目录 |
| `TFOLD_WEIGHTS_DIR` | `$TFOLD_DIR/weights` | tFold 权重目录 |
| `TFOLD_SERVER_URL` | (空) | tFold 推理服务 URL |
| `AF3_DIR` | `~/tools/alphafold3` | AF3 代码目录 |
| `AF3_MODEL_DIR` | `$AF3_DIR/models` | AF3 模型权重目录 |
| `AF3_DB_DIR` | `$AF3_DIR/databases` | AF3 序列数据库目录 |
| `AF3_USE_DOCKER` | `true` | 使用 Docker 模式 |
| `AF3_DOCKER_IMAGE` | `alphafold3` | Docker 镜像名 |
| `AF3_GPU_DEVICE` | `all` | Docker GPU 设备 |
| `PROTEINMPNN_DIR` | `~/tools/ProteinMPNN` | MPNN 代码目录 |
| `BOLTZ_DIR` | `~/tools/boltz` | Boltz 代码目录 |
| `BOLTZ_CACHE` | `~/.boltz` | Boltz 权重和缓存目录 |
| `BOLTZ_MODEL` | `boltz2` | 模型版本 (boltz1/boltz2) |
| `BOLTZ_DEVICE` | `gpu` | 设备 (gpu/cpu) |
| `PYTHONUTF8` | `1` | Windows UTF-8 兼容 |

可将这些变量写入 `.env` 文件（参考 `.env.example`）。

---

## 8. 故障排除

### tFold 找不到权重

```
Error: tFold predict script not found under ~/tools/tfold
```

确认 `TFOLD_DIR` 指向正确路径，且 `projects/tfold_ag/predict.py` 或 `projects/tfold_tcr/predict.py` 存在。

### AF3 Docker GPU 错误

```
Error: could not select device driver "nvidia"
```

需要安装 nvidia-docker2。或改用本地模式 (`AF3_USE_DOCKER=false`)。

### ProteinMPNN CUDA 错误

如果没有 GPU，ProteinMPNN 会自动使用 CPU（较慢但可用）。确保 PyTorch 已正确安装：

```bash
python -c "import torch; print(torch.cuda.is_available())"
```

### Windows 编码问题

始终设置 `PYTHONUTF8=1`：

```bash
set PYTHONUTF8=1  # cmd
$env:PYTHONUTF8="1"  # PowerShell
export PYTHONUTF8=1  # bash/git bash
```
