# AlphaFold3 Linux 安装教程

**适用系统**: Ubuntu 22.04 LTS (推荐) / CentOS 8+ / Rocky Linux 9  
**最后更新**: 2026-04-03  
**AF3 版本**: v3.0.1

---

## 0. 硬件要求

| 项目 | 最低要求 | 推荐配置 |
|------|---------|---------|
| GPU | NVIDIA Ampere+ (CC 8.0+) | A100 80GB / H100 |
| GPU 显存 | 24 GB (RTX 3090/4090) | 80 GB (A100) |
| 系统内存 | 64 GB | 128 GB |
| 磁盘 | 1 TB SSD | 2 TB NVMe SSD |
| CUDA 驱动 | 525+ | 565+ |

> RTX 3090/4090 (24GB) 可以跑小分子 (<2000 tokens)，但 pMHC 三聚体约 800-1200 tokens，24GB 够用。

### 已验证 GPU 列表

| GPU | 显存 | 最大 tokens | 备注 |
|-----|------|------------|------|
| A100 80GB | 80 GB | ~5120 | 官方参考硬件 |
| A100 40GB | 40 GB | ~3000 | 需要 unified memory |
| RTX 4090 | 24 GB | ~2000 | 需要 unified memory + 内存溢出 |
| RTX 3090 | 24 GB | ~1500 | 同上 |
| V100 | 32 GB | ~2000 | 需要特殊 XLA flag (见下文) |

---

## 1. 系统依赖

```bash
# Ubuntu 22.04
sudo apt update && sudo apt install -y \
    git wget curl zstd build-essential cmake \
    gcc g++ zlib1g-dev libbz2-dev liblzma-dev

# 确认 NVIDIA 驱动
nvidia-smi
# 应该显示 Driver Version: 565.x+ 和 CUDA Version: 12.x
```

---

## 方法 A: Docker 安装（官方推荐）

适合生产环境、不想折腾依赖冲突的场景。

### A1. 安装 Docker

```bash
# 添加 Docker 官方源
sudo apt-get install -y ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg \
    -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo "deb [arch=$(dpkg --print-architecture) \
  signed-by=/etc/apt/keyrings/docker.asc] \
  https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io \
    docker-buildx-plugin docker-compose-plugin

# 免 sudo 运行 docker
sudo usermod -aG docker $USER
newgrp docker
```

### A2. 安装 NVIDIA Container Toolkit

```bash
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | \
    sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg

curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# 验证 GPU 透传
docker run --rm --gpus all nvidia/cuda:12.6.0-base-ubuntu22.04 nvidia-smi
```

### A3. 克隆并构建 AF3 镜像

```bash
mkdir -p ~/tools && cd ~/tools
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3

# 构建 Docker 镜像（约 15-30 分钟）
docker build -t alphafold3 -f docker/Dockerfile .

# 如果遇到 "No file descriptors" 错误 (RHEL/Rocky):
# docker build --ulimit nofile=65535:65535 -t alphafold3 -f docker/Dockerfile .
```

### A4. 下载序列数据库（~630 GB）

```bash
# 重要：放在 repo 目录之外！
mkdir -p ~/af3_data/databases
cd ~/tools/alphafold3
bash fetch_databases.sh ~/af3_data/databases

# 下载完成后设置权限
sudo chmod 755 -R ~/af3_data/databases
```

| 数据库 | 解压后大小 |
|--------|-----------|
| BFD (small) | ~272 GB |
| MGnify | ~120 GB |
| UniProt | ~99 GB |
| PDB (mmCIF) | ~78 GB |
| UniRef90 | ~59 GB |
| 其他 (NT, RNACentral 等) | ~2 GB |
| **合计** | **~630 GB** |

### A5. 放置模型权重

1. 申请权重: https://forms.gle/svvpY4u2jsHEwWYS6 （通常 2-3 个工作日）
2. 下载 `af3.bin.zst` 并解压:

```bash
mkdir -p ~/af3_data/models
cd ~/af3_data/models
# 将下载的 af3.bin.zst 放到此目录
unzstd af3.bin.zst
```

### A6. 运行预测

```bash
# 创建输入/输出目录
mkdir -p ~/af_input ~/af_output

# 准备输入 JSON（见下方"输入格式"章节）
# ...

# 运行
docker run -it --rm \
    --volume ~/af_input:/root/af_input \
    --volume ~/af_output:/root/af_output \
    --volume ~/af3_data/models:/root/models \
    --volume ~/af3_data/databases:/root/public_databases \
    --gpus all \
    alphafold3 \
    python run_alphafold.py \
        --json_path=/root/af_input/fold_input.json \
        --model_dir=/root/models \
        --output_dir=/root/af_output
```

---

## 方法 B: Conda 安装（无 Docker）

适合 HPC 集群、需要与 Python pipeline 深度集成的场景。**非官方支持，但社区验证可行。**

### B1. 创建 Conda 环境

```bash
conda create -n af3 python=3.11 -y
conda activate af3

# 防止 pip 安装到用户目录
conda env config vars set PYTHONUSERBASE=intentionally-disabled
conda deactivate && conda activate af3
```

### B2. 安装系统工具（通过 Conda）

```bash
conda install -c conda-forge cmake gcc gxx boost numpy bzip2 zstd git zlib gsl -y
conda install -c bioconda hmmer -y

# 验证 hmmer
jackhmmer -h | head -3
```

### B3. 安装 Python 依赖

```bash
# 核心依赖
pip install \
    pandas==2.2.3 \
    absl-py==2.1.0 \
    chex==0.1.87 \
    dm-haiku==0.0.14 \
    dm-tree==0.1.8 \
    "jax[cuda12]==0.4.34" \
    jaxlib==0.4.34 \
    jaxtyping==0.2.34 \
    jmp==0.0.4 \
    ml-dtypes==0.5.0 \
    numpy==2.1.3 \
    rdkit \
    tabulate==0.9.0

# CUDA 运行时（如果系统没有全局安装 CUDA 12）
pip install \
    nvidia-cublas-cu12 \
    nvidia-cuda-cupti-cu12 \
    nvidia-cuda-nvcc-cu12 \
    nvidia-cuda-runtime-cu12 \
    nvidia-cudnn-cu12 \
    nvidia-cufft-cu12 \
    nvidia-cusolver-cu12 \
    nvidia-cusparse-cu12 \
    nvidia-nccl-cu12 \
    nvidia-nvjitlink-cu12
```

### B4. 安装 AlphaFold3

```bash
cd ~/tools
git clone https://github.com/google-deepmind/alphafold3.git
cd alphafold3
pip install --no-deps .

# 构建数据组件
cd ${CONDA_PREFIX}/bin
./build_data
```

### B5. 验证安装

```bash
python -c "
import jax
print('JAX version:', jax.__version__)
print('GPU devices:', jax.devices())
print('GPU count:', jax.device_count())
"

python -c "from alphafold3 import structure; print('AF3 import OK')"
```

### B6. 下载数据库 + 权重

同方法 A 的 A4 和 A5 步骤。

### B7. 运行预测

```bash
cd ~/tools/alphafold3

python run_alphafold.py \
    --json_path=~/af_input/fold_input.json \
    --model_dir=~/af3_data/models \
    --db_dir=~/af3_data/databases \
    --output_dir=~/af_output
```

---

## 3. pMHC 输入 JSON 格式

这是 Decoy B pipeline 使用的输入格式（对应 `decoy_b/tools/alphafold3.py` 中的 `_build_input_json`）：

```json
{
    "name": "pmhc_GILGFVFTL_HLAA0201",
    "modelSeeds": [1],
    "sequences": [
        {
            "protein": {
                "id": "A",
                "sequence": "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRT"
            }
        },
        {
            "protein": {
                "id": "B",
                "sequence": "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
            }
        },
        {
            "protein": {
                "id": "C",
                "sequence": "GILGFVFTL"
            }
        }
    ],
    "dialect": "alphafold3",
    "version": 2
}
```

- Chain A = HLA 重链 (MHC heavy chain)
- Chain B = β2-microglobulin
- Chain C = Peptide

---

## 4. GPU 显存不够时的优化

```bash
# 方法 1: 开启 unified memory（让 GPU 内存溢出到系统 RAM）
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_CLIENT_MEM_FRACTION=4.0

# 方法 2: 分离 data pipeline 和 inference（批量处理时有效）
# 先跑序列搜索（CPU 密集，不需要 GPU）
python run_alphafold.py \
    --json_path=input.json \
    --db_dir=~/af3_data/databases \
    --output_dir=~/af_output \
    --norun_inference

# 再跑结构推理（GPU 密集）
python run_alphafold.py \
    --json_path=input.json \
    --model_dir=~/af3_data/models \
    --output_dir=~/af_output \
    --norun_data_pipeline

# 方法 3: 启用 JAX 编译缓存（避免重复编译）
python run_alphafold.py \
    --json_path=input.json \
    --jax_compilation_cache_dir=~/af3_cache \
    ...
```

---

## 5. V100 / 旧 GPU 特殊处理

V100 (Compute Capability 7.0) 需要禁用自定义 kernel：

```bash
export XLA_FLAGS="--xla_disable_hlo_passes=custom-kernel-fusion-rewriter"
```

---

## 6. 与 Decoy B Pipeline 集成

安装完成后，配置环境变量让 `decoy_b/tools/alphafold3.py` 自动发现：

```bash
# 添加到 ~/.bashrc 或 .env 文件
export AF3_DIR=~/tools/alphafold3
export AF3_MODEL_DIR=~/af3_data/models
export AF3_DB_DIR=~/af3_data/databases

# Docker 模式
export AF3_USE_DOCKER=true
export AF3_DOCKER_IMAGE=alphafold3

# 或本地 Python 模式
export AF3_USE_DOCKER=false
```

验证集成：

```bash
cd /path/to/decoy_library
python -c "from decoy_b.tools.alphafold3 import check_available; print('AF3 ready:', check_available())"
```

运行完整 Decoy B pipeline：

```bash
PYTHONUTF8=1 python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01
```

---

## 7. 常见问题

| 症状 | 原因 | 解决 |
|------|------|------|
| `Allocator (GPU_0_bfc) ran out of memory` | GPU 显存不足 | 设置 `TF_FORCE_UNIFIED_MEMORY=1` 等（见第4节） |
| `could not select device driver "nvidia"` | 缺少 nvidia-container-toolkit | 安装 NVIDIA Container Toolkit（见 A2） |
| `No file descriptors available` (Docker build) | RHEL 系文件描述符限制 | `docker build --ulimit nofile=65535:65535 ...` |
| Docker build 极慢 / 磁盘爆满 | 数据库放在了 repo 目录里 | 把 databases 和 models 放在 repo **外面** |
| `Segmentation fault` (非 Docker) | JAX import 顺序问题 | 确保 JAX 在其他 C++ 扩展前 import |
| `NVIDIA driver CUDA version < PTX compiler` | 驱动太旧 | 升级到 565+ 并重启 |
| hmmer 相关错误 | jackhmmer 未安装 | `conda install -c bioconda hmmer` 或 `apt install hmmer` |

---

## 8. 参考资源

- 官方仓库: https://github.com/google-deepmind/alphafold3
- 官方安装文档: https://github.com/google-deepmind/alphafold3/blob/main/docs/installation.md
- 社区 Conda 安装指南: https://github.com/Model3DBio/AlphaFold3-Conda-Install
- 权重申请: https://forms.gle/svvpY4u2jsHEwWYS6
- 权重使用条款: 仅限非商业用途
