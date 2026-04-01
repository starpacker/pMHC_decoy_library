#!/usr/bin/env python3
"""
Decoy B Tool Setup Script
==========================
Clones repositories, creates directory structure, installs Python
dependencies, and prepares weight directories for:
  - tFold (TCR-pMHC structure prediction)
  - AlphaFold3 (high-accuracy structure prediction)
  - ProteinMPNN (inverse sequence design)

Usage:
    python scripts/setup_decoy_b.py --all
    python scripts/setup_decoy_b.py --tfold
    python scripts/setup_decoy_b.py --af3
    python scripts/setup_decoy_b.py --mpnn
    python scripts/setup_decoy_b.py --deps-only   # Python deps only

After running this script, place model weights at the paths printed
at the end. See decoy_b/DEPLOYMENT.md for detailed instructions.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

# ── Default Paths ─────────────────────────────────────────────────────────

TOOLS_DIR = Path(os.getenv("DECOY_TOOLS_DIR", Path.home() / "tools"))

TFOLD_DIR = TOOLS_DIR / "tfold"
AF3_DIR = TOOLS_DIR / "alphafold3"
MPNN_DIR = TOOLS_DIR / "ProteinMPNN"

# Git repos
TFOLD_REPO = "https://github.com/TencentAI4S/tfold.git"
AF3_REPO = "https://github.com/google-deepmind/alphafold3.git"
MPNN_REPO = "https://github.com/dauparas/ProteinMPNN.git"

# GitHub proxy for China mainland (comment out if not needed)
GITHUB_PROXY = os.getenv("GITHUB_PROXY", "")  # e.g., "https://ghfast.top/"

IS_WINDOWS = platform.system() == "Windows"


def _proxy_url(url: str) -> str:
    """Apply GitHub proxy if set."""
    if GITHUB_PROXY and "github.com" in url:
        return url.replace("https://github.com/", GITHUB_PROXY)
    return url


def run_cmd(cmd: list[str], cwd: str | None = None, check: bool = True) -> int:
    """Run a command with live output."""
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd)
    if check and result.returncode != 0:
        print(f"  [WARN] Command exited with code {result.returncode}")
    return result.returncode


def check_git() -> bool:
    """Check git is available."""
    try:
        subprocess.run(["git", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


def check_conda() -> bool:
    """Check conda is available."""
    try:
        subprocess.run(["conda", "--version"], capture_output=True, check=True)
        return True
    except (FileNotFoundError, subprocess.CalledProcessError):
        return False


# ── tFold Setup ───────────────────────────────────────────────────────────

def setup_tfold(force: bool = False) -> None:
    """Clone tFold repo and prepare weight directories."""
    print("\n" + "=" * 60)
    print("  Setting up tFold")
    print("=" * 60)

    if TFOLD_DIR.exists() and not force:
        print(f"  tFold already exists at {TFOLD_DIR}")
        print("  Use --force to re-clone")
    else:
        if TFOLD_DIR.exists():
            shutil.rmtree(TFOLD_DIR)
        TOOLS_DIR.mkdir(parents=True, exist_ok=True)
        run_cmd(["git", "clone", "--depth", "1", _proxy_url(TFOLD_REPO), str(TFOLD_DIR)])

    # Create weight directories
    weights_dir = TFOLD_DIR / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)

    # tFold expects these subdirectories
    for subdir in ["esm_ppi_650m", "tfold_tcr_pmhc"]:
        (weights_dir / subdir).mkdir(exist_ok=True)

    print(f"\n  tFold directory: {TFOLD_DIR}")
    print(f"  Weights directory: {weights_dir}")
    print("  Place weights:")
    print(f"    - ESM-PPI model   -> {weights_dir / 'esm_ppi_650m/'}")
    print(f"    - tFold TCR model -> {weights_dir / 'tfold_tcr_pmhc/'}")

    # Install tFold Python dependencies
    tfold_req = TFOLD_DIR / "requirements.txt"
    if tfold_req.exists():
        print("\n  Installing tFold Python dependencies...")
        run_cmd([sys.executable, "-m", "pip", "install", "-r", str(tfold_req)], check=False)

    # Try installing tfold as a package
    tfold_setup = TFOLD_DIR / "setup.py"
    tfold_pyproject = TFOLD_DIR / "pyproject.toml"
    if tfold_setup.exists() or tfold_pyproject.exists():
        print("  Installing tFold as Python package...")
        run_cmd([sys.executable, "-m", "pip", "install", "-e", str(TFOLD_DIR)], check=False)


# ── AlphaFold3 Setup ─────────────────────────────────────────────────────

def setup_af3(force: bool = False) -> None:
    """Clone AF3 repo and prepare directories."""
    print("\n" + "=" * 60)
    print("  Setting up AlphaFold3")
    print("=" * 60)

    if AF3_DIR.exists() and not force:
        print(f"  AlphaFold3 already exists at {AF3_DIR}")
        print("  Use --force to re-clone")
    else:
        if AF3_DIR.exists():
            shutil.rmtree(AF3_DIR)
        TOOLS_DIR.mkdir(parents=True, exist_ok=True)
        run_cmd(["git", "clone", "--depth", "1", _proxy_url(AF3_REPO), str(AF3_DIR)])

    # Create model and database directories
    model_dir = AF3_DIR / "models"
    db_dir = AF3_DIR / "databases"
    model_dir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  AF3 directory: {AF3_DIR}")
    print(f"  Model directory: {model_dir}")
    print(f"  Database directory: {db_dir}")
    print("  Place weights:")
    print(f"    - AF3 model params -> {model_dir}/")
    print("  Note: Request weights at https://forms.gle/svvpY4u2jsHEwWYS6")
    print("  Note: For local (non-Docker) mode, also run:")
    print(f"    cd {AF3_DIR} && python fetch_databases.sh {db_dir}")

    # Check Docker availability for AF3
    try:
        result = subprocess.run(
            ["docker", "--version"], capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            print(f"\n  Docker detected: {result.stdout.strip()}")
            print("  To build AF3 Docker image:")
            print(f"    cd {AF3_DIR}")
            print("    docker build -t alphafold3 -f docker/Dockerfile .")
        else:
            print("\n  Docker not found. AF3 will run in local Python mode.")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("\n  Docker not found. AF3 will run in local Python mode.")
        print("  Set AF3_USE_DOCKER=false in environment.")


# ── ProteinMPNN Setup ────────────────────────────────────────────────────

def setup_mpnn(force: bool = False) -> None:
    """Clone ProteinMPNN repo (weights are included in the repo)."""
    print("\n" + "=" * 60)
    print("  Setting up ProteinMPNN")
    print("=" * 60)

    if MPNN_DIR.exists() and not force:
        print(f"  ProteinMPNN already exists at {MPNN_DIR}")
        print("  Use --force to re-clone")
    else:
        if MPNN_DIR.exists():
            shutil.rmtree(MPNN_DIR)
        TOOLS_DIR.mkdir(parents=True, exist_ok=True)
        run_cmd(["git", "clone", "--depth", "1", _proxy_url(MPNN_REPO), str(MPNN_DIR)])

    # Check if weights are already included
    weights_dir = MPNN_DIR / "vanilla_model_weights"
    if weights_dir.exists() and any(weights_dir.iterdir()):
        print(f"  Weights found at {weights_dir}")
    else:
        weights_dir.mkdir(parents=True, exist_ok=True)
        print(f"  Weights directory created: {weights_dir}")
        print("  ProteinMPNN weights are normally included in the repo.")
        print("  If missing, download from:")
        print("    https://github.com/dauparas/ProteinMPNN/tree/main/vanilla_model_weights")

    # Also check for soluble model weights
    soluble_dir = MPNN_DIR / "soluble_model_weights"
    if not soluble_dir.exists():
        soluble_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n  MPNN directory: {MPNN_DIR}")
    print(f"  Main weights: {weights_dir}")
    print(f"  Run script: {MPNN_DIR / 'protein_mpnn_run.py'}")


# ── Python Dependencies ─────────────────────────────────────────────────

def install_python_deps() -> None:
    """Install shared Python dependencies for Decoy B."""
    print("\n" + "=" * 60)
    print("  Installing Python Dependencies")
    print("=" * 60)

    deps = [
        "numpy",
        "pandas",
        "pyarrow",
        "pydantic>=2.0",
        "biopython",       # For structural comparison (PDB parsing, superimposition)
        "scipy",           # For statistical tests in structure comparison
        "torch",           # Required by tFold and ProteinMPNN
    ]

    for dep in deps:
        print(f"  Installing {dep}...")
        run_cmd(
            [sys.executable, "-m", "pip", "install", dep],
            check=False,
        )

    # Optional: mhcflurry (for HLA filtering in MPNN design branch)
    print("\n  Checking mhcflurry (optional, for MPNN HLA filter)...")
    try:
        import mhcflurry  # noqa: F401
        print("  mhcflurry already installed")
    except ImportError:
        print("  mhcflurry not found. Install with:")
        print("    pip install mhcflurry && mhcflurry-downloads fetch models_class1_pan")


# ── Environment File Generation ──────────────────────────────────────────

def generate_env_file() -> None:
    """Generate .env file with tool paths."""
    project_root = Path(__file__).resolve().parent.parent
    env_path = project_root / ".env"

    content = f"""# Decoy B Tool Paths — auto-generated by setup_decoy_b.py
# Adjust these if you installed tools elsewhere.

# tFold
TFOLD_DIR={TFOLD_DIR}
TFOLD_WEIGHTS_DIR={TFOLD_DIR / 'weights'}

# AlphaFold3
AF3_DIR={AF3_DIR}
AF3_MODEL_DIR={AF3_DIR / 'models'}
AF3_DB_DIR={AF3_DIR / 'databases'}
AF3_USE_DOCKER=false
AF3_DOCKER_IMAGE=alphafold3

# ProteinMPNN
PROTEINMPNN_DIR={MPNN_DIR}

# General
DECOY_TOOLS_DIR={TOOLS_DIR}
PYTHONUTF8=1
"""

    if env_path.exists():
        print(f"\n  .env already exists at {env_path}, writing to .env.decoy_b")
        env_path = project_root / ".env.decoy_b"

    env_path.write_text(content, encoding="utf-8")
    print(f"  Environment file written to {env_path}")


# ── Summary ──────────────────────────────────────────────────────────────

def print_summary() -> None:
    """Print deployment summary."""
    print("\n" + "=" * 60)
    print("  Decoy B Setup Summary")
    print("=" * 60)

    checks = {
        "tFold repo": TFOLD_DIR.exists(),
        "tFold weights dir": (TFOLD_DIR / "weights").exists(),
        "AF3 repo": AF3_DIR.exists(),
        "AF3 models dir": (AF3_DIR / "models").exists(),
        "ProteinMPNN repo": MPNN_DIR.exists(),
        "MPNN weights": (MPNN_DIR / "vanilla_model_weights").exists(),
    }

    for name, ok in checks.items():
        status = "OK" if ok else "MISSING"
        print(f"  [{status:>7}] {name}")

    print("\n  Next steps:")
    print("  1. Place tFold weights in:     ~/tools/tfold/weights/")
    print("  2. Place AF3 weights in:       ~/tools/alphafold3/models/")
    print("  3. (MPNN weights are included in the repo)")
    print("  4. Run verification:           python scripts/verify_decoy_b.py")
    print("  5. Test the pipeline:          python -m decoy_b scan-b --target GILGFVFTL")


# ── Main ─────────────────────────────────────────────────────────────────

def main() -> None:
    global TOOLS_DIR, TFOLD_DIR, AF3_DIR, MPNN_DIR, GITHUB_PROXY

    parser = argparse.ArgumentParser(
        description="Setup external tools for the Decoy B pipeline",
    )
    parser.add_argument("--all", action="store_true", help="Setup all tools")
    parser.add_argument("--tfold", action="store_true", help="Setup tFold only")
    parser.add_argument("--af3", action="store_true", help="Setup AlphaFold3 only")
    parser.add_argument("--mpnn", action="store_true", help="Setup ProteinMPNN only")
    parser.add_argument("--deps-only", action="store_true", help="Install Python deps only")
    parser.add_argument("--force", action="store_true", help="Force re-clone existing repos")
    parser.add_argument(
        "--tools-dir", type=str, default=None,
        help=f"Base directory for tools (default: {TOOLS_DIR})",
    )
    parser.add_argument(
        "--proxy", type=str, default=None,
        help="GitHub proxy URL (e.g., https://ghfast.top/)",
    )

    args = parser.parse_args()

    if args.tools_dir:
        TOOLS_DIR = Path(args.tools_dir)
        TFOLD_DIR = TOOLS_DIR / "tfold"
        AF3_DIR = TOOLS_DIR / "alphafold3"
        MPNN_DIR = TOOLS_DIR / "ProteinMPNN"

    if args.proxy:
        GITHUB_PROXY = args.proxy

    if not check_git():
        print("ERROR: git is not available. Please install git first.")
        sys.exit(1)

    # If no specific flag, default to --all
    if not any([args.all, args.tfold, args.af3, args.mpnn, args.deps_only]):
        args.all = True

    print(f"Tools directory: {TOOLS_DIR}")
    print(f"Platform: {platform.system()} {platform.machine()}")

    if args.deps_only or args.all:
        install_python_deps()

    if args.tfold or args.all:
        setup_tfold(force=args.force)

    if args.af3 or args.all:
        setup_af3(force=args.force)

    if args.mpnn or args.all:
        setup_mpnn(force=args.force)

    generate_env_file()
    print_summary()


if __name__ == "__main__":
    main()
