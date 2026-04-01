#!/usr/bin/env python3
"""
Decoy B Tool Verification Script
==================================
Checks availability and readiness of all external tools required
by the Decoy B structural similarity pipeline.

Usage:
    python scripts/verify_decoy_b.py
    python scripts/verify_decoy_b.py --verbose
    python scripts/verify_decoy_b.py --smoke-test   # Run quick predictions

Exit codes:
    0 — All critical tools available
    1 — Some tools missing (prints what's needed)
"""

from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
from pathlib import Path

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

os.environ.setdefault("PYTHONUTF8", "1")


def _check(name: str, fn, critical: bool = True) -> bool:
    """Run a check and print result."""
    try:
        ok = fn()
        status = "OK" if ok else "MISSING"
        icon = "+" if ok else ("-" if not critical else "!")
        print(f"  [{icon}] {name}: {status}")
        return ok
    except Exception as e:
        icon = "!" if critical else "-"
        print(f"  [{icon}] {name}: ERROR — {e}")
        return False


def check_python_deps() -> dict[str, bool]:
    """Check Python package dependencies."""
    print("\n== Python Dependencies ==")
    deps = {
        "numpy": True,
        "pandas": True,
        "pyarrow": True,
        "pydantic": True,
        "biopython (Bio)": False,  # optional but recommended
        "scipy": False,
        "torch": True,
        "mhcflurry": False,
    }

    results = {}
    for dep, critical in deps.items():
        module_name = dep.split(" ")[0]
        # Special cases
        if module_name == "biopython":
            module_name = "Bio"
        if module_name == "pyarrow":
            module_name = "pyarrow"

        ok = _check(
            dep,
            lambda m=module_name: importlib.import_module(m) is not None,
            critical=critical,
        )
        results[dep] = ok

    return results


def check_tfold() -> bool:
    """Check tFold installation."""
    print("\n== tFold ==")
    from decoy_b.tools.tfold import (
        TFOLD_DIR,
        TFOLD_WEIGHTS_DIR,
        TFOLD_SERVER_URL,
        check_available,
    )

    _check("TFOLD_DIR exists", lambda: TFOLD_DIR.exists())
    _check("Weights directory", lambda: TFOLD_WEIGHTS_DIR.exists())

    # Check for predict script
    for project in ["tfold_ag", "tfold_tcr"]:
        script = TFOLD_DIR / "projects" / project / "predict.py"
        _check(
            f"Predict script ({project})",
            lambda s=script: s.exists(),
            critical=False,
        )

    # Check Python import
    _check(
        "Python package (import tfold)",
        lambda: importlib.import_module("tfold") is not None,
        critical=False,
    )

    if TFOLD_SERVER_URL:
        _check("Server URL", lambda: bool(TFOLD_SERVER_URL), critical=False)

    available = _check("tFold available (any mode)", check_available)
    return available


def check_af3() -> bool:
    """Check AlphaFold3 installation."""
    print("\n== AlphaFold3 ==")
    from decoy_b.tools.alphafold3 import (
        AF3_DIR,
        AF3_DB_DIR,
        AF3_MODEL_DIR,
        AF3_USE_DOCKER,
        check_available,
    )

    _check("AF3_DIR exists", lambda: AF3_DIR.exists())
    _check("Model directory", lambda: AF3_MODEL_DIR.exists())
    _check("Database directory", lambda: AF3_DB_DIR.exists(), critical=False)

    # Check model weights
    has_weights = False
    if AF3_MODEL_DIR.exists():
        weight_files = list(AF3_MODEL_DIR.rglob("*.npz")) + list(AF3_MODEL_DIR.rglob("*.pkl")) + list(AF3_MODEL_DIR.rglob("*.npy"))
        has_weights = len(weight_files) > 0
    _check("Model weights present", lambda: has_weights)

    print(f"  [i] Docker mode: {'enabled' if AF3_USE_DOCKER else 'disabled'}")

    if AF3_USE_DOCKER:
        try:
            import subprocess
            result = subprocess.run(
                ["docker", "images", "-q", "alphafold3"],
                capture_output=True, text=True, timeout=5,
            )
            has_image = bool(result.stdout.strip())
            _check("Docker image 'alphafold3'", lambda: has_image)
        except Exception:
            _check("Docker available", lambda: False)
    else:
        _check(
            "run_alphafold.py",
            lambda: (AF3_DIR / "run_alphafold.py").exists(),
        )

    available = _check("AF3 available", check_available)
    return available


def check_mpnn() -> bool:
    """Check ProteinMPNN installation."""
    print("\n== ProteinMPNN ==")
    from decoy_b.tools.proteinmpnn import (
        PROTEINMPNN_DIR,
        PROTEINMPNN_RUN,
        PROTEINMPNN_WEIGHTS,
        check_available,
    )

    _check("PROTEINMPNN_DIR exists", lambda: PROTEINMPNN_DIR.exists())
    _check("Run script", lambda: PROTEINMPNN_RUN.exists())
    _check("Weights directory", lambda: PROTEINMPNN_WEIGHTS.exists())

    # Check if weights actually contain model files
    if PROTEINMPNN_WEIGHTS.exists():
        weight_files = list(PROTEINMPNN_WEIGHTS.glob("*.pt"))
        _check(
            f"Weight files (.pt): {len(weight_files)} found",
            lambda: len(weight_files) > 0,
        )

    # Check helper scripts
    for script_name in ["parse_multiple_chains.py", "make_fixed_positions_dict.py"]:
        script = PROTEINMPNN_DIR / "helper_scripts" / script_name
        _check(f"Helper: {script_name}", lambda s=script: s.exists(), critical=False)

    available = _check("ProteinMPNN available", check_available)
    return available


def check_hla_backend() -> bool:
    """Check HLA prediction backend (for MPNN filtering)."""
    print("\n== HLA Prediction Backend ==")

    # NetMHCpan
    try:
        from decoy_a.tools.netmhcpan import check_available as nmhc_check
        _check("NetMHCpan", nmhc_check, critical=False)
    except ImportError:
        _check("NetMHCpan (import)", lambda: False, critical=False)

    # mhcflurry
    try:
        from decoy_a.tools.mhcflurry import check_available as mhcf_check
        ok = _check("mhcflurry", mhcf_check, critical=False)
        return ok
    except ImportError:
        _check("mhcflurry (import)", lambda: False, critical=False)
        return False


def check_decoy_b_data() -> None:
    """Check Decoy B data directories."""
    print("\n== Data Directories ==")
    from decoy_a.config import A_DATA_DIR, B_DATA_DIR

    _check("data/decoy_a/", lambda: A_DATA_DIR.exists())
    _check("data/decoy_b/", lambda: B_DATA_DIR.exists(), critical=False)

    # Check for k-mer DB and HLA filtered data (needed by scanner)
    kmer_path = A_DATA_DIR / "human_kmer_db.parquet"
    _check("K-mer database", lambda: kmer_path.exists())

    hla_path = A_DATA_DIR / "hla_filtered_HLA-A0201.parquet"
    _check("HLA-filtered (A*02:01)", lambda: hla_path.exists())

    expr_path = A_DATA_DIR / "gene_expression.parquet"
    _check("Expression database", lambda: expr_path.exists())


def run_smoke_test() -> None:
    """Run a quick smoke test for each tool."""
    print("\n== Smoke Tests ==")

    # Test Atchley factor computation (no external tool needed)
    print("  Testing physicochemical screening...")
    try:
        from decoy_b.scanner import compute_physicochemical_features
        import numpy as np

        target_vec = None
        feat = compute_physicochemical_features("GILGFVFTL")
        print(f"    GILGFVFTL contact residues: {feat.contact_residues}")
        print(f"    Feature vector dim: {len(feat.feature_vector)}")
        print("    [+] Physicochemical computation OK")
    except Exception as e:
        print(f"    [!] Physicochemical test failed: {e}")

    # Test tFold (if available)
    try:
        from decoy_b.tools.tfold import check_available, predict_pmhc
        if check_available():
            print("  Testing tFold prediction...")
            result = predict_pmhc("GILGFVFTL", "HLA-A*02:01")
            if result.success:
                print(f"    PDB: {result.pdb_path}")
                print("    [+] tFold smoke test PASSED")
            else:
                print(f"    [-] tFold smoke test failed: {result.error}")
        else:
            print("  [skip] tFold not available")
    except Exception as e:
        print(f"  [!] tFold smoke test error: {e}")

    # Test ProteinMPNN (if available)
    try:
        from decoy_b.tools.proteinmpnn import check_available
        if check_available():
            print("  [+] ProteinMPNN ready (smoke test requires PDB input)")
        else:
            print("  [skip] ProteinMPNN not available")
    except Exception as e:
        print(f"  [!] ProteinMPNN check error: {e}")

    # Test AF3 (if available)
    try:
        from decoy_b.tools.alphafold3 import check_available
        if check_available():
            print("  [+] AF3 ready (smoke test skipped — too slow)")
        else:
            print("  [skip] AF3 not available")
    except Exception as e:
        print(f"  [!] AF3 check error: {e}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify Decoy B tool installation")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("--smoke-test", action="store_true", help="Run quick predictions")
    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)

    print("=" * 60)
    print("  Decoy B Tool Verification")
    print("=" * 60)

    results = {}

    results["python"] = check_python_deps()
    results["tfold"] = check_tfold()
    results["af3"] = check_af3()
    results["mpnn"] = check_mpnn()
    results["hla"] = check_hla_backend()
    check_decoy_b_data()

    if args.smoke_test:
        run_smoke_test()

    # Summary
    print("\n" + "=" * 60)
    print("  Summary")
    print("=" * 60)

    tool_status = {
        "tFold": results["tfold"],
        "AlphaFold3": results["af3"],
        "ProteinMPNN": results["mpnn"],
        "HLA backend": results["hla"],
    }

    all_ok = True
    for tool, ok in tool_status.items():
        status = "READY" if ok else "NOT READY"
        print(f"  {tool:20s} : {status}")
        if not ok:
            all_ok = False

    if all_ok:
        print("\n  All tools ready! Run:")
        print("    python -m decoy_b scan-b --target GILGFVFTL --hla HLA-A*02:01")
    else:
        print("\n  Some tools are missing. The pipeline will still work with")
        print("  available tools — missing stages will be skipped gracefully.")
        print("  See decoy_b/DEPLOYMENT.md for installation instructions.")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
