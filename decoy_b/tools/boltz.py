"""
Boltz Wrapper — Cross-Validation Structure Prediction
=====================================================
Uses Boltz-2 (MIT-licensed, open-source) for pMHC structure prediction as a
cross-validation tool alongside tFold/AF3. Boltz-2 approaches AlphaFold3
accuracy and provides rich confidence metrics (pLDDT, iPTM, PAE).

Supports three execution modes (tried in order):
    1. Python API  — `import boltz` (pip install from local source)
    2. CLI command  — `boltz predict` (installed via pip install -e .)
    3. Subprocess   — `python -m boltz predict` (from BOLTZ_DIR)

Setup:
    git clone <boltz-repo> ~/tools/boltz
    cd ~/tools/boltz && pip install -e .
    # Weights auto-download to ~/.boltz/ on first run

Environment Variables:
    BOLTZ_DIR      — Root of the Boltz repo   (default: ~/tools/boltz)
    BOLTZ_CACHE    — Cache/weights directory   (default: ~/.boltz)
    BOLTZ_MODEL    — Model variant             (default: boltz2)
    BOLTZ_DEVICE   — Device: gpu or cpu        (default: gpu)

Public API:
    check_available()                              -> bool
    predict_pmhc(peptide, hla_allele, ...)         -> BoltzResult
    predict_pmhc_batch(peptides, hla_allele, ...)  -> list[BoltzResult]
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

log = logging.getLogger(__name__)

# ── Configuration ───────────────────────────────────────────────────────

BOLTZ_DIR = Path(os.getenv("BOLTZ_DIR", os.path.expanduser("~/tools/boltz")))
BOLTZ_CACHE = Path(os.getenv("BOLTZ_CACHE", os.path.expanduser("~/.boltz")))
BOLTZ_MODEL = os.getenv("BOLTZ_MODEL", "boltz2")
BOLTZ_DEVICE = os.getenv("BOLTZ_DEVICE", "gpu")


@dataclass
class BoltzResult:
    """Result from a Boltz structure prediction."""
    peptide: str
    hla_allele: str
    cif_path: Optional[str] = None
    pdb_path: Optional[str] = None
    confidence_score: float = 0.0     # 0.8 * pLDDT + 0.2 * PTM
    plddt: float = 0.0               # predicted LDDT (0-1)
    ptm: float = 0.0                 # predicted TM-score (0-1)
    iptm: float = 0.0                # interface TM-score (0-1)
    complex_plddt: float = 0.0       # average pLDDT over entire complex
    success: bool = False
    error: Optional[str] = None


def check_available() -> bool:
    """Check if Boltz is available via any mode."""
    # Mode 1: Python package
    try:
        import boltz  # noqa: F401
        return True
    except ImportError:
        pass

    # Mode 2: CLI command
    try:
        proc = subprocess.run(
            ["boltz", "--help"],
            capture_output=True, text=True, timeout=10,
        )
        if proc.returncode == 0:
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass

    # Mode 3: Subprocess from BOLTZ_DIR
    main_py = BOLTZ_DIR / "src" / "boltz" / "main.py"
    if main_py.exists():
        return True

    return False


def _safe_name(peptide: str, hla_allele: str) -> str:
    """Generate a filesystem-safe name for output files."""
    return f"pmhc_{peptide}_{hla_allele.replace('*', '').replace(':', '')}"


def _build_yaml_input(
    peptide: str,
    mhc_heavy_seq: str,
    b2m_seq: str,
    yaml_path: Path,
) -> None:
    """
    Build a Boltz YAML input for pMHC complex prediction.

    Chain assignment:
        A = MHC heavy chain
        B = Beta-2-microglobulin
        C = Peptide
    """
    yaml_content = (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        f"      id: A\n"
        f"      sequence: {mhc_heavy_seq}\n"
        "  - protein:\n"
        f"      id: B\n"
        f"      sequence: {b2m_seq}\n"
        "  - protein:\n"
        f"      id: C\n"
        f"      sequence: {peptide}\n"
    )
    yaml_path.write_text(yaml_content, encoding="utf-8")


def _parse_confidence_json(output_dir: Path, stem: str) -> Dict[str, float]:
    """Parse Boltz confidence metrics from output JSON."""
    # Boltz outputs: confidence_<stem>_model_<N>.json
    candidates = []
    for p in output_dir.rglob(f"confidence_{stem}*.json"):
        candidates.append(p)
    for p in output_dir.rglob("confidence_*.json"):
        if p not in candidates:
            candidates.append(p)

    for conf_file in candidates:
        try:
            data = json.loads(conf_file.read_text(encoding="utf-8"))
            return {
                "confidence_score": data.get("confidence_score", 0.0),
                "plddt": data.get("complex_plddt", 0.0),
                "ptm": data.get("ptm", 0.0),
                "iptm": data.get("iptm", data.get("protein_iptm", 0.0)),
                "complex_plddt": data.get("complex_plddt", 0.0),
            }
        except (json.JSONDecodeError, KeyError):
            continue

    return {}


def _find_output_structure(output_dir: Path, stem: str) -> Optional[Path]:
    """Find the output CIF or PDB file from Boltz."""
    # Boltz outputs: <stem>_model_0.cif (or .pdb)
    for suffix in [".cif", ".pdb"]:
        for pattern in [f"{stem}_model*{suffix}", f"*model*{suffix}"]:
            matches = list(output_dir.rglob(pattern))
            if matches:
                return matches[0]

    # Fallback: any structure file
    for suffix in ["*.cif", "*.pdb"]:
        matches = list(output_dir.rglob(suffix))
        if matches:
            return matches[0]

    return None


def _cif_to_pdb(cif_path: Path) -> Optional[Path]:
    """Convert mmCIF to PDB using BioPython (if available)."""
    pdb_path = cif_path.with_suffix(".pdb")
    if pdb_path.exists():
        return pdb_path

    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("pmhc", str(cif_path))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        return pdb_path
    except ImportError:
        log.debug("BioPython not available for CIF->PDB conversion")
        return None
    except Exception as e:
        log.debug("CIF->PDB conversion failed: %s", e)
        return None


# ── Mode 1: Python API ──────────────────────────────────────────────────

def _predict_python_api(
    yaml_path: Path,
    output_dir: Path,
) -> tuple[bool, str]:
    """Run Boltz via Python API (import boltz)."""
    try:
        from boltz.main import predict as boltz_predict

        boltz_predict(
            data=str(yaml_path),
            out_dir=str(output_dir),
            cache=str(BOLTZ_CACHE),
            accelerator=BOLTZ_DEVICE,
            model=BOLTZ_MODEL,
            output_format="mmcif",
            devices=1,
            recycling_steps=3,
            sampling_steps=200,
            diffusion_samples=1,
            override=False,
            use_msa_server=False,
        )
        return True, ""
    except ImportError:
        return False, "boltz package not importable"
    except Exception as e:
        return False, f"Boltz Python API error: {e}"


# ── Mode 2: CLI Command ─────────────────────────────────────────────────

def _predict_cli(
    yaml_path: Path,
    output_dir: Path,
) -> tuple[bool, str]:
    """Run Boltz via `boltz predict` CLI."""
    cmd = [
        "boltz", "predict",
        str(yaml_path),
        "--out_dir", str(output_dir),
        "--cache", str(BOLTZ_CACHE),
        "--accelerator", BOLTZ_DEVICE,
        "--model", BOLTZ_MODEL,
        "--output_format", "mmcif",
        "--devices", "1",
        "--recycling_steps", "3",
        "--sampling_steps", "200",
        "--diffusion_samples", "1",
    ]

    try:
        log.debug("Boltz CLI: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800,
            env={**os.environ, "PYTHONUTF8": "1"},
        )
        if proc.returncode != 0:
            return False, f"Boltz CLI failed (rc={proc.returncode}): {proc.stderr[:500]}"
        return True, ""
    except FileNotFoundError:
        return False, "boltz command not found on PATH"
    except subprocess.TimeoutExpired:
        return False, "Boltz CLI timed out (30 min)"


# ── Mode 3: Subprocess from BOLTZ_DIR ───────────────────────────────────

def _predict_subprocess(
    yaml_path: Path,
    output_dir: Path,
) -> tuple[bool, str]:
    """Run Boltz via python -m boltz from BOLTZ_DIR."""
    main_py = BOLTZ_DIR / "src" / "boltz" / "main.py"
    if not main_py.exists():
        return False, f"Boltz main.py not found at {main_py}"

    cmd = [
        "python", "-m", "boltz", "predict",
        str(yaml_path),
        "--out_dir", str(output_dir),
        "--cache", str(BOLTZ_CACHE),
        "--accelerator", BOLTZ_DEVICE,
        "--model", BOLTZ_MODEL,
        "--output_format", "mmcif",
        "--devices", "1",
        "--recycling_steps", "3",
        "--sampling_steps", "200",
        "--diffusion_samples", "1",
    ]

    try:
        log.debug("Boltz subprocess: %s", " ".join(cmd))
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800,
            cwd=str(BOLTZ_DIR),
            env={
                **os.environ,
                "PYTHONUTF8": "1",
                "PYTHONPATH": str(BOLTZ_DIR / "src"),
            },
        )
        if proc.returncode != 0:
            return False, f"Boltz subprocess failed (rc={proc.returncode}): {proc.stderr[:500]}"
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "Boltz subprocess timed out (30 min)"


# ── Public API ───────────────────────────────────────────────────────────

def predict_pmhc(
    peptide: str,
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
    mhc_heavy_seq: Optional[str] = None,
    b2m_seq: Optional[str] = None,
) -> BoltzResult:
    """
    Predict pMHC complex structure using Boltz-2.

    Tries Python API -> CLI -> Subprocess in order.

    Parameters
    ----------
    peptide : str
        Peptide sequence (8-15 AA).
    hla_allele : str
        HLA allele name (e.g., HLA-A*02:01).
    output_dir : Path, optional
        Directory for output files.
    mhc_heavy_seq : str, optional
        Custom MHC heavy chain sequence (overrides built-in).
    b2m_seq : str, optional
        Custom beta-2-microglobulin sequence.

    Returns
    -------
    BoltzResult
    """
    peptide = peptide.strip().upper()

    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "boltz"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve MHC sequences (reuse from tfold module)
    if mhc_heavy_seq is None or b2m_seq is None:
        from .tfold import _get_mhc_sequences
        mhc_seqs = _get_mhc_sequences(hla_allele)
        if mhc_seqs is None:
            return BoltzResult(
                peptide=peptide, hla_allele=hla_allele,
                error=f"No MHC sequence available for {hla_allele}. "
                      f"Supported alleles are defined in tfold.HLA_SEQUENCES.",
            )
        if mhc_heavy_seq is None:
            mhc_heavy_seq = mhc_seqs["mhc_heavy"]
        if b2m_seq is None:
            b2m_seq = mhc_seqs["b2m"]

    name = _safe_name(peptide, hla_allele)

    # Boltz creates output under: output_dir/boltz_results_<stem>/predictions/<stem>/
    # Check cache: look for existing prediction
    for result_subdir in [
        output_dir / f"boltz_results_{name}" / "predictions" / name,
        output_dir / f"boltz_results_{name}" / "predictions",
        output_dir / name,
        output_dir,
    ]:
        if result_subdir.exists():
            cached = _find_output_structure(result_subdir, name)
            if cached is not None:
                log.debug("Using cached Boltz result: %s", cached)
                conf = _parse_confidence_json(result_subdir, name)
                pdb_path = _cif_to_pdb(cached) if cached.suffix == ".cif" else cached
                return BoltzResult(
                    peptide=peptide, hla_allele=hla_allele,
                    cif_path=str(cached) if cached.suffix == ".cif" else None,
                    pdb_path=str(pdb_path) if pdb_path else None,
                    success=True, **conf,
                )

    # Build YAML input
    input_dir = output_dir / "_boltz_inputs"
    input_dir.mkdir(exist_ok=True)
    yaml_path = input_dir / f"{name}.yaml"
    _build_yaml_input(peptide, mhc_heavy_seq, b2m_seq, yaml_path)

    # Try each mode in order
    modes = [
        ("Python API", _predict_python_api),
        ("CLI", _predict_cli),
        ("Subprocess", _predict_subprocess),
    ]

    last_error = ""
    for mode_name, mode_fn in modes:
        try:
            success, err = mode_fn(yaml_path, output_dir)
            if success:
                # Find output
                # Boltz creates: output_dir/boltz_results_<stem>/predictions/<stem>/
                pred_dir = output_dir / f"boltz_results_{name}" / "predictions" / name
                if not pred_dir.exists():
                    # Try broader search
                    pred_dir = output_dir / f"boltz_results_{name}" / "predictions"
                if not pred_dir.exists():
                    pred_dir = output_dir

                struct_path = _find_output_structure(pred_dir, name)
                if struct_path is None:
                    last_error = f"Boltz {mode_name} succeeded but no structure output found"
                    log.debug(last_error)
                    continue

                conf = _parse_confidence_json(pred_dir, name)
                pdb_path = _cif_to_pdb(struct_path) if struct_path.suffix == ".cif" else struct_path

                return BoltzResult(
                    peptide=peptide, hla_allele=hla_allele,
                    cif_path=str(struct_path) if struct_path.suffix == ".cif" else None,
                    pdb_path=str(pdb_path) if pdb_path else str(struct_path),
                    success=True, **conf,
                )

            last_error = err
            log.debug("Boltz %s failed: %s", mode_name, err)
        except Exception as e:
            last_error = str(e)
            log.debug("Boltz %s error: %s", mode_name, e)

    return BoltzResult(
        peptide=peptide, hla_allele=hla_allele,
        error=f"All Boltz modes failed. Last error: {last_error}",
    )


def predict_pmhc_batch(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
) -> List[BoltzResult]:
    """
    Predict pMHC structures for multiple peptides using Boltz.

    Parameters
    ----------
    peptides : list[str]
        Peptide sequences.
    hla_allele : str
        HLA allele.
    output_dir : Path, optional
        Output directory for structure files.

    Returns
    -------
    list[BoltzResult]
    """
    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "boltz"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, pep in enumerate(peptides):
        log.info(
            "Boltz prediction %d/%d: %s (%s)",
            i + 1, len(peptides), pep, hla_allele,
        )
        result = predict_pmhc(pep, hla_allele, output_dir=output_dir)
        results.append(result)

        if not result.success:
            log.warning("  Failed: %s", result.error)

    n_ok = sum(1 for r in results if r.success)
    log.info("Boltz batch: %d/%d succeeded", n_ok, len(results))
    return results
