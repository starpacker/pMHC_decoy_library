"""
AlphaFold3 Wrapper — High-Accuracy Structure Prediction
========================================================
Uses AlphaFold3 for precise pMHC complex structure prediction.
Intended for the refinement stage (Top 200 candidates) after tFold bulk screening.

Supports two execution modes:
    1. Docker (recommended for production)
    2. Local Python (requires manual AF3 install)

Setup (Docker mode):
    git clone https://github.com/google-deepmind/alphafold3.git
    cd alphafold3 && docker build -t alphafold3 -f docker/Dockerfile .
    # Request weights: https://forms.gle/svvpY4u2jsHEwWYS6
    # Download databases: ./fetch_databases.sh <DB_DIR>

Setup (Local mode):
    git clone https://github.com/google-deepmind/alphafold3.git
    cd alphafold3 && pip install -e .
    # Place model weights in $AF3_MODEL_DIR

Environment Variables:
    AF3_DIR          — Root of the AF3 repo      (default: ~/tools/alphafold3)
    AF3_MODEL_DIR    — Model weights directory    (default: $AF3_DIR/models)
    AF3_DB_DIR       — Sequence databases dir     (default: $AF3_DIR/databases)
    AF3_USE_DOCKER   — Use Docker mode            (default: true)
    AF3_DOCKER_IMAGE — Docker image name          (default: alphafold3)
    AF3_GPU_DEVICE   — GPU device for Docker      (default: all)

Public API:
    check_available()                            -> bool
    predict_pmhc(peptide, hla_allele, ...)       -> AF3Result
    predict_pmhc_batch(peptides, hla_allele, ...) -> list[AF3Result]
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

log = logging.getLogger(__name__)

# ── Configuration ───────────────────────────────────────────────────────
from decoy_a.config import AF3_DIR, AF3_MODEL_DIR, AF3_DB_DIR

AF3_DOCKER_IMAGE = os.getenv("AF3_DOCKER_IMAGE", "alphafold3")
AF3_USE_DOCKER = os.getenv("AF3_USE_DOCKER", "true").lower() in ("true", "1", "yes")
AF3_GPU_DEVICE = os.getenv("AF3_GPU_DEVICE", "all")


@dataclass
class AF3Result:
    """Result from an AlphaFold3 prediction."""
    peptide: str
    hla_allele: str
    cif_path: Optional[str] = None
    pdb_path: Optional[str] = None
    plddt: float = 0.0
    ptm: float = 0.0
    iptm: float = 0.0
    ranking_score: float = 0.0
    success: bool = False
    error: Optional[str] = None


def check_available() -> bool:
    """Check if AlphaFold3 is available (Docker or local)."""
    if AF3_USE_DOCKER:
        try:
            result = subprocess.run(
                ["docker", "images", "-q", AF3_DOCKER_IMAGE],
                capture_output=True, text=True, timeout=10,
            )
            return bool(result.stdout.strip())
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    else:
        # Check for local installation
        run_script = AF3_DIR / "run_alphafold.py"
        if run_script.exists() and AF3_MODEL_DIR.exists():
            return True
        # Also check pip-installed alphafold3
        try:
            import alphafold3  # noqa: F401
            return True
        except ImportError:
            pass
        return False


def _build_input_json(
    peptide: str,
    mhc_heavy_seq: str,
    b2m_seq: str,
    job_name: str,
    seeds: Optional[List[int]] = None,
) -> dict:
    """
    Build AlphaFold3 input JSON for pMHC complex prediction.

    Chain assignment:
        A = MHC heavy chain (alpha chain)
        B = Beta-2-microglobulin
        C = Peptide
    """
    if seeds is None:
        seeds = [1]

    return {
        "name": job_name,
        "modelSeeds": seeds,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": mhc_heavy_seq,
                }
            },
            {
                "protein": {
                    "id": "B",
                    "sequence": b2m_seq,
                }
            },
            {
                "protein": {
                    "id": "C",
                    "sequence": peptide,
                }
            },
        ],
        "dialect": "alphafold3",
        "version": 2,
    }


def _parse_confidence(output_dir: Path, job_name: str) -> Dict[str, float]:
    """Parse AF3 confidence metrics from output JSON."""
    # Try multiple naming conventions
    candidates = [
        output_dir / f"{job_name}_summary_confidences.json",
        output_dir / "summary_confidences.json",
    ]
    # Also search recursively
    for p in output_dir.rglob("*summary_confidences*.json"):
        candidates.append(p)

    for conf_file in candidates:
        if conf_file.exists():
            try:
                data = json.loads(conf_file.read_text(encoding="utf-8"))
                return {
                    "plddt": data.get("mean_plddt", data.get("plddt", 0.0)),
                    "ptm": data.get("ptm", 0.0),
                    "iptm": data.get("iptm", 0.0),
                    "ranking_score": data.get("ranking_score", 0.0),
                }
            except (json.JSONDecodeError, KeyError):
                continue

    return {}


def _find_output_structure(job_dir: Path, job_name: str) -> Optional[Path]:
    """Find the output CIF or PDB file from AF3."""
    # Try expected naming
    for suffix in [".cif", ".pdb"]:
        for pattern in [f"{job_name}_model{suffix}", f"{job_name}{suffix}", f"*model*{suffix}"]:
            matches = list(job_dir.rglob(pattern))
            if matches:
                return matches[0]

    # Fallback: any CIF or PDB
    for suffix in ["*.cif", "*.pdb"]:
        matches = list(job_dir.rglob(suffix))
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


# ── Docker Execution ─────────────────────────────────────────────────────

def _run_docker(
    input_path: Path,
    output_dir: Path,
) -> tuple[bool, str]:
    """Run AF3 via Docker. Returns (success, error_message)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "docker", "run", "--rm",
        "--volume", f"{input_path.parent}:/root/af_input",
        "--volume", f"{output_dir}:/root/af_output",
        "--volume", f"{AF3_MODEL_DIR}:/root/models",
        "--volume", f"{AF3_DB_DIR}:/root/public_databases",
        f"--gpus={AF3_GPU_DEVICE}",
        AF3_DOCKER_IMAGE,
        "python", "run_alphafold.py",
        f"--json_path=/root/af_input/{input_path.name}",
        "--model_dir=/root/models",
        "--output_dir=/root/af_output",
        "--run_data_pipeline=true",
        "--run_inference=true",
    ]

    try:
        log.info("Running AF3 Docker: %s", input_path.stem)
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=3600,  # 1 hour max
        )
        if proc.returncode != 0:
            return False, f"Docker failed (rc={proc.returncode}): {proc.stderr[:500]}"
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "AF3 Docker timed out (60 min)"
    except FileNotFoundError:
        return False, "Docker not found on PATH"


# ── Local Python Execution ───────────────────────────────────────────────

def _run_local(
    input_path: Path,
    output_dir: Path,
) -> tuple[bool, str]:
    """Run AF3 locally (without Docker). Returns (success, error_message)."""
    run_script = AF3_DIR / "run_alphafold.py"
    if not run_script.exists():
        return False, f"AF3 script not found: {run_script}"

    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python", str(run_script),
        f"--json_path={input_path}",
        f"--model_dir={AF3_MODEL_DIR}",
        f"--output_dir={output_dir}",
        "--run_data_pipeline=true",
        "--run_inference=true",
    ]
    if AF3_DB_DIR.exists():
        cmd.append(f"--db_dir={AF3_DB_DIR}")

    try:
        log.info("Running AF3 local: %s", input_path.stem)
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
            timeout=3600,
            env={**os.environ, "PYTHONUTF8": "1"},
        )
        if proc.returncode != 0:
            return False, f"AF3 local failed (rc={proc.returncode}): {proc.stderr[:500]}"
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "AF3 local timed out (60 min)"


# ── Public API ───────────────────────────────────────────────────────────

def predict_pmhc(
    peptide: str,
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
    mhc_heavy_seq: Optional[str] = None,
    b2m_seq: Optional[str] = None,
    seeds: Optional[List[int]] = None,
) -> AF3Result:
    """
    Predict pMHC complex structure using AlphaFold3.

    Parameters
    ----------
    peptide : str
        Peptide sequence (8-15 AA).
    hla_allele : str
        HLA allele.
    output_dir : Path, optional
        Output directory.
    mhc_heavy_seq, b2m_seq : str, optional
        Custom MHC sequences (defaults to built-in from tfold module).
    seeds : list[int], optional
        Random seeds for prediction.

    Returns
    -------
    AF3Result
    """
    peptide = peptide.strip().upper()

    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve MHC sequences (reuse from tfold module)
    if mhc_heavy_seq is None or b2m_seq is None:
        from .tfold import _get_mhc_sequences
        mhc_seqs = _get_mhc_sequences(hla_allele)
        if mhc_seqs is None:
            return AF3Result(
                peptide=peptide, hla_allele=hla_allele,
                error=f"No MHC sequence for {hla_allele}",
            )
        if mhc_heavy_seq is None:
            mhc_heavy_seq = mhc_seqs["mhc_heavy"]
        if b2m_seq is None:
            b2m_seq = mhc_seqs["b2m"]

    job_name = f"pmhc_{peptide}_{hla_allele.replace('*', '').replace(':', '')}"
    job_dir = output_dir / job_name

    # Check cache
    cached_struct = _find_output_structure(job_dir, job_name) if job_dir.exists() else None
    if cached_struct is not None:
        log.debug("Using cached AF3 result: %s", cached_struct)
        conf = _parse_confidence(job_dir, job_name)
        pdb_path = _cif_to_pdb(cached_struct) if cached_struct.suffix == ".cif" else cached_struct
        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            cif_path=str(cached_struct) if cached_struct.suffix == ".cif" else None,
            pdb_path=str(pdb_path) if pdb_path else None,
            success=True, **conf,
        )

    # Build input JSON
    input_json = _build_input_json(
        peptide, mhc_heavy_seq, b2m_seq, job_name, seeds,
    )
    input_dir = output_dir / "_af3_inputs"
    input_dir.mkdir(exist_ok=True)
    input_path = input_dir / f"{job_name}.json"
    input_path.write_text(
        json.dumps(input_json, indent=2), encoding="utf-8",
    )

    # Run AF3
    job_dir.mkdir(parents=True, exist_ok=True)
    if AF3_USE_DOCKER:
        success, err = _run_docker(input_path, job_dir)
    else:
        success, err = _run_local(input_path, job_dir)

    if not success:
        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            error=err,
        )

    # Parse output
    struct_path = _find_output_structure(job_dir, job_name)
    if struct_path is None:
        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            error="AF3 produced no structure output",
        )

    conf = _parse_confidence(job_dir, job_name)
    pdb_path = _cif_to_pdb(struct_path) if struct_path.suffix == ".cif" else struct_path

    return AF3Result(
        peptide=peptide, hla_allele=hla_allele,
        cif_path=str(struct_path) if struct_path.suffix == ".cif" else None,
        pdb_path=str(pdb_path) if pdb_path else str(struct_path),
        success=True, **conf,
    )


def predict_pmhc_batch(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
) -> List[AF3Result]:
    """Predict pMHC structures for multiple peptides."""
    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"

    results = []
    for i, pep in enumerate(peptides):
        log.info("AF3 prediction %d/%d: %s (%s)", i + 1, len(peptides), pep, hla_allele)
        result = predict_pmhc(pep, hla_allele, output_dir=output_dir)
        results.append(result)
        if not result.success:
            log.warning("  Failed: %s", result.error)

    n_ok = sum(1 for r in results if r.success)
    log.info("AF3 batch: %d/%d succeeded", n_ok, len(results))
    return results
