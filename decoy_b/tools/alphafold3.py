"""
AlphaFold3 Wrapper — High-Accuracy pMHC Structure Refinement
=============================================================
Optional refinement stage for Decoy B: uses AlphaFold3 to re-predict the
Top-N candidates identified by tFold, providing higher-accuracy structures
for the final structural similarity comparison.

Supports three execution modes (tried in order):
    1. **Python API** — `pip install alphafold3` (preferred on GPU server)
    2. **Docker**     — official AF3 Docker image
    3. **CLI script** — `python alphafold3/run_alphafold.py`

Setup (Python API — recommended):
    pip install alphafold3
    # Place model weights: af3.bin.zst
    # Set AF3_WEIGHTS_PATH or place in $AF3_MODEL_DIR/af3.bin.zst
    # Download databases: python -m alphafold3.download_databases --db_dir $AF3_DB_DIR

Setup (Docker):
    docker pull us-central1-docker.pkg.dev/alphafold3/alphafold3/alphafold3
    # Or build: git clone ... && docker build -t alphafold3 .

Environment Variables:
    AF3_DIR            — Root of the AF3 repo         (default: /share/liuyutian/alphafold3)
    AF3_MODEL_DIR      — Model weights directory       (default: $AF3_DIR/models)
    AF3_WEIGHTS_PATH   — Path to af3.bin.zst           (default: $AF3_MODEL_DIR/af3.bin.zst)
    AF3_DB_DIR         — Sequence databases dir        (default: $AF3_DIR/databases)
    AF3_USE_DOCKER     — Force Docker mode             (default: false)
    AF3_DOCKER_IMAGE   — Docker image name             (default: alphafold3)
    AF3_GPU_DEVICE     — GPU device for Docker         (default: all)

Chain Convention:
    AF3 outputs chains as A (MHC heavy), B (β2m), C (peptide).
    tFold uses M (MHC heavy), N (β2m), P (peptide).
    The wrapper remaps AF3 chains to M/N/P for compatibility with the
    downstream structure comparison pipeline.

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
from decoy_a.config import AF3_DIR, AF3_MODEL_DIR, AF3_DB_DIR, AF3_WEIGHTS_PATH

AF3_DOCKER_IMAGE = os.getenv("AF3_DOCKER_IMAGE", "alphafold3")
AF3_USE_DOCKER = os.getenv("AF3_USE_DOCKER", "false").lower() in ("true", "1", "yes")
AF3_GPU_DEVICE = os.getenv("AF3_GPU_DEVICE", "all")

# Chain remapping: AF3 (A/B/C) → tFold convention (M/N/P)
_AF3_TO_TFOLD_CHAIN = {"A": "M", "B": "N", "C": "P"}


@dataclass
class AF3Result:
    """Result from an AlphaFold3 prediction."""
    peptide: str
    hla_allele: str
    cif_path: Optional[str] = None
    pdb_path: Optional[str] = None      # PDB with tFold-compatible chain IDs (M/N/P)
    plddt: float = 0.0
    ptm: float = 0.0
    iptm: float = 0.0
    ranking_score: float = 0.0
    success: bool = False
    error: Optional[str] = None


def check_available() -> bool:
    """Check if AlphaFold3 is available (Python API, Docker, or local script)."""
    # Mode 1: Python API
    try:
        from alphafold3 import structure_prediction  # noqa: F401
        if AF3_WEIGHTS_PATH.exists():
            return True
    except ImportError:
        pass

    # Mode 2: Docker
    if AF3_USE_DOCKER:
        try:
            result = subprocess.run(
                ["docker", "images", "-q", AF3_DOCKER_IMAGE],
                capture_output=True, text=True, timeout=10,
            )
            return bool(result.stdout.strip())
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

    # Mode 3: Local script
    run_script = AF3_DIR / "run_alphafold.py"
    if run_script.exists() and AF3_MODEL_DIR.exists():
        return True

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

    Chain assignment (AF3 convention):
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
    candidates = [
        output_dir / f"{job_name}_summary_confidences.json",
        output_dir / "summary_confidences.json",
    ]
    for p in output_dir.rglob("*summary_confidences*.json"):
        candidates.append(p)
    # Also check for confidence in the full output JSON
    for p in output_dir.rglob("*confidences*.json"):
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
    for suffix in [".cif", ".pdb"]:
        for pattern in [f"{job_name}_model{suffix}", f"{job_name}{suffix}", f"*model*{suffix}"]:
            matches = list(job_dir.rglob(pattern))
            if matches:
                return matches[0]

    for suffix in ["*.cif", "*.pdb"]:
        matches = list(job_dir.rglob(suffix))
        if matches:
            return matches[0]

    return None


def _remap_chains_to_tfold(input_path: Path, output_path: Path) -> Path:
    """
    Remap AF3 chain IDs (A/B/C) to tFold convention (M/N/P).

    This ensures AF3-produced PDBs are compatible with the downstream
    structure comparison pipeline in scanner.py, which expects M/N/P.
    """
    lines = input_path.read_text(encoding="utf-8").splitlines()
    remapped = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM", "TER")):
            # PDB column 22 (0-indexed: 21) is the chain ID
            if len(line) > 21:
                old_chain = line[21]
                new_chain = _AF3_TO_TFOLD_CHAIN.get(old_chain, old_chain)
                line = line[:21] + new_chain + line[22:]
        remapped.append(line)
    output_path.write_text("\n".join(remapped), encoding="utf-8")
    return output_path


def _cif_to_pdb(cif_path: Path, remap_chains: bool = True) -> Optional[Path]:
    """Convert mmCIF to PDB using BioPython, with optional chain remapping."""
    pdb_suffix = "_remapped.pdb" if remap_chains else ".pdb"
    pdb_path = cif_path.with_name(cif_path.stem + pdb_suffix)
    if pdb_path.exists():
        return pdb_path

    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("pmhc", str(cif_path))

        if remap_chains:
            for model in structure:
                for chain in model:
                    if chain.id in _AF3_TO_TFOLD_CHAIN:
                        chain.id = _AF3_TO_TFOLD_CHAIN[chain.id]

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        return pdb_path
    except ImportError:
        log.debug("BioPython not available for CIF->PDB conversion")
        return None
    except Exception as e:
        log.warning("CIF->PDB conversion failed: %s", e)
        return None


# ── Mode 1: Python API Execution ────────────────────────────────────────

def _run_python_api(
    input_json: dict,
    output_dir: Path,
    job_name: str,
) -> tuple:
    """Run AF3 via the Python API. Returns (success, error_message)."""
    try:
        from alphafold3.structure_prediction import predict_structure

        output_dir.mkdir(parents=True, exist_ok=True)
        input_path = output_dir / f"{job_name}_input.json"
        input_path.write_text(json.dumps(input_json, indent=2), encoding="utf-8")

        log.info("Running AF3 Python API: %s (weights: %s)", job_name, AF3_WEIGHTS_PATH)

        predict_structure(
            input_json_path=str(input_path),
            model_dir=str(AF3_WEIGHTS_PATH.parent),
            output_dir=str(output_dir),
            db_dir=str(AF3_DB_DIR) if AF3_DB_DIR.exists() else None,
        )
        return True, ""

    except ImportError:
        return False, "alphafold3 Python package not installed"
    except Exception as exc:
        return False, f"AF3 Python API failed: {exc}"


# ── Mode 2: Docker Execution ────────────────────────────────────────────

def _run_docker(
    input_path: Path,
    output_dir: Path,
) -> tuple:
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
            timeout=3600,
        )
        if proc.returncode != 0:
            return False, f"Docker failed (rc={proc.returncode}): {proc.stderr[:500]}"
        return True, ""
    except subprocess.TimeoutExpired:
        return False, "AF3 Docker timed out (60 min)"
    except FileNotFoundError:
        return False, "Docker not found on PATH"


# ── Mode 3: Local Script Execution ──────────────────────────────────────

def _run_local(
    input_path: Path,
    output_dir: Path,
) -> tuple:
    """Run AF3 via local script. Returns (success, error_message)."""
    run_script = AF3_DIR / "run_alphafold.py"
    if not run_script.exists():
        return False, f"AF3 script not found: {run_script}"

    output_dir.mkdir(parents=True, exist_ok=True)

    af3_python = "/home/liuyutian/server/miniconda3/envs/af3/bin/python"
    cmd = [
        af3_python if Path(af3_python).exists() else "python", str(run_script),
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
            env={**os.environ, "PYTHONUTF8": "1", "PYTHONPATH": f"{str(AF3_DIR / 'src')}:{os.environ.get('PYTHONPATH', '')}"},
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

    Tries execution modes in order: Python API → Docker → local script.
    The output PDB has chain IDs remapped to tFold convention (M/N/P)
    for compatibility with the structural comparison pipeline.

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
        With pdb_path pointing to a tFold-compatible (M/N/P chains) PDB.
    """
    peptide = peptide.strip().upper()

    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve MHC sequences
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

    # ── Check cache ──────────────────────────────────────────────────────
    cached_struct = _find_output_structure(job_dir, job_name) if job_dir.exists() else None
    if cached_struct is not None:
        log.debug("Using cached AF3 result: %s", cached_struct)
        conf = _parse_confidence(job_dir, job_name)

        if cached_struct.suffix == ".cif":
            pdb_path = _cif_to_pdb(cached_struct, remap_chains=True)
        else:
            # Remap chains in existing PDB
            remapped = cached_struct.with_name(cached_struct.stem + "_remapped.pdb")
            if remapped.exists():
                pdb_path = remapped
            else:
                pdb_path = _remap_chains_to_tfold(cached_struct, remapped)

        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            cif_path=str(cached_struct) if cached_struct.suffix == ".cif" else None,
            pdb_path=str(pdb_path) if pdb_path else None,
            success=True, **conf,
        )

    # ── Build input JSON ─────────────────────────────────────────────────
    input_json = _build_input_json(
        peptide, mhc_heavy_seq, b2m_seq, job_name, seeds,
    )
    input_dir = output_dir / "_af3_inputs"
    input_dir.mkdir(exist_ok=True)
    input_path = input_dir / f"{job_name}.json"
    input_path.write_text(
        json.dumps(input_json, indent=2), encoding="utf-8",
    )

    # ── Run AF3 (try modes in priority order) ────────────────────────────
    job_dir.mkdir(parents=True, exist_ok=True)
    success = False
    err = "No AF3 execution mode available"

    # Mode 1: Python API
    if not AF3_USE_DOCKER:
        try:
            from alphafold3 import structure_prediction  # noqa: F401
            if AF3_WEIGHTS_PATH.exists():
                success, err = _run_python_api(input_json, job_dir, job_name)
        except ImportError:
            pass

    # Mode 2: Docker
    if not success and AF3_USE_DOCKER:
        success, err = _run_docker(input_path, job_dir)

    # Mode 3: Local script
    if not success:
        run_script = AF3_DIR / "run_alphafold.py"
        if run_script.exists():
            success, err = _run_local(input_path, job_dir)

    if not success:
        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            error=err,
        )

    # ── Parse output ─────────────────────────────────────────────────────
    struct_path = _find_output_structure(job_dir, job_name)
    if struct_path is None:
        return AF3Result(
            peptide=peptide, hla_allele=hla_allele,
            error="AF3 produced no structure output",
        )

    conf = _parse_confidence(job_dir, job_name)

    # Convert to PDB with tFold-compatible chain IDs (M/N/P)
    if struct_path.suffix == ".cif":
        pdb_path = _cif_to_pdb(struct_path, remap_chains=True)
    else:
        remapped = struct_path.with_name(struct_path.stem + "_remapped.pdb")
        pdb_path = _remap_chains_to_tfold(struct_path, remapped)

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
    """
    Predict pMHC structures for multiple peptides.

    Each prediction is run sequentially. Cached results are reused
    automatically (skip already-predicted peptides).
    """
    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"

    results = []
    for i, pep in enumerate(peptides):
        log.info("AF3 prediction %d/%d: %s (%s)", i + 1, len(peptides), pep, hla_allele)
        result = predict_pmhc(pep, hla_allele, output_dir=output_dir)
        results.append(result)
        if result.success:
            log.info("  OK: pLDDT=%.1f, ipTM=%.3f, ranking=%.3f",
                     result.plddt, result.iptm, result.ranking_score)
        else:
            log.warning("  Failed: %s", result.error)

    n_ok = sum(1 for r in results if r.success)
    log.info("AF3 batch: %d/%d succeeded", n_ok, len(results))
    return results
