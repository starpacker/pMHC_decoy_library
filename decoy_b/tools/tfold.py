"""
tFold Wrapper — TCR-pMHC Structure Prediction
==============================================
Uses tFold (TencentAI4S / DP Tech) for fast pMHC and TCR-pMHC complex
structure prediction. tFold is ~25x faster than AlphaFold3 for this task.

Supports three execution modes (tried in order):
    1. Python API  — `import tfold` (pip install)
    2. CLI script  — `python tfold/projects/tfold_ag/predict.py`
    3. Inference server — HTTP POST to tFold server (if deployed)

Setup:
    git clone https://github.com/TencentAI4S/tfold.git ~/tools/tfold
    cd ~/tools/tfold && pip install -e .
    # Download weights from internal model hub or Zenodo
    # Place under ~/tools/tfold/weights/

Environment Variables:
    TFOLD_DIR          — Root of the tFold repo  (default: ~/tools/tfold)
    TFOLD_WEIGHTS_DIR  — Weight directory         (default: $TFOLD_DIR/weights)
    TFOLD_SERVER_URL   — Inference server URL     (optional, e.g. http://localhost:8501)

Public API:
    check_available()                          -> bool
    predict_pmhc(peptide, hla_allele, ...)     -> TFoldResult
    predict_pmhc_batch(peptides, hla_allele, ...) -> list[TFoldResult]
"""

from __future__ import annotations

import csv
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

TFOLD_DIR = Path(os.getenv("TFOLD_DIR", os.path.expanduser("~/tools/tfold")))
TFOLD_WEIGHTS_DIR = Path(os.getenv("TFOLD_WEIGHTS_DIR", str(TFOLD_DIR / "weights")))
TFOLD_SERVER_URL = os.getenv("TFOLD_SERVER_URL", "")

# ── Reference MHC Sequences ────────────────────────────────────────────
# Full extracellular domain sequences for common HLA alleles.
# Used to construct the input complex for tFold prediction.

HLA_SEQUENCES: Dict[str, Dict[str, str]] = {
    "HLA-A*02:01": {
        "mhc_heavy": (
            "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGE"
            "TRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIA"
            "LNEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHA"
            "VSDHEATL"
        ),
        "b2m": (
            "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWS"
            "FYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
        ),
    },
    "HLA-A*01:01": {
        "mhc_heavy": (
            "GSHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRMEPRAPWIEQEGPEYWDR"
            "NTQIFKTNTQTYRENLRIALRYYNQSEAGSHIIQRMYGCDLGPDGRLLRGHDQYAYDGKDYI"
            "ALNEDLRSWTAADMAAQITKRKWEAAHEAEQLRAYLDGTCVEWLRRYLENGKETLQRADPPKTHVTHHP"
            "ISDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQ"
        ),
        "b2m": (
            "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWS"
            "FYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
        ),
    },
    "HLA-A*24:02": {
        "mhc_heavy": (
            "GSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQ"
            "ETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHTIQIMYGCDVGSDGRFLRGYRQDAYDGKDYI"
            "ALNEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDPPKTHMTHHA"
            "VSDHEATL"
        ),
        "b2m": (
            "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWS"
            "FYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
        ),
    },
    "HLA-B*07:02": {
        "mhc_heavy": (
            "GSHSMRYFYTSMSRPGRGEPRFISVGYVDDTQFVRFDSDAASPREEPRAPWIEQEGPEYWDR"
            "NTQIYKAQAQTDRESLRNLRGYYNQSEAGSHTLQWMHGCELGPDGRFLRGYEQFAYDGKDYI"
            "ALKEDLRSWTAADMAAQITKRKWEAARVAEQLRAYLEGLCVESLRRYLENGKETLQRADPPKTHVTHHP"
            "ISDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQ"
        ),
        "b2m": (
            "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWS"
            "FYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
        ),
    },
}


@dataclass
class TFoldResult:
    """Result from a tFold structure prediction."""
    peptide: str
    hla_allele: str
    pdb_path: Optional[str] = None
    confidence: float = 0.0        # pLDDT or equivalent confidence
    success: bool = False
    error: Optional[str] = None


def check_available() -> bool:
    """Check if tFold is available via any mode."""
    # Mode 1: Python package
    try:
        import tfold  # noqa: F401
        return True
    except ImportError:
        pass

    # Mode 2: CLI script (tfold_ag or tfold_tcr project)
    for project in ["tfold_ag", "tfold_tcr"]:
        predict_script = TFOLD_DIR / "projects" / project / "predict.py"
        if predict_script.exists():
            return True

    # Mode 3: Inference server
    if TFOLD_SERVER_URL:
        try:
            import urllib.request
            req = urllib.request.Request(
                TFOLD_SERVER_URL + "/health",
                method="GET",
            )
            with urllib.request.urlopen(req, timeout=3):
                return True
        except Exception:
            pass

    return False


def _get_mhc_sequences(hla_allele: str) -> Optional[Dict[str, str]]:
    """Get MHC heavy chain and B2M sequences for an HLA allele."""
    if hla_allele in HLA_SEQUENCES:
        return HLA_SEQUENCES[hla_allele]
    # Fuzzy match: strip formatting
    normalized = hla_allele.replace("HLA-", "").replace("*", "").replace(":", "")
    for key, val in HLA_SEQUENCES.items():
        if normalized in key.replace("HLA-", "").replace("*", "").replace(":", ""):
            return val
    return None


def _safe_name(peptide: str, hla_allele: str) -> str:
    """Generate a filesystem-safe name for output files."""
    return f"pmhc_{peptide}_{hla_allele.replace('*', '').replace(':', '')}"


# ── Mode 1: Python API ──────────────────────────────────────────────────

def _predict_python_api(
    peptide: str,
    hla_allele: str,
    mhc_heavy_seq: str,
    b2m_seq: str,
    pdb_path: Path,
) -> TFoldResult:
    """Run tFold via the Python API (pip install tfold)."""
    import tfold

    # Try the pMHC-specific predictor first (tfold >= 1.0)
    try:
        predictor = tfold.deploy.pMHCPredictor.from_pretrained()
    except AttributeError:
        # Fallback for older API: TCRpMHCPredictor
        ppi_model = tfold.model.esm_ppi_650m_tcr()
        trunk_model = tfold.model.tfold_tcr_pmhc_trunk()
        predictor = tfold.deploy.TCRpMHCPredictor(ppi_model, trunk_model)

    data = [
        {"id": "M", "sequence": mhc_heavy_seq},
        {"id": "N", "sequence": b2m_seq},
        {"id": "P", "sequence": peptide},
    ]

    predictor.infer_pdb(data, str(pdb_path))

    if pdb_path.exists():
        return TFoldResult(
            peptide=peptide, hla_allele=hla_allele,
            pdb_path=str(pdb_path), success=True,
        )
    return TFoldResult(
        peptide=peptide, hla_allele=hla_allele,
        error="tFold Python API produced no output",
    )


# ── Mode 2: CLI Script ──────────────────────────────────────────────────

def _find_predict_script() -> Optional[Path]:
    """Find the tFold prediction script."""
    for project in ["tfold_ag", "tfold_tcr"]:
        script = TFOLD_DIR / "projects" / project / "predict.py"
        if script.exists():
            return script
    # Also check root-level scripts
    for name in ["predict.py", "infer.py", "run_prediction.py"]:
        script = TFOLD_DIR / name
        if script.exists():
            return script
    return None


def _predict_cli(
    peptide: str,
    hla_allele: str,
    mhc_heavy_seq: str,
    b2m_seq: str,
    pdb_path: Path,
) -> TFoldResult:
    """Run tFold via command-line script."""
    predict_script = _find_predict_script()
    if predict_script is None:
        return TFoldResult(
            peptide=peptide, hla_allele=hla_allele,
            error=f"tFold predict script not found under {TFOLD_DIR}",
        )

    # Write input JSON
    input_data = {
        "sequences": [
            {"id": "M", "sequence": mhc_heavy_seq},
            {"id": "N", "sequence": b2m_seq},
            {"id": "P", "sequence": peptide},
        ]
    }

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8",
    ) as f:
        json.dump(input_data, f)
        input_path = f.name

    try:
        cmd = [
            "python", str(predict_script),
            "--json", input_path,
            "--output", str(pdb_path.parent),
        ]

        # Add model directory if weights exist
        if TFOLD_WEIGHTS_DIR.exists():
            cmd.extend(["--model_dir", str(TFOLD_WEIGHTS_DIR)])

        # Some tFold versions use --mode flag
        if "tfold_ag" in str(predict_script):
            cmd.extend(["--mode", "pmhc"])

        log.debug("tFold CLI: %s", " ".join(cmd))

        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
            env={**os.environ, "PYTHONUTF8": "1"},
        )

        if proc.returncode != 0:
            return TFoldResult(
                peptide=peptide, hla_allele=hla_allele,
                error=f"tFold CLI failed (rc={proc.returncode}): {proc.stderr[:500]}",
            )

        # Find output PDB
        if pdb_path.exists():
            return TFoldResult(
                peptide=peptide, hla_allele=hla_allele,
                pdb_path=str(pdb_path), success=True,
            )

        # tFold may use a different output naming convention
        for p in pdb_path.parent.glob("*.pdb"):
            if peptide.lower() in p.name.lower() or "pmhc" in p.name.lower():
                p.rename(pdb_path)
                return TFoldResult(
                    peptide=peptide, hla_allele=hla_allele,
                    pdb_path=str(pdb_path), success=True,
                )

        return TFoldResult(
            peptide=peptide, hla_allele=hla_allele,
            error="tFold CLI produced no PDB output",
        )
    finally:
        try:
            os.unlink(input_path)
        except OSError:
            pass


# ── Mode 3: Server API ──────────────────────────────────────────────────

def _predict_server(
    peptide: str,
    hla_allele: str,
    mhc_heavy_seq: str,
    b2m_seq: str,
    pdb_path: Path,
) -> TFoldResult:
    """Run tFold via HTTP inference server."""
    import urllib.request

    payload = json.dumps({
        "sequences": [
            {"id": "M", "sequence": mhc_heavy_seq},
            {"id": "N", "sequence": b2m_seq},
            {"id": "P", "sequence": peptide},
        ],
        "output_format": "pdb",
    }).encode("utf-8")

    req = urllib.request.Request(
        TFOLD_SERVER_URL + "/predict",
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    try:
        with urllib.request.urlopen(req, timeout=300) as resp:
            pdb_content = resp.read().decode("utf-8")
            pdb_path.write_text(pdb_content, encoding="utf-8")
            return TFoldResult(
                peptide=peptide, hla_allele=hla_allele,
                pdb_path=str(pdb_path), success=True,
            )
    except Exception as e:
        return TFoldResult(
            peptide=peptide, hla_allele=hla_allele,
            error=f"tFold server request failed: {e}",
        )


# ── Public API ───────────────────────────────────────────────────────────

def predict_pmhc(
    peptide: str,
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
    mhc_heavy_seq: Optional[str] = None,
    b2m_seq: Optional[str] = None,
) -> TFoldResult:
    """
    Predict pMHC complex structure using tFold.

    Tries Python API -> CLI -> Server in order.

    Parameters
    ----------
    peptide : str
        Peptide sequence (8-15 AA).
    hla_allele : str
        HLA allele name (e.g., HLA-A*02:01).
    output_dir : Path, optional
        Directory for output PDB files.
    mhc_heavy_seq : str, optional
        Custom MHC heavy chain sequence (overrides built-in).
    b2m_seq : str, optional
        Custom beta-2-microglobulin sequence.

    Returns
    -------
    TFoldResult
    """
    peptide = peptide.strip().upper()

    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="tfold_"))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve MHC sequences
    if mhc_heavy_seq is None or b2m_seq is None:
        mhc_seqs = _get_mhc_sequences(hla_allele)
        if mhc_seqs is None:
            return TFoldResult(
                peptide=peptide, hla_allele=hla_allele,
                error=f"No MHC sequence available for {hla_allele}. "
                      f"Supported: {', '.join(HLA_SEQUENCES.keys())}",
            )
        if mhc_heavy_seq is None:
            mhc_heavy_seq = mhc_seqs["mhc_heavy"]
        if b2m_seq is None:
            b2m_seq = mhc_seqs["b2m"]

    name = _safe_name(peptide, hla_allele)
    pdb_path = output_dir / f"{name}.pdb"

    # Return cached result
    if pdb_path.exists():
        log.debug("Using cached tFold result: %s", pdb_path)
        return TFoldResult(
            peptide=peptide, hla_allele=hla_allele,
            pdb_path=str(pdb_path), success=True,
        )

    # Try each mode in order
    # Mode 1: Python API
    try:
        result = _predict_python_api(
            peptide, hla_allele, mhc_heavy_seq, b2m_seq, pdb_path,
        )
        if result.success:
            return result
        log.debug("Python API failed: %s", result.error)
    except ImportError:
        log.debug("tFold Python package not available")
    except Exception as e:
        log.debug("Python API error: %s", e)

    # Mode 2: CLI
    try:
        result = _predict_cli(
            peptide, hla_allele, mhc_heavy_seq, b2m_seq, pdb_path,
        )
        if result.success:
            return result
        log.debug("CLI failed: %s", result.error)
    except Exception as e:
        log.debug("CLI error: %s", e)

    # Mode 3: Server
    if TFOLD_SERVER_URL:
        try:
            result = _predict_server(
                peptide, hla_allele, mhc_heavy_seq, b2m_seq, pdb_path,
            )
            if result.success:
                return result
        except Exception as e:
            log.debug("Server error: %s", e)

    return TFoldResult(
        peptide=peptide, hla_allele=hla_allele,
        error="All tFold modes failed. Check TFOLD_DIR or TFOLD_SERVER_URL.",
    )


def predict_pmhc_batch(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    output_dir: Optional[Path] = None,
) -> List[TFoldResult]:
    """
    Predict pMHC structures for multiple peptides.

    Parameters
    ----------
    peptides : list[str]
        Peptide sequences.
    hla_allele : str
        HLA allele.
    output_dir : Path, optional
        Output directory for PDB files.

    Returns
    -------
    list[TFoldResult]
    """
    if output_dir is None:
        from decoy_a.config import B_DATA_DIR
        output_dir = B_DATA_DIR / "pmhc_models" / "tfold"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for i, pep in enumerate(peptides):
        log.info(
            "tFold prediction %d/%d: %s (%s)",
            i + 1, len(peptides), pep, hla_allele,
        )
        result = predict_pmhc(pep, hla_allele, output_dir=output_dir)
        results.append(result)

        if not result.success:
            log.warning("  Failed: %s", result.error)

    n_ok = sum(1 for r in results if r.success)
    log.info("tFold batch: %d/%d succeeded", n_ok, len(results))
    return results
