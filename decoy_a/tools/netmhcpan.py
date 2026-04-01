"""
NetMHCpan 4.1 Wrapper
=====================
Predicts peptide-MHC-I binding using the Eluted Ligand (EL) model.

Setup:
    1. Download from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
    2. Extract and run `./netMHCpan` to verify
    3. Set NETMHCPAN_DIR in config or environment variable

Public API:
    check_available()         -> bool
    predict_binding(peptides, hla) -> list[NetMHCpanResult]
    predict_binding_batch(peptides, hla, batch_size) -> list[NetMHCpanResult]
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

log = logging.getLogger(__name__)

# ── Configuration ───────────────────────────────────────────────────────

# Where NetMHCpan is installed. Override via env or config.
NETMHCPAN_DIR = Path(os.getenv(
    "NETMHCPAN_DIR",
    os.path.expanduser("~/tools/netMHCpan-4.1"),
))
NETMHCPAN_BIN = NETMHCPAN_DIR / "netMHCpan"

# Thresholds
EL_RANK_STRONG = 0.5    # %Rank <= 0.5 = strong binder
EL_RANK_WEAK = 2.0      # %Rank <= 2.0 = weak binder


@dataclass
class NetMHCpanResult:
    """Single peptide-HLA prediction result."""
    sequence: str
    hla_allele: str
    core: str           # 9-residue binding core
    el_score: float     # Raw EL prediction score
    el_rank: float      # Percentile rank (lower = stronger binder)
    ba_score: float     # Binding affinity score (if -BA flag used)
    ba_rank: float      # BA percentile rank
    affinity_nm: float  # Predicted IC50 in nM (if -BA flag used)
    bind_level: str     # "SB", "WB", or ""


def check_available() -> bool:
    """Check if NetMHCpan is installed and executable."""
    if not NETMHCPAN_BIN.exists():
        # Try system PATH
        try:
            result = subprocess.run(
                ["netMHCpan", "-h"],
                capture_output=True, text=True, timeout=10,
            )
            return "NetMHCpan" in result.stdout or "NetMHCpan" in result.stderr
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    try:
        result = subprocess.run(
            [str(NETMHCPAN_BIN), "-h"],
            capture_output=True, text=True, timeout=10,
            env=_build_env(),
        )
        return result.returncode == 0 or "NetMHCpan" in result.stderr
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _build_env() -> dict:
    """Build environment variables for NetMHCpan execution."""
    env = os.environ.copy()
    env["NHOME"] = str(NETMHCPAN_DIR)
    env["TMPDIR"] = tempfile.gettempdir()
    # NetMHCpan needs its own data directory
    data_dir = NETMHCPAN_DIR / "data"
    if data_dir.exists():
        env["NMHOME"] = str(NETMHCPAN_DIR)
    return env


def _get_binary() -> str:
    """Get the NetMHCpan binary path."""
    if NETMHCPAN_BIN.exists():
        return str(NETMHCPAN_BIN)
    # Fall back to system PATH
    return "netMHCpan"


def _format_allele(hla_allele: str) -> str:
    """
    Convert HLA allele to NetMHCpan format.

    NetMHCpan accepts: HLA-A02:01, HLA-A*02:01, or A0201
    """
    # Normalize: HLA-A*02:01 -> HLA-A02:01
    allele = hla_allele.strip()
    allele = allele.replace("*", "")
    return allele


def _parse_output(stdout: str, hla_allele: str) -> List[NetMHCpanResult]:
    """
    Parse NetMHCpan text output into structured results.

    Output format (with -BA flag):
    Pos  HLA  Peptide  Core  Of  Gp  Gl  Ip  Il  Icore  Identity
    Score_EL  %Rank_EL  Score_BA  Affinity(nM)  %Rank_BA  BindLevel
    """
    results = []
    in_data = False

    for line in stdout.splitlines():
        line = line.strip()

        # Detect data section start
        if line.startswith("---"):
            in_data = True
            continue

        if not in_data or not line or line.startswith("#"):
            continue

        # Skip header lines
        if line.startswith("Pos") or line.startswith("Number"):
            continue

        fields = line.split()
        if len(fields) < 13:
            continue

        try:
            # Verify this is a data line (first field should be a number)
            int(fields[0])
        except ValueError:
            continue

        try:
            peptide = fields[2]
            # Verify it's a valid peptide
            if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in peptide):
                continue

            core = fields[3]
            el_score = float(fields[11])
            el_rank = float(fields[12])

            # BA columns (present if -BA was used)
            ba_score = 0.0
            affinity_nm = 0.0
            ba_rank = 0.0
            bind_level = ""

            if len(fields) >= 16:
                ba_score = float(fields[13])
                affinity_nm = float(fields[14])
                ba_rank = float(fields[15])
            if len(fields) >= 17:
                bind_level = fields[16]  # "SB", "WB", or "<=WB"/"<=SB"

            # Normalize bind_level
            if "SB" in bind_level:
                bind_level = "SB"
            elif "WB" in bind_level:
                bind_level = "WB"
            else:
                bind_level = ""

            results.append(NetMHCpanResult(
                sequence=peptide,
                hla_allele=hla_allele,
                core=core,
                el_score=el_score,
                el_rank=el_rank,
                ba_score=ba_score,
                ba_rank=ba_rank,
                affinity_nm=affinity_nm,
                bind_level=bind_level,
            ))

        except (ValueError, IndexError):
            continue

    return results


def predict_binding(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    include_ba: bool = True,
    timeout: int = 3600,
) -> List[NetMHCpanResult]:
    """
    Predict peptide-MHC-I binding for a list of peptides.

    Parameters
    ----------
    peptides : list[str]
        Peptide sequences (8-15 AA).
    hla_allele : str
        HLA allele (e.g., HLA-A*02:01).
    include_ba : bool
        Include binding affinity prediction (-BA flag).
    timeout : int
        Max seconds to wait.

    Returns
    -------
    list[NetMHCpanResult]
        One result per peptide (non-binders included).
    """
    if not peptides:
        return []

    binary = _get_binary()
    allele = _format_allele(hla_allele)

    # Write peptides to temp file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pep", delete=False, dir=tempfile.gettempdir(),
    ) as f:
        for pep in peptides:
            f.write(pep.strip().upper() + "\n")
        pep_path = f.name

    try:
        cmd = [binary, "-p", pep_path, "-a", allele]
        if include_ba:
            cmd.append("-BA")

        log.debug("Running: %s (%d peptides)", " ".join(cmd[:4]), len(peptides))

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=_build_env(),
        )

        if proc.returncode != 0:
            log.error("NetMHCpan failed (rc=%d): %s", proc.returncode, proc.stderr[:500])
            return []

        return _parse_output(proc.stdout, hla_allele)

    except FileNotFoundError:
        log.error(
            "NetMHCpan binary not found at '%s'. "
            "Set NETMHCPAN_DIR env variable or install NetMHCpan 4.1.",
            binary,
        )
        return []
    except subprocess.TimeoutExpired:
        log.error("NetMHCpan timed out after %ds for %d peptides", timeout, len(peptides))
        return []
    finally:
        os.unlink(pep_path)


def predict_binding_batch(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    batch_size: int = 50000,
    include_ba: bool = False,
) -> List[NetMHCpanResult]:
    """
    Predict binding in batches (for large peptide sets).

    Parameters
    ----------
    peptides : list[str]
        All peptide sequences.
    hla_allele : str
        Target HLA allele.
    batch_size : int
        Peptides per batch (default: 50K).
    include_ba : bool
        Include BA prediction (slower).

    Returns
    -------
    list[NetMHCpanResult]
        Combined results from all batches.
    """
    all_results = []
    total = len(peptides)

    for start in range(0, total, batch_size):
        batch = peptides[start:start + batch_size]
        batch_num = start // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size

        log.info(
            "NetMHCpan batch %d/%d: %d peptides (%s)",
            batch_num, total_batches, len(batch), hla_allele,
        )

        results = predict_binding(batch, hla_allele, include_ba=include_ba)
        all_results.extend(results)

        n_binders = sum(1 for r in results if r.el_rank <= EL_RANK_WEAK)
        log.info(
            "  Batch %d: %d/%d binders (SB: %d, WB: %d)",
            batch_num,
            n_binders,
            len(results),
            sum(1 for r in results if r.el_rank <= EL_RANK_STRONG),
            sum(1 for r in results if EL_RANK_STRONG < r.el_rank <= EL_RANK_WEAK),
        )

    return all_results


def filter_binders(
    results: List[NetMHCpanResult],
    rank_threshold: float = EL_RANK_WEAK,
) -> List[NetMHCpanResult]:
    """Keep only peptides with EL %Rank <= threshold."""
    return [r for r in results if r.el_rank <= rank_threshold]
