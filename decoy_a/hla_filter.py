"""
Step 2: MHC Presentation Gate (The HLA Gate)
=============================================
Filters the ~11M human k-mer dictionary down to ~100-200K peptides that
are predicted to be presented by the target HLA allele.

Uses **NetMHCpan 4.1 EL (Eluted Ligand)** mode, which implicitly models
proteasomal cleavage and TAP transport.

Supports two backends:
    1. Local NetMHCpan binary (preferred for production)
    2. IEDB MHC-I prediction API (fallback, rate-limited)

Public API
----------
    run_hla_filter(hla_allele, kmer_df) -> pd.DataFrame
    load_hla_filtered()                 -> pd.DataFrame
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import List, Optional

from .config import (
    AB_DATA_DIR,
    DEFAULT_HLA_ALLELE,
    HLA_FILTERED_PATH,
    NETMHCPAN_BATCH_SIZE,
    NETMHCPAN_CMD,
    NETMHCPAN_EL_RANK_THRESHOLD,
    NETMHCPAN_STRONG_THRESHOLD,
)
from .models import BindingStrength

log = logging.getLogger(__name__)


# ── NetMHCpan Local Backend ──────────────────────────────────────────────

def _check_netmhcpan() -> bool:
    """Check if NetMHCpan is available on PATH."""
    try:
        result = subprocess.run(
            [NETMHCPAN_CMD, "-h"],
            capture_output=True, text=True, timeout=10,
        )
        return result.returncode == 0 or "NetMHCpan" in result.stderr
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def _run_netmhcpan_batch(
    peptides: List[str],
    hla_allele: str,
    batch_id: int = 0,
) -> List[dict]:
    """
    Run NetMHCpan 4.1 in EL mode on a batch of peptides.

    Parameters
    ----------
    peptides : list[str]
        Peptide sequences to predict.
    hla_allele : str
        HLA allele in NetMHCpan format (e.g., HLA-A02:01).
    batch_id : int
        Batch number for logging.

    Returns
    -------
    list[dict] with keys: sequence, el_rank, el_score, binding
    """
    # NetMHCpan expects allele format without "HLA-" prefix in some versions
    allele = hla_allele.replace("HLA-", "").replace("*", "")

    results = []
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".pep", delete=False
    ) as pep_file:
        for pep in peptides:
            pep_file.write(pep + "\n")
        pep_path = pep_file.name

    try:
        cmd = [
            NETMHCPAN_CMD,
            "-p", pep_path,
            "-a", allele,
            "-BA",  # Include both BA and EL predictions
        ]
        log.debug("Running NetMHCpan batch %d: %d peptides", batch_id, len(peptides))

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour max per batch
        )

        if proc.returncode != 0:
            log.error("NetMHCpan batch %d failed: %s", batch_id, proc.stderr[:500])
            return results

        # Parse NetMHCpan output
        for line in proc.stdout.splitlines():
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("-"):
                continue

            fields = line.split()
            # NetMHCpan output format (EL columns):
            # Pos  MHC  Peptide  Core  Of  Gp  Gl  Ip  Il  Icore  Identity
            # Score_EL  %Rank_EL  Score_BA  Aff(nM)  %Rank_BA  BindLevel
            if len(fields) >= 13:
                try:
                    peptide = fields[2]
                    # Verify it's actually a peptide line
                    if not all(c in "ACDEFGHIKLMNPQRSTVWY" for c in peptide):
                        continue
                    el_score = float(fields[11])
                    el_rank = float(fields[12])

                    binding = "Non_Binder"
                    if el_rank <= NETMHCPAN_STRONG_THRESHOLD:
                        binding = "Strong_Binder"
                    elif el_rank <= NETMHCPAN_EL_RANK_THRESHOLD:
                        binding = "Weak_Binder"

                    results.append({
                        "sequence": peptide,
                        "el_score": el_score,
                        "el_rank": el_rank,
                        "binding": binding,
                    })
                except (ValueError, IndexError):
                    continue

    finally:
        os.unlink(pep_path)

    log.info(
        "NetMHCpan batch %d: %d/%d peptides predicted as binders",
        batch_id,
        sum(1 for r in results if r["binding"] != "Non_Binder"),
        len(peptides),
    )
    return results


# ── IEDB API Backend (Fallback) ──────────────────────────────────────────

IEDB_MHC_URL = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"


def _run_iedb_prediction(
    peptides: List[str],
    hla_allele: str,
    method: str = "netmhcpan_el",
) -> List[dict]:
    """
    Use IEDB's MHC-I prediction API as a fallback.

    Note: This is rate-limited and slower than local NetMHCpan.
    Best for small batches (< 1000 peptides).
    """
    import requests

    # IEDB expects allele format like "HLA-A*02:01"
    allele = hla_allele
    if not allele.startswith("HLA-"):
        allele = f"HLA-{allele}"

    results = []
    # IEDB API processes one request at a time
    sequences_str = "\n".join(peptides)

    try:
        resp = requests.post(
            IEDB_MHC_URL,
            data={
                "method": method,
                "sequence_text": sequences_str,
                "allele": allele,
                "length": ",".join(str(len(p)) for p in set(len(p) for p in peptides)),
            },
            timeout=300,
        )
        resp.raise_for_status()

        # Parse IEDB response (TSV format)
        lines = resp.text.strip().splitlines()
        if len(lines) < 2:
            return results

        header = lines[0].split("\t")
        for line in lines[1:]:
            fields = line.split("\t")
            if len(fields) < len(header):
                continue

            row = dict(zip(header, fields))
            try:
                peptide = row.get("peptide", "")
                rank = float(row.get("percentile_rank", row.get("rank", "100")))
                score = float(row.get("score", "0"))

                binding = "Non_Binder"
                if rank <= NETMHCPAN_STRONG_THRESHOLD:
                    binding = "Strong_Binder"
                elif rank <= NETMHCPAN_EL_RANK_THRESHOLD:
                    binding = "Weak_Binder"

                results.append({
                    "sequence": peptide,
                    "el_score": score,
                    "el_rank": rank,
                    "binding": binding,
                })
            except (ValueError, KeyError):
                continue

    except Exception as exc:
        log.error("IEDB prediction failed: %s", exc)

    return results


# ── Unified HLA Filter ──────────────────────────────────────────────────

def run_hla_filter(
    hla_allele: str = DEFAULT_HLA_ALLELE,
    kmer_df=None,
    force: bool = False,
    use_iedb_fallback: bool = True,
) -> "pd.DataFrame":
    """
    Filter k-mers by HLA presentation using NetMHCpan EL mode.

    Parameters
    ----------
    hla_allele : str
        Target HLA allele (e.g., HLA-A*02:01).
    kmer_df : pd.DataFrame, optional
        Pre-loaded k-mer database. If None, loads from disk.
    force : bool
        Re-run even if cached results exist.
    use_iedb_fallback : bool
        Use IEDB API if local NetMHCpan is not available.

    Returns
    -------
    pd.DataFrame
        Filtered k-mers with EL rank and binding strength columns.
    """
    import pandas as pd

    # Check cache
    cache_path = AB_DATA_DIR / f"hla_filtered_{hla_allele.replace('*', '').replace(':', '')}.parquet"
    if cache_path.exists() and not force:
        log.info("Loading cached HLA-filtered k-mers from %s", cache_path)
        return pd.read_parquet(cache_path)

    # Load k-mer database
    if kmer_df is None:
        from .kmer_builder import load_kmer_db
        kmer_df = load_kmer_db()

    all_peptides = kmer_df["sequence"].tolist()
    total = len(all_peptides)
    log.info("Running HLA filter for %s on %d k-mers...", hla_allele, total)

    # Choose backend: NetMHCpan > mhcflurry > IEDB API
    has_netmhcpan = False
    has_mhcflurry = False

    try:
        from .tools.netmhcpan import check_available as _nmhc_available
        has_netmhcpan = _nmhc_available()
    except ImportError:
        has_netmhcpan = _check_netmhcpan()

    if not has_netmhcpan:
        try:
            from .tools.mhcflurry import check_available as _mhcf_available
            has_mhcflurry = _mhcf_available()
        except ImportError:
            pass

    if has_netmhcpan:
        log.info("Using NetMHCpan (tools.netmhcpan wrapper)")
        try:
            from .tools.netmhcpan import predict_binding_batch
            raw_results = predict_binding_batch(
                all_peptides, hla_allele, batch_size=NETMHCPAN_BATCH_SIZE,
            )
            all_results = []
            for r in raw_results:
                binding = "Non_Binder"
                if r.el_rank <= NETMHCPAN_STRONG_THRESHOLD:
                    binding = "Strong_Binder"
                elif r.el_rank <= NETMHCPAN_EL_RANK_THRESHOLD:
                    binding = "Weak_Binder"
                all_results.append({
                    "sequence": r.sequence,
                    "el_score": r.el_score,
                    "el_rank": r.el_rank,
                    "binding": binding,
                })
        except ImportError:
            log.info("Falling back to legacy NetMHCpan calls")
            all_results = []
            for batch_start in range(0, total, NETMHCPAN_BATCH_SIZE):
                batch = all_peptides[batch_start : batch_start + NETMHCPAN_BATCH_SIZE]
                batch_id = batch_start // NETMHCPAN_BATCH_SIZE
                results = _run_netmhcpan_batch(batch, hla_allele, batch_id)
                all_results.extend(results)

    elif has_mhcflurry:
        log.info("Using mhcflurry (local Python predictor)")
        from .tools.mhcflurry import predict_binding_batch as mhcf_batch
        raw_results = mhcf_batch(
            all_peptides, hla_allele, batch_size=NETMHCPAN_BATCH_SIZE,
        )
        all_results = []
        for r in raw_results:
            binding = "Non_Binder"
            if r.el_rank <= NETMHCPAN_STRONG_THRESHOLD:
                binding = "Strong_Binder"
            elif r.el_rank <= NETMHCPAN_EL_RANK_THRESHOLD:
                binding = "Weak_Binder"
            all_results.append({
                "sequence": r.sequence,
                "el_score": r.presentation_score,
                "el_rank": r.el_rank,
                "binding": binding,
            })

    elif use_iedb_fallback:
        log.warning(
            "NetMHCpan not found locally. Using IEDB API fallback "
            "(much slower, suitable for < 10K peptides). "
            "Install NetMHCpan for production use."
        )
        all_results = []
        # IEDB API: small batches with rate limiting
        iedb_batch_size = 500
        for batch_start in range(0, total, iedb_batch_size):
            batch = all_peptides[batch_start : batch_start + iedb_batch_size]
            results = _run_iedb_prediction(batch, hla_allele)
            all_results.extend(results)
            time.sleep(1.0)  # Rate limit
            if batch_start % 5000 == 0:
                log.info("IEDB progress: %d/%d", batch_start, total)
    else:
        raise RuntimeError(
            "No HLA prediction backend available. "
            "Install NetMHCpan 4.1 or enable IEDB fallback."
        )

    if not all_results:
        log.warning("No HLA binding predictions returned!")
        return pd.DataFrame()

    # Build results DataFrame
    results_df = pd.DataFrame(all_results)

    # Filter to binders only (EL %Rank <= threshold)
    binders_df = results_df[results_df["el_rank"] <= NETMHCPAN_EL_RANK_THRESHOLD].copy()

    # Merge back gene/protein annotations from k-mer DB
    binders_df = binders_df.merge(
        kmer_df[["sequence", "source_proteins", "gene_symbols"]],
        on="sequence",
        how="left",
    )

    # Cache results
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    binders_df.to_parquet(cache_path, index=False, engine="pyarrow")

    n_strong = len(binders_df[binders_df["binding"] == "Strong_Binder"])
    n_weak = len(binders_df[binders_df["binding"] == "Weak_Binder"])
    log.info(
        "HLA filter complete for %s: %d -> %d binders "
        "(%d strong, %d weak) [%.1f%% pass rate]",
        hla_allele, total, len(binders_df),
        n_strong, n_weak,
        100 * len(binders_df) / total if total > 0 else 0,
    )

    return binders_df


def load_hla_filtered(hla_allele: str = DEFAULT_HLA_ALLELE) -> "pd.DataFrame":
    """Load cached HLA-filtered k-mers.

    Prefers the presentation-scored file (``_presentation.parquet``) over the
    affinity-only file, since it includes processing_score and
    presentation_percentile for more accurate filtering.
    """
    import pandas as pd

    allele_tag = hla_allele.replace('*', '').replace(':', '')

    df = None

    # Priority 1: presentation-scored file (97万+)
    presentation_path = AB_DATA_DIR / f"hla_filtered_{allele_tag}_presentation.parquet"
    if presentation_path.exists():
        df = pd.read_parquet(presentation_path)
    else:
        # Priority 2: affinity-only file
        cache_path = AB_DATA_DIR / f"hla_filtered_{allele_tag}.parquet"
        if cache_path.exists():
            df = pd.read_parquet(cache_path)
        else:
            # Priority 3: default path
            if HLA_FILTERED_PATH.exists():
                df = pd.read_parquet(HLA_FILTERED_PATH)

    if df is None:
        raise FileNotFoundError(
            f"No HLA-filtered data for {hla_allele}. Run `run_hla_filter()` first."
        )

    # Filter out Non_Binder explicitly as requested
    if "presentation_binding" in df.columns:
        df = df[df["presentation_binding"] != "Non_Binder"].copy()
    elif "binding" in df.columns:
        df = df[df["binding"] != "Non_Binder"].copy()

    return df
