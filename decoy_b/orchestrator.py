"""
Decoy A/B Pipeline Orchestrator
================================
End-to-end coordination of the four-step funnel pipeline:

    Step 1  ->  K-mer dictionary + tissue expression mapping
    Step 2  ->  HLA presentation gate (NetMHCpan EL)
    Step 3A ->  Decoy A: sequence homology scan (Hamming <= 2)
    Step 3B ->  Decoy B: structural / physicochemical similarity scan
    Step 4  ->  Composite risk scoring & final ranking

Supports checkpointing (resume from any step), parallel A/B branching,
and audit-quality logging suitable for IND submission packages.

Public API
----------
    run_full_pipeline(target, hla, ...) -> PipelineResult
    run_step(step_name, state)          -> PipelineState
    load_state() / save_state()         -> PipelineState
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from decoy_a.config import (
    AB_DATA_DIR,
    B_DATA_DIR,
    DECOY_A_OUTPUT,
    DECOY_B_OUTPUT,
    DEFAULT_HLA_ALLELE,
    FINAL_RANKED_OUTPUT,
    FINAL_TOP_N,
    HAMMING_DISTANCE_MAX,
    PHYSICOCHEMICAL_COSINE_THRESHOLD,
    PHYSICOCHEMICAL_TOP_K,
    SURFACE_SIMILARITY_THRESHOLD,
)
from decoy_a.models import (
    DecoyABEntry,
    DecoyAHit,
    DecoyBHit,
    PipelineState,
)

log = logging.getLogger(__name__)

# Checkpoint file for resumable runs
_CHECKPOINT_PATH = B_DATA_DIR / "pipeline_checkpoint.json"


# =====================================================================
#  Pipeline result container
# =====================================================================

@dataclass
class PipelineResult:
    """Immutable summary returned after a full pipeline run."""

    target_sequence: str
    hla_allele: str
    total_kmers: int = 0
    hla_filtered: int = 0
    decoy_a_hits: int = 0
    decoy_b_hits: int = 0
    final_entries: int = 0
    top_risk_score: float = 0.0
    top_risk_peptide: str = ""
    output_path: Path = field(default_factory=lambda: FINAL_RANKED_OUTPUT)
    elapsed_seconds: float = 0.0
    entries: List[DecoyABEntry] = field(default_factory=list)

    def summary(self) -> str:
        sep_thick = "=" * 62
        sep_thin = "-" * 62
        lines = [
            "",
            sep_thick,
            "  Decoy A/B Pipeline  --  Run Summary",
            sep_thick,
            "  Target peptide : {}".format(self.target_sequence),
            "  HLA allele     : {}".format(self.hla_allele),
            sep_thin,
            "  Step 1  K-mers generated      : {:>12,}".format(self.total_kmers),
            "  Step 2  HLA-filtered binders   : {:>12,}".format(self.hla_filtered),
            "  Step 3A Decoy A hits           : {:>12,}".format(self.decoy_a_hits),
            "  Step 3B Decoy B hits           : {:>12,}".format(self.decoy_b_hits),
            "  Step 4  Final ranked entries    : {:>12,}".format(self.final_entries),
            sep_thin,
            "  Highest risk    : {:.4f}  ({})".format(
                self.top_risk_score, self.top_risk_peptide
            ),
            "  Elapsed time    : {:.1f}s".format(self.elapsed_seconds),
            "  Output file     : {}".format(self.output_path),
            sep_thick,
            "",
        ]
        return "\n".join(lines)


# =====================================================================
#  Checkpoint helpers
# =====================================================================

def save_state(state: PipelineState) -> Path:
    """Persist pipeline state to a JSON checkpoint."""
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    _CHECKPOINT_PATH.write_text(
        state.model_dump_json(indent=2), encoding="utf-8",
    )
    log.debug("Checkpoint saved to %s", _CHECKPOINT_PATH)
    return _CHECKPOINT_PATH


def load_state() -> Optional[PipelineState]:
    """Load the most recent checkpoint, or return None."""
    if not _CHECKPOINT_PATH.exists():
        return None
    try:
        return PipelineState.model_validate_json(
            _CHECKPOINT_PATH.read_text(encoding="utf-8"),
        )
    except Exception as exc:
        log.warning("Failed to load checkpoint: %s", exc)
        return None


def clear_state() -> None:
    """Remove the checkpoint file."""
    if _CHECKPOINT_PATH.exists():
        _CHECKPOINT_PATH.unlink()
        log.info("Checkpoint cleared")


# =====================================================================
#  Intermediate result persistence
# =====================================================================

def _save_intermediate(data: list, path: Path) -> None:
    """Write intermediate results as JSON for auditability."""
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )
    log.debug("Intermediate results saved to %s (%d records)", path, len(data))


# =====================================================================
#  Individual step runners
# =====================================================================

def _run_step1(state: PipelineState, force: bool = False) -> PipelineState:
    """Step 1: Build k-mer dictionary and expression database."""
    from decoy_a.kmer_builder import build_expression_database, build_kmer_database

    log.info("== Step 1: Building K-mer dictionary & expression database ==")
    t0 = time.time()

    kmer_path = build_kmer_database(force=force)
    expr_path = build_expression_database(force=force)

    # Read counts for audit
    import pandas as pd

    kmer_df = pd.read_parquet(kmer_path)
    state.total_kmers_generated = len(kmer_df)

    expr_df = pd.read_parquet(expr_path)
    state.genes_with_expression = expr_df["gene_symbol"].nunique()

    log.info(
        "Step 1 complete: %d unique k-mers, %d genes with expression (%.1fs)",
        state.total_kmers_generated,
        state.genes_with_expression,
        time.time() - t0,
    )
    return state


def _run_step2(
    state: PipelineState,
    hla_allele: str,
    force: bool = False,
) -> PipelineState:
    """Step 2: HLA presentation gate."""
    from decoy_a.hla_filter import run_hla_filter

    log.info("== Step 2: HLA presentation gate (%s) ==", hla_allele)
    t0 = time.time()

    filtered_df = run_hla_filter(hla_allele=hla_allele, force=force)
    state.hla_filtered_count = len(filtered_df)

    log.info(
        "Step 2 complete: %d binders (%.1fs)",
        state.hla_filtered_count,
        time.time() - t0,
    )
    return state


def _run_step3a(
    state: PipelineState,
    target_sequence: str,
    hla_allele: str,
    max_hamming: int = HAMMING_DISTANCE_MAX,
) -> PipelineState:
    """Step 3A: Decoy A -- sequence homology scan."""
    from decoy_a.scanner import scan_decoy_a

    log.info("== Step 3A: Decoy A -- Sequence homology scan ==")
    t0 = time.time()

    hits = scan_decoy_a(
        target_sequence=target_sequence,
        hla_allele=hla_allele,
        max_hamming=max_hamming,
    )
    state.decoy_a_hits = hits

    # Persist intermediate results
    _save_intermediate(
        [h.model_dump(mode="json") for h in hits],
        DECOY_A_OUTPUT,
    )

    log.info(
        "Step 3A complete: %d Decoy A hits (%.1fs)",
        len(hits),
        time.time() - t0,
    )
    return state


def _run_step3b(
    state: PipelineState,
    target_sequence: str,
    hla_allele: str,
    run_structural: bool = True,
    run_mpnn: bool = False,
    cosine_threshold: float = PHYSICOCHEMICAL_COSINE_THRESHOLD,
    top_k: int = PHYSICOCHEMICAL_TOP_K,
    surface_threshold: float = SURFACE_SIMILARITY_THRESHOLD,
) -> PipelineState:
    """Step 3B: Decoy B -- structural / physicochemical similarity scan."""
    from .scanner import scan_decoy_b

    log.info("== Step 3B: Decoy B -- Structural similarity scan ==")
    t0 = time.time()

    hits = scan_decoy_b(
        target_sequence=target_sequence,
        hla_allele=hla_allele,
        run_structural=run_structural,
        run_mpnn=run_mpnn,
        cosine_threshold=cosine_threshold,
        top_k=top_k,
        surface_threshold=surface_threshold,
    )
    state.decoy_b_hits = hits

    # Track stage counts
    state.physicochemical_candidates = sum(
        1 for h in hits if h.structural is None
    )
    state.structural_candidates = sum(
        1 for h in hits if h.structural is not None
    )

    # Persist intermediate results
    _save_intermediate(
        [h.model_dump(mode="json") for h in hits],
        DECOY_B_OUTPUT,
    )

    log.info(
        "Step 3B complete: %d Decoy B hits "
        "(physchem-only: %d, with structure: %d) (%.1fs)",
        len(hits),
        state.physicochemical_candidates,
        state.structural_candidates,
        time.time() - t0,
    )
    return state


def _run_step4(
    state: PipelineState,
    top_n: int = FINAL_TOP_N,
) -> PipelineState:
    """Step 4: Composite risk scoring and final ranking."""
    from .risk_scorer import save_ranked_results, score_and_rank

    log.info("== Step 4: Composite risk scoring & final ranking ==")
    t0 = time.time()

    entries = score_and_rank(
        decoy_a_hits=state.decoy_a_hits,
        decoy_b_hits=state.decoy_b_hits,
        top_n=top_n,
    )
    state.final_entries = entries
    state.final_top_n = len(entries)

    save_ranked_results(entries)

    log.info(
        "Step 4 complete: %d final ranked entries (%.1fs)",
        len(entries),
        time.time() - t0,
    )
    return state


# =====================================================================
#  Named step dispatch (for selective re-runs)
# =====================================================================

_STEP_NAMES = {
    "step1": "kmer_and_expression",
    "step2": "hla_filter",
    "step3a": "decoy_a",
    "step3b": "decoy_b",
    "step4": "risk_scoring",
    # Aliases
    "kmer": "kmer_and_expression",
    "hla": "hla_filter",
    "a": "decoy_a",
    "b": "decoy_b",
    "score": "risk_scoring",
    "rank": "risk_scoring",
}


def run_step(
    step_name: str,
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    force: bool = False,
    **kwargs,
) -> PipelineState:
    """
    Run a single named pipeline step.

    Parameters
    ----------
    step_name : str
        One of: step1, step2, step3a, step3b, step4
        (or aliases: kmer, hla, a, b, score)
    target_sequence : str
        Target peptide sequence.
    hla_allele : str
        Target HLA allele.
    force : bool
        Re-run even if cached data exists.

    Returns
    -------
    PipelineState
    """
    canonical = _STEP_NAMES.get(step_name.lower())
    if canonical is None:
        raise ValueError(
            "Unknown step '{}'. Valid steps: {}".format(
                step_name, ", ".join(sorted(set(_STEP_NAMES.values())))
            )
        )

    # Load or create state
    state = load_state()
    if state is None:
        state = PipelineState(
            target_sequence=target_sequence,
            hla_allele=hla_allele,
        )

    if canonical == "kmer_and_expression":
        state = _run_step1(state, force=force)
    elif canonical == "hla_filter":
        state = _run_step2(state, hla_allele, force=force)
    elif canonical == "decoy_a":
        state = _run_step3a(state, target_sequence, hla_allele, **kwargs)
    elif canonical == "decoy_b":
        state = _run_step3b(state, target_sequence, hla_allele, **kwargs)
    elif canonical == "risk_scoring":
        state = _run_step4(state, **kwargs)

    save_state(state)
    return state


# =====================================================================
#  Full pipeline
# =====================================================================

def run_full_pipeline(
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    target_protein: str = "",
    force: bool = False,
    skip_structural: bool = False,
    run_mpnn: bool = False,
    top_n: int = FINAL_TOP_N,
    resume: bool = False,
) -> PipelineResult:
    """
    Execute the complete Decoy A/B funnel pipeline.

    Parameters
    ----------
    target_sequence : str
        The therapeutic target peptide (e.g. a MAGE-A3 epitope).
    hla_allele : str
        Target HLA allele (default: HLA-A*02:01).
    target_protein : str
        Name of the target protein (for reporting).
    force : bool
        Re-build all cached databases and re-run all steps.
    skip_structural : bool
        Skip 3D modeling in Decoy B (use physicochemical only).
    run_mpnn : bool
        Enable MPNN inverse design branch in Decoy B.
    top_n : int
        Number of final ranked entries to output.
    resume : bool
        If True, attempt to resume from the last checkpoint.

    Returns
    -------
    PipelineResult
        Summary of the entire run, including final ranked entries.
    """
    target_sequence = target_sequence.strip().upper()
    wall_start = time.time()

    log.info("=" * 62)
    log.info("  Decoy A/B Pipeline -- Starting full run")
    log.info("  Target : %s  (%s)", target_sequence, target_protein or "unnamed")
    log.info("  HLA    : %s", hla_allele)
    log.info("=" * 62)

    # -- Initialise or resume state --------------------------------
    state = None
    if resume:
        state = load_state()
        if state is not None:
            log.info("Resuming from checkpoint (target=%s)", state.target_sequence)

    if state is None:
        state = PipelineState(
            target_sequence=target_sequence,
            target_protein=target_protein,
            hla_allele=hla_allele,
        )

    # -- Step 1: K-mer dictionary + expression ---------------------
    if state.total_kmers_generated == 0 or force:
        state = _run_step1(state, force=force)
        save_state(state)

    # -- Step 2: HLA presentation gate -----------------------------
    if state.hla_filtered_count == 0 or force:
        state = _run_step2(state, hla_allele, force=force)
        save_state(state)

    # -- Step 3A: Decoy A (sequence homology) ----------------------
    if not state.decoy_a_hits or force:
        state = _run_step3a(state, target_sequence, hla_allele)
        save_state(state)

    # -- Step 3B: Decoy B (structural similarity) ------------------
    if not state.decoy_b_hits or force:
        state = _run_step3b(
            state,
            target_sequence,
            hla_allele,
            run_structural=not skip_structural,
            run_mpnn=run_mpnn,
        )
        save_state(state)

    # -- Step 4: Risk scoring & final ranking ----------------------
    state = _run_step4(state, top_n=top_n)
    save_state(state)

    # -- Build result summary --------------------------------------
    elapsed = time.time() - wall_start

    top_score = 0.0
    top_peptide = ""
    if state.final_entries:
        top_score = state.final_entries[0].total_risk_score
        top_peptide = state.final_entries[0].sequence

    result = PipelineResult(
        target_sequence=target_sequence,
        hla_allele=hla_allele,
        total_kmers=state.total_kmers_generated,
        hla_filtered=state.hla_filtered_count,
        decoy_a_hits=len(state.decoy_a_hits),
        decoy_b_hits=len(state.decoy_b_hits),
        final_entries=len(state.final_entries),
        top_risk_score=top_score,
        top_risk_peptide=top_peptide,
        output_path=FINAL_RANKED_OUTPUT,
        elapsed_seconds=elapsed,
        entries=state.final_entries,
    )

    log.info(result.summary())
    return result
