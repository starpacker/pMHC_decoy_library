"""
Step 4: Composite Risk Scoring & Final Ranking
===============================================
Merges Decoy A and Decoy B hits, applies the composite risk formula,
and outputs the final ranked list of high-risk off-target peptides.

Risk Formula:
    Total Risk = Similarity × (1 / EL_Rank) × Vital_Organ_TPM_Weight

Where:
    - Similarity: sequence similarity (A) or physicochemical/structural (B)
    - EL_Rank:    NetMHCpan EL %Rank (lower = stronger binder = higher risk)
    - TPM_Weight: organ expression weight (Heart/Brain high, Testis low)

Public API
----------
    score_and_rank(decoy_a_hits, decoy_b_hits, ...) -> list[DecoyABEntry]
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

from decoy_a.config import (
    AB_DATA_DIR,
    FINAL_RANKED_OUTPUT,
    FINAL_TOP_N,
    LOW_RELEVANCE_TISSUES,
    ORGAN_TPM_WEIGHT_DEFAULT,
    ORGAN_TPM_WEIGHT_HIGH,
    ORGAN_TPM_WEIGHT_LOW,
    TPM_EXPRESSED_THRESHOLD,
    TPM_HIGH_RISK_THRESHOLD,
    VITAL_ORGANS,
)
from decoy_a.models import (
    DecoyABEntry,
    DecoyAHit,
    DecoyBHit,
    DecoySource,
    MismatchDetail,
    TissueExpression,
)

log = logging.getLogger(__name__)


# ── TPM Weight Calculation ───────────────────────────────────────────────

def compute_tpm_weight(expression: Optional[TissueExpression]) -> float:
    """
    Compute the vital-organ TPM weight multiplier.

    Rules:
        - Heart/Brain TPM > 10  → weight = 10.0
        - Any vital organ TPM > 1 → weight = 1.0 + (max_tpm / 10)
        - Only testis/placenta   → weight = 0.1
        - No expression data     → weight = 1.0 (neutral)
    """
    if expression is None:
        return ORGAN_TPM_WEIGHT_DEFAULT

    max_vital = expression.max_vital_organ_tpm
    category = expression.expression_category

    if category == "high_risk" or max_vital >= TPM_HIGH_RISK_THRESHOLD:
        return ORGAN_TPM_WEIGHT_HIGH
    elif category == "restricted":
        return ORGAN_TPM_WEIGHT_LOW
    elif category == "silent":
        return ORGAN_TPM_WEIGHT_LOW
    elif max_vital >= TPM_EXPRESSED_THRESHOLD:
        return 1.0 + (max_vital / 10.0)
    else:
        return ORGAN_TPM_WEIGHT_DEFAULT


def _get_critical_organs(expression: Optional[TissueExpression]) -> List[str]:
    """Identify which vital organs have high expression."""
    if expression is None:
        return []

    vital_set = {t.lower() for t in VITAL_ORGANS}
    organs = []
    for tissue, tpm in expression.tissue_tpm.items():
        if tissue.lower() in vital_set and tpm >= TPM_EXPRESSED_THRESHOLD:
            organs.append(tissue.title())
    return sorted(organs)


# ── Risk Score Calculation ───────────────────────────────────────────────

def compute_risk_score(
    similarity: float,
    el_rank: float,
    tpm_weight: float,
) -> float:
    """
    Compute the composite risk score.

    Total Risk = Similarity × (1 / EL_Rank) × TPM_Weight

    Parameters
    ----------
    similarity : float
        Sequence (Decoy A) or structural (Decoy B) similarity [0, 1].
    el_rank : float
        NetMHCpan EL %Rank (lower = stronger binder).
    tpm_weight : float
        Vital organ expression weight.

    Returns
    -------
    float
        Composite risk score (higher = more dangerous).
    """
    # Clamp EL rank to avoid division by zero
    el_rank = max(el_rank, 0.01)
    return similarity * (1.0 / el_rank) * tpm_weight


# ── Merge & Rank ─────────────────────────────────────────────────────────

def score_and_rank(
    decoy_a_hits: List[DecoyAHit],
    decoy_b_hits: List[DecoyBHit],
    top_n: int = FINAL_TOP_N,
) -> List[DecoyABEntry]:
    """
    Merge Decoy A and B hits, compute composite risk scores, and rank.

    Parameters
    ----------
    decoy_a_hits : list[DecoyAHit]
        Hits from the sequence homology scan.
    decoy_b_hits : list[DecoyBHit]
        Hits from the structural similarity scan.
    top_n : int
        Number of top-ranked entries to return.

    Returns
    -------
    list[DecoyABEntry]
        Final ranked entries with composite risk scores.
    """
    log.info(
        "Scoring and ranking: %d Decoy A hits + %d Decoy B hits",
        len(decoy_a_hits), len(decoy_b_hits),
    )

    entries: List[DecoyABEntry] = []
    seen_sequences: set = set()
    counter = 1

    # ── Process Decoy A hits ─────────────────────────────────────────────
    for hit in decoy_a_hits:
        tpm_weight = compute_tpm_weight(hit.expression)
        risk = compute_risk_score(hit.similarity_score, hit.el_rank, tpm_weight)
        critical_organs = _get_critical_organs(hit.expression)

        entry = DecoyABEntry(
            decoy_ab_id=f"DAB-{counter:04d}",
            sequence=hit.sequence,
            target_sequence=hit.target_sequence,
            hla_allele=hit.hla_allele,
            source=DecoySource.DECOY_A,
            gene_symbols=hit.gene_symbols,
            source_proteins=hit.source_proteins,
            el_rank=hit.el_rank,
            hamming_distance=hit.hamming_distance,
            sequence_similarity=hit.similarity_score,
            expression=hit.expression,
            vital_organ_tpm_weight=tpm_weight,
            critical_organs=critical_organs,
            total_risk_score=risk,
            mismatches=hit.mismatches,
        )
        entries.append(entry)
        seen_sequences.add(hit.sequence)
        counter += 1

    # ── Process Decoy B hits ─────────────────────────────────────────────
    for hit in decoy_b_hits:
        if hit.sequence in seen_sequences:
            # Peptide already found by Decoy A — upgrade source to BOTH
            for e in entries:
                if e.sequence == hit.sequence:
                    e.source = DecoySource.BOTH
                    e.physicochemical_similarity = hit.similarity_score
                    e.structural = hit.structural
                    # Recalculate risk with the better similarity
                    best_sim = max(e.sequence_similarity, hit.similarity_score)
                    e.total_risk_score = compute_risk_score(
                        best_sim, e.el_rank, e.vital_organ_tpm_weight,
                    )
                    break
            continue

        tpm_weight = compute_tpm_weight(hit.expression)
        risk = compute_risk_score(hit.similarity_score, hit.el_rank, tpm_weight)
        critical_organs = _get_critical_organs(hit.expression)

        entry = DecoyABEntry(
            decoy_ab_id=f"DAB-{counter:04d}",
            sequence=hit.sequence,
            target_sequence=hit.target_sequence,
            hla_allele=hit.hla_allele,
            source=DecoySource.DECOY_B,
            gene_symbols=hit.gene_symbols,
            source_proteins=hit.source_proteins,
            el_rank=hit.el_rank,
            hamming_distance=hit.hamming_distance,
            physicochemical_similarity=hit.physicochemical.cosine_similarity,
            structural_similarity=(
                hit.structural.surface_correlation if hit.structural else 0.0
            ),
            expression=hit.expression,
            vital_organ_tpm_weight=tpm_weight,
            critical_organs=critical_organs,
            total_risk_score=risk,
            structural=hit.structural,
        )
        entries.append(entry)
        seen_sequences.add(hit.sequence)
        counter += 1

    # ── Sort by risk score (descending) and assign ranks ─────────────────
    entries.sort(key=lambda e: -e.total_risk_score)

    # Trim to top-N
    entries = entries[:top_n]

    # Assign ranks and re-number IDs
    for i, entry in enumerate(entries):
        entry.risk_rank = i + 1
        entry.decoy_ab_id = f"DAB-{i + 1:04d}"

    log.info(
        "Final ranking: %d entries (A: %d, B: %d, Both: %d)",
        len(entries),
        sum(1 for e in entries if e.source == DecoySource.DECOY_A),
        sum(1 for e in entries if e.source == DecoySource.DECOY_B),
        sum(1 for e in entries if e.source == DecoySource.BOTH),
    )

    if entries:
        log.info(
            "Top risk score: %.4f (%s, %s)",
            entries[0].total_risk_score,
            entries[0].sequence,
            ", ".join(entries[0].gene_symbols) or "unknown gene",
        )

    return entries


# ── Persistence ──────────────────────────────────────────────────────────

def save_ranked_results(
    entries: List[DecoyABEntry],
    output_path: Optional[Path] = None,
) -> Path:
    """Save the final ranked results to JSON."""
    if output_path is None:
        output_path = FINAL_RANKED_OUTPUT

    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)

    data = [e.model_dump(mode="json") for e in entries]
    output_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    log.info("Saved %d ranked entries to %s", len(entries), output_path)
    return output_path


def load_ranked_results(
    input_path: Optional[Path] = None,
) -> List[DecoyABEntry]:
    """Load ranked results from JSON."""
    if input_path is None:
        input_path = FINAL_RANKED_OUTPUT

    if not input_path.exists():
        raise FileNotFoundError(f"No ranked results at {input_path}")

    data = json.loads(input_path.read_text(encoding="utf-8"))
    return [DecoyABEntry.model_validate(d) for d in data]
                       
