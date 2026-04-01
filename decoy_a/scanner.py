"""
Step 3A: Decoy A — Sequence Homology Scan
==========================================
Identifies off-target peptides within Hamming distance ≤ 2 of the target
peptide, from the HLA-filtered candidate pool (~100-200K peptides).

Key features:
    - Exact Hamming distance calculation (vectorised with NumPy)
    - Anchor vs TCR-contact position annotation
    - Expression-weighted ranking
    - Typically yields 50-200 high-confidence hits

Public API
----------
    scan_decoy_a(target_seq, hla_allele, ...) -> list[DecoyAHit]
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Set

import numpy as np

from .config import (
    ANCHOR_POSITIONS_A0201,
    DEFAULT_HLA_ALLELE,
    HAMMING_DISTANCE_MAX,
    TCR_CONTACT_POSITIONS_9MER,
    VITAL_ORGANS,
)
from .models import DecoyAHit, MismatchDetail, TissueExpression

log = logging.getLogger(__name__)


# ── Core Hamming Distance ────────────────────────────────────────────────

def hamming_distance(s1: str, s2: str) -> int:
    """Compute Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))  # Unequal lengths → max distance
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def hamming_distance_vectorised(
    target: str,
    candidates: np.ndarray,
) -> np.ndarray:
    """
    Compute Hamming distances between a target and many candidates.

    Parameters
    ----------
    target : str
        Target peptide sequence.
    candidates : np.ndarray
        Array of candidate sequences (as bytes or str).

    Returns
    -------
    np.ndarray of int
        Hamming distances.
    """
    target_arr = np.array(list(target), dtype="U1")
    n = len(candidates)
    k = len(target)

    # Convert candidates to a 2D character array
    # Each candidate is split into individual characters
    cand_chars = np.array([list(str(s)) for s in candidates], dtype="U1")

    if cand_chars.shape[1] != k:
        # Length mismatch — fall back to scalar
        return np.array([hamming_distance(target, str(s)) for s in candidates])

    # Vectorised comparison
    mismatches = cand_chars != target_arr
    return mismatches.sum(axis=1)


# ── Mismatch Annotation ─────────────────────────────────────────────────

def _get_anchor_positions(hla_allele: str, peptide_length: int) -> Set[int]:
    """
    Get HLA anchor positions (0-indexed) for a given allele and length.

    For HLA-A*02:01 9-mers: positions 1 (p2) and 8 (p9).
    For other lengths, adjust proportionally.
    """
    if "A*02" in hla_allele:
        if peptide_length == 9:
            return {1, 8}
        elif peptide_length == 10:
            return {1, 9}
        elif peptide_length == 11:
            return {1, 10}
        elif peptide_length == 8:
            return {1, 7}
    # Default: first and last position after p1
    return {1, peptide_length - 1}


def _get_tcr_contact_positions(peptide_length: int) -> List[int]:
    """
    Get TCR-facing contact positions (0-indexed).

    For 9-mers: p4-p8 (indices 3-7).
    Scales proportionally for other lengths.
    """
    if peptide_length == 9:
        return [3, 4, 5, 6, 7]
    elif peptide_length == 10:
        return [3, 4, 5, 6, 7, 8]
    elif peptide_length == 11:
        return [3, 4, 5, 6, 7, 8, 9]
    elif peptide_length == 8:
        return [3, 4, 5, 6]
    # Default: middle 60% of positions
    start = max(2, peptide_length // 4)
    end = min(peptide_length - 1, peptide_length * 3 // 4 + 1)
    return list(range(start, end))


def annotate_mismatches(
    target: str,
    candidate: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
) -> List[MismatchDetail]:
    """
    Annotate each mismatch position between target and candidate.

    Returns
    -------
    list[MismatchDetail]
        One entry per mismatch position.
    """
    pep_len = len(target)
    anchors = _get_anchor_positions(hla_allele, pep_len)
    tcr_contacts = set(_get_tcr_contact_positions(pep_len))

    mismatches = []
    for i, (t_aa, c_aa) in enumerate(zip(target, candidate)):
        if t_aa != c_aa:
            mismatches.append(MismatchDetail(
                position=i,
                target_aa=t_aa,
                candidate_aa=c_aa,
                is_anchor=i in anchors,
                is_tcr_contact=i in tcr_contacts,
            ))
    return mismatches


# ── Main Scan ────────────────────────────────────────────────────────────

def scan_decoy_a(
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    hla_filtered_df=None,
    expr_df=None,
    max_hamming: int = HAMMING_DISTANCE_MAX,
) -> List[DecoyAHit]:
    """
    Scan HLA-filtered peptides for sequence-homologous off-targets.

    Parameters
    ----------
    target_sequence : str
        The therapeutic target peptide (e.g., MAGE-A3 epitope).
    hla_allele : str
        Target HLA allele.
    hla_filtered_df : pd.DataFrame, optional
        HLA-filtered k-mers. If None, loads from cache.
    expr_df : pd.DataFrame, optional
        Expression database. If None, loads from cache.
    max_hamming : int
        Maximum Hamming distance to include (default: 2).

    Returns
    -------
    list[DecoyAHit]
        Hits sorted by (hamming_distance ASC, expression risk DESC).
    """
    import pandas as pd

    target_sequence = target_sequence.strip().upper()
    target_len = len(target_sequence)
    log.info(
        "Decoy A scan: target=%s (len=%d), HLA=%s, max_hamming=%d",
        target_sequence, target_len, hla_allele, max_hamming,
    )

    # Load data
    if hla_filtered_df is None:
        from .hla_filter import load_hla_filtered
        hla_filtered_df = load_hla_filtered(hla_allele)

    if expr_df is None:
        try:
            from .kmer_builder import load_expression_db
            expr_df = load_expression_db()
        except FileNotFoundError:
            log.warning("Expression database not available; skipping expression annotation")
            expr_df = None

    # Filter to same-length peptides only
    same_len = hla_filtered_df[
        hla_filtered_df["sequence"].str.len() == target_len
    ].copy()
    log.info("Candidates of length %d: %d", target_len, len(same_len))

    if same_len.empty:
        log.warning("No candidates of matching length found")
        return []

    # Exclude the target itself
    same_len = same_len[same_len["sequence"] != target_sequence]

    # Vectorised Hamming distance
    candidates = same_len["sequence"].values
    distances = hamming_distance_vectorised(target_sequence, candidates)

    # Filter by max Hamming distance
    mask = distances <= max_hamming
    hits_df = same_len[mask].copy()
    hits_df["hamming_distance"] = distances[mask]

    log.info(
        "Decoy A: %d hits within Hamming distance ≤ %d",
        len(hits_df), max_hamming,
    )

    if hits_df.empty:
        return []

    # Build expression lookup
    expr_lookup: Dict[str, dict] = {}
    if expr_df is not None:
        from .kmer_builder import get_gene_expression
        # Pre-compute for unique genes
        all_genes = set()
        for genes_list in hits_df["gene_symbols"]:
            if hasattr(genes_list, "tolist"):
                genes_list = genes_list.tolist()
            if isinstance(genes_list, (list, tuple)):
                all_genes.update(genes_list)
        for gene in all_genes:
            expr_data = get_gene_expression(gene, expr_df)
            if expr_data:
                expr_lookup[gene] = expr_data

    # Build DecoyAHit objects
    hits: List[DecoyAHit] = []
    for _, row in hits_df.iterrows():
        seq = row["sequence"]
        hd = int(row["hamming_distance"])

        mismatches = annotate_mismatches(target_sequence, seq, hla_allele)
        n_tcr = sum(1 for m in mismatches if m.is_tcr_contact)
        n_anchor = sum(1 for m in mismatches if m.is_anchor)

        # Gene expression — Parquet stores lists as numpy arrays
        gene_symbols = row.get("gene_symbols", [])
        if hasattr(gene_symbols, "tolist"):
            gene_symbols = gene_symbols.tolist()
        elif not isinstance(gene_symbols, list):
            gene_symbols = list(gene_symbols) if gene_symbols is not None else []
        source_proteins = row.get("source_proteins", [])
        if hasattr(source_proteins, "tolist"):
            source_proteins = source_proteins.tolist()
        elif not isinstance(source_proteins, list):
            source_proteins = list(source_proteins) if source_proteins is not None else []

        expression = None
        if gene_symbols and expr_lookup:
            # Use the highest-risk gene for this peptide
            best_expr = None
            for gene in gene_symbols:
                if gene in expr_lookup:
                    expr_data = expr_lookup[gene]
                    if best_expr is None or expr_data["max_vital_organ_tpm"] > best_expr["max_vital_organ_tpm"]:
                        best_expr = expr_data
            if best_expr:
                expression = TissueExpression(**best_expr)

        # Similarity score: inverse of Hamming distance, weighted by position
        # 1-mismatch = 0.9, 2-mismatch = 0.8 (for 9-mer)
        seq_similarity = 1.0 - (hd / target_len)

        hit = DecoyAHit(
            sequence=seq,
            target_sequence=target_sequence,
            hamming_distance=hd,
            mismatches=mismatches,
            n_tcr_contact_mismatches=n_tcr,
            n_anchor_mismatches=n_anchor,
            el_rank=float(row.get("el_rank", 99.0)),
            hla_allele=hla_allele,
            gene_symbols=gene_symbols,
            source_proteins=source_proteins,
            expression=expression,
            similarity_score=seq_similarity,
        )
        hits.append(hit)

    # Sort: fewer mismatches first, then by expression risk (high TPM = high risk)
    hits.sort(key=lambda h: (
        h.hamming_distance,
        -(h.expression.max_vital_organ_tpm if h.expression else 0.0),
        h.el_rank,
    ))

    log.info(
        "Decoy A complete: %d hits (1-mismatch: %d, 2-mismatch: %d)",
        len(hits),
        sum(1 for h in hits if h.hamming_distance == 1),
        sum(1 for h in hits if h.hamming_distance == 2),
    )

    return hits
