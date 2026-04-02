"""
Step 3B: Decoy B — Structural & Physicochemical Similarity Scan
================================================================
Identifies off-target peptides that are **sequentially distant** (Hamming > 2)
but share similar 3D surface conformation and electrostatic charge with the
target peptide — the classic "Titin-like" scenario.

Four-stage approach:
    1. **Atchley factor screening** — Fast physicochemical similarity on
       TCR-contact residues → top 2,000-5,000 candidates.
    2. **tFold bulk prediction** — Batch pMHC structure prediction for
       top candidates (10-30s per peptide, GPU).
    3. **AF3 refinement** — High-accuracy prediction for top 200
       (2-5 min per peptide).
    4. **Structure comparison** — RMSD / TM-score / surface electrostatics.

Optional parallel branch:
    **MPNN inverse design** — Generate sequences that preserve TCR-facing
    surface, then map back to human proteome.

Public API
----------
    scan_decoy_b(target_seq, hla_allele, ...) -> list[DecoyBHit]
    run_mpnn_design(target_pdb, ...) -> list[DecoyBHit]
    compute_physicochemical_features(sequence) -> PhysicochemFeatures
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from decoy_a.config import (
    AB_DATA_DIR,
    B_DATA_DIR,
    ATCHLEY_FACTORS,
    DEFAULT_HLA_ALLELE,
    HAMMING_DISTANCE_MAX,
    PHYSICOCHEMICAL_COSINE_THRESHOLD,
    PHYSICOCHEMICAL_TOP_K,
    SURFACE_SIMILARITY_THRESHOLD,
    VITAL_ORGANS,
)
from decoy_a.models import (
    DecoyBHit,
    PhysicochemFeatures,
    StructuralScore,
    TissueExpression,
)

log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
#  Stage 1: Physicochemical Feature Screening (Atchley Factors)
# ═══════════════════════════════════════════════════════════════════════════

def _get_tcr_contact_residues(sequence: str, peptide_length: int = 0) -> str:
    """
    Extract the TCR-contact residues from a peptide.
    For a 9-mer: positions 4-8 (indices 3-7) = central 5 residues.
    """
    n = len(sequence)
    if n <= 5:
        return sequence
    n_contact = min(5, n - 2)
    start = (n - n_contact) // 2
    return sequence[start : start + n_contact]


def _sequence_to_atchley_vector(sequence: str) -> np.ndarray:
    """Convert amino acid sequence to flattened Atchley factor vector."""
    factors = []
    for aa in sequence:
        if aa in ATCHLEY_FACTORS:
            factors.extend(ATCHLEY_FACTORS[aa])
        else:
            factors.extend([0.0] * 5)
    return np.array(factors, dtype=np.float64)


def compute_physicochemical_features(
    sequence: str,
    target_vector: Optional[np.ndarray] = None,
) -> PhysicochemFeatures:
    """Compute Atchley-based features for a peptide's TCR-contact region."""
    contact = _get_tcr_contact_residues(sequence)
    vec = _sequence_to_atchley_vector(contact)

    cos_sim = 0.0
    if target_vector is not None:
        dot = np.dot(vec, target_vector)
        norm_a = np.linalg.norm(vec)
        norm_b = np.linalg.norm(target_vector)
        if norm_a > 0 and norm_b > 0:
            cos_sim = float(dot / (norm_a * norm_b))

    return PhysicochemFeatures(
        contact_residues=contact,
        feature_vector=vec.tolist(),
        cosine_similarity=cos_sim,
    )


def _batch_cosine_similarity(
    target_vec: np.ndarray,
    candidate_vecs: np.ndarray,
) -> np.ndarray:
    """Vectorised cosine similarity between target and many candidates."""
    target_norm = target_vec / (np.linalg.norm(target_vec) + 1e-12)
    norms = np.linalg.norm(candidate_vecs, axis=1, keepdims=True) + 1e-12
    cand_normed = candidate_vecs / norms
    return cand_normed @ target_norm


def physicochemical_screen(
    target_sequence: str,
    hla_filtered_df,
    min_hamming: int = HAMMING_DISTANCE_MAX + 1,
    cosine_threshold: float = PHYSICOCHEMICAL_COSINE_THRESHOLD,
    top_k: int = PHYSICOCHEMICAL_TOP_K,
) -> "pd.DataFrame":
    """
    Stage 1: Fast Atchley-factor cosine similarity screening.
    Filters out Decoy-A candidates (Hamming <= threshold).
    """
    import pandas as pd
    from decoy_a.scanner import hamming_distance_vectorised

    target_sequence = target_sequence.strip().upper()
    target_len = len(target_sequence)

    log.info(
        "Decoy B physicochemical screen: target=%s, min_hamming=%d, threshold=%.2f",
        target_sequence, min_hamming, cosine_threshold,
    )

    # Accept peptides of length 8-11 (HLA-I binding range), not just same-length.
    # Atchley TCR-contact core is always 5 residues regardless of peptide length,
    # so cross-length cosine similarity is well-defined.
    ALLOWED_LENGTHS = {8, 9, 10, 11}
    candidates_df = hla_filtered_df[
        hla_filtered_df["sequence"].str.len().isin(ALLOWED_LENGTHS)
    ].copy()

    if candidates_df.empty:
        return pd.DataFrame()

    candidates_df = candidates_df[candidates_df["sequence"] != target_sequence]

    # Exclude Decoy A territory (Hamming distance only defined for same-length)
    # For different-length peptides, they are automatically outside Decoy A scope
    cand_seqs = candidates_df["sequence"].values
    cand_lens = np.array([len(s) for s in cand_seqs])

    # Same-length candidates: apply Hamming filter
    same_mask = cand_lens == target_len
    diff_mask = ~same_mask

    hamming_dists = np.full(len(cand_seqs), -1, dtype=np.int32)
    if same_mask.any():
        hamming_dists[same_mask] = hamming_distance_vectorised(
            target_sequence, cand_seqs[same_mask]
        )

    # Keep: same-length with HD >= min_hamming, OR different-length (always HD > length)
    keep_mask = (same_mask & (hamming_dists >= min_hamming)) | diff_mask
    candidates_df = candidates_df[keep_mask].copy()
    candidates_df["hamming_distance"] = hamming_dists[keep_mask]
    # For different-length peptides, set hamming_distance = peptide length (max dissimilarity)
    diff_len_mask = candidates_df["hamming_distance"] < 0
    candidates_df["hamming_distance"] = candidates_df["hamming_distance"].astype(np.int64)
    candidates_df.loc[diff_len_mask, "hamming_distance"] = \
        candidates_df.loc[diff_len_mask, "sequence"].str.len().astype(np.int64)

    log.info(
        "After excluding Decoy A territory: %d candidates "
        "(same-len %d-mer: %d, other lengths: %d)",
        len(candidates_df),
        target_len,
        int((candidates_df["sequence"].str.len() == target_len).sum()),
        int((candidates_df["sequence"].str.len() != target_len).sum()),
    )
    if candidates_df.empty:
        return pd.DataFrame()

    # Batch Atchley cosine similarity
    target_contact = _get_tcr_contact_residues(target_sequence)
    target_vec = _sequence_to_atchley_vector(target_contact)

    candidate_seqs = candidates_df["sequence"].values
    contact_residues = [_get_tcr_contact_residues(str(s)) for s in candidate_seqs]

    dim = len(target_vec)
    feature_matrix = np.zeros((len(candidate_seqs), dim), dtype=np.float64)
    for i, cr in enumerate(contact_residues):
        feature_matrix[i] = _sequence_to_atchley_vector(cr)

    cos_sims = _batch_cosine_similarity(target_vec, feature_matrix)

    candidates_df["cosine_similarity"] = cos_sims
    candidates_df["contact_residues"] = contact_residues

    above_threshold = candidates_df[candidates_df["cosine_similarity"] >= cosine_threshold]
    result = above_threshold.nlargest(top_k, "cosine_similarity")

    log.info(
        "Physicochemical screen: %d -> %d above threshold -> top %d",
        len(candidates_df), len(above_threshold), len(result),
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Stage 2: Structure Prediction (tFold bulk + AF3 refinement)
# ═══════════════════════════════════════════════════════════════════════════

def _predict_structures_tfold(
    peptides: List[str],
    hla_allele: str,
    output_dir: Optional[Path] = None,
) -> Dict[str, Optional[str]]:
    """
    Predict pMHC structures using tFold (bulk stage).
    Returns dict: peptide -> pdb_path (or None if failed).
    """
    from .tools.tfold import check_available, predict_pmhc_batch

    if not check_available():
        log.info("tFold not available; skipping bulk structure prediction")
        return {pep: None for pep in peptides}

    if output_dir is None:
        output_dir = B_DATA_DIR / "pmhc_models" / "tfold"

    results = predict_pmhc_batch(peptides, hla_allele, output_dir)
    return {r.peptide: r.pdb_path for r in results}


def _predict_structures_af3(
    peptides: List[str],
    hla_allele: str,
    output_dir: Optional[Path] = None,
) -> Dict[str, Optional[str]]:
    """
    Predict pMHC structures using AF3 (refinement stage).
    Returns dict: peptide -> pdb_path (or None).
    """
    from .tools.alphafold3 import check_available, predict_pmhc_batch

    if not check_available():
        log.info("AlphaFold3 not available; skipping refinement prediction")
        return {pep: None for pep in peptides}

    if output_dir is None:
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"

    results = predict_pmhc_batch(peptides, hla_allele, output_dir)

    pdb_map = {}
    for r in results:
        # Prefer PDB over CIF
        pdb_map[r.peptide] = r.pdb_path or r.cif_path
    return pdb_map


# ═══════════════════════════════════════════════════════════════════════════
#  Stage 3: Structure Comparison
# ═══════════════════════════════════════════════════════════════════════════

def compute_structure_similarity(
    target_pdb: Optional[str],
    candidate_pdb: Optional[str],
) -> StructuralScore:
    """
    Compare two pMHC structures.

    Computes:
    - Peptide backbone RMSD (after superposition on MHC)
    - B-factor correlation (proxy for surface flexibility)
    - If TMalign available: TM-score

    Parameters
    ----------
    target_pdb, candidate_pdb : str or None
        Paths to PDB files.

    Returns
    -------
    StructuralScore
    """
    if target_pdb is None or candidate_pdb is None:
        return StructuralScore(modeling_tool="none", surface_correlation=0.0)

    # Try BioPython structural comparison
    try:
        from Bio.PDB import PDBParser, Superimposer
        parser = PDBParser(QUIET=True)

        target_struct = parser.get_structure("target", str(target_pdb))
        cand_struct = parser.get_structure("candidate", str(candidate_pdb))

        # ── Separate MHC (chains M+N) and peptide (chain P) CA atoms ──
        # tFold outputs: chain M = MHC heavy, N = β2m, P = peptide
        PEPTIDE_CHAIN = "P"
        MHC_CHAINS = {"M", "N"}

        def _get_cas_by_chains(struct, chain_ids):
            """Extract CA atoms from specified chains."""
            cas = []
            for model in struct:
                for chain in model:
                    if chain.id in chain_ids:
                        for residue in chain:
                            for atom in residue:
                                if atom.get_name() == "CA":
                                    cas.append(atom)
            return cas

        target_mhc_cas = _get_cas_by_chains(target_struct, MHC_CHAINS)
        cand_mhc_cas = _get_cas_by_chains(cand_struct, MHC_CHAINS)
        target_pep_cas = _get_cas_by_chains(target_struct, {PEPTIDE_CHAIN})
        cand_pep_cas = _get_cas_by_chains(cand_struct, {PEPTIDE_CHAIN})

        peptide_rmsd = None
        whole_rmsd = None

        # Step 1: Superimpose on MHC backbone (align the HLA groove)
        if (len(target_mhc_cas) == len(cand_mhc_cas) and len(target_mhc_cas) > 0):
            sup_mhc = Superimposer()
            sup_mhc.set_atoms(target_mhc_cas, cand_mhc_cas)
            whole_rmsd = sup_mhc.rms

            # Apply the MHC-derived rotation to ALL candidate atoms
            sup_mhc.apply(list(cand_struct.get_atoms()))

            # Step 2: Peptide RMSD (no further fitting, just measure after MHC alignment)
            n_tgt = len(target_pep_cas)
            n_cand = len(cand_pep_cas)

            if n_tgt > 0 and n_cand > 0:
                if n_tgt == n_cand:
                    # Same length: full peptide RMSD
                    target_pep_coords = np.array([a.get_vector().get_array() for a in target_pep_cas])
                    cand_pep_coords = np.array([a.get_vector().get_array() for a in cand_pep_cas])
                    diff = target_pep_coords - cand_pep_coords
                    peptide_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
                else:
                    # Cross-length: compare TCR-contact core (central 5 residues)
                    # For a peptide of length n, the TCR-contact core is the central
                    # min(5, n-2) residues. We extract these CA atoms and compute RMSD.
                    def _core_cas(cas_list):
                        n = len(cas_list)
                        n_contact = min(5, n - 2) if n > 5 else n
                        start = (n - n_contact) // 2
                        return cas_list[start : start + n_contact]

                    tgt_core = _core_cas(target_pep_cas)
                    cand_core = _core_cas(cand_pep_cas)

                    if len(tgt_core) == len(cand_core) and len(tgt_core) > 0:
                        tgt_coords = np.array([a.get_vector().get_array() for a in tgt_core])
                        cand_coords = np.array([a.get_vector().get_array() for a in cand_core])
                        diff = tgt_coords - cand_coords
                        peptide_rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

        # B-factor correlation on peptide chain (surface flexibility proxy)
        # For cross-length, use the shorter set for correlation
        n_bf = min(len(target_pep_cas), len(cand_pep_cas))
        target_pep_bf = np.array([a.get_bfactor() for a in target_pep_cas[:n_bf]])
        cand_pep_bf = np.array([a.get_bfactor() for a in cand_pep_cas[:n_bf]])

        corr = 0.0
        if len(target_pep_bf) == len(cand_pep_bf) and len(target_pep_bf) > 1:
            if np.std(target_pep_bf) > 0 and np.std(cand_pep_bf) > 0:
                corr = float(np.corrcoef(target_pep_bf, cand_pep_bf)[0, 1])

        tool = "biopython"
        # Use peptide RMSD for scoring (this is what matters for TCR cross-reactivity)
        rmsd = peptide_rmsd if peptide_rmsd is not None else whole_rmsd
        if rmsd is not None and rmsd < 3.0:
            # Convert RMSD to a 0-1 similarity: RMSD=0 → 1.0, RMSD=3 → 0.0
            rmsd_similarity = max(0.0, 1.0 - rmsd / 3.0)
            corr = max(corr, rmsd_similarity)
            n_tgt = len(target_pep_cas)
            n_cand = len(cand_pep_cas)
            if n_tgt == n_cand:
                tool = "biopython+peptide_rmsd"
            else:
                tool = f"biopython+core_rmsd({n_tgt}v{n_cand})"

        return StructuralScore(
            modeling_tool=tool,
            pdb_path=str(candidate_pdb),
            surface_correlation=max(0.0, corr),
            rmsd=rmsd,
        )

    except ImportError:
        log.debug("BioPython not available for structural comparison")
        return StructuralScore(
            modeling_tool="unavailable",
            pdb_path=str(candidate_pdb),
            surface_correlation=0.0,
        )
    except Exception as exc:
        log.warning("Structural comparison failed: %s", exc)
        return StructuralScore(modeling_tool="error", surface_correlation=0.0)


# ═══════════════════════════════════════════════════════════════════════════
#  MPNN Inverse Design Branch
# ═══════════════════════════════════════════════════════════════════════════

def run_mpnn_design(
    target_pdb: str,
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    peptide_chain_id: str = "P",
    anchor_positions: Optional[List[int]] = None,
    num_designs: int = 1000,
    hla_filtered_df=None,
    expr_df=None,
    el_rank_threshold: float = 2.0,
) -> List[DecoyBHit]:
    """
    MPNN inverse design branch — enhanced two-tier qualification pipeline.

    Tier 1 (Proteome-matched):
        MPNN designs that exactly match a human proteome peptide already
        known to be HLA-presentable.  Highest confidence.

    Tier 2 (HLA-qualified synthetic):
        Designs that do NOT exist in the proteome, but are predicted by
        mhcflurry to be presentable (presentation_percentile ≤ threshold).
        These represent *theoretical* cross-reactivity risks — sequences
        that could arise from somatic mutations, non-canonical ORFs, or
        proteins absent from Swiss-Prot.

    Pipeline
    --------
    1. ProteinMPNN generates ~N designs from the target pMHC backbone,
       fixing HLA anchor positions (p2/p9) and redesigning TCR-contact
       residues (p4-p8).
    2. Each design is checked against the HLA-filtered human proteome.
       → Exact match  → Tier 1 hit (gene/expression annotated).
    3. Non-matched designs are batch-predicted by mhcflurry.
       → presentation_percentile ≤ threshold → Tier 2 hit.
    4. Both tiers are merged, scored, and returned.

    Parameters
    ----------
    target_pdb : str
        Path to target pMHC structure (from tFold/AF3).
    target_sequence : str
        Target peptide sequence.
    hla_allele : str
        Target HLA allele.
    peptide_chain_id : str
        Chain ID of peptide in PDB.
    anchor_positions : list[int], optional
        1-indexed HLA anchor positions to fix.
    num_designs : int
        Number of designs to generate.
    hla_filtered_df : DataFrame, optional
        HLA-filtered k-mer pool for proteome matching.
    expr_df : DataFrame, optional
        Expression database.
    el_rank_threshold : float
        Max presentation_percentile for Tier 2 qualification (default 2.0).

    Returns
    -------
    list[DecoyBHit]
        Qualified designs from both tiers, sorted by similarity_score.
    """
    from .tools.proteinmpnn import check_available, design_peptide

    if not check_available():
        log.warning("ProteinMPNN not available; skipping inverse design")
        return []

    log.info("== MPNN Inverse Design Branch (Enhanced Two-Tier) ==")

    # ── Step 1: Generate designs ────────────────────────────────────────
    designs = design_peptide(
        pdb_path=target_pdb,
        peptide_chain_id=peptide_chain_id,
        anchor_positions=anchor_positions,
        num_designs=num_designs,
    )
    log.info("MPNN generated %d unique designs", len(designs))

    if not designs:
        return []

    # Build a lookup: sequence → MPNNDesign (for score access)
    design_lookup = {d.sequence: d for d in designs}

    # ── Step 2: Load proteome data ──────────────────────────────────────
    import pandas as pd
    if hla_filtered_df is None:
        try:
            from decoy_a.hla_filter import load_hla_filtered
            hla_filtered_df = load_hla_filtered(hla_allele)
        except FileNotFoundError:
            from decoy_a.kmer_builder import load_kmer_db
            hla_filtered_df = load_kmer_db()
            hla_filtered_df["el_rank"] = 1.0

    design_seqs = set(design_lookup.keys())
    proteome_seqs = set(hla_filtered_df["sequence"].values)

    # ── Tier 1: Exact proteome matches ──────────────────────────────────
    tier1_seqs = design_seqs & proteome_seqs
    tier2_candidates = design_seqs - proteome_seqs

    log.info(
        "Proteome match: %d/%d designs found in human proteome (Tier 1)",
        len(tier1_seqs), len(designs),
    )
    log.info(
        "Non-matched designs for HLA qualification: %d (Tier 2 candidates)",
        len(tier2_candidates),
    )

    # ── Tier 2: mhcflurry HLA presentation filter ──────────────────────
    tier2_seqs: Dict[str, float] = {}  # sequence → el_rank
    if tier2_candidates:
        try:
            from decoy_a.tools.mhcflurry import (
                check_available as mhcf_available,
                predict_binding as mhcf_predict,
            )

            if mhcf_available():
                candidate_list = sorted(tier2_candidates)
                log.info(
                    "Running mhcflurry HLA presentation filter on %d "
                    "non-proteome designs (threshold: EL%%Rank ≤ %.1f)...",
                    len(candidate_list), el_rank_threshold,
                )
                results = mhcf_predict(candidate_list, hla_allele)
                for r in results:
                    if r.el_rank <= el_rank_threshold:
                        tier2_seqs[r.sequence] = r.el_rank

                log.info(
                    "mhcflurry Tier 2 filter: %d/%d designs predicted as "
                    "HLA-presentable (EL%%Rank ≤ %.1f)",
                    len(tier2_seqs), len(candidate_list), el_rank_threshold,
                )
            else:
                log.warning("mhcflurry not available; Tier 2 skipped")
        except ImportError:
            log.warning("mhcflurry import failed; Tier 2 skipped")

    # ── Step 3: Build DecoyBHit objects ─────────────────────────────────
    target_contact = _get_tcr_contact_residues(target_sequence)
    target_vec = _sequence_to_atchley_vector(target_contact)

    expr_lookup: Dict[str, dict] = {}
    if expr_df is not None:
        from decoy_a.kmer_builder import get_gene_expression

    from decoy_a.scanner import hamming_distance

    hits: List[DecoyBHit] = []

    # ── Tier 1 hits (proteome-matched) ──────────────────────────────────
    for seq in tier1_seqs:
        row = hla_filtered_df[hla_filtered_df["sequence"] == seq].iloc[0]
        physchem = compute_physicochemical_features(seq, target_vec)

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
        if expr_df is not None and gene_symbols:
            for gene in gene_symbols:
                if gene not in expr_lookup:
                    result = get_gene_expression(gene, expr_df)
                    expr_lookup[gene] = result  # may be None
                ed = expr_lookup.get(gene)
                if ed is not None and (expression is None or ed["max_vital_organ_tpm"] > (expression.max_vital_organ_tpm if expression else 0)):
                    expression = TissueExpression(**ed)

        hd = hamming_distance(target_sequence, seq)

        hit = DecoyBHit(
            sequence=seq,
            target_sequence=target_sequence,
            hamming_distance=hd,
            el_rank=float(row.get("el_rank", 99.0)),
            hla_allele=hla_allele,
            gene_symbols=gene_symbols,
            source_proteins=source_proteins,
            expression=expression,
            physicochemical=physchem,
            structural=None,
            similarity_score=physchem.cosine_similarity,
            mpnn_source="proteome_matched",
        )
        hits.append(hit)

    # ── Tier 2 hits (HLA-qualified synthetic) ───────────────────────────
    for seq, el_rank in tier2_seqs.items():
        physchem = compute_physicochemical_features(seq, target_vec)
        hd = hamming_distance(target_sequence, seq)

        mpnn_design = design_lookup.get(seq)
        mpnn_score = mpnn_design.score if mpnn_design else 0.0

        hit = DecoyBHit(
            sequence=seq,
            target_sequence=target_sequence,
            hamming_distance=hd,
            el_rank=el_rank,
            hla_allele=hla_allele,
            gene_symbols=[],
            source_proteins=[],
            expression=None,
            physicochemical=physchem,
            structural=None,
            similarity_score=physchem.cosine_similarity,
            mpnn_source="hla_qualified_synthetic",
            mpnn_score=mpnn_score,
        )
        hits.append(hit)

    hits.sort(key=lambda h: -h.similarity_score)

    n_t1 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "proteome_matched")
    n_t2 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "hla_qualified_synthetic")
    log.info(
        "MPNN branch complete: %d total hits "
        "(Tier 1 proteome-matched: %d, Tier 2 HLA-qualified: %d)",
        len(hits), n_t1, n_t2,
    )
    return hits


# ═══════════════════════════════════════════════════════════════════════════
#  Main Decoy B Scan (Full Pipeline)
# ═══════════════════════════════════════════════════════════════════════════

def scan_decoy_b(
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    hla_filtered_df=None,
    expr_df=None,
    run_structural: bool = True,
    run_af3_refinement: bool = True,
    run_mpnn: bool = False,
    cosine_threshold: float = PHYSICOCHEMICAL_COSINE_THRESHOLD,
    top_k: int = PHYSICOCHEMICAL_TOP_K,
    af3_top_n: int = 200,
    surface_threshold: float = SURFACE_SIMILARITY_THRESHOLD,
) -> List[DecoyBHit]:
    """
    Full Decoy B pipeline:
        Stage 1: Atchley physicochemical screen → top K
        Stage 2: tFold bulk structure prediction
        Stage 3: AF3 refinement (top N)
        Stage 4: Structure comparison & scoring

    Parameters
    ----------
    target_sequence : str
        Target peptide.
    hla_allele : str
        Target HLA allele.
    hla_filtered_df : DataFrame, optional
        HLA-filtered k-mers.
    expr_df : DataFrame, optional
        Expression database.
    run_structural : bool
        Run structure prediction (tFold).
    run_af3_refinement : bool
        Run AF3 on top candidates.
    run_mpnn : bool
        Run MPNN inverse design branch.
    cosine_threshold : float
        Min Atchley cosine similarity.
    top_k : int
        Candidates for structure prediction.
    af3_top_n : int
        Candidates for AF3 refinement.
    surface_threshold : float
        Min surface similarity for final results.

    Returns
    -------
    list[DecoyBHit]
    """
    import pandas as pd

    target_sequence = target_sequence.strip().upper()
    log.info(
        "Decoy B scan: target=%s, HLA=%s, structural=%s, af3=%s, mpnn=%s",
        target_sequence, hla_allele, run_structural, run_af3_refinement, run_mpnn,
    )

    # Load data if needed
    if hla_filtered_df is None:
        from decoy_a.hla_filter import load_hla_filtered
        hla_filtered_df = load_hla_filtered(hla_allele)

    if expr_df is None:
        try:
            from decoy_a.kmer_builder import load_expression_db
            expr_df = load_expression_db()
        except FileNotFoundError:
            log.warning("Expression database not available")

    # ── Stage 1: Physicochemical Screening ──────────────────────────────
    physchem_df = physicochemical_screen(
        target_sequence, hla_filtered_df,
        cosine_threshold=cosine_threshold, top_k=top_k,
    )

    if physchem_df.empty:
        log.warning("No candidates passed physicochemical screening")
        return []

    log.info("Stage 1 complete: %d candidates", len(physchem_df))

    # ── Stage 2: tFold Bulk Structure Prediction ────────────────────────
    target_pdb_path = None
    candidate_pdbs: Dict[str, Optional[str]] = {}

    if run_structural:
        candidate_peptides = physchem_df["sequence"].tolist()

        # Predict target structure first
        log.info("Stage 2: tFold structure prediction for %d + 1 peptides", len(candidate_peptides))
        all_peptides = [target_sequence] + candidate_peptides
        pdb_map = _predict_structures_tfold(all_peptides, hla_allele)

        target_pdb_path = pdb_map.get(target_sequence)
        candidate_pdbs = {p: pdb_map.get(p) for p in candidate_peptides}

        n_predicted = sum(1 for v in candidate_pdbs.values() if v is not None)
        log.info("Stage 2 complete: %d/%d structures predicted", n_predicted, len(candidate_peptides))

    # ── Stage 3: AF3 Refinement (Top N) ─────────────────────────────────
    if run_af3_refinement and run_structural:
        # Select top candidates by cosine similarity for AF3
        top_candidates = physchem_df.nlargest(af3_top_n, "cosine_similarity")["sequence"].tolist()

        log.info("Stage 3: AF3 refinement for top %d candidates", len(top_candidates))
        af3_pdbs = _predict_structures_af3(top_candidates, hla_allele)

        # AF3 target structure
        if target_pdb_path is None:
            af3_target = _predict_structures_af3([target_sequence], hla_allele)
            target_pdb_path = af3_target.get(target_sequence)

        # Merge: prefer AF3 over tFold
        for pep, pdb in af3_pdbs.items():
            if pdb is not None:
                candidate_pdbs[pep] = pdb

    # ── Stage 4: Build DecoyBHit Objects ────────────────────────────────
    expr_lookup: Dict[str, dict] = {}
    if expr_df is not None:
        from decoy_a.kmer_builder import get_gene_expression
        all_genes = set()
        for genes_list in physchem_df["gene_symbols"]:
            if hasattr(genes_list, "tolist"):
                genes_list = genes_list.tolist()
            if isinstance(genes_list, (list, tuple)):
                all_genes.update(genes_list)
        for gene in all_genes:
            ed = get_gene_expression(gene, expr_df)
            if ed:
                expr_lookup[gene] = ed

    target_contact = _get_tcr_contact_residues(target_sequence)
    target_vec = _sequence_to_atchley_vector(target_contact)

    hits: List[DecoyBHit] = []
    for _, row in physchem_df.iterrows():
        seq = str(row["sequence"])
        cos_sim = float(row["cosine_similarity"])
        hd = int(row.get("hamming_distance", 99))

        physchem = PhysicochemFeatures(
            contact_residues=str(row.get("contact_residues", "")),
            feature_vector=_sequence_to_atchley_vector(
                str(row.get("contact_residues", ""))
            ).tolist(),
            cosine_similarity=cos_sim,
        )

        # Structure comparison
        structural = None
        cand_pdb = candidate_pdbs.get(seq)
        if cand_pdb and target_pdb_path:
            structural = compute_structure_similarity(target_pdb_path, cand_pdb)

        # Expression
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
            best_expr = None
            for gene in gene_symbols:
                if gene in expr_lookup:
                    ed = expr_lookup[gene]
                    if best_expr is None or ed["max_vital_organ_tpm"] > best_expr["max_vital_organ_tpm"]:
                        best_expr = ed
            if best_expr:
                expression = TissueExpression(**best_expr)

        # Combined score: use structural data whenever available, regardless
        # of whether surface_correlation is exactly zero (which can happen
        # from numerical precision, not just missing data)
        if structural is not None:
            combined = 0.4 * cos_sim + 0.6 * structural.surface_correlation
        else:
            combined = cos_sim

        hit = DecoyBHit(
            sequence=seq,
            target_sequence=target_sequence,
            hamming_distance=hd,
            el_rank=float(row.get("el_rank", 99.0)),
            hla_allele=hla_allele,
            gene_symbols=gene_symbols,
            source_proteins=source_proteins,
            expression=expression,
            physicochemical=physchem,
            structural=structural,
            similarity_score=combined,
        )
        hits.append(hit)

    hits.sort(key=lambda h: -h.similarity_score)

    # ── Optional: MPNN Inverse Design Branch ────────────────────────────
    if run_mpnn and target_pdb_path:
        mpnn_hits = run_mpnn_design(
            target_pdb=target_pdb_path,
            target_sequence=target_sequence,
            hla_allele=hla_allele,
            hla_filtered_df=hla_filtered_df,
            expr_df=expr_df,
        )
        # Merge MPNN hits (avoid duplicates)
        existing_seqs = {h.sequence for h in hits}
        for mh in mpnn_hits:
            if mh.sequence not in existing_seqs:
                hits.append(mh)
                existing_seqs.add(mh.sequence)

        hits.sort(key=lambda h: -h.similarity_score)

    n_mpnn_t1 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "proteome_matched")
    n_mpnn_t2 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "hla_qualified_synthetic")
    log.info(
        "Decoy B complete: %d hits (physchem-only: %d, with structure: %d, "
        "MPNN Tier1: %d, MPNN Tier2: %d)",
        len(hits),
        sum(1 for h in hits if h.structural is None and not getattr(h, "mpnn_source", "")),
        sum(1 for h in hits if h.structural is not None),
        n_mpnn_t1, n_mpnn_t2,
    )

    return hits
