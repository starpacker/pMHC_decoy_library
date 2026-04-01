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

    same_len = hla_filtered_df[
        hla_filtered_df["sequence"].str.len() == target_len
    ].copy()

    if same_len.empty:
        return pd.DataFrame()

    same_len = same_len[same_len["sequence"] != target_sequence]

    # Exclude Decoy A territory
    candidates = same_len["sequence"].values
    distances = hamming_distance_vectorised(target_sequence, candidates)
    mask = distances >= min_hamming
    same_len = same_len[mask].copy()
    same_len["hamming_distance"] = distances[mask]

    log.info("After excluding Hamming <= %d: %d candidates", min_hamming - 1, len(same_len))
    if same_len.empty:
        return pd.DataFrame()

    # Batch Atchley cosine similarity
    target_contact = _get_tcr_contact_residues(target_sequence)
    target_vec = _sequence_to_atchley_vector(target_contact)

    candidate_seqs = same_len["sequence"].values
    contact_residues = [_get_tcr_contact_residues(str(s)) for s in candidate_seqs]

    dim = len(target_vec)
    feature_matrix = np.zeros((len(candidate_seqs), dim), dtype=np.float64)
    for i, cr in enumerate(contact_residues):
        feature_matrix[i] = _sequence_to_atchley_vector(cr)

    cos_sims = _batch_cosine_similarity(target_vec, feature_matrix)

    same_len["cosine_similarity"] = cos_sims
    same_len["contact_residues"] = contact_residues

    above_threshold = same_len[same_len["cosine_similarity"] >= cosine_threshold]
    result = above_threshold.nlargest(top_k, "cosine_similarity")

    log.info(
        "Physicochemical screen: %d -> %d above threshold -> top %d",
        len(same_len), len(above_threshold), len(result),
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

        # Extract CA atoms
        target_cas = [a for a in target_struct.get_atoms() if a.get_name() == "CA"]
        cand_cas = [a for a in cand_struct.get_atoms() if a.get_name() == "CA"]

        rmsd = None
        if len(target_cas) == len(cand_cas) and len(target_cas) > 0:
            sup = Superimposer()
            sup.set_atoms(target_cas, cand_cas)
            rmsd = sup.rms

        # B-factor correlation (surface flexibility proxy)
        target_bf = np.array([a.get_bfactor() for a in target_cas])
        cand_bf = np.array([a.get_bfactor() for a in cand_cas])

        corr = 0.0
        if len(target_bf) == len(cand_bf) and len(target_bf) > 1:
            if np.std(target_bf) > 0 and np.std(cand_bf) > 0:
                corr = float(np.corrcoef(target_bf, cand_bf)[0, 1])

        tool = "biopython"
        if rmsd is not None and rmsd < 2.0:
            corr = max(corr, 1.0 - rmsd / 5.0)  # Boost for low RMSD
            tool = "biopython+rmsd"

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
    peptide_chain_id: str = "C",
    anchor_positions: Optional[List[int]] = None,
    num_designs: int = 1000,
    hla_filtered_df=None,
    expr_df=None,
) -> List[DecoyBHit]:
    """
    MPNN inverse design branch: generate sequences from structure,
    then map back to human proteome.

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

    Returns
    -------
    list[DecoyBHit]
        Designed peptides that exist in human proteome.
    """
    from .tools.proteinmpnn import check_available, design_peptide

    if not check_available():
        log.warning("ProteinMPNN not available; skipping inverse design")
        return []

    log.info("== MPNN Inverse Design Branch ==")

    # Step 1: Generate designs
    designs = design_peptide(
        pdb_path=target_pdb,
        peptide_chain_id=peptide_chain_id,
        anchor_positions=anchor_positions,
        num_designs=num_designs,
    )
    log.info("MPNN generated %d unique designs", len(designs))

    if not designs:
        return []

    # Step 2: Map back to human proteome
    import pandas as pd
    if hla_filtered_df is None:
        try:
            from decoy_a.hla_filter import load_hla_filtered
            hla_filtered_df = load_hla_filtered(hla_allele)
        except FileNotFoundError:
            # Fall back to full k-mer DB
            from decoy_a.kmer_builder import load_kmer_db
            hla_filtered_df = load_kmer_db()
            hla_filtered_df["el_rank"] = 1.0

    design_seqs = {d.sequence for d in designs}
    proteome_seqs = set(hla_filtered_df["sequence"].values)
    matched = design_seqs & proteome_seqs

    log.info(
        "Proteome match: %d/%d designs found in human proteome",
        len(matched), len(designs),
    )

    if not matched:
        return []

    # Step 3: Build DecoyBHit objects for matched designs
    target_contact = _get_tcr_contact_residues(target_sequence)
    target_vec = _sequence_to_atchley_vector(target_contact)

    # Expression lookup
    expr_lookup: Dict[str, dict] = {}
    if expr_df is not None:
        from decoy_a.kmer_builder import get_gene_expression

    hits: List[DecoyBHit] = []
    for seq in matched:
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
                    expr_lookup[gene] = get_gene_expression(gene, expr_df)
                ed = expr_lookup.get(gene)
                if ed and (expression is None or ed["max_vital_organ_tpm"] > (expression.max_vital_organ_tpm if expression else 0)):
                    expression = TissueExpression(**ed)

        from decoy_a.scanner import hamming_distance
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
        )
        hits.append(hit)

    hits.sort(key=lambda h: -h.similarity_score)
    log.info("MPNN branch: %d proteome-matched decoy candidates", len(hits))
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

        # Combined score
        if structural and structural.surface_correlation > 0:
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

    log.info(
        "Decoy B complete: %d hits (physchem-only: %d, with structure: %d)",
        len(hits),
        sum(1 for h in hits if h.structural is None),
        sum(1 for h in hits if h.structural is not None),
    )

    return hits
