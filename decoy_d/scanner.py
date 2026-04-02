import logging
from pathlib import Path
from typing import List, Optional, Dict
import pandas as pd

from decoy_a.config import B_DATA_DIR, DEFAULT_HLA_ALLELE
from decoy_a.tools.mhcflurry import check_available as mhcf_available, predict_binding as mhcf_predict
from .tools.proteinmpnn import design_peptide
from .tools.tfold import predict_pmhc_batch

log = logging.getLogger(__name__)

def run_decoy_d(
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    target_pdb: Optional[str] = None,
    num_designs: int = 1000,
    el_rank_threshold: float = 2.0,
    top_k_structures: int = 10,
    output_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Decoy D Pipeline: MPNN + mhcflurry
    1. tFold predict target structure (if not provided).
    2. ProteinMPNN generates designs fixing anchors (p2/p9).
    3. mhcflurry filters designs for HLA presentation.
    4. Top-K designs get 3D structures via tFold.
    """
    if output_dir is None:
        output_dir = Path("data/decoy_d")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Target Structure
    if not target_pdb:
        log.info(f"Predicting target structure for {target_sequence}...")
        res = predict_pmhc_batch([target_sequence], hla_allele, output_dir / "target_pdb")
        target_pdb = res[0].pdb_path
        if not target_pdb:
            raise RuntimeError("Failed to predict target structure with tFold")
    
    # Step 2: ProteinMPNN Design
    log.info(f"Running ProteinMPNN on {target_pdb} for {num_designs} designs...")
    
    # Fix anchor positions based on HLA allele.
    # For HLA-A*02:01: p2 and pΩ (last position) are the canonical anchors.
    # Position numbering is 1-indexed for ProteinMPNN.
    HLA_ANCHOR_MAP = {
        "HLA-A*02:01": {8: [2, 8], 9: [2, 9], 10: [2, 10], 11: [2, 11]},
        "HLA-A*01:01": {9: [2, 9], 10: [2, 10]},
        "HLA-B*07:02": {9: [2, 9], 10: [2, 10]},
    }
    n = len(target_sequence)
    allele_anchors = HLA_ANCHOR_MAP.get(hla_allele, {})
    anchor_positions = allele_anchors.get(n, [2, n])  # fallback: p2 and pΩ
    
    designs = design_peptide(
        pdb_path=target_pdb,
        peptide_chain_id="P",  # tfold uses 'P' for peptide
        anchor_positions=anchor_positions,
        num_designs=num_designs,
    )
    
    log.info(f"MPNN generated {len(designs)} unique designs.")
    if not designs:
        return pd.DataFrame()
        
    design_df = pd.DataFrame([{"sequence": d.sequence, "mpnn_score": d.score} for d in designs])
    
    # Step 3: mhcflurry Filtering
    if not mhcf_available():
        raise RuntimeError("mhcflurry is not available")
        
    log.info("Running mhcflurry presentation filter...")
    candidates = design_df["sequence"].tolist()
    mhcf_results = mhcf_predict(candidates, hla_allele)
    
    mhcf_df = pd.DataFrame([
        {"sequence": r.sequence, "el_rank": r.el_rank, "presentation_score": r.presentation_score} 
        for r in mhcf_results
    ])
    
    df = pd.merge(design_df, mhcf_df, on="sequence")
    
    # Filter by el_rank <= threshold
    filtered_df = df[df["el_rank"] <= el_rank_threshold].copy()
    log.info(f"Filtered {len(filtered_df)} / {len(df)} designs with EL%Rank <= {el_rank_threshold}")
    
    if filtered_df.empty:
        return filtered_df
        
    # Sort by MPNN score (lower is better)
    filtered_df = filtered_df.sort_values("mpnn_score").reset_index(drop=True)
    
    # Step 4: Predict 3D Structures for Top-K
    top_k = filtered_df.head(top_k_structures)
    top_k_seqs = top_k["sequence"].tolist()
    
    log.info(f"Predicting 3D structures for top {len(top_k_seqs)} candidates...")
    tfold_res = predict_pmhc_batch(top_k_seqs, hla_allele, output_dir / "tfold_pdbs")

    # Only include successful predictions (pdb_path is not None)
    pdb_paths = {r.peptide: r.pdb_path for r in tfold_res if r.pdb_path is not None}
    n_failed = sum(1 for r in tfold_res if r.pdb_path is None)
    if n_failed > 0:
        log.warning(f"{n_failed}/{len(tfold_res)} tFold predictions failed — excluded from results")
    filtered_df["pdb_path"] = filtered_df["sequence"].map(pdb_paths)
    
    # Target pdb path
    filtered_df["target_pdb"] = target_pdb
    
    out_csv = output_dir / "decoy_d_results.csv"
    filtered_df.to_csv(out_csv, index=False)
    log.info(f"Decoy D results saved to {out_csv}")
    
    return filtered_df
