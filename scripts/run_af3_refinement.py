#!/usr/bin/env python3
"""
AF3 Refinement Script
=====================
Re-predict the Top-N candidates from Decoy B using AlphaFold3 for
higher-accuracy structural comparison.

Usage:
    # Refine top 10 from existing Decoy B results
    python scripts/run_af3_refinement.py \
        --target GILGFVFTL \
        --hla "HLA-A*02:01" \
        --top-n 10

    # Refine specific sequences
    python scripts/run_af3_refinement.py \
        --target GILGFVFTL \
        --hla "HLA-A*02:01" \
        --sequences LLVGFVFVV LYLGILFAV RALGFLIGL

    # Use existing tFold results as input, refine top 20
    python scripts/run_af3_refinement.py \
        --target GILGFVFTL \
        --hla "HLA-A*02:01" \
        --ranked-json data/decoy_b/final_ranked_decoys.json \
        --top-n 20

    # Custom output directory
    python scripts/run_af3_refinement.py \
        --target GILGFVFTL \
        --hla "HLA-A*02:01" \
        --top-n 10 \
        --output-dir data/decoy_b/af3_refinement

Prerequisites:
    pip install alphafold3 biopython
    # Model weights: set AF3_WEIGHTS_PATH or place af3.bin.zst in $AF3_MODEL_DIR/
    # Databases: python -m alphafold3.download_databases --db_dir $AF3_DB_DIR
"""

import argparse
import json
import logging
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from decoy_a.config import B_DATA_DIR
from decoy_b.tools.alphafold3 import check_available, predict_pmhc_batch, AF3Result
from decoy_b.scanner import compute_structure_similarity

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
log = logging.getLogger("af3_refinement")


def load_candidates_from_ranked(ranked_path: Path, top_n: int) -> list:
    """Load top-N candidate sequences from a Decoy B ranked JSON."""
    data = json.loads(ranked_path.read_text(encoding="utf-8"))
    seqs = []
    for entry in data[:top_n]:
        seq = entry.get("sequence", "")
        if seq:
            seqs.append(seq)
    return seqs


def main():
    parser = argparse.ArgumentParser(
        description="Refine Decoy B candidates with AlphaFold3",
    )
    parser.add_argument("--target", required=True, help="Target peptide sequence")
    parser.add_argument("--hla", default="HLA-A*02:01", help="HLA allele")
    parser.add_argument("--top-n", type=int, default=10, help="Number of top candidates to refine")
    parser.add_argument("--sequences", nargs="+", help="Specific sequences to refine (overrides --top-n)")
    parser.add_argument("--ranked-json", type=str, help="Path to ranked Decoy B JSON")
    parser.add_argument("--output-dir", type=str, help="Output directory for AF3 structures")
    args = parser.parse_args()

    target = args.target.strip().upper()

    # ── Check AF3 availability ───────────────────────────────────────────
    if not check_available():
        log.error(
            "AlphaFold3 is not available.\n"
            "Install: pip install alphafold3\n"
            "Weights: set AF3_WEIGHTS_PATH to point to af3.bin.zst\n"
            "Databases: python -m alphafold3.download_databases --db_dir <path>"
        )
        sys.exit(1)

    log.info("AF3 is available. Starting refinement.")

    # ── Determine candidates ─────────────────────────────────────────────
    if args.sequences:
        candidates = [s.strip().upper() for s in args.sequences]
    elif args.ranked_json:
        ranked_path = Path(args.ranked_json)
        if not ranked_path.exists():
            log.error("Ranked JSON not found: %s", ranked_path)
            sys.exit(1)
        candidates = load_candidates_from_ranked(ranked_path, args.top_n)
    else:
        # Try default path
        default_ranked = B_DATA_DIR / "final_ranked_decoys.json"
        if default_ranked.exists():
            candidates = load_candidates_from_ranked(default_ranked, args.top_n)
        else:
            log.error(
                "No candidates specified. Use --sequences, --ranked-json, "
                "or ensure data/decoy_b/final_ranked_decoys.json exists."
            )
            sys.exit(1)

    log.info("Target: %s (%s)", target, args.hla)
    log.info("Candidates to refine: %d", len(candidates))
    for i, seq in enumerate(candidates):
        log.info("  %2d. %s", i + 1, seq)

    # ── Output directory ─────────────────────────────────────────────────
    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = B_DATA_DIR / "pmhc_models" / "af3"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Predict target structure ─────────────────────────────────────────
    log.info("=" * 60)
    log.info("Step 1: Predicting target structure with AF3...")
    all_peptides = [target] + candidates
    results = predict_pmhc_batch(all_peptides, args.hla, output_dir)

    target_result = results[0]
    candidate_results = results[1:]

    if not target_result.success:
        log.error("Failed to predict target structure: %s", target_result.error)
        sys.exit(1)

    log.info("Target structure: %s (pLDDT=%.1f, ipTM=%.3f)",
             target_result.pdb_path, target_result.plddt, target_result.iptm)

    # ── Structure comparison ─────────────────────────────────────────────
    log.info("=" * 60)
    log.info("Step 2: Computing structural similarity (dual superposition)...")

    comparison_results = []
    for cand_result in candidate_results:
        if not cand_result.success:
            log.warning("  SKIP %s: %s", cand_result.peptide, cand_result.error)
            comparison_results.append({
                "sequence": cand_result.peptide,
                "af3_success": False,
                "error": cand_result.error,
            })
            continue

        struct_score = compute_structure_similarity(
            target_result.pdb_path, cand_result.pdb_path,
        )

        entry = {
            "sequence": cand_result.peptide,
            "af3_success": True,
            "af3_pdb": cand_result.pdb_path,
            "af3_plddt": cand_result.plddt,
            "af3_iptm": cand_result.iptm,
            "af3_ranking_score": cand_result.ranking_score,
            "structural_similarity": struct_score.surface_correlation,
            "rmsd": struct_score.rmsd,
            "modeling_tool": struct_score.modeling_tool,
        }
        comparison_results.append(entry)

        log.info("  %s: similarity=%.4f, RMSD=%.3f, pLDDT=%.1f, ipTM=%.3f",
                 cand_result.peptide,
                 struct_score.surface_correlation,
                 struct_score.rmsd or 0,
                 cand_result.plddt,
                 cand_result.iptm)

    # ── Save results ─────────────────────────────────────────────────────
    output_json = output_dir / f"af3_refinement_{target}.json"
    report = {
        "target": target,
        "hla_allele": args.hla,
        "target_pdb": target_result.pdb_path,
        "target_plddt": target_result.plddt,
        "target_iptm": target_result.iptm,
        "candidates": sorted(
            comparison_results,
            key=lambda x: x.get("structural_similarity", 0),
            reverse=True,
        ),
    }
    output_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    # ── Summary ──────────────────────────────────────────────────────────
    n_ok = sum(1 for r in comparison_results if r.get("af3_success"))
    log.info("=" * 60)
    log.info("AF3 Refinement Complete")
    log.info("  Succeeded: %d/%d", n_ok, len(comparison_results))
    log.info("  Results: %s", output_json)
    if n_ok > 0:
        best = max(comparison_results, key=lambda x: x.get("structural_similarity", 0))
        log.info("  Best match: %s (similarity=%.4f, RMSD=%.3f)",
                 best["sequence"], best["structural_similarity"], best.get("rmsd", 0))


if __name__ == "__main__":
    main()
