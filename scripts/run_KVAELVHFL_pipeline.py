#!/usr/bin/env python3
"""
KVAELVHFL Decoy Library Pipeline
=================================
End-to-end pipeline for peptide KVAELVHFL (HLA-A*02:01).

Generates a summary directory structure identical to GILGFVFTL_summary:
    data/KVAELVHFL_summary/
        Decoy_A/decoy_a_results.json
        Decoy_B/final_ranked_decoys.json
        Decoy_B/3D_structures/  (PDB files)
        Decoy_D/decoy_d_results.csv
        Decoy_D/3D_structures/  (PDB files)
        Decoy_D/target_pdb/     (target PDB)

Key features:
  - Decoy B uses Boltz-2 + tFold cross-validation
  - Similarity computed with TCR-facing surface descriptors (ESP Hodgkin,
    Exposed Shape, rSASA Profile, Exposed Hydrophobicity) — v2 descriptors
    that measure the solvent-exposed surface TCR actually contacts
  - Scoring: 50% RMSD geometry + 50% TCR-facing descriptors
  - Report format matches GILGFVFTL_summary
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Setup ──────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("KVAELVHFL_pipeline")

# ── Parameters ─────────────────────────────────────────────────────────
TARGET = "KVAELVHFL"
HLA = "HLA-A*02:01"

# Output directories
SUMMARY_DIR = PROJECT_ROOT / "data" / f"{TARGET}_summary"
DECOY_A_DIR = SUMMARY_DIR / "Decoy_A"
DECOY_B_DIR = SUMMARY_DIR / "Decoy_B"
DECOY_D_DIR = SUMMARY_DIR / "Decoy_D"

for d in [DECOY_A_DIR, DECOY_B_DIR / "3D_structures", DECOY_D_DIR / "3D_structures", DECOY_D_DIR / "target_pdb"]:
    d.mkdir(parents=True, exist_ok=True)

SEP = "=" * 80


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STEP 1: Decoy A — Sequence Homology Scan                          ║
# ╚══════════════════════════════════════════════════════════════════════╝
def run_decoy_a():
    """Run Decoy A: Hamming-distance-based sequence homology scan."""
    print(f"\n{SEP}")
    print(f"  STEP 1: Decoy A — Sequence Homology Scan for {TARGET}")
    print(SEP)

    from decoy_a.scanner import scan_decoy_a

    t0 = time.time()
    hits = scan_decoy_a(
        target_sequence=TARGET,
        hla_allele=HLA,
        max_hamming=4,  # match GILGFVFTL config
    )
    elapsed = time.time() - t0

    # Convert to JSON-serializable format
    results = []
    for h in hits:
        entry = h.model_dump(mode="json") if hasattr(h, "model_dump") else h.__dict__
        results.append(entry)

    # Save
    output_path = DECOY_A_DIR / "decoy_a_results.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"  Decoy A hits: {len(results)}")
    print(f"  Elapsed: {elapsed:.1f}s")
    print(f"  Output: {output_path}")

    # Print top 10
    if results:
        print(f"\n  Top 10 Decoy A hits:")
        print(f"  {'Sequence':<15} {'HD':>3} {'EL_Rank':>8} {'Genes'}")
        for r in results[:10]:
            genes = ", ".join(r.get("gene_symbols", [])[:3])
            print(f"  {r['sequence']:<15} {r['hamming_distance']:>3} {r.get('el_rank', 0):>8.3f} {genes}")

    print(SEP)
    return results


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STEP 2: Decoy B — Structural Similarity with Boltz+tFold Cross-Val║
# ╚══════════════════════════════════════════════════════════════════════╝
def run_decoy_b():
    """
    Run Decoy B: Structural/physicochemical similarity scan.
    Uses tFold for bulk prediction + Boltz for cross-validation.
    Computes similarity with new interface descriptors.
    """
    print(f"\n{SEP}")
    print(f"  STEP 2: Decoy B — Structural Similarity (Boltz+tFold Cross-Val)")
    print(SEP)

    from decoy_b.scanner import scan_decoy_b

    t0 = time.time()
    hits = scan_decoy_b(
        target_sequence=TARGET,
        hla_allele=HLA,
        run_structural=True,
        run_af3_refinement=False,  # Skip AF3 (slow), rely on tFold+Boltz
        run_boltz_crossval=True,   # Enable Boltz cross-validation
        run_mpnn=False,            # MPNN handled separately in Decoy D
        cosine_threshold=0.70,
        top_k=50,
        boltz_top_n=50,            # Cross-validate top 50
    )
    elapsed = time.time() - t0

    print(f"  Decoy B hits: {len(hits)}")
    print(f"  Elapsed: {elapsed:.1f}s")

    # Build final ranked entries using risk scorer
    from decoy_b.risk_scorer import score_and_rank

    # We need Decoy A hits too for the combined scoring
    decoy_a_path = DECOY_A_DIR / "decoy_a_results.json"
    decoy_a_hits = []
    if decoy_a_path.exists():
        with open(decoy_a_path) as f:
            decoy_a_data = json.load(f)
        from decoy_a.models import DecoyAHit
        for entry in decoy_a_data:
            try:
                decoy_a_hits.append(DecoyAHit.model_validate(entry))
            except Exception:
                pass

    entries = score_and_rank(
        decoy_a_hits=decoy_a_hits,
        decoy_b_hits=hits,
        top_n=50,
    )

    # Convert to JSON-serializable format
    results = []
    for e in entries:
        entry = e.model_dump(mode="json") if hasattr(e, "model_dump") else e.__dict__
        results.append(entry)

    # Save final ranked decoys
    output_path = DECOY_B_DIR / "final_ranked_decoys.json"
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"  Final ranked entries: {len(results)}")
    print(f"  Output: {output_path}")

    # Copy 3D structures to summary directory
    n_copied = 0
    struct_dir = DECOY_B_DIR / "3D_structures"
    for entry in results:
        structural = entry.get("structural", {})
        if structural:
            pdb_path = structural.get("pdb_path")
            if pdb_path and Path(pdb_path).exists():
                dest = struct_dir / Path(pdb_path).name
                if not dest.exists():
                    shutil.copy2(pdb_path, dest)
                    n_copied += 1
            boltz_path = structural.get("boltz_pdb_path")
            if boltz_path and Path(boltz_path).exists():
                dest = struct_dir / f"boltz_{Path(boltz_path).name}"
                if not dest.exists():
                    shutil.copy2(boltz_path, dest)
                    n_copied += 1

    print(f"  3D structures copied: {n_copied}")

    # Print top 10
    if results:
        print(f"\n  Top 10 Decoy B entries:")
        print(f"  {'Rank':>4} {'Sequence':<15} {'HD':>3} {'PhysChem':>8} {'Struct':>8} {'Risk':>8} {'CV_Agree':>8}")
        for r in results[:10]:
            struct = r.get("structural", {}) or {}
            cv = struct.get("cross_validation_agreement", "N/A")
            cv_str = f"{cv:.3f}" if isinstance(cv, (int, float)) else "N/A"
            print(
                f"  {r.get('risk_rank', 0):>4} "
                f"{r['sequence']:<15} "
                f"{r.get('hamming_distance', 0):>3} "
                f"{r.get('physicochemical_similarity', 0):>8.3f} "
                f"{r.get('structural_similarity', 0):>8.3f} "
                f"{r.get('total_risk_score', 0):>8.2f} "
                f"{cv_str:>8}"
            )

    print(SEP)
    return results


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STEP 3: Decoy D — MPNN Inverse Design                             ║
# ╚══════════════════════════════════════════════════════════════════════╝
def run_decoy_d():
    """Run Decoy D: ProteinMPNN inverse design + mhcflurry filter."""
    print(f"\n{SEP}")
    print(f"  STEP 3: Decoy D — MPNN Inverse Design for {TARGET}")
    print(SEP)

    from decoy_d.scanner import run_decoy_d as _run_decoy_d

    output_dir = Path("data/decoy_d") / TARGET
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    try:
        df = _run_decoy_d(
            target_sequence=TARGET,
            hla_allele=HLA,
            num_designs=1000,
            el_rank_threshold=2.0,
            top_k_structures=10,
            output_dir=output_dir,
        )
        elapsed = time.time() - t0

        print(f"  Decoy D designs (filtered): {len(df)}")
        print(f"  Elapsed: {elapsed:.1f}s")

        if not df.empty:
            # Save CSV
            csv_path = DECOY_D_DIR / "decoy_d_results.csv"
            df.to_csv(csv_path, index=False)
            print(f"  Output: {csv_path}")

            # Copy 3D structures
            n_copied = 0
            if "pdb_path" in df.columns:
                for _, row in df.iterrows():
                    pdb = row.get("pdb_path")
                    if pd.notna(pdb) and pdb and Path(str(pdb)).exists():
                        dest = DECOY_D_DIR / "3D_structures" / Path(str(pdb)).name
                        if not dest.exists():
                            shutil.copy2(str(pdb), dest)
                            n_copied += 1

            # Copy target PDB
            if "target_pdb" in df.columns:
                target_pdb = df["target_pdb"].iloc[0]
                if pd.notna(target_pdb) and target_pdb and Path(str(target_pdb)).exists():
                    dest = DECOY_D_DIR / "target_pdb" / Path(str(target_pdb)).name
                    if not dest.exists():
                        shutil.copy2(str(target_pdb), dest)

            print(f"  3D structures copied: {n_copied}")

            # Print top 10
            print(f"\n  Top 10 Decoy D designs:")
            print(f"  {'Sequence':<15} {'MPNN_Score':>10} {'EL_Rank':>8}")
            for _, row in df.head(10).iterrows():
                print(f"  {row['sequence']:<15} {row.get('mpnn_score', 0):>10.4f} {row.get('el_rank', 0):>8.4f}")
        else:
            print("  No designs passed filtering.")

    except Exception as e:
        elapsed = time.time() - t0
        log.warning(f"Decoy D failed: {e}")
        print(f"  Decoy D failed: {e}")
        print(f"  Elapsed: {elapsed:.1f}s")
        df = pd.DataFrame()

    print(SEP)
    return df


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STEP 4: Generate Summary Report                                    ║
# ╚══════════════════════════════════════════════════════════════════════╝
def generate_report(decoy_a_results, decoy_b_results, decoy_d_df):
    """Generate a comprehensive summary report."""
    print(f"\n{SEP}")
    print(f"  STEP 4: Generating Summary Report")
    print(SEP)

    report_lines = []
    report_lines.append(f"{'=' * 72}")
    report_lines.append(f"  pMHC Decoy Library — Summary Report")
    report_lines.append(f"  Target Peptide: {TARGET}")
    report_lines.append(f"  HLA Allele:     {HLA}")
    report_lines.append(f"  Date:           {time.strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append(f"{'=' * 72}")
    report_lines.append("")

    # ── Decoy A Summary ──
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Decoy A — Sequence Homology (Hamming ≤ 4)")
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Total hits: {len(decoy_a_results)}")
    if decoy_a_results:
        hd_counts = {}
        for r in decoy_a_results:
            hd = r.get("hamming_distance", 0)
            hd_counts[hd] = hd_counts.get(hd, 0) + 1
        for hd in sorted(hd_counts):
            report_lines.append(f"    HD={hd}: {hd_counts[hd]} hits")

        # TCR contact vs anchor mismatches
        tcr_mis = [r.get("n_tcr_contact_mismatches", 0) for r in decoy_a_results]
        anc_mis = [r.get("n_anchor_mismatches", 0) for r in decoy_a_results]
        report_lines.append(f"  Avg TCR-contact mismatches: {np.mean(tcr_mis):.2f}")
        report_lines.append(f"  Avg anchor mismatches:      {np.mean(anc_mis):.2f}")
    report_lines.append("")

    # ── Decoy B Summary ──
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Decoy B — Structural Similarity (Boltz-2 + tFold Cross-Validation)")
    report_lines.append(f"  Scoring: 50% RMSD geometry + 50% TCR-facing surface descriptors")
    report_lines.append(f"  Descriptors: ESP Hodgkin (w=0.30), Exposed Shape (w=0.30),")
    report_lines.append(f"               rSASA Profile (w=0.20), Exposed Hydrophobicity (w=0.20)")
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Total ranked entries: {len(decoy_b_results)}")
    if decoy_b_results:
        physchem_sims = [r.get("physicochemical_similarity", 0) for r in decoy_b_results]
        struct_sims = [r.get("structural_similarity", 0) for r in decoy_b_results]
        risk_scores = [r.get("total_risk_score", 0) for r in decoy_b_results]

        report_lines.append(f"  Physicochemical similarity: mean={np.mean(physchem_sims):.3f}, max={np.max(physchem_sims):.3f}")
        report_lines.append(f"  Structural similarity:      mean={np.mean(struct_sims):.3f}, max={np.max(struct_sims):.3f}")
        report_lines.append(f"  Risk score:                 mean={np.mean(risk_scores):.2f}, max={np.max(risk_scores):.2f}")

        # Cross-validation stats
        cv_scores = []
        for r in decoy_b_results:
            s = r.get("structural", {}) or {}
            cv = s.get("cross_validation_agreement")
            if cv is not None:
                cv_scores.append(cv)
        if cv_scores:
            report_lines.append(f"  Cross-validation agreement: mean={np.mean(cv_scores):.3f}, n={len(cv_scores)}")

        # TCR-facing surface descriptor stats (v2)
        iface_keys = [
            ("esp_similarity", "TCR-facing ESP (Hodgkin SI)"),
            ("shape_similarity", "Exposed Shape (SC centroid RMSD)"),
            ("sasa_similarity", "rSASA Profile"),
            ("hydrophobicity_similarity", "Exposed Hydrophobicity"),
            ("interface_combined", "TCR-surface Combined"),
        ]
        for key, label in iface_keys:
            vals = []
            for r in decoy_b_results:
                s = r.get("structural", {}) or {}
                v = s.get(key)
                if v is not None:
                    vals.append(v)
            if vals:
                report_lines.append(f"  {label}: mean={np.mean(vals):.3f}, max={np.max(vals):.3f}, n={len(vals)}")

        # Top 5 highest risk
        report_lines.append(f"\n  Top 5 highest-risk Decoy B entries:")
        report_lines.append(f"  {'Rank':>4} {'Sequence':<15} {'Risk':>8} {'Source'}")
        for r in decoy_b_results[:5]:
            src = r.get("source", "")
            report_lines.append(
                f"  {r.get('risk_rank', 0):>4} {r['sequence']:<15} "
                f"{r.get('total_risk_score', 0):>8.2f} {src}"
            )
    report_lines.append("")

    # ── Decoy D Summary ──
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Decoy D — MPNN Inverse Design")
    report_lines.append(f"{'─' * 72}")
    n_d = len(decoy_d_df) if decoy_d_df is not None and not decoy_d_df.empty else 0
    report_lines.append(f"  Total designs (HLA-filtered): {n_d}")
    if n_d > 0:
        report_lines.append(f"  MPNN score range: {decoy_d_df['mpnn_score'].min():.4f} - {decoy_d_df['mpnn_score'].max():.4f}")
        report_lines.append(f"  EL%Rank range:    {decoy_d_df['el_rank'].min():.4f} - {decoy_d_df['el_rank'].max():.4f}")
    report_lines.append("")

    # ── Output Files ──
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  Output Files")
    report_lines.append(f"{'─' * 72}")
    report_lines.append(f"  {SUMMARY_DIR}/")
    report_lines.append(f"    Decoy_A/decoy_a_results.json")
    report_lines.append(f"    Decoy_B/final_ranked_decoys.json")
    report_lines.append(f"    Decoy_B/3D_structures/")
    report_lines.append(f"    Decoy_D/decoy_d_results.csv")
    report_lines.append(f"    Decoy_D/3D_structures/")
    report_lines.append(f"    Decoy_D/target_pdb/")
    report_lines.append(f"{'=' * 72}")

    report_text = "\n".join(report_lines)
    print(report_text)

    # Save report
    report_path = SUMMARY_DIR / "pipeline_report.txt"
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"\n  Report saved to: {report_path}")

    print(SEP)
    return report_text


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                               ║
# ╚══════════════════════════════════════════════════════════════════════╝
def main():
    wall_start = time.time()

    print(f"\n{'#' * 80}")
    print(f"  pMHC Decoy Library Pipeline — {TARGET} / {HLA}")
    print(f"  Output: {SUMMARY_DIR}")
    print(f"{'#' * 80}")

    # Step 1: Decoy A
    decoy_a_results = run_decoy_a()

    # Step 2: Decoy B (with Boltz+tFold cross-validation)
    decoy_b_results = run_decoy_b()

    # Step 3: Decoy D (MPNN inverse design)
    decoy_d_df = run_decoy_d()

    # Step 4: Summary report
    generate_report(decoy_a_results, decoy_b_results, decoy_d_df)

    wall_elapsed = time.time() - wall_start
    print(f"\n  Total pipeline time: {wall_elapsed:.1f}s ({wall_elapsed/60:.1f} min)")
    print(f"  Summary directory: {SUMMARY_DIR}")


if __name__ == "__main__":
    main()
