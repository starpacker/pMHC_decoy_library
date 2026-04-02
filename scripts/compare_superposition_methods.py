"""
Compare two pMHC structural superposition strategies
=====================================================

Method A (current): Superimpose on MHC Cα → measure peptide RMSD
Method B (new):     Superimpose on peptide Cα → measure MHC groove-contact RMSD

Goal: determine if both methods give consistent rankings for decoy candidates.

Usage:
    python scripts/compare_superposition_methods.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── BioPython imports ────────────────────────────────────────────────────
from Bio.PDB import PDBParser, Superimposer

# ── Paths ────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data" / "GILGFVFTL_summary"
PDB_DIR = DATA / "Decoy_B" / "3D_structures" / "tfold"
TARGET_PDB = PDB_DIR / "pmhc_GILGFVFTL_HLA-A0201.pdb"
FIG_DIR = REPO / "figures"
FIG_DIR.mkdir(exist_ok=True)

# tFold chain convention
PEPTIDE_CHAIN = "P"
MHC_CHAINS = {"M", "N"}
# MHC α1/α2 helices flanking the peptide groove (approximate residue range on chain M)
# These are the residues that directly contact the peptide in class I MHC.
# α1 helix: ~50-85, α2 helix: ~138-175 (HLA-A*02:01 numbering)
GROOVE_HELIX_RANGES = [(50, 85), (138, 175)]


def _get_cas(struct, chain_ids):
    """Extract CA atoms from specified chains."""
    cas = []
    for model in struct:
        for chain in model:
            if chain.id in chain_ids:
                for residue in chain:
                    if residue.id[0] != " ":
                        continue  # skip hetero/water
                    for atom in residue:
                        if atom.get_name() == "CA":
                            cas.append(atom)
    return cas


def _get_groove_cas(struct):
    """Extract CA atoms from MHC groove-contact helices (chain M only).
    These are the α1/α2 helices that line the peptide-binding groove."""
    cas = []
    for model in struct:
        for chain in model:
            if chain.id != "M":
                continue
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resnum = residue.id[1]
                in_groove = any(lo <= resnum <= hi for lo, hi in GROOVE_HELIX_RANGES)
                if in_groove:
                    for atom in residue:
                        if atom.get_name() == "CA":
                            cas.append(atom)
    return cas


def _coords(atoms):
    return np.array([a.get_vector().get_array() for a in atoms])


def _rmsd(coords1, coords2):
    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))


def _core_cas(cas_list, n_core=5):
    """Central TCR-contact core residues."""
    n = len(cas_list)
    if n <= n_core:
        return cas_list
    start = (n - n_core) // 2
    return cas_list[start : start + n_core]


def method_a(target_struct, cand_struct):
    """
    Method A (current): Superimpose on MHC → measure peptide RMSD.
    Returns dict with peptide_rmsd, mhc_fit_rmsd.
    """
    t_mhc = _get_cas(target_struct, MHC_CHAINS)
    c_mhc = _get_cas(cand_struct, MHC_CHAINS)
    t_pep = _get_cas(target_struct, {PEPTIDE_CHAIN})
    c_pep = _get_cas(cand_struct, {PEPTIDE_CHAIN})

    if len(t_mhc) != len(c_mhc) or len(t_mhc) == 0:
        return None

    sup = Superimposer()
    sup.set_atoms(t_mhc, c_mhc)
    mhc_fit_rmsd = sup.rms

    # Apply MHC rotation to ALL candidate atoms
    sup.apply(list(cand_struct.get_atoms()))

    # Peptide RMSD (same-length or core)
    if len(t_pep) == len(c_pep) and len(t_pep) > 0:
        pep_rmsd = _rmsd(_coords(t_pep), _coords(c_pep))
    elif len(t_pep) > 0 and len(c_pep) > 0:
        tc = _core_cas(t_pep)
        cc = _core_cas(c_pep)
        if len(tc) == len(cc):
            pep_rmsd = _rmsd(_coords(tc), _coords(cc))
        else:
            pep_rmsd = None
    else:
        pep_rmsd = None

    # Also measure groove RMSD after MHC superposition (should be ~0 since we fitted on MHC)
    t_groove = _get_groove_cas(target_struct)
    c_groove = _get_groove_cas(cand_struct)
    groove_rmsd = None
    if len(t_groove) == len(c_groove) and len(t_groove) > 0:
        groove_rmsd = _rmsd(_coords(t_groove), _coords(c_groove))

    return {
        "mhc_fit_rmsd": mhc_fit_rmsd,
        "peptide_rmsd": pep_rmsd,
        "groove_rmsd_after_mhc_fit": groove_rmsd,
    }


def method_b(target_struct, cand_struct):
    """
    Method B (new): Superimpose on peptide Cα → measure MHC groove RMSD.

    Rationale: if two peptides adopt the same backbone conformation, and the
    MHC groove also stays in the same place relative to the peptide, then the
    overall pMHC surface presented to TCR is conserved. Measuring the groove
    deviation after peptide alignment captures how much the MHC 'context'
    differs even when the peptide itself is well-aligned.
    """
    t_pep = _get_cas(target_struct, {PEPTIDE_CHAIN})
    c_pep = _get_cas(cand_struct, {PEPTIDE_CHAIN})

    if len(t_pep) != len(c_pep) or len(t_pep) == 0:
        # Cross-length: use core
        t_pep = _core_cas(t_pep)
        c_pep = _core_cas(c_pep)
        if len(t_pep) != len(c_pep) or len(t_pep) == 0:
            return None

    sup = Superimposer()
    sup.set_atoms(t_pep, c_pep)
    pep_fit_rmsd = sup.rms

    # Apply peptide rotation to ALL candidate atoms
    sup.apply(list(cand_struct.get_atoms()))

    # Measure MHC groove deviation
    t_groove = _get_groove_cas(target_struct)
    c_groove = _get_groove_cas(cand_struct)

    groove_rmsd = None
    if len(t_groove) == len(c_groove) and len(t_groove) > 0:
        groove_rmsd = _rmsd(_coords(t_groove), _coords(c_groove))

    # Also measure full MHC deviation
    t_mhc = _get_cas(target_struct, MHC_CHAINS)
    c_mhc = _get_cas(cand_struct, MHC_CHAINS)
    full_mhc_rmsd = None
    if len(t_mhc) == len(c_mhc) and len(t_mhc) > 0:
        full_mhc_rmsd = _rmsd(_coords(t_mhc), _coords(c_mhc))

    return {
        "pep_fit_rmsd": pep_fit_rmsd,
        "groove_rmsd_after_pep_fit": groove_rmsd,
        "full_mhc_rmsd_after_pep_fit": full_mhc_rmsd,
    }


def main():
    parser = PDBParser(QUIET=True)

    if not TARGET_PDB.exists():
        print(f"Target PDB not found: {TARGET_PDB}")
        sys.exit(1)

    # Find all candidate PDBs
    cand_pdbs = sorted(PDB_DIR.glob("pmhc_*_HLA-A0201.pdb"))
    cand_pdbs = [p for p in cand_pdbs if p.name != TARGET_PDB.name]

    print(f"Target: {TARGET_PDB.name}")
    print(f"Candidates: {len(cand_pdbs)}")
    print("=" * 100)

    rows = []
    for cand_pdb in cand_pdbs:
        seq = cand_pdb.stem.split("_")[1]  # pmhc_XXXXX_HLA-A0201 → XXXXX

        # Must re-parse for each candidate because Superimposer modifies coords in-place
        target_struct = parser.get_structure("target", str(TARGET_PDB))
        cand_struct = parser.get_structure("cand", str(cand_pdb))

        res_a = method_a(target_struct, cand_struct)

        # Re-parse again for method B (method A modified cand_struct in-place)
        target_struct = parser.get_structure("target", str(TARGET_PDB))
        cand_struct = parser.get_structure("cand", str(cand_pdb))

        res_b = method_b(target_struct, cand_struct)

        if res_a is None or res_b is None:
            print(f"  SKIP {seq}: superposition failed")
            continue

        row = {"sequence": seq}
        row.update({f"A_{k}": v for k, v in res_a.items()})
        row.update({f"B_{k}": v for k, v in res_b.items()})

        # Convert RMSD → similarity (0-1, cutoff 3Å)
        if res_a["peptide_rmsd"] is not None:
            row["A_pep_similarity"] = max(0.0, 1.0 - res_a["peptide_rmsd"] / 3.0)
        if res_b["groove_rmsd_after_pep_fit"] is not None:
            row["B_groove_similarity"] = max(0.0, 1.0 - res_b["groove_rmsd_after_pep_fit"] / 3.0)

        rows.append(row)

    df = pd.DataFrame(rows)

    # ── Print results table ──────────────────────────────────────────────
    print(f"\n{'Sequence':<14} | {'A: pep RMSD':>11} {'A: pep sim':>10} | "
          f"{'B: pep fit':>10} {'B: groove':>10} {'B: grv sim':>10} | {'Rank A':>6} {'Rank B':>6}")
    print("-" * 100)

    df_sorted = df.sort_values("A_pep_similarity", ascending=False).copy()
    df_sorted["rank_A"] = range(1, len(df_sorted) + 1)
    df_sorted = df_sorted.sort_values("B_groove_similarity", ascending=False)
    rank_b_map = {seq: i + 1 for i, seq in enumerate(df_sorted["sequence"])}
    df_sorted = df_sorted.sort_values("A_pep_similarity", ascending=False)
    df_sorted["rank_B"] = df_sorted["sequence"].map(rank_b_map)

    for _, r in df_sorted.iterrows():
        print(f"{r['sequence']:<14} | "
              f"{r.get('A_peptide_rmsd', 0):>11.4f} {r.get('A_pep_similarity', 0):>10.4f} | "
              f"{r.get('B_pep_fit_rmsd', 0):>10.4f} {r.get('B_groove_rmsd_after_pep_fit', 0):>10.4f} "
              f"{r.get('B_groove_similarity', 0):>10.4f} | "
              f"{int(r['rank_A']):>6} {int(r['rank_B']):>6}")

    # ── Correlation analysis ─────────────────────────────────────────────
    valid = df.dropna(subset=["A_pep_similarity", "B_groove_similarity"])
    if len(valid) > 2:
        from scipy import stats
        pearson_r, pearson_p = stats.pearsonr(valid["A_pep_similarity"], valid["B_groove_similarity"])
        spearman_r, spearman_p = stats.spearmanr(valid["A_pep_similarity"], valid["B_groove_similarity"])
        kendall_tau, kendall_p = stats.kendalltau(valid["A_pep_similarity"], valid["B_groove_similarity"])
        print(f"\n{'='*60}")
        print(f"Correlation between Method A and Method B (n={len(valid)}):")
        print(f"  Pearson  r = {pearson_r:.4f}  (p = {pearson_p:.2e})")
        print(f"  Spearman ρ = {spearman_r:.4f}  (p = {spearman_p:.2e})")
        print(f"  Kendall  τ = {kendall_tau:.4f}  (p = {kendall_p:.2e})")
        print(f"{'='*60}")

    # ── Visualization ────────────────────────────────────────────────────
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        fig.suptitle(
            "Superposition Method Comparison: GILGFVFTL Decoys",
            fontsize=14, fontweight="bold",
        )

        # 1. Scatter: Method A peptide similarity vs Method B groove similarity
        ax = axes[0, 0]
        ax.scatter(valid["A_pep_similarity"], valid["B_groove_similarity"],
                   c="#2196F3", s=50, edgecolors="white", linewidths=0.5, alpha=0.8)
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="y=x")
        ax.set_xlabel("Method A: Peptide Similarity\n(MHC-superposed → peptide RMSD)")
        ax.set_ylabel("Method B: Groove Similarity\n(peptide-superposed → groove RMSD)")
        ax.set_title(f"Pearson r = {pearson_r:.3f}, Spearman ρ = {spearman_r:.3f}")
        ax.legend()
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        # Label outliers (large rank difference)
        df_sorted_copy = df_sorted.copy()
        df_sorted_copy["rank_diff"] = abs(df_sorted_copy["rank_A"] - df_sorted_copy["rank_B"])
        outliers = df_sorted_copy.nlargest(5, "rank_diff")
        for _, r in outliers.iterrows():
            ax.annotate(r["sequence"],
                        (r.get("A_pep_similarity", 0), r.get("B_groove_similarity", 0)),
                        fontsize=6, ha="left", xytext=(3, 3), textcoords="offset points")

        # 2. Scatter: Method A peptide RMSD vs Method B groove RMSD
        ax = axes[0, 1]
        ax.scatter(valid["A_peptide_rmsd"], valid["B_groove_rmsd_after_pep_fit"],
                   c="#FF9800", s=50, edgecolors="white", linewidths=0.5, alpha=0.8)
        ax.set_xlabel("Method A: Peptide RMSD (Å)\n(after MHC superposition)")
        ax.set_ylabel("Method B: Groove RMSD (Å)\n(after peptide superposition)")
        ax.set_title("Raw RMSD Comparison")
        for _, r in outliers.iterrows():
            ax.annotate(r["sequence"],
                        (r.get("A_peptide_rmsd", 0), r.get("B_groove_rmsd_after_pep_fit", 0)),
                        fontsize=6, ha="left", xytext=(3, 3), textcoords="offset points")

        # 3. Rank comparison
        ax = axes[1, 0]
        ax.scatter(df_sorted["rank_A"], df_sorted["rank_B"],
                   c="#4CAF50", s=50, edgecolors="white", linewidths=0.5, alpha=0.8)
        max_rank = max(df_sorted["rank_A"].max(), df_sorted["rank_B"].max())
        ax.plot([1, max_rank], [1, max_rank], "k--", alpha=0.3, label="perfect agreement")
        ax.set_xlabel("Rank by Method A (peptide RMSD)")
        ax.set_ylabel("Rank by Method B (groove RMSD)")
        ax.set_title("Rank Agreement")
        ax.legend()

        # 4. Distribution comparison
        ax = axes[1, 1]
        ax.hist(valid["A_pep_similarity"], bins=15, alpha=0.6, color="#2196F3",
                label="A: peptide sim", edgecolor="white")
        ax.hist(valid["B_groove_similarity"], bins=15, alpha=0.6, color="#FF9800",
                label="B: groove sim", edgecolor="white")
        ax.set_xlabel("Similarity Score (0-1)")
        ax.set_ylabel("Count")
        ax.set_title("Score Distributions")
        ax.legend()

        plt.tight_layout()
        out_path = FIG_DIR / "superposition_method_comparison.png"
        plt.savefig(out_path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"\nFigure saved: {out_path}")

    except ImportError:
        print("\nmatplotlib not available — skipping visualization")

    # ── Save CSV ─────────────────────────────────────────────────────────
    out_csv = FIG_DIR / "superposition_method_comparison.csv"
    df_sorted.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"CSV saved: {out_csv}")


if __name__ == "__main__":
    main()
