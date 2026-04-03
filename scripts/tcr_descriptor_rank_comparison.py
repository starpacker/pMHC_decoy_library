"""
TCR-Facing Descriptor Solo Rank Comparison
===========================================
Compare the new TCR-facing descriptors against RMSD metrics
to verify they provide consistent, meaningful rankings.
"""

import json
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from decoy_b.tools.tcr_surface_descriptors import (
    compute_tcr_facing_descriptors,
    compute_esp_hodgkin_similarity,
    compute_sasa_similarity,
    compute_hydrophobicity_similarity,
    compute_shape_similarity,
)
from Bio.PDB import PDBParser, Superimposer

STRUCT_DIR = PROJECT_ROOT / "data" / "GILGFVFTL_summary" / "Decoy_B" / "3D_structures"
TARGET_PDB = STRUCT_DIR / "pmhc_GILGFVFTL_HLA-A0201.pdb"


def _get_cas_by_chains(struct, chain_ids):
    cas = []
    for model in struct:
        for chain in model:
            if chain.id in chain_ids:
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    for atom in res:
                        if atom.get_name() == "CA":
                            cas.append(atom)
        break
    return cas


def calc_rmsd(a, b):
    diff = a - b
    return float(np.sqrt((diff ** 2).sum() / len(a)))


def main():
    parser = PDBParser(QUIET=True)
    target_seq = "GILGFVFTL"

    all_pdbs = sorted(STRUCT_DIR.glob("pmhc_*_HLA-A0201.pdb"))
    candidate_pdbs = [p for p in all_pdbs if target_seq not in p.name]
    print(f"Candidates: {len(candidate_pdbs)}")

    # Compute target descriptors
    print("Computing target TCR-facing descriptors...")
    desc_target = compute_tcr_facing_descriptors(str(TARGET_PDB))

    results = []
    for i, pdb_path in enumerate(candidate_pdbs):
        seq = pdb_path.stem.split("_")[1]
        print(f"  [{i+1}/{len(candidate_pdbs)}] {seq}...", end="", flush=True)

        row = {"sequence": seq}

        # TCR-facing descriptors
        desc_cand = compute_tcr_facing_descriptors(str(pdb_path))

        if desc_target.esp_profile is not None and desc_cand.esp_profile is not None:
            row["tcr_esp"] = compute_esp_hodgkin_similarity(
                desc_target.esp_profile, desc_cand.esp_profile)
        else:
            row["tcr_esp"] = None

        if desc_target.sasa_profile is not None and desc_cand.sasa_profile is not None:
            row["tcr_sasa"] = compute_sasa_similarity(
                desc_target.sasa_profile, desc_cand.sasa_profile)
        else:
            row["tcr_sasa"] = None

        if (desc_target.hydrophobicity_profile is not None
                and desc_cand.hydrophobicity_profile is not None):
            row["tcr_hydro"] = compute_hydrophobicity_similarity(
                desc_target.hydrophobicity_profile, desc_cand.hydrophobicity_profile)
        else:
            row["tcr_hydro"] = None

        if (desc_target.exposed_sidechain_coords is not None
                and desc_cand.exposed_sidechain_coords is not None):
            row["tcr_shape"] = compute_shape_similarity(
                desc_target.exposed_sidechain_coords,
                desc_cand.exposed_sidechain_coords)
        else:
            row["tcr_shape"] = None

        # RMSD Method A (for comparison)
        try:
            tgt = parser.get_structure(f"t{i}", str(TARGET_PDB))
            cand = parser.get_structure(f"c{i}", str(pdb_path))
            t_mhc = _get_cas_by_chains(tgt, {"M", "N"})
            c_mhc = _get_cas_by_chains(cand, {"M", "N"})
            t_pep = _get_cas_by_chains(tgt, {"P"})
            c_pep = _get_cas_by_chains(cand, {"P"})

            if len(t_mhc) == len(c_mhc) and len(t_mhc) > 0:
                sup = Superimposer()
                sup.set_atoms(t_mhc, c_mhc)
                sup.apply(list(cand.get_atoms()))
                if len(t_pep) == len(c_pep) and len(t_pep) > 0:
                    t_coords = np.array([a.get_vector().get_array() for a in t_pep])
                    c_coords = np.array([a.get_vector().get_array() for a in c_pep])
                    rmsd = calc_rmsd(t_coords, c_coords)
                    row["sim_A"] = max(0.0, 1.0 - rmsd / 3.0)
                else:
                    row["sim_A"] = 0.0
            else:
                row["sim_A"] = 0.0
        except Exception:
            row["sim_A"] = 0.0

        # Combined TCR-facing score
        scores = {}
        w = {"esp": 0.30, "shape": 0.30, "sasa": 0.20, "hydrophobicity": 0.20}
        if row.get("tcr_esp") is not None:
            scores["esp"] = row["tcr_esp"]
        if row.get("tcr_sasa") is not None:
            scores["sasa"] = row["tcr_sasa"]
        if row.get("tcr_hydro") is not None:
            scores["hydrophobicity"] = row["tcr_hydro"]
        if row.get("tcr_shape") is not None:
            scores["shape"] = row["tcr_shape"]

        if scores:
            aw = {k: w[k] for k in scores}
            ws = sum(aw.values())
            row["tcr_combined"] = sum((aw[k] / ws) * scores[k] for k in scores)
        else:
            row["tcr_combined"] = 0.0

        results.append(row)
        print(f" esp={row.get('tcr_esp','N/A'):.3f}"
              f" sasa={row.get('tcr_sasa','N/A'):.3f}"
              f" shape={row.get('tcr_shape','N/A'):.3f}"
              f" combined={row['tcr_combined']:.3f}")

    # Ranking & correlation
    print("\n" + "=" * 90)
    print("SPEARMAN RANK CORRELATION MATRIX")
    print("=" * 90)

    from scipy.stats import spearmanr

    metrics = ["tcr_esp", "tcr_sasa", "tcr_hydro", "tcr_shape", "tcr_combined", "sim_A"]
    ranks = {}
    for m in metrics:
        vals = [r.get(m) if r.get(m) is not None else -999 for r in results]
        order = np.argsort(vals)[::-1]
        rank_arr = np.zeros(len(vals), dtype=int)
        for rank, idx in enumerate(order, 1):
            rank_arr[idx] = rank
        ranks[m] = rank_arr

    header = f"{'':>16}"
    for m in metrics:
        header += f"{m[:14]:>15}"
    print(header)

    for m1 in metrics:
        line = f"{m1[:16]:>16}"
        for m2 in metrics:
            rho, _ = spearmanr(ranks[m1], ranks[m2])
            line += f"{rho:>15.3f}"
        print(line)

    # Top 15 by tcr_combined
    print(f"\n{'=' * 90}")
    print("TOP 15 BY TCR-FACING COMBINED SCORE")
    print(f"{'=' * 90}")
    sorted_idx = np.argsort(ranks["tcr_combined"])
    print(f"{'Rank':>4}  {'Sequence':<12} {'TCR_comb':>9} {'ESP':>6} {'SASA':>6} "
          f"{'Hydro':>6} {'Shape':>6} {'sim_A':>6}")
    for pos, idx in enumerate(sorted_idx[:15]):
        r = results[idx]
        print(f"{pos+1:>4}  {r['sequence']:<12} {r['tcr_combined']:>9.4f} "
              f"{r.get('tcr_esp',0):>6.3f} {r.get('tcr_sasa',0):>6.3f} "
              f"{r.get('tcr_hydro',0):>6.3f} {r.get('tcr_shape',0):>6.3f} "
              f"{r['sim_A']:>6.3f}")

    # Save
    out_path = STRUCT_DIR.parent / "tcr_descriptor_solo_ranks.json"
    output = []
    for i, r in enumerate(results):
        entry = dict(r)
        for m in metrics:
            entry[f"rank_{m}"] = int(ranks[m][i])
        output.append(entry)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved to: {out_path}")


if __name__ == "__main__":
    main()
