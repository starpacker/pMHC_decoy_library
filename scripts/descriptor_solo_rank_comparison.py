"""
Descriptor Solo Rank Comparison
================================
For each pMHC decoy candidate, compute individual descriptor scores
against the GILGFVFTL target, rank by each descriptor independently,
and compare rank consistency (Spearman correlation matrix).

Available descriptors on this machine (BioPython + PyTorch only):
  1. Atchley cosine similarity (sequence-level)
  2. RMSD Method A: MHC-superposed peptide RMSD → similarity
  3. RMSD Method B: peptide-superposed groove RMSD → similarity
  4. B-factor correlation
  5. ESP Coulomb proxy (BioPython)
  6. PeSTo contact proxy (BioPython)
  7. Current combined score (surface_correlation from JSON)
"""

import json
import sys
from pathlib import Path

import numpy as np

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from Bio.PDB import PDBParser, Superimposer
from decoy_a.config import ATCHLEY_FACTORS
from decoy_b.tools.interface_descriptors import (
    compute_esp_vector,
    compute_esp_similarity,
    compute_pesto_embedding,
    compute_pesto_similarity,
)

# ── Constants ──────────────────────────────────────────────────────────
STRUCT_DIR = PROJECT_ROOT / "data" / "GILGFVFTL_summary" / "Decoy_B" / "3D_structures"
JSON_PATH = STRUCT_DIR.parent / "final_ranked_decoys.json"
TARGET_PDB = STRUCT_DIR / "pmhc_GILGFVFTL_HLA-A0201.pdb"
PEPTIDE_CHAIN = "P"
MHC_CHAINS = {"M", "N"}


# ── Helpers from scanner.py ────────────────────────────────────────────
def _get_tcr_contact_residues(sequence: str) -> str:
    n = len(sequence)
    if n <= 5:
        return sequence
    n_contact = min(5, n - 2)
    start = (n - n_contact) // 2
    return sequence[start: start + n_contact]


def _sequence_to_atchley_vector(sequence: str) -> np.ndarray:
    factors = []
    for aa in sequence:
        if aa in ATCHLEY_FACTORS:
            factors.extend(ATCHLEY_FACTORS[aa])
        else:
            factors.extend([0.0] * 5)
    return np.array(factors, dtype=np.float64)


def cosine_sim(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.dot(a, b) / (na * nb))


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


def _ca_coords(atoms):
    return np.array([a.get_vector().get_array() for a in atoms])


def _core_cas(cas):
    """Central residues for cross-length comparison."""
    n = len(cas)
    if n <= 5:
        return cas
    start = (n - 5) // 2
    return cas[start: start + 5]


GROOVE_RANGES = [(50, 85), (138, 175)]

def _get_groove_cas(struct):
    cas = []
    for model in struct:
        for chain in model:
            if chain.id != "M":
                continue
            for res in chain:
                if res.id[0] != " ":
                    continue
                resnum = res.id[1]
                in_groove = any(lo <= resnum <= hi for lo, hi in GROOVE_RANGES)
                if in_groove:
                    for atom in res:
                        if atom.get_name() == "CA":
                            cas.append(atom)
        break
    return cas


def calc_rmsd(coords_a, coords_b):
    diff = coords_a - coords_b
    return float(np.sqrt((diff ** 2).sum() / len(coords_a)))


# ── Main ───────────────────────────────────────────────────────────────
def main():
    parser = PDBParser(QUIET=True)

    # Load existing JSON for the 7 Decoy_B entries
    with open(JSON_PATH) as f:
        all_data = json.load(f)

    decoy_b_entries = [d for d in all_data if d["source"] == "Decoy_B_Structural_Similarity"]
    print(f"Decoy B entries from JSON: {len(decoy_b_entries)}")

    # Also find ALL candidate PDBs (not just the 7 in JSON)
    target_seq = "GILGFVFTL"
    all_pdbs = sorted(STRUCT_DIR.glob("pmhc_*_HLA-A0201.pdb"))
    candidate_pdbs = [p for p in all_pdbs if target_seq not in p.name]
    print(f"Total candidate PDB files: {len(candidate_pdbs)}")

    # Extract peptide sequence from PDB filename
    def seq_from_pdb(p):
        name = p.stem  # pmhc_XXXXX_HLA-A0201
        return name.split("_")[1]

    # ── Compute target descriptors ─────────────────────────────────────
    target_contact = _get_tcr_contact_residues(target_seq)
    target_atchley = _sequence_to_atchley_vector(target_contact)

    target_struct_a = parser.get_structure("tgt_a", str(TARGET_PDB))
    target_mhc_cas = _get_cas_by_chains(target_struct_a, MHC_CHAINS)
    target_pep_cas = _get_cas_by_chains(target_struct_a, {PEPTIDE_CHAIN})
    target_pep_bf = np.array([a.get_bfactor() for a in target_pep_cas])

    print("Computing ESP vector for target...")
    target_esp = compute_esp_vector(str(TARGET_PDB))
    print("Computing PeSTo embedding for target...")
    target_pesto = compute_pesto_embedding(str(TARGET_PDB))

    print(f"Target ESP: {'OK' if target_esp is not None else 'FAILED'} "
          f"(dim={len(target_esp) if target_esp is not None else 0})")
    print(f"Target PeSTo: {'OK' if target_pesto is not None else 'FAILED'} "
          f"(dim={len(target_pesto) if target_pesto is not None else 0})")

    # ── Compute per-candidate scores ───────────────────────────────────
    results = []

    for i, pdb_path in enumerate(candidate_pdbs):
        seq = seq_from_pdb(pdb_path)
        print(f"  [{i+1}/{len(candidate_pdbs)}] {seq}...", end="", flush=True)

        row = {"sequence": seq, "pdb": pdb_path.name}

        # 1. Atchley cosine
        contact = _get_tcr_contact_residues(seq)
        vec = _sequence_to_atchley_vector(contact)
        row["atchley_cosine"] = cosine_sim(target_atchley, vec)

        # 2-3. Dual superposition RMSD
        try:
            cand_struct = parser.get_structure(f"c{i}", str(pdb_path))
            cand_mhc_cas = _get_cas_by_chains(cand_struct, MHC_CHAINS)
            cand_pep_cas = _get_cas_by_chains(cand_struct, {PEPTIDE_CHAIN})

            # Method A: MHC superpose → peptide RMSD
            pep_rmsd = None
            if len(target_mhc_cas) == len(cand_mhc_cas) and len(target_mhc_cas) > 0:
                # Re-parse target each time (superimposer modifies coords)
                tgt_a = parser.get_structure(f"ta{i}", str(TARGET_PDB))
                c_a = parser.get_structure(f"ca{i}", str(pdb_path))
                t_mhc = _get_cas_by_chains(tgt_a, MHC_CHAINS)
                c_mhc = _get_cas_by_chains(c_a, MHC_CHAINS)
                t_pep = _get_cas_by_chains(tgt_a, {PEPTIDE_CHAIN})
                c_pep = _get_cas_by_chains(c_a, {PEPTIDE_CHAIN})

                sup = Superimposer()
                sup.set_atoms(t_mhc, c_mhc)
                sup.apply(list(c_a.get_atoms()))

                if len(t_pep) == len(c_pep) and len(t_pep) > 0:
                    pep_rmsd = calc_rmsd(_ca_coords(t_pep), _ca_coords(c_pep))
                elif len(t_pep) > 0 and len(c_pep) > 0:
                    tc = _core_cas(t_pep)
                    cc = _core_cas(c_pep)
                    if len(tc) == len(cc):
                        pep_rmsd = calc_rmsd(_ca_coords(tc), _ca_coords(cc))

            row["peptide_rmsd"] = pep_rmsd
            row["sim_A"] = max(0.0, 1.0 - pep_rmsd / 3.0) if pep_rmsd is not None else 0.0

            # Method B: peptide superpose → groove RMSD
            groove_rmsd = None
            tgt_b = parser.get_structure(f"tb{i}", str(TARGET_PDB))
            c_b = parser.get_structure(f"cb{i}", str(pdb_path))
            t_pep_b = _get_cas_by_chains(tgt_b, {PEPTIDE_CHAIN})
            c_pep_b = _get_cas_by_chains(c_b, {PEPTIDE_CHAIN})

            t_fit = t_pep_b
            c_fit = c_pep_b
            if len(t_fit) != len(c_fit):
                t_fit = _core_cas(t_fit)
                c_fit = _core_cas(c_fit)

            if len(t_fit) == len(c_fit) and len(t_fit) > 0:
                sup_b = Superimposer()
                sup_b.set_atoms(t_fit, c_fit)
                sup_b.apply(list(c_b.get_atoms()))

                t_groove = _get_groove_cas(tgt_b)
                c_groove = _get_groove_cas(c_b)
                if len(t_groove) == len(c_groove) and len(t_groove) > 0:
                    groove_rmsd = calc_rmsd(_ca_coords(t_groove), _ca_coords(c_groove))

            row["groove_rmsd"] = groove_rmsd
            row["sim_B"] = max(0.0, 1.0 - groove_rmsd / 3.0) if groove_rmsd is not None else 0.0

            # 4. B-factor correlation
            n_bf = min(len(target_pep_cas), len(cand_pep_cas))
            if n_bf > 1:
                t_bf = np.array([a.get_bfactor() for a in target_pep_cas[:n_bf]])
                c_bf = np.array([a.get_bfactor() for a in cand_pep_cas[:n_bf]])
                if np.std(t_bf) > 0 and np.std(c_bf) > 0:
                    row["bf_corr"] = float(np.corrcoef(t_bf, c_bf)[0, 1])
                else:
                    row["bf_corr"] = 0.0
            else:
                row["bf_corr"] = 0.0

        except Exception as exc:
            print(f" STRUCT_ERR({exc})", end="")
            row["peptide_rmsd"] = None
            row["sim_A"] = 0.0
            row["groove_rmsd"] = None
            row["sim_B"] = 0.0
            row["bf_corr"] = 0.0

        # 5. ESP Coulomb proxy
        try:
            cand_esp = compute_esp_vector(str(pdb_path))
            if target_esp is not None and cand_esp is not None:
                row["esp_sim"] = compute_esp_similarity(target_esp, cand_esp)
            else:
                row["esp_sim"] = None
        except Exception:
            row["esp_sim"] = None

        # 6. PeSTo contact proxy
        try:
            cand_pesto = compute_pesto_embedding(str(pdb_path))
            if target_pesto is not None and cand_pesto is not None:
                row["pesto_sim"] = compute_pesto_similarity(target_pesto, cand_pesto)
            else:
                row["pesto_sim"] = None
        except Exception:
            row["pesto_sim"] = None

        # 7. Combined RMSD geo score (same formula as scanner.py)
        sa, sb, bf = row["sim_A"], row["sim_B"], row["bf_corr"]
        if sa > 0 and sb > 0:
            rmsd_geo = 0.5 * sa + 0.5 * sb
        elif sa > 0:
            rmsd_geo = sa
        elif sb > 0:
            rmsd_geo = sb
        else:
            rmsd_geo = 0.0
        row["rmsd_geo"] = max(rmsd_geo, bf)

        results.append(row)
        print(f" OK (atch={row['atchley_cosine']:.3f} simA={row['sim_A']:.3f} simB={row['sim_B']:.3f}"
              f" esp={row.get('esp_sim','N/A')} pesto={row.get('pesto_sim','N/A')})")

    # ── Ranking comparison ─────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("SOLO RANKING COMPARISON")
    print("=" * 100)

    # Metrics to rank by (higher = more similar = higher risk)
    metrics = ["atchley_cosine", "sim_A", "sim_B", "bf_corr", "rmsd_geo"]
    if any(r.get("esp_sim") is not None for r in results):
        metrics.append("esp_sim")
    if any(r.get("pesto_sim") is not None for r in results):
        metrics.append("pesto_sim")

    # Compute ranks for each metric
    ranks = {}
    for m in metrics:
        vals = []
        for r in results:
            v = r.get(m)
            vals.append(v if v is not None else -999)
        # Rank descending (highest score = rank 1)
        order = np.argsort(vals)[::-1]
        rank_arr = np.zeros(len(vals), dtype=int)
        for rank, idx in enumerate(order, 1):
            rank_arr[idx] = rank
        ranks[m] = rank_arr

    # Print top-20 by each metric
    for m in metrics:
        print(f"\n--- Top 20 by {m} ---")
        sorted_idx = np.argsort(ranks[m])
        print(f"{'Rank':>4}  {'Sequence':<12}  {m:>14}  {'RMSD_geo':>10}  {'Atchley':>10}")
        for pos, idx in enumerate(sorted_idx[:20]):
            r = results[idx]
            val = r.get(m)
            val_str = f"{val:.4f}" if val is not None else "N/A"
            print(f"{pos+1:>4}  {r['sequence']:<12}  {val_str:>14}  "
                  f"{r['rmsd_geo']:.4f}  {r['atchley_cosine']:.4f}")

    # Spearman rank correlation matrix
    from scipy.stats import spearmanr

    print(f"\n{'=' * 80}")
    print("SPEARMAN RANK CORRELATION MATRIX")
    print(f"{'=' * 80}")

    # Header
    print(f"{'':>16}", end="")
    for m in metrics:
        label = m[:12]
        print(f"{label:>14}", end="")
    print()

    for m1 in metrics:
        print(f"{m1[:16]:>16}", end="")
        for m2 in metrics:
            rho, _ = spearmanr(ranks[m1], ranks[m2])
            print(f"{rho:>14.3f}", end="")
        print()

    # ── Outlier detection: who ranks very differently? ──────────────────
    print(f"\n{'=' * 80}")
    print("RANK DISAGREEMENT ANALYSIS (|rank_A - rank_B| > 15)")
    print(f"{'=' * 80}")

    for i, r in enumerate(results):
        disagreements = []
        for m1 in metrics:
            for m2 in metrics:
                if m1 >= m2:
                    continue
                diff = abs(int(ranks[m1][i]) - int(ranks[m2][i]))
                if diff > 15:
                    disagreements.append(
                        f"{m1}=#{ranks[m1][i]} vs {m2}=#{ranks[m2][i]} (Δ={diff})"
                    )
        if disagreements:
            print(f"\n  {r['sequence']}:")
            for d in disagreements:
                print(f"    {d}")

    # ── Save full results to JSON ──────────────────────────────────────
    output = []
    for i, r in enumerate(results):
        entry = dict(r)
        entry.pop("pdb", None)
        for m in metrics:
            entry[f"rank_{m}"] = int(ranks[m][i])
        output.append(entry)

    out_path = STRUCT_DIR.parent / "descriptor_solo_ranks.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nFull results saved to: {out_path}")


if __name__ == "__main__":
    main()
