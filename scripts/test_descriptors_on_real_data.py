"""
Test interface descriptors on real GILGFVFTL pMHC structures.

Since PLIP/FreeSASA/PRODIGY/APBS/PeSTo require C extensions or GPU libraries
not available on Windows, this script implements pure-Python equivalents
using BioPython:

1. Contact Fingerprint (PLIP proxy):
   - H-bond proxy: N/O pairs < 3.5 A across interface
   - Hydrophobic: C-C pairs < 4.5 A
   - Salt bridge: charged atom pairs < 5.5 A
   - Aromatic: aromatic ring C pairs < 5.5 A

2. BSA proxy (FreeSASA proxy):
   - Counts atoms "buried" at interface (< 4.0 A from partner chain)
   - Normalizes by total atom count

3. Contact Affinity (PRODIGY proxy):
   - Classifies interfacial contacts by residue polarity
   - Charged-Charged, Charged-Polar, Charged-Apolar, Polar-Polar, etc.
   - Uses PRODIGY's linear model coefficients

4. Coulomb ESP (APBS proxy):
   - Computes electrostatic potential at interface peptide atoms
   - Sum of q_j/r_ij from nearby MHC atoms using AMBER-like charges

5. Interface Propensity (PeSTo proxy):
   - Per-residue contact density, physicochemical type, proximity
   - Sigmoid-based interface probability estimate

All run on the 51 tFold PDB structures from GILGFVFTL_summary.
"""

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

# ── Config ──────────────────────────────────────────────────────────────

PDB_DIR = Path(r"C:\Users\30670\AppData\Local\Temp\pMHC_decoy_library"
               r"\data\GILGFVFTL_summary\Decoy_B\3D_structures")
RESULTS_JSON = Path(r"C:\Users\30670\AppData\Local\Temp\pMHC_decoy_library"
                    r"\data\GILGFVFTL_summary\Decoy_B\final_ranked_decoys.json")
TARGET_SEQ = "GILGFVFTL"
TARGET_PDB = PDB_DIR / f"pmhc_{TARGET_SEQ}_HLA-A0201.pdb"

# Residue polarity classification (PRODIGY-style)
CHARGED = {"ARG", "LYS", "ASP", "GLU", "HIS"}
POLAR = {"SER", "THR", "ASN", "GLN", "TYR", "CYS", "TRP"}
APOLAR = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "GLY"}

AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}

parser = PDBParser(QUIET=True)


# ═══════════════════════════════════════════════════════════════════════
#  1. Contact Fingerprint (PLIP proxy)
# ═══════════════════════════════════════════════════════════════════════

def compute_contact_fingerprint(struct):
    """Extract non-covalent interaction counts at peptide-MHC interface."""
    pep_atoms = []
    mhc_atoms = []
    for chain in struct.get_chains():
        if chain.id == "P":
            pep_atoms = list(chain.get_atoms())
        elif chain.id == "M":
            mhc_atoms = list(chain.get_atoms())

    if not pep_atoms or not mhc_atoms:
        return None

    pep_coords = np.array([a.get_vector().get_array() for a in pep_atoms])
    mhc_coords = np.array([a.get_vector().get_array() for a in mhc_atoms])

    hbond = 0
    hydrophobic = 0
    salt_bridge = 0
    aromatic = 0
    total_contacts = 0

    # Per-residue interaction map (peptide residue -> types)
    residue_map = defaultdict(list)

    for i, pa in enumerate(pep_atoms):
        pa_coord = pep_coords[i]
        dists = np.linalg.norm(mhc_coords - pa_coord, axis=1)
        pa_elem = pa.element.strip() if pa.element else pa.get_name()[0]
        pa_resname = pa.get_parent().get_resname()
        pa_resnum = pa.get_parent().get_id()[1]

        for j in np.where(dists < 5.5)[0]:
            d = dists[j]
            ma = mhc_atoms[j]
            ma_elem = ma.element.strip() if ma.element else ma.get_name()[0]
            ma_resname = ma.get_parent().get_resname()

            if d < 4.5:
                total_contacts += 1

            # H-bond proxy: N/O pairs < 3.5 A
            if d < 3.5 and pa_elem in ("N", "O") and ma_elem in ("N", "O"):
                hbond += 1
                residue_map[pa_resnum].append("hbond")

            # Hydrophobic: C-C pairs < 4.5 A
            elif d < 4.5 and pa_elem == "C" and ma_elem == "C":
                hydrophobic += 1
                residue_map[pa_resnum].append("hydrophobic")

            # Salt bridge: charged residue atom pairs < 5.5 A
            elif d < 5.5:
                if (pa_resname in CHARGED and ma_resname in CHARGED
                        and pa_elem in ("N", "O") and ma_elem in ("N", "O")):
                    salt_bridge += 1
                    residue_map[pa_resnum].append("salt_bridge")

            # Aromatic: aromatic C pairs < 5.5 A
            if (d < 5.5 and pa_resname in AROMATIC_RES
                    and ma_resname in AROMATIC_RES
                    and pa_elem == "C" and ma_elem == "C"):
                aromatic += 1
                residue_map[pa_resnum].append("aromatic")

    return {
        "hbond": hbond,
        "hydrophobic": hydrophobic,
        "salt_bridge": salt_bridge,
        "aromatic": aromatic,
        "total_contacts": total_contacts,
        "residue_map": dict(residue_map),
    }


def tanimoto(vec_a, vec_b):
    """Generalized Tanimoto on count vectors."""
    a = np.array(vec_a, dtype=float)
    b = np.array(vec_b, dtype=float)
    dot = float(np.dot(a, b))
    denom = float(np.dot(a, a) + np.dot(b, b) - dot)
    if denom < 1e-12:
        return 1.0 if np.dot(a, a) < 1e-12 else 0.0
    return max(0.0, min(1.0, dot / denom))


# ═══════════════════════════════════════════════════════════════════════
#  2. BSA proxy (FreeSASA proxy)
# ═══════════════════════════════════════════════════════════════════════

def compute_bsa_proxy(struct, cutoff=4.0):
    """
    Count atoms buried at interface as BSA proxy.
    An atom is 'buried' if it has a partner atom from the other chain within cutoff.
    Returns: n_buried_peptide, n_buried_mhc, n_total_peptide, n_total_mhc
    """
    pep_atoms = []
    mhc_atoms = []
    for chain in struct.get_chains():
        if chain.id == "P":
            pep_atoms = list(chain.get_atoms())
        elif chain.id == "M":
            mhc_atoms = list(chain.get_atoms())

    if not pep_atoms or not mhc_atoms:
        return None

    pep_coords = np.array([a.get_vector().get_array() for a in pep_atoms])
    mhc_coords = np.array([a.get_vector().get_array() for a in mhc_atoms])

    # Peptide atoms buried
    n_pep_buried = 0
    for i in range(len(pep_coords)):
        dists = np.linalg.norm(mhc_coords - pep_coords[i], axis=1)
        if np.any(dists < cutoff):
            n_pep_buried += 1

    # MHC atoms buried
    n_mhc_buried = 0
    for j in range(len(mhc_coords)):
        dists = np.linalg.norm(pep_coords - mhc_coords[j], axis=1)
        if np.any(dists < cutoff):
            n_mhc_buried += 1

    bsa_proxy = n_pep_buried + n_mhc_buried
    return {
        "bsa_proxy": bsa_proxy,
        "pep_buried": n_pep_buried,
        "mhc_buried": n_mhc_buried,
        "pep_total": len(pep_atoms),
        "mhc_total": len(mhc_atoms),
        "pep_burial_frac": n_pep_buried / max(len(pep_atoms), 1),
    }


def bsa_similarity(bsa_a, bsa_b, max_diff=150.0):
    """BSA proxy similarity."""
    diff = abs(bsa_a - bsa_b)
    return max(0.0, 1.0 - diff / max_diff)


# ═══════════════════════════════════════════════════════════════════════
#  3. Contact Affinity (PRODIGY proxy)
# ═══════════════════════════════════════════════════════════════════════

def classify_residue(resname):
    if resname in CHARGED:
        return "C"
    elif resname in POLAR:
        return "P"
    elif resname in APOLAR:
        return "A"
    return "A"


def compute_prodigy_proxy(struct, cutoff=8.5):
    """
    PRODIGY-like contact classification at peptide-MHC interface.
    Counts contacts by residue polarity pair type.
    Uses simplified PRODIGY linear model to estimate dG.
    """
    pep_residues = {}
    mhc_residues = {}
    for chain in struct.get_chains():
        if chain.id == "P":
            for res in chain.get_residues():
                if res.get_id()[0] == " ":
                    ca = res["CA"] if "CA" in res else None
                    if ca:
                        pep_residues[res.get_id()[1]] = {
                            "resname": res.get_resname(),
                            "class": classify_residue(res.get_resname()),
                            "coord": ca.get_vector().get_array(),
                        }
        elif chain.id == "M":
            for res in chain.get_residues():
                if res.get_id()[0] == " ":
                    ca = res["CA"] if "CA" in res else None
                    if ca:
                        mhc_residues[res.get_id()[1]] = {
                            "resname": res.get_resname(),
                            "class": classify_residue(res.get_resname()),
                            "coord": ca.get_vector().get_array(),
                        }

    if not pep_residues or not mhc_residues:
        return None

    # Contact matrix by polarity class
    contact_types = defaultdict(int)
    total_contacts = 0

    for p_id, p_info in pep_residues.items():
        for m_id, m_info in mhc_residues.items():
            d = np.linalg.norm(p_info["coord"] - m_info["coord"])
            if d < cutoff:
                pair = tuple(sorted([p_info["class"], m_info["class"]]))
                contact_types[pair] += 1
                total_contacts += 1

    # Simplified PRODIGY linear model (approximate coefficients)
    # Real PRODIGY: dG = -0.09459*CC + -0.10007*CP + 0.19577*CA
    #              + -0.22671*PP + 0.18681*PA + 0.13810*AA + ...
    # We use simplified version
    cc = contact_types.get(("C", "C"), 0)
    cp = contact_types.get(("C", "P"), 0)
    ca = contact_types.get(("A", "C"), 0)
    pp = contact_types.get(("P", "P"), 0)
    pa = contact_types.get(("A", "P"), 0)
    aa = contact_types.get(("A", "A"), 0)

    # Simplified dG estimate (kcal/mol, more negative = stronger)
    dg_est = (-0.09 * cc - 0.10 * cp + 0.20 * ca
              - 0.23 * pp + 0.19 * pa + 0.14 * aa
              - 15.0)  # offset to get into typical pMHC range

    return {
        "total_contacts": total_contacts,
        "CC": cc, "CP": cp, "CA": ca,
        "PP": pp, "PA": pa, "AA": aa,
        "dg_estimate": dg_est,
        "contact_types": dict(contact_types),
    }


def prodigy_similarity(dg_a, dg_b, max_diff=5.0):
    diff = abs(dg_a - dg_b)
    return max(0.0, 1.0 - diff / max_diff)


# ═══════════════════════════════════════════════════════════════════════
#  4. Coulomb ESP proxy (APBS proxy)
# ═══════════════════════════════════════════════════════════════════════

# Simplified partial charges by heavy-atom element
_PARTIAL_CHARGES = {"N": -0.4, "O": -0.5, "C": 0.1, "S": -0.1}
_CHARGED_ATOMS = {
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5, ("ARG", "NE"): 0.3,
    ("LYS", "NZ"): 0.7,
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("HIS", "ND1"): 0.3, ("HIS", "NE2"): 0.3,
}


def compute_esp_proxy(struct, cutoff=5.0):
    """Coulomb electrostatic potential at peptide interface atoms."""
    model = list(struct.get_models())[0]
    pep_atoms = []  # (coord, charge, resnum)
    mhc_atoms = []  # (coord, charge)

    for chain in model:
        for res in chain:
            if res.id[0] != " ":
                continue
            resname = res.get_resname()
            for atom in res:
                if atom.element == "H":
                    continue
                coord = atom.get_vector().get_array()
                aname = atom.get_name()
                q = _CHARGED_ATOMS.get((resname, aname))
                if q is None:
                    q = _PARTIAL_CHARGES.get(atom.element, 0.0)

                if chain.id == "P":
                    pep_atoms.append((coord, q, res.id[1]))
                elif chain.id in ("M", "N"):
                    mhc_atoms.append((coord, q))

    if not pep_atoms or not mhc_atoms:
        return None

    mhc_coords = np.array([a[0] for a in mhc_atoms])
    mhc_charges = np.array([a[1] for a in mhc_atoms])

    esp_values = []
    for pc, pq, resnum in pep_atoms:
        dists = np.linalg.norm(mhc_coords - pc, axis=1)
        mask = dists < cutoff
        if not np.any(mask):
            continue
        r = np.maximum(dists[mask], 0.5)
        esp = float(np.sum(mhc_charges[mask] / r))
        esp_values.append(esp)

    return np.array(esp_values) if esp_values else None


def esp_similarity(esp_a, esp_b):
    """Pearson correlation mapped to [0, 1]."""
    if esp_a is None or esp_b is None:
        return 0.0
    n = min(len(esp_a), len(esp_b))
    if n < 3:
        return 0.0
    a, b = esp_a[:n], esp_b[:n]
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 1.0 if (np.std(a) < 1e-12 and np.std(b) < 1e-12) else 0.0
    corr = float(np.corrcoef(a, b)[0, 1])
    return max(0.0, min(1.0, (corr + 1.0) / 2.0))


# ═══════════════════════════════════════════════════════════════════════
#  5. Interface Propensity proxy (PeSTo proxy)
# ═══════════════════════════════════════════════════════════════════════

_RESIDUE_PROPS = {
    "ALA": [1, 0, 0, 0, 0], "VAL": [1, 0, 0, 0, 0], "LEU": [1, 0, 0, 0, 0],
    "ILE": [1, 0, 0, 0, 0], "MET": [1, 0, 0, 0, 0], "PRO": [1, 0, 0, 0, 0],
    "GLY": [0, 0, 0, 0, 0],
    "SER": [0, 1, 0, 0, 0], "THR": [0, 1, 0, 0, 0], "CYS": [0, 1, 0, 0, 0],
    "ASN": [0, 1, 0, 0, 0], "GLN": [0, 1, 0, 0, 0],
    "ARG": [0, 0, 1, 0, 0], "LYS": [0, 0, 1, 0, 0], "HIS": [0, 0, 1, 0, 0],
    "ASP": [0, 0, 0, 1, 0], "GLU": [0, 0, 0, 1, 0],
    "PHE": [1, 0, 0, 0, 1], "TYR": [0, 1, 0, 0, 1], "TRP": [1, 0, 0, 0, 1],
}
MAX_CONTACTS = 20.0


def compute_pesto_proxy(struct, contact_cutoff=5.0):
    """Contact-based interface propensity vector for peptide residues."""
    model = list(struct.get_models())[0]
    mhc_coords = []
    for chain in model:
        if chain.id in ("M", "N"):
            for res in chain:
                if res.id[0] != " ":
                    continue
                for atom in res:
                    if atom.element != "H":
                        mhc_coords.append(atom.get_vector().get_array())

    if not mhc_coords:
        return None
    mhc_arr = np.array(mhc_coords)

    pep_chain = None
    for chain in model:
        if chain.id == "P":
            pep_chain = chain
            break
    if pep_chain is None:
        return None

    residue_features = []
    for res in pep_chain:
        if res.id[0] != " ":
            continue
        resname = res.get_resname()
        props = _RESIDUE_PROPS.get(resname, [0, 0, 0, 0, 0])

        res_atoms = [a for a in res if a.element != "H"]
        if not res_atoms:
            residue_features.append(props + [0.0, 0.0, 0.0])
            continue

        res_coords = np.array([a.get_vector().get_array() for a in res_atoms])
        n_contacts = 0
        min_dist = float("inf")
        for rc in res_coords:
            dists = np.linalg.norm(mhc_arr - rc, axis=1)
            n_contacts += int(np.sum(dists < contact_cutoff))
            min_dist = min(min_dist, float(dists.min()))

        contact_density = min(n_contacts / MAX_CONTACTS, 1.0)
        interface_prob = 1.0 / (1.0 + np.exp(-5.0 * (contact_density - 0.3)))
        proximity = max(0.0, 1.0 - min_dist / contact_cutoff)

        residue_features.append(props + [contact_density, interface_prob, proximity])

    if not residue_features:
        return None

    feat_matrix = np.array(residue_features, dtype=float)
    mean_feat = feat_matrix.mean(axis=0)
    probs = feat_matrix[:, 6]  # interface_prob column

    # Fixed-length descriptor
    max_len = 15
    padded_probs = np.zeros(max_len)
    n = min(len(probs), max_len)
    padded_probs[:n] = probs[:n]
    padded_contacts = np.zeros(max_len)
    padded_contacts[:n] = feat_matrix[:n, 5]
    padded_proximity = np.zeros(max_len)
    padded_proximity[:n] = feat_matrix[:n, 7]

    return np.concatenate([padded_probs, padded_contacts, padded_proximity, mean_feat])


def pesto_similarity(emb_a, emb_b):
    """Cosine similarity mapped to [0, 1]."""
    if emb_a is None or emb_b is None:
        return 0.0
    norm_a = np.linalg.norm(emb_a)
    norm_b = np.linalg.norm(emb_b)
    if norm_a < 1e-12 or norm_b < 1e-12:
        return 1.0 if (norm_a < 1e-12 and norm_b < 1e-12) else 0.0
    cosine = float(np.dot(emb_a, emb_b) / (norm_a * norm_b))
    return max(0.0, min(1.0, (cosine + 1.0) / 2.0))


# ═══════════════════════════════════════════════════════════════════════
#  6. RMSD-based score (existing baseline)
# ═══════════════════════════════════════════════════════════════════════

def compute_rmsd_score(target_struct, cand_struct):
    """Existing RMSD + B-factor scoring."""
    from Bio.PDB import Superimposer

    def get_chain_cas(struct, chain_id):
        for chain in struct.get_chains():
            if chain.id == chain_id:
                return [a for a in chain.get_atoms() if a.get_name() == "CA"]
        return []

    tgt_cas = get_chain_cas(target_struct, "P")
    cand_cas = get_chain_cas(cand_struct, "P")

    if len(tgt_cas) != len(cand_cas) or len(tgt_cas) == 0:
        return {"rmsd": None, "bf_corr": 0, "rmsd_sim": 0, "score": 0}

    sup = Superimposer()
    sup.set_atoms(tgt_cas, cand_cas)
    rmsd = sup.rms

    tgt_bf = np.array([a.get_bfactor() for a in tgt_cas])
    cand_bf = np.array([a.get_bfactor() for a in cand_cas])

    bf_corr = 0.0
    if np.std(tgt_bf) > 0 and np.std(cand_bf) > 0:
        bf_corr = float(np.corrcoef(tgt_bf, cand_bf)[0, 1])

    rmsd_sim = max(0.0, 1.0 - rmsd / 3.0) if rmsd < 3.0 else 0.0
    score = max(rmsd_sim, bf_corr)

    return {"rmsd": rmsd, "bf_corr": bf_corr, "rmsd_sim": rmsd_sim, "score": score}


# ═══════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    pdbs = sorted(PDB_DIR.glob("*.pdb"))
    print(f"Found {len(pdbs)} PDB structures")
    print(f"Target: {TARGET_SEQ} ({TARGET_PDB.name})")
    print()

    if not TARGET_PDB.exists():
        print(f"ERROR: Target PDB not found: {TARGET_PDB}")
        return

    target_struct = parser.get_structure("target", str(TARGET_PDB))

    # Compute target descriptors
    t0 = time.time()
    tgt_fp = compute_contact_fingerprint(target_struct)
    tgt_bsa = compute_bsa_proxy(target_struct)
    tgt_prodigy = compute_prodigy_proxy(target_struct)
    tgt_esp = compute_esp_proxy(target_struct)
    tgt_pesto = compute_pesto_proxy(target_struct)
    t_target = time.time() - t0

    print("=" * 90)
    print(f"  TARGET: {TARGET_SEQ}  (computed in {t_target:.3f}s)")
    print("=" * 90)
    print(f"  Contact fingerprint: HB={tgt_fp['hbond']}, Hydro={tgt_fp['hydrophobic']}, "
          f"SaltBr={tgt_fp['salt_bridge']}, Arom={tgt_fp['aromatic']}, Total={tgt_fp['total_contacts']}")
    print(f"  BSA proxy: {tgt_bsa['bsa_proxy']} buried atoms "
          f"(pep={tgt_bsa['pep_buried']}/{tgt_bsa['pep_total']}, "
          f"mhc={tgt_bsa['mhc_buried']}/{tgt_bsa['mhc_total']})")
    print(f"  PRODIGY proxy: dG={tgt_prodigy['dg_estimate']:.2f} kcal/mol, "
          f"{tgt_prodigy['total_contacts']} CA contacts "
          f"(CC={tgt_prodigy['CC']}, CP={tgt_prodigy['CP']}, CA={tgt_prodigy['CA']}, "
          f"PP={tgt_prodigy['PP']}, PA={tgt_prodigy['PA']}, AA={tgt_prodigy['AA']})")
    if tgt_esp is not None:
        print(f"  ESP proxy: {len(tgt_esp)} interface atoms, range [{tgt_esp.min():.3f}, {tgt_esp.max():.3f}]")
    if tgt_pesto is not None:
        probs = tgt_pesto[:15]  # first 15 = padded interface probs
        n_res = int(np.sum(probs > 0))
        print(f"  PeSTo proxy: {n_res} residues, mean interface_prob={probs[probs > 0].mean():.3f}")

    # Per-residue interactions for target
    print(f"\n  Per-residue interactions (peptide):")
    for pos in sorted(tgt_fp["residue_map"].keys()):
        types = tgt_fp["residue_map"][pos]
        from collections import Counter
        tc = Counter(types)
        aa = TARGET_SEQ[pos - 1] if 1 <= pos <= len(TARGET_SEQ) else "?"
        print(f"    p{pos} ({aa}): {dict(tc)}")

    # ── Compare all 9-mer candidates ───────────────────────────────────
    tgt_fp_vec = [tgt_fp["hbond"], tgt_fp["hydrophobic"],
                  tgt_fp["salt_bridge"], tgt_fp["aromatic"]]

    results = []
    for pdb_path in pdbs:
        seq = pdb_path.stem.split("_")[1]
        if seq == TARGET_SEQ:
            continue

        cand_struct = parser.get_structure("c", str(pdb_path))

        # Check same length
        cand_pep_cas = []
        for chain in cand_struct.get_chains():
            if chain.id == "P":
                cand_pep_cas = [a for a in chain.get_atoms() if a.get_name() == "CA"]
        if len(cand_pep_cas) != len(TARGET_SEQ):
            continue

        t0 = time.time()

        # 1. RMSD baseline
        rmsd_res = compute_rmsd_score(target_struct, cand_struct)

        # 2. Contact fingerprint
        cand_fp = compute_contact_fingerprint(cand_struct)
        cand_fp_vec = [cand_fp["hbond"], cand_fp["hydrophobic"],
                       cand_fp["salt_bridge"], cand_fp["aromatic"]]
        fp_tanimoto = tanimoto(tgt_fp_vec, cand_fp_vec)

        # 3. BSA proxy
        cand_bsa = compute_bsa_proxy(cand_struct)
        bsa_sim = bsa_similarity(tgt_bsa["bsa_proxy"], cand_bsa["bsa_proxy"])

        # 4. PRODIGY proxy
        cand_prodigy = compute_prodigy_proxy(cand_struct)
        prod_sim = prodigy_similarity(tgt_prodigy["dg_estimate"],
                                       cand_prodigy["dg_estimate"])

        # 5. ESP proxy (APBS)
        cand_esp = compute_esp_proxy(cand_struct)
        esp_sim = esp_similarity(tgt_esp, cand_esp)

        # 6. PeSTo proxy
        cand_pesto = compute_pesto_proxy(cand_struct)
        pesto_sim = pesto_similarity(tgt_pesto, cand_pesto)

        # 7. Combined interface score (5-descriptor adaptive weighted)
        # Weights: PLIP=0.25, BSA=0.10, PRODIGY=0.20, ESP=0.25, PeSTo=0.20
        active_scores = {}
        active_weights = {}
        for name, sim, w in [
            ("plip", fp_tanimoto, 0.25), ("bsa", bsa_sim, 0.10),
            ("prodigy", prod_sim, 0.20), ("esp", esp_sim, 0.25),
            ("pesto", pesto_sim, 0.20),
        ]:
            if sim > 0 or name in ("plip", "bsa", "prodigy"):  # always include first 3
                active_scores[name] = sim
                active_weights[name] = w

        w_sum = sum(active_weights.values())
        interface_combined = sum(
            (active_weights[k] / w_sum) * active_scores[k]
            for k in active_scores
        )

        # Old 3-descriptor combined (for comparison)
        old_interface = 0.40 * fp_tanimoto + 0.25 * bsa_sim + 0.35 * prod_sim

        # 8. New score: 50/50 blend
        new_score = 0.50 * rmsd_res["score"] + 0.50 * interface_combined

        dt = time.time() - t0

        results.append({
            "seq": seq,
            "rmsd": rmsd_res["rmsd"],
            "bf_corr": rmsd_res["bf_corr"],
            "old_score": rmsd_res["score"],
            "fp_tanimoto": fp_tanimoto,
            "bsa_sim": bsa_sim,
            "prodigy_sim": prod_sim,
            "esp_sim": esp_sim,
            "pesto_sim": pesto_sim,
            "old_interface": old_interface,
            "interface_combined": interface_combined,
            "new_score": new_score,
            "delta": new_score - rmsd_res["score"],
            "time": dt,
        })

    results.sort(key=lambda r: -r["new_score"])

    # ── Print comparison table ─────────────────────────────────────────
    print("\n")
    print("=" * 150)
    print(f"  COMPARISON: {len(results)} candidates (sorted by NEW score, 5-descriptor)")
    print("=" * 150)
    print(f"{'Rank':>4} {'Peptide':<12} {'RMSD':>6} {'OLD':>7} | "
          f"{'FP_Tan':>7} {'BSA':>7} {'PRODY':>7} {'ESP':>7} {'PeSTo':>7} {'IF_CMB':>7} | "
          f"{'NEW':>7} {'Delta':>7}")
    print("-" * 150)

    for i, r in enumerate(results):
        rmsd_s = f"{r['rmsd']:.3f}" if r["rmsd"] else "N/A"
        print(f"{i+1:>4} {r['seq']:<12} {rmsd_s:>6} {r['old_score']:>7.4f} | "
              f"{r['fp_tanimoto']:>7.4f} {r['bsa_sim']:>7.4f} {r['prodigy_sim']:>7.4f} "
              f"{r['esp_sim']:>7.4f} {r['pesto_sim']:>7.4f} "
              f"{r['interface_combined']:>7.4f} | "
              f"{r['new_score']:>7.4f} {r['delta']:>+7.4f}")

    # ── Statistics ─────────────────────────────────────────────────────
    old_scores = [r["old_score"] for r in results]
    new_scores = [r["new_score"] for r in results]
    fp_sims = [r["fp_tanimoto"] for r in results]
    bsa_sims = [r["bsa_sim"] for r in results]
    prod_sims = [r["prodigy_sim"] for r in results]
    esp_sims = [r["esp_sim"] for r in results]
    pesto_sims = [r["pesto_sim"] for r in results]
    old_iface = [r["old_interface"] for r in results]
    deltas = [r["delta"] for r in results]
    times = [r["time"] for r in results]

    print("\n")
    print("=" * 90)
    print("  STATISTICS (5-descriptor)")
    print("=" * 90)
    print(f"  {'Metric':<30} {'Mean':>8} {'Std':>8} {'Min':>8} {'Max':>8}")
    print(f"  {'-'*30} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
    for name, vals in [
        ("Old score (RMSD+BF)", old_scores),
        ("FP Tanimoto", fp_sims),
        ("BSA similarity", bsa_sims),
        ("PRODIGY similarity", prod_sims),
        ("ESP similarity (NEW)", esp_sims),
        ("PeSTo similarity (NEW)", pesto_sims),
        ("Old interface (3-desc)", old_iface),
        ("New interface (5-desc)", [r["interface_combined"] for r in results]),
        ("New score (blended)", new_scores),
        ("Delta (new - old)", deltas),
    ]:
        print(f"  {name:<30} {np.mean(vals):>8.4f} {np.std(vals):>8.4f} "
              f"{min(vals):>8.4f} {max(vals):>8.4f}")

    print(f"\n  Avg time per candidate: {np.mean(times):.3f}s")
    print(f"  Total time: {sum(times):.1f}s for {len(results)} candidates")

    # ── Rank changes ───────────────────────────────────────────────────
    old_ranked = sorted(results, key=lambda r: -r["old_score"])
    new_ranked = sorted(results, key=lambda r: -r["new_score"])

    old_rank = {r["seq"]: i + 1 for i, r in enumerate(old_ranked)}
    new_rank = {r["seq"]: i + 1 for i, r in enumerate(new_ranked)}

    print("\n")
    print("=" * 90)
    print("  RANK CHANGES (top 15)")
    print("=" * 90)
    print(f"  {'Peptide':<12} {'Old_Rank':>9} {'New_Rank':>9} {'Change':>8} {'Old_Score':>10} {'New_Score':>10}")
    print(f"  {'-'*12} {'-'*9} {'-'*9} {'-'*8} {'-'*10} {'-'*10}")

    for r in new_ranked[:15]:
        seq = r["seq"]
        change = old_rank[seq] - new_rank[seq]
        arrow = f"+{change}" if change > 0 else str(change)
        print(f"  {seq:<12} {old_rank[seq]:>9} {new_rank[seq]:>9} {arrow:>8} "
              f"{r['old_score']:>10.4f} {r['new_score']:>10.4f}")

    # Biggest movers
    movers = sorted(results, key=lambda r: abs(old_rank[r["seq"]] - new_rank[r["seq"]]), reverse=True)
    print(f"\n  Biggest rank changes:")
    for r in movers[:5]:
        seq = r["seq"]
        change = old_rank[seq] - new_rank[seq]
        direction = "UP" if change > 0 else "DOWN"
        print(f"    {seq}: {old_rank[seq]} -> {new_rank[seq]} ({direction} {abs(change)} places)")

    # ── Discrimination improvement ─────────────────────────────────────
    print("\n")
    print("=" * 90)
    print("  DISCRIMINATION ANALYSIS")
    print("=" * 90)
    print(f"  Old score std: {np.std(old_scores):.4f}  (higher = better discrimination)")
    print(f"  New score std: {np.std(new_scores):.4f}")
    print(f"  Improvement:   {(np.std(new_scores)/np.std(old_scores) - 1)*100:+.1f}%")
    print()
    print(f"  Old score range: [{min(old_scores):.4f}, {max(old_scores):.4f}] "
          f"(spread = {max(old_scores)-min(old_scores):.4f})")
    print(f"  New score range: [{min(new_scores):.4f}, {max(new_scores):.4f}] "
          f"(spread = {max(new_scores)-min(new_scores):.4f})")

    # Score > 0.9 count
    print(f"\n  Candidates with score > 0.90:")
    print(f"    Old: {sum(1 for s in old_scores if s > 0.9)}/{len(old_scores)}")
    print(f"    New: {sum(1 for s in new_scores if s > 0.9)}/{len(new_scores)}")
    print(f"  Candidates with score > 0.85:")
    print(f"    Old: {sum(1 for s in old_scores if s > 0.85)}/{len(old_scores)}")
    print(f"    New: {sum(1 for s in new_scores if s > 0.85)}/{len(new_scores)}")

    # ── Correlation between descriptors ────────────────────────────────
    print("\n")
    print("=" * 110)
    print("  INTER-DESCRIPTOR CORRELATIONS (6x6)")
    print("=" * 110)
    labels = ["RMSD_score", "FP_Tan", "BSA_sim", "PRODIGY", "ESP", "PeSTo"]
    data = np.array([old_scores, fp_sims, bsa_sims, prod_sims, esp_sims, pesto_sims])
    corr = np.corrcoef(data)
    print(f"  {'':>12}", end="")
    for l in labels:
        print(f" {l:>10}", end="")
    print()
    for i, l in enumerate(labels):
        print(f"  {l:>12}", end="")
        for j in range(len(labels)):
            print(f" {corr[i,j]:>10.4f}", end="")
        print()

    print("\n  Key correlations:")
    pairs = [
        (0, 1, "RMSD vs PLIP"), (0, 3, "RMSD vs PRODIGY"),
        (0, 4, "RMSD vs ESP"), (0, 5, "RMSD vs PeSTo"),
        (4, 1, "ESP vs PLIP"), (4, 3, "ESP vs PRODIGY"),
        (5, 1, "PeSTo vs PLIP"), (5, 4, "PeSTo vs ESP"),
    ]
    for i, j, label in pairs:
        r = corr[i, j]
        if abs(r) < 0.3:
            interp = "LOW (orthogonal)"
        elif abs(r) < 0.7:
            interp = "MODERATE (partial overlap)"
        else:
            interp = "HIGH (redundant)"
        print(f"    {label:<20} r={r:+.3f}  {interp}")


if __name__ == "__main__":
    main()
