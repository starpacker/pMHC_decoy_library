"""
TCR-facing surface descriptors for pMHC structures.
====================================================
Computes descriptors on the **solvent-exposed, TCR-accessible surface**
of pMHC complexes — the surface that TCR actually sees — rather than
the peptide-MHC binding groove interface.

Four descriptors:
1. **TCR-facing ESP**     — Coulomb potential at exposed peptide atoms
                            (what the TCR "feels" electrostatically)
2. **rSASA profile**      — Per-position relative solvent-accessible surface area
                            (which positions are exposed to TCR)
3. **Exposed hydrophobicity** — SASA-weighted Kyte-Doolittle profile
                            (hydrophobic patches on TCR-facing surface)
4. **Exposed shape**      — Coordinates of exposed side-chain atoms after
                            MHC superposition (3D shape of what TCR sees)

All pairwise similarity functions return normalized [0, 1] scores.

Usage:
    from decoy_b.tools.tcr_surface_descriptors import compute_tcr_facing_similarity
    sim = compute_tcr_facing_similarity(target_pdb, candidate_pdb)
    # sim.esp_similarity, sim.sasa_similarity, sim.hydrophobicity_similarity,
    # sim.shape_similarity, sim.combined
"""

from __future__ import annotations

import logging
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# Chain convention from tFold / AF3
DEFAULT_PEPTIDE_CHAIN = "P"
DEFAULT_MHC_CHAINS = ("M", "N")

# Maximum peptide length for fixed-length vectors (pad/truncate)
MAX_PEPTIDE_LEN = 15

# Kyte-Doolittle hydrophobicity scale
KYTE_DOOLITTLE = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
    "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
    "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
    "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

# Partial charges for Coulomb proxy (by element)
PARTIAL_CHARGES = {
    "N": -0.4, "O": -0.5, "C": 0.1, "S": -0.1,
}
# Charged residue corrections
CHARGED_ATOMS = {
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5, ("ARG", "NE"): 0.3,
    ("LYS", "NZ"): 0.7,
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("HIS", "ND1"): 0.3, ("HIS", "NE2"): 0.3,
}


# ======================================================================
#  Data Structures
# ======================================================================

class TCRFacingDescriptors:
    """All TCR-facing descriptors for a single pMHC structure."""

    __slots__ = (
        "sasa_profile", "esp_profile", "hydrophobicity_profile",
        "exposed_sidechain_coords", "peptide_length",
    )

    def __init__(
        self,
        sasa_profile: Optional[np.ndarray] = None,
        esp_profile: Optional[np.ndarray] = None,
        hydrophobicity_profile: Optional[np.ndarray] = None,
        exposed_sidechain_coords: Optional[np.ndarray] = None,
        peptide_length: int = 0,
    ):
        self.sasa_profile = sasa_profile
        self.esp_profile = esp_profile
        self.hydrophobicity_profile = hydrophobicity_profile
        self.exposed_sidechain_coords = exposed_sidechain_coords
        self.peptide_length = peptide_length


class TCRFacingSimilarity:
    """Pairwise similarity between target and candidate TCR-facing surfaces."""

    __slots__ = (
        "esp_similarity", "sasa_similarity",
        "hydrophobicity_similarity", "shape_similarity", "combined",
    )

    def __init__(
        self,
        esp_similarity: float = 0.0,
        sasa_similarity: float = 0.0,
        hydrophobicity_similarity: float = 0.0,
        shape_similarity: float = 0.0,
        combined: float = 0.0,
    ):
        self.esp_similarity = esp_similarity
        self.sasa_similarity = sasa_similarity
        self.hydrophobicity_similarity = hydrophobicity_similarity
        self.shape_similarity = shape_similarity
        self.combined = combined

    def to_dict(self) -> dict:
        return {
            "esp_similarity": self.esp_similarity,
            "sasa_similarity": self.sasa_similarity,
            "hydrophobicity_similarity": self.hydrophobicity_similarity,
            "shape_similarity": self.shape_similarity,
            "tcr_combined": self.combined,
        }


# ======================================================================
#  1. Per-Residue SASA (rSASA profile)
# ======================================================================

def compute_sasa_profile(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[np.ndarray]:
    """
    Compute per-position relative SASA (rSASA) for the peptide.

    rSASA(i) = SASA_in_complex(i) / SASA_free(i)
      - 0.0 = fully buried in groove (anchor residue)
      - ~0.4-0.8 = exposed to solvent (TCR-contact residue)

    Returns a vector of length = peptide residue count.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError:
        log.warning("BioPython SASA not available")
        return None

    try:
        parser = PDBParser(QUIET=True)
        sr = ShrakeRupley()

        # --- Complex SASA (all chains present) ---
        struct_complex = parser.get_structure("complex", pdb_path)
        model_complex = struct_complex[0]
        sr.compute(model_complex, level="A")

        pep_chain_complex = None
        for ch in model_complex:
            if ch.id == peptide_chain:
                pep_chain_complex = ch
                break
        if pep_chain_complex is None:
            log.warning("Peptide chain %s not found in %s", peptide_chain, pdb_path)
            return None

        complex_sasa_per_res = []
        for res in pep_chain_complex:
            if res.id[0] != " ":
                continue
            complex_sasa_per_res.append(sum(a.sasa for a in res))

        # --- Free peptide SASA (peptide alone) ---
        pdb_lines = Path(pdb_path).read_text().splitlines()
        pep_lines = [
            l for l in pdb_lines
            if l.startswith(("ATOM", "HETATM")) and len(l) > 21
            and l[21] == peptide_chain
        ]
        if not pep_lines:
            return None

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".pdb", delete=False
        ) as f:
            f.write("\n".join(pep_lines) + "\nEND\n")
            tmp_path = f.name

        try:
            struct_free = parser.get_structure("free", tmp_path)
            model_free = struct_free[0]
            sr.compute(model_free, level="A")

            free_sasa_per_res = []
            for ch in model_free:
                for res in ch:
                    if res.id[0] != " ":
                        continue
                    free_sasa_per_res.append(sum(a.sasa for a in res))
        finally:
            Path(tmp_path).unlink(missing_ok=True)

        if len(complex_sasa_per_res) != len(free_sasa_per_res):
            log.warning("SASA residue count mismatch: complex=%d free=%d",
                        len(complex_sasa_per_res), len(free_sasa_per_res))
            return None

        rsasa = np.array([
            c / max(f, 0.01)
            for c, f in zip(complex_sasa_per_res, free_sasa_per_res)
        ])

        log.debug(
            "rSASA %s: %s",
            Path(pdb_path).name,
            " ".join(f"{v:.2f}" for v in rsasa),
        )
        return rsasa

    except Exception as exc:
        log.warning("SASA computation failed for %s: %s", pdb_path, exc)
        return None


def compute_sasa_similarity(
    sasa_target: np.ndarray,
    sasa_candidate: np.ndarray,
) -> float:
    """
    Cosine similarity between two rSASA profiles.

    Captures whether the same peptide positions are exposed to TCR.
    """
    # Align lengths (pad shorter with 0)
    n = max(len(sasa_target), len(sasa_candidate))
    a = np.zeros(n)
    b = np.zeros(n)
    a[:len(sasa_target)] = sasa_target
    b[:len(sasa_candidate)] = sasa_candidate

    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.clip(np.dot(a, b) / (na * nb), 0.0, 1.0))


# ======================================================================
#  2. TCR-Facing Electrostatic Potential (Coulomb proxy)
# ======================================================================

def compute_tcr_facing_esp(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    cutoff: float = 8.0,
) -> Optional[np.ndarray]:
    """
    Compute electrostatic potential at TCR-exposed peptide atoms.

    For each peptide atom with SASA > 0 in the pMHC complex, compute:
        ESP_i = sum_j(q_j / r_ij) for all non-hydrogen atoms j within cutoff

    Then average per residue position to get a fixed-length profile.
    This is what TCR "feels" electrostatically at each peptide position.

    Returns a per-position ESP vector, or None.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError:
        return None

    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("esp", pdb_path)
        model = struct[0]

        # Compute atom-level SASA to identify exposed atoms
        sr = ShrakeRupley()
        sr.compute(model, level="A")

        mhc_set = set(mhc_chains)

        # Collect all non-H atoms with their coordinates and charges
        all_atoms = []  # (coord, charge)
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
                    q = CHARGED_ATOMS.get((resname, aname))
                    if q is None:
                        q = PARTIAL_CHARGES.get(atom.element, 0.0)
                    all_atoms.append((coord, q))

        if not all_atoms:
            return None

        all_coords = np.array([a[0] for a in all_atoms])
        all_charges = np.array([a[1] for a in all_atoms])

        # Find peptide chain and compute per-position ESP at exposed atoms
        pep_chain = None
        for ch in model:
            if ch.id == peptide_chain:
                pep_chain = ch
                break
        if pep_chain is None:
            return None

        per_position_esp = []
        for res in pep_chain:
            if res.id[0] != " ":
                continue

            # Collect exposed atoms for this residue
            exposed_coords = []
            for atom in res:
                if atom.element == "H":
                    continue
                if atom.sasa > 0.1:  # solvent-exposed
                    exposed_coords.append(atom.get_vector().get_array())

            if not exposed_coords:
                per_position_esp.append(0.0)
                continue

            # For each exposed peptide atom, compute Coulomb sum from all other atoms
            esp_sum = 0.0
            n_atoms = 0
            for ec in exposed_coords:
                dists = np.linalg.norm(all_coords - ec, axis=1)
                mask = (dists < cutoff) & (dists > 0.5)  # exclude self
                if np.any(mask):
                    r = np.maximum(dists[mask], 0.5)
                    esp_sum += float(np.sum(all_charges[mask] / r))
                    n_atoms += 1

            if n_atoms > 0:
                per_position_esp.append(esp_sum / n_atoms)
            else:
                per_position_esp.append(0.0)

        esp_arr = np.array(per_position_esp)

        log.debug(
            "TCR-facing ESP %s: %s",
            Path(pdb_path).name,
            " ".join(f"{v:.3f}" for v in esp_arr),
        )
        return esp_arr

    except Exception as exc:
        log.warning("TCR-facing ESP failed for %s: %s", pdb_path, exc)
        return None


def compute_esp_hodgkin_similarity(
    esp_target: np.ndarray,
    esp_candidate: np.ndarray,
) -> float:
    """
    Hodgkin Similarity Index between two ESP profiles.

    SI = 2 * dot(A, B) / (dot(A, A) + dot(B, B))

    Range: [-1, 1]. 1.0 = identical, -1.0 = opposite charge.
    We map to [0, 1] via (SI + 1) / 2 for consistency with other descriptors.

    This is the standard metric in PIPSA/MatchTope literature.
    """
    # Align lengths
    n = max(len(esp_target), len(esp_candidate))
    a = np.zeros(n)
    b = np.zeros(n)
    a[:len(esp_target)] = esp_target
    b[:len(esp_candidate)] = esp_candidate

    denom = float(np.dot(a, a) + np.dot(b, b))
    if denom < 1e-12:
        return 1.0  # both zero = identical (no charge)

    hodgkin_si = 2.0 * float(np.dot(a, b)) / denom
    # Map [-1, 1] → [0, 1]
    return float(np.clip((hodgkin_si + 1.0) / 2.0, 0.0, 1.0))


# ======================================================================
#  3. Exposed Hydrophobicity Profile
# ======================================================================

def compute_exposed_hydrophobicity(
    pdb_path: str,
    sasa_profile: Optional[np.ndarray] = None,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[np.ndarray]:
    """
    SASA-weighted Kyte-Doolittle hydrophobicity at each peptide position.

    hydro(i) = rSASA(i) * KD(residue_i)

    Captures the hydrophobic character of TCR-facing residues.
    Buried residues (rSASA ≈ 0) contribute nothing — they don't face TCR.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        return None

    try:
        if sasa_profile is None:
            sasa_profile = compute_sasa_profile(pdb_path, peptide_chain, mhc_chains)
        if sasa_profile is None:
            return None

        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("hydro", pdb_path)
        model = struct[0]

        pep_chain = None
        for ch in model:
            if ch.id == peptide_chain:
                pep_chain = ch
                break
        if pep_chain is None:
            return None

        residues = [r for r in pep_chain if r.id[0] == " "]
        if len(residues) != len(sasa_profile):
            return None

        hydro_profile = np.array([
            sasa_profile[i] * KYTE_DOOLITTLE.get(res.get_resname(), 0.0)
            for i, res in enumerate(residues)
        ])

        log.debug(
            "Exposed hydrophobicity %s: %s",
            Path(pdb_path).name,
            " ".join(f"{v:.2f}" for v in hydro_profile),
        )
        return hydro_profile

    except Exception as exc:
        log.warning("Hydrophobicity computation failed for %s: %s", pdb_path, exc)
        return None


def compute_hydrophobicity_similarity(
    hydro_target: np.ndarray,
    hydro_candidate: np.ndarray,
) -> float:
    """
    Cosine similarity between SASA-weighted hydrophobicity profiles.

    Captures whether the same hydrophobic/hydrophilic patches face the TCR.
    """
    n = max(len(hydro_target), len(hydro_candidate))
    a = np.zeros(n)
    b = np.zeros(n)
    a[:len(hydro_target)] = hydro_target
    b[:len(hydro_candidate)] = hydro_candidate

    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    return float(np.clip(np.dot(a, b) / (na * nb), 0.0, 1.0))


# ======================================================================
#  4. Exposed Side-Chain Shape (post-superposition)
# ======================================================================

def compute_exposed_sidechain_coords(
    pdb_path: str,
    sasa_profile: Optional[np.ndarray] = None,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    sasa_threshold: float = 0.15,
) -> Optional[np.ndarray]:
    """
    Extract 3D coordinates of exposed peptide side-chain heavy atoms.

    For each residue with rSASA >= sasa_threshold (TCR-facing):
      - Collect all heavy atoms with SASA > 0 in the complex
      - Store their coordinates

    Returns per-position centroid coordinates (N x 3) for exposed residues.
    These can be compared after MHC superposition to capture side-chain
    orientation differences that backbone RMSD misses.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError:
        return None

    try:
        if sasa_profile is None:
            sasa_profile = compute_sasa_profile(pdb_path, peptide_chain, mhc_chains)
        if sasa_profile is None:
            return None

        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("shape", pdb_path)
        model = struct[0]

        sr = ShrakeRupley()
        sr.compute(model, level="A")

        pep_chain_obj = None
        for ch in model:
            if ch.id == peptide_chain:
                pep_chain_obj = ch
                break
        if pep_chain_obj is None:
            return None

        residues = [r for r in pep_chain_obj if r.id[0] == " "]
        if len(residues) != len(sasa_profile):
            return None

        centroids = []
        for i, res in enumerate(residues):
            if sasa_profile[i] < sasa_threshold:
                continue  # skip buried residues

            # Collect coordinates of exposed heavy atoms
            exposed_coords = []
            for atom in res:
                if atom.element == "H":
                    continue
                if atom.sasa > 0.1:
                    exposed_coords.append(atom.get_vector().get_array())

            if exposed_coords:
                centroid = np.mean(exposed_coords, axis=0)
                centroids.append(centroid)

        if not centroids:
            return None

        return np.array(centroids)

    except Exception as exc:
        log.warning("Exposed shape computation failed for %s: %s", pdb_path, exc)
        return None


def compute_shape_similarity(
    coords_target: np.ndarray,
    coords_candidate: np.ndarray,
    max_rmsd: float = 3.0,
) -> float:
    """
    Shape similarity based on exposed side-chain centroid RMSD.

    Both coordinate sets should be in MHC-superposed frames (done by
    the dual-superposition RMSD in scanner.py Method A).

    Uses Procrustes alignment on the exposed centroids to handle
    minor frame differences, then computes RMSD.

    Returns 1.0 - RMSD/max_rmsd, clamped to [0, 1].
    """
    n_t = len(coords_target)
    n_c = len(coords_candidate)

    if n_t == 0 or n_c == 0:
        return 0.0

    # If different number of exposed residues, use smaller set (core matching)
    n = min(n_t, n_c)
    # Center-align (take central n residues from each)
    t_start = (n_t - n) // 2
    c_start = (n_c - n) // 2
    ct = coords_target[t_start:t_start + n].copy()
    cc = coords_candidate[c_start:c_start + n].copy()

    # Center both
    ct -= ct.mean(axis=0)
    cc -= cc.mean(axis=0)

    # Optimal rotation via SVD (Kabsch algorithm)
    H = cc.T @ ct
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign_matrix @ U.T
    cc_rotated = (R @ cc.T).T

    # RMSD
    diff = ct - cc_rotated
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))

    sim = max(0.0, 1.0 - rmsd / max_rmsd)

    log.debug("Exposed shape RMSD=%.3f (n=%d exposed residues) → sim=%.3f",
              rmsd, n, sim)
    return sim


# ======================================================================
#  Composite: Compute All & Compare
# ======================================================================

def compute_tcr_facing_descriptors(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> TCRFacingDescriptors:
    """
    Compute all TCR-facing descriptors for a single pMHC structure.
    """
    sasa = compute_sasa_profile(pdb_path, peptide_chain, mhc_chains)
    esp = compute_tcr_facing_esp(pdb_path, peptide_chain, mhc_chains)
    hydro = compute_exposed_hydrophobicity(pdb_path, sasa, peptide_chain, mhc_chains)
    shape = compute_exposed_sidechain_coords(pdb_path, sasa, peptide_chain, mhc_chains)

    pep_len = len(sasa) if sasa is not None else 0

    return TCRFacingDescriptors(
        sasa_profile=sasa,
        esp_profile=esp,
        hydrophobicity_profile=hydro,
        exposed_sidechain_coords=shape,
        peptide_length=pep_len,
    )


def compute_tcr_facing_similarity(
    target_pdb: str,
    candidate_pdb: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    weights: Optional[Dict[str, float]] = None,
) -> TCRFacingSimilarity:
    """
    Compute pairwise TCR-facing surface similarity between two pMHC complexes.

    Default weights (4-descriptor):
        {"esp": 0.30, "shape": 0.30, "sasa": 0.20, "hydrophobicity": 0.20}
    - ESP (0.30):  electrostatic complementarity — primary TCR recognition driver
    - Shape (0.30): exposed side-chain geometry — determines TCR fit
    - SASA (0.20):  exposure pattern — which positions face TCR
    - Hydrophobicity (0.20): surface character — drives binding energetics

    If a descriptor fails, its weight is redistributed to successful ones.
    """
    if weights is None:
        weights = {
            "esp": 0.30, "shape": 0.30,
            "sasa": 0.20, "hydrophobicity": 0.20,
        }

    desc_t = compute_tcr_facing_descriptors(target_pdb, peptide_chain, mhc_chains)
    desc_c = compute_tcr_facing_descriptors(candidate_pdb, peptide_chain, mhc_chains)

    # ESP similarity (Hodgkin SI)
    esp_sim = 0.0
    esp_ok = False
    if desc_t.esp_profile is not None and desc_c.esp_profile is not None:
        esp_sim = compute_esp_hodgkin_similarity(desc_t.esp_profile, desc_c.esp_profile)
        esp_ok = True

    # SASA profile similarity
    sasa_sim = 0.0
    sasa_ok = False
    if desc_t.sasa_profile is not None and desc_c.sasa_profile is not None:
        sasa_sim = compute_sasa_similarity(desc_t.sasa_profile, desc_c.sasa_profile)
        sasa_ok = True

    # Hydrophobicity similarity
    hydro_sim = 0.0
    hydro_ok = False
    if (desc_t.hydrophobicity_profile is not None
            and desc_c.hydrophobicity_profile is not None):
        hydro_sim = compute_hydrophobicity_similarity(
            desc_t.hydrophobicity_profile, desc_c.hydrophobicity_profile,
        )
        hydro_ok = True

    # Exposed shape similarity
    shape_sim = 0.0
    shape_ok = False
    if (desc_t.exposed_sidechain_coords is not None
            and desc_c.exposed_sidechain_coords is not None):
        shape_sim = compute_shape_similarity(
            desc_t.exposed_sidechain_coords, desc_c.exposed_sidechain_coords,
        )
        shape_ok = True

    # Adaptive weighting
    active = {}
    if esp_ok:
        active["esp"] = weights.get("esp", 0.30)
    if sasa_ok:
        active["sasa"] = weights.get("sasa", 0.20)
    if hydro_ok:
        active["hydrophobicity"] = weights.get("hydrophobicity", 0.20)
    if shape_ok:
        active["shape"] = weights.get("shape", 0.30)

    if not active:
        return TCRFacingSimilarity(
            esp_similarity=esp_sim,
            sasa_similarity=sasa_sim,
            hydrophobicity_similarity=hydro_sim,
            shape_similarity=shape_sim,
            combined=0.0,
        )

    total_w = sum(active.values())
    norm = {k: v / total_w for k, v in active.items()}

    combined = (
        norm.get("esp", 0.0) * esp_sim
        + norm.get("sasa", 0.0) * sasa_sim
        + norm.get("hydrophobicity", 0.0) * hydro_sim
        + norm.get("shape", 0.0) * shape_sim
    )

    log.debug(
        "TCR-facing similarity: ESP=%.3f%s SASA=%.3f%s Hydro=%.3f%s "
        "Shape=%.3f%s → combined=%.3f",
        esp_sim, "" if esp_ok else "(N/A)",
        sasa_sim, "" if sasa_ok else "(N/A)",
        hydro_sim, "" if hydro_ok else "(N/A)",
        shape_sim, "" if shape_ok else "(N/A)",
        combined,
    )

    return TCRFacingSimilarity(
        esp_similarity=esp_sim,
        sasa_similarity=sasa_sim,
        hydrophobicity_similarity=hydro_sim,
        shape_similarity=shape_sim,
        combined=combined,
    )
