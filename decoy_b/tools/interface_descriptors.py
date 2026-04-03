"""
Interface descriptor computation for pMHC structures.

Computes three classes of descriptors from PDB structures:
1. PLIP: non-covalent interaction fingerprint (H-bonds, hydrophobic, salt bridges, pi-stacking)
2. FreeSASA: buried surface area (BSA / SASA) at the peptide-MHC interface
3. PRODIGY: contact-based binding affinity prediction (dG)

All pairwise similarity functions return normalized [0, 1] scores.

Usage:
    from decoy_b.tools.interface_descriptors import compute_interface_similarity
    sim = compute_interface_similarity(target_pdb, candidate_pdb)
    # sim.plip_tanimoto, sim.bsa_similarity, sim.prodigy_similarity, sim.combined
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# tFold / AF3 output chain convention
DEFAULT_PEPTIDE_CHAIN = "P"
DEFAULT_MHC_CHAINS = ("M", "N")


# ======================================================================
#  Data Structures
# ======================================================================

class InteractionFingerprint:
    """Vectorized non-covalent interaction profile from PLIP.

    Stores both global interaction counts and a per-residue interaction map
    for computing Tanimoto similarity between two pMHC interfaces.
    """

    __slots__ = (
        "hbond_count", "hydrophobic_count", "salt_bridge_count",
        "pi_stacking_count", "pi_cation_count", "residue_interactions",
    )

    def __init__(
        self,
        hbond_count: int = 0,
        hydrophobic_count: int = 0,
        salt_bridge_count: int = 0,
        pi_stacking_count: int = 0,
        pi_cation_count: int = 0,
        residue_interactions: Optional[Dict[int, List[str]]] = None,
    ):
        self.hbond_count = hbond_count
        self.hydrophobic_count = hydrophobic_count
        self.salt_bridge_count = salt_bridge_count
        self.pi_stacking_count = pi_stacking_count
        self.pi_cation_count = pi_cation_count
        self.residue_interactions = residue_interactions or {}

    def to_vector(self) -> np.ndarray:
        """Convert to a 5-element count vector for Tanimoto computation."""
        return np.array([
            self.hbond_count,
            self.hydrophobic_count,
            self.salt_bridge_count,
            self.pi_stacking_count,
            self.pi_cation_count,
        ], dtype=float)

    def to_dict(self) -> dict:
        return {
            "hbond_count": self.hbond_count,
            "hydrophobic_count": self.hydrophobic_count,
            "salt_bridge_count": self.salt_bridge_count,
            "pi_stacking_count": self.pi_stacking_count,
            "pi_cation_count": self.pi_cation_count,
        }


class InterfaceDescriptors:
    """All interface descriptors computed for a single pMHC structure."""

    __slots__ = (
        "fingerprint", "bsa_total", "bsa_polar", "bsa_apolar",
        "prodigy_dg", "prodigy_kd",
    )

    def __init__(
        self,
        fingerprint: Optional[InteractionFingerprint] = None,
        bsa_total: Optional[float] = None,
        bsa_polar: Optional[float] = None,
        bsa_apolar: Optional[float] = None,
        prodigy_dg: Optional[float] = None,
        prodigy_kd: Optional[float] = None,
    ):
        self.fingerprint = fingerprint
        self.bsa_total = bsa_total
        self.bsa_polar = bsa_polar
        self.bsa_apolar = bsa_apolar
        self.prodigy_dg = prodigy_dg
        self.prodigy_kd = prodigy_kd


class InterfaceSimilarity:
    """Pairwise similarity between target and candidate interface descriptors."""

    __slots__ = (
        "plip_tanimoto", "bsa_similarity", "prodigy_similarity", "combined",
    )

    def __init__(
        self,
        plip_tanimoto: float = 0.0,
        bsa_similarity: float = 0.0,
        prodigy_similarity: float = 0.0,
        combined: float = 0.0,
    ):
        self.plip_tanimoto = plip_tanimoto
        self.bsa_similarity = bsa_similarity
        self.prodigy_similarity = prodigy_similarity
        self.combined = combined

    def to_dict(self) -> dict:
        return {
            "plip_tanimoto": self.plip_tanimoto,
            "bsa_similarity": self.bsa_similarity,
            "prodigy_similarity": self.prodigy_similarity,
            "combined": self.combined,
        }


# ======================================================================
#  1. PLIP: Non-Covalent Interaction Fingerprint
# ======================================================================

def compute_plip_fingerprint(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[InteractionFingerprint]:
    """
    Run PLIP on a pMHC PDB to extract non-covalent interactions
    at the peptide-MHC interface.

    Returns InteractionFingerprint or None if PLIP is unavailable.
    """
    try:
        from plip.structure.preparation import PDBComplex
    except ImportError:
        log.warning("PLIP not installed. Run: conda install -c conda-forge openbabel && pip install plip")
        return None

    try:
        mol = PDBComplex()
        mol.load_pdb(pdb_path)
        mol.analyze()

        hbond_count = 0
        hydrophobic_count = 0
        salt_bridge_count = 0
        pi_stacking_count = 0
        pi_cation_count = 0
        residue_interactions: Dict[int, List[str]] = {}

        mhc_set = set(mhc_chains)

        for bsid, interaction in mol.interaction_sets.items():
            # Hydrogen bonds
            for hb in interaction.hbonds_pdon + interaction.hbonds_ldon:
                if _crosses_interface(hb, peptide_chain, mhc_set):
                    hbond_count += 1
                    _record_peptide_residue(residue_interactions, hb, peptide_chain, "hbond")

            # Hydrophobic contacts
            for hc in interaction.hydrophobic_contacts:
                if _crosses_interface(hc, peptide_chain, mhc_set):
                    hydrophobic_count += 1
                    _record_peptide_residue(residue_interactions, hc, peptide_chain, "hydrophobic")

            # Salt bridges
            for sb in interaction.saltbridge_lneg + interaction.saltbridge_pneg:
                if _crosses_interface(sb, peptide_chain, mhc_set):
                    salt_bridge_count += 1
                    _record_peptide_residue(residue_interactions, sb, peptide_chain, "salt_bridge")

            # pi-stacking
            for ps in interaction.pistacking:
                if _crosses_interface(ps, peptide_chain, mhc_set):
                    pi_stacking_count += 1
                    _record_peptide_residue(residue_interactions, ps, peptide_chain, "pi_stacking")

            # Cation-pi
            for pc in interaction.pication_laro + interaction.pication_paro:
                if _crosses_interface(pc, peptide_chain, mhc_set):
                    pi_cation_count += 1
                    _record_peptide_residue(residue_interactions, pc, peptide_chain, "pi_cation")

        fp = InteractionFingerprint(
            hbond_count=hbond_count,
            hydrophobic_count=hydrophobic_count,
            salt_bridge_count=salt_bridge_count,
            pi_stacking_count=pi_stacking_count,
            pi_cation_count=pi_cation_count,
            residue_interactions=residue_interactions,
        )

        log.debug(
            "PLIP %s: %d hbonds, %d hydrophobic, %d salt bridges, %d pi-stack, %d pi-cation",
            Path(pdb_path).name,
            hbond_count, hydrophobic_count, salt_bridge_count,
            pi_stacking_count, pi_cation_count,
        )
        return fp

    except Exception as exc:
        log.warning("PLIP analysis failed for %s: %s", pdb_path, exc)
        return None


def compute_plip_tanimoto(
    fp_target: InteractionFingerprint,
    fp_candidate: InteractionFingerprint,
) -> float:
    """
    Generalized Tanimoto similarity between two interaction fingerprints.

    Uses count vectors [hbond, hydrophobic, salt_bridge, pi_stack, pi_cation]:
        T(A, B) = A.B / (|A|^2 + |B|^2 - A.B)

    Returns float in [0, 1]. Returns 1.0 if both vectors are zero.
    """
    vec_t = fp_target.to_vector()
    vec_c = fp_candidate.to_vector()

    dot = float(np.dot(vec_t, vec_c))
    norm_t = float(np.dot(vec_t, vec_t))
    norm_c = float(np.dot(vec_c, vec_c))
    denom = norm_t + norm_c - dot

    if denom < 1e-12:
        return 1.0 if (norm_t < 1e-12 and norm_c < 1e-12) else 0.0

    return max(0.0, min(1.0, dot / denom))


# ======================================================================
#  2. FreeSASA: Buried Surface Area (BSA)
# ======================================================================

def compute_bsa(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[float]:
    """
    Compute buried surface area (BSA) at the peptide-MHC interface.

    BSA = SASA(peptide alone) + SASA(MHC alone) - SASA(complex)

    Returns BSA in A^2, or None if computation fails.
    """
    try:
        import freesasa
    except ImportError:
        log.warning("FreeSASA not installed. Run: pip install freesasa")
        return None

    try:
        pdb_lines = Path(pdb_path).read_text().splitlines()
        atom_lines = [l for l in pdb_lines if l.startswith(("ATOM", "HETATM"))]

        peptide_lines = [l for l in atom_lines if len(l) > 21 and l[21] == peptide_chain]
        mhc_lines = [l for l in atom_lines if len(l) > 21 and l[21] in set(mhc_chains)]

        if not peptide_lines or not mhc_lines:
            log.warning("Could not find peptide chain %s or MHC chains %s in %s",
                        peptide_chain, mhc_chains, pdb_path)
            return None

        def _sasa_from_lines(lines: List[str]) -> float:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
                f.write("\n".join(lines) + "\nEND\n")
                tmp_path = f.name
            try:
                structure = freesasa.Structure(tmp_path)
                result = freesasa.calc(structure)
                return result.totalArea()
            finally:
                Path(tmp_path).unlink(missing_ok=True)

        sasa_complex = _sasa_from_lines(peptide_lines + mhc_lines)
        sasa_peptide = _sasa_from_lines(peptide_lines)
        sasa_mhc = _sasa_from_lines(mhc_lines)

        bsa = sasa_peptide + sasa_mhc - sasa_complex

        log.debug("FreeSASA %s: peptide=%.1f, mhc=%.1f, complex=%.1f, BSA=%.1f A^2",
                  Path(pdb_path).name, sasa_peptide, sasa_mhc, sasa_complex, bsa)

        return max(0.0, bsa)

    except Exception as exc:
        log.warning("FreeSASA computation failed for %s: %s", pdb_path, exc)
        return None


def compute_bsa_similarity(
    bsa_target: float,
    bsa_candidate: float,
    max_bsa_diff: float = 800.0,
) -> float:
    """
    BSA similarity: 1.0 when identical, 0.0 when difference >= max_bsa_diff.

    For 9-mer pMHC, BSA typically ranges 800-1500 A^2.
    """
    diff = abs(bsa_target - bsa_candidate)
    return max(0.0, 1.0 - diff / max_bsa_diff)


# ======================================================================
#  3. PRODIGY: Contact-Based Binding Affinity
# ======================================================================

def compute_prodigy_affinity(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chain: str = "M",
) -> Optional[float]:
    """
    Predict binding dG (kcal/mol) using PRODIGY.

    Returns predicted dG (negative = favorable binding), or None.
    """
    # Attempt 1: Python API
    try:
        from prodigy.predict import predict_IC
        result = predict_IC(pdb_path, [peptide_chain], [mhc_chain])
        if result and hasattr(result, "dg_predicted"):
            dg = float(result.dg_predicted)
            log.debug("PRODIGY (API) %s: dG = %.2f kcal/mol", Path(pdb_path).name, dg)
            return dg
    except (ImportError, AttributeError):
        pass
    except Exception as exc:
        log.debug("PRODIGY Python API failed: %s, trying CLI", exc)

    # Attempt 2: CLI fallback
    try:
        result = subprocess.run(
            ["prodigy", pdb_path, "--selection", peptide_chain, mhc_chain],
            capture_output=True, text=True, timeout=30,
        )
        for line in result.stdout.splitlines():
            if "binding affinity" in line.lower() or "predicted" in line.lower():
                match = re.search(r"(-?\d+\.?\d*)\s*kcal", line)
                if match:
                    dg = float(match.group(1))
                    log.debug("PRODIGY (CLI) %s: dG = %.2f kcal/mol", Path(pdb_path).name, dg)
                    return dg
        # Try parsing the last numeric value as dG
        for line in reversed(result.stdout.splitlines()):
            match = re.search(r"(-?\d+\.?\d+)", line)
            if match:
                dg = float(match.group(1))
                if -30.0 < dg < 0.0:
                    log.debug("PRODIGY (CLI-fallback) %s: dG = %.2f kcal/mol", Path(pdb_path).name, dg)
                    return dg
    except FileNotFoundError:
        log.warning("PRODIGY CLI not found in PATH. Run: pip install prodigy-prot")
    except subprocess.TimeoutExpired:
        log.warning("PRODIGY timed out for %s", pdb_path)
    except Exception as exc:
        log.warning("PRODIGY failed for %s: %s", pdb_path, exc)

    return None


def compute_prodigy_similarity(
    dg_target: float,
    dg_candidate: float,
    max_dg_diff: float = 5.0,
) -> float:
    """
    dG similarity: candidates with similar predicted binding affinity
    to the target are more likely to compete for TCR binding.

    For pMHC-I, dG typically ranges -5 to -15 kcal/mol.
    max_dg_diff=5.0 kcal/mol as normalization range.

    Returns float in [0, 1].
    """
    diff = abs(dg_target - dg_candidate)
    return max(0.0, 1.0 - diff / max_dg_diff)


# ======================================================================
#  Composite: All Descriptors
# ======================================================================

def compute_all_descriptors(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> InterfaceDescriptors:
    """
    Compute all interface descriptors for a single pMHC structure.
    Runs PLIP, FreeSASA, and PRODIGY in sequence.
    """
    fp = compute_plip_fingerprint(pdb_path, peptide_chain, mhc_chains)
    bsa = compute_bsa(pdb_path, peptide_chain, mhc_chains)
    dg = compute_prodigy_affinity(pdb_path, peptide_chain, mhc_chains[0])

    return InterfaceDescriptors(
        fingerprint=fp,
        bsa_total=bsa,
        prodigy_dg=dg,
    )


def compute_interface_similarity(
    target_pdb: str,
    candidate_pdb: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    weights: Optional[Dict[str, float]] = None,
) -> InterfaceSimilarity:
    """
    Compute pairwise interface similarity between target and candidate pMHC.

    Runs all three descriptors on both structures, computes per-descriptor
    similarity, and returns a weighted combination.

    Default weights: {"plip": 0.40, "bsa": 0.25, "prodigy": 0.35}
    - PLIP (0.40): non-covalent interaction pattern is the most direct TCR recognition determinant
    - PRODIGY (0.35): binding affinity reflects thermodynamic stability
    - BSA (0.25): interface area is necessary but not sufficient

    If a descriptor fails, its weight is redistributed to successful ones.
    """
    if weights is None:
        weights = {"plip": 0.40, "bsa": 0.25, "prodigy": 0.35}

    desc_target = compute_all_descriptors(target_pdb, peptide_chain, mhc_chains)
    desc_candidate = compute_all_descriptors(candidate_pdb, peptide_chain, mhc_chains)

    # PLIP similarity
    plip_sim = 0.0
    plip_ok = False
    if desc_target.fingerprint is not None and desc_candidate.fingerprint is not None:
        plip_sim = compute_plip_tanimoto(desc_target.fingerprint, desc_candidate.fingerprint)
        plip_ok = True

    # BSA similarity
    bsa_sim = 0.0
    bsa_ok = False
    if desc_target.bsa_total is not None and desc_candidate.bsa_total is not None:
        bsa_sim = compute_bsa_similarity(desc_target.bsa_total, desc_candidate.bsa_total)
        bsa_ok = True

    # PRODIGY similarity
    prodigy_sim = 0.0
    prodigy_ok = False
    if desc_target.prodigy_dg is not None and desc_candidate.prodigy_dg is not None:
        prodigy_sim = compute_prodigy_similarity(desc_target.prodigy_dg, desc_candidate.prodigy_dg)
        prodigy_ok = True

    # Adaptive weighting: redistribute failed descriptor weights
    active = {}
    if plip_ok:
        active["plip"] = weights["plip"]
    if bsa_ok:
        active["bsa"] = weights["bsa"]
    if prodigy_ok:
        active["prodigy"] = weights["prodigy"]

    if not active:
        return InterfaceSimilarity(
            plip_tanimoto=plip_sim,
            bsa_similarity=bsa_sim,
            prodigy_similarity=prodigy_sim,
            combined=0.0,
        )

    # Normalize weights to sum to 1.0
    total_weight = sum(active.values())
    norm = {k: v / total_weight for k, v in active.items()}

    combined = (
        norm.get("plip", 0.0) * plip_sim
        + norm.get("bsa", 0.0) * bsa_sim
        + norm.get("prodigy", 0.0) * prodigy_sim
    )

    log.debug(
        "Interface similarity: PLIP=%.3f%s BSA=%.3f%s PRODIGY=%.3f%s -> combined=%.3f",
        plip_sim, "" if plip_ok else "(N/A)",
        bsa_sim, "" if bsa_ok else "(N/A)",
        prodigy_sim, "" if prodigy_ok else "(N/A)",
        combined,
    )

    return InterfaceSimilarity(
        plip_tanimoto=plip_sim,
        bsa_similarity=bsa_sim,
        prodigy_similarity=prodigy_sim,
        combined=combined,
    )


# ======================================================================
#  Internal Helpers
# ======================================================================

def _crosses_interface(interaction, peptide_chain: str, mhc_chains: set) -> bool:
    """
    Check if a PLIP interaction crosses the peptide-MHC interface.

    A cross-interface interaction has one partner in the peptide chain
    and the other in an MHC chain.
    """
    try:
        chains = set()

        for attr_chain in ("reschain", "reschain_l"):
            val = getattr(interaction, attr_chain, None)
            if val:
                chains.add(val)

        for attr_chain in ("chain", "chain_l"):
            val = getattr(interaction, attr_chain, None)
            if val:
                chains.add(val)

        has_pep = peptide_chain in chains
        has_mhc = bool(chains & mhc_chains)
        return has_pep and has_mhc

    except Exception:
        return False


def _record_peptide_residue(
    residue_interactions: Dict[int, List[str]],
    interaction,
    peptide_chain: str,
    interaction_type: str,
) -> None:
    """Record which peptide residue is involved in an interaction."""
    try:
        for attr_chain, attr_resnr in [
            ("reschain", "resnr"),
            ("reschain_l", "resnr_l"),
            ("chain", "resnr"),
        ]:
            chain_val = getattr(interaction, attr_chain, None)
            resnr_val = getattr(interaction, attr_resnr, None)
            if chain_val == peptide_chain and resnr_val is not None:
                resnum = int(resnr_val)
                if resnum not in residue_interactions:
                    residue_interactions[resnum] = []
                residue_interactions[resnum].append(interaction_type)
                return
    except Exception:
        pass
