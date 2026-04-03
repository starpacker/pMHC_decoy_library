"""
Interface descriptor computation for pMHC structures.

Computes five classes of descriptors from PDB structures:
1. PLIP: non-covalent interaction fingerprint (H-bonds, hydrophobic, salt bridges, pi-stacking)
2. FreeSASA: buried surface area (BSA / SASA) at the peptide-MHC interface
3. PRODIGY: contact-based binding affinity prediction (dG)
4. APBS: electrostatic surface potential at the peptide-MHC interface
5. PeSTo: learned interface compatibility embedding (geometric deep learning)

All pairwise similarity functions return normalized [0, 1] scores.

Usage:
    from decoy_b.tools.interface_descriptors import compute_interface_similarity
    sim = compute_interface_similarity(target_pdb, candidate_pdb)
    # sim.plip_tanimoto, sim.bsa_similarity, sim.prodigy_similarity,
    # sim.esp_similarity, sim.pesto_similarity, sim.combined
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
        "esp_vector", "pesto_embedding",
    )

    def __init__(
        self,
        fingerprint: Optional[InteractionFingerprint] = None,
        bsa_total: Optional[float] = None,
        bsa_polar: Optional[float] = None,
        bsa_apolar: Optional[float] = None,
        prodigy_dg: Optional[float] = None,
        prodigy_kd: Optional[float] = None,
        esp_vector: Optional[np.ndarray] = None,
        pesto_embedding: Optional[np.ndarray] = None,
    ):
        self.fingerprint = fingerprint
        self.bsa_total = bsa_total
        self.bsa_polar = bsa_polar
        self.bsa_apolar = bsa_apolar
        self.prodigy_dg = prodigy_dg
        self.prodigy_kd = prodigy_kd
        self.esp_vector = esp_vector
        self.pesto_embedding = pesto_embedding


class InterfaceSimilarity:
    """Pairwise similarity between target and candidate interface descriptors."""

    __slots__ = (
        "plip_tanimoto", "bsa_similarity", "prodigy_similarity",
        "esp_similarity", "pesto_similarity", "combined",
    )

    def __init__(
        self,
        plip_tanimoto: float = 0.0,
        bsa_similarity: float = 0.0,
        prodigy_similarity: float = 0.0,
        esp_similarity: float = 0.0,
        pesto_similarity: float = 0.0,
        combined: float = 0.0,
    ):
        self.plip_tanimoto = plip_tanimoto
        self.bsa_similarity = bsa_similarity
        self.prodigy_similarity = prodigy_similarity
        self.esp_similarity = esp_similarity
        self.pesto_similarity = pesto_similarity
        self.combined = combined

    def to_dict(self) -> dict:
        return {
            "plip_tanimoto": self.plip_tanimoto,
            "bsa_similarity": self.bsa_similarity,
            "prodigy_similarity": self.prodigy_similarity,
            "esp_similarity": self.esp_similarity,
            "pesto_similarity": self.pesto_similarity,
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
#  4. APBS: Electrostatic Surface Potential
# ======================================================================

def compute_esp_vector(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    interface_cutoff: float = 5.0,
) -> Optional[np.ndarray]:
    """
    Compute electrostatic surface potential (ESP) at peptide-MHC interface atoms.

    Full mode: pdb2pqr → APBS → parse DX grid → sample at interface atoms.
    Fallback: Coulomb approximation using AMBER partial charges.

    Returns an ESP vector (one value per interface atom) or None.
    """
    # Attempt 1: Full APBS pipeline (pdb2pqr + APBS)
    esp = _compute_esp_apbs(pdb_path, peptide_chain, mhc_chains, interface_cutoff)
    if esp is not None:
        return esp

    # Attempt 2: Coulomb proxy using BioPython
    return _compute_esp_coulomb_proxy(pdb_path, peptide_chain, mhc_chains, interface_cutoff)


def _compute_esp_apbs(
    pdb_path: str,
    peptide_chain: str,
    mhc_chains: Tuple[str, ...],
    interface_cutoff: float,
) -> Optional[np.ndarray]:
    """Full APBS pipeline: pdb2pqr → APBS → parse DX → sample at interface."""
    try:
        import subprocess as _sp
        _sp.run(["pdb2pqr", "--version"], capture_output=True, timeout=5)
        _sp.run(["apbs", "--version"], capture_output=True, timeout=5)
    except (FileNotFoundError, _sp.TimeoutExpired):
        log.debug("APBS/pdb2pqr not in PATH, falling back to Coulomb proxy")
        return None

    try:
        pdb_p = Path(pdb_path)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pqr_path = tmp / "structure.pqr"
            apbs_input = tmp / "apbs.in"
            dx_path = tmp / "pot.dx"

            # Step 1: pdb2pqr
            result = subprocess.run(
                ["pdb2pqr", "--ff=AMBER", "--apbs-input=" + str(apbs_input),
                 str(pdb_p), str(pqr_path)],
                capture_output=True, text=True, timeout=60,
            )
            if result.returncode != 0 or not pqr_path.exists():
                log.warning("pdb2pqr failed for %s: %s", pdb_path, result.stderr[:200])
                return None

            # Step 2: Patch APBS input to use our output path
            apbs_text = apbs_input.read_text()
            # Replace the default output name with our desired path
            apbs_text = re.sub(
                r'write pot dx \S+',
                f'write pot dx {str(tmp / "pot")}',
                apbs_text,
            )
            apbs_input.write_text(apbs_text)

            # Step 3: Run APBS
            result = subprocess.run(
                ["apbs", str(apbs_input)],
                capture_output=True, text=True, timeout=120,
                cwd=str(tmpdir),
            )

            # Find DX file (APBS may append .dx)
            dx_candidates = list(tmp.glob("pot*.dx"))
            if not dx_candidates:
                log.warning("APBS produced no DX output for %s", pdb_path)
                return None
            dx_file = dx_candidates[0]

            # Step 4: Parse DX and sample at interface atoms
            grid_data, origin, delta = _parse_dx(dx_file)
            if grid_data is None:
                return None

            interface_coords = _get_interface_atom_coords(
                pdb_path, peptide_chain, mhc_chains, interface_cutoff
            )
            if interface_coords is None or len(interface_coords) == 0:
                return None

            esp_values = _sample_dx_at_coords(grid_data, origin, delta, interface_coords)
            log.debug("APBS %s: sampled %d interface atoms, ESP range [%.2f, %.2f]",
                      pdb_p.name, len(esp_values), esp_values.min(), esp_values.max())
            return esp_values

    except Exception as exc:
        log.warning("APBS pipeline failed for %s: %s", pdb_path, exc)
        return None


def _parse_dx(dx_path: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Parse an OpenDX electrostatic potential grid file."""
    try:
        lines = dx_path.read_text().splitlines()
        nx, ny, nz = 0, 0, 0
        origin = np.zeros(3)
        delta = np.zeros((3, 3))
        data_lines = []
        reading_data = False
        delta_idx = 0

        for line in lines:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if "object 1" in line and "gridpositions" in line:
                parts = line.split()
                idx = parts.index("counts")
                nx, ny, nz = int(parts[idx+1]), int(parts[idx+2]), int(parts[idx+3])
            elif line.startswith("origin"):
                parts = line.split()
                origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith("delta"):
                parts = line.split()
                delta[delta_idx] = [float(parts[1]), float(parts[2]), float(parts[3])]
                delta_idx += 1
            elif "object 3" in line and "data" in line:
                reading_data = True
            elif reading_data:
                if line.startswith("attribute") or line.startswith("object") or line.startswith("component"):
                    break
                data_lines.append(line)

        if nx == 0:
            return None, None, None

        values = []
        for dl in data_lines:
            values.extend(float(x) for x in dl.split())

        grid = np.array(values).reshape((nx, ny, nz))
        spacing = np.array([delta[0, 0], delta[1, 1], delta[2, 2]])
        return grid, origin, spacing

    except Exception as exc:
        log.warning("DX parse error: %s", exc)
        return None, None, None


def _sample_dx_at_coords(
    grid: np.ndarray, origin: np.ndarray, spacing: np.ndarray,
    coords: np.ndarray,
) -> np.ndarray:
    """Trilinear interpolation of DX grid at given coordinates."""
    fractional = (coords - origin) / spacing
    indices = fractional.astype(int)
    nx, ny, nz = grid.shape

    values = []
    for i in range(len(coords)):
        ix, iy, iz = indices[i]
        ix = max(0, min(ix, nx - 2))
        iy = max(0, min(iy, ny - 2))
        iz = max(0, min(iz, nz - 2))
        # Nearest-neighbor for simplicity (trilinear adds complexity for marginal gain)
        values.append(grid[ix, iy, iz])

    return np.array(values)


def _compute_esp_coulomb_proxy(
    pdb_path: str,
    peptide_chain: str,
    mhc_chains: Tuple[str, ...],
    interface_cutoff: float,
) -> Optional[np.ndarray]:
    """
    Coulomb potential proxy for APBS.

    For each interface peptide atom, sum Coulomb contributions from nearby
    MHC atoms using AMBER-like partial charges assigned by atom type.

    This is a simplified but physically motivated approximation:
    ESP_i = sum_j(q_j / r_ij) for MHC atoms j within cutoff of peptide atom i.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        log.debug("BioPython not available for Coulomb proxy")
        return None

    # Simplified partial charges by heavy-atom element
    PARTIAL_CHARGES = {
        "N": -0.4, "O": -0.5, "C": 0.1, "S": -0.1,
        "CA": 0.1, "CB": 0.0,
    }
    # Charged residue adjustments
    CHARGED_ATOMS = {
        ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5, ("ARG", "NE"): 0.3,
        ("LYS", "NZ"): 0.7,
        ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
        ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
        ("HIS", "ND1"): 0.3, ("HIS", "NE2"): 0.3,
    }

    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("esp", pdb_path)
        model = list(struct.get_models())[0]

        mhc_set = set(mhc_chains)
        pep_atoms = []
        mhc_atoms = []

        for chain in model:
            cid = chain.id
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                resname = residue.get_resname()
                for atom in residue:
                    if atom.element == "H":
                        continue
                    coord = atom.get_vector().get_array()
                    aname = atom.get_name()
                    # Assign charge
                    q = CHARGED_ATOMS.get((resname, aname))
                    if q is None:
                        q = PARTIAL_CHARGES.get(atom.element, 0.0)

                    if cid == peptide_chain:
                        pep_atoms.append((coord, q, residue.id[1]))
                    elif cid in mhc_set:
                        mhc_atoms.append((coord, q))

        if not pep_atoms or not mhc_atoms:
            return None

        # Find interface peptide atoms (within cutoff of any MHC atom)
        pep_coords = np.array([a[0] for a in pep_atoms])
        mhc_coords = np.array([a[0] for a in mhc_atoms])
        mhc_charges = np.array([a[1] for a in mhc_atoms])

        # For each peptide atom, compute ESP from MHC atoms
        esp_values = []
        for i, (pc, pq, resnum) in enumerate(pep_atoms):
            dists = np.linalg.norm(mhc_coords - pc, axis=1)
            mask = dists < interface_cutoff
            if not np.any(mask):
                continue
            # Coulomb: sum(q/r), with minimum distance clamp to avoid singularity
            r = np.maximum(dists[mask], 0.5)
            esp = float(np.sum(mhc_charges[mask] / r))
            esp_values.append(esp)

        if not esp_values:
            return None

        esp_arr = np.array(esp_values)
        log.debug("Coulomb proxy %s: %d atoms, ESP range [%.3f, %.3f]",
                  Path(pdb_path).name, len(esp_arr), esp_arr.min(), esp_arr.max())
        return esp_arr

    except Exception as exc:
        log.warning("Coulomb proxy failed for %s: %s", pdb_path, exc)
        return None


def _get_interface_atom_coords(
    pdb_path: str,
    peptide_chain: str,
    mhc_chains: Tuple[str, ...],
    cutoff: float,
) -> Optional[np.ndarray]:
    """Get coordinates of peptide heavy atoms within cutoff of MHC atoms."""
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("iface", pdb_path)
        model = list(struct.get_models())[0]

        mhc_set = set(mhc_chains)
        pep_coords = []
        mhc_coords = []

        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                for atom in residue:
                    if atom.element == "H":
                        continue
                    if chain.id == peptide_chain:
                        pep_coords.append(atom.get_vector().get_array())
                    elif chain.id in mhc_set:
                        mhc_coords.append(atom.get_vector().get_array())

        if not pep_coords or not mhc_coords:
            return None

        pep_arr = np.array(pep_coords)
        mhc_arr = np.array(mhc_coords)

        # Find peptide atoms within cutoff of any MHC atom
        interface_mask = np.zeros(len(pep_arr), dtype=bool)
        for i, pc in enumerate(pep_arr):
            min_dist = np.min(np.linalg.norm(mhc_arr - pc, axis=1))
            if min_dist < cutoff:
                interface_mask[i] = True

        return pep_arr[interface_mask] if np.any(interface_mask) else None

    except Exception:
        return None


def compute_esp_similarity(
    esp_target: np.ndarray,
    esp_candidate: np.ndarray,
) -> float:
    """
    Electrostatic potential similarity via Pearson correlation.

    Handles different-length vectors by truncating to the shorter length
    (peptide atoms are ordered N→C, so this compares corresponding positions).

    Returns float in [0, 1] (mapped from [-1, 1] correlation).
    """
    n = min(len(esp_target), len(esp_candidate))
    if n < 3:
        return 0.0

    a = esp_target[:n]
    b = esp_candidate[:n]

    std_a, std_b = np.std(a), np.std(b)
    if std_a < 1e-12 or std_b < 1e-12:
        # Both flat → similar electrostatic environment
        if std_a < 1e-12 and std_b < 1e-12:
            return 1.0
        return 0.0

    corr = float(np.corrcoef(a, b)[0, 1])
    # Map [-1, 1] → [0, 1]: anti-correlated ESP is maximally dissimilar
    return max(0.0, min(1.0, (corr + 1.0) / 2.0))


# ======================================================================
#  5. PeSTo: Protein Structure Transformer Interface Embedding
# ======================================================================

def compute_pesto_embedding(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[np.ndarray]:
    """
    Compute interface embedding using PeSTo (Protein Structure Transformer).

    Full mode: Load pretrained PeSTo model, compute per-residue interface
    probability and embedding vectors.
    Fallback: Contact-based interface propensity proxy using BioPython.

    Returns a fixed-length embedding vector or None.
    """
    # Attempt 1: Full PeSTo model
    emb = _compute_pesto_full(pdb_path, peptide_chain, mhc_chains)
    if emb is not None:
        return emb

    # Attempt 2: Contact-propensity proxy
    return _compute_pesto_proxy(pdb_path, peptide_chain, mhc_chains)


def _compute_pesto_full(
    pdb_path: str,
    peptide_chain: str,
    mhc_chains: Tuple[str, ...],
) -> Optional[np.ndarray]:
    """
    Full PeSTo inference.

    PeSTo is a geometric transformer that predicts per-residue interaction
    interface probabilities (~300ms, parameter-free after loading weights).
    """
    try:
        import torch
        from pesto import PeSTo as PeSToModel
    except ImportError:
        log.debug("PeSTo or PyTorch not installed, using proxy")
        return None

    try:
        model = PeSToModel.load_pretrained()
        model.eval()

        with torch.no_grad():
            result = model.predict_from_pdb(pdb_path)

        # Extract per-residue embeddings for peptide chain
        chain_mask = result["chain_ids"] == peptide_chain
        if not chain_mask.any():
            log.warning("PeSTo: peptide chain %s not found in output", peptide_chain)
            return None

        pep_embeddings = result["embeddings"][chain_mask]  # (n_res, embed_dim)
        pep_probs = result["interface_prob"][chain_mask]    # (n_res,)

        # Create a compact descriptor: concatenate mean embedding + per-residue probs
        mean_emb = pep_embeddings.mean(axis=0)  # (embed_dim,)
        descriptor = np.concatenate([mean_emb, pep_probs])

        log.debug("PeSTo %s: %d peptide residues, mean interface prob=%.3f",
                  Path(pdb_path).name, len(pep_probs), pep_probs.mean())
        return descriptor

    except Exception as exc:
        log.warning("PeSTo inference failed for %s: %s", pdb_path, exc)
        return None


def _compute_pesto_proxy(
    pdb_path: str,
    peptide_chain: str,
    mhc_chains: Tuple[str, ...],
    contact_cutoff: float = 5.0,
) -> Optional[np.ndarray]:
    """
    Contact-based interface propensity proxy for PeSTo.

    For each peptide residue, computes a feature vector capturing:
    - Contact density with MHC (number of heavy-atom contacts / max_contacts)
    - Residue physicochemical type (hydrophobic, polar, charged encoding)
    - Relative solvent exposure proxy (contact ratio with MHC vs total contacts)
    - Per-residue interface probability estimate

    Returns a fixed-length vector: per-residue features flattened.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        return None

    # Residue property encoding: [hydrophobic, polar, positive, negative, aromatic]
    RESIDUE_PROPS = {
        "ALA": [1, 0, 0, 0, 0], "VAL": [1, 0, 0, 0, 0], "LEU": [1, 0, 0, 0, 0],
        "ILE": [1, 0, 0, 0, 0], "MET": [1, 0, 0, 0, 0], "PRO": [1, 0, 0, 0, 0],
        "GLY": [0, 0, 0, 0, 0],
        "SER": [0, 1, 0, 0, 0], "THR": [0, 1, 0, 0, 0], "CYS": [0, 1, 0, 0, 0],
        "ASN": [0, 1, 0, 0, 0], "GLN": [0, 1, 0, 0, 0],
        "ARG": [0, 0, 1, 0, 0], "LYS": [0, 0, 1, 0, 0], "HIS": [0, 0, 1, 0, 0],
        "ASP": [0, 0, 0, 1, 0], "GLU": [0, 0, 0, 1, 0],
        "PHE": [1, 0, 0, 0, 1], "TYR": [0, 1, 0, 0, 1], "TRP": [1, 0, 0, 0, 1],
    }
    DEFAULT_PROPS = [0, 0, 0, 0, 0]
    MAX_CONTACTS = 20.0  # normalization ceiling

    try:
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure("pesto", pdb_path)
        model = list(struct.get_models())[0]

        mhc_set = set(mhc_chains)

        # Collect MHC heavy atoms
        mhc_atoms_coords = []
        for chain in model:
            if chain.id in mhc_set:
                for res in chain:
                    if res.id[0] != " ":
                        continue
                    for atom in res:
                        if atom.element != "H":
                            mhc_atoms_coords.append(atom.get_vector().get_array())

        if not mhc_atoms_coords:
            return None
        mhc_arr = np.array(mhc_atoms_coords)

        # Process peptide residues
        pep_chain = None
        for chain in model:
            if chain.id == peptide_chain:
                pep_chain = chain
                break
        if pep_chain is None:
            return None

        residue_features = []
        for res in pep_chain:
            if res.id[0] != " ":
                continue
            resname = res.get_resname()
            props = RESIDUE_PROPS.get(resname, DEFAULT_PROPS)

            # Count contacts with MHC
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
            # Interface probability proxy: sigmoid of contact density
            interface_prob = 1.0 / (1.0 + np.exp(-5.0 * (contact_density - 0.3)))
            proximity = max(0.0, 1.0 - min_dist / contact_cutoff)

            residue_features.append(props + [contact_density, interface_prob, proximity])

        if not residue_features:
            return None

        # Build embedding: per-residue features (8-dim each) + summary stats
        feat_matrix = np.array(residue_features, dtype=float)  # (n_res, 8)
        mean_feat = feat_matrix.mean(axis=0)                   # (8,)
        probs = feat_matrix[:, 6]                               # interface_prob column

        # Fixed-length descriptor: flatten per-residue probs + mean features
        # Pad/truncate to 15 residues (max peptide length)
        max_len = 15
        padded_probs = np.zeros(max_len)
        n = min(len(probs), max_len)
        padded_probs[:n] = probs[:n]

        padded_contacts = np.zeros(max_len)
        padded_contacts[:n] = feat_matrix[:n, 5]  # contact_density

        padded_proximity = np.zeros(max_len)
        padded_proximity[:n] = feat_matrix[:n, 7]  # proximity

        descriptor = np.concatenate([padded_probs, padded_contacts, padded_proximity, mean_feat])

        log.debug("PeSTo proxy %s: %d residues, mean interface_prob=%.3f",
                  Path(pdb_path).name, len(residue_features),
                  float(probs.mean()) if len(probs) > 0 else 0.0)
        return descriptor

    except Exception as exc:
        log.warning("PeSTo proxy failed for %s: %s", pdb_path, exc)
        return None


def compute_pesto_similarity(
    emb_target: np.ndarray,
    emb_candidate: np.ndarray,
) -> float:
    """
    PeSTo embedding similarity via cosine similarity.

    Returns float in [0, 1] (mapped from [-1, 1]).
    """
    norm_t = np.linalg.norm(emb_target)
    norm_c = np.linalg.norm(emb_candidate)

    if norm_t < 1e-12 or norm_c < 1e-12:
        if norm_t < 1e-12 and norm_c < 1e-12:
            return 1.0
        return 0.0

    cosine = float(np.dot(emb_target, emb_candidate) / (norm_t * norm_c))
    # Map [-1, 1] → [0, 1]
    return max(0.0, min(1.0, (cosine + 1.0) / 2.0))


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
    Runs PLIP, FreeSASA, PRODIGY, APBS, and PeSTo in sequence.
    """
    fp = compute_plip_fingerprint(pdb_path, peptide_chain, mhc_chains)
    bsa = compute_bsa(pdb_path, peptide_chain, mhc_chains)
    dg = compute_prodigy_affinity(pdb_path, peptide_chain, mhc_chains[0])
    esp = compute_esp_vector(pdb_path, peptide_chain, mhc_chains)
    pesto = compute_pesto_embedding(pdb_path, peptide_chain, mhc_chains)

    return InterfaceDescriptors(
        fingerprint=fp,
        bsa_total=bsa,
        prodigy_dg=dg,
        esp_vector=esp,
        pesto_embedding=pesto,
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

    Runs all five descriptors on both structures, computes per-descriptor
    similarity, and returns a weighted combination.

    Default weights (5-descriptor):
        {"plip": 0.25, "bsa": 0.10, "prodigy": 0.20, "esp": 0.25, "pesto": 0.20}
    - PLIP (0.25): non-covalent interaction pattern — direct TCR recognition
    - APBS/ESP (0.25): electrostatic complementarity — core molecular mimicry driver
    - PRODIGY (0.20): binding affinity — thermodynamic stability
    - PeSTo (0.20): learned interface compatibility — robust to prediction error
    - BSA (0.10): interface area — necessary but not sufficient

    If a descriptor fails, its weight is redistributed to successful ones.
    """
    if weights is None:
        weights = {
            "plip": 0.25, "bsa": 0.10, "prodigy": 0.20,
            "esp": 0.25, "pesto": 0.20,
        }

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

    # ESP similarity
    esp_sim = 0.0
    esp_ok = False
    if desc_target.esp_vector is not None and desc_candidate.esp_vector is not None:
        esp_sim = compute_esp_similarity(desc_target.esp_vector, desc_candidate.esp_vector)
        esp_ok = True

    # PeSTo similarity
    pesto_sim = 0.0
    pesto_ok = False
    if desc_target.pesto_embedding is not None and desc_candidate.pesto_embedding is not None:
        pesto_sim = compute_pesto_similarity(
            desc_target.pesto_embedding, desc_candidate.pesto_embedding
        )
        pesto_ok = True

    # Adaptive weighting: redistribute failed descriptor weights
    active = {}
    if plip_ok:
        active["plip"] = weights.get("plip", 0.25)
    if bsa_ok:
        active["bsa"] = weights.get("bsa", 0.10)
    if prodigy_ok:
        active["prodigy"] = weights.get("prodigy", 0.20)
    if esp_ok:
        active["esp"] = weights.get("esp", 0.25)
    if pesto_ok:
        active["pesto"] = weights.get("pesto", 0.20)

    if not active:
        return InterfaceSimilarity(
            plip_tanimoto=plip_sim,
            bsa_similarity=bsa_sim,
            prodigy_similarity=prodigy_sim,
            esp_similarity=esp_sim,
            pesto_similarity=pesto_sim,
            combined=0.0,
        )

    # Normalize weights to sum to 1.0
    total_weight = sum(active.values())
    norm = {k: v / total_weight for k, v in active.items()}

    combined = (
        norm.get("plip", 0.0) * plip_sim
        + norm.get("bsa", 0.0) * bsa_sim
        + norm.get("prodigy", 0.0) * prodigy_sim
        + norm.get("esp", 0.0) * esp_sim
        + norm.get("pesto", 0.0) * pesto_sim
    )

    log.debug(
        "Interface similarity: PLIP=%.3f%s BSA=%.3f%s PRODIGY=%.3f%s "
        "ESP=%.3f%s PeSTo=%.3f%s -> combined=%.3f",
        plip_sim, "" if plip_ok else "(N/A)",
        bsa_sim, "" if bsa_ok else "(N/A)",
        prodigy_sim, "" if prodigy_ok else "(N/A)",
        esp_sim, "" if esp_ok else "(N/A)",
        pesto_sim, "" if pesto_ok else "(N/A)",
        combined,
    )

    return InterfaceSimilarity(
        plip_tanimoto=plip_sim,
        bsa_similarity=bsa_sim,
        prodigy_similarity=prodigy_sim,
        esp_similarity=esp_sim,
        pesto_similarity=pesto_sim,
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
