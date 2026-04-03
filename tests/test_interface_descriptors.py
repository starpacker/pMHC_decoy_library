"""
Tests for the interface descriptor module (Phase 1).

Uses existing tFold PDB outputs in data/decoy_b/ for integration tests.
Unit tests for similarity functions use synthetic data.
"""

import pytest
import numpy as np
from pathlib import Path
import sys

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from decoy_b.tools.interface_descriptors import (
    InteractionFingerprint,
    InterfaceSimilarity,
    compute_plip_fingerprint,
    compute_plip_tanimoto,
    compute_bsa,
    compute_bsa_similarity,
    compute_prodigy_affinity,
    compute_prodigy_similarity,
    compute_interface_similarity,
    compute_all_descriptors,
)


# -- Paths --

DATA_DIR = PROJECT_ROOT / "data" / "decoy_b"
STEPWISE_DIR = DATA_DIR / "stepwise_demo" / "tfold_pdbs"


def _find_pdbs() -> list:
    """Find available PDB files for testing."""
    pdbs = []
    for d in [STEPWISE_DIR, DATA_DIR]:
        if d.exists():
            pdbs.extend(sorted(d.rglob("*.pdb")))
    return pdbs


# -- Unit Tests (no PDB needed) --

class TestTanimoto:
    """Test Tanimoto similarity computation on synthetic fingerprints."""

    def test_identical_fingerprints(self):
        fp = InteractionFingerprint(
            hbond_count=5, hydrophobic_count=10, salt_bridge_count=2,
        )
        assert compute_plip_tanimoto(fp, fp) == pytest.approx(1.0)

    def test_zero_fingerprints(self):
        fp1 = InteractionFingerprint()  # all zeros
        fp2 = InteractionFingerprint()
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(1.0)

    def test_disjoint_fingerprints(self):
        fp1 = InteractionFingerprint(hbond_count=10)
        fp2 = InteractionFingerprint(hydrophobic_count=10)
        # Tanimoto of orthogonal vectors = 0
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(0.0)

    def test_similar_fingerprints(self):
        fp1 = InteractionFingerprint(hbond_count=5, hydrophobic_count=10)
        fp2 = InteractionFingerprint(hbond_count=4, hydrophobic_count=9)
        sim = compute_plip_tanimoto(fp1, fp2)
        assert 0.5 < sim < 1.0

    def test_symmetry(self):
        fp1 = InteractionFingerprint(hbond_count=3, salt_bridge_count=1)
        fp2 = InteractionFingerprint(hbond_count=5, hydrophobic_count=2)
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(
            compute_plip_tanimoto(fp2, fp1)
        )

    def test_to_vector_shape(self):
        fp = InteractionFingerprint(hbond_count=3, pi_cation_count=1)
        vec = fp.to_vector()
        assert vec.shape == (5,)
        assert vec[0] == 3.0  # hbond
        assert vec[4] == 1.0  # pi_cation

    def test_to_dict(self):
        fp = InteractionFingerprint(hbond_count=2, hydrophobic_count=7)
        d = fp.to_dict()
        assert d["hbond_count"] == 2
        assert d["hydrophobic_count"] == 7
        assert d["salt_bridge_count"] == 0


class TestBSASimilarity:
    """Test BSA similarity computation."""

    def test_identical(self):
        assert compute_bsa_similarity(1000.0, 1000.0) == pytest.approx(1.0)

    def test_max_diff(self):
        assert compute_bsa_similarity(1000.0, 1800.0) == pytest.approx(0.0)

    def test_mid_range(self):
        sim = compute_bsa_similarity(1000.0, 1400.0)
        assert sim == pytest.approx(0.5)

    def test_symmetry(self):
        assert compute_bsa_similarity(800.0, 1200.0) == pytest.approx(
            compute_bsa_similarity(1200.0, 800.0)
        )

    def test_custom_max_diff(self):
        sim = compute_bsa_similarity(1000.0, 1500.0, max_bsa_diff=1000.0)
        assert sim == pytest.approx(0.5)

    def test_beyond_max_clamps_to_zero(self):
        sim = compute_bsa_similarity(100.0, 2000.0)
        assert sim == pytest.approx(0.0)


class TestProdigySimilarity:
    """Test PRODIGY dG similarity computation."""

    def test_identical(self):
        assert compute_prodigy_similarity(-10.0, -10.0) == pytest.approx(1.0)

    def test_max_diff(self):
        assert compute_prodigy_similarity(-10.0, -15.0) == pytest.approx(0.0)

    def test_mid_range(self):
        sim = compute_prodigy_similarity(-10.0, -12.5)
        assert sim == pytest.approx(0.5)

    def test_symmetry(self):
        assert compute_prodigy_similarity(-8.0, -12.0) == pytest.approx(
            compute_prodigy_similarity(-12.0, -8.0)
        )

    def test_custom_max_diff(self):
        sim = compute_prodigy_similarity(-10.0, -12.0, max_dg_diff=4.0)
        assert sim == pytest.approx(0.5)


class TestInterfaceSimilarity:
    """Test InterfaceSimilarity data class."""

    def test_to_dict(self):
        sim = InterfaceSimilarity(
            plip_tanimoto=0.8, bsa_similarity=0.6,
            prodigy_similarity=0.7, combined=0.72,
        )
        d = sim.to_dict()
        assert d["plip_tanimoto"] == 0.8
        assert d["combined"] == 0.72


class TestModelBackwardCompatibility:
    """Test that StructuralScore and DecoyABEntry remain backward-compatible."""

    def test_structural_score_old_style(self):
        from decoy_a.models import StructuralScore
        s = StructuralScore(modeling_tool="test", surface_correlation=0.5)
        assert s.plip_tanimoto is None
        assert s.interface_combined is None

    def test_structural_score_new_style(self):
        from decoy_a.models import StructuralScore
        s = StructuralScore(
            modeling_tool="biopython+dual_superposition+interface",
            surface_correlation=0.8,
            rmsd=1.2,
            plip_tanimoto=0.75,
            bsa_similarity=0.6,
            prodigy_similarity=0.8,
            interface_combined=0.72,
        )
        assert s.plip_tanimoto == 0.75
        d = s.model_dump(mode="json")
        assert "plip_tanimoto" in d
        assert "interface_combined" in d

    def test_decoy_ab_entry_new_fields(self):
        from decoy_a.models import DecoyABEntry
        entry = DecoyABEntry(
            decoy_ab_id="DAB-0001",
            sequence="GILGFVFTL",
            target_sequence="GILGFVFTL",
            hla_allele="HLA-A*02:01",
            source="Decoy_B_Structural_Similarity",
            el_rank=0.5,
            hamming_distance=0,
            structural_similarity=0.85,
            plip_tanimoto=0.75,
            bsa_similarity=0.6,
            prodigy_similarity=0.8,
            interface_combined=0.72,
            total_risk_score=1.5,
        )
        d = entry.model_dump(mode="json")
        assert d["plip_tanimoto"] == 0.75
        assert d["interface_combined"] == 0.72


# -- Integration Tests (need PDB files) --

class TestPLIPIntegration:
    """Test PLIP fingerprint extraction on real tFold PDBs."""

    def test_fingerprint_extraction(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        fp = compute_plip_fingerprint(str(pdbs[0]))
        if fp is None:
            pytest.skip("PLIP not available (missing OpenBabel?)")

        assert fp.hbond_count >= 0
        assert fp.hydrophobic_count >= 0
        total = fp.hbond_count + fp.hydrophobic_count + fp.salt_bridge_count
        assert total > 0, "Expected at least some interactions in pMHC"

    def test_self_tanimoto_is_one(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        fp = compute_plip_fingerprint(str(pdbs[0]))
        if fp is None:
            pytest.skip("PLIP not available")

        assert compute_plip_tanimoto(fp, fp) == pytest.approx(1.0)


class TestFreeSASAIntegration:
    """Test FreeSASA BSA computation on real tFold PDBs."""

    def test_bsa_positive(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        bsa = compute_bsa(str(pdbs[0]))
        if bsa is None:
            pytest.skip("FreeSASA not available")

        assert bsa > 0, f"Expected positive BSA, got {bsa}"
        assert bsa < 3000, f"BSA unreasonably large: {bsa}"


class TestPRODIGYIntegration:
    """Test PRODIGY affinity prediction on real tFold PDBs."""

    def test_prodigy_returns_negative_dg(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        dg = compute_prodigy_affinity(str(pdbs[0]))
        if dg is None:
            pytest.skip("PRODIGY not available")

        assert dg < 0, f"Expected negative dG, got {dg}"


class TestEndToEnd:
    """End-to-end test: full interface similarity between two pMHCs."""

    def test_interface_similarity(self):
        pdbs = _find_pdbs()
        if len(pdbs) < 2:
            pytest.skip("Need at least 2 PDB files")

        sim = compute_interface_similarity(str(pdbs[0]), str(pdbs[1]))

        assert 0.0 <= sim.combined <= 1.0
        assert 0.0 <= sim.plip_tanimoto <= 1.0
        assert 0.0 <= sim.bsa_similarity <= 1.0
        assert 0.0 <= sim.prodigy_similarity <= 1.0

    def test_self_similarity_is_one(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        sim = compute_interface_similarity(str(pdbs[0]), str(pdbs[0]))

        if sim.plip_tanimoto is not None:
            assert sim.plip_tanimoto == pytest.approx(1.0)
        if sim.bsa_similarity is not None:
            assert sim.bsa_similarity == pytest.approx(1.0)
        if sim.prodigy_similarity is not None:
            assert sim.prodigy_similarity == pytest.approx(1.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
