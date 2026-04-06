"""Quick test to verify PLIP deployment on a pMHC structure.

PLIP is designed for protein-ligand interactions. In pMHC structures,
the peptide is encoded as ATOM (protein), so we convert it to HETATM
to make PLIP treat it as a "ligand" and detect its MHC interactions.
"""
import sys
import tempfile
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))

from plip.structure.preparation import PDBComplex


def prepare_pdb_for_plip(pdb_path: str, peptide_chain: str = "P") -> str:
    """Convert peptide chain ATOM -> HETATM so PLIP sees it as a ligand."""
    tmp = tempfile.NamedTemporaryFile(
        suffix=".pdb", delete=False, mode="w", encoding="utf-8",
    )
    with open(pdb_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ATOM") and len(line) >= 22:
                chain_id = line[21]
                if chain_id == peptide_chain:
                    tmp.write("HETATM" + line[6:])
                    continue
            tmp.write(line)
    tmp.close()
    return tmp.name


# Test on GILGFVFTL target
pdb_path = str(project_root / "data" / "GILGFVFTL_summary" / "Decoy_B"
               / "3D_structures" / "pmhc_GILGFVFTL_HLA-A0201.pdb")

print(f"Testing PLIP on: {Path(pdb_path).name}")
print("Step 1: Converting peptide chain P to HETATM...")
plip_pdb = prepare_pdb_for_plip(pdb_path)

print("Step 2: Running PLIP analysis...")
mol = PDBComplex()
mol.load_pdb(plip_pdb)
mol.analyze()

print(f"Binding sites found: {len(mol.interaction_sets)}")
for bsid, inter in mol.interaction_sets.items():
    n_hbond = len(inter.hbonds_pdon) + len(inter.hbonds_ldon)
    n_hydro = len(inter.hydrophobic_contacts)
    n_salt = len(inter.saltbridge_lneg) + len(inter.saltbridge_pneg)
    n_pi = len(inter.pistacking)
    n_pication = len(inter.pication_laro) + len(inter.pication_paro)
    print(f"  {bsid}:")
    print(f"    H-bonds:      {n_hbond}")
    print(f"    Hydrophobic:   {n_hydro}")
    print(f"    Salt bridges:  {n_salt}")
    print(f"    Pi-stacking:   {n_pi}")
    print(f"    Pi-cation:     {n_pication}")

# Also test on a candidate for Tanimoto comparison
cand_pdb = str(project_root / "data" / "GILGFVFTL_summary" / "Decoy_B"
               / "3D_structures" / "pmhc_LLVGFVFVV_HLA-A0201.pdb")
if Path(cand_pdb).exists():
    print(f"\nStep 3: Comparing with candidate {Path(cand_pdb).stem}...")
    plip_pdb2 = prepare_pdb_for_plip(cand_pdb)
    mol2 = PDBComplex()
    mol2.load_pdb(plip_pdb2)
    mol2.analyze()
    for bsid, inter in mol2.interaction_sets.items():
        n_hbond = len(inter.hbonds_pdon) + len(inter.hbonds_ldon)
        n_hydro = len(inter.hydrophobic_contacts)
        print(f"  {bsid}: H-bonds={n_hbond}, Hydrophobic={n_hydro}")
    import os
    os.unlink(plip_pdb2)

import os
os.unlink(plip_pdb)

total_sites = len(mol.interaction_sets)
if total_sites > 0:
    print(f"\nPLIP deployment test PASSED! ({total_sites} binding site(s) detected)")
else:
    print("\nWARNING: No binding sites found. PLIP may need further configuration.")
