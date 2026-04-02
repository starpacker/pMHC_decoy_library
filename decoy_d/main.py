import argparse
import logging
from pathlib import Path
from .scanner import run_decoy_d

def main():
    parser = argparse.ArgumentParser(description="Decoy D Pipeline (MPNN + mhcflurry)")
    parser.add_argument("--target", required=True, help="Target sequence")
    parser.add_argument("--hla", default="HLA-A*02:01", help="HLA Allele")
    parser.add_argument("--target-pdb", help="Optional pre-computed target PDB")
    parser.add_argument("--designs", type=int, default=1000, help="Number of MPNN designs")
    parser.add_argument("--top-k", type=int, default=10, help="Top K to predict 3D structures for")
    args = parser.parse_args()
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    
    output_dir = Path(f"data/decoy_d/{args.target}")
    df = run_decoy_d(
        target_sequence=args.target,
        hla_allele=args.hla,
        target_pdb=args.target_pdb,
        num_designs=args.designs,
        top_k_structures=args.top_k,
        output_dir=output_dir,
    )
    
    print("\n================== DECOY D RESULTS ==================")
    print(df.head(args.top_k).to_string())
    print("=====================================================")
    
if __name__ == "__main__":
    main()
