#!/usr/bin/env python3
"""
Build a small demo HLA-filtered dataset for testing Decoy B pipeline.

This creates a synthetic but realistic HLA-filtered k-mer parquet file
containing ~10,000 9-mer peptides with EL rank scores, suitable for
testing the full Decoy B pipeline without requiring the full 41M k-mer
database or hours of mhcflurry filtering.

Usage:
    python scripts/build_demo_data.py
"""

import os
import sys
import random
import string

# Force mhcflurry data to /share/liuyutian
os.environ.setdefault("MHCFLURRY_DATA_DIR", "/share/liuyutian/mhcflurry_data/4")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from pathlib import Path

from decoy_a.config import AB_DATA_DIR, DEFAULT_HLA_ALLELE

AA = list("ACDEFGHIKLMNPQRSTVWY")

# Known immunogenic peptides and their variants for realistic testing
SEED_PEPTIDES = [
    # Target: EVDPIGHLY (MAGE-A3, HLA-A*01:01)
    "EVDPIGHLY",
    # Flu M1 peptide (HLA-A*02:01)
    "GILGFVFTL",
    # Some known cross-reactive peptides
    "GLLGFVGTL",  # TAP2
    "GILLFLFTL",  # PKD1L1
    "EVDPIGHVY",  # MAGEA6
    "EVGPIFHLY",  # FGD5
    # More realistic peptides
    "YLQPRTFLL",  # SARS-CoV-2 Spike
    "KLPDDFTGCV",  # HBV
    "RMFPNAPYL",  # WT1
    "SLYNTVATL",  # HIV Gag
]


def random_peptide(length=9):
    """Generate a random peptide of given length."""
    return "".join(random.choices(AA, k=length))


def mutate_peptide(peptide, n_mutations=1):
    """Create a mutant of a peptide with n random mutations."""
    pep = list(peptide)
    positions = random.sample(range(len(pep)), min(n_mutations, len(pep)))
    for pos in positions:
        new_aa = random.choice([a for a in AA if a != pep[pos]])
        pep[pos] = new_aa
    return "".join(pep)


def build_demo_hla_filtered(
    n_peptides: int = 10000,
    hla_allele: str = DEFAULT_HLA_ALLELE,
) -> Path:
    """Build a demo HLA-filtered parquet file."""
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)

    peptides = []
    el_ranks = []
    gene_symbols_list = []
    source_proteins_list = []

    # Add seed peptides and their variants
    fake_genes = [
        "TAP2", "PKD1L1", "MAGEA6", "FGD5", "MAGEA12", "MAGEA8",
        "EPS8L2", "TITIN", "ABCB1", "BRCA1", "TP53", "EGFR",
        "MYC", "KRAS", "PTEN", "RB1", "APC", "VHL", "WT1", "NF1",
        "CDH1", "SMAD4", "MLH1", "MSH2", "BRCA2", "ATM", "CHEK2",
        "PALB2", "RAD51", "FANCA", "FANCC", "FANCG", "NBN", "MRE11",
    ]
    fake_proteins = [f"P{i:05d}" for i in range(len(fake_genes))]

    # Generate variants of seed peptides (Hamming 1-5)
    for seed in SEED_PEPTIDES:
        for hd in range(1, 6):
            for _ in range(20):
                variant = mutate_peptide(seed, hd)
                if variant not in peptides and variant != seed:
                    peptides.append(variant)
                    # Realistic EL rank distribution
                    el_ranks.append(round(random.uniform(0.01, 2.0), 3))
                    gi = random.randint(0, len(fake_genes) - 1)
                    gene_symbols_list.append([fake_genes[gi]])
                    source_proteins_list.append([fake_proteins[gi]])

    # Fill remaining with random peptides
    while len(peptides) < n_peptides:
        length = random.choice([8, 9, 10, 11])
        pep = random_peptide(length)
        if pep not in peptides:
            peptides.append(pep)
            el_ranks.append(round(random.uniform(0.01, 2.0), 3))
            gi = random.randint(0, len(fake_genes) - 1)
            gene_symbols_list.append([fake_genes[gi]])
            source_proteins_list.append([fake_proteins[gi]])

    df = pd.DataFrame({
        "sequence": peptides[:n_peptides],
        "el_rank": el_ranks[:n_peptides],
        "gene_symbols": gene_symbols_list[:n_peptides],
        "source_proteins": source_proteins_list[:n_peptides],
    })

    # Save as parquet
    allele_tag = hla_allele.replace("*", "").replace(":", "")
    cache_path = AB_DATA_DIR / f"hla_filtered_{allele_tag}.parquet"
    df.to_parquet(cache_path, index=False)

    print(f"Demo HLA-filtered data saved: {cache_path}")
    print(f"  Peptides: {len(df)}")
    print(f"  9-mers: {len(df[df['sequence'].str.len() == 9])}")
    print(f"  Strong binders (EL < 0.5): {len(df[df['el_rank'] <= 0.5])}")
    print(f"  Weak binders (EL < 2.0): {len(df[df['el_rank'] <= 2.0])}")

    return cache_path


def build_demo_expression() -> Path:
    """Build a demo expression database."""
    from decoy_a.config import EXPRESSION_DB_PATH, VITAL_ORGANS

    if EXPRESSION_DB_PATH.exists():
        print(f"Expression database already exists: {EXPRESSION_DB_PATH}")
        return EXPRESSION_DB_PATH

    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)

    genes = [
        "TAP2", "PKD1L1", "MAGEA6", "FGD5", "MAGEA12", "MAGEA8",
        "EPS8L2", "TITIN", "ABCB1", "BRCA1", "TP53", "EGFR",
        "MYC", "KRAS", "PTEN", "RB1", "APC", "VHL", "WT1", "NF1",
    ]

    tissues = VITAL_ORGANS + ["testis", "placenta", "skin", "bone marrow", "thymus"]

    rows = []
    for gene in genes:
        max_vital = 0.0
        for tissue in tissues:
            tpm = round(random.uniform(0, 50), 1)
            if tissue in VITAL_ORGANS:
                max_vital = max(max_vital, tpm)

            if max_vital >= 10:
                category = "high_risk"
            elif max_vital >= 1:
                category = "expressed"
            else:
                category = "silent"

            rows.append({
                "gene_symbol": gene,
                "tissue": tissue,
                "tpm": tpm,
                "max_vital_organ_tpm": max_vital,
                "vital_organ_expressed": max_vital >= 1.0,
                "expression_category": category,
            })

    df = pd.DataFrame(rows)
    df.to_parquet(EXPRESSION_DB_PATH, index=False)
    print(f"Demo expression database saved: {EXPRESSION_DB_PATH}")
    return EXPRESSION_DB_PATH


if __name__ == "__main__":
    random.seed(42)
    np.random.seed(42)

    print("Building demo HLA-filtered data...")
    build_demo_hla_filtered()

    print("\nBuilding demo expression database...")
    build_demo_expression()

    print("\nDone! You can now run:")
    print("  python -m decoy_b scan-b --target GILGFVFTL --skip-structural")
    print("  python -m decoy_b scan-b --target EVDPIGHLY --hla HLA-A*01:01 --skip-structural")
    print("  python -m decoy_b run --target GILGFVFTL --hla HLA-A*02:01 --skip-structural --top-n 50")
