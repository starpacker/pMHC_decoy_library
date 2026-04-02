#!/usr/bin/env python3
"""
Decoy A — Standalone Runner
============================
End-to-end Decoy A pipeline: from Swiss-Prot download to final ranked decoy list.

Supports three modes:
  1. Full mode:  Download proteome → build k-mers → HLA filter → Hamming scan
  2. Local mode: Use pre-built k-mer / HLA-filtered Parquet
  3. Demo mode:  Built-in known cross-reactive peptides for quick validation

Usage:
    # Demo (immediate, no downloads)
    python -m decoy_ab.run_decoy_a --target EVDPIGHLY --demo

    # Full pipeline (requires NetMHCpan or IEDB fallback)
    python -m decoy_ab.run_decoy_a --target EVDPIGHLY --hla HLA-A*02:01

    # From custom FASTA/CSV peptide pool
    python -m decoy_ab.run_decoy_a --target EVDPIGHLY --pool my_peptides.csv

    # Step-by-step (build databases first, scan later)
    python -m decoy_ab.run_decoy_a --step build-kmer
    python -m decoy_ab.run_decoy_a --step hla-filter --hla HLA-A*02:01
    python -m decoy_ab.run_decoy_a --target EVDPIGHLY --hla HLA-A*02:01
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from .config import (
    AB_DATA_DIR,
    DECOY_A_OUTPUT,
    DEFAULT_HLA_ALLELE,
    HAMMING_DISTANCE_MAX,
    VITAL_ORGANS,
)
from .scanner import annotate_mismatches, hamming_distance_vectorised, scan_decoy_a
from .models import DecoyAHit, TissueExpression

log = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════════
#  Demo peptide pool — known cross-reactive peptides for validation
# ═════════════════════════════════════════════════════════════════════════

# These are well-documented pMHC off-target cases from literature.
# Used for quick validation that the Hamming scan works correctly.
DEMO_PEPTIDES = [
    # MAGE-A3 / HLA-A*01:01 cross-reactivity cluster
    {"sequence": "EVDPIGHLY", "gene_symbols": ["MAGEA3"], "source_proteins": ["P43357"], "el_rank": 0.15},
    {"sequence": "ESDPIVAQY", "gene_symbols": ["TTN"],    "source_proteins": ["Q8WZ42"], "el_rank": 0.32},
    {"sequence": "EVDPIGHVY", "gene_symbols": ["EPS8L2"], "source_proteins": ["Q9H6S3"], "el_rank": 0.28},
    {"sequence": "EVDPIGHAY", "gene_symbols": ["DEMO1"],  "source_proteins": ["DEMO"],   "el_rank": 0.40},
    {"sequence": "EVDPIGHLF", "gene_symbols": ["DEMO2"],  "source_proteins": ["DEMO"],   "el_rank": 0.55},
    {"sequence": "EVDPIGHLS", "gene_symbols": ["DEMO3"],  "source_proteins": ["DEMO"],   "el_rank": 0.60},
    {"sequence": "EVDAIGHLY", "gene_symbols": ["DEMO4"],  "source_proteins": ["DEMO"],   "el_rank": 0.70},
    {"sequence": "EVDPIAHLY", "gene_symbols": ["DEMO5"],  "source_proteins": ["DEMO"],   "el_rank": 0.45},
    {"sequence": "AVDPIGHLY", "gene_symbols": ["DEMO6"],  "source_proteins": ["DEMO"],   "el_rank": 0.50},
    {"sequence": "EVDPIGKLY", "gene_symbols": ["DEMO7"],  "source_proteins": ["DEMO"],   "el_rank": 0.35},
    {"sequence": "EVDPNGHLY", "gene_symbols": ["DEMO8"],  "source_proteins": ["DEMO"],   "el_rank": 0.65},
    # 2-mismatch variants
    {"sequence": "AVDPIGHVY", "gene_symbols": ["DEMO9"],  "source_proteins": ["DEMO"],   "el_rank": 0.80},
    {"sequence": "EVDAIAHLY", "gene_symbols": ["DEMO10"], "source_proteins": ["DEMO"],   "el_rank": 0.90},
    {"sequence": "EVDPKGHLS", "gene_symbols": ["DEMO11"], "source_proteins": ["DEMO"],   "el_rank": 0.75},
    {"sequence": "AVDPIGHLF", "gene_symbols": ["DEMO12"], "source_proteins": ["DEMO"],   "el_rank": 0.85},
    # >2 mismatch (should NOT be caught by Decoy A)
    {"sequence": "AADAIAHLY", "gene_symbols": ["FAR1"],   "source_proteins": ["DEMO"],   "el_rank": 1.20},
    {"sequence": "QRTPIGHLY", "gene_symbols": ["FAR2"],   "source_proteins": ["DEMO"],   "el_rank": 1.50},
    # Melanoma / HLA-A*02:01 cluster (MART-1/Melan-A)
    {"sequence": "ELAGIGILTV", "gene_symbols": ["MLANA"], "source_proteins": ["Q16655"], "el_rank": 0.10},
    {"sequence": "EAAGIGILTV", "gene_symbols": ["MLANA"], "source_proteins": ["Q16655"], "el_rank": 0.12},
    {"sequence": "ELAGIGALTV", "gene_symbols": ["DEMO13"], "source_proteins": ["DEMO"],  "el_rank": 0.50},
    {"sequence": "ELAGIAILTV", "gene_symbols": ["DEMO14"], "source_proteins": ["DEMO"],  "el_rank": 0.55},
    # NY-ESO-1 / HLA-A*02:01 cluster
    {"sequence": "SLLMWITQC", "gene_symbols": ["CTAG1B"], "source_proteins": ["P78358"], "el_rank": 0.20},
    {"sequence": "SLLMWITQA", "gene_symbols": ["DEMO15"], "source_proteins": ["DEMO"],   "el_rank": 0.45},
    {"sequence": "SLLMWIVQC", "gene_symbols": ["DEMO16"], "source_proteins": ["DEMO"],   "el_rank": 0.50},
    # WT1 / HLA-A*02:01 cluster
    {"sequence": "RMFPNAPYL", "gene_symbols": ["WT1"],    "source_proteins": ["P19544"], "el_rank": 0.08},
    {"sequence": "RMFPNAPYL", "gene_symbols": ["WT1"],    "source_proteins": ["P19544"], "el_rank": 0.08},  # dup (tests dedup)
    {"sequence": "RMFPNAPYA", "gene_symbols": ["DEMO17"], "source_proteins": ["DEMO"],   "el_rank": 0.50},
    {"sequence": "RMFANAPYL", "gene_symbols": ["DEMO18"], "source_proteins": ["DEMO"],   "el_rank": 0.60},
]


def _build_demo_df() -> pd.DataFrame:
    """Build a demo HLA-filtered DataFrame from known peptides."""
    rows = []
    seen = set()
    for p in DEMO_PEPTIDES:
        if p["sequence"] in seen:
            continue
        seen.add(p["sequence"])
        rows.append(p)
    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════
#  Load peptide pool from custom file
# ═════════════════════════════════════════════════════════════════════════

def load_peptide_pool(path: str) -> pd.DataFrame:
    """
    Load a peptide pool from CSV, TSV, or FASTA.

    CSV/TSV must have a 'sequence' column. Optional: gene_symbols, source_proteins, el_rank.
    FASTA: headers parsed for gene info, el_rank defaults to 1.0.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Peptide pool not found: {path}")

    suffix = p.suffix.lower()
    if suffix in (".csv", ".tsv"):
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(p, sep=sep)
        if "sequence" not in df.columns:
            raise ValueError("CSV/TSV must have a 'sequence' column")
        if "el_rank" not in df.columns:
            df["el_rank"] = 1.0
        if "gene_symbols" not in df.columns:
            df["gene_symbols"] = [[] for _ in range(len(df))]
        if "source_proteins" not in df.columns:
            df["source_proteins"] = [[] for _ in range(len(df))]
        return df

    if suffix in (".fasta", ".fa", ".faa"):
        sequences = []
        current_gene = ""
        current_seq_parts: list[str] = []

        text = p.read_text(encoding="utf-8")
        for line in text.splitlines():
            if line.startswith(">"):
                if current_seq_parts:
                    seq = "".join(current_seq_parts).strip().upper()
                    if 8 <= len(seq) <= 15:
                        sequences.append({
                            "sequence": seq,
                            "gene_symbols": [current_gene] if current_gene else [],
                            "source_proteins": [],
                            "el_rank": 1.0,
                        })
                current_gene = line[1:].split()[0] if len(line) > 1 else ""
                current_seq_parts = []
            else:
                current_seq_parts.append(line.strip())

        if current_seq_parts:
            seq = "".join(current_seq_parts).strip().upper()
            if 8 <= len(seq) <= 15:
                sequences.append({
                    "sequence": seq,
                    "gene_symbols": [current_gene] if current_gene else [],
                    "source_proteins": [],
                    "el_rank": 1.0,
                })

        if not sequences:
            raise ValueError("No valid peptides (8-15 AA) found in FASTA")
        return pd.DataFrame(sequences)

    raise ValueError(f"Unsupported file format: {suffix}. Use .csv, .tsv, or .fasta")


# ═════════════════════════════════════════════════════════════════════════
#  Result formatting
# ═════════════════════════════════════════════════════════════════════════

def format_results_table(
    target: str,
    hits: List[DecoyAHit],
    hla_allele: str,
    max_hamming: int,
) -> str:
    """Format Decoy A results as a readable table."""
    lines = []
    lines.append("")
    lines.append("=" * 90)
    lines.append("  Decoy A Results — Sequence Homology Scan")
    lines.append("=" * 90)
    lines.append(f"  Target peptide : {target}")
    lines.append(f"  HLA allele     : {hla_allele}")
    lines.append(f"  Max Hamming    : {max_hamming}")
    lines.append(f"  Total hits     : {len(hits)}")

    n1 = sum(1 for h in hits if h.hamming_distance == 1)
    n2 = sum(1 for h in hits if h.hamming_distance == 2)
    lines.append(f"    1-mismatch   : {n1}")
    lines.append(f"    2-mismatch   : {n2}")
    lines.append("-" * 90)

    if not hits:
        lines.append(f"  No hits found within Hamming distance <= {max_hamming}")
        lines.append("=" * 90)
        return "\n".join(lines)

    # Header
    lines.append(
        f"  {'#':<4} {'Sequence':<15} {'HD':>3} {'EL%Rank':>8} "
        f"{'Genes':<12} {'Mismatches':<25} {'TCR_mis':>7} {'Anc_mis':>7}"
    )
    lines.append("  " + "-" * 86)

    for i, h in enumerate(hits, 1):
        genes = ",".join(h.gene_symbols[:2]) if h.gene_symbols else "-"
        # Format mismatch details
        mm_parts = []
        for m in h.mismatches:
            flag = ""
            if m.is_tcr_contact:
                flag = "*"  # TCR contact
            elif m.is_anchor:
                flag = "^"  # Anchor
            mm_parts.append(f"p{m.position+1}:{m.target_aa}>{m.candidate_aa}{flag}")
        mm_str = " ".join(mm_parts) if mm_parts else "-"

        lines.append(
            f"  {i:<4} {h.sequence:<15} {h.hamming_distance:>3} {h.el_rank:>8.2f} "
            f"{genes:<12} {mm_str:<25} {h.n_tcr_contact_mismatches:>7} {h.n_anchor_mismatches:>7}"
        )

    lines.append("")
    lines.append("  Legend: * = TCR contact position, ^ = HLA anchor position")
    lines.append("=" * 90)
    return "\n".join(lines)


def save_results(
    hits: List[DecoyAHit],
    output_path: Optional[Path] = None,
) -> Path:
    """Save Decoy A results to JSON."""
    if output_path is None:
        output_path = DECOY_A_OUTPUT
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    data = [h.model_dump(mode="json") for h in hits]
    output_path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return output_path


# ═════════════════════════════════════════════════════════════════════════
#  Step runners
# ═════════════════════════════════════════════════════════════════════════

def step_build_kmer(force: bool = False) -> None:
    """Step 1: Download Swiss-Prot and build k-mer + expression databases."""
    from .kmer_builder import build_expression_database, build_kmer_database

    print("=" * 60)
    print("  Step 1: Building K-mer Dictionary + Expression Database")
    print("=" * 60)

    t0 = time.time()
    kmer_path = build_kmer_database(force=force)
    df = pd.read_parquet(kmer_path)
    print(f"  K-mer database: {len(df):,} unique peptides")
    print(f"  Saved to: {kmer_path}")

    expr_path = build_expression_database(force=force)
    edf = pd.read_parquet(expr_path)
    n_genes = edf["gene_symbol"].nunique()
    print(f"  Expression database: {n_genes:,} genes")
    print(f"  Saved to: {expr_path}")
    print(f"  Elapsed: {time.time() - t0:.1f}s")
    print("=" * 60)


def step_hla_filter(hla_allele: str = DEFAULT_HLA_ALLELE, force: bool = False) -> None:
    """Step 2: HLA presentation gate."""
    from .hla_filter import run_hla_filter

    print("=" * 60)
    print(f"  Step 2: HLA Presentation Gate ({hla_allele})")
    print("=" * 60)

    t0 = time.time()
    df = run_hla_filter(hla_allele=hla_allele, force=force)
    n_strong = len(df[df["binding"] == "Strong_Binder"]) if "binding" in df.columns else 0
    n_weak = len(df[df["binding"] == "Weak_Binder"]) if "binding" in df.columns else 0
    print(f"  Total binders: {len(df):,}")
    print(f"    Strong (%Rank <= 0.5): {n_strong:,}")
    print(f"    Weak   (%Rank <= 2.0): {n_weak:,}")
    print(f"  Elapsed: {time.time() - t0:.1f}s")
    print("=" * 60)


def run_decoy_a(
    target_sequence: str,
    hla_allele: str = DEFAULT_HLA_ALLELE,
    max_hamming: int = HAMMING_DISTANCE_MAX,
    pool_path: Optional[str] = None,
    demo: bool = False,
    output: Optional[str] = None,
) -> List[DecoyAHit]:
    """
    Run the Decoy A pipeline.

    Parameters
    ----------
    target_sequence : str
        Target peptide to scan against.
    hla_allele : str
        HLA allele context.
    max_hamming : int
        Maximum Hamming distance.
    pool_path : str, optional
        Path to custom peptide pool (CSV/TSV/FASTA).
    demo : bool
        Use built-in demo peptide pool.
    output : str, optional
        Output JSON path.
    """
    target_sequence = target_sequence.strip().upper()
    print("")
    print("=" * 90)
    print("  Decoy A — Sequence Homology Scan")
    print("=" * 90)
    print(f"  Target       : {target_sequence}")
    print(f"  HLA          : {hla_allele}")
    print(f"  Max Hamming  : {max_hamming}")

    t0 = time.time()

    if demo:
        print("  Mode         : DEMO (built-in peptide pool)")
        hla_filtered_df = _build_demo_df()
        print(f"  Pool size    : {len(hla_filtered_df)} peptides")
    elif pool_path:
        print(f"  Mode         : Custom pool ({pool_path})")
        hla_filtered_df = load_peptide_pool(pool_path)
        print(f"  Pool size    : {len(hla_filtered_df)} peptides")
    else:
        print("  Mode         : Full pipeline (pre-built HLA-filtered database)")
        hla_filtered_df = None  # scan_decoy_a will load from cache

    print("-" * 90)

    # Run the scan
    hits = scan_decoy_a(
        target_sequence=target_sequence,
        hla_allele=hla_allele,
        hla_filtered_df=hla_filtered_df,
        max_hamming=max_hamming,
    )

    # Display results
    table = format_results_table(target_sequence, hits, hla_allele, max_hamming)
    print(table)

    # Save results
    out_path = Path(output) if output else DECOY_A_OUTPUT
    saved = save_results(hits, out_path)
    print(f"\n  Results saved to: {saved}")
    print(f"  Elapsed: {time.time() - t0:.1f}s")

    return hits


# ═════════════════════════════════════════════════════════════════════════
#  CLI
# ═════════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="decoy_a",
        description="Decoy A — Sequence Homology Off-Target Screening",
    )
    parser.add_argument("-v", "--verbose", action="store_true")

    sub = parser.add_subparsers(dest="command")

    # -- scan (default behavior) ------------------------------------
    p_scan = sub.add_parser("scan", help="Run Decoy A Hamming scan")
    p_scan.add_argument("--target", required=True, help="Target peptide sequence")
    p_scan.add_argument("--hla", default=DEFAULT_HLA_ALLELE, help="HLA allele")
    p_scan.add_argument("--max-hamming", type=int, default=HAMMING_DISTANCE_MAX)
    p_scan.add_argument("--pool", help="Custom peptide pool (CSV/TSV/FASTA)")
    p_scan.add_argument("--demo", action="store_true", help="Use demo peptide pool")
    p_scan.add_argument("--output", "-o", help="Output JSON path")

    # -- build-kmer -------------------------------------------------
    p_kmer = sub.add_parser("build-kmer", help="Step 1: Build k-mer + expression DB")
    p_kmer.add_argument("--force", action="store_true")

    # -- hla-filter -------------------------------------------------
    p_hla = sub.add_parser("hla-filter", help="Step 2: HLA presentation gate")
    p_hla.add_argument("--hla", default=DEFAULT_HLA_ALLELE)
    p_hla.add_argument("--force", action="store_true")

    # -- demo (shortcut) --------------------------------------------
    p_demo = sub.add_parser("demo", help="Quick demo with built-in peptides")
    p_demo.add_argument("--target", default="EVDPIGHLY", help="Target peptide")
    p_demo.add_argument("--hla", default="HLA-A*01:01")

    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.command == "build-kmer":
        step_build_kmer(force=args.force)
    elif args.command == "hla-filter":
        step_hla_filter(hla_allele=args.hla, force=args.force)
    elif args.command == "demo":
        run_decoy_a(
            target_sequence=args.target,
            hla_allele=args.hla,
            demo=True,
        )
    elif args.command == "scan":
        run_decoy_a(
            target_sequence=args.target,
            hla_allele=args.hla,
            max_hamming=args.max_hamming,
            pool_path=args.pool,
            demo=args.demo,
            output=args.output,
        )
    else:
        # No subcommand: if --target given, do scan; else show help
        # Re-parse with scan defaults
        parser.print_help()
        print("\nQuick start:")
        print("  python -m decoy_ab.run_decoy_a demo")
        print("  python -m decoy_ab.run_decoy_a scan --target EVDPIGHLY --demo")
        print("  python -m decoy_ab.run_decoy_a scan --target EVDPIGHLY --pool peptides.csv")
        print("  python -m decoy_ab.run_decoy_a build-kmer")
        print("  python -m decoy_ab.run_decoy_a hla-filter --hla HLA-A*02:01")


if __name__ == "__main__":
    main()
