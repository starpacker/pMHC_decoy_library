#!/usr/bin/env python3
"""
Decoy A/B Pipeline CLI
======================
Command-line interface for the computational off-target screening pipeline.

Usage examples:

    # Full pipeline run
    python -m decoy_ab run --target EVDPIGHLY --hla HLA-A*02:01

    # Run individual steps
    python -m decoy_ab build-kmer
    python -m decoy_ab hla-filter --hla HLA-A*02:01
    python -m decoy_ab scan-a --target EVDPIGHLY
    python -m decoy_ab scan-b --target EVDPIGHLY
    python -m decoy_ab score

    # Show results
    python -m decoy_ab show
    python -m decoy_ab stats
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


# =====================================================================
#  Subcommands
# =====================================================================

def cmd_run(args: argparse.Namespace) -> None:
    """Run the full Decoy A/B pipeline end-to-end."""
    from .orchestrator import run_full_pipeline

    result = run_full_pipeline(
        target_sequence=args.target,
        hla_allele=args.hla,
        target_protein=args.protein or "",
        force=args.force,
        skip_structural=args.skip_structural,
        run_mpnn=args.mpnn,
        top_n=args.top_n,
        resume=args.resume,
    )
    print(result.summary())


def cmd_build_kmer(args: argparse.Namespace) -> None:
    """Build the human k-mer dictionary and expression database."""
    from decoy_a.kmer_builder import build_expression_database, build_kmer_database

    print("Building k-mer database...")
    kmer_path = build_kmer_database(force=args.force)
    print("K-mer database: {}".format(kmer_path))

    print("Building expression database...")
    expr_path = build_expression_database(force=args.force)
    print("Expression database: {}".format(expr_path))
    print("Done.")


def cmd_hla_filter(args: argparse.Namespace) -> None:
    """Run the HLA presentation gate."""
    from decoy_a.hla_filter import run_hla_filter

    print("Running HLA filter for {}...".format(args.hla))
    df = run_hla_filter(hla_allele=args.hla, force=args.force)
    print("HLA filter complete: {} binders".format(len(df)))


def cmd_scan_a(args: argparse.Namespace) -> None:
    """Run Decoy A sequence homology scan."""
    from decoy_a.scanner import scan_decoy_a

    print("Running Decoy A scan for target: {}".format(args.target))
    hits = scan_decoy_a(
        target_sequence=args.target,
        hla_allele=args.hla,
        max_hamming=args.max_hamming,
    )
    print("Decoy A: {} hits found".format(len(hits)))
    if hits:
        print("\nTop 10 hits:")
        header = "{:<15} {:>4} {:>8} {:>10} {}".format(
            "Sequence", "HD", "EL_Rank", "Genes", "TCR_mis"
        )
        print(header)
        print("-" * 60)
        for h in hits[:10]:
            genes = ",".join(h.gene_symbols[:2]) if h.gene_symbols else "-"
            row = "{:<15} {:>4} {:>8.2f} {:>10} {}".format(
                h.sequence,
                h.hamming_distance,
                h.el_rank,
                genes,
                h.n_tcr_contact_mismatches,
            )
            print(row)


def cmd_scan_b(args: argparse.Namespace) -> None:
    """Run Decoy B structural similarity scan."""
    from .scanner import scan_decoy_b

    print("Running Decoy B scan for target: {}".format(args.target))
    hits = scan_decoy_b(
        target_sequence=args.target,
        hla_allele=args.hla,
        run_structural=not args.skip_structural,
        run_mpnn=args.mpnn,
    )
    print("Decoy B: {} hits found".format(len(hits)))

    # Separate MPNN tier stats
    n_mpnn_t1 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "proteome_matched")
    n_mpnn_t2 = sum(1 for h in hits if getattr(h, "mpnn_source", "") == "hla_qualified_synthetic")
    if n_mpnn_t1 or n_mpnn_t2:
        print("  MPNN Tier 1 (proteome-matched): {}".format(n_mpnn_t1))
        print("  MPNN Tier 2 (HLA-qualified):    {}".format(n_mpnn_t2))
    if hits:
        print("\nTop 10 hits:")
        header = "{:<15} {:>8} {:>10} {}".format(
            "Sequence", "CosSim", "EL_Rank", "Genes"
        )
        print(header)
        print("-" * 55)
        for h in hits[:10]:
            genes = ",".join(h.gene_symbols[:2]) if h.gene_symbols else "-"
            row = "{:<15} {:>8.3f} {:>10.2f} {}".format(
                h.sequence,
                h.physicochemical.cosine_similarity,
                h.el_rank,
                genes,
            )
            print(row)


def cmd_score(args: argparse.Namespace) -> None:
    """Run risk scoring on existing A/B results."""
    from .orchestrator import load_state, run_step

    state = load_state()
    if state is None:
        print("No pipeline state found. Run the full pipeline first.")
        sys.exit(1)

    if not state.decoy_a_hits and not state.decoy_b_hits:
        print("No Decoy A/B hits found in state. Run scan-a and scan-b first.")
        sys.exit(1)

    run_step(
        "score",
        target_sequence=state.target_sequence,
        hla_allele=state.hla_allele,
        top_n=args.top_n,
    )
    print("Scoring complete. Results saved.")


def cmd_show(args: argparse.Namespace) -> None:
    """Display the final ranked results."""
    from .risk_scorer import load_ranked_results

    try:
        entries = load_ranked_results()
    except FileNotFoundError:
        print("No ranked results found. Run the pipeline first.")
        sys.exit(1)

    if args.format == "json":
        data = [e.model_dump(mode="json") for e in entries]
        print(json.dumps(data, indent=2, ensure_ascii=False))
        return

    # Table format
    header = "{:<10} {:<15} {:>4} {:>8} {:>8} {:<8} {:>10} {}".format(
        "ID", "Sequence", "HD", "EL_Rank", "Risk", "Source", "TPM_Wt", "Genes"
    )
    print(header)
    print("-" * 100)

    for e in entries:
        genes = ",".join(e.gene_symbols[:2]) if e.gene_symbols else "-"
        src = "A" if "Sequence" in e.source.value else ("B" if "Structural" in e.source.value else "A+B")
        row = "{:<10} {:<15} {:>4} {:>8.2f} {:>8.4f} {:<8} {:>10.1f} {}".format(
            e.decoy_ab_id,
            e.sequence,
            e.hamming_distance,
            e.el_rank,
            e.total_risk_score,
            src,
            e.vital_organ_tpm_weight,
            genes,
        )
        print(row)

    print("\nTotal: {} entries".format(len(entries)))


def cmd_stats(args: argparse.Namespace) -> None:
    """Show pipeline statistics."""
    from .orchestrator import load_state
    from .risk_scorer import load_ranked_results

    state = load_state()
    if state is None:
        print("No pipeline state found.")
        sys.exit(1)

    print("Decoy A/B Pipeline Statistics")
    print("=" * 40)
    print("Target peptide  : {}".format(state.target_sequence))
    print("Target protein  : {}".format(state.target_protein or "-"))
    print("HLA allele      : {}".format(state.hla_allele))
    print("")
    print("Step 1: K-mers generated   : {:,}".format(state.total_kmers_generated))
    print("Step 1: Genes w/ expression: {:,}".format(state.genes_with_expression))
    print("Step 2: HLA-filtered       : {:,}".format(state.hla_filtered_count))
    print("Step 3A: Decoy A hits      : {:,}".format(len(state.decoy_a_hits)))
    print("Step 3B: Decoy B hits      : {:,}".format(len(state.decoy_b_hits)))
    print("  - Physchem only          : {:,}".format(state.physicochemical_candidates))
    print("  - With structure         : {:,}".format(state.structural_candidates))
    print("Step 4: Final entries      : {:,}".format(state.final_top_n))

    try:
        entries = load_ranked_results()
        if entries:
            from decoy_a.models import DecoySource
            n_a = sum(1 for e in entries if e.source == DecoySource.DECOY_A)
            n_b = sum(1 for e in entries if e.source == DecoySource.DECOY_B)
            n_both = sum(1 for e in entries if e.source == DecoySource.BOTH)
            print("")
            print("Final results breakdown:")
            print("  Decoy A only  : {}".format(n_a))
            print("  Decoy B only  : {}".format(n_b))
            print("  Both A and B  : {}".format(n_both))
            print("  Top risk score: {:.4f} ({})".format(
                entries[0].total_risk_score, entries[0].sequence
            ))
    except FileNotFoundError:
        pass


def cmd_clear(args: argparse.Namespace) -> None:
    """Clear pipeline checkpoint."""
    from .orchestrator import clear_state
    clear_state()
    print("Pipeline checkpoint cleared.")


# =====================================================================
#  Main parser
# =====================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="decoy_ab",
        description="Decoy A/B -- Computational Off-Target Screening Pipeline",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    sub = parser.add_subparsers(dest="command", help="Available commands")

    # -- run (full pipeline) ----------------------------------------
    p_run = sub.add_parser("run", help="Run the full Decoy A/B pipeline")
    p_run.add_argument("--target", required=True, help="Target peptide sequence")
    p_run.add_argument("--hla", default="HLA-A*02:01", help="HLA allele")
    p_run.add_argument("--protein", default="", help="Target protein name")
    p_run.add_argument("--force", action="store_true", help="Force re-run all steps")
    p_run.add_argument("--skip-structural", action="store_true", help="Skip 3D modeling")
    p_run.add_argument("--mpnn", action="store_true", help="Enable MPNN inverse design branch")
    p_run.add_argument("--top-n", type=int, default=100, help="Number of top results")
    p_run.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    p_run.set_defaults(func=cmd_run)

    # -- build-kmer -------------------------------------------------
    p_kmer = sub.add_parser("build-kmer", help="Build k-mer + expression databases")
    p_kmer.add_argument("--force", action="store_true", help="Force rebuild")
    p_kmer.set_defaults(func=cmd_build_kmer)

    # -- hla-filter -------------------------------------------------
    p_hla = sub.add_parser("hla-filter", help="Run HLA presentation gate")
    p_hla.add_argument("--hla", default="HLA-A*02:01", help="HLA allele")
    p_hla.add_argument("--force", action="store_true", help="Force re-run")
    p_hla.set_defaults(func=cmd_hla_filter)

    # -- scan-a -----------------------------------------------------
    p_a = sub.add_parser("scan-a", help="Run Decoy A sequence homology scan")
    p_a.add_argument("--target", required=True, help="Target peptide sequence")
    p_a.add_argument("--hla", default="HLA-A*02:01", help="HLA allele")
    p_a.add_argument("--max-hamming", type=int, default=2, help="Max Hamming distance")
    p_a.set_defaults(func=cmd_scan_a)

    # -- scan-b -----------------------------------------------------
    p_b = sub.add_parser("scan-b", help="Run Decoy B structural similarity scan")
    p_b.add_argument("--target", required=True, help="Target peptide sequence")
    p_b.add_argument("--hla", default="HLA-A*02:01", help="HLA allele")
    p_b.add_argument("--skip-structural", action="store_true", help="Skip 3D modeling")
    p_b.add_argument("--mpnn", action="store_true", help="Enable MPNN inverse design branch")
    p_b.set_defaults(func=cmd_scan_b)

    # -- score ------------------------------------------------------
    p_sc = sub.add_parser("score", help="Run risk scoring on existing results")
    p_sc.add_argument("--top-n", type=int, default=100, help="Number of top results")
    p_sc.set_defaults(func=cmd_score)

    # -- show -------------------------------------------------------
    p_sh = sub.add_parser("show", help="Display final ranked results")
    p_sh.add_argument("--format", choices=["json", "table"], default="table")
    p_sh.set_defaults(func=cmd_show)

    # -- stats ------------------------------------------------------
    p_st = sub.add_parser("stats", help="Show pipeline statistics")
    p_st.set_defaults(func=cmd_stats)

    # -- clear ------------------------------------------------------
    p_cl = sub.add_parser("clear", help="Clear pipeline checkpoint")
    p_cl.set_defaults(func=cmd_clear)

    args = parser.parse_args()
    setup_logging(args.verbose)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
