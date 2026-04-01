#!/usr/bin/env python3
"""Decoy C Library CLI Entry Point."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from decoy_c.models import DecoyLibrary
from decoy_c.orchestrator import (
    load_library,
    run_cold_start,
    run_pipeline,
    save_library,
)
from decoy_c.seed_data import build_seed_library
from decoy_c.validator import validate_batch


def setup_logging(verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def cmd_cold_start(args):
    """Initialize the library with seed data and PubMed discovery."""
    if args.seed_only:
        lib = build_seed_library()
        if not args.no_validate:
            logging.info("Validating seed entries...")
            lib.entries = validate_batch(lib.entries)
        save_library(lib)
        n = len(lib.entries)
        print("Seed library created with " + str(n) + " entries.")
    else:
        lib = run_cold_start(validate=not args.no_validate)
        n = len(lib.entries)
        print("Cold start complete: " + str(n) + " entries in library.")


def cmd_search(args):
    """Search PubMed and extract decoys from results."""
    lib = run_pipeline(query=args.query, validate=not args.no_validate)
    n = len(lib.entries)
    print("Search complete: " + str(n) + " total entries in library.")


def cmd_fetch(args):
    """Fetch and process specific PMIDs."""
    lib = run_pipeline(pmids=args.pmids, validate=not args.no_validate)
    n = len(lib.entries)
    print("Fetch complete: " + str(n) + " total entries in library.")


def cmd_validate(args):
    """Re-validate all entries in the existing library."""
    lib = load_library()
    if not lib.entries:
        print("Library is empty. Run cold-start first.")
        return
    print("Re-validating " + str(len(lib.entries)) + " entries...")
    lib.entries = validate_batch(lib.entries)
    save_library(lib)

    statuses = {}
    for e in lib.entries:
        s = e.validation_flags.get("overall_status", "UNKNOWN")
        statuses[s] = statuses.get(s, 0) + 1
    print("Validation summary:")
    for status, count in sorted(statuses.items()):
        print("  " + status + ": " + str(count))


def cmd_show(args):
    """Display library contents."""
    lib = load_library()
    if not lib.entries:
        print("Library is empty. Run cold-start first.")
        return

    if args.format == "json":
        print(lib.model_dump_json(indent=2))
    else:
        header = "{:<10} {:<15} {:<15} {:<10} {:<30} {:<12} {}".format(
            "ID", "Sequence", "HLA", "Gene", "Level", "Status", "Source"
        )
        print(header)
        print("-" * 160)
        for e in lib.entries:
            seq = e.peptide_info.decoy_sequence
            hla = e.peptide_info.hla_allele
            gene = e.peptide_info.gene_symbol
            lvl = e.risk_profile.evidence_level.value
            status = e.validation_flags.get("overall_status", "-")
            src = e.source.citation if e.source and e.source.citation else "-"
            row = "{:<10} {:<15} {:<15} {:<10} {:<30} {:<12} {}".format(
                e.decoy_id, seq, hla, gene, lvl, status, src
            )
            print(row)
        print("\nTotal: " + str(len(lib.entries)) + " entries")


def cmd_stats(args):
    """Show library statistics."""
    lib = load_library()
    if not lib.entries:
        print("Library is empty.")
        return

    levels = {}
    for e in lib.entries:
        lv = e.risk_profile.evidence_level.value
        levels[lv] = levels.get(lv, 0) + 1

    print("Decoy C Library Statistics")
    print("==========================")
    print("Total entries: " + str(len(lib.entries)))
    print("\nBy evidence level:")
    for lv, cnt in sorted(levels.items()):
        print("  " + lv + ": " + str(cnt))

    ms_true = sum(1 for e in lib.entries if e.experimental_evidence.mass_spec_confirmed is True)
    ms_false = sum(1 for e in lib.entries if e.experimental_evidence.mass_spec_confirmed is False)
    ms_none = sum(1 for e in lib.entries if e.experimental_evidence.mass_spec_confirmed is None)
    print("\nMass-spec confirmation:")
    print("  Confirmed: " + str(ms_true))
    print("  Not confirmed: " + str(ms_false))
    print("  Unknown: " + str(ms_none))

    genes = set(e.peptide_info.gene_symbol for e in lib.entries)
    print("\nUnique genes: " + str(len(genes)))
    print("Genes: " + ", ".join(sorted(genes)))


def main():
    parser = argparse.ArgumentParser(
        prog="decoy_library",
        description="Decoy C Library — Cold-Start Collection System",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging")
    sub = parser.add_subparsers(dest="command", help="Available commands")

    # cold-start
    p_cs = sub.add_parser("cold-start", help="Initialize library with seed + PubMed discovery")
    p_cs.add_argument("--no-validate", action="store_true", help="Skip validation")
    p_cs.add_argument("--seed-only", action="store_true", help="Only seed data, no PubMed")
    p_cs.set_defaults(func=cmd_cold_start)

    # search
    p_s = sub.add_parser("search", help="Search PubMed and extract decoys")
    p_s.add_argument("--query", required=True, help="PubMed search query")
    p_s.add_argument("--no-validate", action="store_true")
    p_s.set_defaults(func=cmd_search)

    # fetch
    p_f = sub.add_parser("fetch", help="Fetch specific PMIDs")
    p_f.add_argument("--pmids", nargs="+", required=True, help="PMID list")
    p_f.add_argument("--no-validate", action="store_true")
    p_f.set_defaults(func=cmd_fetch)

    # validate
    p_v = sub.add_parser("validate", help="Re-validate all entries")
    p_v.set_defaults(func=cmd_validate)

    # show
    p_sh = sub.add_parser("show", help="Display library contents")
    p_sh.add_argument("--format", choices=["json", "table"], default="table")
    p_sh.set_defaults(func=cmd_show)

    # stats
    p_st = sub.add_parser("stats", help="Show library statistics")
    p_st.set_defaults(func=cmd_stats)

    args = parser.parse_args()
    setup_logging(args.verbose)

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
