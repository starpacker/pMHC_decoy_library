#!/usr/bin/env python3
"""
Unified CLI entry point for running individual decoy strategies.

Usage:
    python run_decoy.py GILGFVFTL d                  # Run Decoy D only
    python run_decoy.py GILGFVFTL a                  # Run Decoy A only
    python run_decoy.py GILGFVFTL b                  # Run Decoy B only
    python run_decoy.py GILGFVFTL a b d              # Run A, B, D
    python run_decoy.py GILGFVFTL d --hla HLA-A*01:01
    python run_decoy.py GILGFVFTL d --designs 2000 --top-k 20
    python run_decoy.py GILGFVFTL a b d --skip-structural  # B without structure prediction
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("PYTHONUTF8", "1")

PROJECT_ROOT = Path(__file__).resolve().parent


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


# ── Strategy runners ─────────────────────────────────────────────────────

def run_strategy_a(target: str, hla: str, **kwargs) -> None:
    """Decoy A: sequence homology scan."""
    from decoy_b.scanner import scan_decoy_a

    log = logging.getLogger("run_decoy.A")
    log.info("=== Decoy A: Sequence Homology Scan ===")
    log.info("Target: %s | HLA: %s", target, hla)

    hits = scan_decoy_a(target, hla)
    out_dir = PROJECT_ROOT / "data" / "decoy_a" / target
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "decoy_a_results.json"
    results = [h.model_dump() if hasattr(h, "model_dump") else h.dict() for h in hits]
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    log.info("Decoy A: %d hits -> %s", len(hits), out_file)


def run_strategy_b(target: str, hla: str, *, skip_structural: bool = False, **kwargs) -> None:
    """Decoy B: physicochemical + structural similarity scan."""
    from decoy_b.scanner import scan_decoy_b

    log = logging.getLogger("run_decoy.B")
    log.info("=== Decoy B: Structural Similarity Scan ===")
    log.info("Target: %s | HLA: %s | skip_structural: %s", target, hla, skip_structural)

    hits = scan_decoy_b(target, hla, run_structural=not skip_structural)
    out_dir = PROJECT_ROOT / "data" / "decoy_b" / target
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "decoy_b_results.json"
    results = [h.model_dump() if hasattr(h, "model_dump") else h.dict() for h in hits]
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)
    log.info("Decoy B: %d hits -> %s", len(hits), out_file)


def run_strategy_d(
    target: str,
    hla: str,
    *,
    designs: int = 1000,
    top_k: int = 10,
    target_pdb: str | None = None,
    **kwargs,
) -> None:
    """Decoy D: ProteinMPNN inverse design + mhcflurry filtering."""
    from decoy_d.scanner import run_decoy_d

    log = logging.getLogger("run_decoy.D")
    log.info("=== Decoy D: MPNN Inverse Design ===")
    log.info("Target: %s | HLA: %s | designs: %d | top_k: %d", target, hla, designs, top_k)

    output_dir = PROJECT_ROOT / "data" / "decoy_d" / target
    df = run_decoy_d(
        target_sequence=target,
        hla_allele=hla,
        target_pdb=target_pdb,
        num_designs=designs,
        top_k_structures=top_k,
        output_dir=output_dir,
    )

    log.info("Decoy D: %d candidates passed filter", len(df))
    if not df.empty:
        print(f"\n{'='*60}")
        print(f" Decoy D Results: {target} ({hla})")
        print(f"{'='*60}")
        print(df.head(top_k).to_string(index=False))
        print(f"{'='*60}")
        print(f"Full results: {output_dir / 'decoy_d_results.csv'}")


STRATEGY_MAP = {
    "a": run_strategy_a,
    "b": run_strategy_b,
    "d": run_strategy_d,
}


# ── CLI ──────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Unified CLI for running decoy strategies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  python run_decoy.py GILGFVFTL d
  python run_decoy.py GILGFVFTL a b d
  python run_decoy.py NLVPMVATV d --hla HLA-A*02:01 --designs 2000
  python run_decoy.py EVDPIGHLY d --hla HLA-A*01:01
""",
    )
    parser.add_argument("target", help="Target peptide sequence (e.g. GILGFVFTL)")
    parser.add_argument(
        "strategies",
        nargs="+",
        choices=["a", "b", "c", "d"],
        help="Decoy strategies to run (a/b/c/d)",
    )
    parser.add_argument("--hla", default="HLA-A*02:01", help="HLA allele (default: HLA-A*02:01)")
    parser.add_argument("--designs", type=int, default=1000, help="[D] Number of MPNN designs (default: 1000)")
    parser.add_argument("--top-k", type=int, default=10, help="[D] Top K structures to predict (default: 10)")
    parser.add_argument("--target-pdb", help="[D] Pre-computed target PDB path")
    parser.add_argument("--skip-structural", action="store_true", help="[B] Skip structural prediction stage")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose (DEBUG) logging")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    setup_logging(args.verbose)
    log = logging.getLogger("run_decoy")

    target = args.target.upper()
    strategies = sorted(set(args.strategies))

    log.info("Target: %s | HLA: %s | Strategies: %s", target, args.hla, strategies)

    for s in strategies:
        if s == "c":
            log.warning("Decoy C (literature mining) is not yet integrated in this CLI — skipping")
            continue

        runner = STRATEGY_MAP.get(s)
        if runner is None:
            log.error("Unknown strategy: %s", s)
            continue

        runner(
            target=target,
            hla=args.hla,
            designs=args.designs,
            top_k=args.top_k,
            target_pdb=args.target_pdb,
            skip_structural=args.skip_structural,
        )

    log.info("All requested strategies completed.")


if __name__ == "__main__":
    main()
