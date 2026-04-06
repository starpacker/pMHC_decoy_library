#!/usr/bin/env python3
"""
Generate Decoy D detail figure for a given target peptide.

Reads:  data/decoy_d/{TARGET}/decoy_d_results.csv
Writes: figures/{TARGET}_decoy_d_detail.png

Usage:
    python scripts/visualize_decoy_d_detail.py GILGFVFTL
    python scripts/visualize_decoy_d_detail.py --all   # All targets in candidate_targets.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "data"
FIGURES_DIR = REPO / "figures"

C_D = "#4CAF50"
C_BG = "#F5F5F5"


def plot_decoy_d_detail(target: str) -> Path | None:
    """Generate the 3-panel Decoy D detail figure for *target*.

    Returns the output path, or None if no data found.
    """
    csv_path = DATA_DIR / "decoy_d" / target / "decoy_d_results.csv"
    if not csv_path.exists():
        print(f"[SKIP] {target}: {csv_path} not found")
        return None

    df = pd.read_csv(csv_path)
    # Exclude the target itself if present
    df = df[df["sequence"] != target].copy()
    if df.empty:
        print(f"[SKIP] {target}: no candidate rows after filtering")
        return None

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.facecolor": C_BG,
        "figure.facecolor": "white",
        "axes.grid": True,
        "grid.alpha": 0.3,
    })

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle(
        f"Decoy D (MPNN Inverse Design) — {target}",
        fontsize=15, fontweight="bold",
    )

    # ── Panel a: Top 15 by MPNN score ────────────────────────────────────
    ax = axes[0]
    top15 = df.nsmallest(15, "mpnn_score")
    y_pos = range(len(top15))
    ax.barh(y_pos, top15["mpnn_score"], color=C_D, edgecolor="white", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top15["sequence"], fontfamily="monospace", fontsize=9)
    ax.set_xlabel("MPNN Score")
    ax.set_title("Top 15 by MPNN Score", fontweight="bold")
    ax.invert_yaxis()
    for i, (_, row) in enumerate(top15.iterrows()):
        ax.text(
            row["mpnn_score"] + 0.01, i,
            f'EL={row["el_rank"]:.2f}', va="center", fontsize=7,
        )

    # ── Panel b: Top 15 by EL%Rank ──────────────────────────────────────
    ax = axes[1]
    top15_el = df.nsmallest(15, "el_rank")
    y_pos2 = range(len(top15_el))
    ax.barh(y_pos2, top15_el["el_rank"], color="#66BB6A", edgecolor="white", height=0.7)
    ax.set_yticks(y_pos2)
    ax.set_yticklabels(top15_el["sequence"], fontfamily="monospace", fontsize=9)
    ax.set_xlabel("EL %Rank")
    ax.set_title("Top 15 by MHC Binding", fontweight="bold")
    ax.invert_yaxis()
    for i, (_, row) in enumerate(top15_el.iterrows()):
        ax.text(
            row["el_rank"] + 0.005, i,
            f'MPNN={row["mpnn_score"]:.2f}', va="center", fontsize=7,
        )

    # ── Panel c: Pareto front (MPNN score vs EL%Rank) ────────────────────
    ax = axes[2]
    ax.scatter(
        df["mpnn_score"], df["el_rank"],
        c=C_D, s=60, edgecolors="black", linewidths=0.5, zorder=3,
    )

    # Compute Pareto-optimal front
    pareto = []
    sorted_by_mpnn = df.sort_values("mpnn_score")
    best_el = float("inf")
    for _, row in sorted_by_mpnn.iterrows():
        if row["el_rank"] < best_el:
            pareto.append(row)
            best_el = row["el_rank"]

    if pareto:
        pareto_df = pd.DataFrame(pareto)
        ax.scatter(
            pareto_df["mpnn_score"], pareto_df["el_rank"],
            c="red", s=100, marker="D", edgecolors="black", linewidths=0.8,
            zorder=5, label=f"Pareto front ({len(pareto)})",
        )
        ax.plot(
            pareto_df["mpnn_score"], pareto_df["el_rank"],
            color="red", linestyle="--", alpha=0.5, zorder=4,
        )
        for _, row in pareto_df.iterrows():
            ax.annotate(
                row["sequence"],
                (row["mpnn_score"], row["el_rank"]),
                fontsize=6.5, ha="left",
                xytext=(5, 5), textcoords="offset points",
            )

    ax.set_xlabel("MPNN Score (lower = better design)")
    ax.set_ylabel("EL %Rank (lower = stronger binder)")
    ax.set_title("Pareto Front: Design vs Binding", fontweight="bold")
    ax.legend(fontsize=9)

    plt.tight_layout()
    FIGURES_DIR.mkdir(exist_ok=True)
    out_path = FIGURES_DIR / f"{target}_decoy_d_detail.png"
    plt.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[OK] {target} -> {out_path}")
    return out_path


def load_all_targets() -> list[dict]:
    """Load targets from candidate_targets.json."""
    json_path = DATA_DIR / "candidate_targets.json"
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    targets = []
    # Existing targets (default to HLA-A*02:01)
    for seq in data.get("existing_targets", []):
        targets.append({"sequence": seq, "hla_allele": "HLA-A*02:01"})
    # Proposed targets
    for entry in data.get("proposed_targets", []):
        targets.append({
            "sequence": entry["sequence"],
            "hla_allele": entry.get("hla_allele", "HLA-A*02:01"),
        })
    return targets


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate Decoy D detail figures",
    )
    parser.add_argument(
        "targets", nargs="*",
        help="Target peptide sequence(s). Use --all for all candidates.",
    )
    parser.add_argument(
        "--all", action="store_true",
        help="Generate for all targets in candidate_targets.json",
    )
    args = parser.parse_args()

    if args.all:
        all_targets = load_all_targets()
        seqs = [t["sequence"] for t in all_targets]
    elif args.targets:
        seqs = [t.upper() for t in args.targets]
    else:
        parser.print_help()
        sys.exit(1)

    generated = 0
    for seq in seqs:
        result = plot_decoy_d_detail(seq)
        if result:
            generated += 1

    print(f"\nDone: {generated}/{len(seqs)} figures generated in {FIGURES_DIR}")


if __name__ == "__main__":
    main()
