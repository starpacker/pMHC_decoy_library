"""
Visualize Interface Descriptor Results for Decoy B

Generates figures showing how the new interface descriptors (PLIP, FreeSASA, PRODIGY)
affect the scoring and ranking of Decoy B candidates.

Usage (on server, after running full pipeline):
    cd /share/liuyutian/pMHC_decoy_library
    python scripts/visualize_interface_descriptors.py

    # Or with a specific results file:
    python scripts/visualize_interface_descriptors.py --results data/decoy_b/final_ranked_decoys.json

Output:
    figures/interface_descriptor_analysis.png
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

OUT = PROJECT_ROOT / "figures"
OUT.mkdir(exist_ok=True)

# ── Color palette ───────────────────────────────────────────────────────
C_PLIP = "#E91E63"
C_BSA = "#2196F3"
C_PROD = "#FF9800"
C_COMB = "#4CAF50"
C_RMSD = "#9C27B0"
C_BG = "#F5F5F5"


def load_results(path: Path) -> List[dict]:
    """Load ranked decoy results."""
    data = json.loads(path.read_text(encoding="utf-8"))
    return data


def plot_interface_analysis(entries: List[dict], target: str = "GILGFVFTL"):
    """Generate comprehensive interface descriptor visualization."""
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.facecolor": C_BG,
        "figure.facecolor": "white",
        "axes.grid": True,
        "grid.alpha": 0.3,
    })

    # Filter entries with interface data
    has_interface = [e for e in entries if e.get("interface_combined") is not None]
    has_struct = [e for e in entries if e.get("structural") and e["structural"].get("rmsd") is not None]

    fig = plt.figure(figsize=(22, 16))
    fig.suptitle(
        f"Interface Descriptor Analysis — {target} Decoys (n={len(entries)})",
        fontsize=16, fontweight="bold", y=0.98,
    )
    gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.35,
                           top=0.93, bottom=0.06, left=0.06, right=0.96)

    # ── 1a: Overview — descriptor availability ─────────────────────────────
    ax1 = fig.add_subplot(gs[0, 0])
    n_plip = sum(1 for e in entries if e.get("plip_tanimoto") is not None)
    n_bsa = sum(1 for e in entries if e.get("bsa_similarity") is not None)
    n_prod = sum(1 for e in entries if e.get("prodigy_similarity") is not None)
    n_comb = sum(1 for e in entries if e.get("interface_combined") is not None)
    n_total = len(entries)

    labels = ["PLIP", "FreeSASA", "PRODIGY", "Combined"]
    vals = [n_plip, n_bsa, n_prod, n_comb]
    colors = [C_PLIP, C_BSA, C_PROD, C_COMB]
    bars = ax1.bar(labels, vals, color=colors, edgecolor="white", linewidth=1.5)
    ax1.axhline(n_total, color="gray", linestyle="--", alpha=0.5, label=f"Total entries: {n_total}")
    for bar, v in zip(bars, vals):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                 f"{v}/{n_total}", ha="center", va="bottom", fontweight="bold", fontsize=10)
    ax1.set_ylabel("Entries with Score")
    ax1.set_title("Descriptor Availability", fontweight="bold")
    ax1.legend(fontsize=9)

    # ── 1b: Interface combined vs structural_similarity scatter ────────────
    ax2 = fig.add_subplot(gs[0, 1])
    if has_interface:
        x = [e["structural_similarity"] for e in has_interface]
        y = [e["interface_combined"] for e in has_interface]
        risk = [e["total_risk_score"] for e in has_interface]
        sc = ax2.scatter(x, y, c=risk, cmap="YlOrRd", s=60, edgecolors="black",
                        linewidths=0.5, zorder=3)
        plt.colorbar(sc, ax=ax2, label="Risk Score")
        ax2.set_xlabel("Structural Similarity (surface_correlation)")
        ax2.set_ylabel("Interface Combined Score")
        ax2.set_title("Structural vs Interface Scores", fontweight="bold")
        # Diagonal
        ax2.plot([0, 1], [0, 1], "k--", alpha=0.3)
    else:
        ax2.text(0.5, 0.5, "No interface data", ha="center", va="center", transform=ax2.transAxes)
        ax2.set_title("Structural vs Interface Scores", fontweight="bold")

    # ── 1c: Per-descriptor distribution (violin or box) ───────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    plot_data = []
    plot_labels = []
    plot_colors = []
    if n_plip > 0:
        plot_data.append([e["plip_tanimoto"] for e in entries if e.get("plip_tanimoto") is not None])
        plot_labels.append("PLIP")
        plot_colors.append(C_PLIP)
    if n_bsa > 0:
        plot_data.append([e["bsa_similarity"] for e in entries if e.get("bsa_similarity") is not None])
        plot_labels.append("BSA")
        plot_colors.append(C_BSA)
    if n_prod > 0:
        plot_data.append([e["prodigy_similarity"] for e in entries if e.get("prodigy_similarity") is not None])
        plot_labels.append("PRODIGY")
        plot_colors.append(C_PROD)
    if n_comb > 0:
        plot_data.append([e["interface_combined"] for e in entries if e.get("interface_combined") is not None])
        plot_labels.append("Combined")
        plot_colors.append(C_COMB)

    if plot_data:
        parts = ax3.violinplot(plot_data, positions=range(len(plot_data)),
                               showmedians=True, showextrema=True)
        for i, pc in enumerate(parts["bodies"]):
            pc.set_facecolor(plot_colors[i])
            pc.set_alpha(0.7)
        parts["cmedians"].set_color("black")
        ax3.set_xticks(range(len(plot_labels)))
        ax3.set_xticklabels(plot_labels)
        ax3.set_ylabel("Score [0-1]")
        ax3.set_title("Descriptor Score Distributions", fontweight="bold")
    else:
        ax3.text(0.5, 0.5, "No descriptor data", ha="center", va="center", transform=ax3.transAxes)
        ax3.set_title("Descriptor Score Distributions", fontweight="bold")

    # ── 2a: Top 20 radar chart ────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0], polar=True)
    top20_interface = [e for e in has_interface if e.get("plip_tanimoto") is not None][:20]
    if top20_interface:
        categories = ["PLIP", "BSA", "PRODIGY", "Structural\nSim", "PhysChem\nSim"]
        N = len(categories)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]

        for i, e in enumerate(top20_interface[:5]):
            vals = [
                e.get("plip_tanimoto", 0) or 0,
                e.get("bsa_similarity", 0) or 0,
                e.get("prodigy_similarity", 0) or 0,
                e.get("structural_similarity", 0) or 0,
                e.get("physicochemical_similarity", 0) or 0,
            ]
            vals += vals[:1]
            alpha = 0.6 - i * 0.08
            ax4.plot(angles, vals, 'o-', linewidth=1.5, alpha=max(0.3, alpha),
                    label=f"#{i+1} {e['sequence']}")
            ax4.fill(angles, vals, alpha=0.05)

        ax4.set_xticks(angles[:-1])
        ax4.set_xticklabels(categories, fontsize=8)
        ax4.set_ylim(0, 1)
        ax4.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.3, 1.1))
        ax4.set_title("Top-5 Descriptor Profiles", fontweight="bold", pad=20)
    else:
        ax4.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax4.transAxes)

    # ── 2b: PLIP vs BSA scatter ───────────────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    both_pb = [e for e in entries
               if e.get("plip_tanimoto") is not None and e.get("bsa_similarity") is not None]
    if both_pb:
        x = [e["plip_tanimoto"] for e in both_pb]
        y = [e["bsa_similarity"] for e in both_pb]
        risk = [e["total_risk_score"] for e in both_pb]
        sc = ax5.scatter(x, y, c=risk, cmap="YlOrRd", s=60,
                        edgecolors="black", linewidths=0.5, zorder=3)
        plt.colorbar(sc, ax=ax5, label="Risk Score")
        ax5.set_xlabel("PLIP Tanimoto")
        ax5.set_ylabel("BSA Similarity")
        ax5.set_title("PLIP vs BSA (colored by risk)", fontweight="bold")
    else:
        ax5.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax5.transAxes)
        ax5.set_title("PLIP vs BSA", fontweight="bold")

    # ── 2c: Risk score: old formula vs new contribution ───────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    if has_interface:
        # Show how interface_combined shifts the ranking
        ranks_by_struct = sorted(has_interface, key=lambda e: -e.get("structural_similarity", 0))
        ranks_by_risk = sorted(has_interface, key=lambda e: -e["total_risk_score"])

        struct_rank = {e["sequence"]: i for i, e in enumerate(ranks_by_struct)}
        risk_rank = {e["sequence"]: i for i, e in enumerate(ranks_by_risk)}

        seqs = list(struct_rank.keys())
        x_vals = [struct_rank[s] for s in seqs]
        y_vals = [risk_rank[s] for s in seqs]
        ax6.scatter(x_vals, y_vals, c=C_COMB, s=40, alpha=0.7, edgecolors="black", linewidths=0.3)
        ax6.plot([0, len(seqs)], [0, len(seqs)], "k--", alpha=0.3, label="No change")
        ax6.set_xlabel("Rank by Structural Similarity")
        ax6.set_ylabel("Rank by Risk Score")
        ax6.set_title("Rank Shift: Structural vs Final", fontweight="bold")
        ax6.legend(fontsize=9)
    else:
        ax6.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax6.transAxes)
        ax6.set_title("Rank Shift Analysis", fontweight="bold")

    # ── 3a: Top 20 horizontal bar with descriptor breakdown ───────────────
    ax7 = fig.add_subplot(gs[2, :2])
    top20 = entries[:20]
    if top20:
        y_pos = np.arange(len(top20))
        seq_labels = [f"#{i+1} {e['sequence']}" for i, e in enumerate(top20)]

        # Stacked bar: physicochemical + structural + interface
        phys = [e.get("physicochemical_similarity", 0) or 0 for e in top20]
        struct = [e.get("structural_similarity", 0) or 0 for e in top20]
        iface = [e.get("interface_combined", 0) or 0 for e in top20]

        ax7.barh(y_pos, phys, color="#90CAF9", label="PhysChem", edgecolor="white", height=0.7)
        ax7.barh(y_pos, struct, left=phys, color=C_RMSD, label="Structural", edgecolor="white", height=0.7)
        ax7.barh(y_pos, iface, left=[p+s for p, s in zip(phys, struct)],
                color=C_COMB, label="Interface", edgecolor="white", height=0.7)

        ax7.set_yticks(y_pos)
        ax7.set_yticklabels(seq_labels, fontfamily="monospace", fontsize=8)
        ax7.invert_yaxis()
        ax7.set_xlabel("Score Components")
        ax7.set_title("Top 20: Score Breakdown", fontweight="bold")
        ax7.legend(fontsize=9, loc="lower right")

    # ── 3b: Summary statistics ────────────────────────────────────────────
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.axis("off")
    summary_text = [
        f"Total entries: {len(entries)}",
        f"With interface data: {len(has_interface)}",
        f"With structural data: {len(has_struct)}",
        "",
        "Interface descriptor coverage:",
        f"  PLIP:    {n_plip}/{n_total} ({100*n_plip/max(n_total,1):.0f}%)",
        f"  BSA:     {n_bsa}/{n_total} ({100*n_bsa/max(n_total,1):.0f}%)",
        f"  PRODIGY: {n_prod}/{n_total} ({100*n_prod/max(n_total,1):.0f}%)",
        "",
    ]

    if has_interface:
        ic_vals = [e["interface_combined"] for e in has_interface]
        summary_text.extend([
            "Interface combined score:",
            f"  Mean: {np.mean(ic_vals):.3f}",
            f"  Std:  {np.std(ic_vals):.3f}",
            f"  Range: [{min(ic_vals):.3f}, {max(ic_vals):.3f}]",
        ])

    ax8.text(0.05, 0.95, "\n".join(summary_text),
             transform=ax8.transAxes, fontsize=10, fontfamily="monospace",
             verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))
    ax8.set_title("Summary", fontweight="bold")

    plt.savefig(OUT / "interface_descriptor_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {OUT / 'interface_descriptor_analysis.png'}")


def main():
    parser = argparse.ArgumentParser(description="Visualize interface descriptor results")
    parser.add_argument("--results", type=Path, default=None,
                       help="Path to final_ranked_decoys.json")
    args = parser.parse_args()

    # Try to find results file
    candidates = [
        args.results,
        PROJECT_ROOT / "data" / "decoy_b" / "final_ranked_decoys.json",
        Path("/share/liuyutian/pMHC_decoy_library/data/decoy_b/final_ranked_decoys.json"),
    ]

    results_path = None
    for c in candidates:
        if c and c.exists():
            results_path = c
            break

    if results_path is None:
        print("No results file found. Run the full Decoy B pipeline first.")
        print("Expected: data/decoy_b/final_ranked_decoys.json")
        sys.exit(1)

    entries = load_results(results_path)
    print(f"Loaded {len(entries)} entries from {results_path}")

    plot_interface_analysis(entries)


if __name__ == "__main__":
    main()
