"""
Presentation Score Re-scoring Visualization
============================================
Compare affinity-only vs presentation-based classification.
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "decoy_a"
OUT_DIR = PROJECT_ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)

PARQUET = DATA_DIR / "hla_filtered_HLA-A0201_presentation.parquet"


def main():
    print("Loading presentation-scored data...")
    df = pd.read_parquet(PARQUET)
    print(f"  {len(df):,} entries loaded")

    fig = plt.figure(figsize=(20, 14))
    fig.suptitle(
        "Presentation Score Model — Affinity vs Presentation Reclassification\n"
        f"HLA-A*02:01  |  N = {len(df):,}  |  mhcflurry 2.2.0 Class1PresentationPredictor",
        fontsize=16, fontweight="bold", y=0.98,
    )

    # Layout: 3 rows × 2 cols
    gs = fig.add_gridspec(3, 2, hspace=0.38, wspace=0.3,
                          left=0.07, right=0.95, top=0.92, bottom=0.05)

    # ── Panel A: Sankey-style reclassification ───────────────────────────
    ax = fig.add_subplot(gs[0, :])
    _plot_reclassification(ax, df)

    # ── Panel B: Processing score distribution ───────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    _plot_processing_dist(ax, df)

    # ── Panel C: Presentation score distribution ─────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    _plot_presentation_dist(ax, df)

    # ── Panel D: Affinity %Rank vs Presentation %Rank scatter ────────────
    ax = fig.add_subplot(gs[2, 0])
    _plot_rank_comparison(ax, df)

    # ── Panel E: Processing score by old binding class ───────────────────
    ax = fig.add_subplot(gs[2, 1])
    _plot_processing_by_class(ax, df)

    out = OUT_DIR / "decoy_a_presentation_rescore.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def _plot_reclassification(ax, df):
    """Panel A: Waterfall / alluvial showing old → new classification."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    ax.set_title("A. Reclassification: Affinity-Only → Presentation Model",
                 fontsize=13, fontweight="bold", loc="left")

    # Old counts
    old_sb = (df["binding"] == "Strong_Binder").sum()
    old_wb = (df["binding"] == "Weak_Binder").sum()
    # New counts
    new_sb = (df["presentation_binding"] == "Strong_Binder").sum()
    new_wb = (df["presentation_binding"] == "Weak_Binder").sum()
    new_nb = (df["presentation_binding"] == "Non_Binder").sum()

    # Cross-tab
    sb_to_sb = ((df["binding"] == "Strong_Binder") & (df["presentation_binding"] == "Strong_Binder")).sum()
    sb_to_wb = ((df["binding"] == "Strong_Binder") & (df["presentation_binding"] == "Weak_Binder")).sum()
    wb_to_sb = ((df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Strong_Binder")).sum()
    wb_to_wb = ((df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Weak_Binder")).sum()
    wb_to_nb = ((df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Non_Binder")).sum()

    total = len(df)

    # Left column: Old
    left_x = 0.5
    box_w = 2.8
    # Old Strong
    h_sb_old = old_sb / total * 3.2
    y_sb_old = 3.5 - h_sb_old
    ax.add_patch(mpatches.FancyBboxPatch(
        (left_x, y_sb_old), box_w, h_sb_old,
        boxstyle="round,pad=0.05", facecolor="#C44E52", alpha=0.85,
        edgecolor="white", linewidth=1.5))
    ax.text(left_x + box_w / 2, y_sb_old + h_sb_old / 2,
            f"Strong Binder\n{old_sb:,}\n({old_sb/total:.1%})",
            ha="center", va="center", fontsize=10, fontweight="bold", color="white")

    # Old Weak
    h_wb_old = old_wb / total * 3.2
    y_wb_old = y_sb_old - h_wb_old - 0.1
    ax.add_patch(mpatches.FancyBboxPatch(
        (left_x, y_wb_old), box_w, h_wb_old,
        boxstyle="round,pad=0.05", facecolor="#DD8452", alpha=0.85,
        edgecolor="white", linewidth=1.5))
    ax.text(left_x + box_w / 2, y_wb_old + h_wb_old / 2,
            f"Weak Binder\n{old_wb:,}\n({old_wb/total:.1%})",
            ha="center", va="center", fontsize=10, fontweight="bold", color="white")

    # Right column: New
    right_x = 6.5
    # New Strong
    h_sb_new = new_sb / total * 3.2
    y_sb_new = 3.5 - h_sb_new
    ax.add_patch(mpatches.FancyBboxPatch(
        (right_x, y_sb_new), box_w, h_sb_new,
        boxstyle="round,pad=0.05", facecolor="#C44E52", alpha=0.85,
        edgecolor="white", linewidth=1.5))
    ax.text(right_x + box_w / 2, y_sb_new + h_sb_new / 2,
            f"Strong Binder\n{new_sb:,}\n({new_sb/total:.1%})",
            ha="center", va="center", fontsize=10, fontweight="bold", color="white")

    # New Weak
    h_wb_new = new_wb / total * 3.2
    y_wb_new = y_sb_new - h_wb_new - 0.08
    ax.add_patch(mpatches.FancyBboxPatch(
        (right_x, y_wb_new), box_w, h_wb_new,
        boxstyle="round,pad=0.05", facecolor="#DD8452", alpha=0.85,
        edgecolor="white", linewidth=1.5))
    ax.text(right_x + box_w / 2, y_wb_new + h_wb_new / 2,
            f"Weak Binder\n{new_wb:,}\n({new_wb/total:.1%})",
            ha="center", va="center", fontsize=10, fontweight="bold", color="white")

    # New Non-Binder
    h_nb_new = new_nb / total * 3.2
    y_nb_new = y_wb_new - h_nb_new - 0.08
    ax.add_patch(mpatches.FancyBboxPatch(
        (right_x, y_nb_new), box_w, h_nb_new,
        boxstyle="round,pad=0.05", facecolor="#bdc3c7", alpha=0.85,
        edgecolor="white", linewidth=1.5))
    ax.text(right_x + box_w / 2, y_nb_new + h_nb_new / 2,
            f"Non-Binder\n{new_nb:,}\n({new_nb/total:.1%})",
            ha="center", va="center", fontsize=10, fontweight="bold", color="#555")

    # Labels
    ax.text(left_x + box_w / 2, 3.8, "Affinity Only", ha="center",
            fontsize=12, fontweight="bold", color="#2c3e50")
    ax.text(right_x + box_w / 2, 3.8, "Presentation Model", ha="center",
            fontsize=12, fontweight="bold", color="#2c3e50")

    # Flow arrows with counts
    mid_x = 5.0
    flows = [
        (y_sb_old + h_sb_old * 0.7, y_sb_new + h_sb_new * 0.7, sb_to_sb, "kept", "#C44E52"),
        (y_sb_old + h_sb_old * 0.3, y_wb_new + h_wb_new * 0.8, sb_to_wb, "46K downgraded", "#e67e22"),
        (y_wb_old + h_wb_old * 0.8, y_sb_new + h_sb_new * 0.2, wb_to_sb, "92K upgraded", "#27ae60"),
        (y_wb_old + h_wb_old * 0.5, y_wb_new + h_wb_new * 0.5, wb_to_wb, "kept", "#DD8452"),
        (y_wb_old + h_wb_old * 0.15, y_nb_new + h_nb_new * 0.5, wb_to_nb, "308K removed", "#95a5a6"),
    ]

    for y_from, y_to, count, label, color in flows:
        if count < 10000:
            continue
        lw = max(0.8, min(4, count / 100000))
        ax.annotate("", xy=(right_x - 0.1, y_to), xytext=(left_x + box_w + 0.1, y_from),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, alpha=0.5,
                                    connectionstyle="arc3,rad=0.1"))
        if label != "kept":
            ax.text(mid_x, (y_from + y_to) / 2, label,
                    ha="center", va="center", fontsize=9, color=color,
                    fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=color, alpha=0.9))


def _plot_processing_dist(ax, df):
    """Panel B: Processing score histogram by old binding class."""
    for label, color in [("Strong_Binder", "#C44E52"), ("Weak_Binder", "#DD8452")]:
        vals = df.loc[df["binding"] == label, "processing_score"]
        ax.hist(vals, bins=80, alpha=0.6, color=color, label=label.replace("_", " "),
                edgecolor="white", linewidth=0.3)

    ax.axvline(df["processing_score"].median(), color="#333", ls="--", lw=1,
               label=f'Median = {df["processing_score"].median():.3f}')
    ax.set_xlabel("Processing Score", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("B. Antigen Processing Score Distribution", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f}k"))


def _plot_presentation_dist(ax, df):
    """Panel C: Presentation score histogram colored by new classification."""
    for label, color in [("Strong_Binder", "#C44E52"), ("Weak_Binder", "#DD8452"),
                          ("Non_Binder", "#bdc3c7")]:
        vals = df.loc[df["presentation_binding"] == label, "presentation_score"]
        ax.hist(vals, bins=80, alpha=0.6, color=color,
                label=f'{label.replace("_", " ")} ({len(vals):,})',
                edgecolor="white", linewidth=0.3)

    ax.set_xlabel("Presentation Score", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("C. Presentation Score Distribution (New Classification)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f}k"))


def _plot_rank_comparison(ax, df):
    """Panel D: Affinity %Rank vs Presentation %Rank density."""
    # Subsample for plotting speed
    rng = np.random.default_rng(42)
    idx = rng.choice(len(df), size=min(50000, len(df)), replace=False)
    sub = df.iloc[idx]

    color_map = {"Strong_Binder": "#C44E52", "Weak_Binder": "#DD8452", "Non_Binder": "#bdc3c7"}
    for label in ["Non_Binder", "Weak_Binder", "Strong_Binder"]:
        mask = sub["presentation_binding"] == label
        ax.scatter(sub.loc[mask, "el_rank"], sub.loc[mask, "presentation_percentile"],
                   c=color_map[label], s=3, alpha=0.3, label=label.replace("_", " "),
                   rasterized=True)

    # Reference lines
    ax.axhline(2.0, color="#888", ls=":", lw=0.8, alpha=0.5)
    ax.axvline(2.0, color="#888", ls=":", lw=0.8, alpha=0.5)
    ax.axhline(0.5, color="#C44E52", ls=":", lw=0.8, alpha=0.5)
    ax.axvline(0.5, color="#C44E52", ls=":", lw=0.8, alpha=0.5)
    ax.plot([0, 4.5], [0, 4.5], color="#aaa", ls="--", lw=0.8, alpha=0.5)

    ax.set_xlabel("Affinity %Rank (old)", fontsize=11)
    ax.set_ylabel("Presentation %Rank (new)", fontsize=11)
    ax.set_title("D. Affinity vs Presentation %Rank (50K subsample)",
                 fontsize=13, fontweight="bold")
    ax.legend(fontsize=9, markerscale=4)
    ax.set_xlim(-0.05, 2.15)
    ax.set_ylim(-0.05, 4.5)


def _plot_processing_by_class(ax, df):
    """Panel E: Processing score violin by reclassification outcome."""
    # Define groups
    groups = {
        "SB → SB\n(kept)": df[(df["binding"] == "Strong_Binder") & (df["presentation_binding"] == "Strong_Binder")],
        "SB → WB\n(downgraded)": df[(df["binding"] == "Strong_Binder") & (df["presentation_binding"] == "Weak_Binder")],
        "WB → SB\n(upgraded)": df[(df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Strong_Binder")],
        "WB → WB\n(kept)": df[(df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Weak_Binder")],
        "WB → NB\n(removed)": df[(df["binding"] == "Weak_Binder") & (df["presentation_binding"] == "Non_Binder")],
    }
    colors = ["#C44E52", "#e67e22", "#27ae60", "#DD8452", "#bdc3c7"]

    # Subsample each group for violin
    rng = np.random.default_rng(42)
    data = []
    for g in groups.values():
        n = min(20000, len(g))
        data.append(g["processing_score"].values[rng.choice(len(g), n, replace=False)])

    positions = range(1, len(groups) + 1)
    vp = ax.violinplot(data, positions=positions, showmedians=True, showextrema=False, widths=0.7)
    for i, body in enumerate(vp["bodies"]):
        body.set_facecolor(colors[i])
        body.set_alpha(0.7)
    vp["cmedians"].set_color("#333")
    vp["cmedians"].set_linewidth(1.5)

    ax.set_xticks(list(positions))
    ax.set_xticklabels(list(groups.keys()), fontsize=9)
    ax.set_ylabel("Processing Score", fontsize=11)
    ax.set_title("E. Processing Score by Reclassification Outcome",
                 fontsize=13, fontweight="bold")

    # Add median labels
    for i, (name, g) in enumerate(groups.items()):
        med = g["processing_score"].median()
        ax.text(i + 1, med + 0.03, f"{med:.3f}", ha="center", fontsize=8, color="#333")


if __name__ == "__main__":
    main()
