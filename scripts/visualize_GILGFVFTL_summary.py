"""
Visualization of GILGFVFTL Decoy Library Results (Decoy A / B / D)
Target peptide: GILGFVFTL (Influenza M1, HLA-A*02:01)
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from collections import Counter

# ── paths ────────────────────────────────────────────────────────────────
REPO = Path(r"C:\Users\30670\AppData\Local\Temp\pMHC_decoy_library")
DATA = REPO / "data" / "GILGFVFTL_summary"
OUT  = Path(r"C:\Users\30670\Desktop\decoy_library\figures")
OUT.mkdir(exist_ok=True)

# ── load data ────────────────────────────────────────────────────────────
with open(DATA / "Decoy_A" / "decoy_a_results.json") as f:
    decoy_a = json.load(f)

with open(DATA / "Decoy_B" / "final_ranked_decoys.json") as f:
    decoy_b = json.load(f)

decoy_d = pd.read_csv(DATA / "Decoy_D" / "decoy_d_results.csv")
# Remove the target itself from Decoy D
decoy_d = decoy_d[decoy_d["sequence"] != "GILGFVFTL"].copy()

TARGET = "GILGFVFTL"

# ── color palette ────────────────────────────────────────────────────────
C_A = "#2196F3"   # blue
C_B = "#FF9800"   # orange
C_D = "#4CAF50"   # green
C_BG = "#F5F5F5"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.facecolor": C_BG,
    "figure.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
})


# ======================================================================
# Figure 1: Overview Dashboard
# ======================================================================
fig = plt.figure(figsize=(20, 16))
fig.suptitle(
    f"Decoy Library Summary for {TARGET} (HLA-A*02:01)",
    fontsize=18, fontweight="bold", y=0.98,
)
gs = gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.35,
                       top=0.93, bottom=0.06, left=0.07, right=0.95)

# ── 1a: Pipeline overview (bar chart of counts) ─────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
counts = [len(decoy_a), len(decoy_b), len(decoy_d)]
labels = ["Decoy A\n(Sequence\nHomology)", "Decoy B\n(Structural\nSimilarity)", "Decoy D\n(MPNN\nInverse Design)"]
colors = [C_A, C_B, C_D]
bars = ax1.bar(labels, counts, color=colors, edgecolor="white", linewidth=1.5, width=0.6)
for bar, c in zip(bars, counts):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
             str(c), ha="center", va="bottom", fontweight="bold", fontsize=13)
ax1.set_ylabel("Number of Decoy Candidates")
ax1.set_title("Pipeline Output Counts", fontweight="bold")
ax1.set_ylim(0, max(counts) * 1.15)

# ── 1b: Decoy A - Hamming Distance Distribution ─────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
hd_a = [d["hamming_distance"] for d in decoy_a]
hd_counts = Counter(hd_a)
hd_x = sorted(hd_counts.keys())
hd_y = [hd_counts[x] for x in hd_x]
bars2 = ax2.bar(hd_x, hd_y, color=C_A, edgecolor="white", linewidth=1.2)
for bar, c in zip(bars2, hd_y):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
             str(c), ha="center", va="bottom", fontsize=10)
ax2.set_xlabel("Hamming Distance to Target")
ax2.set_ylabel("Count")
ax2.set_title("Decoy A: Hamming Distance", fontweight="bold")
ax2.set_xticks(hd_x)

# ── 1c: Decoy A - TCR Contact vs Anchor Mismatches ──────────────────────
ax3 = fig.add_subplot(gs[0, 2])
tcr_mis = [d["n_tcr_contact_mismatches"] for d in decoy_a]
anc_mis = [d["n_anchor_mismatches"] for d in decoy_a]
# 2D histogram as heatmap
max_tcr, max_anc = max(tcr_mis), max(anc_mis)
heatmap = np.zeros((max_anc + 1, max_tcr + 1))
for t, a in zip(tcr_mis, anc_mis):
    heatmap[a, t] += 1
im = ax3.imshow(heatmap, cmap="Blues", aspect="auto", origin="lower")
for i in range(heatmap.shape[0]):
    for j in range(heatmap.shape[1]):
        if heatmap[i, j] > 0:
            ax3.text(j, i, f"{int(heatmap[i,j])}", ha="center", va="center",
                     fontsize=9, color="white" if heatmap[i,j] > heatmap.max()*0.5 else "black")
ax3.set_xlabel("TCR Contact Mismatches")
ax3.set_ylabel("Anchor Mismatches")
ax3.set_title("Decoy A: Mismatch Profile", fontweight="bold")
ax3.set_xticks(range(max_tcr + 1))
ax3.set_yticks(range(max_anc + 1))
plt.colorbar(im, ax=ax3, label="Count", shrink=0.8)

# ── 1d: Decoy A - EL Rank Distribution ──────────────────────────────────
ax4 = fig.add_subplot(gs[1, 0])
el_ranks_a = [d["el_rank"] for d in decoy_a]
ax4.hist(el_ranks_a, bins=40, color=C_A, edgecolor="white", alpha=0.85)
ax4.axvline(0.5, color="red", linestyle="--", linewidth=1.5, label="SB/WB cutoff (0.5)")
ax4.axvline(2.0, color="darkred", linestyle="--", linewidth=1.5, label="WB/NB cutoff (2.0)")
ax4.set_xlabel("EL %Rank")
ax4.set_ylabel("Count")
ax4.set_title("Decoy A: MHC Binding Strength", fontweight="bold")
ax4.legend(fontsize=8)

# ── 1e: Decoy B - Physicochemical vs Structural Similarity ──────────────
ax5 = fig.add_subplot(gs[1, 1])
b_phys = [d["physicochemical_similarity"] for d in decoy_b]
b_struct = [d["structural_similarity"] for d in decoy_b]
b_source = [d["source"] for d in decoy_b]
colors_b = [C_B if "Structural" in s else C_A for s in b_source]
ax5.scatter(b_phys, b_struct, c=colors_b, s=60, edgecolors="white", linewidths=0.5, alpha=0.85, zorder=3)
ax5.set_xlabel("Physicochemical Similarity")
ax5.set_ylabel("Structural Similarity (Surface Correlation)")
ax5.set_title("Decoy B: Phys. vs Structural", fontweight="bold")
# Add legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=C_B, markersize=10, label=f'Structural ({sum(1 for s in b_source if "Structural" in s)})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=C_A, markersize=10, label=f'Sequence ({sum(1 for s in b_source if "Sequence" in s)})')
]
ax5.legend(handles=legend_elements, fontsize=8, loc="lower left")

# ── 1f: Decoy B - Risk Score Distribution ────────────────────────────────
ax6 = fig.add_subplot(gs[1, 2])
risk_scores = [d["total_risk_score"] for d in decoy_b]
colors_risk = [C_B if "Structural" in d["source"] else C_A for d in decoy_b]
ax6.barh(range(len(risk_scores)), risk_scores, color=colors_risk, edgecolor="white", height=0.8)
ax6.set_xlabel("Total Risk Score")
ax6.set_ylabel("Rank")
ax6.set_title("Decoy B: Risk Score Ranking (Top 50)", fontweight="bold")
ax6.set_yticks(range(0, len(risk_scores), 5))
ax6.set_yticklabels([str(i+1) for i in range(0, len(risk_scores), 5)])
ax6.invert_yaxis()

# ── 1g: Decoy D - MPNN Score vs EL Rank ─────────────────────────────────
ax7 = fig.add_subplot(gs[2, 0])
ax7.scatter(decoy_d["mpnn_score"], decoy_d["el_rank"],
            c=C_D, s=50, edgecolors="white", linewidths=0.5, alpha=0.8, zorder=3)
ax7.set_xlabel("MPNN Score (lower = better)")
ax7.set_ylabel("EL %Rank (lower = stronger binder)")
ax7.set_title("Decoy D: MPNN Score vs MHC Binding", fontweight="bold")
# Label top candidates
for _, row in decoy_d.head(5).iterrows():
    ax7.annotate(row["sequence"], (row["mpnn_score"], row["el_rank"]),
                 fontsize=6.5, ha="left", va="bottom",
                 xytext=(4, 4), textcoords="offset points")

# ── 1h: Decoy D - MPNN Score Distribution ───────────────────────────────
ax8 = fig.add_subplot(gs[2, 1])
ax8.hist(decoy_d["mpnn_score"], bins=20, color=C_D, edgecolor="white", alpha=0.85)
ax8.set_xlabel("MPNN Score")
ax8.set_ylabel("Count")
ax8.set_title("Decoy D: MPNN Score Distribution", fontweight="bold")
ax8.axvline(decoy_d["mpnn_score"].median(), color="darkgreen", linestyle="--",
            label=f'Median: {decoy_d["mpnn_score"].median():.2f}')
ax8.legend(fontsize=9)

# ── 1i: Cross-pipeline comparison ───────────────────────────────────────
ax9 = fig.add_subplot(gs[2, 2])
# Compare EL rank distributions across pipelines
el_a = [d["el_rank"] for d in decoy_a]
el_b = [d["el_rank"] for d in decoy_b]
el_d = decoy_d["el_rank"].tolist()
parts = ax9.violinplot([el_a, el_b, el_d], positions=[1, 2, 3], showmedians=True, showextrema=True)
for i, pc in enumerate(parts["bodies"]):
    pc.set_facecolor([C_A, C_B, C_D][i])
    pc.set_alpha(0.7)
parts["cmedians"].set_color("black")
parts["cmins"].set_color("gray")
parts["cmaxes"].set_color("gray")
parts["cbars"].set_color("gray")
ax9.set_xticks([1, 2, 3])
ax9.set_xticklabels(["Decoy A\n(n=591)", "Decoy B\n(n=50)", "Decoy D\n(n=42)"])
ax9.set_ylabel("EL %Rank")
ax9.set_title("MHC Binding Comparison", fontweight="bold")

plt.savefig(OUT / "GILGFVFTL_dashboard.png", dpi=180, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT / 'GILGFVFTL_dashboard.png'}")


# ======================================================================
# Figure 2: Sequence Logo / Position-wise Analysis
# ======================================================================
fig2, axes2 = plt.subplots(2, 2, figsize=(18, 12))
fig2.suptitle(f"Position-wise Analysis for {TARGET} Decoys", fontsize=16, fontweight="bold")

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

# ── 2a: Decoy A - Position-wise AA frequency (top 20 by risk) ───────────
ax = axes2[0, 0]
# Sort Decoy A by el_rank (best binders first), take top 50
sorted_a = sorted(decoy_a, key=lambda d: d["el_rank"])[:50]
seqs_a = [d["sequence"] for d in sorted_a]
freq_a = np.zeros((len(AA_ORDER), len(TARGET)))
for seq in seqs_a:
    for j, aa in enumerate(seq):
        if aa in AA_ORDER:
            freq_a[AA_ORDER.index(aa), j] += 1
freq_a /= max(len(seqs_a), 1)
im2a = ax.imshow(freq_a, cmap="Blues", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(TARGET)))
ax.set_xticklabels([f"p{i+1}\n({TARGET[i]})" for i in range(len(TARGET))])
ax.set_yticks(range(len(AA_ORDER)))
ax.set_yticklabels(list(AA_ORDER), fontsize=8)
ax.set_title("Decoy A Top-50 Binders: AA Frequency", fontweight="bold")
ax.set_xlabel("Position (target residue)")
ax.set_ylabel("Amino Acid")
plt.colorbar(im2a, ax=ax, label="Frequency", shrink=0.8)
# Mark target residues
for j, aa in enumerate(TARGET):
    idx = AA_ORDER.index(aa)
    ax.plot(j, idx, marker="*", color="red", markersize=10, zorder=5)

# ── 2b: Decoy D - Position-wise AA frequency ────────────────────────────
ax = axes2[0, 1]
seqs_d = decoy_d["sequence"].tolist()
freq_d = np.zeros((len(AA_ORDER), len(TARGET)))
for seq in seqs_d:
    for j, aa in enumerate(seq):
        if j < len(TARGET) and aa in AA_ORDER:
            freq_d[AA_ORDER.index(aa), j] += 1
freq_d /= max(len(seqs_d), 1)
im2b = ax.imshow(freq_d, cmap="Greens", aspect="auto", vmin=0, vmax=1)
ax.set_xticks(range(len(TARGET)))
ax.set_xticklabels([f"p{i+1}\n({TARGET[i]})" for i in range(len(TARGET))])
ax.set_yticks(range(len(AA_ORDER)))
ax.set_yticklabels(list(AA_ORDER), fontsize=8)
ax.set_title("Decoy D (MPNN): AA Frequency", fontweight="bold")
ax.set_xlabel("Position (target residue)")
ax.set_ylabel("Amino Acid")
plt.colorbar(im2b, ax=ax, label="Frequency", shrink=0.8)
for j, aa in enumerate(TARGET):
    idx = AA_ORDER.index(aa)
    ax.plot(j, idx, marker="*", color="red", markersize=10, zorder=5)

# ── 2c: Decoy B - Structural metrics scatter ────────────────────────────
ax = axes2[1, 0]
rmsds = [d["structural"]["rmsd"] for d in decoy_b if d.get("structural") and d["structural"].get("rmsd")]
surf_corrs = [d["structural"]["surface_correlation"] for d in decoy_b if d.get("structural") and d["structural"].get("surface_correlation")]
el_b2 = [d["el_rank"] for d in decoy_b if d.get("structural") and d["structural"].get("rmsd")]
if rmsds and surf_corrs:
    sc = ax.scatter(rmsds, surf_corrs, c=el_b2, cmap="YlOrRd_r", s=70,
                    edgecolors="black", linewidths=0.5, zorder=3)
    plt.colorbar(sc, ax=ax, label="EL %Rank")
    ax.set_xlabel("Peptide RMSD (A)")
    ax.set_ylabel("Surface Correlation")
    ax.set_title("Decoy B: Structure Quality vs Similarity", fontweight="bold")
    # Label a few
    b_with_struct = [d for d in decoy_b if d.get("structural") and d["structural"].get("rmsd")]
    for d in b_with_struct[:5]:
        ax.annotate(d["sequence"], (d["structural"]["rmsd"], d["structural"]["surface_correlation"]),
                    fontsize=6, ha="left", xytext=(3, 3), textcoords="offset points")

# ── 2d: Cross-pipeline sequence diversity ────────────────────────────────
ax = axes2[1, 1]
def avg_pairwise_identity(seqs, max_pairs=500):
    """Average pairwise identity for a set of same-length sequences."""
    if len(seqs) < 2:
        return 0
    n = len(seqs[0])
    ids = []
    np.random.seed(42)
    indices = np.random.choice(len(seqs), size=min(len(seqs), max_pairs), replace=False)
    for i in range(len(indices)):
        for j in range(i+1, min(i+10, len(indices))):
            s1, s2 = seqs[indices[i]], seqs[indices[j]]
            if len(s1) == len(s2):
                ids.append(sum(a == b for a, b in zip(s1, s2)) / len(s1))
    return np.mean(ids) if ids else 0

# Position-wise conservation
pos_cons_a = []
pos_cons_d = []
for p in range(len(TARGET)):
    aa_a = [s["sequence"][p] for s in sorted_a if len(s["sequence"]) == len(TARGET)]
    aa_d = [s[p] for s in seqs_d if len(s) == len(TARGET)]
    if aa_a:
        most_common_a = Counter(aa_a).most_common(1)[0][1] / len(aa_a)
        pos_cons_a.append(most_common_a)
    else:
        pos_cons_a.append(0)
    if aa_d:
        most_common_d = Counter(aa_d).most_common(1)[0][1] / len(aa_d)
        pos_cons_d.append(most_common_d)
    else:
        pos_cons_d.append(0)

x = np.arange(len(TARGET))
width = 0.35
bars_a = ax.bar(x - width/2, pos_cons_a, width, color=C_A, label="Decoy A (Top 50)", edgecolor="white")
bars_d = ax.bar(x + width/2, pos_cons_d, width, color=C_D, label="Decoy D (MPNN)", edgecolor="white")
ax.set_xticks(x)
ax.set_xticklabels([f"p{i+1}\n{TARGET[i]}" for i in range(len(TARGET))])
ax.set_ylabel("Most Frequent AA Fraction")
ax.set_title("Position-wise Conservation", fontweight="bold")
ax.legend(fontsize=9)
ax.set_ylim(0, 1.1)
# Mark anchor positions
for pos in [1, 8]:  # 0-indexed p2, p9
    ax.axvspan(pos - 0.45, pos + 0.45, alpha=0.1, color="red")
ax.text(1, 1.03, "anchor", ha="center", fontsize=7, color="red")
ax.text(8, 1.03, "anchor", ha="center", fontsize=7, color="red")

plt.tight_layout()
plt.savefig(OUT / "GILGFVFTL_position_analysis.png", dpi=180, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT / 'GILGFVFTL_position_analysis.png'}")


# ======================================================================
# Figure 3: Decoy D Detail
# ======================================================================
fig3, axes3 = plt.subplots(1, 3, figsize=(18, 6))
fig3.suptitle("Decoy D (MPNN Inverse Design) Detailed Analysis", fontsize=15, fontweight="bold")

# ── 3a: Top 15 candidates table-like horizontal bar ─────────────────────
ax = axes3[0]
top15 = decoy_d.nsmallest(15, "mpnn_score")
y_pos = range(len(top15))
ax.barh(y_pos, top15["mpnn_score"], color=C_D, edgecolor="white", height=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(top15["sequence"], fontfamily="monospace", fontsize=9)
ax.set_xlabel("MPNN Score")
ax.set_title("Top 15 by MPNN Score", fontweight="bold")
ax.invert_yaxis()
for i, (_, row) in enumerate(top15.iterrows()):
    ax.text(row["mpnn_score"] + 0.01, i, f'EL={row["el_rank"]:.2f}', va="center", fontsize=7)

# ── 3b: EL Rank sorted ──────────────────────────────────────────────────
ax = axes3[1]
top15_el = decoy_d.nsmallest(15, "el_rank")
y_pos2 = range(len(top15_el))
ax.barh(y_pos2, top15_el["el_rank"], color="#66BB6A", edgecolor="white", height=0.7)
ax.set_yticks(y_pos2)
ax.set_yticklabels(top15_el["sequence"], fontfamily="monospace", fontsize=9)
ax.set_xlabel("EL %Rank")
ax.set_title("Top 15 by MHC Binding", fontweight="bold")
ax.invert_yaxis()
for i, (_, row) in enumerate(top15_el.iterrows()):
    ax.text(row["el_rank"] + 0.005, i, f'MPNN={row["mpnn_score"]:.2f}', va="center", fontsize=7)

# ── 3c: MPNN score vs rank (Pareto front) ───────────────────────────────
ax = axes3[2]
ax.scatter(decoy_d["mpnn_score"], decoy_d["el_rank"], c=C_D, s=60,
           edgecolors="black", linewidths=0.5, zorder=3)
# Highlight Pareto-optimal candidates
pareto = []
sorted_by_mpnn = decoy_d.sort_values("mpnn_score")
best_el = float("inf")
for _, row in sorted_by_mpnn.iterrows():
    if row["el_rank"] < best_el:
        pareto.append(row)
        best_el = row["el_rank"]
if pareto:
    pareto_df = pd.DataFrame(pareto)
    ax.scatter(pareto_df["mpnn_score"], pareto_df["el_rank"],
               c="red", s=100, marker="D", edgecolors="black", linewidths=0.8, zorder=5,
               label=f"Pareto front ({len(pareto)})")
    ax.plot(pareto_df["mpnn_score"], pareto_df["el_rank"],
            color="red", linestyle="--", alpha=0.5, zorder=4)
    for _, row in pareto_df.iterrows():
        ax.annotate(row["sequence"], (row["mpnn_score"], row["el_rank"]),
                    fontsize=6.5, ha="left", xytext=(5, 5), textcoords="offset points")
ax.set_xlabel("MPNN Score (lower = better design)")
ax.set_ylabel("EL %Rank (lower = stronger binder)")
ax.set_title("Pareto Front: Design vs Binding", fontweight="bold")
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig(OUT / "GILGFVFTL_decoy_d_detail.png", dpi=180, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT / 'GILGFVFTL_decoy_d_detail.png'}")

print("\nAll figures saved to:", OUT)
