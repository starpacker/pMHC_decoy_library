"""
Decoy A Pipeline Visualization
===============================
Part 1: Matplotlib — Overall HLA-filtered database statistics
Part 2: Plotly/HTML — Interactive case studies for individual scan results
"""

import json
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd

# ── Paths ───────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "decoy_a"
OUT_DIR = PROJECT_ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)

HLA_PARQUET = DATA_DIR / "hla_filtered_HLA-A0201.parquet"
EXPRESSION_PARQUET = DATA_DIR / "gene_expression.parquet"
RESULTS_JSON = DATA_DIR / "decoy_a_results.json"

VITAL_ORGANS = [
    "heart muscle", "cerebral cortex", "lung",
    "liver", "kidney", "colon", "small intestine",
]


# =========================================================================
# Part 1:  Matplotlib — Overall Statistics
# =========================================================================

def plot_overview(df: pd.DataFrame):
    """4-panel overview figure of the HLA-filtered database."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Decoy A · HLA-A*02:01 Filtered Database Overview\n"
        f"(N = {len(df):,} binders from 41.7 M human k-mers, pass rate 3.07%)",
        fontsize=14, fontweight="bold", y=0.98,
    )

    # ── Panel A: EL%Rank distribution ────────────────────────────────────
    ax = axes[0, 0]
    ax.hist(df["el_rank"], bins=100, color="#4C72B0", edgecolor="white",
            linewidth=0.3, alpha=0.85)
    ax.axvline(0.5, color="#C44E52", ls="--", lw=1.5, label="Strong ≤ 0.5")
    ax.axvline(2.0, color="#DD8452", ls="--", lw=1.5, label="Weak ≤ 2.0")
    ax.set_xlabel("EL %Rank")
    ax.set_ylabel("Count")
    ax.set_title("A. EL %Rank Distribution")
    ax.legend(fontsize=9)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f}k"))

    # ── Panel B: Binding category pie ────────────────────────────────────
    ax = axes[0, 1]
    counts = df["binding"].value_counts()
    colors = {"Strong_Binder": "#C44E52", "Weak_Binder": "#DD8452"}
    labels = [f"{k.replace('_', ' ')}\n({v:,})" for k, v in counts.items()]
    wedges, texts, autotexts = ax.pie(
        counts.values, labels=labels,
        colors=[colors.get(k, "#999") for k in counts.index],
        autopct="%1.1f%%", startangle=90, textprops={"fontsize": 10},
    )
    ax.set_title("B. Binding Strength")

    # ── Panel C: Peptide length distribution ─────────────────────────────
    ax = axes[1, 0]
    df["pep_len"] = df["sequence"].str.len()
    len_counts = df["pep_len"].value_counts().sort_index()
    bars = ax.bar(len_counts.index, len_counts.values,
                  color=["#55A868", "#4C72B0", "#8172B2", "#CCB974"],
                  edgecolor="white", width=0.6)
    for bar, val in zip(bars, len_counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{val:,}", ha="center", va="bottom", fontsize=9)
    ax.set_xlabel("Peptide Length (aa)")
    ax.set_ylabel("Count")
    ax.set_title("C. Peptide Length Distribution")
    ax.set_xticks(len_counts.index)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e3:.0f}k"))

    # ── Panel D: EL%Rank by peptide length (violin) ──────────────────────
    ax = axes[1, 1]
    lengths = sorted(df["pep_len"].unique())
    data_by_len = [df.loc[df["pep_len"] == l, "el_rank"].values for l in lengths]
    vp = ax.violinplot(data_by_len, positions=lengths, showmedians=True,
                       showextrema=False, widths=0.5)
    for body in vp["bodies"]:
        body.set_facecolor("#8172B2")
        body.set_alpha(0.7)
    vp["cmedians"].set_color("#C44E52")
    ax.set_xlabel("Peptide Length (aa)")
    ax.set_ylabel("EL %Rank")
    ax.set_title("D. EL %Rank by Peptide Length")
    ax.set_xticks(lengths)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = OUT_DIR / "decoy_a_overview.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_funnel():
    """Pipeline diagram: generation step then filtering funnel."""
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(7, 7.6, "Decoy A Pipeline — Generation & Filtering",
            ha="center", va="center", fontsize=16, fontweight="bold",
            color="#2c3e50")

    # ── Step 1: Generation (not a filter — expansion) ────────────────────
    gen_stages = [
        ("Swiss-Prot\n20,431 proteins", "#4C72B0", 3.5),
        ("Sliding Window\n41.7 M unique k-mers", "#55A868", 5.0),
    ]
    y_gen = 6.2
    ax.add_patch(mpatches.FancyBboxPatch(
        (1.5, y_gen - 0.45), 3.5, 0.9,
        boxstyle="round,pad=0.12", facecolor="#4C72B0",
        edgecolor="white", linewidth=2, alpha=0.85))
    ax.text(3.25, y_gen, "Swiss-Prot Human\n20,431 proteins",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="white")

    ax.annotate("", xy=(6.2, y_gen), xytext=(5.0, y_gen),
                arrowprops=dict(arrowstyle="-|>", color="#7f8c8d", lw=2))
    ax.text(5.6, y_gen + 0.3, "8-11 mer\nsliding window",
            ha="center", va="center", fontsize=9, color="#555")

    ax.add_patch(mpatches.FancyBboxPatch(
        (6.2, y_gen - 0.45), 5.0, 0.9,
        boxstyle="round,pad=0.12", facecolor="#55A868",
        edgecolor="white", linewidth=2, alpha=0.85))
    ax.text(8.7, y_gen, "Unique K-mers\n41,700,000 peptide candidates",
            ha="center", va="center", fontsize=12, fontweight="bold",
            color="white")

    # ── Divider ──────────────────────────────────────────────────────────
    ax.plot([1, 13], [5.2, 5.2], color="#ddd", lw=1, ls="--")
    ax.text(7, 5.35, "FILTERING FUNNEL", ha="center", fontsize=10,
            color="#aaa", fontweight="bold")

    # ── Funnel stages (these actually filter down) ───────────────────────
    funnel = [
        ("HLA-A*02:01 Binders  (EL%Rank \u2264 2.0)", 1_280_221, "#DD8452"),
        ("Strong Binders  (EL%Rank \u2264 0.5)", 292_479, "#C44E52"),
        ("Decoy Hits  (Hamming \u2264 2)", 200, "#8B0000"),
    ]
    base_count = 41_700_000
    y_positions = [4.3, 3.0, 1.7]
    max_w = 10.0
    min_w = 4.0

    prev_count = base_count
    for i, ((label, count, color), y) in enumerate(zip(funnel, y_positions)):
        # Width proportional to log scale
        ratio = np.log10(count + 1) / np.log10(base_count + 1)
        w = max(max_w * ratio, min_w)
        cx = 7.0
        rect = mpatches.FancyBboxPatch(
            (cx - w / 2, y - 0.4), w, 0.8,
            boxstyle="round,pad=0.12", facecolor=color,
            edgecolor="white", linewidth=2, alpha=0.85,
        )
        ax.add_patch(rect)
        ax.text(cx, y, f"{label}\n{count:,}",
                ha="center", va="center", fontsize=11,
                fontweight="bold", color="white")

        # Pass rate from previous stage
        pct = count / prev_count * 100
        if i == 0:
            pct_text = f"{pct:.2f}% pass"
            from_label = f"of 41.7M"
        else:
            pct_text = f"{pct:.2f}% pass"
            from_label = f"of {prev_count:,}"

        # Arrow between stages
        if i == 0:
            ax.annotate("", xy=(cx, y + 0.4), xytext=(cx, 5.2),
                        arrowprops=dict(arrowstyle="-|>", color="#999", lw=1.5))
            ax.text(cx + max_w / 2 + 0.3, y + 0.7, f"{pct:.2f}% pass rate",
                    ha="left", va="center", fontsize=10, color="#888",
                    fontweight="bold")
        else:
            ax.annotate("", xy=(cx, y + 0.4), xytext=(cx, y_positions[i - 1] - 0.4),
                        arrowprops=dict(arrowstyle="-|>", color="#999", lw=1.5))
            ax.text(cx + max_w / 2 + 0.3, (y + y_positions[i - 1]) / 2,
                    f"{pct:.2f}% pass rate",
                    ha="left", va="center", fontsize=10, color="#888",
                    fontweight="bold")
        prev_count = count

    plt.tight_layout()
    out = OUT_DIR / "decoy_a_funnel.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_mismatch_position_heatmap(hits: list[dict]):
    """Heatmap: which positions get mutated across all hits."""
    if not hits:
        return
    pep_len = len(hits[0]["target_sequence"])
    target = hits[0]["target_sequence"]

    pos_counts = {"TCR contact": np.zeros(pep_len),
                  "Anchor": np.zeros(pep_len),
                  "Other": np.zeros(pep_len)}

    for h in hits:
        for m in h["mismatches"]:
            p = m["position"]
            if m["is_tcr_contact"]:
                pos_counts["TCR contact"][p] += 1
            elif m["is_anchor"]:
                pos_counts["Anchor"][p] += 1
            else:
                pos_counts["Other"][p] += 1

    fig, ax = plt.subplots(figsize=(10, 3.5))
    x = np.arange(pep_len)
    width = 0.25
    colors = {"TCR contact": "#C44E52", "Anchor": "#4C72B0", "Other": "#55A868"}
    for i, (label, vals) in enumerate(pos_counts.items()):
        ax.bar(x + i * width, vals, width, label=label, color=colors[label],
               edgecolor="white")

    ax.set_xticks(x + width)
    ax.set_xticklabels([f"p{i+1}\n{target[i]}" for i in range(pep_len)], fontsize=10)
    ax.set_ylabel("Mismatch Count")
    ax.set_title(
        f"Mismatch Position Distribution — Target: {target}  (N={len(hits)} hits)",
        fontsize=12, fontweight="bold",
    )
    ax.legend()
    plt.tight_layout()
    out = OUT_DIR / "decoy_a_mismatch_positions.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_scatter_elrank_vs_hamming(hits: list[dict]):
    """Scatter plot: EL%Rank vs Hamming distance, colored by mismatch type."""
    if not hits:
        return
    target = hits[0]["target_sequence"]
    fig, ax = plt.subplots(figsize=(8, 5))

    for h in hits:
        color = "#C44E52" if h["n_anchor_mismatches"] > 0 else (
            "#DD8452" if h["n_tcr_contact_mismatches"] > 0 else "#55A868"
        )
        marker = "^" if h["n_anchor_mismatches"] > 0 else "o"
        ax.scatter(h["hamming_distance"], h["el_rank"], c=color, marker=marker,
                   s=80, edgecolors="white", linewidth=0.5, zorder=3)
        ax.annotate(h["sequence"], (h["hamming_distance"], h["el_rank"]),
                    fontsize=7, xytext=(5, 3), textcoords="offset points")

    ax.axhline(0.5, color="#C44E52", ls=":", lw=1, alpha=0.5, label="Strong ≤ 0.5")
    ax.axhline(2.0, color="#DD8452", ls=":", lw=1, alpha=0.5, label="Weak ≤ 2.0")

    # legend handles
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="^", color="w", markerfacecolor="#C44E52",
               markersize=10, label="Anchor mismatch"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#DD8452",
               markersize=10, label="TCR contact mismatch only"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#55A868",
               markersize=10, label="Non-contact mismatch only"),
    ]
    ax.legend(handles=handles, fontsize=9)
    ax.set_xlabel("Hamming Distance")
    ax.set_ylabel("EL %Rank")
    ax.set_title(f"EL %Rank vs Hamming Distance — Target: {target}",
                 fontsize=12, fontweight="bold")
    ax.set_xticks([0, 1, 2])
    plt.tight_layout()
    out = OUT_DIR / "decoy_a_scatter.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


# =========================================================================
# Part 2:  Plotly / HTML — Interactive Case Studies
# =========================================================================

def generate_interactive_html(hits: list[dict]):
    """Generate a single interactive HTML report with two case studies."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    target = hits[0]["target_sequence"]

    # ── Pick two example hits ────────────────────────────────────────────
    # Case 1: the one with real expression data (EPS8L2)
    case1 = next((h for h in hits if h.get("expression")), hits[0])
    # Case 2: an anchor-mismatch hit
    case2 = next((h for h in hits if h["n_anchor_mismatches"] > 0), hits[-1])

    html_parts = []
    html_parts.append(_html_header(target, len(hits)))

    # ── Case Study 1: High-risk with expression ─────────────────────────
    html_parts.append(_case_section(case1, 1, target))

    # ── Case Study 2: Anchor mismatch ───────────────────────────────────
    html_parts.append(_case_section(case2, 2, target))

    # ── Close ────────────────────────────────────────────────────────────
    html_parts.append("</div></body></html>")

    out = OUT_DIR / "decoy_a_interactive.html"
    out.write_text("\n".join(html_parts), encoding="utf-8")
    print(f"[plotly/html] Saved: {out}")


def _html_header(target: str, n_hits: int) -> str:
    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>Decoy A — Interactive Case Studies</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
  body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: #f8f9fa;
         margin: 0; padding: 20px; color: #333; }}
  .container {{ max-width: 1100px; margin: 0 auto; }}
  h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
  h2 {{ color: #2c3e50; margin-top: 40px; }}
  .case-card {{ background: white; border-radius: 12px; padding: 24px;
                margin: 20px 0; box-shadow: 0 2px 12px rgba(0,0,0,0.08); }}
  .seq-align {{ font-family: 'Consolas', 'Courier New', monospace; font-size: 18px;
                line-height: 2.2; background: #1a1a2e; color: #e0e0e0;
                padding: 20px 24px; border-radius: 8px; letter-spacing: 4px; }}
  .seq-align .match {{ color: #a0a0a0; }}
  .seq-align .mismatch-tcr {{ color: #ff6b6b; font-weight: bold;
                               background: rgba(255,107,107,0.15);
                               border-radius: 3px; padding: 2px 1px; }}
  .seq-align .mismatch-anchor {{ color: #ffd93d; font-weight: bold;
                                  background: rgba(255,217,61,0.15);
                                  border-radius: 3px; padding: 2px 1px; }}
  .seq-align .mismatch-other {{ color: #6bcb77; font-weight: bold;
                                 background: rgba(107,203,119,0.15);
                                 border-radius: 3px; padding: 2px 1px; }}
  .label-row {{ color: #888; font-size: 13px; letter-spacing: 4px; }}
  .meta-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 12px; margin: 16px 0; }}
  .meta-item {{ background: #f0f4f8; border-radius: 8px; padding: 12px 16px; }}
  .meta-item .val {{ font-size: 22px; font-weight: bold; color: #2c3e50; }}
  .meta-item .lbl {{ font-size: 12px; color: #7f8c8d; text-transform: uppercase; }}
  .risk-high {{ border-left: 4px solid #e74c3c; }}
  .risk-low {{ border-left: 4px solid #27ae60; }}
  .legend {{ font-size: 13px; color: #666; margin-top: 8px; }}
  .legend span {{ padding: 2px 8px; border-radius: 3px; margin: 0 4px; }}
  .legend .tcr {{ background: rgba(255,107,107,0.2); color: #ff6b6b; }}
  .legend .anchor {{ background: rgba(255,217,61,0.2); color: #c89b00; }}
  .legend .other {{ background: rgba(107,203,119,0.2); color: #27ae60; }}
</style>
</head><body><div class="container">
<h1>Decoy A — Interactive Case Studies</h1>
<p>Target: <strong>{target}</strong> &nbsp;|&nbsp; Total hits: <strong>{n_hits}</strong>
&nbsp;|&nbsp; HLA: HLA-A*01:01</p>
<div class="legend">
  Mismatch types:
  <span class="tcr">TCR Contact</span>
  <span class="anchor">Anchor</span>
  <span class="other">Other</span>
</div>
"""


def _case_section(hit: dict, case_num: int, target: str) -> str:
    """Build HTML + Plotly section for one case study."""
    import plotly.graph_objects as go

    seq = hit["sequence"]
    mismatches = {m["position"]: m for m in hit["mismatches"]}
    gene = hit["gene_symbols"][0] if hit["gene_symbols"] else "Unknown"
    expr = hit.get("expression")
    risk = expr["expression_category"] if expr else "no data"
    risk_class = "risk-high" if risk == "high_risk" else "risk-low"

    # ── Sequence alignment HTML ──────────────────────────────────────────
    pos_labels = "".join(f"p{i+1:<3d}" for i in range(len(target)))
    target_html = ""
    decoy_html = ""
    for i in range(len(target)):
        if i in mismatches:
            m = mismatches[i]
            cls = ("mismatch-tcr" if m["is_tcr_contact"]
                   else "mismatch-anchor" if m["is_anchor"]
                   else "mismatch-other")
            target_html += f'<span class="{cls}">{target[i]}</span>'
            decoy_html += f'<span class="{cls}">{seq[i]}</span>'
        else:
            target_html += f'<span class="match">{target[i]}</span>'
            decoy_html += f'<span class="match">{seq[i]}</span>'

    section = f"""
<div class="case-card {risk_class}">
  <h2>Case {case_num}: {seq} — Gene: {gene}</h2>
  <div class="seq-align">
    <div class="label-row">{"".join(f"p{i+1:<3d}" for i in range(len(target)))}</div>
    <div>Target:&nbsp; {target_html}</div>
    <div>Decoy:&nbsp;&nbsp; {decoy_html}</div>
  </div>
  <div class="meta-grid">
    <div class="meta-item"><div class="val">{hit['hamming_distance']}</div><div class="lbl">Hamming Distance</div></div>
    <div class="meta-item"><div class="val">{hit['el_rank']:.2f}</div><div class="lbl">EL %Rank</div></div>
    <div class="meta-item"><div class="val">{hit['n_tcr_contact_mismatches']}</div><div class="lbl">TCR Contact Mismatches</div></div>
    <div class="meta-item"><div class="val">{hit['n_anchor_mismatches']}</div><div class="lbl">Anchor Mismatches</div></div>
    <div class="meta-item"><div class="val">{hit['similarity_score']:.1%}</div><div class="lbl">Sequence Similarity</div></div>
    <div class="meta-item"><div class="val" style="font-size:16px">{risk.replace('_', ' ').title()}</div><div class="lbl">Expression Risk</div></div>
  </div>
"""

    # ── Expression heatmap (Plotly) ──────────────────────────────────────
    if expr and expr.get("tissue_tpm"):
        tpm = expr["tissue_tpm"]
        tissues = sorted(tpm.keys(), key=lambda t: tpm[t], reverse=True)
        values = [tpm[t] for t in tissues]
        colors = []
        for t in tissues:
            if t in VITAL_ORGANS:
                colors.append("#e74c3c" if tpm[t] >= 10 else "#f39c12")
            else:
                colors.append("#3498db")

        fig = go.Figure(go.Bar(
            y=tissues, x=values, orientation="h",
            marker=dict(color=colors, line=dict(width=0.5, color="white")),
            hovertemplate="%{y}: %{x:.1f} TPM<extra></extra>",
        ))
        fig.add_vline(x=10, line_dash="dash", line_color="#e74c3c",
                      annotation_text="High Risk (10 TPM)",
                      annotation_position="top right",
                      annotation_font_size=10)
        fig.update_layout(
            title=dict(text=f"Tissue Expression — {gene}", font_size=14),
            xaxis_title="TPM (nTPM consensus)",
            height=max(600, len(tissues) * 16),
            width=700,
            margin=dict(l=150, r=30, t=50, b=40),
            yaxis=dict(autorange="reversed"),
            plot_bgcolor="white",
        )
        div_id = f"expr_plot_{case_num}"
        plot_json = fig.to_json()
        section += f"""
  <div id="{div_id}" style="margin-top: 16px;"></div>
  <script>
  (function() {{
    var spec = {plot_json};
    Plotly.newPlot('{div_id}', spec.data, spec.layout);
  }})();
  </script>
"""
    else:
        section += """
  <p style="color: #999; font-style: italic; margin-top: 16px;">
    No tissue expression data available (demo entry)
  </p>
"""

    section += "</div>"
    return section


# =========================================================================
# Main
# =========================================================================

def main():
    print("=" * 60)
    print("  Decoy A Pipeline — Visualization Suite")
    print("=" * 60)

    # Load data
    print("\nLoading HLA-filtered database...")
    df = pd.read_parquet(HLA_PARQUET)
    print(f"  {len(df):,} binders loaded")

    print("Loading scan results...")
    with open(RESULTS_JSON, encoding="utf-8") as f:
        hits = json.load(f)
    print(f"  {len(hits)} hits loaded")

    # Part 1: Matplotlib
    print("\n── Part 1: Matplotlib Overview Figures ──")
    plot_overview(df)
    plot_funnel()
    plot_mismatch_position_heatmap(hits)
    plot_scatter_elrank_vs_hamming(hits)

    # Part 2: Plotly HTML
    print("\n── Part 2: Plotly Interactive Case Studies ──")
    generate_interactive_html(hits)

    print(f"\nAll figures saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
