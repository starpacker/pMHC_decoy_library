"""
Decoy C Library Visualization
===============================
Part 1: Matplotlib — Pipeline overview, library statistics, filtering logic
Part 2: Plotly/HTML — Interactive case study (DC-0001: MAGE-A3 / Titin)
"""

import json
from collections import Counter
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Paths ───────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "decoy_c"
OUT_DIR = PROJECT_ROOT / "figures"
OUT_DIR.mkdir(exist_ok=True)

LIBRARY_JSON = DATA_DIR / "decoy_library.json"


def load_library() -> list[dict]:
    with open(LIBRARY_JSON, encoding="utf-8") as f:
        lib = json.load(f)
    return lib.get("entries", lib) if isinstance(lib, dict) else lib


# =========================================================================
# Part 1:  Matplotlib — Pipeline & Library Statistics
# =========================================================================

def plot_pipeline_flowchart():
    """Show the Decoy C pipeline as a clear step-by-step flowchart."""
    fig, ax = plt.subplots(figsize=(20, 13))
    ax.set_xlim(0, 20)
    ax.set_ylim(-0.5, 12.5)
    ax.axis("off")

    # ── Title ────────────────────────────────────────────────────────────
    ax.text(10, 11.4, "Decoy C Pipeline — Literature-Mined Cross-Reactivity Database",
            ha="center", va="center", fontsize=20, fontweight="bold", color="#2c3e50")

    # ── Stage boxes ──────────────────────────────────────────────────────
    stages = [
        {
            "x": 0.8, "y": 8.5, "w": 4.5, "h": 2.4,
            "color": "#3498db", "title": "Stage 1: Discovery",
            "subtitle": "PubMed",
            "lines": [
                '8 keyword queries',
                'e.g. "TCR" AND "cross-reactivity"',
                'NCBI E-utilities API',
                '~200 papers retrieved',
            ],
        },
        {
            "x": 6.8, "y": 8.5, "w": 4.8, "h": 2.4,
            "color": "#9b59b6", "title": "Stage 2: LLM Extraction",
            "subtitle": "GPT-4o",
            "lines": [
                'Title + Abstract  →  JSON',
                'Structured DecoyEntry output',
                'Temperature = 0.1 (strict)',
                '~250 raw entries',
            ],
        },
        {
            "x": 13.2, "y": 8.5, "w": 5.8, "h": 2.4,
            "color": "#e67e22", "title": "Stage 3: Validation",
            "subtitle": "UniProt + IEDB",
            "lines": [
                'Gene → UniProt (human, reviewed)',
                'Peptide → IEDB epitope DB',
                'Mass-spec confirmation check',
                '→ VALIDATED / PARTIAL / REVIEW',
            ],
        },
        {
            "x": 7.0, "y": 3.5, "w": 6.0, "h": 2.4,
            "color": "#27ae60", "title": "Stage 4: Library Curation",
            "subtitle": "JSON DB",
            "lines": [
                'Deduplicate by peptide sequence',
                'Assign DC-NNNN sequential IDs',
                '4-tier evidence classification',
                '250 curated entries (current)',
            ],
        },
    ]

    for s in stages:
        rect = mpatches.FancyBboxPatch(
            (s["x"], s["y"] - s["h"] / 2), s["w"], s["h"],
            boxstyle="round,pad=0.15", facecolor=s["color"],
            edgecolor="white", linewidth=2.5, alpha=0.9,
        )
        ax.add_patch(rect)
        # Title
        ax.text(s["x"] + s["w"] / 2, s["y"] + s["h"] / 2 - 0.22,
                s["title"], ha="center", va="center",
                fontsize=14, fontweight="bold", color="white")
        # Subtitle badge
        ax.text(s["x"] + s["w"] / 2, s["y"] + s["h"] / 2 - 0.6,
                s["subtitle"], ha="center", va="center",
                fontsize=10, color="white", alpha=0.7,
                fontstyle="italic")
        # Detail lines
        for i, line in enumerate(s["lines"]):
            ax.text(s["x"] + 0.35, s["y"] - 0.15 - i * 0.35,
                    line, ha="left", va="center",
                    fontsize=11, color="white", alpha=0.95)

    # ── Arrows between stages ────────────────────────────────────────────
    arrow_kw = dict(arrowstyle="-|>", color="#7f8c8d", lw=2.5,
                    connectionstyle="arc3,rad=0")
    # Stage 1 → Stage 2
    ax.annotate("", xy=(6.8, 8.5), xytext=(5.3, 8.5),
                arrowprops=arrow_kw)
    ax.text(6.05, 8.9, "papers", ha="center", fontsize=11, color="#555",
            fontweight="bold")

    # Stage 2 → Stage 3
    ax.annotate("", xy=(13.2, 8.5), xytext=(11.6, 8.5),
                arrowprops=arrow_kw)
    ax.text(12.4, 8.9, "raw entries", ha="center", fontsize=11, color="#555",
            fontweight="bold")

    # Stage 3 → Stage 4  (down-arrow)
    ax.annotate("", xy=(10.5, 5.7), xytext=(16.1, 7.3),
                arrowprops=dict(arrowstyle="-|>", color="#7f8c8d", lw=2.5,
                                connectionstyle="arc3,rad=-0.3"))
    ax.text(14.0, 6.6, "validated", ha="center", fontsize=11, color="#555",
            fontweight="bold")

    # ── Seed data input arrow ────────────────────────────────────────────
    seed_rect = mpatches.FancyBboxPatch(
        (1.0, 3.1), 4.2, 1.6,
        boxstyle="round,pad=0.15", facecolor="#95a5a6",
        edgecolor="white", linewidth=2.5, alpha=0.85,
    )
    ax.add_patch(seed_rect)
    ax.text(3.1, 4.25, "Seed Data", ha="center", va="center",
            fontsize=14, fontweight="bold", color="white")
    ax.text(3.1, 3.55, "51 hand-curated entries\nfrom key publications",
            ha="center", va="center", fontsize=11, color="white", alpha=0.9)

    ax.annotate("", xy=(7.0, 4.0), xytext=(5.2, 4.0),
                arrowprops=dict(arrowstyle="-|>", color="#7f8c8d", lw=2.5))
    ax.text(6.1, 4.4, "merge", ha="center", fontsize=11, color="#555",
            fontweight="bold")

    # ── Evidence level legend (bottom) ───────────────────────────────────
    ev_colors = {
        "Level 1\nClinical Fatal": "#c0392b",
        "Level 2\nIn Vitro Confirmed": "#e74c3c",
        "Level 3\nHT Screened": "#f39c12",
        "Level 4\nIn Silico": "#3498db",
    }
    y_legend = 0.8
    ax.text(10, 1.8, "Evidence Classification (4-Tier)",
            ha="center", fontsize=14, fontweight="bold", color="#2c3e50")
    for i, (label, col) in enumerate(ev_colors.items()):
        x = 1.5 + i * 4.5
        rect = mpatches.FancyBboxPatch(
            (x, y_legend - 0.5), 3.8, 1.1,
            boxstyle="round,pad=0.1", facecolor=col,
            edgecolor="white", linewidth=2, alpha=0.85,
        )
        ax.add_patch(rect)
        ax.text(x + 1.9, y_legend + 0.05, label,
                ha="center", va="center", fontsize=12,
                fontweight="bold", color="white")

    plt.tight_layout()
    out = OUT_DIR / "decoy_c_pipeline.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_library_overview(entries: list[dict]):
    """4-panel overview of the library contents."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Decoy C Library Overview — {len(entries)} Curated Entries",
        fontsize=14, fontweight="bold", y=0.98,
    )

    # ── Panel A: Evidence Level Distribution ─────────────────────────────
    ax = axes[0, 0]
    ev_counter = Counter()
    for e in entries:
        lv = e.get("risk_profile", {}).get("evidence_level", "Unknown")
        ev_counter[lv] += 1

    ev_order = [
        "Level_1_Clinical_Fatal",
        "Level_2_In_Vitro_Confirmed",
        "Level_3_High_Throughput_Screened",
        "Level_4_In_Silico_High_Risk",
        "Unknown",
    ]
    ev_labels = ["L1: Fatal", "L2: In Vitro", "L3: HT Screen", "L4: In Silico", "Unknown"]
    ev_colors = ["#c0392b", "#e74c3c", "#f39c12", "#3498db", "#bdc3c7"]
    ev_vals = [ev_counter.get(k, 0) for k in ev_order]

    bars = ax.barh(ev_labels[::-1], ev_vals[::-1],
                   color=ev_colors[::-1], edgecolor="white", height=0.6)
    for bar, val in zip(bars, ev_vals[::-1]):
        if val > 0:
            ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height() / 2,
                    str(val), ha="left", va="center", fontsize=10, fontweight="bold")
    ax.set_xlabel("Count")
    ax.set_title("A. Evidence Level Distribution")

    # ── Panel B: Validation Status ───────────────────────────────────────
    ax = axes[0, 1]
    val_counter = Counter()
    for e in entries:
        vs = e.get("validation_flags", {}).get("overall_status", "Unknown")
        val_counter[vs] += 1

    val_order = ["VALIDATED", "PARTIAL", "NEEDS_REVIEW", "Unknown"]
    val_colors_map = {"VALIDATED": "#27ae60", "PARTIAL": "#f39c12",
                      "NEEDS_REVIEW": "#e74c3c", "Unknown": "#bdc3c7"}
    val_vals = [val_counter.get(k, 0) for k in val_order]
    val_labels_fmt = [f"{k}\n({v})" for k, v in zip(val_order, val_vals) if v > 0]
    val_vals_nz = [v for v in val_vals if v > 0]
    val_cols_nz = [val_colors_map[k] for k, v in zip(val_order, val_vals) if v > 0]

    ax.pie(val_vals_nz, labels=val_labels_fmt, colors=val_cols_nz,
           autopct="%1.1f%%", startangle=90, textprops={"fontsize": 9})
    ax.set_title("B. Validation Status")

    # ── Panel C: HLA Allele Distribution (top 10) ────────────────────────
    ax = axes[1, 0]
    hla_counter = Counter()
    for e in entries:
        h = e.get("peptide_info", {}).get("hla_allele", "Unknown")
        hla_counter[h] += 1

    top_hla = hla_counter.most_common(10)
    hla_names = [h[0] for h in top_hla][::-1]
    hla_vals = [h[1] for h in top_hla][::-1]
    cmap = plt.cm.Set2(np.linspace(0, 1, len(hla_names)))
    bars = ax.barh(hla_names, hla_vals, color=cmap, edgecolor="white", height=0.6)
    for bar, val in zip(bars, hla_vals):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height() / 2,
                str(val), ha="left", va="center", fontsize=9, fontweight="bold")
    ax.set_xlabel("Count")
    ax.set_title("C. HLA Allele Distribution (Top 10)")

    # ── Panel D: Critical Organs Affected (top 10) ───────────────────────
    ax = axes[1, 1]
    organ_counter = Counter()
    for e in entries:
        for o in e.get("risk_profile", {}).get("critical_organs_affected", []):
            organ_counter[o] += 1

    top_organs = organ_counter.most_common(10)
    organ_names = [o[0] for o in top_organs][::-1]
    organ_vals = [o[1] for o in top_organs][::-1]

    vital_set = {"Heart", "Brain", "Lung", "Liver", "Kidney"}
    organ_colors = ["#e74c3c" if o in vital_set else "#3498db" for o in organ_names]
    bars = ax.barh(organ_names, organ_vals, color=organ_colors,
                   edgecolor="white", height=0.6)
    for bar, val in zip(bars, organ_vals):
        ax.text(bar.get_width() + 0.2, bar.get_y() + bar.get_height() / 2,
                str(val), ha="left", va="center", fontsize=9, fontweight="bold")
    ax.set_xlabel("Count")
    ax.set_title("D. Critical Organs Affected (Top 10)")

    # Legend for organ colors
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#e74c3c",
               markersize=10, label="Vital organ"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#3498db",
               markersize=10, label="Other organ"),
    ]
    ax.legend(handles=handles, fontsize=8, loc="lower right")

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = OUT_DIR / "decoy_c_overview.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_filtering_sankey(entries: list[dict]):
    """Illustrate the filtering/enrichment flow as a horizontal funnel."""
    fig, ax = plt.subplots(figsize=(18, 8))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(9, 7.5, "Decoy C — Filtering & Enrichment Flow",
            ha="center", va="center", fontsize=18, fontweight="bold",
            color="#2c3e50")

    # Count stats
    total = len(entries)
    known_ev = sum(1 for e in entries
                   if e.get("risk_profile", {}).get("evidence_level", "Unknown") != "Unknown")
    validated = sum(1 for e in entries
                    if e.get("validation_flags", {}).get("overall_status") == "VALIDATED")
    ms_conf = sum(1 for e in entries
                  if e.get("experimental_evidence", {}).get("mass_spec_confirmed"))

    # Stages as blocks (left to right)
    blocks = [
        {"x": 0.3, "y": 4.2, "w": 3.2, "h": 3.8,
         "color": "#3498db", "alpha": 0.85,
         "title": "PubMed\nDiscovery",
         "detail": f"~200 papers\nscreened"},
        {"x": 4.0, "y": 4.2, "w": 3.2, "h": 3.2,
         "color": "#9b59b6", "alpha": 0.85,
         "title": "LLM\nExtraction",
         "detail": f"{total} entries\nextracted"},
        {"x": 7.7, "y": 4.2, "w": 3.2, "h": 2.5,
         "color": "#e67e22", "alpha": 0.85,
         "title": "UniProt\nValidation",
         "detail": f"{known_ev} with known\nevidence level"},
        {"x": 11.4, "y": 4.2, "w": 3.0, "h": 1.9,
         "color": "#27ae60", "alpha": 0.85,
         "title": "IEDB\nConfirmed",
         "detail": f"{validated} VALIDATED\n{ms_conf} mass-spec"},
        {"x": 14.8, "y": 4.2, "w": 2.8, "h": 1.3,
         "color": "#c0392b", "alpha": 0.85,
         "title": "Clinical\nEvidence",
         "detail": "5 fatal cases"},
    ]

    for b in blocks:
        cy = b["y"]
        rect = mpatches.FancyBboxPatch(
            (b["x"], cy - b["h"] / 2), b["w"], b["h"],
            boxstyle="round,pad=0.18", facecolor=b["color"],
            edgecolor="white", linewidth=2.5, alpha=b["alpha"],
        )
        ax.add_patch(rect)
        ax.text(b["x"] + b["w"] / 2, cy + 0.25, b["title"],
                ha="center", va="center", fontsize=14,
                fontweight="bold", color="white")
        ax.text(b["x"] + b["w"] / 2, cy - 0.5, b["detail"],
                ha="center", va="center", fontsize=11, color="white", alpha=0.95)

    # Arrows
    arrow_xs = [(3.5, 4.0), (7.2, 7.7), (10.9, 11.4), (14.4, 14.8)]
    for x_from, x_to in arrow_xs:
        ax.annotate("", xy=(x_to, 4.2), xytext=(x_from, 4.2),
                    arrowprops=dict(arrowstyle="-|>", color="#7f8c8d", lw=2.5))

    # Bottom annotations — key filtering criteria
    criteria = [
        (1.9, "Keywords:\nTCR, cross-reactivity\ntoxicity, off-target"),
        (5.6, "Extraction:\nSequence (8-15 AA)\nHLA, Gene, Organ\nEvidence level"),
        (9.3, "Filter:\ngene in UniProt\norganism = Human\nreviewed = true"),
        (12.9, "Confirm:\npeptide in IEDB\nmass-spec assays\nHLA match"),
        (16.2, "Classify:\n4-tier evidence\norgan risk\nexpression"),
    ]
    for x, text in criteria:
        ax.text(x, 1.2, text, ha="center", va="center", fontsize=10,
                color="#444",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                          edgecolor="#ccc", alpha=0.9))

    plt.tight_layout()
    out = OUT_DIR / "decoy_c_filtering.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


def plot_assay_and_year(entries: list[dict]):
    """Two-panel: assay type frequency + publication year timeline."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # ── Left: Top assay types ────────────────────────────────────────────
    ax = axes[0]
    assay_counter = Counter()
    for e in entries:
        for a in e.get("experimental_evidence", {}).get("assays_performed", []):
            clean = a.split("|")[0].strip()
            assay_counter[clean] += 1

    top_assays = assay_counter.most_common(12)
    names = [a[0] for a in top_assays][::-1]
    vals = [a[1] for a in top_assays][::-1]
    ax.barh(names, vals, color="#8e44ad", edgecolor="white", height=0.6, alpha=0.85)
    for i, (n, v) in enumerate(zip(names, vals)):
        ax.text(v + 1, i, str(v), ha="left", va="center", fontsize=9)
    ax.set_xlabel("Occurrences across entries")
    ax.set_title("A. Top Experimental Assay Types", fontsize=12, fontweight="bold")

    # ── Right: Publication year timeline ─────────────────────────────────
    ax = axes[1]
    year_counter = Counter()
    for e in entries:
        s = e.get("source", {})
        if s and s.get("year"):
            year_counter[s["year"]] += 1

    if year_counter:
        years = sorted(year_counter.keys())
        counts = [year_counter[y] for y in years]
        ax.bar(years, counts, color="#2980b9", edgecolor="white", width=0.7, alpha=0.85)
        # Cumulative line
        ax2 = ax.twinx()
        cumulative = np.cumsum(counts)
        ax2.plot(years, cumulative, color="#e74c3c", lw=2, marker="o", markersize=4)
        ax2.set_ylabel("Cumulative entries", color="#e74c3c")
        ax2.tick_params(axis="y", labelcolor="#e74c3c")
    ax.set_xlabel("Publication Year")
    ax.set_ylabel("Entries")
    ax.set_title("B. Publication Year Distribution", fontsize=12, fontweight="bold")

    plt.tight_layout()
    out = OUT_DIR / "decoy_c_assay_year.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[matplotlib] Saved: {out}")


# =========================================================================
# Part 2:  Plotly / HTML — Interactive Case Study
# =========================================================================

def generate_interactive_html(entries: list[dict]):
    """Interactive HTML report for DC-0001 (MAGE-A3/Titin) case study."""
    import plotly.graph_objects as go

    # Find DC-0001
    case = next((e for e in entries if e["decoy_id"] == "DC-0001"), entries[0])

    # Also collect all entries related to MAGE-A3 target
    mage_entries = [e for e in entries
                    if (e.get("discovery_context", {}).get("original_target_protein") or "")
                    .upper().startswith("MAGE")]

    html_parts = [_html_header()]

    # ── Section 1: Case card ─────────────────────────────────────────────
    html_parts.append(_case_card(case))

    # ── Section 2: MAGE-A3 related entries interactive table ─────────────
    html_parts.append(_mage_network_section(mage_entries))

    # ── Section 3: Evidence level sunburst for full library ──────────────
    html_parts.append(_evidence_sunburst(entries))

    html_parts.append("</div></body></html>")

    out = OUT_DIR / "decoy_c_interactive.html"
    out.write_text("\n".join(html_parts), encoding="utf-8")
    print(f"[plotly/html] Saved: {out}")


def _html_header() -> str:
    return """<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8">
<title>Decoy C — Interactive Case Study</title>
<script src="https://cdn.plot.ly/plotly-2.35.0.min.js"></script>
<style>
  body { font-family: 'Segoe UI', system-ui, sans-serif; background: #f8f9fa;
         margin: 0; padding: 20px; color: #333; }
  .container { max-width: 1200px; margin: 0 auto; }
  h1 { color: #2c3e50; border-bottom: 3px solid #9b59b6; padding-bottom: 10px; }
  h2 { color: #2c3e50; margin-top: 40px; }
  .case-card { background: white; border-radius: 12px; padding: 28px;
               margin: 20px 0; box-shadow: 0 2px 12px rgba(0,0,0,0.08);
               border-left: 5px solid #c0392b; }
  .meta-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
               gap: 12px; margin: 16px 0; }
  .meta-item { background: #f0f4f8; border-radius: 8px; padding: 14px 18px; }
  .meta-item .val { font-size: 20px; font-weight: bold; color: #2c3e50; }
  .meta-item .lbl { font-size: 11px; color: #7f8c8d; text-transform: uppercase;
                     margin-top: 2px; }
  .badge { display: inline-block; padding: 4px 12px; border-radius: 20px;
           font-size: 12px; font-weight: bold; color: white; margin: 2px 4px; }
  .badge-fatal { background: #c0392b; }
  .badge-validated { background: #27ae60; }
  .badge-ms { background: #8e44ad; }
  .badge-assay { background: #2980b9; }
  .seq-display { font-family: 'Consolas', monospace; font-size: 28px;
                  letter-spacing: 6px; background: #1a1a2e; color: #ffd93d;
                  padding: 16px 24px; border-radius: 8px; display: inline-block;
                  margin: 10px 0; }
  .provenance { background: #fef9e7; border-radius: 8px; padding: 16px;
                margin-top: 16px; border-left: 4px solid #f39c12; }
  .provenance a { color: #2980b9; text-decoration: none; }
  .provenance a:hover { text-decoration: underline; }
  table { width: 100%; border-collapse: collapse; margin: 12px 0; }
  th { background: #2c3e50; color: white; padding: 10px 12px; text-align: left;
       font-size: 12px; }
  td { padding: 8px 12px; border-bottom: 1px solid #eee; font-size: 13px; }
  tr:hover { background: #f5f5f5; }
  .section-card { background: white; border-radius: 12px; padding: 24px;
                  margin: 20px 0; box-shadow: 0 2px 12px rgba(0,0,0,0.08); }
</style>
</head><body><div class="container">
<h1>Decoy C Library — Interactive Case Study</h1>
<p>Literature-mined cross-reactivity database with experimental evidence classification</p>
"""


def _case_card(case: dict) -> str:
    pi = case["peptide_info"]
    dc = case.get("discovery_context", {})
    rp = case["risk_profile"]
    ee = case.get("experimental_evidence", {})
    prov = case.get("provenance", {})
    src = case.get("source", {})
    vf = case.get("validation_flags", {})

    # Format PMIDs as links
    pmid_links = " ".join(
        f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" target="_blank">PMID:{p}</a>'
        for p in prov.get("pmid", [])[:6]
    )
    if len(prov.get("pmid", [])) > 6:
        pmid_links += f" ... +{len(prov['pmid']) - 6} more"

    # Format clinical trial IDs
    trial_links = " ".join(
        f'<a href="https://clinicaltrials.gov/study/{t}" target="_blank">{t}</a>'
        for t in prov.get("clinical_trial_id", [])
    )

    # Format assays as badges
    key_assays = ["Clinical Trial", "Cytotoxicity", "IFN-gamma", "IFNg release",
                  "dissociation constant KD", "3D structure", "ligand presentation"]
    assay_badges = ""
    for a in ee.get("assays_performed", []):
        clean = a.split("|")[0].strip()
        if clean in key_assays:
            assay_badges += f'<span class="badge badge-assay">{clean}</span>'

    target_seq = dc.get("original_target_sequence", "N/A")
    target_display = ""
    if target_seq != "N/A" and len(target_seq) <= 15:
        # Highlight differences between target and decoy
        decoy_seq = pi["decoy_sequence"]
        for i in range(min(len(target_seq), len(decoy_seq))):
            if target_seq[i] != decoy_seq[i]:
                target_display += f'<span style="color:#e74c3c;font-weight:bold">{target_seq[i]}</span>'
            else:
                target_display += target_seq[i]
        target_display = f'<span style="font-family:Consolas;letter-spacing:4px;font-size:18px">{target_display}</span>'
    else:
        target_display = target_seq

    decoy_display = ""
    if target_seq != "N/A" and len(target_seq) <= 15:
        decoy_seq = pi["decoy_sequence"]
        for i in range(min(len(target_seq), len(decoy_seq))):
            if target_seq[i] != decoy_seq[i]:
                decoy_display += f'<span style="color:#ffd93d;font-weight:bold">{decoy_seq[i]}</span>'
            else:
                decoy_display += f'<span style="opacity:0.6">{decoy_seq[i]}</span>'
        decoy_display = f'<span style="font-family:Consolas;letter-spacing:4px;font-size:18px">{decoy_display}</span>'

    return f"""
<h2>Case Study: {case['decoy_id']} — The MAGE-A3 / Titin Cardiac Fatality</h2>
<div class="case-card">
  <div style="display:flex; align-items:center; gap:12px; margin-bottom:12px;">
    <span class="badge badge-fatal">Level 1: Clinical Fatal</span>
    <span class="badge badge-validated">VALIDATED</span>
    <span class="badge badge-ms">Mass-Spec Confirmed</span>
  </div>

  <div style="background:#1a1a2e; padding:18px 24px; border-radius:8px; margin:12px 0;
              font-family:Consolas; color:#e0e0e0; font-size:15px; line-height:2.4;">
    <div style="color:#888; font-size:12px;">
      &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
      {"&nbsp;&nbsp;".join(f"p{i+1}" for i in range(len(pi['decoy_sequence'])))}
    </div>
    <div>Target (MAGE-A3):&nbsp; {target_display}</div>
    <div>Decoy&nbsp; (Titin):&nbsp;&nbsp;&nbsp; {decoy_display}</div>
  </div>

  <div class="meta-grid">
    <div class="meta-item">
      <div class="val">{pi['decoy_sequence']}</div>
      <div class="lbl">Decoy Peptide</div>
    </div>
    <div class="meta-item">
      <div class="val">{pi['hla_allele']}</div>
      <div class="lbl">HLA Allele</div>
    </div>
    <div class="meta-item">
      <div class="val">{pi['source_protein']} ({pi['gene_symbol']})</div>
      <div class="lbl">Source Protein</div>
    </div>
    <div class="meta-item">
      <div class="val">{pi.get('uniprot_id', 'N/A')}</div>
      <div class="lbl">UniProt ID</div>
    </div>
    <div class="meta-item">
      <div class="val">{dc.get('original_target_protein', 'N/A')}</div>
      <div class="lbl">Original Target</div>
    </div>
    <div class="meta-item">
      <div class="val">{dc.get('tcr_name_or_id', 'N/A')}</div>
      <div class="lbl">TCR Clone</div>
    </div>
    <div class="meta-item">
      <div class="val">{', '.join(rp.get('critical_organs_affected', []))}</div>
      <div class="lbl">Organs Affected</div>
    </div>
    <div class="meta-item">
      <div class="val">{(rp.get('expression_pattern') or 'N/A').replace('_', ' ')}</div>
      <div class="lbl">Expression Pattern</div>
    </div>
  </div>

  <div style="margin-top:12px;">
    <strong>Key Assays:</strong><br>{assay_badges}
  </div>

  <div class="provenance">
    <strong>Evidence Summary:</strong> {prov.get('evidence_summary', 'N/A')}<br><br>
    <strong>Primary Source:</strong> {src.get('citation', 'N/A')}
    {f' — <a href="https://doi.org/{src["doi"]}" target="_blank">DOI</a>' if src.get('doi') else ''}<br>
    <strong>PMIDs:</strong> {pmid_links}<br>
    <strong>Clinical Trials:</strong> {trial_links or 'N/A'}
  </div>
</div>
"""


def _mage_network_section(mage_entries: list[dict]) -> str:
    """Table of all MAGE-A3 related decoy entries."""
    if not mage_entries:
        return ""

    rows = ""
    for e in mage_entries:
        pi = e["peptide_info"]
        rp = e.get("risk_profile", {})
        ee = e.get("experimental_evidence", {})
        vf = e.get("validation_flags", {})
        ev = rp.get("evidence_level", "Unknown")

        ev_color = {
            "Level_1_Clinical_Fatal": "#c0392b",
            "Level_2_In_Vitro_Confirmed": "#e74c3c",
            "Level_3_High_Throughput_Screened": "#f39c12",
            "Level_4_In_Silico_High_Risk": "#3498db",
        }.get(ev, "#bdc3c7")

        ev_short = ev.replace("Level_", "L").split("_")[0] + ev.split("_")[0][-1:]
        ev_label = {"Level_1_Clinical_Fatal": "L1 Fatal",
                    "Level_2_In_Vitro_Confirmed": "L2 In Vitro",
                    "Level_3_High_Throughput_Screened": "L3 HT Screen",
                    "Level_4_In_Silico_High_Risk": "L4 In Silico"}.get(ev, "Unknown")

        ms = "Yes" if ee.get("mass_spec_confirmed") else "No"
        status = vf.get("overall_status", "?")
        organs = ", ".join(rp.get("critical_organs_affected", [])[:3])

        rows += f"""<tr>
  <td><strong>{e['decoy_id']}</strong></td>
  <td style="font-family:Consolas; letter-spacing:2px">{pi['decoy_sequence']}</td>
  <td>{pi['hla_allele']}</td>
  <td>{pi['gene_symbol']}</td>
  <td><span style="color:{ev_color};font-weight:bold">{ev_label}</span></td>
  <td>{organs}</td>
  <td>{ms}</td>
  <td>{status}</td>
</tr>"""

    return f"""
<div class="section-card">
  <h2>MAGE-A3 Cross-Reactivity Network ({len(mage_entries)} entries)</h2>
  <p>All decoy peptides discovered in context of MAGE-A3 targeted TCR therapies:</p>
  <table>
    <tr>
      <th>ID</th><th>Peptide</th><th>HLA</th><th>Gene</th>
      <th>Evidence</th><th>Organs</th><th>Mass-Spec</th><th>Status</th>
    </tr>
    {rows}
  </table>
</div>
"""


def _evidence_sunburst(entries: list[dict]) -> str:
    """Plotly sunburst: Evidence Level → Expression Pattern → Validation."""
    import plotly.graph_objects as go

    labels, parents, values, colors = [], [], [], []

    # Root
    labels.append("All Entries")
    parents.append("")
    values.append(len(entries))
    colors.append("#2c3e50")

    # Level 1: Evidence levels
    ev_map = {
        "Level_1_Clinical_Fatal": ("L1: Fatal", "#c0392b"),
        "Level_2_In_Vitro_Confirmed": ("L2: In Vitro", "#e74c3c"),
        "Level_3_High_Throughput_Screened": ("L3: HT Screen", "#f39c12"),
        "Level_4_In_Silico_High_Risk": ("L4: In Silico", "#3498db"),
        "Unknown": ("Unknown", "#bdc3c7"),
    }

    from collections import defaultdict
    ev_groups = defaultdict(list)
    for e in entries:
        ev = e.get("risk_profile", {}).get("evidence_level", "Unknown")
        ev_groups[ev].append(e)

    for ev_key, group in ev_groups.items():
        ev_label, ev_color = ev_map.get(ev_key, (ev_key, "#999"))
        labels.append(ev_label)
        parents.append("All Entries")
        values.append(len(group))
        colors.append(ev_color)

        # Level 2: Validation status within each evidence level
        val_groups = defaultdict(int)
        for e in group:
            vs = e.get("validation_flags", {}).get("overall_status", "Unknown")
            val_groups[vs] += 1

        val_color_map = {"VALIDATED": "#27ae60", "PARTIAL": "#f39c12",
                         "NEEDS_REVIEW": "#e74c3c", "Unknown": "#bdc3c7"}
        for vs, count in val_groups.items():
            child_label = f"{ev_label} / {vs}"
            labels.append(child_label)
            parents.append(ev_label)
            values.append(count)
            colors.append(val_color_map.get(vs, "#999"))

    fig = go.Figure(go.Sunburst(
        labels=labels, parents=parents, values=values,
        marker=dict(colors=colors),
        branchvalues="total",
        hovertemplate="<b>%{label}</b><br>Count: %{value}<extra></extra>",
        textinfo="label+value",
        insidetextorientation="radial",
    ))
    fig.update_layout(
        title=dict(text="Evidence Level → Validation Status (Sunburst)",
                   font_size=14),
        width=650, height=550,
        margin=dict(t=50, l=10, r=10, b=10),
    )

    plot_json = fig.to_json()
    return f"""
<div class="section-card">
  <h2>Library Composition — Evidence & Validation</h2>
  <p>Interactive sunburst: click to drill down from evidence level to validation status.</p>
  <div id="sunburst_plot" style="margin: 0 auto; width: 650px;"></div>
  <script>
  (function() {{
    var spec = {plot_json};
    Plotly.newPlot('sunburst_plot', spec.data, spec.layout);
  }})();
  </script>
</div>
"""


# =========================================================================
# Main
# =========================================================================

def main():
    print("=" * 60)
    print("  Decoy C Library — Visualization Suite")
    print("=" * 60)

    print("\nLoading library...")
    entries = load_library()
    print(f"  {len(entries)} entries loaded")

    # Part 1: Matplotlib
    print("\n-- Part 1: Matplotlib Figures --")
    plot_pipeline_flowchart()
    plot_library_overview(entries)
    plot_filtering_sankey(entries)
    plot_assay_and_year(entries)

    # Part 2: Plotly HTML
    print("\n-- Part 2: Plotly Interactive Case Study --")
    generate_interactive_html(entries)

    print(f"\nAll figures saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
