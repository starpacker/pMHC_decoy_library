"""
Benchmark: Interface Descriptor Methods for Decoy B

Tests each descriptor (PLIP, FreeSASA, PRODIGY) individually and in combination
on available pMHC PDB structures. Compares:
  1. Per-method availability and runtime
  2. Per-method score distributions
  3. Effect on final surface_correlation ranking
  4. Old vs new formula ranking comparison

Usage (on server):
    cd /share/liuyutian/pMHC_decoy_library
    python scripts/benchmark_interface_descriptors.py

    # Or point to specific PDB directory:
    python scripts/benchmark_interface_descriptors.py --pdb-dir data/decoy_b/pmhc_models/tfold

Output:
    figures/benchmark_interface_descriptors.png
    figures/benchmark_report.txt
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
log = logging.getLogger("benchmark")

OUT = PROJECT_ROOT / "figures"
OUT.mkdir(exist_ok=True)


# ── Find PDB files ──────────────────────────────────────────────────────

def find_pdbs(pdb_dir: Optional[Path] = None) -> List[Path]:
    """Find pMHC PDB files for benchmarking."""
    search_dirs = []
    if pdb_dir:
        search_dirs.append(pdb_dir)
    else:
        data_b = PROJECT_ROOT / "data" / "decoy_b"
        search_dirs.extend([
            data_b / "pmhc_models" / "tfold",
            data_b / "pmhc_models" / "af3",
            data_b / "stepwise_demo" / "tfold_pdbs",
            data_b,
        ])
        # Also check temp location
        temp_data = Path(r"C:\Users\30670\AppData\Local\Temp\pMHC_decoy_library\data\decoy_b")
        if temp_data.exists():
            search_dirs.append(temp_data)
        # Server location
        server_data = Path("/share/liuyutian/pMHC_decoy_library/data/decoy_b")
        if server_data.exists():
            search_dirs.extend([
                server_data / "pmhc_models" / "tfold",
                server_data / "pmhc_models" / "af3",
                server_data,
            ])

    pdbs = []
    seen = set()
    for d in search_dirs:
        if d.exists():
            for p in sorted(d.rglob("*.pdb")):
                if p.name not in seen:
                    pdbs.append(p)
                    seen.add(p.name)
    return pdbs


# ── Benchmark each method ───────────────────────────────────────────────

def benchmark_plip(pdbs: List[Path]) -> Dict:
    """Benchmark PLIP on all PDBs."""
    from decoy_b.tools.interface_descriptors import compute_plip_fingerprint

    results = {"available": False, "timings": [], "scores": [], "failures": 0}

    for pdb in pdbs:
        t0 = time.time()
        fp = compute_plip_fingerprint(str(pdb))
        dt = time.time() - t0

        if fp is None:
            if not results["available"]:
                log.warning("PLIP not available")
                return results
            results["failures"] += 1
        else:
            results["available"] = True
            results["timings"].append(dt)
            total = fp.hbond_count + fp.hydrophobic_count + fp.salt_bridge_count + fp.pi_stacking_count + fp.pi_cation_count
            results["scores"].append({
                "pdb": pdb.name,
                "hbond": fp.hbond_count,
                "hydrophobic": fp.hydrophobic_count,
                "salt_bridge": fp.salt_bridge_count,
                "pi_stacking": fp.pi_stacking_count,
                "pi_cation": fp.pi_cation_count,
                "total_interactions": total,
            })

    return results


def benchmark_freesasa(pdbs: List[Path]) -> Dict:
    """Benchmark FreeSASA on all PDBs."""
    from decoy_b.tools.interface_descriptors import compute_bsa

    results = {"available": False, "timings": [], "scores": [], "failures": 0}

    for pdb in pdbs:
        t0 = time.time()
        bsa = compute_bsa(str(pdb))
        dt = time.time() - t0

        if bsa is None:
            if not results["available"]:
                log.warning("FreeSASA not available")
                return results
            results["failures"] += 1
        else:
            results["available"] = True
            results["timings"].append(dt)
            results["scores"].append({
                "pdb": pdb.name,
                "bsa": bsa,
            })

    return results


def benchmark_prodigy(pdbs: List[Path]) -> Dict:
    """Benchmark PRODIGY on all PDBs."""
    from decoy_b.tools.interface_descriptors import compute_prodigy_affinity

    results = {"available": False, "timings": [], "scores": [], "failures": 0}

    for pdb in pdbs:
        t0 = time.time()
        dg = compute_prodigy_affinity(str(pdb))
        dt = time.time() - t0

        if dg is None:
            if not results["available"]:
                log.warning("PRODIGY not available")
                return results
            results["failures"] += 1
        else:
            results["available"] = True
            results["timings"].append(dt)
            results["scores"].append({
                "pdb": pdb.name,
                "dg": dg,
            })

    return results


def benchmark_combined(pdbs: List[Path]) -> Dict:
    """Benchmark combined interface similarity (all pairs)."""
    from decoy_b.tools.interface_descriptors import compute_interface_similarity

    if len(pdbs) < 2:
        return {"available": False, "pairs": []}

    results = {"available": False, "timings": [], "pairs": []}

    # Use first PDB as target, compare against all others
    target = pdbs[0]

    for cand in pdbs[1:]:
        t0 = time.time()
        sim = compute_interface_similarity(str(target), str(cand))
        dt = time.time() - t0

        results["available"] = True
        results["timings"].append(dt)
        results["pairs"].append({
            "target": target.name,
            "candidate": cand.name,
            "plip_tanimoto": sim.plip_tanimoto,
            "bsa_similarity": sim.bsa_similarity,
            "prodigy_similarity": sim.prodigy_similarity,
            "combined": sim.combined,
        })

    return results


def benchmark_old_vs_new(pdbs: List[Path]) -> Dict:
    """Compare old scoring (RMSD-only) vs new (RMSD + interface)."""
    from decoy_b.scanner import compute_structure_similarity

    if len(pdbs) < 2:
        return {"comparisons": []}

    target = pdbs[0]
    comparisons = []

    for cand in pdbs[1:]:
        result = compute_structure_similarity(str(target), str(cand))
        comparisons.append({
            "target": target.name,
            "candidate": cand.name,
            "surface_correlation": result.surface_correlation,
            "rmsd": result.rmsd,
            "plip_tanimoto": result.plip_tanimoto,
            "bsa_similarity": result.bsa_similarity,
            "prodigy_similarity": result.prodigy_similarity,
            "interface_combined": result.interface_combined,
            "modeling_tool": result.modeling_tool,
        })

    return {"comparisons": comparisons}


# ── Report ──────────────────────────────────────────────────────────────

def generate_report(
    pdbs: List[Path],
    plip_res: Dict,
    sasa_res: Dict,
    prodigy_res: Dict,
    combined_res: Dict,
    old_vs_new: Dict,
) -> str:
    """Generate a text report."""
    lines = []
    lines.append("=" * 70)
    lines.append("  Interface Descriptor Benchmark Report")
    lines.append("=" * 70)
    lines.append(f"\nTotal PDB files found: {len(pdbs)}")
    if pdbs:
        lines.append(f"PDB directory: {pdbs[0].parent}")
        lines.append(f"Sample PDBs: {', '.join(p.name for p in pdbs[:5])}")

    # PLIP
    lines.append(f"\n{'─' * 50}")
    lines.append("1. PLIP (Non-Covalent Interaction Fingerprint)")
    lines.append(f"{'─' * 50}")
    if plip_res["available"]:
        lines.append(f"  Status: AVAILABLE")
        lines.append(f"  Structures analyzed: {len(plip_res['scores'])}")
        lines.append(f"  Failures: {plip_res['failures']}")
        if plip_res["timings"]:
            lines.append(f"  Avg time/structure: {np.mean(plip_res['timings']):.3f}s")
        if plip_res["scores"]:
            totals = [s["total_interactions"] for s in plip_res["scores"]]
            lines.append(f"  Interactions per structure: mean={np.mean(totals):.1f}, "
                        f"min={min(totals)}, max={max(totals)}")
            hbonds = [s["hbond"] for s in plip_res["scores"]]
            hydro = [s["hydrophobic"] for s in plip_res["scores"]]
            lines.append(f"    H-bonds: mean={np.mean(hbonds):.1f}")
            lines.append(f"    Hydrophobic: mean={np.mean(hydro):.1f}")
            lines.append(f"    Salt bridges: mean={np.mean([s['salt_bridge'] for s in plip_res['scores']]):.1f}")
    else:
        lines.append("  Status: NOT AVAILABLE (install: conda install openbabel && pip install plip)")

    # FreeSASA
    lines.append(f"\n{'─' * 50}")
    lines.append("2. FreeSASA (Buried Surface Area)")
    lines.append(f"{'─' * 50}")
    if sasa_res["available"]:
        lines.append(f"  Status: AVAILABLE")
        lines.append(f"  Structures analyzed: {len(sasa_res['scores'])}")
        lines.append(f"  Failures: {sasa_res['failures']}")
        if sasa_res["timings"]:
            lines.append(f"  Avg time/structure: {np.mean(sasa_res['timings']):.3f}s")
        if sasa_res["scores"]:
            bsas = [s["bsa"] for s in sasa_res["scores"]]
            lines.append(f"  BSA range: {min(bsas):.1f} - {max(bsas):.1f} A^2 "
                        f"(mean={np.mean(bsas):.1f})")
    else:
        lines.append("  Status: NOT AVAILABLE (install: pip install freesasa)")

    # PRODIGY
    lines.append(f"\n{'─' * 50}")
    lines.append("3. PRODIGY (Binding Affinity Prediction)")
    lines.append(f"{'─' * 50}")
    if prodigy_res["available"]:
        lines.append(f"  Status: AVAILABLE")
        lines.append(f"  Structures analyzed: {len(prodigy_res['scores'])}")
        lines.append(f"  Failures: {prodigy_res['failures']}")
        if prodigy_res["timings"]:
            lines.append(f"  Avg time/structure: {np.mean(prodigy_res['timings']):.3f}s")
        if prodigy_res["scores"]:
            dgs = [s["dg"] for s in prodigy_res["scores"]]
            lines.append(f"  dG range: {min(dgs):.2f} to {max(dgs):.2f} kcal/mol "
                        f"(mean={np.mean(dgs):.2f})")
    else:
        lines.append("  Status: NOT AVAILABLE (install: pip install prodigy-prot)")

    # Combined similarity
    lines.append(f"\n{'─' * 50}")
    lines.append("4. Combined Interface Similarity (Pairwise)")
    lines.append(f"{'─' * 50}")
    if combined_res["available"] and combined_res["pairs"]:
        lines.append(f"  Pairs compared: {len(combined_res['pairs'])}")
        if combined_res["timings"]:
            lines.append(f"  Avg time/pair: {np.mean(combined_res['timings']):.3f}s")
        combos = [p["combined"] for p in combined_res["pairs"]]
        lines.append(f"  Combined similarity: mean={np.mean(combos):.3f}, "
                    f"min={min(combos):.3f}, max={max(combos):.3f}")
        lines.append("\n  Per-pair breakdown:")
        lines.append(f"  {'Target':<25} {'Candidate':<25} {'PLIP':>6} {'BSA':>6} {'PRODY':>6} {'Combo':>6}")
        for p in combined_res["pairs"]:
            lines.append(
                f"  {p['target']:<25} {p['candidate']:<25} "
                f"{p['plip_tanimoto']:6.3f} {p['bsa_similarity']:6.3f} "
                f"{p['prodigy_similarity']:6.3f} {p['combined']:6.3f}"
            )
    else:
        lines.append("  No pairs available for comparison")

    # Old vs New scoring
    lines.append(f"\n{'─' * 50}")
    lines.append("5. Old vs New Scoring Comparison")
    lines.append(f"{'─' * 50}")
    if old_vs_new["comparisons"]:
        lines.append(f"  {'Candidate':<25} {'RMSD':>6} {'OldCorr':>8} {'NewCorr':>8} {'IntComb':>8} {'Tool'}")
        for c in old_vs_new["comparisons"]:
            rmsd_str = f"{c['rmsd']:.2f}" if c['rmsd'] is not None else "N/A"
            ic_str = f"{c['interface_combined']:.3f}" if c['interface_combined'] is not None else "N/A"
            lines.append(
                f"  {c['candidate']:<25} {rmsd_str:>6} "
                f"{c['surface_correlation']:8.4f}   ---     {ic_str:>8} {c['modeling_tool']}"
            )
    else:
        lines.append("  No comparisons available")

    # Summary
    lines.append(f"\n{'=' * 70}")
    lines.append("  Summary")
    lines.append(f"{'=' * 70}")
    n_avail = sum(1 for r in [plip_res, sasa_res, prodigy_res] if r["available"])
    lines.append(f"  Methods available: {n_avail}/3")
    methods = []
    if plip_res["available"]: methods.append("PLIP")
    if sasa_res["available"]: methods.append("FreeSASA")
    if prodigy_res["available"]: methods.append("PRODIGY")
    lines.append(f"  Active methods: {', '.join(methods) if methods else 'None'}")
    if not methods:
        lines.append("\n  ACTION REQUIRED: Install dependencies on server:")
        lines.append("    conda install -c conda-forge openbabel -y")
        lines.append("    pip install plip freesasa prodigy-prot")

    return "\n".join(lines)


# ── Visualization ───────────────────────────────────────────────────────

def plot_benchmark(
    plip_res: Dict,
    sasa_res: Dict,
    prodigy_res: Dict,
    combined_res: Dict,
):
    """Generate benchmark visualization."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        log.warning("matplotlib not available, skipping visualization")
        return

    fig = plt.figure(figsize=(18, 12))
    fig.suptitle("Interface Descriptor Benchmark", fontsize=16, fontweight="bold", y=0.98)
    gs = gridspec.GridSpec(2, 3, hspace=0.35, wspace=0.35,
                           top=0.92, bottom=0.08, left=0.07, right=0.95)

    C_PLIP = "#E91E63"
    C_BSA = "#2196F3"
    C_PROD = "#FF9800"
    C_COMB = "#4CAF50"

    plt.rcParams.update({
        "font.size": 10,
        "axes.facecolor": "#F5F5F5",
        "figure.facecolor": "white",
        "axes.grid": True,
        "grid.alpha": 0.3,
    })

    # 1. Method availability & timing
    ax1 = fig.add_subplot(gs[0, 0])
    methods = ["PLIP", "FreeSASA", "PRODIGY"]
    avail = [plip_res["available"], sasa_res["available"], prodigy_res["available"]]
    colors = [C_PLIP, C_BSA, C_PROD]
    timings = []
    for res in [plip_res, sasa_res, prodigy_res]:
        if res["timings"]:
            timings.append(np.mean(res["timings"]))
        else:
            timings.append(0)

    bars = ax1.bar(methods, timings, color=colors, edgecolor="white", linewidth=1.5)
    for bar, a, t in zip(bars, avail, timings):
        label = f"{t:.3f}s" if a else "N/A"
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                 label, ha="center", va="bottom", fontweight="bold")
    ax1.set_ylabel("Avg Time per Structure (s)")
    ax1.set_title("Method Runtime", fontweight="bold")

    # 2. PLIP interaction type distribution
    ax2 = fig.add_subplot(gs[0, 1])
    if plip_res["available"] and plip_res["scores"]:
        types = ["hbond", "hydrophobic", "salt_bridge", "pi_stacking", "pi_cation"]
        type_means = [np.mean([s[t] for s in plip_res["scores"]]) for t in types]
        ax2.barh(types, type_means, color=C_PLIP, edgecolor="white")
        ax2.set_xlabel("Mean Count per Structure")
        ax2.set_title("PLIP: Interaction Types", fontweight="bold")
    else:
        ax2.text(0.5, 0.5, "PLIP not available", ha="center", va="center", transform=ax2.transAxes)
        ax2.set_title("PLIP: Interaction Types", fontweight="bold")

    # 3. BSA distribution
    ax3 = fig.add_subplot(gs[0, 2])
    if sasa_res["available"] and sasa_res["scores"]:
        bsas = [s["bsa"] for s in sasa_res["scores"]]
        ax3.hist(bsas, bins=15, color=C_BSA, edgecolor="white", alpha=0.85)
        ax3.axvline(np.mean(bsas), color="darkblue", linestyle="--",
                    label=f"Mean: {np.mean(bsas):.0f} A^2")
        ax3.set_xlabel("BSA (A^2)")
        ax3.set_ylabel("Count")
        ax3.legend(fontsize=9)
        ax3.set_title("FreeSASA: BSA Distribution", fontweight="bold")
    else:
        ax3.text(0.5, 0.5, "FreeSASA not available", ha="center", va="center", transform=ax3.transAxes)
        ax3.set_title("FreeSASA: BSA Distribution", fontweight="bold")

    # 4. PRODIGY dG distribution
    ax4 = fig.add_subplot(gs[1, 0])
    if prodigy_res["available"] and prodigy_res["scores"]:
        dgs = [s["dg"] for s in prodigy_res["scores"]]
        ax4.hist(dgs, bins=15, color=C_PROD, edgecolor="white", alpha=0.85)
        ax4.axvline(np.mean(dgs), color="darkorange", linestyle="--",
                    label=f"Mean: {np.mean(dgs):.2f} kcal/mol")
        ax4.set_xlabel("dG (kcal/mol)")
        ax4.set_ylabel("Count")
        ax4.legend(fontsize=9)
        ax4.set_title("PRODIGY: Binding Affinity", fontweight="bold")
    else:
        ax4.text(0.5, 0.5, "PRODIGY not available", ha="center", va="center", transform=ax4.transAxes)
        ax4.set_title("PRODIGY: Binding Affinity", fontweight="bold")

    # 5. Combined similarity scatter
    ax5 = fig.add_subplot(gs[1, 1])
    if combined_res["available"] and combined_res["pairs"]:
        plip_vals = [p["plip_tanimoto"] for p in combined_res["pairs"]]
        bsa_vals = [p["bsa_similarity"] for p in combined_res["pairs"]]
        prod_vals = [p["prodigy_similarity"] for p in combined_res["pairs"]]
        combo_vals = [p["combined"] for p in combined_res["pairs"]]

        x = range(len(combined_res["pairs"]))
        ax5.scatter(x, plip_vals, c=C_PLIP, label="PLIP", s=40, zorder=3)
        ax5.scatter(x, bsa_vals, c=C_BSA, label="BSA", s=40, zorder=3)
        ax5.scatter(x, prod_vals, c=C_PROD, label="PRODIGY", s=40, zorder=3)
        ax5.plot(x, combo_vals, c=C_COMB, linewidth=2, label="Combined", zorder=4)
        ax5.set_xlabel("Candidate Index")
        ax5.set_ylabel("Similarity [0-1]")
        ax5.legend(fontsize=8, loc="lower left")
        ax5.set_title("Per-Pair Similarity Scores", fontweight="bold")
    else:
        ax5.text(0.5, 0.5, "No pairs available", ha="center", va="center", transform=ax5.transAxes)
        ax5.set_title("Per-Pair Similarity Scores", fontweight="bold")

    # 6. Radar chart: descriptor contribution
    ax6 = fig.add_subplot(gs[1, 2], polar=True)
    if combined_res["available"] and combined_res["pairs"]:
        categories = ["PLIP\nTanimoto", "BSA\nSimilarity", "PRODIGY\nSimilarity"]
        N = len(categories)
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]

        # Mean values
        mean_vals = [
            np.mean([p["plip_tanimoto"] for p in combined_res["pairs"]]),
            np.mean([p["bsa_similarity"] for p in combined_res["pairs"]]),
            np.mean([p["prodigy_similarity"] for p in combined_res["pairs"]]),
        ]
        mean_vals += mean_vals[:1]

        ax6.plot(angles, mean_vals, 'o-', linewidth=2, color=C_COMB)
        ax6.fill(angles, mean_vals, alpha=0.25, color=C_COMB)
        ax6.set_xticks(angles[:-1])
        ax6.set_xticklabels(categories, fontsize=9)
        ax6.set_ylim(0, 1)
        ax6.set_title("Mean Descriptor Profile", fontweight="bold", pad=20)
    else:
        ax6.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax6.transAxes)

    plt.savefig(OUT / "benchmark_interface_descriptors.png", dpi=150, bbox_inches="tight")
    plt.close()
    log.info("Saved: %s", OUT / "benchmark_interface_descriptors.png")


# ── Main ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Benchmark interface descriptors")
    parser.add_argument("--pdb-dir", type=Path, default=None, help="Directory with PDB files")
    parser.add_argument("--max-pdbs", type=int, default=20, help="Max PDBs to benchmark")
    args = parser.parse_args()

    pdbs = find_pdbs(args.pdb_dir)
    if not pdbs:
        log.error("No PDB files found. Run Decoy B pipeline first or specify --pdb-dir")
        sys.exit(1)

    pdbs = pdbs[:args.max_pdbs]
    log.info("Found %d PDB files for benchmarking", len(pdbs))

    # Run benchmarks
    log.info("=" * 50)
    log.info("Benchmarking PLIP...")
    plip_res = benchmark_plip(pdbs)

    log.info("Benchmarking FreeSASA...")
    sasa_res = benchmark_freesasa(pdbs)

    log.info("Benchmarking PRODIGY...")
    prodigy_res = benchmark_prodigy(pdbs)

    log.info("Benchmarking combined similarity...")
    combined_res = benchmark_combined(pdbs)

    log.info("Comparing old vs new scoring...")
    old_vs_new = benchmark_old_vs_new(pdbs)

    # Generate report
    report = generate_report(pdbs, plip_res, sasa_res, prodigy_res, combined_res, old_vs_new)
    print(report)

    report_path = OUT / "benchmark_report.txt"
    report_path.write_text(report, encoding="utf-8")
    log.info("Report saved to: %s", report_path)

    # Generate visualization
    plot_benchmark(plip_res, sasa_res, prodigy_res, combined_res)

    # Save raw results as JSON for further analysis
    raw = {
        "n_pdbs": len(pdbs),
        "pdb_names": [p.name for p in pdbs],
        "plip": {"available": plip_res["available"], "n_scores": len(plip_res["scores"]), "failures": plip_res["failures"]},
        "freesasa": {"available": sasa_res["available"], "n_scores": len(sasa_res["scores"]), "failures": sasa_res["failures"]},
        "prodigy": {"available": prodigy_res["available"], "n_scores": len(prodigy_res["scores"]), "failures": prodigy_res["failures"]},
        "combined_pairs": combined_res.get("pairs", []),
        "old_vs_new": old_vs_new.get("comparisons", []),
    }
    json_path = OUT / "benchmark_raw_results.json"
    json_path.write_text(json.dumps(raw, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info("Raw results saved to: %s", json_path)


if __name__ == "__main__":
    main()
