#!/usr/bin/env python3
"""
Decoy B — Step-by-Step Pipeline Demo
=====================================
Target: GILGFVFTL (Influenza M1 / HLA-A*02:01)

This script runs each stage of Decoy B individually, saves all intermediate
files, and prints exactly what was produced at every step.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Setup ──────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
os.chdir(PROJECT_ROOT)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("decoy_b_stepwise")

# ── Parameters ─────────────────────────────────────────────────────────
TARGET = "GILGFVFTL"
HLA = "HLA-A*02:01"
DEMO_DIR = PROJECT_ROOT / "data" / "decoy_b" / "stepwise_demo"
DEMO_DIR.mkdir(parents=True, exist_ok=True)

SEPARATOR = "=" * 80


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STAGE 0 — 输入数据概览                                              ║
# ╚══════════════════════════════════════════════════════════════════════╝
def stage0_overview():
    print(f"\n{SEPARATOR}")
    print("  STAGE 0 — 输入数据概览")
    print(SEPARATOR)
    print(f"  靶标肽段 (Target):  {TARGET}")
    print(f"  HLA 型别:           {HLA}")
    print(f"  肽段长度:           {len(TARGET)}")
    print(f"  TCR 接触核心 (p4-p8): {TARGET[3:8]}")
    print(f"  HLA 锚定位 (p2,p9):  {TARGET[1]}, {TARGET[8]}")
    print()

    from decoy_a.hla_filter import load_hla_filtered
    hla_df = load_hla_filtered(HLA)
    n_total = len(hla_df)
    lens = hla_df["sequence"].str.len()
    print(f"  HLA-A*02:01 可呈递肽段总数: {n_total:,}")
    for l in sorted(lens.unique()):
        n = int((lens == l).sum())
        print(f"    {l}-mer: {n:>10,}  ({n/n_total*100:5.1f}%)")
    print(f"  数据文件: data/decoy_a/hla_filtered_HLA-A0201_presentation.parquet")
    print(SEPARATOR)
    return hla_df


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STAGE 1 — Atchley 理化初筛 (CPU, 秒级)                              ║
# ╚══════════════════════════════════════════════════════════════════════╝
def stage1_atchley(hla_df):
    print(f"\n{SEPARATOR}")
    print("  STAGE 1 — Atchley 理化因子初筛")
    print(SEPARATOR)

    from decoy_b.scanner import (
        physicochemical_screen,
        _get_tcr_contact_residues,
        _sequence_to_atchley_vector,
    )
    from decoy_a.scanner import hamming_distance_vectorised

    t0 = time.time()

    # 1a. 统计各长度候选数
    ALLOWED_LENGTHS = {8, 9, 10, 11}
    all_cands = hla_df[hla_df["sequence"].str.len().isin(ALLOWED_LENGTHS)].copy()
    all_cands = all_cands[all_cands["sequence"] != TARGET]
    cand_lens = all_cands["sequence"].str.len()

    # Same-length: exclude Decoy A territory (HD <= 2)
    same_len_mask = cand_lens == len(TARGET)
    same_len_seqs = all_cands.loc[same_len_mask, "sequence"].values
    distances = hamming_distance_vectorised(TARGET, same_len_seqs)
    n_decoy_a_territory = int((distances <= 2).sum())
    n_same_len_b = int((distances >= 3).sum())
    n_diff_len = int((~same_len_mask).sum())

    print(f"\n  [1a] 候选池统计 (HLA-I 8-11mer)")
    print(f"       总候选 (8-11mer):       {len(all_cands):,}")
    for l in sorted(ALLOWED_LENGTHS):
        n = int((cand_lens == l).sum())
        tag = " ← 靶标同长度" if l == len(TARGET) else ""
        print(f"         {l}-mer: {n:>10,}{tag}")
    print(f"       同长度属于 Decoy A (HD≤2): {n_decoy_a_territory}")
    print(f"       进入 Decoy B 候选池:        {n_same_len_b + n_diff_len:,}")
    print(f"         同长度 (HD≥3): {n_same_len_b:,}")
    print(f"         异长度:        {n_diff_len:,}")

    # 1b. Atchley 向量化 + 余弦相似度
    physchem_df = physicochemical_screen(
        TARGET, hla_df,
        min_hamming=3,
        cosine_threshold=0.70,
        top_k=5000,
    )

    elapsed = time.time() - t0
    print(f"\n  [1b] Atchley 因子余弦相似度筛选")
    print(f"       TCR 接触核心: {_get_tcr_contact_residues(TARGET)}")
    print(f"       Atchley 向量维度: {len(_sequence_to_atchley_vector(_get_tcr_contact_residues(TARGET)))}")
    print(f"       余弦相似度阈值: ≥ 0.70")
    print(f"       通过阈值的候选数: {len(physchem_df)}")
    print(f"       耗时: {elapsed:.1f}s")

    # Save intermediate
    stage1_file = DEMO_DIR / "stage1_atchley_candidates.csv"
    save_cols = ["sequence", "cosine_similarity", "hamming_distance", "el_rank", "gene_symbols"]
    save_df = physchem_df[[c for c in save_cols if c in physchem_df.columns]].copy()
    save_df.to_csv(stage1_file, index=False)
    print(f"\n  ✅ 中间文件已保存: {stage1_file.relative_to(PROJECT_ROOT)}")

    # Show top 20
    print(f"\n  Stage 1 筛选结果 — Top 20 (共 {len(physchem_df)} 条):")
    print(f"  {'#':>3}  {'Sequence':<14} {'Len':>3} {'HD':>3} {'CosSim':>7} {'EL%Rank':>8} {'Genes':<20}")
    print(f"  {'-'*75}")
    for i, (_, row) in enumerate(physchem_df.head(20).iterrows()):
        seq = row["sequence"]
        seq_len = len(seq)
        hd = row.get("hamming_distance", seq_len)
        cos = row["cosine_similarity"]
        elr = row.get("el_rank", 99.0)
        genes = row.get("gene_symbols", [])
        if hasattr(genes, "tolist"):
            genes = genes.tolist()
        gene_str = ",".join(genes[:3]) if isinstance(genes, list) else str(genes)
        hd_str = str(int(hd)) if seq_len == len(TARGET) else f"×{seq_len}"
        print(f"  {i+1:>3}  {seq:<14} {seq_len:>3} {hd_str:>3} {cos:>7.3f} {elr:>8.2f} {gene_str:<20}")

    print(SEPARATOR)
    return physchem_df


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STAGE 2 — tFold 批量结构预测 (GPU, 分钟级)                           ║
# ╚══════════════════════════════════════════════════════════════════════╝
def stage2_tfold(physchem_df):
    print(f"\n{SEPARATOR}")
    print("  STAGE 2 — tFold 批量 pMHC 结构预测")
    print(SEPARATOR)

    from decoy_b.tools.tfold import check_available, predict_pmhc

    available = check_available()
    print(f"\n  tFold 可用性: {'✅ 可用' if available else '❌ 不可用'}")

    tfold_dir = DEMO_DIR / "tfold_pdbs"
    tfold_dir.mkdir(parents=True, exist_ok=True)

    # We predict the target + top 10 candidates for demo
    top_seqs = physchem_df.head(10)["sequence"].tolist()
    all_seqs = [TARGET] + top_seqs

    results = {}
    if available:
        print(f"  预测 {len(all_seqs)} 个 pMHC 结构 (1 靶标 + {len(top_seqs)} 候选)...")
        t0 = time.time()
        for seq in all_seqs:
            print(f"    → 预测 {seq} ...", end=" ", flush=True)
            try:
                res = predict_pmhc(seq, HLA, tfold_dir)
                results[seq] = res
                if res.success:
                    print(f"✅ {Path(res.pdb_path).name} (confidence={res.confidence:.2f})")
                else:
                    print(f"❌ {res.error}")
            except Exception as e:
                print(f"❌ Error: {e}")
                results[seq] = None
        elapsed = time.time() - t0
        print(f"\n  耗时: {elapsed:.1f}s")
    else:
        print("  ⚠️  tFold 不可用，跳过结构预测。")
        print("  （在实际部署中，此步骤会为每个候选生成一个 .pdb 文件）")

    # List generated files
    pdb_files = sorted(tfold_dir.glob("*.pdb"))
    print(f"\n  生成的 PDB 文件 ({len(pdb_files)} 个):")
    for f in pdb_files:
        size = f.stat().st_size
        print(f"    📄 {f.name}  ({size:,} bytes)")

    if not pdb_files:
        print("    （无 PDB 文件生成）")

    print(f"\n  ✅ PDB 文件目录: {tfold_dir.relative_to(PROJECT_ROOT)}")
    print(SEPARATOR)
    return results, tfold_dir


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STAGE 3 — 结构比对与评分                                            ║
# ╚══════════════════════════════════════════════════════════════════════╝
def stage3_structure_comparison(physchem_df, tfold_results, tfold_dir):
    print(f"\n{SEPARATOR}")
    print("  STAGE 3 — 结构比对与综合评分")
    print(SEPARATOR)

    from decoy_b.scanner import compute_structure_similarity

    target_pdb = None
    if tfold_results and TARGET in tfold_results:
        r = tfold_results[TARGET]
        if r and r.success:
            target_pdb = r.pdb_path

    comparison_results = []
    top_seqs = physchem_df.head(10)["sequence"].tolist()

    if target_pdb:
        print(f"\n  靶标 PDB: {Path(target_pdb).name}")
        print(f"  ⚠️  RMSD 计算: MHC 对齐后，同长度→全肽段 RMSD，异长度→TCR 接触核心 RMSD")
        print(f"\n  {'#':>3}  {'Sequence':<14} {'Len':>3} {'HD':>3} {'CosSim':>7} {'PepRMSD':>8} {'SurfCorr':>8} {'Combined':>9} {'Tool':<30}")
        print(f"  {'-'*100}")

        for i, seq in enumerate(top_seqs):
            cos_sim = float(physchem_df[physchem_df["sequence"] == seq]["cosine_similarity"].iloc[0])
            seq_len = len(seq)
            hd = int(physchem_df[physchem_df["sequence"] == seq]["hamming_distance"].iloc[0])

            cand_pdb = None
            if seq in tfold_results and tfold_results[seq] and tfold_results[seq].success:
                cand_pdb = tfold_results[seq].pdb_path

            score = compute_structure_similarity(target_pdb, cand_pdb)

            if score.surface_correlation > 0:
                combined = 0.4 * cos_sim + 0.6 * score.surface_correlation
            else:
                combined = cos_sim

            rmsd_str = f"{score.rmsd:.3f}" if score.rmsd is not None else "N/A"
            corr_str = f"{score.surface_correlation:.3f}"
            hd_str = str(hd) if seq_len == len(TARGET) else f"×{seq_len}"

            print(f"  {i+1:>3}  {seq:<14} {seq_len:>3} {hd_str:>3} {cos_sim:>7.3f} {rmsd_str:>8} {corr_str:>8} {combined:>9.3f} {score.modeling_tool:<30}")

            comparison_results.append({
                "sequence": seq,
                "hamming_distance": hd,
                "cosine_similarity": cos_sim,
                "rmsd": score.rmsd,
                "bfactor_correlation": score.surface_correlation,
                "combined_score": combined,
                "modeling_tool": score.modeling_tool,
                "pdb_path": str(cand_pdb) if cand_pdb else None,
            })
    else:
        print("\n  ⚠️  靶标 PDB 不可用，仅使用 Atchley 理化相似度作为最终评分")
        for i, seq in enumerate(top_seqs):
            cos_sim = float(physchem_df[physchem_df["sequence"] == seq]["cosine_similarity"].iloc[0])
            hd = hamming(TARGET, seq)
            comparison_results.append({
                "sequence": seq,
                "hamming_distance": hd,
                "cosine_similarity": cos_sim,
                "rmsd": None,
                "bfactor_correlation": None,
                "combined_score": cos_sim,
                "modeling_tool": "physchem_only",
                "pdb_path": None,
            })

    # Save
    stage3_file = DEMO_DIR / "stage3_structure_comparison.json"
    with open(stage3_file, "w") as f:
        json.dump(comparison_results, f, indent=2)
    print(f"\n  ✅ 中间文件已保存: {stage3_file.relative_to(PROJECT_ROOT)}")
    print(SEPARATOR)
    return comparison_results


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  STAGE 4 — 最终 Decoy B 输出汇总                                     ║
# ╚══════════════════════════════════════════════════════════════════════╝
def stage4_final_output(physchem_df, comparison_results):
    print(f"\n{SEPARATOR}")
    print("  STAGE 4 — 最终 Decoy B 输出汇总")
    print(SEPARATOR)

    from decoy_a.config import ATCHLEY_FACTORS
    from decoy_b.scanner import _get_tcr_contact_residues

    # Build final table
    final_entries = []
    for cr in comparison_results:
        seq = cr["sequence"]
        row = physchem_df[physchem_df["sequence"] == seq]
        genes = []
        if not row.empty:
            g = row.iloc[0].get("gene_symbols", [])
            if hasattr(g, "tolist"):
                g = g.tolist()
            genes = g if isinstance(g, list) else []

        entry = {
            "rank": len(final_entries) + 1,
            "sequence": seq,
            "target": TARGET,
            "hamming_distance": cr["hamming_distance"],
            "tcr_contact_target": _get_tcr_contact_residues(TARGET),
            "tcr_contact_decoy": _get_tcr_contact_residues(seq),
            "cosine_similarity": cr["cosine_similarity"],
            "rmsd": cr["rmsd"],
            "combined_score": cr["combined_score"],
            "el_rank": float(row.iloc[0].get("el_rank", 99.0)) if not row.empty else 99.0,
            "gene_symbols": genes,
            "hla_allele": HLA,
            "source": "decoy_b",
        }
        final_entries.append(entry)

    # Sort by combined score
    final_entries.sort(key=lambda e: -e["combined_score"])
    for i, e in enumerate(final_entries):
        e["rank"] = i + 1

    # Save final output
    final_file = DEMO_DIR / "final_decoy_b_results.json"
    with open(final_file, "w") as f:
        json.dump(final_entries, f, indent=2, ensure_ascii=False)

    print(f"\n  最终 Decoy B 结果 (靶标: {TARGET} / {HLA}):")
    print(f"\n  {'Rank':>4}  {'Sequence':<14} {'Len':>3} {'HD':>3} {'TCR核心(靶标)':>14} {'TCR核心(Decoy)':>14} {'CosSim':>7} {'Score':>7} {'EL%Rank':>8} {'Genes':<20}")
    print(f"  {'-'*110}")
    for e in final_entries:
        gene_str = ",".join(e["gene_symbols"][:3]) if e["gene_symbols"] else "?"
        seq = e["sequence"]
        seq_len = len(seq)
        hd = e["hamming_distance"]
        hd_str = str(hd) if seq_len == len(TARGET) else f"×{seq_len}"
        print(f"  {e['rank']:>4}  {seq:<14} {seq_len:>3} {hd_str:>3} {e['tcr_contact_target']:>14} {e['tcr_contact_decoy']:>14} {e['cosine_similarity']:>7.3f} {e['combined_score']:>7.3f} {e['el_rank']:>8.2f} {gene_str:<20}")

    print(f"\n  ✅ 最终结果已保存: {final_file.relative_to(PROJECT_ROOT)}")
    print(SEPARATOR)
    return final_entries


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  文件清单汇总                                                        ║
# ╚══════════════════════════════════════════════════════════════════════╝
def print_file_manifest():
    print(f"\n{SEPARATOR}")
    print("  📁 生成的全部文件清单")
    print(SEPARATOR)

    for root, dirs, files in os.walk(DEMO_DIR):
        level = len(Path(root).relative_to(DEMO_DIR).parts)
        indent = "  " + "  " * level
        folder = Path(root).name
        if level == 0:
            print(f"  {DEMO_DIR.relative_to(PROJECT_ROOT)}/")
        else:
            print(f"{indent}{folder}/")
        for f in sorted(files):
            fp = Path(root) / f
            size = fp.stat().st_size
            if size > 1024 * 1024:
                size_str = f"{size / 1024 / 1024:.1f} MB"
            elif size > 1024:
                size_str = f"{size / 1024:.1f} KB"
            else:
                size_str = f"{size} B"
            print(f"{indent}  📄 {f}  ({size_str})")

    print(SEPARATOR)


# ╔══════════════════════════════════════════════════════════════════════╗
# ║  MAIN                                                              ║
# ╚══════════════════════════════════════════════════════════════════════╝
if __name__ == "__main__":
    print("\n" + "█" * 80)
    print("  Decoy B — 逐步管线演示")
    print(f"  靶标: {TARGET} / {HLA}")
    print("█" * 80)

    hla_df = stage0_overview()
    physchem_df = stage1_atchley(hla_df)
    tfold_results, tfold_dir = stage2_tfold(physchem_df)
    comparison_results = stage3_structure_comparison(physchem_df, tfold_results, tfold_dir)
    final_entries = stage4_final_output(physchem_df, comparison_results)
    print_file_manifest()

    print(f"\n{'█' * 80}")
    print(f"  ✅ Decoy B 逐步演示完成！")
    print(f"  所有中间文件保存在: {DEMO_DIR.relative_to(PROJECT_ROOT)}")
    print(f"{'█' * 80}\n")
