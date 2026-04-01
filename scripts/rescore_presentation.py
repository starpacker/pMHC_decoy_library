"""
Re-score HLA-filtered binders with mhcflurry Presentation Score Model.

The original filtering used Class1AffinityPredictor (binding affinity only).
This script adds presentation_score and processing_score columns using
Class1PresentationPredictor, which integrates:
  - Binding affinity (MHC-peptide interaction)
  - Antigen processing (proteasomal cleavage + TAP transport)

Usage:
    PYTHONUTF8=1 python scripts/rescore_presentation.py
"""

import logging
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

os.environ.setdefault("PYTHONUTF8", "1")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "decoy_a"
INPUT_PATH = DATA_DIR / "hla_filtered_HLA-A0201.parquet"
OUTPUT_PATH = DATA_DIR / "hla_filtered_HLA-A0201_presentation.parquet"

HLA_ALLELE = "HLA-A0201"
BATCH_SIZE = 50_000


def main():
    log.info("=" * 60)
    log.info("  Presentation Score Re-scoring")
    log.info("=" * 60)

    # Load existing binders
    log.info("Loading %s ...", INPUT_PATH.name)
    df = pd.read_parquet(INPUT_PATH)
    log.info("  %d binders loaded", len(df))
    log.info("  Existing columns: %s", list(df.columns))

    # Load presentation predictor
    log.info("Loading mhcflurry Class1PresentationPredictor ...")
    from mhcflurry import Class1PresentationPredictor
    predictor = Class1PresentationPredictor.load()
    log.info("  Predictor loaded (supports %d alleles)", len(predictor.supported_alleles))

    peptides = df["sequence"].tolist()
    total = len(peptides)
    n_batches = (total + BATCH_SIZE - 1) // BATCH_SIZE

    # Pre-allocate result arrays
    processing_scores = np.full(total, np.nan, dtype=np.float32)
    presentation_scores = np.full(total, np.nan, dtype=np.float32)
    presentation_pctls = np.full(total, np.nan, dtype=np.float32)
    affinities = np.full(total, np.nan, dtype=np.float32)

    t_start = time.time()

    for batch_idx in range(n_batches):
        start = batch_idx * BATCH_SIZE
        end = min(start + BATCH_SIZE, total)
        batch_peps = peptides[start:end]

        t0 = time.time()
        result_df = predictor.predict(
            peptides=batch_peps,
            alleles=[HLA_ALLELE],
            verbose=0,
        )
        elapsed = time.time() - t0

        # Store results
        processing_scores[start:end] = result_df["processing_score"].values.astype(np.float32)
        presentation_scores[start:end] = result_df["presentation_score"].values.astype(np.float32)
        presentation_pctls[start:end] = result_df["presentation_percentile"].values.astype(np.float32)
        affinities[start:end] = result_df["affinity"].values.astype(np.float32)

        speed = len(batch_peps) / elapsed
        eta = (total - end) / speed if speed > 0 else 0
        log.info(
            "  Batch %d/%d: %d-%d (%d peps) in %.1fs (%.0f pep/s, ETA %.1f min)",
            batch_idx + 1, n_batches, start, end, len(batch_peps),
            elapsed, speed, eta / 60,
        )

    total_time = time.time() - t_start
    log.info("Total scoring time: %.1f min (%.0f pep/s avg)",
             total_time / 60, total / total_time)

    # Add new columns to dataframe
    df["processing_score"] = processing_scores
    df["presentation_score"] = presentation_scores
    df["presentation_percentile"] = presentation_pctls
    df["affinity_nM"] = affinities

    # Reclassify binding using presentation_percentile instead of affinity
    def classify_presentation(row):
        pctl = row["presentation_percentile"]
        if pctl <= 0.5:
            return "Strong_Binder"
        elif pctl <= 2.0:
            return "Weak_Binder"
        else:
            return "Non_Binder"

    df["presentation_binding"] = df.apply(classify_presentation, axis=1)

    # Stats
    log.info("")
    log.info("=" * 60)
    log.info("  Results Summary")
    log.info("=" * 60)

    log.info("")
    log.info("Old classification (affinity-based el_rank):")
    log.info("  %s", df["binding"].value_counts().to_dict())

    log.info("")
    log.info("New classification (presentation_percentile):")
    log.info("  %s", df["presentation_binding"].value_counts().to_dict())

    # Cross-tabulation
    ct = pd.crosstab(df["binding"], df["presentation_binding"], margins=True)
    log.info("")
    log.info("Cross-tabulation (old vs new):")
    log.info("\n%s", ct.to_string())

    log.info("")
    log.info("Presentation score stats:")
    log.info("\n%s", df["presentation_score"].describe().to_string())

    log.info("")
    log.info("Processing score stats:")
    log.info("\n%s", df["processing_score"].describe().to_string())

    # How many binders get reclassified?
    upgraded = ((df["binding"] == "Weak_Binder") &
                (df["presentation_binding"] == "Strong_Binder")).sum()
    downgraded = ((df["binding"] == "Strong_Binder") &
                  (df["presentation_binding"] != "Strong_Binder")).sum()
    lost = (df["presentation_binding"] == "Non_Binder").sum()

    log.info("")
    log.info("Reclassification impact:")
    log.info("  Weak → Strong (upgraded): %d", upgraded)
    log.info("  Strong → Weak/Non (downgraded): %d", downgraded)
    log.info("  Now Non_Binder (lost from pool): %d", lost)
    log.info("  Remaining binders (presentation): %d",
             (df["presentation_binding"] != "Non_Binder").sum())

    # Save enriched file
    log.info("")
    log.info("Saving to %s ...", OUTPUT_PATH.name)
    df.to_parquet(OUTPUT_PATH, index=False)
    log.info("  Done. File size: %.1f MB", OUTPUT_PATH.stat().st_size / 1e6)


if __name__ == "__main__":
    main()
