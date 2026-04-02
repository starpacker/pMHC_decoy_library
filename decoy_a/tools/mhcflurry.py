"""
MHCflurry Wrapper — Local MHC-I Binding Prediction
===================================================
Pure-Python alternative to NetMHCpan. No license required.
Uses deep learning models trained on mass spec eluted ligand data.

Setup:
    pip install mhcflurry
    mhcflurry-downloads fetch models_class1_presentation
    mhcflurry-downloads fetch models_class1_pan

Public API:
    check_available()         -> bool
    predict_binding(peptides, hla) -> list[MHCflurryResult]
    predict_binding_batch(peptides, hla, batch_size) -> list[MHCflurryResult]
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import List, Optional

log = logging.getLogger(__name__)

# Ensure UTF-8 on Windows (mhcgnomes YAML parser requirement)
os.environ.setdefault("PYTHONUTF8", "1")

# Force mhcflurry data to /share/liuyutian (not home directory)
os.environ.setdefault("MHCFLURRY_DATA_DIR", "/share/liuyutian/mhcflurry_data/4")

# Thresholds (same as NetMHCpan convention)
EL_RANK_STRONG = 0.5
EL_RANK_WEAK = 2.0

# Lazy-loaded predictor singleton
_predictor = None


@dataclass
class MHCflurryResult:
    """Single peptide-HLA prediction result."""
    sequence: str
    hla_allele: str
    affinity_nm: float       # Predicted IC50 in nM (lower = stronger)
    affinity_percentile: float  # %Rank for affinity (lower = stronger)
    processing_score: float  # Antigen processing score
    presentation_score: float  # Combined presentation score
    presentation_percentile: float  # %Rank for presentation

    @property
    def el_rank(self) -> float:
        """EL-equivalent %Rank (use presentation_percentile)."""
        return self.presentation_percentile

    @property
    def bind_level(self) -> str:
        if self.el_rank <= EL_RANK_STRONG:
            return "SB"
        elif self.el_rank <= EL_RANK_WEAK:
            return "WB"
        return ""


def check_available() -> bool:
    """Check if mhcflurry is installed and pan-allele models are downloaded."""
    try:
        import mhcflurry
        from mhcflurry import Class1AffinityPredictor
        # Try loading — will fail if models not downloaded
        Class1AffinityPredictor.load()
        return True
    except (ImportError, Exception):
        return False


def _get_predictor():
    """Get or create the singleton predictor instance (Class1AffinityPredictor)."""
    global _predictor
    if _predictor is None:
        try:
            from mhcflurry import Class1AffinityPredictor
            log.info("Loading mhcflurry Class1AffinityPredictor...")
            _predictor = Class1AffinityPredictor.load()
            log.info("mhcflurry predictor loaded")
        except Exception as e:
            log.error("Failed to load mhcflurry: %s", e)
            raise
    return _predictor


def _normalize_allele(hla_allele: str) -> str:
    """Normalize HLA allele to mhcflurry format: HLA-A*02:01 → HLA-A0201."""
    allele = hla_allele.strip()
    # Remove '*' and ':' for mhcflurry: HLA-A*02:01 → HLA-A0201
    allele = allele.replace("*", "").replace(":", "")
    return allele


def predict_binding(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
) -> List[MHCflurryResult]:
    """
    Predict peptide-MHC-I binding for a list of peptides.

    Parameters
    ----------
    peptides : list[str]
        Peptide sequences (8-15 AA).
    hla_allele : str
        HLA allele (e.g., HLA-A*02:01).

    Returns
    -------
    list[MHCflurryResult]
        One result per peptide.
    """
    if not peptides:
        return []

    predictor = _get_predictor()
    allele = _normalize_allele(hla_allele)

    alleles = [allele] * len(peptides)

    try:
        df = predictor.predict_to_dataframe(
            peptides=peptides,
            allele=allele,
        )

        results = []
        for _, row in df.iterrows():
            pct = float(row.get("prediction_percentile", 100))
            aff = float(row.get("prediction", 50000))
            results.append(MHCflurryResult(
                sequence=str(row.get("peptide", "")),
                hla_allele=hla_allele,
                affinity_nm=aff,
                affinity_percentile=pct,
                processing_score=0.0,
                presentation_score=0.0,
                presentation_percentile=pct,  # Use affinity %Rank as EL proxy
            ))
        return results

    except Exception as e:
        log.error("mhcflurry prediction failed: %s", e)
        return []


def predict_binding_batch(
    peptides: List[str],
    hla_allele: str = "HLA-A*02:01",
    batch_size: int = 100000,
) -> List[MHCflurryResult]:
    """
    Predict binding in batches for large peptide sets.

    Parameters
    ----------
    peptides : list[str]
        All peptide sequences.
    hla_allele : str
        Target HLA allele.
    batch_size : int
        Peptides per batch.

    Returns
    -------
    list[MHCflurryResult]
    """
    all_results = []
    total = len(peptides)

    for start in range(0, total, batch_size):
        batch = peptides[start:start + batch_size]
        batch_num = start // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size

        log.info(
            "mhcflurry batch %d/%d: %d peptides (%s)",
            batch_num, total_batches, len(batch), hla_allele,
        )

        results = predict_binding(batch, hla_allele)
        all_results.extend(results)

        n_binders = sum(1 for r in results if r.el_rank <= EL_RANK_WEAK)
        log.info(
            "  Batch %d: %d/%d binders (SB: %d, WB: %d)",
            batch_num, n_binders, len(results),
            sum(1 for r in results if r.el_rank <= EL_RANK_STRONG),
            sum(1 for r in results if EL_RANK_STRONG < r.el_rank <= EL_RANK_WEAK),
        )

    return all_results


def filter_binders(
    results: List[MHCflurryResult],
    rank_threshold: float = EL_RANK_WEAK,
) -> List[MHCflurryResult]:
    """Keep only peptides with presentation %Rank <= threshold."""
    return [r for r in results if r.el_rank <= rank_threshold]
