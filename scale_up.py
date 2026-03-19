#!/usr/bin/env python3
"""
Scale-Up Generator for Decoy C Library
=======================================

Automatically discovers, extracts, validates and deduplicates decoy peptide
entries until the library reaches a user-specified target size N.

Usage
-----
    python scale_up.py --required_num 10000
    python scale_up.py --required_num 500 --no-validate --batch-size 50
    python scale_up.py --required_num 10000 --resume          # continue previous run
    python scale_up.py --required_num 10000 --strategy all    # use all query strategies

How it works
------------
1. Load (or cold-start) existing library from ``data/decoy_library.json``
2. Compute ``gap = required_num − len(library)``
3. Generate PubMed queries via multiple *expansion strategies*
4. For each batch of PMIDs:
   a. Fetch papers
   b. Extract decoy entries (LLM or rule-based)
   c. Validate against UniProt / IEDB
   d. Deduplicate by peptide sequence
   e. Save checkpoint
5. Repeat until gap ≤ 0 or all strategies exhausted
6. Print final statistics
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import random
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# ── Ensure the parent package is importable ─────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from decoy_library.config import LIBRARY_JSON, DATA_DIR
from decoy_library.fetcher import (
    COLD_START_QUERIES,
    PaperRecord,
    fetch_abstract,
    fetch_multiple,
    search_pubmed,
)
from decoy_library.models import DecoyEntry, DecoyLibrary
from decoy_library.orchestrator import (
    load_library,
    process_papers,
    run_cold_start,
    save_library,
)
from decoy_library.seed_data import build_seed_library
from decoy_library.validator import validate_batch

# Multi-source fetchers
try:
    from decoy_library.multi_source_fetcher import (
        search_arxiv,
        search_biorxiv,
        search_clinical_trials,
        search_semantic_scholar,
        search_europe_pmc,
        search_all_sources,
        PaperRecord as MultiSourcePaperRecord,
        ARXIV_QUERIES,
        BIORXIV_QUERIES,
        CLINICAL_TRIALS_QUERIES,
        SEMANTIC_SCHOLAR_QUERIES,
    )
    MULTI_SOURCE_AVAILABLE = True
except ImportError:
    MULTI_SOURCE_AVAILABLE = False

log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────
# 1.  Query Expansion Strategies
# ─────────────────────────────────────────────────────────────────────────

# Strategy A: Domain-specific keyword combos (broad immunotherapy safety)
STRATEGY_A_QUERIES = [
    # ── BLOCK 1: Clinical toxicity reports ──
    '"TCR" AND "cross-reactivity" AND ("fatal" OR "toxicity" OR "cardiac")',
    '"adoptive cell therapy" AND "off-target" AND "peptide"',
    '"engineered T cell" AND "cross-reactivity" AND "safety"',
    '"affinity-enhanced TCR" AND "off-target"',
    '"MAGE-A3" AND "TCR" AND ("Titin" OR "cross-reactive")',
    '"T cell therapy" AND "serious adverse event" AND "peptide"',
    '"adoptive immunotherapy" AND "death" AND "cross-reactivity"',
    '"TCR gene therapy" AND "neurotoxicity"',
    '"TCR" AND "on-target off-tumor" AND "normal tissue"',
    # ── BLOCK 2: Peptide-MHC binding and TCR specificity ──
    '"peptide-MHC" AND "TCR" AND "cross-reactivity"',
    '"T-cell receptor" AND "off-target" AND "immunotherapy"',
    '"TCR-T" AND ("adverse event" OR "toxicity" OR "safety")',
    '"CAR-T" AND "on-target off-tumor" AND "peptide"',
    '"pMHC multimer" AND "cross-reactive" AND "T cell"',
    '"peptide exchange" AND "HLA" AND "T cell recognition"',
    '"TCR degeneracy" AND "peptide"',
    '"polyspecific" AND "T cell receptor" AND "peptide"',
    '"TCR promiscuity" AND "peptide" AND "MHC"',
    # ── BLOCK 3: High-throughput screening ──
    '"X-scan" AND "TCR" AND "peptide"',
    '"alanine scanning" AND "TCR" AND "cross-reactivity"',
    '"yeast display" AND "pMHC" AND "TCR"',
    '"peptide library" AND "TCR" AND "screening"',
    '"combinatorial peptide library" AND "T cell"',
    '"positional scanning" AND "peptide" AND "T cell"',
    '"random peptide library" AND "CTL" AND "HLA"',
    '"phage display" AND "peptide" AND "MHC"',
    '"peptide microarray" AND "T cell" AND "HLA"',
    '"deep mutational scanning" AND "peptide" AND "MHC"',
    # ── BLOCK 4: Computational prediction ──
    '"ARDitox" OR ("EpiTox" AND "TCR" AND "safety")',
    '"neoantigen" AND "TCR" AND "cross-reactivity"',
    '"molecular mimicry" AND "TCR" AND "peptide"',
    '"NetMHCpan" AND "peptide" AND "binding"',
    '"MHCflurry" AND "peptide" AND "presentation"',
    '"pVACseq" AND "neoantigen" AND "peptide"',
    '"TCR-pMHC" AND "binding prediction" AND "deep learning"',
    '"peptide binding affinity" AND "HLA class I" AND "prediction"',
    '"immunogenicity prediction" AND "peptide" AND "CTL"',
    '"cross-reactivity prediction" AND "TCR" AND "computational"',
    # ── BLOCK 5: Cancer-testis antigens & shared antigens ──
    '"cancer-testis antigen" AND "TCR" AND "peptide"',
    '"NY-ESO-1" AND "TCR" AND "cross-reactivity"',
    '"MAGE" AND "TCR" AND "peptide" AND "HLA"',
    '"PRAME" AND "T cell" AND "HLA"',
    '"WT1" AND "TCR" AND "peptide"',
    '"SSX2" AND "CTL" AND "peptide" AND "HLA"',
    '"LAGE-1" AND "T cell" AND "epitope"',
    '"cancer testis" AND "epitope" AND "HLA-A*02"',
    '"shared tumor antigen" AND "T cell" AND "peptide"',
    # ── BLOCK 6: Self-antigen / autoimmunity / tissue cross-reactivity ──
    '"tumor-associated antigen" AND "TCR" AND "normal tissue"',
    '"autoimmune" AND "TCR" AND "cross-reactive" AND "peptide"',
    '"self-antigen" AND "T cell" AND "peptide" AND "toxicity"',
    '"thymic selection" AND "cross-reactive" AND "peptide"',
    '"autoimmune encephalomyelitis" AND "molecular mimicry" AND "peptide"',
    '"myelin" AND "T cell" AND "cross-reactivity" AND "peptide"',
    '"type 1 diabetes" AND "T cell" AND "peptide" AND "HLA"',
    '"celiac disease" AND "T cell" AND "peptide" AND "HLA-DQ"',
    '"rheumatoid arthritis" AND "T cell" AND "peptide" AND "HLA-DR"',
    '"vitiligo" AND "T cell" AND "melanocyte" AND "peptide"',
    '"graft versus host" AND "minor histocompatibility" AND "peptide"',
    # ── BLOCK 7: MHC class I specific — major alleles ──
    '"HLA-A*02:01" AND "peptide" AND "TCR" AND "cross-reactivity"',
    '"HLA-A*01:01" AND "peptide" AND "TCR" AND "toxicity"',
    '"HLA-B*07:02" AND "peptide" AND "T cell"',
    '"HLA-A*11:01" AND "peptide" AND "T cell"',
    '"HLA-A*24:02" AND "peptide" AND "T cell"',
    '"HLA-A*03:01" AND "peptide" AND "T cell" AND "epitope"',
    '"HLA-B*08:01" AND "peptide" AND "CTL"',
    '"HLA-B*35:01" AND "peptide" AND "T cell"',
    '"HLA-B*44:02" AND "peptide" AND "CTL"',
    '"HLA-C*07:02" AND "peptide" AND "T cell"',
    '"HLA-A*68:01" AND "peptide" AND "epitope"',
    '"HLA-B*15:01" AND "peptide" AND "T cell"',
    '"HLA-B*27:05" AND "peptide" AND "spondylitis"',
    '"HLA-B*57:01" AND "peptide" AND "HIV"',
    '"HLA-A*02:01" AND "9mer" AND "binding"',
    '"HLA-A*02:01" AND "epitope" AND "mass spectrometry"',
    # ── BLOCK 8: Specific proteins of interest ──
    '"gp100" AND "TCR" AND ("melanocyte" OR "cross-reactivity")',
    '"tyrosinase" AND "TCR" AND "HLA"',
    '"CEA" AND "TCR" AND "colitis"',
    '"mesothelin" AND "TCR" AND "T cell"',
    '"AFP" AND "TCR" AND "hepatocellular"',
    '"HER2" AND "T cell" AND "cardiac"',
    '"EGFR" AND "T cell" AND "peptide"',
    '"survivin" AND "CTL" AND "peptide"',
    '"PSMA" AND "T cell" AND "peptide"',
    '"claudin" AND "T cell" AND "peptide" AND "HLA"',
    '"GPC3" AND "T cell" AND "peptide" AND "HLA"',
    '"MUC1" AND "CTL" AND "epitope" AND "HLA"',
    '"HPV E6 E7" AND "CTL" AND "peptide" AND "HLA"',
    '"EBV LMP" AND "T cell" AND "peptide" AND "HLA"',
    '"CMV pp65" AND "T cell" AND "peptide" AND "HLA"',
    '"HTLV Tax" AND "CTL" AND "peptide"',
    '"HIV gag" AND "CTL" AND "peptide" AND "HLA"',
    '"influenza" AND "CTL" AND "peptide" AND "HLA" AND "epitope"',
    '"hepatitis B" AND "CTL" AND "peptide" AND "HLA"',
    '"SARS-CoV-2" AND "T cell" AND "epitope" AND "HLA"',
    # ── BLOCK 9: Eluted peptides / immunopeptidome (HIGH YIELD — peptide sequences) ──
    '"immunopeptidome" AND "HLA" AND "ligand"',
    '"HLA ligandome" AND "mass spectrometry"',
    '"eluted peptide" AND "HLA" AND "tumor"',
    '"MHC peptidome" AND "healthy tissue"',
    '"immunopeptidome" AND "melanoma" AND "HLA"',
    '"immunopeptidome" AND "lung cancer" AND "HLA"',
    '"immunopeptidome" AND "breast cancer" AND "HLA"',
    '"immunopeptidome" AND "leukemia" AND "HLA"',
    '"immunopeptidome" AND "glioblastoma" AND "HLA"',
    '"immunopeptidome" AND "ovarian cancer" AND "HLA"',
    '"immunopeptidome" AND "colon cancer" AND "HLA"',
    '"immunopeptidome" AND "renal cell" AND "HLA"',
    '"immunopeptidome" AND "pancreatic" AND "HLA"',
    '"HLA ligand atlas" AND "tissue"',
    '"MHC ligand" AND "mass spectrometry" AND "normal tissue"',
    '"MHC I" AND "eluted" AND "peptidome"',
    '"antigen presentation" AND "mass spectrometry" AND "HLA" AND "peptide"',
    '"tumor-specific antigen" AND "mass spectrometry" AND "HLA"',
    '"surface HLA" AND "peptide" AND "proteomics"',
    # ── BLOCK 10: Epitope databases & systematic epitope mapping ──
    '"T cell epitope" AND "database" AND "HLA"',
    '"IEDB" AND "T cell" AND "epitope" AND "HLA"',
    '"systematic epitope mapping" AND "T cell"',
    '"epitope mapping" AND "overlapping peptide" AND "T cell"',
    '"tetramer staining" AND "peptide" AND "HLA" AND "specific"',
    '"ELISPOT" AND "peptide" AND "HLA" AND "T cell"',
    '"intracellular cytokine staining" AND "peptide" AND "HLA"',
    '"IFN-gamma" AND "peptide" AND "T cell" AND "HLA"',
    '"cytotoxicity assay" AND "peptide" AND "HLA" AND "CTL"',
    # ── BLOCK 11: Structural TCR-pMHC studies ──
    '"TCR-pMHC" AND "crystal structure" AND "peptide"',
    '"T cell receptor" AND "structure" AND "cross-reactivity" AND "peptide"',
    '"TCR" AND "structural basis" AND "peptide recognition"',
    '"MHC class I" AND "peptide" AND "crystal structure" AND "TCR"',
    '"peptide conformation" AND "TCR recognition" AND "HLA"',
    '"altered peptide ligand" AND "TCR" AND "cross-reactivity"',
    '"hotspot" AND "TCR" AND "peptide" AND "binding"',
    # ── BLOCK 12: Neoantigen / personalized immunotherapy (peptide-rich) ──
    '"neoantigen" AND "peptide" AND "HLA" AND "patient"',
    '"personalized vaccine" AND "neoantigen" AND "peptide" AND "HLA"',
    '"neoantigen" AND "mass spectrometry" AND "validated"',
    '"tumor mutation burden" AND "neoantigen" AND "peptide"',
    '"neoantigen" AND "immunopeptidome" AND "validation"',
    '"mutant peptide" AND "HLA" AND "T cell" AND "recognition"',
    '"neoepitope" AND "HLA" AND "binding" AND "T cell"',
    '"somatic mutation" AND "peptide" AND "HLA" AND "CTL"',
    # ── BLOCK 13: Peptide vaccine studies (peptide sequences always listed) ──
    '"peptide vaccine" AND "cancer" AND "HLA" AND "clinical trial"',
    '"peptide vaccine" AND "melanoma" AND "CTL"',
    '"peptide vaccine" AND "HLA-A*02:01" AND "clinical"',
    '"multi-peptide vaccine" AND "cancer" AND "T cell"',
    '"long peptide vaccine" AND "T cell" AND "HLA"',
    '"synthetic long peptide" AND "immunotherapy" AND "T cell"',
    '"minimal epitope" AND "CTL" AND "HLA" AND "tumor"',
    # ── BLOCK 14: Minor histocompatibility antigens (transplant cross-reactivity) ──
    '"minor histocompatibility antigen" AND "peptide" AND "HLA"',
    '"minor H antigen" AND "T cell" AND "graft versus host"',
    '"HA-1" AND "HLA" AND "peptide" AND "T cell"',
    '"mismatch" AND "HLA" AND "peptide" AND "alloimmune"',
    '"alloreactive T cell" AND "peptide" AND "HLA"',
]

# Strategy B: Gene-focused queries — one per gene family
STRATEGY_B_GENES = [
    # Classical decoy / cross-reactivity genes
    "TTN", "MLANA", "PMEL", "TYR", "MAGEA3", "MAGEA4", "MAGEA1",
    # Cancer-testis antigens
    "CTAG1B", "CTAG2", "SSX2", "SSX1", "SSX4", "LAGE1", "BAGE", "BAGE2",
    "GAGE", "SAGE1", "XAGE1", "PAGE4", "SPANX", "SPANXB1",
    "MAGEA10", "MAGEA12", "MAGEB2", "MAGEC1", "MAGEC2",
    "CT45A1", "CT47A1", "CT83", "ADAM2", "Boris", "TDRD1",
    # Tumor-associated antigens
    "WT1", "PRAME", "AFP", "ERBB2", "MUC1", "TERT", "MSLN",
    "CEACAM5", "BIRC5", "ACPP", "FOLH1", "GPC3",
    # Melanoma differentiation antigens
    "SLC45A2", "TYRP1", "DCT", "MC1R", "MITF", "gp100",
    # Oncogenes / neoantigen-relevant
    "KRAS", "TP53", "BRAF", "CDK4", "EGFR", "PIK3CA", "NRAS",
    "HRAS", "IDH1", "IDH2", "ALK", "ROS1", "RET", "MET",
    "CTNNB1", "FBXW7", "SMAD4", "APC", "VHL", "PTEN",
    "NPM1", "FLT3", "JAK2", "CALR", "BCR-ABL",
    # Immuno-oncology surface targets
    "CLEC12A", "CD33", "FOLR1", "NECTIN4", "CLDN18", "CLDN6",
    "DLL3", "CD276", "ROR1", "CD19", "CD20", "CD22", "BCMA",
    "GPNMB", "TROP2", "HER3", "FGFR2", "FGFR3",
    # Viral antigens (rich in known epitopes)
    "HBV", "HPV16", "HPV18", "EBV", "CMV", "HTLV",
    "HIV", "HCV", "influenza", "SARS-CoV-2",
    "EBNA1", "EBNA3", "LMP1", "LMP2",
    # Minor histocompatibility antigens
    "HA-1", "HA-2", "HB-1", "UGT2B17", "HMHA1",
    # Autoimmune-relevant self-antigens
    "MBP", "MOG", "PLP1", "GAD65", "IAPP", "insulin",
    "thyroglobulin", "CYP21A2", "TSHR",
    "collagen II", "aggrecan", "fibrillin",
    "desmoglein", "aquaporin-4",
]

# Strategy C: MeSH-term based systematic expansion
STRATEGY_C_MESH = [
    '"Receptors, Antigen, T-Cell"[MeSH] AND "Cross Reactions"[MeSH] AND "Peptides"[MeSH]',
    '"Immunotherapy, Adoptive"[MeSH] AND "Drug-Related Side Effects and Adverse Reactions"[MeSH]',
    '"T-Lymphocytes, Cytotoxic"[MeSH] AND "Molecular Mimicry"[MeSH]',
    '"HLA-A Antigens"[MeSH] AND "Epitopes, T-Lymphocyte"[MeSH] AND "Cross Reactions"[MeSH]',
    '"Neoplasm Proteins"[MeSH] AND "T-Lymphocytes, Cytotoxic"[MeSH] AND "HLA Antigens"[MeSH]',
    '"Cancer Vaccines"[MeSH] AND "Autoimmunity"[MeSH] AND "Peptides"[MeSH]',
    '"Receptors, Chimeric Antigen"[MeSH] AND "Off-Target Effects"[MeSH]',
    '"Melanoma"[MeSH] AND "T-Lymphocytes, Cytotoxic"[MeSH] AND "Antigens, Neoplasm"[MeSH]',
    '"Leukemia"[MeSH] AND "Immunotherapy, Adoptive"[MeSH] AND "Antigens, Differentiation"[MeSH]',
    # Additional MeSH — immunopeptidome, epitope discovery, vaccines
    '"Proteomics"[MeSH] AND "Histocompatibility Antigens Class I"[MeSH] AND "Peptides"[MeSH]',
    '"Mass Spectrometry"[MeSH] AND "HLA Antigens"[MeSH] AND "Peptides"[MeSH]',
    '"Epitopes, T-Lymphocyte"[MeSH] AND "Peptide Mapping"[MeSH]',
    '"Cancer Vaccines"[MeSH] AND "Peptides"[MeSH] AND "Clinical Trials as Topic"[MeSH]',
    '"Antigens, Neoplasm"[MeSH] AND "Peptides"[MeSH] AND "T-Lymphocytes, Cytotoxic"[MeSH]',
    '"Minor Histocompatibility Antigens"[MeSH] AND "Graft vs Host Disease"[MeSH]',
    '"Autoimmune Diseases"[MeSH] AND "Molecular Mimicry"[MeSH] AND "Peptides"[MeSH]',
    '"Neoantigen"[MeSH] AND "Immunogenicity, Vaccine"[MeSH]',
    '"Antigen Presentation"[MeSH] AND "Peptide Fragments"[MeSH] AND "HLA Antigens"[MeSH]',
    '"T-Cell Antigen Receptor Specificity"[MeSH] AND "Peptides"[MeSH]',
    '"Viral Vaccines"[MeSH] AND "Epitopes, T-Lymphocyte"[MeSH] AND "HLA Antigens"[MeSH]',
]

# Strategy D: Year-range expansion — systematically sweep publication years
def generate_year_queries(base_queries: List[str], year_start: int = 1990, year_end: int = 2026) -> List[str]:
    """Add year filters to base queries to pull different result sets."""
    expanded = []
    for q in base_queries[:8]:  # use a subset of base queries
        for y_start in range(year_start, year_end, 5):
            y_end = min(y_start + 4, year_end)
            expanded.append(f'({q}) AND ("{y_start}"[PDAT] : "{y_end}"[PDAT])')
    return expanded


# Strategy E: LLM-guided query generation
def generate_llm_queries(existing_entries: List[DecoyEntry], n_queries: int = 20) -> List[str]:
    """
    Use LLM to propose novel PubMed queries based on patterns in existing data.
    Falls back to empty list if no LLM available.
    """
    from decoy_library.config import OPENAI_API_KEY, OPENAI_BASE_URL, OPENAI_MODEL
    if not OPENAI_API_KEY:
        log.warning("No LLM API key — skipping LLM-guided query generation")
        return []

    try:
        from openai import OpenAI
    except ImportError:
        return []

    # Summarise what we already have
    genes = sorted(set(e.peptide_info.gene_symbol for e in existing_entries))
    hla_alleles = sorted(set(e.peptide_info.hla_allele for e in existing_entries))
    proteins = sorted(set(e.peptide_info.source_protein for e in existing_entries))[:30]

    prompt = (
        "You are a biomedical literature search expert.\n"
        "Our Decoy Peptide Library already covers these genes: " + ", ".join(genes[:40]) + "\n"
        "HLA alleles: " + ", ".join(hla_alleles) + "\n"
        "Proteins: " + ", ".join(proteins) + "\n\n"
        f"Generate {n_queries} diverse PubMed search queries to find NEW papers about:\n"
        "- TCR cross-reactivity with off-target peptides\n"
        "- Immunotherapy toxicity from peptide mimicry\n"
        "- HLA-restricted peptides that could be decoys\n"
        "- Immunopeptidome / HLA-ligandome studies\n"
        "- Neoantigen cross-reactivity\n"
        "- Peptide-MHC binding with unexpected TCR recognition\n\n"
        "Focus on genes/HLA alleles NOT yet in our library.\n"
        "Output ONLY a JSON array of query strings. No markdown.\n"
    )

    kwargs = {"api_key": OPENAI_API_KEY}
    if OPENAI_BASE_URL:
        kwargs["base_url"] = OPENAI_BASE_URL
    client = OpenAI(**kwargs)

    try:
        for attempt in range(3):  # Retry up to 3 times for malformed JSON
            resp = client.chat.completions.create(
                model=OPENAI_MODEL,
                messages=[
                    {"role": "system", "content": "You are a PubMed search query generator. Output only a valid JSON array of strings, nothing else."},
                    {"role": "user", "content": prompt},
                ],
                temperature=0.8,
                max_tokens=2048,
            )
            raw = resp.choices[0].message.content or "[]"
            raw = raw.strip()
            if raw.startswith("```"):
                import re
                raw = re.sub(r"^```(?:json)?\s*\n?", "", raw)
                raw = re.sub(r"\n?```\s*$", "", raw)
            try:
                queries = json.loads(raw)
                if isinstance(queries, list) and len(queries) > 0:
                    log.info("LLM generated %d new search queries (attempt %d)", len(queries), attempt + 1)
                    return [str(q) for q in queries][:n_queries]
            except json.JSONDecodeError as je:
                log.warning("LLM JSON parse attempt %d failed: %s", attempt + 1, je)
                continue
        log.warning("LLM query generation: all 3 attempts returned invalid JSON")
    except Exception as exc:
        log.warning("LLM query generation failed: %s", exc)
    return []


# Strategy F: Multi-source queries (arXiv, bioRxiv, ClinicalTrials, etc.)
MULTI_SOURCE_QUERIES = {
    "arxiv": [
        "TCR T cell receptor cross-reactivity",
        "immunotherapy off-target toxicity peptide",
        "peptide MHC binding prediction safety",
        "T cell epitope cross-recognition",
        "CAR-T cell therapy adverse events",
        "neoantigen TCR recognition",
        "HLA peptide binding affinity",
        "molecular mimicry autoimmunity TCR",
        "immunopeptidome mass spectrometry HLA",
        "deep learning peptide MHC binding",
        "TCR repertoire specificity prediction",
        "pMHC tetramer epitope discovery",
    ],
    "biorxiv": [
        "TCR cross-reactivity off-target",
        "T cell receptor peptide specificity",
        "immunotherapy cardiac toxicity",
        "HLA ligandome tumor",
        "adoptive cell therapy safety",
        "peptide vaccine cross-reactive",
        "immunopeptidome cancer",
        "neoantigen peptide HLA",
        "MHC peptide elution mass spectrometry",
        "T cell epitope mapping",
        "TCR-pMHC structure",
        "peptide binding affinity HLA",
    ],
    "medrxiv": [
        "TCR therapy clinical trial",
        "CAR-T adverse events",
        "immunotherapy toxicity report",
        "T cell therapy safety profile",
        "peptide vaccine clinical",
        "neoantigen vaccine cancer",
        "adoptive cell therapy outcome",
        "immune checkpoint inhibitor autoimmune",
    ],
    "clinicaltrials": [
        "TCR T cell receptor therapy cancer",
        "MAGE-A3 TCR adoptive",
        "NY-ESO-1 TCR therapy",
        "adoptive cell transfer solid tumor",
        "engineered T cell immunotherapy",
        "affinity enhanced TCR",
        "peptide vaccine melanoma",
        "HLA-A2 restricted TCR",
        "neoantigen vaccine personalized",
        "peptide vaccine cancer phase",
        "WT1 peptide vaccine",
        "gp100 peptide melanoma",
        "survivin peptide cancer",
        "multi-epitope vaccine cancer",
        "PRAME TCR therapy",
    ],
    "semantic_scholar": [
        "TCR cross-reactivity off-target peptide toxicity",
        "T cell receptor cardiac toxicity MAGE",
        "HLA-restricted peptide decoy mimicry",
        "immunotherapy autoimmune cross-reactive",
        "adoptive T cell therapy safety failure",
        "immunopeptidome HLA ligand mass spectrometry",
        "neoantigen peptide HLA binding validation",
        "cancer testis antigen T cell epitope",
        "peptide vaccine cancer clinical trial HLA",
        "minor histocompatibility antigen peptide",
        "molecular mimicry autoimmune disease peptide",
        "TCR-pMHC structural cross-reactivity",
        "HLA ligandome tissue-specific peptide",
        "viral epitope CTL HLA peptide",
        "tumor antigen T cell recognition peptide",
    ],
    "europepmc": [
        "TCR cross-reactivity peptide",
        "T cell therapy toxicity",
        "HLA peptide immunotherapy",
        "off-target T cell recognition",
        "adoptive cell therapy adverse",
        "immunopeptidome tumor HLA ligand",
        "neoantigen mass spectrometry validation",
        "peptide vaccine HLA-A*02:01",
        "cancer testis antigen CTL epitope",
        "minor histocompatibility antigen graft",
        "molecular mimicry T cell autoimmune",
        "viral peptide HLA T cell epitope",
        "eluted peptide HLA mass spectrometry",
        "TCR-pMHC crystal structure peptide",
        "tumor-associated antigen peptide HLA",
    ],
}


def search_multi_sources_for_papers(
    lib: DecoyLibrary,
    max_per_source: int = 30,
    processed_ids: Set[str] = None
) -> List:
    """
    Search multiple non-PubMed sources for relevant papers.
    
    Uses fail-fast logic: after 2 consecutive empty results for a source,
    skip the remaining queries for that source (likely unreachable).
    
    Returns list of (source, query, papers) tuples.
    """
    if not MULTI_SOURCE_AVAILABLE:
        log.warning("Multi-source fetcher not available")
        return []
    
    if processed_ids is None:
        processed_ids = set()
    
    results = []
    source_empty_streak: Dict[str, int] = {}  # consecutive-empty count
    SKIP_THRESHOLD = 2                         # skip source after N empties
    
    for source, queries in MULTI_SOURCE_QUERIES.items():
        for query in queries:
            query_key = f"{source}:{query}"
            if query_key in processed_ids:
                continue
            
            # Skip source that appears unreachable / empty
            if source_empty_streak.get(source, 0) >= SKIP_THRESHOLD:
                log.info("  ⏭️  Skipping remaining '%s' queries (likely unreachable)", source)
                break
            
            try:
                if source == "arxiv":
                    papers = search_arxiv(query, max_per_source)
                elif source == "biorxiv":
                    papers = search_biorxiv(query, max_per_source, "biorxiv")
                elif source == "medrxiv":
                    papers = search_biorxiv(query, max_per_source, "medrxiv")
                elif source == "clinicaltrials":
                    papers = search_clinical_trials(query, max_per_source)
                elif source == "semantic_scholar":
                    papers = search_semantic_scholar(query, max_per_source)
                elif source == "europepmc":
                    papers = search_europe_pmc(query, max_per_source)
                else:
                    continue
                
                if papers:
                    results.append((source, query, papers))
                    source_empty_streak[source] = 0          # reset on success
                else:
                    source_empty_streak[source] = source_empty_streak.get(source, 0) + 1
                    
            except Exception as exc:
                log.warning("Error searching %s with '%s': %s", source, query[:30], exc)
                source_empty_streak[source] = source_empty_streak.get(source, 0) + 1
    
    return results


# ─────────────────────────────────────────────────────────────────────────
# 2.  Scale-Up Orchestrator
# ─────────────────────────────────────────────────────────────────────────

CHECKPOINT_FILE = DATA_DIR / "scale_up_checkpoint.json"


def _load_checkpoint() -> Dict:
    """Load progress checkpoint."""
    if CHECKPOINT_FILE.exists():
        return json.loads(CHECKPOINT_FILE.read_text())
    return {
        "processed_pmids": [],
        "processed_queries": [],
        "processed_multi_source": [],  # Track multi-source queries
        "started_at": datetime.now().isoformat(),
        "rounds": 0,
    }


def _save_checkpoint(ckpt: Dict):
    """Persist checkpoint to disk."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINT_FILE.write_text(json.dumps(ckpt, indent=2))


def _unique_sequences(lib: DecoyLibrary) -> Set[str]:
    return {e.peptide_info.decoy_sequence for e in lib.entries}


def _quality_filter(lib: DecoyLibrary) -> DecoyLibrary:
    """
    Remove low-quality entries that shouldn't count toward target N.
    Keeps entries that have at least ONE of:
        - non-UNKNOWN gene symbol
        - valid HLA allele (not HLA-A*00:xx)
        - evidence summary length > 20
    """
    good = []
    dropped = 0
    for e in lib.entries:
        gene_ok = e.peptide_info.gene_symbol not in ("UNKNOWN", "", None)
        hla_ok = "00:xx" not in e.peptide_info.hla_allele
        summary_ok = len(e.provenance.evidence_summary) > 20
        if gene_ok or hla_ok or summary_ok:
            good.append(e)
        else:
            dropped += 1
    if dropped:
        log.info("Quality filter dropped %d low-quality entries", dropped)
    lib.entries = good
    return lib


def scale_up(
    required_num: int,
    batch_size: int = 30,
    max_rounds: int = 500,
    validate: bool = True,
    resume: bool = True,
    strategy: str = "all",
    retmax_per_query: int = 40,
    save_every: int = 5,
    quality_filter: bool = True,
    multi_source: bool = False,
    skip_pubmed: bool = False,
    infinite: bool = False,
) -> DecoyLibrary:
    """
    Main scale-up loop.

    Parameters
    ----------
    required_num : int
        Target number of unique decoy entries.
    batch_size : int
        Number of PMIDs to process per round.
    max_rounds : int
        Safety limit on total rounds.
    validate : bool
        Whether to validate against UniProt/IEDB.
    resume : bool
        Whether to resume from checkpoint.
    strategy : str
        Query strategy: 'keyword', 'gene', 'mesh', 'year', 'llm', 'all'
    retmax_per_query : int
        Max results per PubMed query.
    save_every : int
        Save checkpoint every N rounds.
    quality_filter : bool
        Whether to apply quality filtering.
    multi_source : bool
        Whether to search additional sources (arXiv, bioRxiv, etc.)
    skip_pubmed : bool
        Whether to skip PubMed queries.
    infinite : bool
        Whether to loop indefinitely, generating new queries when exhausted.

    Returns
    -------
    DecoyLibrary
    """
    # ── Step 0: Load or initialise ──────────────────────────────────
    lib = load_library()
    if not lib.entries:
        log.info("📦 Library empty — running cold start first...")
        lib = run_cold_start(validate=validate)

    # Resume checkpoint
    ckpt = _load_checkpoint() if resume else {
        "processed_pmids": [],
        "processed_queries": [],
        "processed_multi_source": [],
        "started_at": datetime.now().isoformat(),
        "rounds": 0,
    }
    processed_pmids: Set[str] = set(ckpt.get("processed_pmids", []))
    processed_queries: Set[str] = set(ckpt.get("processed_queries", []))
    processed_multi_source: Set[str] = set(ckpt.get("processed_multi_source", []))

    current_count = len(lib.entries)
    log.info("=" * 60)
    log.info("🚀 SCALE-UP START")
    log.info("   Target:  %d entries", required_num)
    log.info("   Current: %d entries (%d unique sequences)", current_count, len(_unique_sequences(lib)))
    log.info("   Gap:     %d", max(0, required_num - current_count))
    log.info("   Strategy: %s", strategy)
    log.info("   Multi-source: %s", "enabled" if multi_source else "disabled")
    log.info("   Skip PubMed: %s", "yes" if skip_pubmed else "no")
    log.info("=" * 60)

    if current_count >= required_num:
        log.info("✅ Already have %d ≥ %d — nothing to do!", current_count, required_num)
        return lib

    # ── Step 1: Build query pool ────────────────────────────────────
    query_pool: List[str] = []

    if strategy in ("keyword", "all"):
        query_pool.extend(STRATEGY_A_QUERIES)
        log.info("  📋 Strategy A (keyword): %d queries", len(STRATEGY_A_QUERIES))

    if strategy in ("gene", "all"):
        gene_queries = [
            f'("{g}" AND "TCR" AND "peptide" AND "HLA")' for g in STRATEGY_B_GENES
        ]
        query_pool.extend(gene_queries)
        log.info("  🧬 Strategy B (gene): %d queries", len(gene_queries))

    if strategy in ("mesh", "all"):
        query_pool.extend(STRATEGY_C_MESH)
        log.info("  🏷️  Strategy C (MeSH): %d queries", len(STRATEGY_C_MESH))

    if strategy in ("year", "all"):
        year_queries = generate_year_queries(STRATEGY_A_QUERIES[:6])
        query_pool.extend(year_queries)
        log.info("  📅 Strategy D (year): %d queries", len(year_queries))

    if strategy in ("llm", "all"):
        llm_queries = generate_llm_queries(lib.entries)
        query_pool.extend(llm_queries)
        log.info("  🤖 Strategy E (LLM): %d queries", len(llm_queries))

    # Shuffle to get diversity early
    random.shuffle(query_pool)

    # Filter out already-processed queries
    query_pool = [q for q in query_pool if q not in processed_queries]
    log.info("  📊 Total query pool: %d (after removing %d already-processed)",
             len(query_pool), len(processed_queries))

    # ── Step 2: Main PubMed loop ──────────────────────────────────────
    round_num = ckpt.get("rounds", 0)
    total_new = 0
    failed_queries = 0
    empty_queries = 0
    multi_source_papers_processed = 0

    # Skip PubMed queries if --no-pubmed is set
    if skip_pubmed:
        log.info("⏭️  Skipping PubMed queries (--no-pubmed)")
        query_pool = []

    for qi, query in enumerate(query_pool):
        if len(lib.entries) >= required_num:
            break
        if round_num >= max_rounds:
            log.warning("⚠️  Hit max_rounds=%d limit", max_rounds)
            break

        round_num += 1
        gap = required_num - len(lib.entries)

        log.info("")
        log.info("━" * 60)
        log.info("🔄 Round %d | Query %d/%d | Have: %d | Need: %d more",
                 round_num, qi + 1, len(query_pool), len(lib.entries), gap)
        log.info("   Query: %s", query[:80])
        log.info("━" * 60)

        # Search PubMed
        try:
            pmids = search_pubmed(query, retmax=retmax_per_query)
        except Exception as exc:
            log.error("❌ PubMed search failed: %s", exc)
            failed_queries += 1
            processed_queries.add(query)
            continue

        # Filter out already-processed PMIDs
        new_pmids = [p for p in pmids if p not in processed_pmids]
        if not new_pmids:
            log.info("   ⏭️  No new PMIDs (all %d already processed)", len(pmids))
            empty_queries += 1
            processed_queries.add(query)
            continue

        # Take a batch
        batch = new_pmids[:batch_size]
        log.info("   📄 Found %d PMIDs, %d new, processing %d",
                 len(pmids), len(new_pmids), len(batch))

        # Fetch papers
        try:
            papers = fetch_multiple(batch)
        except Exception as exc:
            log.error("❌ Fetch failed: %s", exc)
            processed_pmids.update(batch)
            continue

        # Process: extract + validate + deduplicate
        before = len(lib.entries)
        try:
            new_entries = process_papers(papers, lib, validate=validate)
        except Exception as exc:
            log.error("❌ Processing failed: %s", exc)
            new_entries = []

        after = len(lib.entries)
        added = after - before
        total_new += added

        processed_pmids.update(batch)
        processed_queries.add(query)

        log.info("   ✅ Added %d new entries (extracted %d, deduped %d)",
                 added, len(new_entries), len(new_entries) - added)
        log.info("   📊 Library: %d / %d (%.1f%%)",
                 len(lib.entries), required_num,
                 100.0 * len(lib.entries) / required_num)

        # Checkpoint
        if round_num % save_every == 0:
            save_library(lib)
            ckpt.update({
                "processed_pmids": list(processed_pmids),
                "processed_queries": list(processed_queries),
                "rounds": round_num,
                "library_size": len(lib.entries),
                "last_saved": datetime.now().isoformat(),
            })
            _save_checkpoint(ckpt)
            log.info("   💾 Checkpoint saved (round %d)", round_num)

    # ── Step 2.5: Multi-source search ───────────────────────────────
    if multi_source and len(lib.entries) < required_num and MULTI_SOURCE_AVAILABLE:
        log.info("")
        log.info("=" * 60)
        log.info("🌐 MULTI-SOURCE SEARCH")
        log.info("=" * 60)
        
        ms_results = search_multi_sources_for_papers(
            lib, max_per_source=30, processed_ids=processed_multi_source
        )
        
        log.info("  Found %d source/query combinations with papers", len(ms_results))
        
        for source, query, papers in ms_results:
            if len(lib.entries) >= required_num:
                break
            
            query_key = f"{source}:{query}"
            if query_key in processed_multi_source:
                continue
            
            round_num += 1
            gap = required_num - len(lib.entries)
            
            log.info("")
            log.info("━" * 60)
            log.info("🔄 Round %d | Source: %s | Have: %d | Need: %d more",
                     round_num, source.upper(), len(lib.entries), gap)
            log.info("   Query: %s", query[:60])
            log.info("   Papers: %d", len(papers))
            log.info("━" * 60)
            
            # Convert multi-source PaperRecords → fetcher.PaperRecord
            # so the existing process_papers pipeline works unchanged.
            converted: List[PaperRecord] = []
            for p in papers[:batch_size]:
                pr = PaperRecord(
                    pmid=p.source_id,           # use source_id as pseudo-PMID
                    title=p.title or "",
                    abstract=p.abstract or "",
                    authors=p.authors if p.authors else [],
                    journal=f"[{source}] {p.doi}" if p.doi else f"[{source}]",
                    year=str(p.year) if p.year else "",
                    doi=p.doi or "",
                    full_text=p.full_text or "",
                )
                converted.append(pr)
            
            if not converted:
                processed_multi_source.add(query_key)
                continue
            
            # Feed through the standard pipeline
            before = len(lib.entries)
            try:
                new_entries = process_papers(converted, lib, validate=validate)
            except Exception as exc:
                log.error("❌ Multi-source processing failed for %s: %s", source, exc)
                new_entries = []
            
            after = len(lib.entries)
            added = after - before
            total_new += added
            multi_source_papers_processed += len(converted)
            processed_multi_source.add(query_key)
            
            log.info("   ✅ Added %d new entries from %s (processed %d papers)",
                     added, source, len(converted))
            log.info("   📊 Library: %d / %d (%.1f%%)",
                     len(lib.entries), required_num,
                     100.0 * len(lib.entries) / required_num)
            
            # Checkpoint
            if round_num % save_every == 0:
                save_library(lib)
                ckpt.update({
                    "processed_pmids": list(processed_pmids),
                    "processed_queries": list(processed_queries),
                    "processed_multi_source": list(processed_multi_source),
                    "rounds": round_num,
                    "library_size": len(lib.entries),
                    "last_saved": datetime.now().isoformat(),
                })
                _save_checkpoint(ckpt)
                log.info("   💾 Checkpoint saved (round %d)", round_num)

    # ── Step 3: Quality filter ──────────────────────────────────────
    if quality_filter:
        before_qf = len(lib.entries)
        lib = _quality_filter(lib)
        log.info("Quality filter: %d → %d entries", before_qf, len(lib.entries))

    # ── Step 4: Final save ──────────────────────────────────────────
    save_library(lib)
    ckpt.update({
        "processed_pmids": list(processed_pmids),
        "processed_queries": list(processed_queries),
        "processed_multi_source": list(processed_multi_source),
        "rounds": round_num,
        "library_size": len(lib.entries),
        "completed_at": datetime.now().isoformat(),
    })
    _save_checkpoint(ckpt)

    # ── Step 5: Summary ─────────────────────────────────────────────
    unique_seqs = _unique_sequences(lib)
    log.info("")
    log.info("=" * 60)
    log.info("🏁 SCALE-UP COMPLETE")
    log.info("=" * 60)
    log.info("   Total entries:     %d", len(lib.entries))
    log.info("   Unique sequences:  %d", len(unique_seqs))
    log.info("   Target:            %d", required_num)
    log.info("   Achieved:          %.1f%%", 100.0 * len(lib.entries) / required_num)
    log.info("   Rounds executed:   %d", round_num)
    log.info("   New entries added: %d", total_new)
    log.info("   Failed queries:    %d", failed_queries)
    log.info("   Empty queries:     %d", empty_queries)
    log.info("   PMIDs processed:   %d", len(processed_pmids))
    if multi_source:
        log.info("   Multi-source papers: %d", multi_source_papers_processed)
        log.info("   Multi-source queries: %d", len(processed_multi_source))

    if len(lib.entries) < required_num:
        shortfall = required_num - len(lib.entries)
        log.warning("")
        log.warning("⚠️  Fell short by %d entries.", shortfall)
        log.warning("   Suggestions:")
        log.warning("   1. Run again (--resume) — new LLM queries may find more")
        log.warning("   2. Add custom queries: python -m decoy_library search --query '...'")
        log.warning("   3. Lower quality bar: --no-quality-filter")
        log.warning("   4. Increase retmax: --retmax 100")

    # Evidence level distribution
    levels: Dict[str, int] = {}
    for e in lib.entries:
        lv = e.risk_profile.evidence_level.value
        levels[lv] = levels.get(lv, 0) + 1
    log.info("")
    log.info("Evidence level distribution:")
    for lv, cnt in sorted(levels.items()):
        log.info("   %s: %d", lv, cnt)

    return lib


# ─────────────────────────────────────────────────────────────────────────
# 3.  CLI Entry Point
# ─────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        prog="scale_up",
        description="Scale up the Decoy C Library to N unique high-quality entries",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scale_up.py --required_num 500
  python scale_up.py --required_num 10000 --strategy all --batch-size 50
  python scale_up.py --required_num 10000 --resume
  python scale_up.py --required_num 200  --no-validate --strategy keyword
        """,
    )

    parser.add_argument(
        "--required_num", "-n",
        type=int, required=True,
        help="Target number of unique decoy entries to generate",
    )
    parser.add_argument(
        "--batch-size", "-b",
        type=int, default=30,
        help="Number of PMIDs to process per round (default: 30)",
    )
    parser.add_argument(
        "--max-rounds",
        type=int, default=500,
        help="Maximum number of query rounds (default: 500)",
    )
    parser.add_argument(
        "--strategy", "-s",
        choices=["keyword", "gene", "mesh", "year", "llm", "all"],
        default="all",
        help="Query expansion strategy (default: all)",
    )
    parser.add_argument(
        "--retmax",
        type=int, default=40,
        help="Max results per PubMed query (default: 40)",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip UniProt/IEDB validation (faster but less quality)",
    )
    parser.add_argument(
        "--no-quality-filter",
        action="store_true",
        help="Disable quality filtering on final output",
    )
    parser.add_argument(
        "--resume",
        action="store_true", default=True,
        help="Resume from last checkpoint (default: True)",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Ignore checkpoint and start fresh",
    )
    parser.add_argument(
        "--save-every",
        type=int, default=5,
        help="Save checkpoint every N rounds (default: 5)",
    )
    parser.add_argument(
        "--multi-source",
        action="store_true",
        help="Enable multi-source search (arXiv, bioRxiv, ClinicalTrials, etc.)",
    )
    parser.add_argument(
        "--no-pubmed",
        action="store_true",
        help="Skip PubMed queries (only use multi-source if --multi-source is set)",
    )
    parser.add_argument(
        "--infinite",
        action="store_true",
        help="Run indefinitely, generating new queries when exhausted",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging",
    )

    args = parser.parse_args()

    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    print(f"\n{'='*60}")
    print(f"  Decoy C Library — Scale-Up Generator")
    print(f"  Target: {args.required_num} unique entries")
    print(f"  Strategy: {args.strategy}")
    print(f"  Validate: {not args.no_validate}")
    print(f"  Multi-source: {args.multi_source}")
    print(f"{'='*60}\n")

    cycle = 0
    while True:
        cycle += 1
        if cycle > 1:
            print(f"\n{'🔄'*20}")
            print(f"  INFINITE MODE — Cycle {cycle}")
            print(f"  Clearing query checkpoint for fresh LLM queries...")
            print(f"{'🔄'*20}\n")
            # Clear query-related checkpoint data but keep processed PMIDs
            # so we don't re-process the same papers.
            ckpt = _load_checkpoint()
            ckpt["processed_queries"] = []
            ckpt["processed_multi_source"] = []
            _save_checkpoint(ckpt)
            time.sleep(10)  # Brief pause between cycles

        lib = scale_up(
            required_num=args.required_num,
            batch_size=args.batch_size,
            max_rounds=args.max_rounds,
            validate=not args.no_validate,
            resume=True if cycle > 1 else not args.fresh,  # First cycle respects --fresh, later always resume
            strategy=args.strategy,
            retmax_per_query=args.retmax,
            save_every=args.save_every,
            quality_filter=not args.no_quality_filter,
            multi_source=args.multi_source,
            skip_pubmed=args.no_pubmed,
            infinite=args.infinite,
        )

        print(f"\n{'='*60}")
        print(f"  Cycle {cycle} complete: {len(lib.entries)} entries")
        print(f"  📁 Saved to: {LIBRARY_JSON}")
        print(f"{'='*60}\n")

        # Check if we've reached the target
        if len(lib.entries) >= args.required_num:
            print(f"  ✅ TARGET REACHED! {len(lib.entries)} ≥ {args.required_num}")
            break

        # If not infinite mode, stop after one cycle
        if not args.infinite:
            print(f"  ✅ Final library: {len(lib.entries)} entries")
            break

        print(f"  📊 Progress: {len(lib.entries)} / {args.required_num}")
        print(f"  🔄 Infinite mode: starting cycle {cycle + 1}...")


if __name__ == "__main__":
    main()
