#!/usr/bin/env python3
"""
IEDB Bulk Miner
===============
Systematically mines the Immune Epitope Database (IEDB) as a **primary data
source** for the Decoy C Library.

IEDB contains ~1 million T-cell epitope records with experimentally validated
peptide-MHC-TCR interactions.  This module queries the IEDB T-cell assay and
epitope APIs to discover decoy-relevant entries that literature search alone
would miss.

Mining strategies
-----------------
1. **T-cell assay search** — query for assays involving cross-reactivity,
   autoimmune self-antigens, and off-target recognition.
2. **Source-organism filter** — mine all human-self-antigen epitopes that
   bind common HLA alleles (potential off-target pool).
3. **Protein-specific search** — for known dangerous proteins (TTN, MAGEA3,
   cardiac, neural) pull ALL known epitopes.
4. **MHC-elution / mass-spec confirmed** — highest-confidence naturally
   presented peptides.

Public API
----------
    mine_iedb(lib, strategy="all") -> list[DecoyEntry]
"""

from __future__ import annotations

import logging
import re
import time
from typing import Any, Dict, List, Optional, Set

import requests

from .config import IEDB_DELAY_SEC
from .models import (
    DecoyEntry,
    DecoyLibrary,
    DiscoveryContext,
    EvidenceLevel,
    ExperimentalEvidence,
    PeptideInfo,
    Provenance,
    RiskProfile,
)

log = logging.getLogger(__name__)

# ── IEDB API endpoints ────────────────────────────────────────────────
IEDB_TCELL = "https://query-api.iedb.org/tcell_search"
IEDB_EPITOPE = "https://query-api.iedb.org/epitope_search"

# ── Rate limiting ─────────────────────────────────────────────────────
_last_call = 0.0


def _rate_limit():
    global _last_call
    elapsed = time.time() - _last_call
    if elapsed < IEDB_DELAY_SEC:
        time.sleep(IEDB_DELAY_SEC - elapsed)
    _last_call = time.time()


# ── Helper: query IEDB T-cell assay API ───────────────────────────────

def _query_tcell(params: Dict[str, str], max_results: int = 500) -> List[Dict]:
    """
    Query the IEDB T-cell assay search API.

    Returns a list of result dicts (each representing one T-cell assay record).
    """
    _rate_limit()
    try:
        resp = requests.get(IEDB_TCELL, params=params, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, list):
                return data[:max_results]
        else:
            log.warning("IEDB T-cell API returned %d for params %s", resp.status_code, params)
    except Exception as exc:
        log.warning("IEDB T-cell query failed: %s", exc)
    return []


def _query_epitope(params: Dict[str, str], max_results: int = 500) -> List[Dict]:
    """Query the IEDB epitope search API."""
    _rate_limit()
    try:
        resp = requests.get(IEDB_EPITOPE, params=params, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            if isinstance(data, list):
                return data[:max_results]
        else:
            log.warning("IEDB epitope API returned %d for params %s", resp.status_code, params)
    except Exception as exc:
        log.warning("IEDB epitope query failed: %s", exc)
    return []


# ── Conversion: IEDB record → DecoyEntry ─────────────────────────────

_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")
_HLA_RE = re.compile(r"HLA-[A-Z]+\d?\*\d{2}:\d{2}")


def _normalise_hla(raw: str) -> Optional[str]:
    """Try to extract a standard HLA-X*NN:NN allele from IEDB allele strings."""
    if not raw:
        return None
    m = _HLA_RE.search(raw)
    if m:
        return m.group(0)
    # Try low-res: "HLA-A2" → "HLA-A*02:01"
    m2 = re.match(r"HLA-([A-Z])(\d+)", raw)
    if m2:
        return f"HLA-{m2.group(1)}*{int(m2.group(2)):02d}:01"
    return None


def _determine_evidence_level(record: Dict) -> EvidenceLevel:
    """
    Map IEDB assay metadata to our 4-tier evidence level.

    - Mass-spec (elution) confirmed → Level 2 (in vitro confirmed)
    - Cytotoxicity / killing assay → Level 2
    - Binding / ELISPOT / tetramer → Level 3
    - Default → Level 3 (high-throughput screened)
    """
    assay_type = str(record.get("assay_type", "")).lower()
    qualitative = str(record.get("qualitative_measure", "")).lower()

    if "cytotoxicity" in assay_type or "killing" in assay_type or "51cr" in assay_type:
        return EvidenceLevel.LEVEL_2
    if "elution" in assay_type or "mass spec" in assay_type:
        return EvidenceLevel.LEVEL_2
    if "ifn" in assay_type or "elispot" in assay_type or "tetramer" in assay_type:
        return EvidenceLevel.LEVEL_3
    return EvidenceLevel.LEVEL_3


_PROTEIN_TO_GENE: Dict[str, str] = {
    # Cardiac / muscle
    "titin": "TTN", "Titin": "TTN", "myosin-6": "MYH6", "myosin-7": "MYH7",
    "troponin T": "TNNT2", "troponin I": "TNNI3", "desmin": "DES",
    "ryanodine receptor 2": "RYR2", "dystrophin": "DMD",
    "Actin, alpha cardiac": "ACTC1", "Tropomyosin": "TPM1",
    # Melanocyte differentiation
    "Melanocyte protein PMEL": "PMEL", "gp100": "PMEL",
    "Melan-A": "MLANA", "MART-1": "MLANA",
    "Tyrosinase": "TYR", "L-dopachrome tautomerase": "DCT",
    "5,6-dihydroxyindole-2-carboxylic acid oxidase": "TYRP1",
    # Cancer-testis
    "Melanoma-associated antigen 3": "MAGEA3", "MAGE-A3": "MAGEA3",
    "Melanoma-associated antigen 4": "MAGEA4",
    "Melanoma-associated antigen 1": "MAGEA1",
    "Cancer/testis antigen 1": "CTAG1B", "NY-ESO-1": "CTAG1B",
    "Protein SSX2": "SSX2", "PRAME": "PRAME", "WT1": "WT1",
    # Neural
    "Glial fibrillary acidic protein": "GFAP", "SRY-box 2": "SOX2",
    "Enolase 2": "ENO2", "Recoverin": "RCVRN",
    "Glypican-3": "GPC3", "Glypican 3": "GPC3",
    # Common self-antigens
    "Myelin basic protein": "MBP", "Myelin-oligodendrocyte glycoprotein": "MOG",
    "Glutamate decarboxylase": "GAD2", "Insulin": "INS",
    "Thyroglobulin": "TG", "Desmoglein-3": "DSG3",
    "Aquaporin-4": "AQP4", "Collagen alpha-1(II)": "COL2A1",
    # Signaling / oncology
    "GTPase KRas": "KRAS", "Epidermal growth factor receptor": "EGFR",
    "Receptor tyrosine-protein kinase erbB-2": "ERBB2", "HER2": "ERBB2",
    "Mucin-1": "MUC1", "Alpha-fetoprotein": "AFP",
    "Prostate-specific antigen": "KLK3", "Mesothelin": "MSLN",
    "Telomerase reverse transcriptase": "TERT",
    "Carcinoembryonic antigen": "CEACAM5", "Survivin": "BIRC5",
    # Housekeeping
    "Actin, cytoplasmic 1": "ACTB", "Lamin-B1": "LMNB1",
    "Heat shock 70 kDa protein": "HSPA1A", "Heat shock protein": "HSP90AA1",
    "Vimentin": "VIM", "Claudin-18": "CLDN18",
}


def _protein_to_gene(protein_name: str) -> str:
    """Map IEDB protein/antigen names to gene symbols."""
    # Direct match
    if protein_name in _PROTEIN_TO_GENE:
        return _PROTEIN_TO_GENE[protein_name]
    # Case-insensitive match
    lower = protein_name.lower()
    for k, v in _PROTEIN_TO_GENE.items():
        if k.lower() == lower:
            return v
    # Partial match (protein name starts with known key)
    for k, v in _PROTEIN_TO_GENE.items():
        if lower.startswith(k.lower()):
            return v
    # If short enough and looks like a gene symbol already, keep it
    if len(protein_name) <= 10 and protein_name.isalnum():
        return protein_name.upper()
    # Fallback: truncate
    return protein_name[:20]


def _record_to_entry(record: Dict) -> Optional[DecoyEntry]:
    """
    Convert a single IEDB T-cell/epitope record to a DecoyEntry.
    Returns None if the record doesn't meet minimum quality criteria.
    """
    # --- Extract peptide sequence ---
    seq = (record.get("linear_sequence")
           or record.get("epitope_linear_sequence")
           or record.get("description", ""))
    if isinstance(seq, str):
        seq = seq.strip().upper()
    else:
        return None

    if not seq or len(seq) < 8 or len(seq) > 15:
        return None
    if not _AA_RE.match(seq):
        return None

    # --- Extract HLA allele ---
    # tcell_search has singular fields, epitope_search has plural (list) fields
    allele_raw = (record.get("mhc_allele_name")
                  or record.get("mhc_restriction")
                  or record.get("allele_name")
                  or "")
    if isinstance(allele_raw, list):
        allele_raw = allele_raw[0] if allele_raw else ""
    # Handle plural field from epitope_search
    if not allele_raw:
        allele_names = record.get("mhc_allele_names") or []
        if allele_names:
            allele_raw = allele_names[0] if isinstance(allele_names, list) else str(allele_names)
    hla = _normalise_hla(str(allele_raw))
    if not hla:
        return None

    # --- Source protein / gene ---
    source_antigen = record.get("antigen_name") or record.get("parent_source_antigen_name") or ""
    if not source_antigen:
        # Plural field from epitope_search
        names = record.get("parent_source_antigen_names") or []
        source_antigen = names[0] if names else "Unknown"
    source_antigen = str(source_antigen)
    # Clean UniProt prefix: "Titin (UniProt:Q8WZ42)" → "Titin"
    source_antigen_clean = re.sub(r"\s*\(UniProt:[A-Z0-9]+\)\s*$", "", source_antigen).strip()
    if not source_antigen_clean:
        source_antigen_clean = source_antigen

    source_organism = record.get("source_organism_name") or ""
    if not source_organism:
        org_names = record.get("source_organism_names") or []
        source_organism = org_names[0] if org_names else ""
    source_organism = str(source_organism)

    # --- Extract UniProt accession ---
    uniprot_id = None
    antigen_iri = record.get("parent_source_antigen_iri") or ""
    if isinstance(antigen_iri, str) and antigen_iri.startswith("UNIPROT:"):
        uniprot_id = antigen_iri.split(":", 1)[1].strip()
    if not uniprot_id:
        # Try extracting from antigen name: "Titin (UniProt:Q8WZ42)"
        m_up = re.search(r"\(UniProt:([A-Z0-9]+)\)", source_antigen)
        if m_up:
            uniprot_id = m_up.group(1)

    # --- Gene symbol ---
    gene_symbol = str(record.get("gene_symbol") or "UNKNOWN")
    if gene_symbol == "UNKNOWN":
        gene_symbol = _protein_to_gene(source_antigen_clean)

    # --- Assay / evidence info ---
    assay_type = str(record.get("assay_type", "IEDB record"))
    qualitative = str(record.get("qualitative_measure", ""))
    pmids = []
    ref_id = record.get("pubmed_id") or record.get("reference_id")
    if ref_id:
        pmids = [str(ref_id)]

    evidence_level = _determine_evidence_level(record)

    # --- Mass-spec ---
    ms_confirmed = None
    if "elution" in assay_type.lower() or "mass spec" in assay_type.lower():
        ms_confirmed = True

    # --- Build entry ---
    try:
        entry = DecoyEntry(
            decoy_id="DC-0000",
            peptide_info=PeptideInfo(
                decoy_sequence=seq,
                hla_allele=hla,
                source_protein=source_antigen_clean,
                gene_symbol=gene_symbol if gene_symbol != "UNKNOWN" else source_antigen_clean[:30],
                uniprot_id=uniprot_id,
            ),
            discovery_context=DiscoveryContext(
                original_target_sequence=None,
                original_target_protein=None,
                tcr_name_or_id=str(record.get("receptor_name", "")) or None,
            ),
            risk_profile=RiskProfile(
                evidence_level=evidence_level,
                critical_organs_affected=[],
            ),
            experimental_evidence=ExperimentalEvidence(
                mass_spec_confirmed=ms_confirmed,
                assays_performed=[assay_type] if assay_type else [],
                cross_reactivity_affinity=qualitative if qualitative else None,
            ),
            provenance=Provenance(
                pmid=pmids,
                evidence_summary=f"IEDB record: {assay_type}. Organism: {source_organism}. "
                                 f"Antigen: {source_antigen_clean}.",
            ),
            thought_process=None,
            validation_flags={
                "extraction_method": "iedb_miner",
                "iedb_match": "FOUND",
                "source_text_match": "N/A",
            },
        )
        return entry
    except Exception as exc:
        log.debug("Failed to build DecoyEntry from IEDB record: %s", exc)
        return None


# ── Mining strategies ─────────────────────────────────────────────────

# Proteins known to cause or risk off-target TCR toxicity
DANGER_PROTEINS = [
    # Cardiac / muscle (fatal cross-reactivity risk)
    "titin", "TTN", "MYH6", "MYH7", "TNNT2", "TNNI3", "MYL2",
    "MYBPC3", "DES", "RYR2", "ACTC1", "TPM1", "ACTN2", "DMD",
    # Melanocyte differentiation antigens
    "MLANA", "MART-1", "gp100", "PMEL", "TYR", "tyrosinase",
    "DCT", "TYRP1", "SLC45A2",
    # Cancer-testis with known cross-reactivity
    "MAGEA3", "MAGEA4", "MAGEA1", "CTAG1B", "NY-ESO-1", "SSX2",
    "PRAME", "WT1", "LAGE-1",
    # Neural antigens
    "GFAP", "SOX2", "SOX10", "OLIG2", "ENO2", "MAP2",
    "recoverin", "CRMP5", "HuD",
    # Autoimmune targets
    "MBP", "MOG", "PLP1", "GAD65", "insulin", "thyroglobulin",
    "desmoglein", "aquaporin-4", "collagen II",
]

# HLA alleles most commonly used in TCR-T therapy
COMMON_HLA = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*35:01",
    "HLA-B*44:02", "HLA-B*57:01", "HLA-C*07:02",
]


def _mine_strategy_protein(existing_seqs: Set[str]) -> List[DecoyEntry]:
    """
    Strategy 1: For each known danger protein, pull ALL T-cell assay records.
    Uses PostgREST ``ilike`` operator for case-insensitive pattern matching.
    """
    entries = []
    for protein in DANGER_PROTEINS:
        log.info("IEDB mining protein: %s", protein)
        # tcell_search uses singular field names
        results = _query_tcell({
            "parent_source_antigen_name": f"ilike.*{protein}*",
            "limit": "500",
        })
        if not results:
            # epitope_search uses plural field names wrapped in arrays
            results = _query_epitope({
                "linear_sequence": f"neq.",  # non-empty
                "linear_sequence_length": "gte.8",
            })
            # Filter client-side for antigen name
            results = [r for r in results
                       if protein.lower() in str(r.get("parent_source_antigen_names", "")).lower()
                       or protein.lower() in str(r.get("curated_source_antigens", "")).lower()]

        for r in results:
            entry = _record_to_entry(r)
            if entry and entry.peptide_info.decoy_sequence not in existing_seqs:
                entries.append(entry)
                existing_seqs.add(entry.peptide_info.decoy_sequence)

    log.info("Strategy protein: %d new entries from %d danger proteins",
             len(entries), len(DANGER_PROTEINS))
    return entries


def _mine_strategy_selfantigen(existing_seqs: Set[str]) -> List[DecoyEntry]:
    """
    Strategy 2: Mine human self-antigen epitopes that bind common HLA alleles.
    These are the pool of potential off-target peptides for TCR-T therapies.
    """
    entries = []
    for hla in COMMON_HLA:
        log.info("IEDB mining self-antigens for %s", hla)
        results = _query_tcell({
            "mhc_allele_name": f"eq.{hla}",
            "source_organism_name": f"ilike.*Homo sapiens*",
            "qualitative_measure": "eq.Positive",
            "limit": "500",
        })
        for r in results:
            entry = _record_to_entry(r)
            if entry and entry.peptide_info.decoy_sequence not in existing_seqs:
                entries.append(entry)
                existing_seqs.add(entry.peptide_info.decoy_sequence)

    log.info("Strategy self-antigen: %d new entries across %d HLA alleles",
             len(entries), len(COMMON_HLA))
    return entries


def _mine_strategy_massspec(existing_seqs: Set[str]) -> List[DecoyEntry]:
    """
    Strategy 3: Mine mass-spec / MHC-elution confirmed peptides from human
    source organisms.  These are the highest-confidence naturally presented
    peptides.
    """
    entries = []
    log.info("IEDB mining mass-spec confirmed human peptides")

    # Query for human source epitopes with elution data
    for hla in COMMON_HLA[:5]:  # Top 5 most common HLAs
        log.info("  Mass-spec mining for %s", hla)
        results = _query_epitope({
            "source_organism_names": f"ilike.*Homo sapiens*",
            "mhc_allele_names": f"ilike.*{hla.replace('*', '%')}*",
            "linear_sequence_length": "gte.8",
            "limit": "200",
        })

        for r in results:
            # Only keep if there's elution evidence
            elution_ids = r.get("elution_ids") or r.get("elution_iris") or []
            if not elution_ids:
                continue

            entry = _record_to_entry(r)
            if entry and entry.peptide_info.decoy_sequence not in existing_seqs:
                entry.experimental_evidence.mass_spec_confirmed = True
                entries.append(entry)
                existing_seqs.add(entry.peptide_info.decoy_sequence)

    log.info("Strategy mass-spec: %d new entries", len(entries))
    return entries


def _mine_strategy_crossreactive(existing_seqs: Set[str]) -> List[DecoyEntry]:
    """
    Strategy 4: Search for T-cell assays explicitly annotated as
    cross-reactive or involving autoimmune diseases.
    """
    entries = []
    cross_reactive_queries = [
        {"disease_names": "ilike.*autoimmun*"},
        {"disease_names": "ilike.*myocarditis*"},
        {"disease_names": "ilike.*encephalomyelitis*"},
        {"disease_names": "ilike.*type 1 diabetes*"},
        {"disease_names": "ilike.*vitiligo*"},
        {"disease_names": "ilike.*graft versus host*"},
        {"disease_names": "ilike.*multiple sclerosis*"},
        {"disease_names": "ilike.*rheumatoid arthritis*"},
        {"disease_names": "ilike.*celiac*"},
    ]

    for params in cross_reactive_queries:
        log.info("IEDB mining cross-reactive: %s", params)
        # Use epitope_search for disease-annotated records (plural fields)
        results = _query_epitope({**params, "limit": "200"})
        for r in results:
            # Only keep human-HLA restricted
            allele_names = r.get("mhc_allele_names") or []
            has_hla = any("HLA" in str(a) for a in allele_names)
            if not has_hla:
                continue
            entry = _record_to_entry(r)
            if entry and entry.peptide_info.decoy_sequence not in existing_seqs:
                entries.append(entry)
                existing_seqs.add(entry.peptide_info.decoy_sequence)

    # Also search tcell_search for Homo sapiens source with autoimmune diseases
    for disease in ["autoimmune", "myocarditis", "diabetes", "vitiligo", "multiple sclerosis"]:
        log.info("IEDB mining tcell cross-reactive: %s", disease)
        results = _query_tcell({
            "disease_names": f"ilike.*{disease}*",
            "mhc_allele_name": "ilike.*HLA*",
            "qualitative_measure": "eq.Positive",
            "limit": "200",
        })
        for r in results:
            entry = _record_to_entry(r)
            if entry and entry.peptide_info.decoy_sequence not in existing_seqs:
                entries.append(entry)
                existing_seqs.add(entry.peptide_info.decoy_sequence)

    log.info("Strategy cross-reactive: %d new entries", len(entries))
    return entries


# ── Main entry point ──────────────────────────────────────────────────

def mine_iedb(
    lib: DecoyLibrary,
    strategy: str = "all",
    validate: bool = True,
) -> List[DecoyEntry]:
    """
    Mine IEDB for decoy-relevant peptide entries.

    Parameters
    ----------
    lib : DecoyLibrary
        Existing library (used for deduplication).
    strategy : str
        Mining strategy: 'protein', 'selfantigen', 'massspec',
        'crossreactive', or 'all'.
    validate : bool
        Whether to run UniProt validation on mined entries.

    Returns
    -------
    list[DecoyEntry]
        New entries ready to be added to the library.
    """
    existing_seqs = {e.peptide_info.decoy_sequence for e in lib.entries}
    all_entries: List[DecoyEntry] = []

    log.info("=" * 60)
    log.info("IEDB MINER START (strategy=%s)", strategy)
    log.info("Existing library: %d entries, %d unique sequences",
             len(lib.entries), len(existing_seqs))
    log.info("=" * 60)

    if strategy in ("protein", "all"):
        all_entries.extend(_mine_strategy_protein(existing_seqs))

    if strategy in ("selfantigen", "all"):
        all_entries.extend(_mine_strategy_selfantigen(existing_seqs))

    if strategy in ("massspec", "all"):
        all_entries.extend(_mine_strategy_massspec(existing_seqs))

    if strategy in ("crossreactive", "all"):
        all_entries.extend(_mine_strategy_crossreactive(existing_seqs))

    log.info("IEDB miner total: %d new candidate entries", len(all_entries))

    # Validate if requested
    if validate and all_entries:
        from .validator import validate_batch
        log.info("Validating %d IEDB-mined entries...", len(all_entries))
        all_entries = validate_batch(all_entries)

    return all_entries
