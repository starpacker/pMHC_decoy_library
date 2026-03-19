"""
Validator Agent
===============
Cross-checks extracted DecoyEntry records against external databases:

    1. UniProt - verifies gene_symbol <-> uniprot_id mapping
    2. IEDB   - checks if the peptide exists in the Immune Epitope Database

Public API
----------
    validate_entry(entry)   -> DecoyEntry  (with validation_flags set)
    validate_batch(entries)  -> list[DecoyEntry]
"""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, List, Optional

import requests

from .config import (
    IEDB_DELAY_SEC,
    IEDB_EPITOPE_SEARCH,
    UNIPROT_DELAY_SEC,
    UNIPROT_SEARCH,
)
from .models import DecoyEntry

log = logging.getLogger(__name__)


def _query_uniprot(gene_symbol: str, organism_id: int = 9606) -> Optional[Dict[str, Any]]:
    """Search UniProt for a human protein by gene symbol."""
    time.sleep(UNIPROT_DELAY_SEC)
    params = {
        "query": f"gene:{gene_symbol} AND organism_id:{organism_id}",
        "fields": "accession,gene_names,protein_name",
        "format": "json",
        "size": "3",
    }
    try:
        resp = requests.get(UNIPROT_SEARCH, params=params, timeout=30)
        resp.raise_for_status()
        results = resp.json().get("results", [])
        for r in results:
            if r.get("entryType", "").startswith("UniProtKB reviewed"):
                return r
        return results[0] if results else None
    except Exception as exc:
        log.warning("UniProt query failed for %s: %s", gene_symbol, exc)
        return None


def _query_iedb(sequence: str) -> Optional[Dict[str, Any]]:
    """Search IEDB for a peptide by exact linear sequence."""
    time.sleep(IEDB_DELAY_SEC)
    params = {"linear_sequence": f"eq.{sequence}"}
    try:
        resp = requests.get(IEDB_EPITOPE_SEARCH, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list) and len(data) > 0:
            return data[0]
        return None
    except Exception as exc:
        log.warning("IEDB query failed for %s: %s", sequence, exc)
        return None


def validate_uniprot(entry: DecoyEntry) -> DecoyEntry:
    """Validate and enrich UniProt fields of a DecoyEntry."""
    gene = entry.peptide_info.gene_symbol
    claimed_id = entry.peptide_info.uniprot_id

    if gene in ("UNKNOWN", "", None):
        entry.validation_flags["uniprot_match"] = "SKIPPED_NO_GENE"
        return entry

    hit = _query_uniprot(gene)
    if hit is None:
        entry.validation_flags["uniprot_match"] = "NOT_FOUND"
        log.warning("UniProt: no result for gene %s", gene)
        return entry

    found_accession = hit.get("primaryAccession", "")
    found_protein = ""
    prot_desc = hit.get("proteinDescription", {})
    rec_name = prot_desc.get("recommendedName", {})
    if rec_name:
        found_protein = rec_name.get("fullName", {}).get("value", "")

    found_genes = []
    for g in hit.get("genes", []):
        gn = g.get("geneName", {}).get("value", "")
        if gn:
            found_genes.append(gn)

    gene_matches = gene.upper() in [g.upper() for g in found_genes]

    if not claimed_id:
        entry.peptide_info.uniprot_id = found_accession
        if found_protein and entry.peptide_info.source_protein.startswith("Unknown"):
            entry.peptide_info.source_protein = found_protein
        entry.validation_flags["uniprot_match"] = "ENRICHED"
        log.info("UniProt: enriched %s -> %s (%s)", gene, found_accession, found_protein)
    elif claimed_id == found_accession:
        entry.validation_flags["uniprot_match"] = "OK"
        log.info("UniProt: confirmed %s = %s", gene, found_accession)
    else:
        entry.validation_flags["uniprot_match"] = f"MISMATCH_expected_{found_accession}"
        log.warning("UniProt: MISMATCH for %s — claimed %s, found %s",
                    gene, claimed_id, found_accession)

    if not gene_matches:
        entry.validation_flags["uniprot_gene_warning"] = (
            f"Gene symbol '{gene}' not in UniProt gene list: {found_genes}"
        )

    return entry


def validate_iedb(entry: DecoyEntry) -> DecoyEntry:
    """Check if the decoy peptide exists in IEDB and enrich metadata."""
    seq = entry.peptide_info.decoy_sequence

    hit = _query_iedb(seq)
    if hit is None:
        entry.validation_flags["iedb_match"] = "NOT_FOUND"
        log.info("IEDB: no record for %s", seq)
        return entry

    entry.validation_flags["iedb_match"] = "FOUND"
    structure_id = hit.get("structure_id", "")
    entry.validation_flags["iedb_structure_id"] = str(structure_id)

    # Enrich HLA allele info
    allele_names = hit.get("mhc_allele_names", [])
    if allele_names:
        entry.validation_flags["iedb_hla_alleles"] = ", ".join(str(a) for a in allele_names if a)

    # Check mass-spec / elution evidence
    elution_ids = hit.get("elution_ids") or []
    if elution_ids and entry.experimental_evidence.mass_spec_confirmed is None:
        entry.experimental_evidence.mass_spec_confirmed = True
        entry.validation_flags["iedb_ms_confirmed"] = f"{len(elution_ids)} elution assays"
        log.info("IEDB: mass-spec confirmed for %s (%d elution assays)", seq, len(elution_ids))

    # Enrich assay information
    assay_names = hit.get("assay_names") or []
    if assay_names:
        existing = set(entry.experimental_evidence.assays_performed)
        for name in assay_names:
            if name and name not in existing:
                entry.experimental_evidence.assays_performed.append(name)

    # Check TCR receptor info
    receptor_names = hit.get("receptor_names") or []
    if receptor_names and not entry.discovery_context.tcr_name_or_id:
        entry.discovery_context.tcr_name_or_id = receptor_names[0]
        entry.validation_flags["iedb_tcr_enriched"] = ", ".join(str(r) for r in receptor_names if r)

    # Check source antigen
    parent_names = hit.get("parent_source_antigen_names") or []
    if parent_names:
        entry.validation_flags["iedb_source_antigen"] = ", ".join(str(p) for p in parent_names if p)

    # Enrich PubMed IDs from IEDB
    iedb_pmids = hit.get("pubmed_ids") or []
    if iedb_pmids:
        existing_pmids = set(entry.provenance.pmid)
        for pmid in iedb_pmids:
            if pmid and str(pmid) not in existing_pmids:
                entry.provenance.pmid.append(str(pmid))

    log.info("IEDB: found record for %s (structure_id=%s)", seq, structure_id)
    return entry


def validate_entry(entry: DecoyEntry) -> DecoyEntry:
    """
    Run all validation checks on a single DecoyEntry.

    Parameters
    ----------
    entry : DecoyEntry
        Entry to validate (modified in-place and returned).

    Returns
    -------
    DecoyEntry
        Same entry with validation_flags populated.
    """
    entry = validate_uniprot(entry)
    entry = validate_iedb(entry)

    # Compute overall validation status
    flags = entry.validation_flags
    issues = []
    if flags.get("uniprot_match", "").startswith("MISMATCH"):
        issues.append("uniprot_mismatch")
    if flags.get("uniprot_match") == "NOT_FOUND":
        issues.append("uniprot_not_found")
    if flags.get("iedb_match") == "NOT_FOUND":
        issues.append("iedb_not_found")

    if not issues:
        flags["overall_status"] = "VALIDATED"
    elif "uniprot_mismatch" in issues:
        flags["overall_status"] = "NEEDS_REVIEW"
    else:
        flags["overall_status"] = "PARTIAL"

    return entry


def validate_batch(entries: List[DecoyEntry]) -> List[DecoyEntry]:
    """Validate a list of DecoyEntry records."""
    validated = []
    for i, entry in enumerate(entries):
        log.info("Validating entry %d/%d: %s", i + 1, len(entries),
                 entry.peptide_info.decoy_sequence)
        try:
            validated.append(validate_entry(entry))
        except Exception as exc:
            log.error("Validation failed for %s: %s",
                      entry.peptide_info.decoy_sequence, exc)
            entry.validation_flags["overall_status"] = "ERROR"
            entry.validation_flags["error"] = str(exc)
            validated.append(entry)
    return validated
