"""
Orchestrator
============
Coordinates the Fetcher -> Extractor -> Validator pipeline and manages
the persistent DecoyLibrary JSON database.

Public API
----------
    run_cold_start()                   -> DecoyLibrary
    run_pipeline(pmids=..., query=...) -> DecoyLibrary
    load_library()                     -> DecoyLibrary
    save_library(lib)                  -> None
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import List, Optional

from .config import DATA_DIR, LIBRARY_JSON
from .extractor import extract_from_paper
from .fetcher import (
    COLD_START_QUERIES,
    PaperRecord,
    fetch_abstract,
    fetch_multiple,
    search_pubmed,
)
from .models import DecoyEntry, DecoyLibrary
from .seed_data import build_seed_library
from .validator import hard_reject, validate_batch, validate_entry

log = logging.getLogger(__name__)


# -- Persistence -----------------------------------------------------------

def load_library() -> DecoyLibrary:
    """Load the library from disk, or return an empty one."""
    if LIBRARY_JSON.exists():
        raw = LIBRARY_JSON.read_text(encoding="utf-8")
        return DecoyLibrary.model_validate_json(raw)
    return DecoyLibrary()


def save_library(lib: DecoyLibrary) -> Path:
    """Persist the library to JSON on disk."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    data = lib.model_dump_json(indent=2)
    LIBRARY_JSON.write_text(data, encoding="utf-8")
    log.info("Saved library with %d entries to %s", len(lib.entries), LIBRARY_JSON)
    return LIBRARY_JSON


# -- Pipeline steps --------------------------------------------------------

def _assign_ids(lib: DecoyLibrary, entries: List[DecoyEntry]) -> List[DecoyEntry]:
    """Assign sequential DC-NNNN IDs to new entries."""
    current_max = 0
    for e in lib.entries:
        try:
            num = int(e.decoy_id.split("-")[1])
            if num > current_max:
                current_max = num
        except (IndexError, ValueError):
            pass

    for entry in entries:
        if entry.decoy_id == "DC-0000":
            current_max += 1
            entry.decoy_id = f"DC-{current_max:04d}"
    return entries


def process_papers(
    papers: List[PaperRecord],
    lib: DecoyLibrary,
    validate: bool = True,
    consensus: bool = False,
) -> List[DecoyEntry]:
    """
    Run Extractor + Validator + Hard-Rejection on a batch of papers and
    merge surviving entries into the library.

    Parameters
    ----------
    papers : list[PaperRecord]
        Papers to process.
    lib : DecoyLibrary
        Existing library (entries are appended in-place).
    validate : bool
        Whether to run UniProt/IEDB/protein-containment validation.
    consensus : bool
        Whether to use dual-extraction consensus mode.

    Returns list of newly added entries.
    """
    all_new = []
    total_rejected = 0

    for paper in papers:
        log.info("Processing PMID %s: %s", paper.pmid, paper.title[:60])

        # Extract (with optional consensus mode and source-text verification)
        extracted = extract_from_paper(paper, consensus=consensus)
        if not extracted:
            log.info("  No decoy entries extracted from PMID %s", paper.pmid)
            continue

        log.info("  Extracted %d candidate entries", len(extracted))

        # Assign IDs
        extracted = _assign_ids(lib, extracted)

        # Validate
        if validate:
            extracted = validate_batch(extracted)

        # Hard-rejection gate + deduplicate and add
        added = 0
        rejected = 0
        for entry in extracted:
            if validate and hard_reject(entry):
                rejected += 1
                log.info("  REJECTED: %s (%s) — %s",
                         entry.peptide_info.decoy_sequence,
                         entry.peptide_info.gene_symbol,
                         entry.validation_flags.get("overall_status", "?"))
                continue

            if lib.add_entry(entry, deduplicate=True):
                all_new.append(entry)
                added += 1
            else:
                log.info("  Skipped duplicate: %s", entry.peptide_info.decoy_sequence)

        total_rejected += rejected
        log.info("  Added %d, rejected %d from PMID %s", added, rejected, paper.pmid)

    if total_rejected:
        log.info("Total hard-rejected across batch: %d entries", total_rejected)

    return all_new


# -- High-level workflows --------------------------------------------------

def run_cold_start(validate: bool = True) -> DecoyLibrary:
    """
    Initialize the library with seed data + automated PubMed discovery.

    Steps:
        1. Load seed data (hand-curated high-confidence entries).
        2. Validate seed entries against UniProt + IEDB.
        3. Search PubMed with cold-start queries.
        4. Fetch and process discovered papers.
        5. Save the resulting library.

    Returns
    -------
    DecoyLibrary
        The initialized library.
    """
    log.info("=== COLD START: Initializing Decoy C Library ===")

    # Step 1: Seed data
    lib = build_seed_library()
    log.info("Loaded %d seed entries", len(lib.entries))

    # Step 2: Validate seed entries
    if validate:
        log.info("Validating seed entries...")
        lib.entries = validate_batch(lib.entries)

    # Step 3-4: Discover and process papers
    all_pmids = set()
    # Known important PMIDs to always include
    known_pmids = {"23926201", "23863783", "24475783", "26457759", "19451549", "21282551"}
    all_pmids.update(known_pmids)

    for query in COLD_START_QUERIES:
        try:
            pmids = search_pubmed(query, retmax=10)
            all_pmids.update(pmids)
            log.info("Query '%s' returned %d PMIDs", query[:40], len(pmids))
        except Exception as exc:
            log.error("Search failed for '%s': %s", query[:40], exc)

    # Remove PMIDs already in seed data
    seed_pmids = set()
    for e in lib.entries:
        seed_pmids.update(e.provenance.pmid)
    new_pmids = all_pmids - seed_pmids

    log.info("Discovered %d total PMIDs, %d new after dedup", len(all_pmids), len(new_pmids))

    if new_pmids:
        papers = fetch_multiple(list(new_pmids)[:30])  # Limit to 30 for cold start
        new_entries = process_papers(papers, lib, validate=validate)
        log.info("Added %d new entries from PubMed discovery", len(new_entries))

    # Step 5: Save
    save_library(lib)
    log.info("=== COLD START COMPLETE: %d total entries ===", len(lib.entries))
    return lib


def run_pipeline(
    pmids: Optional[List[str]] = None,
    query: Optional[str] = None,
    validate: bool = True,
) -> DecoyLibrary:
    """
    Run the extraction pipeline on specific PMIDs or a PubMed query.

    Parameters
    ----------
    pmids : list[str], optional
        Specific PMIDs to process.
    query : str, optional
        PubMed search query to discover papers.
    validate : bool
        Whether to run validation.

    Returns
    -------
    DecoyLibrary
        Updated library.
    """
    lib = load_library()
    if not lib.entries:
        log.info("Library is empty - running cold start first")
        return run_cold_start(validate=validate)

    target_pmids = set()
    if pmids:
        target_pmids.update(pmids)
    if query:
        found = search_pubmed(query, retmax=20)
        target_pmids.update(found)

    if not target_pmids:
        log.warning("No PMIDs to process")
        return lib

    papers = fetch_multiple(list(target_pmids))
    new_entries = process_papers(papers, lib, validate=validate)
    
    if new_entries:
        save_library(lib)
        log.info("Pipeline complete: added %d new entries (total: %d)",
                 len(new_entries), len(lib.entries))
    else:
        log.info("Pipeline complete: no new entries added")

    return lib
