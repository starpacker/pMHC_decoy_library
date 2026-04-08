#!/usr/bin/env python3
"""
Post-processing: fix gene_symbol / uniprot_id for IEDB-mined entries,
then re-run protein containment validation.

This addresses the issue where IEDB entries had full protein names
as gene_symbol (e.g., "Glypican-3") instead of proper gene symbols
(e.g., "GPC3"), and uniprot_id was None despite IEDB providing it.

Usage:
    python -m decoy_c.fix_and_revalidate
"""

from __future__ import annotations

import json
import logging
import re
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from decoy_c.config import DATA_DIR, LIBRARY_JSON
from decoy_c.iedb_miner import _protein_to_gene
from decoy_c.models import DecoyEntry, DecoyLibrary
from decoy_c.orchestrator import save_library
from decoy_c.validator import verify_peptide_in_protein

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("fix_revalidate")


def fix_gene_and_uniprot(lib: DecoyLibrary) -> int:
    """
    For each entry missing uniprot_id, try to extract it from
    validation_flags['iedb_source_antigen'] or provenance fields.
    Also fix gene_symbol using the protein-to-gene mapping.
    """
    fixed = 0
    for entry in lib.entries:
        pi = entry.peptide_info
        flags = entry.validation_flags

        # Already has uniprot_id and confirmed protein_containment? Skip.
        if pi.uniprot_id and flags.get("protein_containment") == "CONFIRMED":
            continue

        # Try to extract uniprot_id from iedb_source_antigen field
        # Format: "Titin (UniProt:Q8WZ42)" or "Actin, cytoplasmic 1"
        source_antigen = flags.get("iedb_source_antigen", "") or pi.source_protein or ""
        if not pi.uniprot_id:
            m = re.search(r"\(UniProt:([A-Z0-9]+)\)", source_antigen)
            if m:
                pi.uniprot_id = m.group(1)
                fixed += 1
                log.info("Fixed uniprot_id for %s: %s (from iedb_source_antigen)",
                         pi.decoy_sequence, pi.uniprot_id)

        # Try from provenance evidence_summary
        if not pi.uniprot_id:
            m = re.search(r"UniProt:([A-Z0-9]+)", entry.provenance.evidence_summary)
            if m:
                pi.uniprot_id = m.group(1)
                fixed += 1

        # Fix gene_symbol if it looks like a full protein name
        if pi.gene_symbol and len(pi.gene_symbol) > 10:
            new_gene = _protein_to_gene(pi.gene_symbol)
            if new_gene != pi.gene_symbol:
                log.info("Fixed gene_symbol: %s -> %s", pi.gene_symbol, new_gene)
                pi.gene_symbol = new_gene

        # Also fix source_protein gene_symbol
        if pi.gene_symbol and len(pi.gene_symbol) > 10:
            new_gene = _protein_to_gene(pi.source_protein)
            if new_gene != pi.source_protein and len(new_gene) <= 10:
                pi.gene_symbol = new_gene

    return fixed


def revalidate_protein_containment(lib: DecoyLibrary) -> dict:
    """Re-run protein containment for entries that now have uniprot_id."""
    stats = {"checked": 0, "confirmed": 0, "not_found": 0, "error": 0, "skipped": 0}

    for entry in lib.entries:
        flags = entry.validation_flags
        uid = entry.peptide_info.uniprot_id

        # Skip if already confirmed
        if flags.get("protein_containment") == "CONFIRMED":
            continue

        # Skip if no uniprot_id
        if not uid:
            stats["skipped"] += 1
            continue

        stats["checked"] += 1
        entry = verify_peptide_in_protein(entry)

        result = flags.get("protein_containment", "ERROR")
        if result == "CONFIRMED":
            stats["confirmed"] += 1
        elif result == "NOT_FOUND":
            stats["not_found"] += 1
        else:
            stats["error"] += 1

        if stats["checked"] % 20 == 0:
            log.info("  Progress: %d checked, %d confirmed, %d not_found",
                     stats["checked"], stats["confirmed"], stats["not_found"])

    return stats


def drop_unverified_iedb_entries(lib: DecoyLibrary) -> int:
    """
    Drop IEDB-mined entries where protein_containment is NOT_FOUND.
    These are peptides that IEDB lists but do NOT appear in the claimed protein.
    """
    before = len(lib.entries)
    kept = []
    dropped = []

    for entry in lib.entries:
        flags = entry.validation_flags
        if (flags.get("extraction_method") == "iedb_miner"
                and flags.get("protein_containment") == "NOT_FOUND"):
            dropped.append(entry)
            log.warning("Dropping %s: protein_containment=NOT_FOUND (uniprot=%s)",
                        entry.peptide_info.decoy_sequence,
                        entry.peptide_info.uniprot_id)
        else:
            kept.append(entry)

    lib.entries = kept

    if dropped:
        dropped_path = DATA_DIR / "dropped_protein_containment.json"
        dropped_data = {
            "description": "IEDB entries dropped because peptide not found in claimed protein",
            "count": len(dropped),
            "entries": [json.loads(e.model_dump_json()) for e in dropped],
        }
        dropped_path.write_text(json.dumps(dropped_data, indent=2, ensure_ascii=False), encoding="utf-8")
        log.info("Saved %d dropped entries to %s", len(dropped), dropped_path)

    return before - len(lib.entries)


def print_summary(lib: DecoyLibrary):
    entries = lib.entries
    log.info("")
    log.info("=" * 60)
    log.info("FINAL LIBRARY SUMMARY")
    log.info("=" * 60)
    log.info("Total entries: %d", len(entries))

    seqs = set(e.peptide_info.decoy_sequence for e in entries)
    log.info("Unique sequences: %d", len(seqs))

    protein = Counter(e.validation_flags.get("protein_containment", "NOT_CHECKED") for e in entries)
    log.info("\nProtein containment:")
    for p, cnt in sorted(protein.items()):
        log.info("  %-15s %d", p, cnt)

    iedb = Counter(e.validation_flags.get("iedb_match", "NOT_CHECKED") for e in entries)
    log.info("\nIEDB match:")
    for i, cnt in sorted(iedb.items()):
        log.info("  %-15s %d", i, cnt)

    methods = Counter(e.validation_flags.get("extraction_method", "llm") for e in entries)
    log.info("\nExtraction methods:")
    for m, cnt in sorted(methods.items()):
        log.info("  %-25s %d", m, cnt)

    ms = sum(1 for e in entries if e.experimental_evidence.mass_spec_confirmed)
    log.info("\nMass-spec confirmed: %d", ms)

    hlas = set(e.peptide_info.hla_allele for e in entries)
    genes = set(e.peptide_info.gene_symbol for e in entries)
    log.info("Unique HLA alleles: %d", len(hlas))
    log.info("Unique genes: %d", len(genes))


def main():
    start = time.time()

    log.info("Loading library from %s", LIBRARY_JSON)
    raw = json.loads(LIBRARY_JSON.read_text(encoding="utf-8"))
    lib = DecoyLibrary(entries=[DecoyEntry.model_validate(e) for e in raw.get("entries", [])])
    log.info("Loaded %d entries", len(lib.entries))

    # Phase 1: Fix gene_symbol and uniprot_id
    log.info("")
    log.info("=" * 60)
    log.info("PHASE 1: FIX GENE_SYMBOL AND UNIPROT_ID")
    log.info("=" * 60)
    fixed = fix_gene_and_uniprot(lib)
    log.info("Fixed %d entries", fixed)

    # Phase 2: Re-validate protein containment
    log.info("")
    log.info("=" * 60)
    log.info("PHASE 2: RE-VALIDATE PROTEIN CONTAINMENT")
    log.info("=" * 60)
    stats = revalidate_protein_containment(lib)
    log.info("Protein containment re-validation: %s", stats)

    # Phase 3: Drop IEDB entries that fail protein containment
    log.info("")
    log.info("=" * 60)
    log.info("PHASE 3: DROP UNVERIFIED ENTRIES")
    log.info("=" * 60)
    dropped = drop_unverified_iedb_entries(lib)
    log.info("Dropped %d entries with protein_containment=NOT_FOUND", dropped)

    # Save
    save_library(lib)
    log.info("Saved library: %d entries", len(lib.entries))

    print_summary(lib)

    elapsed = time.time() - start
    log.info("\nTotal time: %.1f minutes (%.0f seconds)", elapsed / 60, elapsed)


if __name__ == "__main__":
    main()
