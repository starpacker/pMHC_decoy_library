#!/usr/bin/env python3
"""
Decoy C — Schema Migration + Triple Revalidation + IEDB Mining
===============================================================
1. Migrate 145 old-schema entries (from generate_50.py) to standard schema
2. Re-validate ALL entries through triple-validation pipeline
3. Hard-reject entries with zero evidence
4. Mine IEDB for new high-confidence entries
5. Save the cleaned + expanded library

Usage:
    python -m decoy_c.revalidate_and_mine
"""

from __future__ import annotations

import json
import logging
import re
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from decoy_c.config import DATA_DIR, LIBRARY_JSON
from decoy_c.models import (
    DecoyEntry,
    DecoyLibrary,
    DiscoveryContext,
    EvidenceLevel,
    ExperimentalEvidence,
    PeptideInfo,
    Provenance,
    RiskProfile,
)
from decoy_c.orchestrator import save_library
from decoy_c.validator import (
    hard_reject,
    validate_entry,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("revalidate")


# ── Phase 0: Load raw JSON + migrate schema ──────────────────────────

def _migrate_risk_level(old_level: str) -> EvidenceLevel:
    """Map old risk_level (high/medium/low) to new EvidenceLevel."""
    mapping = {
        "critical": EvidenceLevel.LEVEL_1,
        "high": EvidenceLevel.LEVEL_2,
        "medium": EvidenceLevel.LEVEL_3,
        "low": EvidenceLevel.LEVEL_4,
    }
    return mapping.get(old_level, EvidenceLevel.LEVEL_4)


def _infer_gene_symbol(source_protein: str) -> str:
    """Best-effort gene symbol from source protein name."""
    known = {
        "MBP": "MBP", "MOG": "MOG", "PLP": "PLP1",
        "GAD65": "GAD2", "GAD67": "GAD1", "insulin": "INS",
        "thyroglobulin": "TG", "TPO": "TPO",
        "TTN": "TTN", "titin": "TTN", "Titin": "TTN",
        "collagen": "COL2A1", "desmoglein": "DSG3",
        "aquaporin": "AQP4", "aquaporin-4": "AQP4",
        "MAGE-A3": "MAGEA3", "MAGEA3": "MAGEA3",
        "MAGE-A4": "MAGEA4", "MAGE-A1": "MAGEA1",
        "NY-ESO-1": "CTAG1B", "CTAG1B": "CTAG1B",
        "gp100": "PMEL", "PMEL": "PMEL",
        "MLANA": "MLANA", "MART-1": "MLANA", "Melan-A": "MLANA",
        "TYR": "TYR", "tyrosinase": "TYR",
        "WT1": "WT1", "PRAME": "PRAME", "SSX2": "SSX2",
        "survivin": "BIRC5", "CEA": "CEACAM5",
        "HER2": "ERBB2", "AFP": "AFP", "PSA": "KLK3",
        "MUC1": "MUC1", "TERT": "TERT",
        "EGFR": "EGFR", "mesothelin": "MSLN",
        "DCT": "DCT", "TYRP1": "TYRP1",
        "SOX2": "SOX2", "SOX10": "SOX10",
        "GFAP": "GFAP", "ENO2": "ENO2",
        "recoverin": "RCVRN", "CRMP5": "DPYSL5",
        "MYH6": "MYH6", "MYH7": "MYH7",
        "TNNT2": "TNNT2", "TNNI3": "TNNI3",
    }
    # Direct match
    if source_protein in known:
        return known[source_protein]
    # Case-insensitive
    for k, v in known.items():
        if k.lower() == source_protein.lower():
            return v
    # First word match
    first_word = source_protein.split()[0].split("/")[0].split("(")[0].strip()
    if first_word in known:
        return known[first_word]
    return source_protein.upper()[:10]  # fallback


def load_and_migrate() -> DecoyLibrary:
    """Load the raw JSON, migrate old-schema entries to new schema."""
    raw = json.loads(LIBRARY_JSON.read_text(encoding="utf-8"))
    raw_entries = raw.get("entries", [])

    log.info("Raw JSON: %d entries", len(raw_entries))

    migrated_entries = []
    new_schema_count = 0
    old_schema_count = 0
    failed_count = 0

    for i, e in enumerate(raw_entries):
        pi = e.get("peptide_info", {})

        if "hla_allele" in pi and "gene_symbol" in pi:
            # ── New schema — try direct Pydantic parse ────────────
            try:
                entry = DecoyEntry.model_validate(e)
                new_schema_count += 1
                migrated_entries.append(entry)
                continue
            except Exception as exc:
                log.warning("Entry %d: new-schema parse failed: %s", i, exc)

        if "hla_restriction" in pi:
            # ── Old schema — migrate ──────────────────────────────
            old_schema_count += 1
            try:
                seq = pi.get("decoy_sequence", "")
                hla = pi.get("hla_restriction", "HLA-A*02:01")
                source_protein = pi.get("source_protein", "Unknown")
                gene = _infer_gene_symbol(source_protein)

                old_risk = e.get("risk_profile", {})
                evidence_level = _migrate_risk_level(old_risk.get("risk_level", "low"))

                old_ctx = e.get("discovery_context", {})
                old_prov = e.get("provenance", {})
                old_exp = e.get("experimental_evidence", {})

                entry = DecoyEntry(
                    decoy_id=e.get("decoy_id", "DC-0000"),
                    peptide_info=PeptideInfo(
                        decoy_sequence=seq,
                        hla_allele=hla,
                        source_protein=source_protein,
                        gene_symbol=gene,
                        uniprot_id=None,
                    ),
                    discovery_context=DiscoveryContext(
                        original_target_sequence=None,
                        original_target_protein=None,
                        tcr_name_or_id=None,
                    ),
                    risk_profile=RiskProfile(
                        evidence_level=evidence_level,
                        critical_organs_affected=[],
                    ),
                    experimental_evidence=ExperimentalEvidence(
                        mass_spec_confirmed=None,
                        assays_performed=[],
                    ),
                    provenance=Provenance(
                        pmid=[],
                        evidence_summary=old_ctx.get("description", "Migrated from LLM-batch generation"),
                    ),
                    thought_process=None,
                    validation_flags={
                        "extraction_method": "llm_batch_migrated",
                        "migration_source": "generate_50",
                        "needs_validation": "true",
                    },
                )
                migrated_entries.append(entry)
            except Exception as exc:
                log.warning("Entry %d: migration failed for %s: %s",
                           i, pi.get("decoy_sequence", "?"), exc)
                failed_count += 1
        else:
            log.warning("Entry %d: unknown schema, keys=%s", i, list(pi.keys()))
            failed_count += 1

    log.info("Migration results:")
    log.info("  New-schema (kept as-is): %d", new_schema_count)
    log.info("  Old-schema (migrated):   %d", old_schema_count)
    log.info("  Failed:                  %d", failed_count)
    log.info("  Total migrated:          %d", len(migrated_entries))

    lib = DecoyLibrary(entries=migrated_entries)
    return lib


# ── Phase 1: Triple revalidation ─────────────────────────────────────

def revalidate_existing(lib: DecoyLibrary) -> DecoyLibrary:
    """Re-validate all entries through triple-validation + hard-reject."""
    log.info("")
    log.info("=" * 60)
    log.info("PHASE 1: TRIPLE REVALIDATION OF %d ENTRIES", len(lib.entries))
    log.info("=" * 60)

    validated = []
    rejected = []
    errors = []

    for i, entry in enumerate(lib.entries):
        seq = entry.peptide_info.decoy_sequence
        gene = entry.peptide_info.gene_symbol
        log.info("[%d/%d] %s (%s)", i + 1, len(lib.entries), seq, gene)

        try:
            flags = entry.validation_flags

            # Check if already fully validated with our new triple checks
            already_done = (
                flags.get("protein_containment") in ("CONFIRMED", "NOT_FOUND")
                and flags.get("overall_status") in ("VALIDATED", "PARTIAL", "NEEDS_REVIEW", "REJECTED")
                and flags.get("iedb_match") in ("FOUND", "NOT_FOUND")
                and flags.get("uniprot_match") not in (None, "")
            )

            if already_done:
                log.info("  Already triple-validated, reusing existing flags")
            else:
                entry = validate_entry(entry)
                flags = entry.validation_flags

            # Apply hard rejection
            if hard_reject(entry):
                rejected.append(entry)
                log.warning("  REJECTED: %s | protein=%s | iedb=%s | source=%s | overall=%s",
                           seq,
                           flags.get("protein_containment", "?"),
                           flags.get("iedb_match", "?"),
                           flags.get("source_text_match", "?"),
                           flags.get("overall_status", "?"))
            else:
                validated.append(entry)
                log.info("  OK: protein=%s | iedb=%s | overall=%s",
                        flags.get("protein_containment", "?"),
                        flags.get("iedb_match", "?"),
                        flags.get("overall_status", "?"))

        except Exception as exc:
            log.error("  ERROR validating %s: %s", seq, exc)
            entry.validation_flags["overall_status"] = "ERROR"
            entry.validation_flags["error"] = str(exc)
            # Keep errors for now, they can be re-tried
            errors.append(entry)

    log.info("")
    log.info("=" * 40)
    log.info("RE-VALIDATION RESULTS:")
    log.info("  Passed:   %d", len(validated))
    log.info("  Rejected: %d", len(rejected))
    log.info("  Errors:   %d", len(errors))
    log.info("=" * 40)

    # Save rejected for audit
    if rejected:
        rejected_path = DATA_DIR / "rejected_by_revalidation.json"
        rejected_data = {
            "description": "Entries rejected during triple-validation revalidation",
            "date": datetime.now().isoformat(),
            "count": len(rejected),
            "entries": [json.loads(e.model_dump_json()) for e in rejected],
        }
        rejected_path.write_text(json.dumps(rejected_data, indent=2, ensure_ascii=False), encoding="utf-8")
        log.info("  Rejected entries saved to %s", rejected_path)

    lib.entries = validated + errors
    return lib


# ── Phase 2: IEDB mining ─────────────────────────────────────────────

def mine_iedb_entries(lib: DecoyLibrary) -> DecoyLibrary:
    """Mine IEDB for new entries and add to library."""
    log.info("")
    log.info("=" * 60)
    log.info("PHASE 2: IEDB SYSTEMATIC MINING")
    log.info("=" * 60)

    try:
        from decoy_c.iedb_miner import mine_iedb
    except ImportError as exc:
        log.error("Cannot import iedb_miner: %s", exc)
        return lib

    new_entries = mine_iedb(lib, strategy="all", validate=True)

    # Assign IDs and add
    added = 0
    rejected_protein = 0
    current_max = 0
    for e in lib.entries:
        try:
            num = int(e.decoy_id.split("-")[1])
            if num > current_max:
                current_max = num
        except (IndexError, ValueError):
            pass

    for entry in new_entries:
        # Extra check: reject if protein containment failed
        if entry.validation_flags.get("protein_containment") == "NOT_FOUND":
            rejected_protein += 1
            continue

        current_max += 1
        entry.decoy_id = f"DC-{current_max:04d}"

        if lib.add_entry(entry, deduplicate=True):
            added += 1

    log.info("")
    log.info("IEDB MINING RESULTS:")
    log.info("  Candidates: %d", len(new_entries))
    log.info("  Added:      %d", added)
    log.info("  Rejected (protein containment): %d", rejected_protein)

    return lib


# ── Summary ───────────────────────────────────────────────────────────

def print_summary(lib: DecoyLibrary):
    entries = lib.entries
    log.info("")
    log.info("=" * 60)
    log.info("FINAL LIBRARY SUMMARY")
    log.info("=" * 60)
    log.info("Total entries: %d", len(entries))

    seqs = set(e.peptide_info.decoy_sequence for e in entries)
    log.info("Unique sequences: %d", len(seqs))

    statuses = Counter(e.validation_flags.get("overall_status", "?") for e in entries)
    log.info("\nValidation status:")
    for s, cnt in sorted(statuses.items()):
        log.info("  %-15s %d", s, cnt)

    levels = Counter(e.risk_profile.evidence_level.value for e in entries)
    log.info("\nEvidence levels:")
    for lv, cnt in sorted(levels.items()):
        log.info("  %-35s %d", lv, cnt)

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

    hlas = set(e.peptide_info.hla_allele for e in entries)
    genes = set(e.peptide_info.gene_symbol for e in entries)
    ms = sum(1 for e in entries if e.experimental_evidence.mass_spec_confirmed)
    log.info("\nUnique HLA alleles: %d", len(hlas))
    log.info("Unique genes: %d", len(genes))
    log.info("Mass-spec confirmed: %d", ms)


# ── Main ──────────────────────────────────────────────────────────────

def main():
    start = time.time()

    # Backup
    if LIBRARY_JSON.exists():
        backup = DATA_DIR / f"decoy_library_pre_revalidation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        import shutil
        shutil.copy2(LIBRARY_JSON, backup)
        log.info("Backup: %s", backup)

    # Phase 0: Migrate
    lib = load_and_migrate()
    save_library(lib)
    log.info("Post-migration: %d entries saved", len(lib.entries))

    # Phase 1: Triple revalidation
    lib = revalidate_existing(lib)
    save_library(lib)
    log.info("Post-revalidation: %d entries saved", len(lib.entries))

    # Phase 2: IEDB mining
    lib = mine_iedb_entries(lib)
    save_library(lib)
    log.info("Post-IEDB-mining: %d entries saved", len(lib.entries))

    # Summary
    print_summary(lib)

    elapsed = time.time() - start
    log.info("\nTotal time: %.1f minutes (%.0f seconds)", elapsed / 60, elapsed)


if __name__ == "__main__":
    main()
