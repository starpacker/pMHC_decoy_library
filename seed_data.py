"""
Seed Data for Cold Start
========================
Loads curated decoy entries from seed_entries.json.
Covers 50+ documented TCR cross-reactivity / off-target peptide cases
from published literature, each with full source attribution.
"""

import json
import os
from typing import List

from .models import DecoyEntry, DecoyLibrary, Source


_SEED_FILE = os.path.join(os.path.dirname(__file__), "data", "seed_entries.json")


def _load_seed_entries() -> List[DecoyEntry]:
    """Load seed entries from the JSON data file."""
    with open(_SEED_FILE, "r") as f:
        raw = json.load(f)

    entries = []
    for item in raw:
        # Build Source if present
        src = None
        if "source" in item and item["source"]:
            src = Source(**item["source"])

        entry = DecoyEntry(
            decoy_id=item["decoy_id"],
            peptide_info=item["peptide_info"],
            discovery_context=item.get("discovery_context", {}),
            risk_profile=item["risk_profile"],
            experimental_evidence=item.get("experimental_evidence", {}),
            provenance=item.get("provenance", {}),
            source=src,
            thought_process=item.get("thought_process"),
        )
        entries.append(entry)
    return entries


def build_seed_library() -> DecoyLibrary:
    """Construct the initial DecoyLibrary with curated seed entries."""
    entries = _load_seed_entries()
    lib = DecoyLibrary()
    for entry in entries:
        lib.add_entry(entry)
    return lib
