#!/usr/bin/env python3
"""Split decoy_library.json into strict off-target pMHC vs broad reference vs removed."""
import json
from collections import Counter
from datetime import datetime
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "decoy_c"
LIB = DATA / "decoy_library.json"

lib = json.loads(LIB.read_text())
entries = lib["entries"]

strict = []    # 严格脱靶毒性 pMHC
broad = []     # 广义交叉反应参考
removed = []   # 不符合的 (neoantigen etc.)

for e in entries:
    cat = e.get("discovery_context", {}).get("category", "unknown")
    if cat == "unknown":
        strict.append(e)
    elif cat in ("cardiac_muscle", "autoimmune_self", "tissue_specific"):
        strict.append(e)
    elif cat == "neoantigen":
        removed.append(e)
    else:
        broad.append(e)

print(f"=== Classification Results ===")
print(f"STRICT off-target toxicity pMHC: {len(strict)}")
print(f"BROAD cross-reactivity reference: {len(broad)}")
print(f"REMOVED (not off-target): {len(removed)}")
print(f"Total: {len(strict) + len(broad) + len(removed)}")

# Breakdown
for label, lst in [("Strict", strict), ("Broad", broad), ("Removed", removed)]:
    cats = Counter(e.get("discovery_context", {}).get("category", "unknown") for e in lst)
    print(f"\n{label} breakdown:")
    for cat, cnt in cats.most_common():
        print(f"  {cat}: {cnt}")

# --- Save strict to decoy_library.json ---
now = datetime.now().isoformat()
strict_lib = {
    "version": lib.get("version", "0.3.0"),
    "metadata": {
        "description": "Strict off-target toxicity pMHC targets for TCR-T safety screening",
        "total_entries": len(strict),
        "last_updated": now,
        "categories": ["seed_curated", "cardiac_muscle", "autoimmune_self", "tissue_specific"],
        "note": "Only clinically known or high-confidence off-target toxicity pMHC. Broad cross-reactivity reference epitopes are in cross_reactivity_reference.json."
    },
    "entries": strict,
}

# Backup original
backup = DATA / "decoy_library_full_backup.json"
backup.write_text(json.dumps(lib, indent=2, ensure_ascii=False))
print(f"\nBackup saved: {backup} ({len(entries)} entries)")

# Write strict
LIB.write_text(json.dumps(strict_lib, indent=2, ensure_ascii=False))
print(f"Strict library saved: {LIB} ({len(strict)} entries)")

# --- Save broad to cross_reactivity_reference.json ---
broad_lib = {
    "version": "0.3.0",
    "metadata": {
        "description": "Broad cross-reactivity reference epitopes (viral, cancer, drug hypersensitivity, mHAg, infectious disease)",
        "total_entries": len(broad),
        "last_updated": now,
        "categories": list(set(e.get("discovery_context", {}).get("category", "unknown") for e in broad)),
        "note": "These are known T-cell epitopes that could potentially cause cross-reactivity but lack direct evidence of off-target toxicity in TCR-T therapy."
    },
    "entries": broad,
}
broad_path = DATA / "cross_reactivity_reference.json"
broad_path.write_text(json.dumps(broad_lib, indent=2, ensure_ascii=False))
print(f"Broad reference saved: {broad_path} ({len(broad)} entries)")

# --- Save removed to removed_entries.json (for audit) ---
removed_lib = {
    "version": "0.3.0",
    "metadata": {
        "description": "Removed entries - not off-target toxicity targets (e.g. neoantigens are therapeutic targets)",
        "total_entries": len(removed),
        "last_updated": now,
        "reason": "Neoantigens (KRAS G12D, TP53, BRAF V600E etc.) are therapeutic targets, not off-target toxicity pMHC"
    },
    "entries": removed,
}
removed_path = DATA / "removed_entries.json"
removed_path.write_text(json.dumps(removed_lib, indent=2, ensure_ascii=False))
print(f"Removed entries saved: {removed_path} ({len(removed)} entries)")

print("\n=== DONE ===")
print(f"decoy_library.json: {len(strict)} strict off-target pMHC")
print(f"cross_reactivity_reference.json: {len(broad)} broad reference epitopes")
print(f"removed_entries.json: {len(removed)} removed (neoantigens)")
