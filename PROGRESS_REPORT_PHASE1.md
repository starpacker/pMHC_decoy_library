# Phase 1: Interface Descriptor Upgrade — Progress Report

**Date**: 2026-04-03
**Status**: Code implementation complete, pending server deployment & validation

---

## Overview

Upgraded Decoy B scoring from pure RMSD-based to RMSD + interface descriptor fusion.

### Scoring Formula Change

```
OLD:  surface_correlation = max(avg(sim_a, sim_b), bf_corr)

NEW:  rmsd_geo           = max(avg(sim_a, sim_b), bf_corr)        # unchanged
      interface_combined = adaptive_weighted_avg(
          0.40 * plip_tanimoto,       # non-covalent interactions
          0.25 * bsa_similarity,      # buried surface area
          0.35 * prodigy_similarity,  # binding affinity
      )
      surface_correlation = 0.50 * rmsd_geo + 0.50 * interface_combined
```

Upstream formulas unchanged:
- `similarity_score = 0.4 * cos_sim + 0.6 * surface_correlation`
- `risk = similarity_score * (1/EL_Rank) * TPM_Weight`

---

## Completed Tasks

### Task 1: Interface Descriptors Module
- **File**: `decoy_b/tools/interface_descriptors.py` (NEW, ~430 lines)
- Three descriptor classes:
  - **PLIP**: Non-covalent interaction fingerprint (H-bond, hydrophobic, salt bridge, pi-stacking, pi-cation) + Tanimoto similarity
  - **FreeSASA**: Buried surface area at peptide-MHC interface + normalized similarity
  - **PRODIGY**: Contact-based binding affinity (dG) prediction + similarity
- Adaptive weighting: if a descriptor fails, its weight redistributes to successful ones
- Data structures: `InteractionFingerprint`, `InterfaceDescriptors`, `InterfaceSimilarity`

### Task 2: StructuralScore Model Extension
- **File**: `decoy_a/models.py` — `StructuralScore` class
- Added 8 Optional fields: `plip_tanimoto`, `bsa_target`, `bsa_candidate`, `bsa_similarity`, `prodigy_dg_target`, `prodigy_dg_candidate`, `prodigy_similarity`, `interface_combined`
- All backward-compatible (Optional + default None)

### Task 3: Scanner Integration
- **File**: `decoy_b/scanner.py` — `compute_structure_similarity()`
- After RMSD calculation, calls interface descriptors module
- Blends RMSD geometric score (50%) with interface combined (50%)
- Gracefully degrades: if descriptors unavailable, falls back to pure RMSD
- Tool label appends `+interface` when descriptors are active

### Task 4: Output Format & Risk Scorer Update
- **File**: `decoy_a/models.py` — `DecoyABEntry` class
  - Added 4 fields: `plip_tanimoto`, `bsa_similarity`, `prodigy_similarity`, `interface_combined`
- **File**: `decoy_b/risk_scorer.py`
  - Passes interface descriptor fields from `StructuralScore` to `DecoyABEntry`

### Task 5: Unit Tests
- **File**: `tests/test_interface_descriptors.py` (NEW, 28 tests)
- **Results**: 22/22 unit tests PASSED
  - `TestTanimoto` (7 tests): identical, zero, disjoint, similar, symmetry, vector shape, dict
  - `TestBSASimilarity` (6 tests): identical, max_diff, mid_range, symmetry, custom, clamp
  - `TestProdigySimilarity` (5 tests): identical, max_diff, mid_range, symmetry, custom
  - `TestInterfaceSimilarity` (1 test): data class to_dict
  - `TestModelBackwardCompatibility` (3 tests): old/new StructuralScore, DecoyABEntry
- 6 integration tests (skipped locally — need PDB files + external tools on server)

### Task 6: Benchmark Script
- **File**: `scripts/benchmark_interface_descriptors.py` (NEW)
- Tests each method individually: timing, score distributions, failure rates
- Combined pairwise similarity comparison
- Old vs new scoring comparison
- Generates: `figures/benchmark_interface_descriptors.png` + `figures/benchmark_report.txt`

### Task 7: Visualization Script
- **File**: `scripts/visualize_interface_descriptors.py` (NEW)
- 9-panel figure: descriptor availability, score distributions, radar profiles,
  PLIP vs BSA scatter, rank shift analysis, top-20 score breakdown

---

## File Change Summary

| File | Action | Lines |
|------|--------|-------|
| `decoy_b/tools/interface_descriptors.py` | NEW | ~430 |
| `decoy_a/models.py` | MODIFIED | +20 (StructuralScore + DecoyABEntry) |
| `decoy_b/scanner.py` | MODIFIED | +60 (interface descriptor integration) |
| `decoy_b/risk_scorer.py` | MODIFIED | +8 (pass new fields) |
| `tests/test_interface_descriptors.py` | NEW | ~250 |
| `scripts/benchmark_interface_descriptors.py` | NEW | ~320 |
| `scripts/visualize_interface_descriptors.py` | NEW | ~280 |

---

## Pending: Server Deployment

### Step 1: Install dependencies
```bash
cd /share/liuyutian/pMHC_decoy_library
conda install -c conda-forge openbabel -y
pip install plip freesasa prodigy-prot
```

### Step 2: Copy code changes
Copy the following files from local to server:
- `decoy_b/tools/interface_descriptors.py`
- `decoy_a/models.py` (modified)
- `decoy_b/scanner.py` (modified)
- `decoy_b/risk_scorer.py` (modified)
- `tests/test_interface_descriptors.py`
- `scripts/benchmark_interface_descriptors.py`
- `scripts/visualize_interface_descriptors.py`

### Step 3: Run benchmark
```bash
python scripts/benchmark_interface_descriptors.py
```

### Step 4: Run full pipeline regression
```bash
# Backup old results
cp data/decoy_b/final_ranked_decoys.json data/decoy_b/final_ranked_decoys_BACKUP.json

# Run pipeline
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# Visualize
python scripts/visualize_interface_descriptors.py
```

### Step 5: Compare old vs new
- Top 20 overlap should be >50%
- New fields should be populated in output JSON
- `modeling_tool` should show `+interface` suffix

---

## Design Decisions

1. **50/50 weight split** (RMSD vs interface): Conservative starting point. Can be tuned after collecting benchmark data.
2. **Adaptive weighting**: If PLIP fails but BSA/PRODIGY succeed, their weights redistribute proportionally. Ensures graceful degradation.
3. **PLIP 0.40 > PRODIGY 0.35 > BSA 0.25**: Non-covalent interaction patterns are the most direct TCR recognition determinant. BSA is necessary but not sufficient.
4. **Backward compatibility**: All new fields are Optional with None defaults. Existing code paths unaffected.

---

# Phase 2: Boltz-2 Cross-Validation Integration

**Date**: 2026-04-03
**Status**: Code implementation complete, pending server deployment & Boltz weights

---

## Overview

Integrated Boltz-2 as an independent structure prediction engine for cross-validation of tFold/AF3 predictions. This addresses the key reliability concern: a single model's structural similarity score may have systematic bias. By requiring agreement between two independent prediction engines, we significantly increase confidence in structural similarity assessments.

### Scoring Formula Change

```
Phase 1:  similarity_score = 0.4 * cos_sim + 0.6 * surface_correlation

Phase 2:  base_score       = 0.4 * cos_sim + 0.6 * surface_correlation
          cv_boost          = 0.10 * cross_validation_agreement    # [0, 0.10]
          similarity_score  = min(1.0, base_score + cv_boost)

Where:
  cross_validation_agreement = max(0, 1 - RMSD_primary_vs_boltz / 4.0)
  # RMSD < 0.5A → agreement ~0.87, RMSD ~2A → agreement ~0.50, RMSD > 4A → 0
```

### Why Cross-Validation?

- **tFold** is fast (~25x faster than AF3) but specialized for TCR-pMHC; may miss unusual conformations
- **Boltz-2** is a general-purpose structure predictor (MIT open-source, near-AF3 accuracy)
- Two independent models agreeing on a similar structure = high confidence
- Two models disagreeing = flag for manual review or deprioritize

---

## Completed Tasks

### Task 1: Boltz Wrapper Module
- **File**: `decoy_b/tools/boltz.py` (NEW, ~370 lines)
- Three execution modes with automatic fallback:
  1. **Python API** — `from boltz.main import predict` (pip-installed)
  2. **CLI command** — `boltz predict <yaml>` (PATH-installed)
  3. **Subprocess** — `python -m boltz predict` (from BOLTZ_DIR)
- YAML input builder for pMHC three-chain complex (A=MHC, B=B2M, C=peptide)
- Confidence metrics parser (confidence_score, pLDDT, PTM, iPTM)
- Output structure finder (mmCIF/PDB) with Boltz naming convention
- CIF-to-PDB conversion (BioPython fallback)
- Result caching: skips re-prediction if output already exists
- Environment variables: `BOLTZ_DIR`, `BOLTZ_CACHE`, `BOLTZ_MODEL`, `BOLTZ_DEVICE`

### Task 2: Cross-Validation Agreement Function
- **File**: `decoy_b/scanner.py` — `_compute_cross_validation_agreement()`
- Loads primary (tFold/AF3) and Boltz PDB structures
- Superimposes CA atoms using BioPython
- Computes backbone RMSD → converts to agreement score [0, 1]
- Graceful degradation if BioPython or structures unavailable

### Task 3: Scanner Pipeline Extension
- **File**: `decoy_b/scanner.py`
- Pipeline expanded from 4 stages to 5 stages:
  - Stage 4 (NEW): Boltz cross-validation for top N candidates
  - Stage 5: Structure comparison (was Stage 4) now includes CV data
- New parameters: `run_boltz_crossval` (default True), `boltz_top_n` (default 200)
- Hit building enriched with Boltz confidence and agreement scores
- CV boost applied to combined similarity score (+10% max)
- Summary logging includes cross-validated hit count

### Task 4: Data Model Extension
- **File**: `decoy_a/models.py` — `StructuralScore` class
- Added 5 new Optional fields:
  - `boltz_pdb_path`: Path to Boltz-predicted structure
  - `boltz_confidence`: Boltz confidence score [0-1]
  - `boltz_iptm`: Boltz interface TM-score [0-1]
  - `boltz_rmsd`: RMSD between primary and Boltz predictions
  - `cross_validation_agreement`: Agreement score [0-1]
- All backward-compatible (Optional + default None)

### Task 5: Documentation Updates
- **README.md**: Updated Decoy B section with 5-stage pipeline, Boltz cross-validation details table, updated scoring formula
- **PROGRESS_REPORT_PHASE1.md**: Added Phase 2 section (this document)

---

## File Change Summary

| File | Action | Lines |
|------|--------|-------|
| `decoy_b/tools/boltz.py` | NEW | ~370 |
| `decoy_a/models.py` | MODIFIED | +15 (StructuralScore cross-validation fields) |
| `decoy_b/scanner.py` | MODIFIED | +90 (Boltz stage, cross-validation, CV boost) |
| `README.md` | MODIFIED | Updated Decoy B section |

---

## Pending: Server Deployment

### Step 1: Install Boltz
```bash
# Clone Boltz source (already at ~/tools/boltz or copy from local)
cd ~/tools
git clone <boltz-repo> boltz
cd boltz && pip install -e .

# Weights auto-download on first run to ~/.boltz/
# Alternatively, set BOLTZ_CACHE to a custom directory
export BOLTZ_CACHE=/share/liuyutian/boltz_cache
```

### Step 2: Copy code changes
Copy the following files from local to server:
- `decoy_b/tools/boltz.py` (NEW)
- `decoy_a/models.py` (modified)
- `decoy_b/scanner.py` (modified)
- `README.md` (modified)

### Step 3: Verify Boltz installation
```bash
python -c "from decoy_b.tools.boltz import check_available; print(check_available())"
# Should print True

# Quick test
boltz predict examples/prot.yaml --out_dir /tmp/boltz_test --devices 1
```

### Step 4: Run full pipeline with cross-validation
```bash
# Run with Boltz cross-validation enabled (default)
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# Or disable Boltz if weights not available yet
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --no-boltz
```

### Step 5: Validate cross-validation
- Check for `cross_validation_agreement` field in output JSON
- `modeling_tool` should show `+interface` suffix
- High-agreement candidates (>0.7) should rank higher than low-agreement ones

---

## Design Decisions (Phase 2)

1. **+10% CV boost cap**: Conservative — enough to differentiate but not dominate. Can be tuned after validation.
2. **Linear agreement mapping**: `agreement = 1 - RMSD/4.0`. Simple and interpretable. RMSD > 4A (poor agreement) maps to 0.
3. **Boltz-2 over Boltz-1**: Boltz-2 is the latest model with near-AF3 accuracy and joint structure/affinity prediction.
4. **Top 200 for CV**: Same as AF3 refinement count. Balances compute cost vs coverage.
5. **Three-mode fallback**: Python API → CLI → Subprocess ensures Boltz works regardless of installation method.
6. **Graceful degradation**: If Boltz unavailable, pipeline proceeds without CV boost (no regression).
