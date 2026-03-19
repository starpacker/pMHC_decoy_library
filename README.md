# Decoy C Library — Cold-Start Collection System

A multi-agent, multi-source pipeline for collecting, extracting, validating, and curating a standardized database of off-target / cross-reactive peptide "decoys" that have caused (or are predicted to cause) toxicity in TCR-based immunotherapies.

## Overview

The Decoy C Library addresses a critical safety need in T-cell receptor (TCR) immunotherapy: maintaining a curated blacklist of peptides known to cause dangerous cross-reactivity. When an engineered TCR binds an unintended peptide on healthy tissue, the consequences can be fatal — as demonstrated by the MAGE-A3/Titin case that killed patients via cardiac toxicity.

## Architecture

The system implements a three-agent pipeline with multi-source fetching:

```
┌─────────────────────────────────────────────────┐
│              DATA SOURCES                        │
│  ┌─────────┐ ┌──────────┐ ┌───────────────────┐ │
│  │ PubMed/ │ │ Semantic │ │ arXiv / bioRxiv / │ │
│  │ PMC     │ │ Scholar  │ │ medRxiv / EuropePMC│ │
│  └────┬────┘ └────┬─────┘ └────────┬──────────┘ │
│       │           │                │             │
│  ┌────┴───┐  ┌────┴────────┐  ┌───┴──────────┐  │
│  │Fetcher │  │Multi-Source │  │ClinicalTrials│  │
│  │Agent   │  │Fetcher     │  │.gov Fetcher  │  │
│  └────┬───┘  └────┬────────┘  └───┬──────────┘  │
└───────┼───────────┼───────────────┼──────────────┘
        └───────────┼───────────────┘
                    ▼
          ┌──────────────┐    ┌──────────────┐
          │  Extractor   │───>│  Validator   │
          │   Agent      │    │   Agent      │
          └──────────────┘    └──────────────┘
           LLM or Rules        UniProt + IEDB
           Structured JSON      Cross-check
```

### Fetcher Agent (`fetcher.py`)
- Searches PubMed via NCBI E-utilities (ESearch + EFetch)
- Retrieves abstracts, metadata, and MeSH terms
- Fetches PMC Open Access full text when available
- Pre-built cold-start search queries for TCR cross-reactivity literature

### Multi-Source Fetcher (`multi_source_fetcher.py`)
- **arXiv**: Searches quantitative biology and ML categories via the arXiv API
- **bioRxiv / medRxiv**: Date-range retrieval with local keyword filtering via bioRxiv API
- **ClinicalTrials.gov**: Searches trial protocols via the ClinicalTrials.gov v2 API
- **Semantic Scholar**: Bulk and standard endpoint search with automatic 429 retry
- **Europe PMC**: Searches the European PMC index (includes preprints and patents)
- Returns unified `PaperRecord` objects compatible with the existing extraction pipeline

### Extractor Agent (`extractor.py`)
- **LLM mode**: Uses OpenAI-compatible API with a specialized system prompt
- **Rule-based fallback**: Pattern matching for peptide sequences, HLA alleles, and safety-relevant keywords
- Produces structured `DecoyEntry` JSON following the standardized schema
- Includes `thought_process` field with inline citations (anti-hallucination)

### Validator Agent (`validator.py`)
- **UniProt validation**: Confirms gene_symbol ↔ uniprot_id mapping
- **IEDB validation**: Checks if peptide exists in the Immune Epitope Database
- Enriches entries with additional assay data, PubMed IDs, and TCR receptor info
- Sets validation flags: `VALIDATED`, `PARTIAL`, `NEEDS_REVIEW`, `ERROR`

## Schema

Each entry follows the standardized Decoy C Library schema:

| Section | Fields | Description |
|---------|--------|-------------|
| `peptide_info` | decoy_sequence, hla_allele, source_protein, gene_symbol, uniprot_id | Core peptide identity |
| `discovery_context` | original_target_sequence, original_target_protein, tcr_name_or_id | How the decoy was found |
| `risk_profile` | evidence_level, critical_organs_affected, expression_pattern, computational_safety_score | Danger assessment |
| `experimental_evidence` | mass_spec_confirmed, assays_performed, cross_reactivity_affinity | Hard experimental data |
| `provenance` | pmid, clinical_trial_id, evidence_summary | Literature traceability |
| `source` | title, authors, journal, year, pmid, doi, url, citation | Paper/report from which the entry was extracted |

### Evidence Levels (4-tier classification)

| Level | Name | Description |
|-------|------|-------------|
| 1 | `Level_1_Clinical_Fatal` | Patient SAE, death, or FDA clinical hold |
| 2 | `Level_2_In_Vitro_Confirmed` | In vitro cytotoxicity against healthy cells |
| 3 | `Level_3_High_Throughput_Screened` | X-scan/yeast display hits, no killing data |
| 4 | `Level_4_In_Silico_High_Risk` | Computational prediction only (ARDitox, EpiTox) |

## Quick Start

### 1. Cold Start with Seed Data Only (fastest, no PubMed search)

```bash
python -m decoy_library cold-start --seed-only
```

### 2. Full Cold Start (seed + PubMed discovery + validation)

```bash
python -m decoy_library cold-start
```

### 3. View the Library

```bash
python -m decoy_library show              # Table format
python -m decoy_library show --format json # Full JSON
python -m decoy_library stats             # Statistics
```

### 4. Add Papers by PMID

```bash
python -m decoy_library fetch --pmids 23926201 23863783
```

### 5. Search PubMed for New Papers

```bash
python -m decoy_library search --query "TCR cross-reactivity cardiac toxicity"
```

### 6. Re-validate All Entries

```bash
python -m decoy_library validate
```

### 7. Scale Up to N Entries (Automated Discovery)

```bash
# Generate 10,000 unique decoy entries automatically
python scale_up.py --required_num 10000

# Fast mode (skip validation)
python scale_up.py --required_num 500 --no-validate

# Use specific strategy
python scale_up.py --required_num 1000 --strategy keyword

# Resume from checkpoint
python scale_up.py --required_num 10000 --resume
```

---

## 🔬 Scale-Up Data Discovery System

The `scale_up.py` module implements an automated, scalable data discovery pipeline that can generate **any target number N** of unique, high-quality decoy peptide entries.

### How It Works

```
┌────────────────────────────────────────────────────────────────────────────┐
│                         SCALE-UP DATA DISCOVERY                            │
├────────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌────────────┐  │
│  │  Strategy   │───>│  PubMed     │───>│   Fetcher   │───>│ Extractor  │  │
│  │  Generator  │    │   Search    │    │   Agent     │    │   Agent    │  │
│  └─────────────┘    └─────────────┘    └─────────────┘    └────────────┘  │
│        │                  │                  │                  │         │
│        │                  │                  │                  ▼         │
│        │                  │                  │          ┌────────────┐    │
│        │                  │                  │          │ Validator  │    │
│        │                  │                  │          │   Agent    │    │
│        │                  │                  │          └────────────┘    │
│        │                  │                  │                  │         │
│        │                  │                  │                  ▼         │
│        │                  │                  │          ┌────────────┐    │
│        │                  │                  └─────────>│ Deduplicator│   │
│        │                  │                             │ (by seq)   │    │
│        │                  │                             └────────────┘    │
│        │                  │                                    │          │
│        │                  │                                    ▼          │
│        │                  │                             ┌────────────┐    │
│        │                  │                             │  Quality   │    │
│        │                  │                             │  Filter    │    │
│        │                  │                             └────────────┘    │
│        │                  │                                    │          │
│        │                  │                                    ▼          │
│        │                  │    ┌────────────────────────────────────┐     │
│        │                  └───>│  Checkpoint: data/scale_up_checkpoint │  │
│        │                       │  Library:    data/decoy_library.json  │  │
│        │                       └────────────────────────────────────┘     │
│        │                                          │                       │
│        │                                          ▼                       │
│        │                                 ┌────────────────┐               │
│        │                                 │ len(lib) >= N? │               │
│        │                                 └────────────────┘               │
│        │                                    │NO      │YES                 │
│        │◄───────────────────────────────────┘        │                    │
│        │  (next query)                               ▼                    │
│        │                                        ✅ DONE                   │
│        └──────────────────────────────────────────────────────────────────┤
└────────────────────────────────────────────────────────────────────────────┘
```

### 6 Query Expansion Strategies

The system uses **6 complementary strategies** to discover diverse papers containing decoy peptides (Strategies A–E target PubMed; Strategy F searches additional sources):

#### Strategy A: Keyword-Based Queries (149 queries)

Domain-specific boolean queries targeting TCR cross-reactivity literature:

```python
# Clinical toxicity reports
'"TCR" AND "cross-reactivity" AND ("fatal" OR "toxicity" OR "cardiac")'
'"affinity-enhanced TCR" AND "off-target"'
'"MAGE-A3" AND "TCR" AND ("Titin" OR "cross-reactive")'

# Peptide-MHC binding and specificity
'"peptide-MHC" AND "TCR" AND "cross-reactivity"'
'"T-cell receptor" AND "off-target" AND "immunotherapy"'

# High-throughput screening
'"X-scan" AND "TCR" AND "peptide"'
'"alanine scanning" AND "TCR" AND "cross-reactivity"'
'"yeast display" AND "pMHC" AND "TCR"'

# Immunopeptidome / HLA ligandome
'"immunopeptidome" AND "HLA" AND "ligand"'
'"HLA ligandome" AND "mass spectrometry"'
'"eluted peptide" AND "HLA" AND "tumor"'
```

#### Strategy B: Gene-Focused Queries (126 queries)

One query per high-risk gene family known to be involved in TCR cross-reactivity:

```python
GENES = [
    # Classical decoy genes
    "TTN", "MLANA", "PMEL", "TYR", "MAGEA3", "MAGEA4",
    
    # Cancer-testis antigens
    "CTAG1B", "NY-ESO-1", "PRAME", "SSX2", "GAGE", "BAGE",
    
    # Tumor-associated antigens
    "WT1", "AFP", "ERBB2", "CEACAM5", "MUC1", "MSLN",
    
    # Oncogenes with neoantigen potential
    "KRAS", "BRAF", "TP53", "PIK3CA", "IDH1",
    
    # Viral antigens
    "HBV", "HPV16", "EBV", "CMV", "HTLV",
    ...  # 126 genes total
]

# Generated queries:
# '("TTN" AND "TCR" AND "peptide" AND "HLA")'
# '("MAGEA3" AND "TCR" AND "peptide" AND "HLA")'
# ...
```

#### Strategy C: MeSH Term-Based Queries (20 queries)

Systematic searches using controlled Medical Subject Headings:

```python
'"Receptors, Antigen, T-Cell"[MeSH] AND "Cross Reactions"[MeSH] AND "Peptides"[MeSH]'
'"Immunotherapy, Adoptive"[MeSH] AND "Drug-Related Side Effects and Adverse Reactions"[MeSH]'
'"T-Lymphocytes, Cytotoxic"[MeSH] AND "Molecular Mimicry"[MeSH]'
'"HLA-A Antigens"[MeSH] AND "Epitopes, T-Lymphocyte"[MeSH] AND "Cross Reactions"[MeSH]'
```

#### Strategy D: Year-Range Expansion (48 queries)

Temporal slicing to retrieve papers from different publication periods:

```python
# Takes base queries and adds year filters:
'("TCR" AND "cross-reactivity" AND "toxicity") AND ("1990"[PDAT] : "1994"[PDAT])'
'("TCR" AND "cross-reactivity" AND "toxicity") AND ("1995"[PDAT] : "1999"[PDAT])'
'("TCR" AND "cross-reactivity" AND "toxicity") AND ("2000"[PDAT] : "2004"[PDAT])'
# ... up to 2026
```

#### Strategy E: LLM-Guided Query Generation (Dynamic)

Uses GPT-4 to propose novel search queries based on gaps in the existing library:

```python
# LLM analyzes existing entries and generates queries like:
'"HLA-B*35:01" AND "peptide" AND "T cell" AND "cross-reactivity"'  # Underrepresented HLA
'"SPANXB1" AND "CTL" AND "epitope"'  # Novel cancer-testis antigen
'"cardiac myosin" AND "TCR" AND "autoimmune"'  # Cardiac cross-reactivity
```

#### Strategy F: Multi-Source Queries (enabled with `--multi-source`)

Searches 6 non-PubMed sources with domain-specific queries for each:

| Source | Queries | Description |
|--------|---------|-------------|
| arXiv | 12 | Quantitative biology & ML papers on TCR/pMHC |
| bioRxiv | 12 | Preprints on TCR specificity and immunotherapy safety |
| medRxiv | 8 | Clinical preprints on TCR therapy and CAR-T adverse events |
| ClinicalTrials.gov | 15 | TCR therapy trials (MAGE-A3, NY-ESO-1, adoptive cell transfer) |
| Semantic Scholar | 15 | Cross-reactivity, cardiac toxicity, decoy mimicry literature |
| Europe PMC | 15 | European index including preprints and patents |

Each source uses fail-fast logic: after 2 consecutive empty results, remaining queries for that source are skipped to avoid wasting time on unreachable APIs.

### Infinite Mode (`--infinite`)

When `--infinite` is set, the scale-up process loops indefinitely:
- After each cycle, the query checkpoint is cleared (but processed PMIDs are preserved)
- LLM generates fresh queries each cycle based on the updated library
- Multi-source queries are re-tried each cycle
- A convenience script `scale_up_infinite.sh` wraps this for background execution:

```bash
nohup bash scale_up_infinite.sh 1000 > scale_up_infinite.log 2>&1 &
```

### Deduplication Logic

Entries are deduplicated by **peptide sequence** (case-insensitive):

```python
class DecoyLibrary:
    def add_entry(self, entry: DecoyEntry, deduplicate: bool = True) -> bool:
        """Add entry; returns False if duplicate sequence already exists."""
        if deduplicate and self.find_by_sequence(entry.peptide_info.decoy_sequence):
            return False  # Reject duplicate
        self.entries.append(entry)
        return True
```

### Quality Filtering

Low-quality entries are filtered out before counting toward target N:

```python
def _quality_filter(lib: DecoyLibrary) -> DecoyLibrary:
    """
    Keep entries that have at least ONE of:
        - Non-UNKNOWN gene symbol
        - Valid HLA allele (not HLA-A*00:xx placeholder)
        - Evidence summary length > 20 characters
    """
```

### Checkpoint & Resumability

Progress is saved to `data/scale_up_checkpoint.json`:

```json
{
  "processed_pmids": ["23926201", "23863783", ...],
  "processed_queries": ["\"TCR\" AND \"cross-reactivity\"...", ...],
  "processed_multi_source": ["arxiv:TCR pMHC cross-reactivity...", ...],
  "started_at": "2026-03-17T10:30:00",
  "rounds": 45,
  "library_size": 1523,
  "last_saved": "2026-03-17T11:15:00"
}
```

Run with `--resume` (default) to continue from the last checkpoint, or `--fresh` to start over.

### CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `--required_num, -n` | Target number of unique entries | (required) |
| `--batch-size, -b` | PMIDs to process per round | 30 |
| `--max-rounds` | Safety limit on total rounds | 500 |
| `--strategy, -s` | Query strategy: `keyword`, `gene`, `mesh`, `year`, `llm`, `all` | `all` |
| `--retmax` | Max results per PubMed query | 40 |
| `--no-validate` | Skip UniProt/IEDB validation | False |
| `--no-quality-filter` | Disable quality filtering | False |
| `--resume` | Resume from checkpoint | True |
| `--fresh` | Ignore checkpoint, start fresh | False |
| `--save-every` | Save checkpoint every N rounds | 5 |
| `--multi-source` | Enable multi-source search (arXiv, bioRxiv, etc.) | False |
| `--no-pubmed` | Skip PubMed queries (use only multi-source) | False |
| `--infinite` | Loop indefinitely, regenerating queries each cycle | False |
| `-v, --verbose` | Enable debug logging | False |

### Example Output

```
============================================================
🚀 SCALE-UP START
   Target:  10000 entries
   Current: 51 entries (51 unique sequences)
   Gap:     9949
   Strategy: all
============================================================
  📋 Strategy A (keyword): 149 queries
  🧬 Strategy B (gene): 126 queries
  🏷️  Strategy C (MeSH): 20 queries
  📅 Strategy D (year): 48 queries
  🤖 Strategy E (LLM): 20 queries
  📊 Total query pool: 363 queries

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
🔄 Round 1 | Query 1/181 | Have: 51 | Need: 9949 more
   Query: "TCR" AND "cross-reactivity" AND ("fatal" OR "toxicity"...
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
   📄 Found 40 PMIDs, 38 new, processing 30
   ✅ Added 12 new entries (extracted 15, deduped 3)
   📊 Library: 63 / 10000 (0.6%)
   💾 Checkpoint saved (round 5)
...

============================================================
🏁 SCALE-UP COMPLETE
============================================================
   Total entries:     10247
   Unique sequences:  10247
   Target:            10000
   Achieved:          102.5%
   Rounds executed:   312
   New entries added: 10196
   PMIDs processed:   8934

Evidence level distribution:
   Level_1_Clinical_Fatal: 89
   Level_2_In_Vitro_Confirmed: 1456
   Level_3_High_Throughput_Screened: 5823
   Level_4_In_Silico_High_Risk: 2879
```

---

## Seed Data (51 curated entries)

The library ships with 51 hand-curated high-confidence entries covering the most well-documented TCR cross-reactivity cases:

| ID | Peptide | Gene | HLA | Level | Case |
|----|---------|------|-----|-------|------|
| DC-0001 | ESDPIVAQY | TTN | HLA-A*01:01 | Level 1 | MAGE-A3 TCR → Titin cardiac death |
| DC-0002 | EVDPIGHVY | EPS8L2 | HLA-A*01:01 | Level 1 | MAGE-A3 TCR → Brain toxicity |
| DC-0003 | ELAGIGILTV | MLANA | HLA-A*02:01 | Level 2 | DMF5 TCR → Melanocyte destruction |
| DC-0004 | SLLMWITQV | CTAG1B | HLA-A*02:01 | Level 3 | 1G4 c259 → NY-ESO-1 X-scan screened |
| DC-0005 | IMIGVLVGV | CEACAM5 | HLA-A*02:01 | Level 1 | Anti-CEA TCR → Severe colitis |
| ... | ... | ... | ... | ... | (46 more entries) |

See `data/seed_entries.json` for the complete list.

## Configuration

Environment variables:

| Variable | Description | Default |
|----------|-------------|---------|
| `OPENAI_API_KEY` | OpenAI-compatible API key for LLM extraction | Pre-configured (DP Tech internal gateway) |
| `OPENAI_MODEL` | LLM model name | `gpt-4o` |
| `OPENAI_BASE_URL` | Custom API endpoint | `https://ai-gateway-internal.dp.tech/v1` |
| `NCBI_EMAIL` | Email for NCBI E-utilities | `decoy_library@example.com` |
| `NCBI_API_KEY` | NCBI API key (faster rate limits) | (none) |
| `SEMANTIC_SCHOLAR_API_KEY` | Semantic Scholar API key (higher rate limits) | (none — defined in `multi_source_fetcher.py`) |

## Dependencies

- Python 3.9+
- pydantic >= 2.0
- requests
- httpx (for multi-source fetcher async HTTP)
- openai (optional, for LLM extraction)

## File Structure

```
decoy_library/
├── __init__.py              # Package metadata (__version__ = "0.1.0")
├── __main__.py              # python -m entry point
├── config.py                # Configuration constants & API endpoints
├── models.py                # Pydantic v2 schema models
├── fetcher.py               # Fetcher Agent (PubMed/PMC via NCBI E-utilities)
├── multi_source_fetcher.py  # ⭐ Multi-source Fetcher (arXiv, bioRxiv, Semantic Scholar, etc.)
├── extractor.py             # Extractor Agent (LLM + rule-based fallback)
├── validator.py             # Validator Agent (UniProt + IEDB cross-check)
├── orchestrator.py          # Pipeline coordination & library I/O
├── seed_data.py             # Curated cold-start entries (writes seed_entries.json)
├── main.py                  # CLI entry point (argparse subcommands)
├── scale_up.py              # ⭐ Scale-up data discovery system
├── scale_up_infinite.sh     # Convenience wrapper for infinite mode
├── generate_50.py           # Utility: generate initial 50-entry set
├── Decoy_AB.md              # Additional documentation (Decoy A/B design notes)
├── README.md                # This file
├── scale_up_*.log           # Log files from scale-up runs (auto-generated)
└── data/
    ├── decoy_library.json       # Persistent database
    ├── seed_entries.json        # 51 curated seed entries
    └── scale_up_checkpoint.json # Resume checkpoint (auto-generated)
```

## API Usage (Programmatic)

```python
from decoy_library.orchestrator import load_library, run_cold_start
from decoy_library.fetcher import fetch_abstract
from decoy_library.validator import validate_entry

# Initialize
lib = run_cold_start()

# Or load existing
lib = load_library()

# Find a peptide
entry = lib.find_by_sequence("ESDPIVAQY")
print(entry.risk_profile.evidence_level)  # Level_1_Clinical_Fatal

# Fetch and validate a new paper
from decoy_library.extractor import extract_from_paper
paper = fetch_abstract("23926201")
entries = extract_from_paper(paper)
for e in entries:
    e = validate_entry(e)
```

## References

- Cameron et al. (2013) Sci Transl Med. PMID: 23926201
- Linette et al. (2013) Blood. PMID: 23863783
- Morgan et al. (2013) J Immunother. PMID: 24475783
- Raman et al. (2016) J Biol Chem. PMID: 26457759
- Parkhurst et al. (2011) Mol Ther. PMID: 21282551
