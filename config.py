"""
Configuration constants and helpers for the Decoy C Library pipeline.
"""

import os
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data"
LIBRARY_JSON = DATA_DIR / "decoy_library.json"

# ── NCBI E-utilities ─────────────────────────────────────────────────────
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_ESEARCH = f"{NCBI_BASE}/esearch.fcgi"
NCBI_EFETCH = f"{NCBI_BASE}/efetch.fcgi"
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "decoy_library@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")

# ── UniProt REST API ─────────────────────────────────────────────────────
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"

# ── IEDB Query API ───────────────────────────────────────────────────────
IEDB_EPITOPE_SEARCH = "https://query-api.iedb.org/epitope_search"

# ── LLM (OpenAI-compatible) ─────────────────────────────────────────────
# Default: Use DP Tech internal API gateway with gpt-4o (most stable)
# Override via environment variables if needed
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "sk-Zj3a7RQDVCXr-Axg-0gtkg")
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o")
OPENAI_BASE_URL = os.getenv("OPENAI_BASE_URL", "https://ai-gateway-internal.dp.tech/v1")

# ── Alternative Models (tested and stable on DP Tech gateway) ────────────
# Ranked by speed and reliability:
#   1. gpt-4o           - 1.48s avg, 100% success (RECOMMENDED)
#   2. gemini-2.5-flash - 2.60s avg, 100% success (fast alternative)
#   3. claude-sonnet-4-5 - 3.19s avg, 100% success (Claude option)
#   4. gpt-5.2          - 3.99s avg, 100% success (slower but capable)

# ── Rate-limiting ────────────────────────────────────────────────────────
NCBI_DELAY_SEC = 0.34  # NCBI asks ≤3 req/sec without an API key
UNIPROT_DELAY_SEC = 0.5
IEDB_DELAY_SEC = 1.0
