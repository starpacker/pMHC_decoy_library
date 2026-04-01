"""
Fetcher Agent
=============
Retrieves paper abstracts (and full-text when available) from PubMed / PMC
using NCBI E-utilities.  Also supports keyword-based discovery searches.

Public API
----------
    fetch_abstract(pmid)       -> PaperRecord
    search_pubmed(query, max)  -> list[str]           (list of PMIDs)
    fetch_multiple(pmids)      -> list[PaperRecord]
"""

from __future__ import annotations

import logging
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Optional

import requests

from .config import (
    NCBI_API_KEY,
    NCBI_DELAY_SEC,
    NCBI_EFETCH,
    NCBI_EMAIL,
    NCBI_ESEARCH,
)

log = logging.getLogger(__name__)


# ── Data container ───────────────────────────────────────────────────────

@dataclass
class PaperRecord:
    """Lightweight container for data pulled from PubMed."""
    pmid: str
    title: str = ""
    abstract: str = ""
    authors: List[str] = field(default_factory=list)
    journal: str = ""
    year: str = ""
    doi: str = ""
    pmc_id: str = ""
    full_text: str = ""  # populated only when PMC OA full text is available
    mesh_terms: List[str] = field(default_factory=list)


# ── Internal helpers ─────────────────────────────────────────────────────

def _ncbi_params(**extra) -> dict:
    """Base params shared by all E-utility calls."""
    params = {"email": NCBI_EMAIL, "tool": "decoy_library"}
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    params.update(extra)
    return params


def _rate_limit():
    time.sleep(NCBI_DELAY_SEC)


def _text(el: Optional[ET.Element], default: str = "") -> str:
    """Safe text extraction from an XML element."""
    if el is None:
        return default
    # Collect all text recursively (handles mixed-content like <i>, <sup>)
    return "".join(el.itertext()).strip() or default


# ── Public functions ─────────────────────────────────────────────────────

def search_pubmed(query: str, retmax: int = 20) -> List[str]:
    """
    Run an ESearch against PubMed and return a list of PMIDs.

    Parameters
    ----------
    query : str
        PubMed search query string (supports MeSH, boolean operators, etc.)
    retmax : int
        Maximum number of results to return.

    Returns
    -------
    list[str]
        PMIDs matching the query.
    """
    params = _ncbi_params(db="pubmed", term=query, retmax=retmax, retmode="json")
    _rate_limit()
    resp = requests.get(NCBI_ESEARCH, params=params, timeout=30)
    resp.raise_for_status()
    data = resp.json()
    pmids = data.get("esearchresult", {}).get("idlist", [])
    log.info("PubMed search '%s' → %d hits (returning %d)", query, 
             int(data.get("esearchresult", {}).get("count", 0)), len(pmids))
    return pmids


def fetch_abstract(pmid: str) -> PaperRecord:
    """
    Fetch a single PubMed record by PMID and parse into a PaperRecord.

    Parameters
    ----------
    pmid : str
        PubMed identifier.

    Returns
    -------
    PaperRecord
    """
    params = _ncbi_params(db="pubmed", id=pmid, rettype="xml", retmode="xml")
    _rate_limit()
    resp = requests.get(NCBI_EFETCH, params=params, timeout=30)
    resp.raise_for_status()

    root = ET.fromstring(resp.content)
    article = root.find(".//PubmedArticle")
    if article is None:
        log.warning("No PubmedArticle found for PMID %s", pmid)
        return PaperRecord(pmid=pmid)

    # ── Title ────────────────────────────────────────────────────────
    title_el = article.find(".//ArticleTitle")
    title = _text(title_el)

    # ── Abstract ─────────────────────────────────────────────────────
    abstract_parts = []
    for abs_text in article.findall(".//AbstractText"):
        label = abs_text.get("Label", "")
        txt = _text(abs_text)
        if label:
            abstract_parts.append(f"{label}: {txt}")
        else:
            abstract_parts.append(txt)
    abstract = "\n".join(abstract_parts)

    # ── Authors ──────────────────────────────────────────────────────
    authors = []
    for author in article.findall(".//Author"):
        last = _text(author.find("LastName"))
        fore = _text(author.find("ForeName"))
        if last:
            authors.append(f"{fore} {last}".strip())

    # ── Journal / Year ───────────────────────────────────────────────
    journal = _text(article.find(".//Journal/Title"))
    year_el = article.find(".//PubDate/Year")
    year = _text(year_el)
    if not year:
        medline_date = _text(article.find(".//PubDate/MedlineDate"))
        if medline_date:
            year = medline_date[:4]

    # ── DOI ──────────────────────────────────────────────────────────
    doi = ""
    for aid in article.findall(".//ArticleId"):
        if aid.get("IdType") == "doi":
            doi = _text(aid)
            break

    # ── PMC ID ───────────────────────────────────────────────────────
    pmc_id = ""
    for aid in article.findall(".//ArticleId"):
        if aid.get("IdType") == "pmc":
            pmc_id = _text(aid)
            break

    # ── MeSH terms ───────────────────────────────────────────────────
    mesh_terms = []
    for mh in article.findall(".//MeshHeading/DescriptorName"):
        mesh_terms.append(_text(mh))

    record = PaperRecord(
        pmid=pmid,
        title=title,
        abstract=abstract,
        authors=authors,
        journal=journal,
        year=year,
        doi=doi,
        pmc_id=pmc_id,
        mesh_terms=mesh_terms,
    )

    # ── Try fetching PMC full text if available ──────────────────────
    if pmc_id:
        record.full_text = _try_fetch_pmc_text(pmc_id)

    log.info("Fetched PMID %s: %s (%s)", pmid, title[:60], year)
    return record


def _try_fetch_pmc_text(pmc_id: str) -> str:
    """
    Attempt to retrieve plain-text body from PMC Open Access.
    Returns empty string on failure (many articles are paywalled).
    """
    try:
        _rate_limit()
        url = f"{NCBI_EFETCH}"
        params = _ncbi_params(db="pmc", id=pmc_id, rettype="xml", retmode="xml")
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()
        root = ET.fromstring(resp.content)
        # Collect all <p> text under <body>
        paragraphs = []
        for p in root.iter("p"):
            txt = _text(p)
            if txt:
                paragraphs.append(txt)
        full_text = "\n\n".join(paragraphs)
        if full_text:
            log.info("PMC full text retrieved for %s (%d chars)", pmc_id, len(full_text))
        return full_text
    except Exception as exc:
        log.debug("PMC full-text fetch failed for %s: %s", pmc_id, exc)
        return ""


def fetch_multiple(pmids: List[str]) -> List[PaperRecord]:
    """Convenience: fetch abstracts for a list of PMIDs."""
    records = []
    for pmid in pmids:
        try:
            records.append(fetch_abstract(pmid))
        except Exception as exc:
            log.error("Failed to fetch PMID %s: %s", pmid, exc)
            records.append(PaperRecord(pmid=pmid))
    return records


# ── Discovery search queries (pre-built for cold start) ─────────────────

COLD_START_QUERIES = [
    # Direct clinical toxicity cases
    '"TCR" AND "cross-reactivity" AND ("fatal" OR "toxicity" OR "cardiac")',
    '"affinity-enhanced TCR" AND "off-target"',
    '"MAGE-A3" AND "TCR" AND ("Titin" OR "cross-reactive")',
    # Broader TCR safety
    '"adoptive cell therapy" AND "off-target" AND "peptide"',
    '"engineered T cell" AND "cross-reactivity" AND "safety"',
    # Computational prediction
    '"ARDitox" OR "EpiTox" AND "TCR" AND "safety"',
    # X-scan / peptide scanning
    '"X-scan" AND "TCR" AND "peptide"',
    '"alanine scanning" AND "TCR" AND "cross-reactivity"',
]
