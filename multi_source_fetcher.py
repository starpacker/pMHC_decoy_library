#!/usr/bin/env python3
"""
Multi-Source Fetcher
====================
Fetches papers from multiple sources beyond PubMed:
- arXiv (via arXiv API)
- bioRxiv / medRxiv (via bioRxiv API)
- ClinicalTrials.gov (via ClinicalTrials API)
- Semantic Scholar (via Semantic Scholar API)
- Europe PMC (via Europe PMC API)

Each source returns standardized PaperRecord objects.
"""

from __future__ import annotations

import logging
import os
import re
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any
from urllib.parse import quote_plus

import requests

log = logging.getLogger(__name__)

# Rate limiting delays (seconds)
ARXIV_DELAY = 3.0  # arXiv asks for 3 seconds between requests
BIORXIV_DELAY = 1.0
CLINICAL_TRIALS_DELAY = 0.5
SEMANTIC_SCHOLAR_DELAY = 5.0  # increased to avoid 429
EUROPEPMC_DELAY = 0.5

# Optional API keys (set via environment variables)
SEMANTIC_SCHOLAR_API_KEY = os.environ.get("SEMANTIC_SCHOLAR_API_KEY", "")


@dataclass
class PaperRecord:
    """Unified paper record from any source."""
    source: str  # 'pubmed', 'arxiv', 'biorxiv', 'clinicaltrials', 'semantic_scholar', 'europepmc'
    source_id: str  # Source-specific ID (PMID, arXiv ID, DOI, NCT ID, etc.)
    title: str = ""
    abstract: str = ""
    authors: List[str] = field(default_factory=list)
    year: str = ""
    doi: str = ""
    full_text: str = ""
    url: str = ""
    keywords: List[str] = field(default_factory=list)
    
    # For compatibility with existing code
    @property
    def pmid(self) -> str:
        """Return source_id for backward compatibility."""
        return self.source_id


# ══════════════════════════════════════════════════════════════════════════
# arXiv Fetcher
# ══════════════════════════════════════════════════════════════════════════

ARXIV_API = "https://export.arxiv.org/api/query"

def search_arxiv(query: str, max_results: int = 50) -> List[PaperRecord]:
    """
    Search arXiv for papers.
    
    Parameters
    ----------
    query : str
        Search query (supports arXiv query syntax)
    max_results : int
        Maximum number of results
        
    Returns
    -------
    List[PaperRecord]
    """
    time.sleep(ARXIV_DELAY)
    
    # Build arXiv query: restrict to bio/ML categories for relevance
    # Use AND between key terms in abstract; restrict to biology categories
    # arXiv categories: q-bio.* (quantitative biology), cs.LG, cs.AI
    terms = query.strip().split()
    # Join terms with AND for abstract search
    abs_query = " AND ".join(f"abs:{t}" for t in terms if len(t) > 2)
    arxiv_query = f'({abs_query}) AND (cat:q-bio* OR cat:cs.LG OR cat:cs.AI)'
    
    params = {
        "search_query": arxiv_query,
        "start": 0,
        "max_results": max_results,
        "sortBy": "relevance",
        "sortOrder": "descending"
    }
    
    try:
        resp = requests.get(ARXIV_API, params=params, timeout=30)
        resp.raise_for_status()
    except Exception as exc:
        log.warning("arXiv search failed: %s", exc)
        return []
    
    # Parse Atom feed
    papers = []
    try:
        root = ET.fromstring(resp.content)
        ns = {"atom": "http://www.w3.org/2005/Atom"}
        
        for entry in root.findall("atom:entry", ns):
            arxiv_id = entry.find("atom:id", ns)
            arxiv_id = arxiv_id.text.split("/abs/")[-1] if arxiv_id is not None else ""
            
            title = entry.find("atom:title", ns)
            title = title.text.strip().replace("\n", " ") if title is not None else ""
            
            abstract = entry.find("atom:summary", ns)
            abstract = abstract.text.strip().replace("\n", " ") if abstract is not None else ""
            
            authors = []
            for author in entry.findall("atom:author", ns):
                name = author.find("atom:name", ns)
                if name is not None:
                    authors.append(name.text)
            
            published = entry.find("atom:published", ns)
            year = published.text[:4] if published is not None else ""
            
            # Get DOI if available
            doi = ""
            for link in entry.findall("atom:link", ns):
                if link.get("title") == "doi":
                    doi = link.get("href", "").replace("http://dx.doi.org/", "")
            
            paper = PaperRecord(
                source="arxiv",
                source_id=arxiv_id,
                title=title,
                abstract=abstract,
                authors=authors,
                year=year,
                doi=doi,
                url=f"https://arxiv.org/abs/{arxiv_id}"
            )
            papers.append(paper)
            
    except Exception as exc:
        log.warning("Failed to parse arXiv response: %s", exc)
    
    log.info("arXiv search '%s' → %d papers", query[:50], len(papers))
    return papers


# ══════════════════════════════════════════════════════════════════════════
# bioRxiv / medRxiv Fetcher
# ══════════════════════════════════════════════════════════════════════════

BIORXIV_API = "https://api.biorxiv.org/details/biorxiv"
MEDRXIV_API = "https://api.biorxiv.org/details/medrxiv"

def search_biorxiv(query: str, max_results: int = 50, server: str = "biorxiv") -> List[PaperRecord]:
    """
    Search bioRxiv or medRxiv for papers.
    
    Note: bioRxiv API doesn't support full-text search directly.
    We use date-based retrieval and filter locally.
    
    Parameters
    ----------
    query : str
        Search terms to filter results
    max_results : int
        Maximum number of results
    server : str
        'biorxiv' or 'medrxiv'
        
    Returns
    -------
    List[PaperRecord]
    """
    time.sleep(BIORXIV_DELAY)
    
    base_url = BIORXIV_API if server == "biorxiv" else MEDRXIV_API
    
    # bioRxiv API requires date range, get recent papers
    # Format: /details/{server}/{interval}/{cursor}
    # interval: YYYY-MM-DD/YYYY-MM-DD
    from datetime import datetime, timedelta
    end_date = datetime.now().strftime("%Y-%m-%d")
    start_date = (datetime.now() - timedelta(days=365)).strftime("%Y-%m-%d")  # Last year
    
    url = f"{base_url}/{start_date}/{end_date}/0"
    
    try:
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        log.warning("%s search failed: %s", server, exc)
        return []
    
    papers = []
    query_lower = query.lower()
    # Filter out common stop words to avoid overly broad matching
    stop_words = {"and", "or", "the", "a", "an", "in", "of", "for", "to", "with", "by", "on", "at", "is", "are"}
    query_terms = [t for t in query_lower.split() if t not in stop_words and len(t) > 2]
    # Require at least 2 query terms to match (prevents overly broad results)
    min_match = min(2, len(query_terms))
    
    for item in data.get("collection", []):
        # Filter by query terms in title or abstract
        title = item.get("title", "")
        abstract = item.get("abstract", "")
        combined = (title + " " + abstract).lower()
        
        # Count how many key terms match
        matches = sum(1 for term in query_terms if term in combined)
        if matches < min_match:
            continue
        
        paper = PaperRecord(
            source=server,
            source_id=item.get("doi", ""),
            title=title,
            abstract=abstract,
            authors=item.get("authors", "").split("; "),
            year=item.get("date", "")[:4],
            doi=item.get("doi", ""),
            url=f"https://www.{server}.org/content/{item.get('doi', '')}"
        )
        papers.append(paper)
        
        if len(papers) >= max_results:
            break
    
    log.info("%s search '%s' → %d papers", server, query[:50], len(papers))
    return papers


# ══════════════════════════════════════════════════════════════════════════
# ClinicalTrials.gov Fetcher
# ══════════════════════════════════════════════════════════════════════════

CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"

def search_clinical_trials(query: str, max_results: int = 50) -> List[PaperRecord]:
    """
    Search ClinicalTrials.gov for relevant trials.
    
    Parameters
    ----------
    query : str
        Search query
    max_results : int
        Maximum number of results
        
    Returns
    -------
    List[PaperRecord]
    """
    time.sleep(CLINICAL_TRIALS_DELAY)
    
    params = {
        "query.term": query,
        "pageSize": min(max_results, 100),
        "format": "json"
    }
    
    try:
        resp = requests.get(CLINICALTRIALS_API, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        log.warning("ClinicalTrials.gov search failed: %s", exc)
        return []
    
    papers = []
    
    for study in data.get("studies", []):
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        desc_module = protocol.get("descriptionModule", {})
        
        nct_id = id_module.get("nctId", "")
        title = id_module.get("officialTitle", "") or id_module.get("briefTitle", "")
        
        # Combine brief and detailed descriptions
        brief_summary = desc_module.get("briefSummary", "")
        detailed_desc = desc_module.get("detailedDescription", "")
        abstract = f"{brief_summary}\n\n{detailed_desc}".strip()
        
        # Get conditions and interventions as keywords
        conditions = protocol.get("conditionsModule", {}).get("conditions", [])
        interventions = []
        for i in protocol.get("armsInterventionsModule", {}).get("interventions", []):
            interventions.append(i.get("name", ""))
        
        paper = PaperRecord(
            source="clinicaltrials",
            source_id=nct_id,
            title=title,
            abstract=abstract,
            year=protocol.get("statusModule", {}).get("startDateStruct", {}).get("date", "")[:4],
            url=f"https://clinicaltrials.gov/study/{nct_id}",
            keywords=conditions + interventions
        )
        papers.append(paper)
    
    log.info("ClinicalTrials search '%s' → %d trials", query[:50], len(papers))
    return papers


# ══════════════════════════════════════════════════════════════════════════
# Semantic Scholar Fetcher
# ══════════════════════════════════════════════════════════════════════════

SEMANTIC_SCHOLAR_API = "https://api.semanticscholar.org/graph/v1/paper/search"
SEMANTIC_SCHOLAR_BULK_API = "https://api.semanticscholar.org/graph/v1/paper/search/bulk"

def search_semantic_scholar(query: str, max_results: int = 50) -> List[PaperRecord]:
    """
    Search Semantic Scholar for papers.
    
    Uses the bulk endpoint (more lenient rate limits) with fallback
    to the standard endpoint.
    
    Parameters
    ----------
    query : str
        Search query
    max_results : int
        Maximum number of results
        
    Returns
    -------
    List[PaperRecord]
    """
    time.sleep(SEMANTIC_SCHOLAR_DELAY)
    
    headers = {}
    if SEMANTIC_SCHOLAR_API_KEY:
        headers["x-api-key"] = SEMANTIC_SCHOLAR_API_KEY
    
    data = None
    
    # Try bulk endpoint first (more lenient rate limits)
    for endpoint in [SEMANTIC_SCHOLAR_BULK_API, SEMANTIC_SCHOLAR_API]:
        params = {
            "query": query,
            "limit" if endpoint == SEMANTIC_SCHOLAR_API else "limit": min(max_results, 100),
            "fields": "paperId,title,abstract,authors,year,externalIds,url"
        }
        
        for attempt in range(2):
            try:
                resp = requests.get(endpoint, params=params,
                                    headers=headers, timeout=30)
                if resp.status_code == 429:
                    wait = (attempt + 1) * 10
                    log.warning("Semantic Scholar 429 at %s, waiting %ds (attempt %d/2)",
                                "bulk" if "bulk" in endpoint else "standard",
                                wait, attempt + 1)
                    time.sleep(wait)
                    continue
                resp.raise_for_status()
                data = resp.json()
                break
            except requests.exceptions.HTTPError:
                continue
            except Exception as exc:
                log.warning("Semantic Scholar search failed at %s: %s",
                            "bulk" if "bulk" in endpoint else "standard", exc)
                break
        
        if data and data.get("data"):
            break
    
    if not data or not data.get("data"):
        log.warning("Semantic Scholar: no results for '%s'", query[:50])
        return []
    
    papers = []
    
    for item in data.get("data", []):
        external_ids = item.get("externalIds") or {}
        
        paper = PaperRecord(
            source="semantic_scholar",
            source_id=item.get("paperId", ""),
            title=item.get("title", "") or "",
            abstract=item.get("abstract", "") or "",
            authors=[a.get("name", "") for a in (item.get("authors") or [])],
            year=str(item.get("year", "") or ""),
            doi=external_ids.get("DOI", "") or "",
            url=item.get("url", "") or ""
        )
        papers.append(paper)
    
    log.info("Semantic Scholar search '%s' → %d papers", query[:50], len(papers))
    return papers


# ══════════════════════════════════════════════════════════════════════════
# Europe PMC Fetcher
# ══════════════════════════════════════════════════════════════════════════

EUROPEPMC_API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

def search_europe_pmc(query: str, max_results: int = 50) -> List[PaperRecord]:
    """
    Search Europe PMC for papers (includes preprints and patents).
    
    Parameters
    ----------
    query : str
        Search query
    max_results : int
        Maximum number of results
        
    Returns
    -------
    List[PaperRecord]
    """
    time.sleep(EUROPEPMC_DELAY)
    
    params = {
        "query": query,
        "format": "json",
        "pageSize": min(max_results, 100),
        "resultType": "core"  # Include full abstract
    }
    
    try:
        resp = requests.get(EUROPEPMC_API, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        log.warning("Europe PMC search failed: %s", exc)
        return []
    
    papers = []
    
    for item in data.get("resultList", {}).get("result", []):
        source_type = item.get("source", "MED").lower()
        
        paper = PaperRecord(
            source=f"europepmc_{source_type}",
            source_id=item.get("pmid", "") or item.get("id", ""),
            title=item.get("title", ""),
            abstract=item.get("abstractText", "") or "",
            authors=[f"{a.get('firstName', '')} {a.get('lastName', '')}".strip() 
                     for a in item.get("authorList", {}).get("author", [])],
            year=item.get("pubYear", ""),
            doi=item.get("doi", ""),
            url=item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url", "")
        )
        papers.append(paper)
    
    log.info("Europe PMC search '%s' → %d papers", query[:50], len(papers))
    return papers


# ══════════════════════════════════════════════════════════════════════════
# Unified Multi-Source Search
# ══════════════════════════════════════════════════════════════════════════

def search_all_sources(
    query: str,
    max_per_source: int = 30,
    sources: List[str] = None
) -> List[PaperRecord]:
    """
    Search multiple sources for papers.
    
    Parameters
    ----------
    query : str
        Search query
    max_per_source : int
        Maximum results per source
    sources : List[str]
        Sources to search. Default: all sources.
        Options: 'arxiv', 'biorxiv', 'medrxiv', 'clinicaltrials', 
                 'semantic_scholar', 'europepmc'
        
    Returns
    -------
    List[PaperRecord]
    """
    if sources is None:
        sources = ['arxiv', 'biorxiv', 'medrxiv', 'clinicaltrials', 
                   'semantic_scholar', 'europepmc']
    
    all_papers = []
    seen_ids = set()
    
    for source in sources:
        try:
            if source == 'arxiv':
                papers = search_arxiv(query, max_per_source)
            elif source == 'biorxiv':
                papers = search_biorxiv(query, max_per_source, 'biorxiv')
            elif source == 'medrxiv':
                papers = search_biorxiv(query, max_per_source, 'medrxiv')
            elif source == 'clinicaltrials':
                papers = search_clinical_trials(query, max_per_source)
            elif source == 'semantic_scholar':
                papers = search_semantic_scholar(query, max_per_source)
            elif source == 'europepmc':
                papers = search_europe_pmc(query, max_per_source)
            else:
                log.warning("Unknown source: %s", source)
                continue
            
            # Deduplicate by DOI or source_id
            for paper in papers:
                key = paper.doi if paper.doi else f"{paper.source}:{paper.source_id}"
                if key not in seen_ids:
                    seen_ids.add(key)
                    all_papers.append(paper)
                    
        except Exception as exc:
            log.error("Error searching %s: %s", source, exc)
    
    log.info("Multi-source search '%s' → %d total papers from %d sources",
             query[:50], len(all_papers), len(sources))
    
    return all_papers


# ══════════════════════════════════════════════════════════════════════════
# TCR/Immunotherapy-Specific Queries for Each Source
# ══════════════════════════════════════════════════════════════════════════

ARXIV_QUERIES = [
    "TCR T cell receptor cross-reactivity",
    "immunotherapy off-target toxicity",
    "peptide MHC binding prediction",
    "T cell epitope prediction",
    "CAR-T cell safety",
    "neoantigen prediction",
    "HLA peptide binding",
]

BIORXIV_QUERIES = [
    "TCR cross-reactivity",
    "T cell receptor off-target",
    "immunotherapy safety",
    "peptide HLA",
    "CAR-T toxicity",
    "adoptive cell therapy",
]

CLINICAL_TRIALS_QUERIES = [
    "TCR T cell receptor therapy",
    "MAGE-A3 TCR",
    "NY-ESO-1 TCR",
    "adoptive cell therapy solid tumor",
    "affinity enhanced TCR",
    "engineered T cell cancer",
]

SEMANTIC_SCHOLAR_QUERIES = [
    "TCR cross-reactivity off-target peptide",
    "T cell receptor cardiac toxicity",
    "MAGE-A3 Titin cross-reactive",
    "HLA-restricted peptide decoy",
    "TCR molecular mimicry autoimmunity",
]
