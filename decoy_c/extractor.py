"""
Extractor Agent
===============
Uses an LLM (OpenAI-compatible API) to convert raw paper text into
structured DecoyEntry records following the Decoy C Library schema.

When no LLM API key is available, falls back to a rule-based heuristic
extractor that can still parse well-known cases.

Public API
----------
    extract_from_paper(paper: PaperRecord) -> list[DecoyEntry]
"""

from __future__ import annotations

import json
import logging
import re
from typing import Any, Dict, List, Optional

from .config import OPENAI_API_KEY, OPENAI_BASE_URL, OPENAI_MODEL
from .fetcher import PaperRecord
from .models import (
    ComputationalSafetyScore,
    DecoyEntry,
    DiscoveryContext,
    EvidenceLevel,
    ExperimentalEvidence,
    ExpressionPattern,
    PeptideInfo,
    Provenance,
    RiskProfile,
)

log = logging.getLogger(__name__)

# ── System prompt for the LLM ───────────────────────────────────────────

SYSTEM_PROMPT = (
    "You are a highly specialised biomedical data-extraction agent.\n"
    "Your task is to extract structured 'Decoy' entries from immunotherapy papers.\n\n"
    "A 'Decoy' is an off-target peptide that a T-cell receptor (TCR) cross-reacts\n"
    "with, causing (or risking) unintended toxicity against healthy tissue.\n\n"
    "For every distinct decoy peptide you find in the provided text, output a JSON\n"
    "object following this EXACT schema:\n\n"
    "{\n"
    '  "decoy_id": "PLACEHOLDER",\n'
    '  "peptide_info": {\n'
    '    "decoy_sequence": "<8-15 uppercase AA letters>",\n'
    '    "hla_allele": "<HLA-X*NN:NN format>",\n'
    '    "source_protein": "<protein name>",\n'
    '    "gene_symbol": "<HGNC symbol>",\n'
    '    "uniprot_id": "<UniProt accession or null>"\n'
    "  },\n"
    '  "discovery_context": {\n'
    '    "original_target_sequence": "<intended target peptide or null>",\n'
    '    "original_target_protein": "<intended target protein or null>",\n'
    '    "tcr_name_or_id": "<TCR clone name or null>"\n'
    "  },\n"
    '  "risk_profile": {\n'
    '    "evidence_level": "<Level_1_Clinical_Fatal | Level_2_In_Vitro_Confirmed | Level_3_High_Throughput_Screened | Level_4_In_Silico_High_Risk>",\n'
    '    "critical_organs_affected": ["<organ1>", ...],\n'
    '    "expression_pattern": "<Ubiquitous | Essential_Organ_Specific | Restricted | null>",\n'
    '    "computational_safety_score": {"ARDitox_score": null, "other_scores": {}}\n'
    "  },\n"
    '  "experimental_evidence": {\n'
    '    "mass_spec_confirmed": <true|false|null>,\n'
    '    "assays_performed": ["<assay1>", ...],\n'
    '    "cross_reactivity_affinity": "<EC50/KD value or qualitative or null>"\n'
    "  },\n"
    '  "provenance": {\n'
    '    "pmid": ["<PMID>"],\n'
    '    "clinical_trial_id": ["<NCTxxxxxxxx>"],\n'
    '    "evidence_summary": "<1-2 sentence summary>"\n'
    "  },\n"
    '  "thought_process": "<Your reasoning with inline quotes from the text>"\n'
    "}\n\n"
    "STRICT RULES:\n"
    "1. decoy_sequence MUST be pure uppercase amino-acid letters (8-15 chars).\n"
    "2. hla_allele MUST be in HLA-X*NN:NN format. If only low-res (e.g. HLA-A2), use HLA-A*02:xx.\n"
    "3. evidence_level MUST be exactly one of the four enum values.\n"
    "4. expression_pattern MUST be one of: Ubiquitous, Essential_Organ_Specific, Restricted, or null.\n"
    "5. thought_process MUST quote original sentences from the paper to justify each field.\n"
    "6. If information is not available, use null (not empty string).\n"
    "7. Output ONLY a JSON array of decoy objects. No markdown, no commentary.\n"
    "8. If no decoy peptides are found, output an empty array: []\n"
)


# ── LLM-based extraction ────────────────────────────────────────────────

def _call_llm(paper_text: str) -> str:
    """Call OpenAI-compatible chat completion and return raw text."""
    try:
        from openai import OpenAI
    except ImportError:
        raise RuntimeError("openai package is required for LLM extraction")

    kwargs: Dict[str, Any] = {"api_key": OPENAI_API_KEY}
    if OPENAI_BASE_URL:
        kwargs["base_url"] = OPENAI_BASE_URL

    client = OpenAI(**kwargs)

    # Truncate to ~12k tokens worth of text to stay within context
    max_chars = 48_000
    if len(paper_text) > max_chars:
        paper_text = paper_text[:max_chars] + "\n\n[TEXT TRUNCATED]"

    response = client.chat.completions.create(
        model=OPENAI_MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"Extract all decoy peptide entries from this paper:\n\n{paper_text}"},
        ],
        temperature=0.1,
        max_tokens=4096,
    )
    return response.choices[0].message.content or "[]"


def _parse_llm_response(raw: str, pmid: str) -> List[Dict]:
    """Parse the LLM JSON response, handling markdown fencing etc."""
    # Strip markdown code fencing if present
    raw = raw.strip()
    if raw.startswith("```"):
        raw = re.sub(r"^```(?:json)?\s*\n?", "", raw)
        raw = re.sub(r"\n?```\s*$", "", raw)
    
    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        log.error("LLM returned invalid JSON for PMID %s: %s", pmid, exc)
        return []

    if isinstance(data, dict):
        data = [data]
    if not isinstance(data, list):
        log.error("LLM returned unexpected type for PMID %s: %s", pmid, type(data))
        return []
    return data


def _dict_to_entry(d: Dict, pmid: str) -> Optional[DecoyEntry]:
    """Convert a raw dict from LLM output to a validated DecoyEntry."""
    try:
        # Ensure pmid is in provenance
        prov = d.get("provenance", {})
        if pmid and pmid not in prov.get("pmid", []):
            prov.setdefault("pmid", []).append(pmid)
        d["provenance"] = prov

        # ── Fix decoy_id: LLM may return various formats ────────────────
        # Accept: "DC-0000", "Decoy_1", "PLACEHOLDER", None, "", etc.
        # Always normalize to "DC-0000" (will be reassigned by orchestrator)
        decoy_id = d.get("decoy_id", "")
        if not decoy_id or not re.match(r"^DC-\d{4,}$", str(decoy_id)):
            d["decoy_id"] = "DC-0000"

        # ── Fix HLA allele format ───────────────────────────────────────
        # LLM may return "HLA-A*01:xx" or "HLA-A1" - normalize to valid format
        pep_info = d.get("peptide_info", {})
        hla = pep_info.get("hla_allele", "")
        if hla:
            # Replace :xx with :01 (placeholder resolution)
            hla = re.sub(r":xx$", ":01", hla)
            # If no colon, try to add proper format
            if "*" in hla and ":" not in hla:
                # e.g., HLA-A*01 -> HLA-A*01:01
                hla = hla + ":01"
            pep_info["hla_allele"] = hla
            d["peptide_info"] = pep_info

        # ── Fix critical_organs_affected: ensure list of strings ────────
        risk = d.get("risk_profile", {})
        organs = risk.get("critical_organs_affected", [])
        if organs:
            # Capitalize organ names for consistency
            risk["critical_organs_affected"] = [
                o.capitalize() if isinstance(o, str) else str(o) 
                for o in organs
            ]
            d["risk_profile"] = risk

        entry = DecoyEntry.model_validate(d)
        return entry
    except Exception as exc:
        log.warning("Failed to validate LLM extraction for PMID %s: %s\nRaw: %s",
                    pmid, exc, json.dumps(d, indent=2)[:500])
        return None


def extract_with_llm(paper: PaperRecord) -> List[DecoyEntry]:
    """
    Extract decoy entries from a paper using an LLM.

    Parameters
    ----------
    paper : PaperRecord
        Paper data from the Fetcher agent.

    Returns
    -------
    list[DecoyEntry]
        Validated decoy entries (may be empty).
    """
    # Build input text: prefer full text, fall back to abstract
    text_parts = [f"Title: {paper.title}", f"PMID: {paper.pmid}"]
    if paper.full_text:
        text_parts.append(f"Full Text:\n{paper.full_text}")
    elif paper.abstract:
        text_parts.append(f"Abstract:\n{paper.abstract}")
    else:
        log.warning("No text content for PMID %s – skipping LLM extraction", paper.pmid)
        return []

    paper_text = "\n\n".join(text_parts)
    raw = _call_llm(paper_text)
    dicts = _parse_llm_response(raw, paper.pmid)

    entries = []
    for d in dicts:
        entry = _dict_to_entry(d, paper.pmid)
        if entry:
            entries.append(entry)

    log.info("LLM extracted %d decoy entries from PMID %s", len(entries), paper.pmid)
    return entries


# ── Rule-based fallback extractor ────────────────────────────────────────

# Patterns for common amino acid sequences in text (8-11 uppercase letters)
_AA_PATTERN = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY]{8,15})\b")
# HLA allele pattern
_HLA_PATTERN = re.compile(r"HLA-[A-Z]+\d?\*\d{2}:\d{2}")
# NCT trial ID
_NCT_PATTERN = re.compile(r"NCT\d{7,8}")


def extract_rule_based(paper: PaperRecord) -> List[DecoyEntry]:
    """
    Simple rule-based extraction for when no LLM is available.
    
    This will only catch obvious cases where peptide sequences, HLA alleles,
    and key terms (off-target, cross-reactive, toxicity) co-occur.
    
    Returns
    -------
    list[DecoyEntry]
        Best-effort extraction (likely incomplete).
    """
    text = f"{paper.title}\n{paper.abstract}\n{paper.full_text}"
    text_lower = text.lower()

    # Only process papers that look relevant
    relevance_terms = ["off-target", "cross-react", "decoy", "toxicity", "adverse event",
                       "cross-recognition", "safety", "cardiac", "fatal"]
    if not any(term in text_lower for term in relevance_terms):
        return []

    # Find peptide sequences
    aa_matches = _AA_PATTERN.findall(text)
    hla_matches = _HLA_PATTERN.findall(text)
    nct_matches = _NCT_PATTERN.findall(text)

    if not aa_matches:
        return []

    # Determine evidence level from keywords
    if any(w in text_lower for w in ["fatal", "death", "died", "cardiogenic shock",
                                       "serious adverse event", "sae", "clinical hold"]):
        evidence_level = EvidenceLevel.LEVEL_1
    elif any(w in text_lower for w in ["cytotoxicity", "cell killing", "chromium release",
                                         "in vitro"]):
        evidence_level = EvidenceLevel.LEVEL_2
    elif any(w in text_lower for w in ["x-scan", "yeast display", "high-throughput"]):
        evidence_level = EvidenceLevel.LEVEL_3
    else:
        evidence_level = EvidenceLevel.LEVEL_4

    # Determine affected organs
    organs = []
    organ_map = {
        "heart": "Heart", "cardiac": "Heart", "cardiomyocyte": "Heart",
        "brain": "Brain", "neural": "Brain", "lung": "Lungs",
        "liver": "Liver", "kidney": "Kidney", "muscle": "Skeletal Muscle",
    }
    for keyword, organ in organ_map.items():
        if keyword in text_lower and organ not in organs:
            organs.append(organ)

    # Build entries — we can't reliably pair sequences with their context,
    # so we emit one entry per unique peptide found
    entries = []
    seen_seqs = set()
    for seq in aa_matches:
        if seq in seen_seqs:
            continue
        seen_seqs.add(seq)

        # Skip very common words that happen to be all-caps AA letters
        if seq in ("METHODS", "RESULTS", "ABSTRACT", "RESEARCH", "CLINICAL",
                    "PATIENTS", "THERAPY", "PEPTIDE", "PROTEIN", "DESIGN"):
            continue

        entry = DecoyEntry(
            decoy_id="DC-0000",
            peptide_info=PeptideInfo(
                decoy_sequence=seq,
                hla_allele=hla_matches[0] if hla_matches else "HLA-A*00:xx",
                source_protein="Unknown (rule-based extraction)",
                gene_symbol="UNKNOWN",
                uniprot_id=None,
            ),
            discovery_context=DiscoveryContext(
                tcr_name_or_id=None,
            ),
            risk_profile=RiskProfile(
                evidence_level=evidence_level,
                critical_organs_affected=organs,
            ),
            experimental_evidence=ExperimentalEvidence(
                mass_spec_confirmed=None,
                assays_performed=["rule-based extraction — needs manual review"],
            ),
            provenance=Provenance(
                pmid=[paper.pmid] if paper.pmid else [],
                clinical_trial_id=nct_matches,
                evidence_summary=f"Auto-extracted from PMID {paper.pmid} by rule-based parser. Needs human review.",
            ),
            thought_process="Rule-based extraction: peptide sequence found in relevant paper. Manual curation required.",
            validation_flags={"extraction_method": "rule_based", "needs_review": "true"},
        )
        entries.append(entry)

    log.info("Rule-based extractor found %d candidate sequences in PMID %s", len(entries), paper.pmid)
    return entries


# ── Post-extraction source verification ────────────────────────────────

def verify_peptide_in_source(entries: List[DecoyEntry], paper_text: str) -> List[DecoyEntry]:
    """
    Check whether each extracted peptide sequence actually appears in the
    source paper text.  Sets ``validation_flags["source_text_match"]`` to
    ``FOUND`` or ``NOT_FOUND``.

    This is the primary defence against LLM hallucination of sequences.
    """
    text_upper = paper_text.upper()
    for entry in entries:
        seq = entry.peptide_info.decoy_sequence.upper()
        if seq in text_upper:
            entry.validation_flags["source_text_match"] = "FOUND"
        else:
            # Also try with common delimiters stripped (spaces, dashes)
            text_clean = re.sub(r"[\s\-–—·.]", "", text_upper)
            if seq in text_clean:
                entry.validation_flags["source_text_match"] = "FOUND"
            else:
                entry.validation_flags["source_text_match"] = "NOT_FOUND"
                log.warning(
                    "Peptide %s NOT found in source text of PMID %s — possible hallucination",
                    seq, entry.provenance.pmid,
                )
    return entries


# ── Dual-extraction consensus ──────────────────────────────────────────

def _extract_consensus(paper: PaperRecord) -> List[DecoyEntry]:
    """
    Run LLM extraction twice at different temperatures and keep only
    peptide sequences that appear in BOTH passes (consensus filtering).

    This drastically reduces hallucinated sequences at the cost of 2× API
    calls.  Sequences unique to one pass are logged for review.
    """
    text_parts = [f"Title: {paper.title}", f"PMID: {paper.pmid}"]
    if paper.full_text:
        text_parts.append(f"Full Text:\n{paper.full_text}")
    elif paper.abstract:
        text_parts.append(f"Abstract:\n{paper.abstract}")
    else:
        return []

    paper_text = "\n\n".join(text_parts)

    # --- Pass 1: temperature=0.1 (deterministic) ---
    raw1 = _call_llm(paper_text)
    dicts1 = _parse_llm_response(raw1, paper.pmid)

    # --- Pass 2: temperature=0.4 (slightly creative) ---
    raw2 = _call_llm_t2(paper_text)
    dicts2 = _parse_llm_response(raw2, paper.pmid)

    # Build sequence sets
    seqs1 = {d.get("peptide_info", {}).get("decoy_sequence", "").upper()
             for d in dicts1}
    seqs2 = {d.get("peptide_info", {}).get("decoy_sequence", "").upper()
             for d in dicts2}

    consensus_seqs = seqs1 & seqs2
    only_pass1 = seqs1 - seqs2
    only_pass2 = seqs2 - seqs1

    if only_pass1:
        log.info("Consensus: pass-1-only sequences (dropped): %s", only_pass1)
    if only_pass2:
        log.info("Consensus: pass-2-only sequences (dropped): %s", only_pass2)

    # Keep dicts from pass 1 that are in consensus (pass 1 is lower temp → higher fidelity)
    kept = [d for d in dicts1
            if d.get("peptide_info", {}).get("decoy_sequence", "").upper() in consensus_seqs]
    # Also add any pass-2-only that appear in consensus (shouldn't happen, but safe)
    kept_seqs = {d.get("peptide_info", {}).get("decoy_sequence", "").upper() for d in kept}
    for d in dicts2:
        seq = d.get("peptide_info", {}).get("decoy_sequence", "").upper()
        if seq in consensus_seqs and seq not in kept_seqs:
            kept.append(d)
            kept_seqs.add(seq)

    entries = []
    for d in kept:
        entry = _dict_to_entry(d, paper.pmid)
        if entry:
            entry.validation_flags["extraction_method"] = "dual_consensus"
            entries.append(entry)

    log.info("Consensus: %d pass1, %d pass2, %d consensus → %d valid entries",
             len(dicts1), len(dicts2), len(consensus_seqs), len(entries))
    return entries


def _call_llm_t2(paper_text: str) -> str:
    """Second LLM call at temperature=0.4 for consensus extraction."""
    try:
        from openai import OpenAI
    except ImportError:
        raise RuntimeError("openai package is required for LLM extraction")

    from .config import OPENAI_API_KEY, OPENAI_BASE_URL, OPENAI_MODEL

    kwargs: Dict[str, Any] = {"api_key": OPENAI_API_KEY}
    if OPENAI_BASE_URL:
        kwargs["base_url"] = OPENAI_BASE_URL

    client = OpenAI(**kwargs)

    max_chars = 48_000
    if len(paper_text) > max_chars:
        paper_text = paper_text[:max_chars] + "\n\n[TEXT TRUNCATED]"

    response = client.chat.completions.create(
        model=OPENAI_MODEL,
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": f"Extract all decoy peptide entries from this paper:\n\n{paper_text}"},
        ],
        temperature=0.4,
        max_tokens=4096,
    )
    return response.choices[0].message.content or "[]"


# ── Unified extraction entry point ──────────────────────────────────────

def extract_from_paper(
    paper: PaperRecord,
    consensus: bool = False,
) -> List[DecoyEntry]:
    """
    Extract decoy entries from a paper, using LLM if available,
    otherwise falling back to rule-based extraction.

    After extraction, every entry is checked against the source text
    to flag possible LLM hallucinations.

    Parameters
    ----------
    paper : PaperRecord
        Paper data from the Fetcher agent.
    consensus : bool
        If True, run dual-extraction consensus mode (2× API cost but
        much lower hallucination rate).

    Returns
    -------
    list[DecoyEntry]
        Extracted and validated entries.
    """
    # ── Try LLM first, fall back to rule-based extraction ──────────────
    entries: List[DecoyEntry] = []

    if not OPENAI_API_KEY:
        log.warning(
            "OPENAI_API_KEY not set — falling back to rule-based extraction "
            "for PMID %s (results will be lower quality)",
            paper.pmid,
        )
        entries = extract_rule_based(paper)
    else:
        try:
            if consensus:
                entries = _extract_consensus(paper)
            else:
                entries = extract_with_llm(paper)
        except Exception as exc:
            log.warning(
                "LLM extraction failed for PMID %s (%s) — "
                "falling back to rule-based extraction",
                paper.pmid, exc,
            )
            entries = extract_rule_based(paper)

    # ── Source-text verification ────────────────────────────────────────
    if entries:
        source_text = f"{paper.title}\n{paper.abstract}\n{paper.full_text}"
        entries = verify_peptide_in_source(entries, source_text)

    return entries
