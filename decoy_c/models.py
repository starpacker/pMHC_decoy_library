"""
Pydantic v2 models for the Decoy C Library standardised schema.

Every field mirrors the JSON schema provided in the specification, with
strict enum validation for evidence_level and expression_pattern.
"""

from __future__ import annotations

import re
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


# ── Enums ────────────────────────────────────────────────────────────────

class EvidenceLevel(str, Enum):
    """Risk evidence classification (must be exactly one of these four)."""
    LEVEL_1 = "Level_1_Clinical_Fatal"
    LEVEL_2 = "Level_2_In_Vitro_Confirmed"
    LEVEL_3 = "Level_3_High_Throughput_Screened"
    LEVEL_4 = "Level_4_In_Silico_High_Risk"


class ExpressionPattern(str, Enum):
    """How broadly the source protein is expressed in the human body."""
    UBIQUITOUS = "Ubiquitous"
    ESSENTIAL_ORGAN_SPECIFIC = "Essential_Organ_Specific"
    RESTRICTED = "Restricted"


# ── Sub-models ───────────────────────────────────────────────────────────

class ComputationalSafetyScore(BaseModel):
    """Algorithmic toxicity / safety scores."""
    ARDitox_score: Optional[float] = Field(
        None, description="ARDitox predicted safety score (lower = safer)"
    )
    other_scores: Dict[str, float] = Field(
        default_factory=dict,
        description="Additional tool scores, e.g. EpiTox",
    )


class PeptideInfo(BaseModel):
    """Core identity of the decoy peptide."""
    decoy_sequence: str = Field(
        ..., min_length=8, max_length=15,
        description="Off-target peptide sequence (8-15 AA)",
    )
    hla_allele: str = Field(
        ..., description="HLA allele in high-resolution format, e.g. HLA-A*02:01"
    )
    source_protein: str = Field(..., description="Full protein name, e.g. Titin")
    gene_symbol: str = Field(..., description="HGNC gene symbol, e.g. TTN")
    uniprot_id: Optional[str] = Field(
        None, description="UniProt accession, e.g. Q8WZ42"
    )

    @field_validator("decoy_sequence")
    @classmethod
    def uppercase_aa(cls, v: str) -> str:
        v = v.strip().upper()
        if not re.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", v):
            raise ValueError(f"Invalid amino-acid sequence: {v}")
        return v

    @field_validator("hla_allele")
    @classmethod
    def normalise_hla(cls, v: str) -> str:
        v = v.strip()
        # Accept HLA-X*NN:NN or HLA-X*NN:xx (low-res placeholder)
        if not re.match(r"^HLA-[A-Z]+\d?\*\d{2}:\w{2,3}$", v):
            m = re.match(r"^HLA-([A-Z]+)(\d+)$", v)
            if m:
                return f"HLA-{m.group(1)}*{int(m.group(2)):02d}:xx"
        return v


class DiscoveryContext(BaseModel):
    """How the decoy was originally discovered."""
    original_target_sequence: Optional[str] = Field(
        None, description="Intended target peptide sequence"
    )
    original_target_protein: Optional[str] = Field(
        None, description="Intended target protein name, e.g. MAGE-A3"
    )
    tcr_name_or_id: Optional[str] = Field(
        None, description="TCR clone name / identifier, e.g. a3a TCR, 1G4, DMF5"
    )


class RiskProfile(BaseModel):
    """Danger assessment for the decoy."""
    evidence_level: EvidenceLevel = Field(
        ..., description="Toxicity evidence level (strict 4-tier)"
    )
    critical_organs_affected: List[str] = Field(
        default_factory=list,
        description="Organs directly affected, e.g. Heart, Brain; or Systemic",
    )
    expression_pattern: Optional[ExpressionPattern] = Field(
        None, description="Source-protein expression breadth"
    )
    computational_safety_score: ComputationalSafetyScore = Field(
        default_factory=ComputationalSafetyScore,
    )


class ExperimentalEvidence(BaseModel):
    """Hard experimental data backing the entry."""
    mass_spec_confirmed: Optional[bool] = Field(
        None,
        description="Was the peptide detected by MS on cell-surface HLA?",
    )
    assays_performed: List[str] = Field(
        default_factory=list,
        description="Names of assays: X-scan, Cr-release, SPR, etc.",
    )
    cross_reactivity_affinity: Optional[str] = Field(
        None,
        description="Affinity metric, e.g. EC50 = 1.2 uM or KD = 50 nM",
    )


class Provenance(BaseModel):
    """Literature / trial traceability."""
    pmid: List[str] = Field(default_factory=list, description="PubMed IDs")
    clinical_trial_id: List[str] = Field(
        default_factory=list, description="ClinicalTrials.gov NCT IDs"
    )
    evidence_summary: str = Field(
        "", description="1-2 sentence English summary of why this entry exists"
    )


class Source(BaseModel):
    """Tracks exactly which paper / article / report this entry was extracted from."""
    title: str = Field(
        ..., description="Title of the paper, article, or report"
    )
    authors: Optional[str] = Field(
        None, description="Author list, e.g. 'Cameron BJ, Gerry AB, et al.'"
    )
    journal: Optional[str] = Field(
        None, description="Journal or publication name"
    )
    year: Optional[int] = Field(
        None, description="Publication year"
    )
    pmid: Optional[str] = Field(
        None, description="PubMed ID of the source paper"
    )
    doi: Optional[str] = Field(
        None, description="DOI of the source paper"
    )
    url: Optional[str] = Field(
        None, description="Direct URL to the paper or report"
    )
    citation: Optional[str] = Field(
        None,
        description="Full formatted citation string, e.g. 'Cameron et al., Sci Transl Med, 2013'"
    )


# ── Top-level entry ─────────────────────────────────────────────────────

class DecoyEntry(BaseModel):
    """
    A single record in the Decoy C Library.

    Maps 1-to-1 with the standardised JSON schema.
    """
    decoy_id: str = Field(..., pattern=r"^DC-\d{4,}$", description="Unique ID, e.g. DC-0001")
    peptide_info: PeptideInfo
    discovery_context: DiscoveryContext = Field(default_factory=DiscoveryContext)
    risk_profile: RiskProfile
    experimental_evidence: ExperimentalEvidence = Field(
        default_factory=ExperimentalEvidence,
    )
    provenance: Provenance = Field(default_factory=Provenance)
    source: Optional[Source] = Field(
        None,
        description="The paper / article / report from which this entry was extracted",
    )

    # ── Optional extraction metadata (not in the public schema) ──────
    thought_process: Optional[str] = Field(
        None,
        description="Agent's reasoning chain with inline citations (anti-hallucination)",
    )
    validation_flags: Dict[str, str] = Field(
        default_factory=dict,
        description="Flags set by the Validator agent, e.g. uniprot_match: OK",
    )


class DecoyLibrary(BaseModel):
    """The full collection – serialised to / from decoy_library.json."""
    entries: List[DecoyEntry] = Field(default_factory=list)
    version: str = "0.1.0"

    # ── helpers ───────────────────────────────────────────────────────
    @property
    def next_id(self) -> str:
        """Return the next sequential DC-NNNN id."""
        if not self.entries:
            return "DC-0001"
        max_num = max(int(e.decoy_id.split("-")[1]) for e in self.entries)
        return f"DC-{max_num + 1:04d}"

    def find_by_sequence(self, seq: str) -> Optional[DecoyEntry]:
        seq = seq.strip().upper()
        for e in self.entries:
            if e.peptide_info.decoy_sequence == seq:
                return e
        return None

    def add_entry(self, entry: DecoyEntry, deduplicate: bool = True) -> bool:
        """Add entry; returns False if duplicate sequence already exists."""
        if deduplicate and self.find_by_sequence(entry.peptide_info.decoy_sequence):
            return False
        self.entries.append(entry)
        return True
