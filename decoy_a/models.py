"""
Data models for the Decoy A/B computational pipeline.

These models represent the intermediate and final data structures used
throughout the funnel-based screening workflow.  They are intentionally
separate from the Decoy C ``models.py`` to keep concerns clean, but the
final ``DecoyABEntry`` can be converted into a Decoy C ``DecoyEntry``
for unified library storage.
"""

from __future__ import annotations

import re
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


# ── Enums ────────────────────────────────────────────────────────────────

class DecoySource(str, Enum):
    """Which branch of the A/B pipeline produced this hit."""
    DECOY_A = "Decoy_A_Sequence_Homology"
    DECOY_B = "Decoy_B_Structural_Similarity"
    BOTH = "Both_A_and_B"


class BindingStrength(str, Enum):
    """NetMHCpan EL binding classification."""
    STRONG = "Strong_Binder"
    WEAK = "Weak_Binder"
    NON_BINDER = "Non_Binder"


# ── K-mer record (Step 1 output) ────────────────────────────────────────

class KmerRecord(BaseModel):
    """A single peptide k-mer derived from the human proteome."""
    sequence: str = Field(..., min_length=8, max_length=15)
    length: int = Field(..., ge=8, le=15)
    source_proteins: List[str] = Field(
        default_factory=list,
        description="UniProt accessions of proteins containing this k-mer",
    )
    gene_symbols: List[str] = Field(
        default_factory=list,
        description="HGNC gene symbols of source proteins",
    )

    @field_validator("sequence")
    @classmethod
    def uppercase_aa(cls, v: str) -> str:
        v = v.strip().upper()
        if not re.match(r"^[ACDEFGHIKLMNPQRSTVWY]+$", v):
            raise ValueError(f"Invalid amino-acid sequence: {v}")
        return v


# ── Expression profile (Step 1 output) ──────────────────────────────────

class TissueExpression(BaseModel):
    """RNA expression profile for a gene across vital organs."""
    gene_symbol: str
    tissue_tpm: Dict[str, float] = Field(
        default_factory=dict,
        description="Mapping of tissue name -> TPM value",
    )
    max_vital_organ_tpm: float = Field(
        0.0, description="Highest TPM across vital organs"
    )
    vital_organ_expressed: bool = Field(
        False, description="True if TPM > threshold in any vital organ"
    )
    expression_category: str = Field(
        "unknown",
        description="high_risk | expressed | silent | restricted",
    )


# ── HLA-filtered peptide (Step 2 output) ────────────────────────────────

class HLABindingResult(BaseModel):
    """NetMHCpan EL prediction result for a single peptide-HLA pair."""
    sequence: str = Field(..., min_length=8, max_length=15)
    hla_allele: str
    el_rank: float = Field(..., description="NetMHCpan EL %Rank")
    el_score: float = Field(0.0, description="Raw EL score")
    binding_strength: BindingStrength = BindingStrength.NON_BINDER
    gene_symbols: List[str] = Field(default_factory=list)
    source_proteins: List[str] = Field(default_factory=list)


# ── Decoy A hit (Step 3A output) ────────────────────────────────────────

class MismatchDetail(BaseModel):
    """Describes a single amino-acid mismatch between target and candidate."""
    position: int = Field(..., description="0-indexed position in the peptide")
    target_aa: str = Field(..., max_length=1)
    candidate_aa: str = Field(..., max_length=1)
    is_anchor: bool = Field(False, description="True if position is an HLA anchor")
    is_tcr_contact: bool = Field(False, description="True if position faces the TCR")


class DecoyAHit(BaseModel):
    """A Decoy A candidate: sequence-homologous peptide within Hamming distance."""
    sequence: str = Field(..., min_length=8, max_length=15)
    target_sequence: str = Field(..., description="The target peptide being screened")
    hamming_distance: int = Field(..., ge=0, le=5)
    mismatches: List[MismatchDetail] = Field(default_factory=list)
    n_tcr_contact_mismatches: int = Field(0)
    n_anchor_mismatches: int = Field(0)
    el_rank: float = Field(..., description="NetMHCpan EL %Rank")
    hla_allele: str = Field(...)
    gene_symbols: List[str] = Field(default_factory=list)
    source_proteins: List[str] = Field(default_factory=list)
    expression: Optional[TissueExpression] = None
    similarity_score: float = Field(
        0.0, description="Sequence similarity score (higher = more similar)"
    )


# ── Decoy B hit (Step 3B output) ────────────────────────────────────────

class PhysicochemFeatures(BaseModel):
    """Atchley/Kidera factor vector for TCR-contact residues."""
    contact_residues: str = Field(
        ..., description="Amino acids at TCR contact positions"
    )
    feature_vector: List[float] = Field(
        default_factory=list,
        description="Flattened physicochemical feature vector",
    )
    cosine_similarity: float = Field(
        0.0, description="Cosine similarity to target peptide features"
    )


class StructuralScore(BaseModel):
    """3D structural and electrostatic similarity assessment."""
    modeling_tool: str = Field("", description="Tool used: biopython+dual_superposition | none | unavailable | error")
    pdb_path: Optional[str] = Field(None, description="Path to modeled pMHC PDB")
    surface_correlation: float = Field(
        0.0, description="Combined structural similarity score [0-1]"
    )
    rmsd: Optional[float] = Field(None, description="Peptide backbone RMSD to target pMHC (A)")
    electrostatic_fingerprint_distance: Optional[float] = Field(
        None, description="Distance in electrostatic fingerprint space"
    )
    # Phase 1: Interface descriptors
    plip_tanimoto: Optional[float] = Field(
        None, description="PLIP non-covalent interaction Tanimoto similarity [0-1]"
    )
    bsa_target: Optional[float] = Field(
        None, description="Target pMHC buried surface area (A^2)"
    )
    bsa_candidate: Optional[float] = Field(
        None, description="Candidate pMHC buried surface area (A^2)"
    )
    bsa_similarity: Optional[float] = Field(
        None, description="BSA similarity [0-1]"
    )
    prodigy_dg_target: Optional[float] = Field(
        None, description="Target PRODIGY dG prediction (kcal/mol)"
    )
    prodigy_dg_candidate: Optional[float] = Field(
        None, description="Candidate PRODIGY dG prediction (kcal/mol)"
    )
    prodigy_similarity: Optional[float] = Field(
        None, description="PRODIGY dG similarity [0-1]"
    )
    interface_combined: Optional[float] = Field(
        None, description="Weighted combination of all interface descriptors [0-1]"
    )
    # Phase 2: APBS + PeSTo descriptors
    esp_similarity: Optional[float] = Field(
        None, description="APBS electrostatic potential similarity [0-1]"
    )
    pesto_similarity: Optional[float] = Field(
        None, description="PeSTo interface embedding similarity [0-1]"
    )
    # Cross-validation (Boltz)
    boltz_pdb_path: Optional[str] = Field(
        None, description="Path to Boltz-predicted pMHC structure"
    )
    boltz_confidence: Optional[float] = Field(
        None, description="Boltz confidence score (0.8*pLDDT + 0.2*PTM) [0-1]"
    )
    boltz_iptm: Optional[float] = Field(
        None, description="Boltz interface TM-score [0-1]"
    )
    boltz_rmsd: Optional[float] = Field(
        None, description="Peptide backbone RMSD between Boltz and primary model (A)"
    )
    cross_validation_agreement: Optional[float] = Field(
        None, description="Agreement score between primary (tFold/AF3) and Boltz predictions [0-1]"
    )


class DecoyBHit(BaseModel):
    """A Decoy B candidate: structurally / physicochemically similar peptide."""
    sequence: str = Field(..., min_length=8, max_length=15)
    target_sequence: str = Field(..., description="The target peptide being screened")
    hamming_distance: int = Field(..., ge=0)
    el_rank: float = Field(..., description="NetMHCpan EL %Rank")
    hla_allele: str = Field(...)
    gene_symbols: List[str] = Field(default_factory=list)
    source_proteins: List[str] = Field(default_factory=list)
    expression: Optional[TissueExpression] = None
    physicochemical: PhysicochemFeatures
    structural: Optional[StructuralScore] = None
    similarity_score: float = Field(
        0.0, description="Combined physicochemical + structural similarity"
    )
    # MPNN inverse design metadata
    mpnn_source: str = Field(
        "",
        description=(
            "Origin of this hit in the MPNN branch: "
            "'proteome_matched' (Tier 1: exact human proteome match) | "
            "'hla_qualified_synthetic' (Tier 2: no proteome match but "
            "mhcflurry predicts HLA-presentable) | '' (not from MPNN)"
        ),
    )
    mpnn_score: float = Field(
        0.0,
        description=(
            "ProteinMPNN negative log-probability score "
            "(lower = higher confidence that the sequence folds correctly)"
        ),
    )


# ── Final ranked entry (Step 4 output) ──────────────────────────────────

class DecoyABEntry(BaseModel):
    """
    Final ranked entry combining all pipeline signals.

    This is the primary output of the Decoy A/B pipeline — one record
    per high-risk off-target peptide candidate.
    """
    decoy_ab_id: str = Field(
        ..., pattern=r"^DAB-\d{4,}$",
        description="Unique ID for Decoy A/B entries, e.g. DAB-0001",
    )
    sequence: str = Field(..., min_length=8, max_length=15)
    target_sequence: str
    hla_allele: str
    source: DecoySource
    gene_symbols: List[str] = Field(default_factory=list)
    source_proteins: List[str] = Field(default_factory=list)

    # Funnel metrics
    el_rank: float
    hamming_distance: int
    sequence_similarity: float = Field(0.0)
    physicochemical_similarity: float = Field(0.0)
    structural_similarity: float = Field(0.0)

    # Interface descriptor details
    plip_tanimoto: Optional[float] = Field(None, description="PLIP interaction Tanimoto [0-1]")
    bsa_similarity: Optional[float] = Field(None, description="BSA similarity [0-1]")
    prodigy_similarity: Optional[float] = Field(None, description="PRODIGY dG similarity [0-1]")
    esp_similarity: Optional[float] = Field(None, description="APBS electrostatic potential similarity [0-1]")
    pesto_similarity: Optional[float] = Field(None, description="PeSTo interface embedding similarity [0-1]")
    interface_combined: Optional[float] = Field(None, description="Interface descriptor combined [0-1]")

    # Expression risk
    expression: Optional[TissueExpression] = None
    vital_organ_tpm_weight: float = Field(1.0)
    critical_organs: List[str] = Field(default_factory=list)

    # Composite risk
    total_risk_score: float = Field(
        0.0, description="Final composite risk score (higher = more dangerous)"
    )
    risk_rank: int = Field(0, description="1-indexed rank in final output")

    # Mismatch annotation (from Decoy A)
    mismatches: List[MismatchDetail] = Field(default_factory=list)

    # Structural annotation (from Decoy B)
    structural: Optional[StructuralScore] = None

    @field_validator("sequence")
    @classmethod
    def uppercase_aa(cls, v: str) -> str:
        return v.strip().upper()


# ── Pipeline state container ─────────────────────────────────────────────

class PipelineState(BaseModel):
    """
    Tracks the full state of a Decoy A/B pipeline run.

    Useful for checkpointing, resuming, and auditing.
    """
    target_sequence: str
    target_protein: str = ""
    hla_allele: str
    pipeline_version: str = "0.1.0"

    # Step 1
    total_kmers_generated: int = 0
    genes_with_expression: int = 0

    # Step 2
    hla_filtered_count: int = 0

    # Step 3A
    decoy_a_hits: List[DecoyAHit] = Field(default_factory=list)

    # Step 3B
    physicochemical_candidates: int = 0
    structural_candidates: int = 0
    decoy_b_hits: List[DecoyBHit] = Field(default_factory=list)

    # Step 4
    final_entries: List[DecoyABEntry] = Field(default_factory=list)
    final_top_n: int = 0
