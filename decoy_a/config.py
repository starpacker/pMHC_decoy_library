"""
Configuration for the Decoy A/B computational pipeline.

Defines paths, thresholds, and external tool settings used across
the funnel-based off-target screening workflow.
"""

from __future__ import annotations

import os
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"

# Decoy A data (shared base + A-specific outputs)
A_DATA_DIR = DATA_DIR / "decoy_a"

# Decoy B data
B_DATA_DIR = DATA_DIR / "decoy_b"

# Legacy alias (code that still references AB_DATA_DIR)
AB_DATA_DIR = A_DATA_DIR

# Pre-built assets (generated once, reused across projects)
KMER_DB_PATH = A_DATA_DIR / "human_kmer_db.parquet"
EXPRESSION_DB_PATH = A_DATA_DIR / "gene_expression.parquet"
HLA_FILTERED_PATH = A_DATA_DIR / "hla_filtered_kmers.parquet"

# Pipeline outputs
DECOY_A_OUTPUT = A_DATA_DIR / "decoy_a_results.json"
DECOY_B_OUTPUT = B_DATA_DIR / "decoy_b_results.json"
FINAL_RANKED_OUTPUT = B_DATA_DIR / "final_ranked_decoys.json"

# ── UniProt Settings ─────────────────────────────────────────────────────
UNIPROT_PROTEOME_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
)
UNIPROT_PROTEOME_PARAMS = {
    "query": "organism_id:9606 AND reviewed:true",
    "format": "fasta",
}
# Local cache of the Swiss-Prot FASTA
SWISSPROT_FASTA = AB_DATA_DIR / "uniprot_sprot_human.fasta"

# ── K-mer Generation ─────────────────────────────────────────────────────
KMER_LENGTHS = [8, 9, 10, 11]  # Standard pMHC-I peptide lengths

# ── Tissue Expression (GTEx / HPA) ──────────────────────────────────────
# Human Protein Atlas downloadable data
HPA_NORMAL_TISSUE_URL = (
    "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
)
HPA_RNA_TISSUE_URL = (
    "https://v23.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
)
HPA_NORMAL_TISSUE_TSV = AB_DATA_DIR / "normal_tissue.tsv"
HPA_RNA_TISSUE_TSV = AB_DATA_DIR / "rna_tissue_consensus.tsv"

# Critical organs for toxicity assessment (high weight)
VITAL_ORGANS = [
    "heart muscle",
    "cerebral cortex",
    "lung",
    "liver",
    "kidney",
    "colon",
    "small intestine",
]

# Organs with low clinical relevance for off-target (low weight)
LOW_RELEVANCE_TISSUES = [
    "testis",
    "placenta",
]

# Expression thresholds
TPM_HIGH_RISK_THRESHOLD = 10.0   # TPM > 10 in vital organ = high weight
TPM_EXPRESSED_THRESHOLD = 1.0    # TPM > 1 = gene is expressed
TPM_SILENT_THRESHOLD = 1.0       # TPM < 1 in ALL vital organs = low weight

# ── MHC / HLA Prediction ────────────────────────────────────────────────
# NetMHCpan 4.1 EL mode settings
NETMHCPAN_CMD = os.getenv("NETMHCPAN_CMD", "netMHCpan")
NETMHCPAN_EL_RANK_THRESHOLD = 2.0   # %Rank <= 2.0 (weak binder + strong)
NETMHCPAN_STRONG_THRESHOLD = 0.5    # %Rank <= 0.5 for strong binders
NETMHCPAN_BATCH_SIZE = 50000        # Peptides per batch call

# Default HLA allele (most studied)
DEFAULT_HLA_ALLELE = "HLA-A*02:01"

# ── Decoy A: Sequence Homology ──────────────────────────────────────────
HAMMING_DISTANCE_MAX = 2

# Anchor positions for HLA-A*02:01 (0-indexed): p2 and p9 (for 9-mers)
# Mutations at anchor positions that still pass EL threshold are included
ANCHOR_POSITIONS_A0201 = {1, 8}  # 0-indexed: position 2 and 9

# TCR contact residues (0-indexed, for 9-mers): p4-p8
TCR_CONTACT_POSITIONS_9MER = [3, 4, 5, 6, 7]

# ── Decoy B: Structural / Physicochemical Similarity ────────────────────
# Atchley factors for the 20 standard amino acids
# Each AA maps to 5 factors: [polarity, secondary_structure, molecular_size,
#                              codon_diversity, electrostatic_charge]
ATCHLEY_FACTORS: dict[str, list[float]] = {
    "A": [-0.591, -1.302, -0.733,  1.570, -0.146],
    "C": [-1.343,  0.465, -0.862, -1.020, -0.255],
    "D": [ 1.050,  0.302, -3.656, -0.259, -3.242],
    "E": [ 1.357, -1.453,  1.477,  0.113, -0.837],
    "F": [-1.006, -0.590,  1.891, -0.397,  0.412],
    "G": [-0.384,  1.652,  1.330,  1.045,  2.064],
    "H": [ 0.336, -0.417, -1.673, -1.474, -0.078],
    "I": [-1.239, -0.547,  2.131,  0.393,  0.816],
    "K": [ 1.831, -0.561,  0.533, -0.277,  1.648],
    "L": [-1.019, -0.987, -1.505,  1.266, -0.912],
    "M": [-0.663, -1.524,  2.219, -1.005,  1.212],
    "N": [ 0.945,  0.828,  1.299, -0.169,  0.933],
    "P": [ 0.189,  2.081, -1.628,  0.421, -1.392],
    "Q": [ 0.931, -0.179, -3.005, -0.503, -1.853],
    "R": [ 1.538, -0.055,  1.502,  0.440,  2.897],
    "S": [-0.228,  1.399, -4.760,  0.670, -2.647],
    "T": [-0.032,  0.326,  2.213,  0.908,  1.313],
    "V": [-1.337, -0.279, -0.544,  1.242, -1.262],
    "W": [-0.595,  0.009,  0.672, -2.128, -0.184],
    "Y": [ 0.260,  0.830,  3.097, -0.838,  1.512],
}

# Physicochemical similarity threshold (cosine similarity)
PHYSICOCHEMICAL_COSINE_THRESHOLD = 0.70
PHYSICOCHEMICAL_TOP_K = 5000  # Keep top-K for structural modeling

# Surface electrostatic similarity threshold
SURFACE_SIMILARITY_THRESHOLD = 0.75

# ── Decoy B: 3D Modeling Tools ──────────────────────────────────────────
# tFold (preferred for bulk pMHC prediction)
TFOLD_DIR = Path(os.getenv("TFOLD_DIR", os.path.expanduser("~/tools/tfold")))
TFOLD_WEIGHTS_DIR = Path(os.getenv("TFOLD_WEIGHTS_DIR", str(TFOLD_DIR / "weights")))

# AlphaFold3 (for high-accuracy refinement)
AF3_DIR = Path(os.getenv("AF3_DIR", os.path.expanduser("~/tools/alphafold3")))
AF3_MODEL_DIR = Path(os.getenv("AF3_MODEL_DIR", str(AF3_DIR / "models")))
AF3_DB_DIR = Path(os.getenv("AF3_DB_DIR", str(AF3_DIR / "databases")))

# ProteinMPNN (for inverse sequence design)
PROTEINMPNN_DIR = Path(os.getenv("PROTEINMPNN_DIR", os.path.expanduser("~/tools/ProteinMPNN")))

# Legacy (PANDORA/APE-Gen — still supported as fallback)
MODELING_TOOL = os.getenv("PMHC_MODELING_TOOL", "tfold")  # tfold | af3 | pandora | apegen
PANDORA_CMD = os.getenv("PANDORA_CMD", "pandora")
APEGEN_CMD = os.getenv("APEGEN_CMD", "APE-Gen")

# ── Final Risk Scoring ──────────────────────────────────────────────────
# Weight multipliers for vital organ expression
ORGAN_TPM_WEIGHT_HIGH = 10.0   # Heart/Brain TPM > threshold
ORGAN_TPM_WEIGHT_LOW = 0.1     # Testis/Placenta only
ORGAN_TPM_WEIGHT_DEFAULT = 1.0

# Final output size
FINAL_TOP_N = 100
