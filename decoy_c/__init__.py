"""
Decoy C Library — Cold-Start Collection System
================================================
A multi-agent pipeline for collecting, extracting, validating, and curating
a standardized database of off-target / cross-reactive peptide "decoys"
that have caused (or are predicted to cause) toxicity in TCR-based
immunotherapies.

Agents:
    Fetcher   – retrieves papers from PubMed / PMC via NCBI E-utilities.
    Extractor – uses an LLM to produce structured DecoyEntry JSON from text.
    Validator – cross-checks entries against UniProt & IEDB public APIs.
"""

__version__ = "0.1.0"
