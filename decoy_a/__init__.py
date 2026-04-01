"""
Decoy A/B Library — Computational Off-Target Screening Pipeline
================================================================
A funnel-based computational pipeline for identifying potential TCR
off-target peptides through two complementary strategies:

    Decoy A — Sequence homology (Hamming distance ≤ 2)
    Decoy B — Structural & physicochemical similarity (3D/electrostatic)

Pipeline Architecture (Funnel Strategy):
    Step 1: Human K-mer Dictionary + Tissue Expression Mapping
    Step 2: MHC Presentation Gate (NetMHCpan EL)
    Step 3A: Sequence-driven scan (Decoy A)
    Step 3B: Structure-driven scan (Decoy B)
    Step 4: Composite Risk Scoring & Final Ranking
"""

__version__ = "0.1.0"
