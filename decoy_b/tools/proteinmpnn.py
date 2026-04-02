"""
ProteinMPNN Wrapper — Inverse Sequence Design
==============================================
Uses ProteinMPNN (Baker Lab) to design peptide sequences that maintain
the structural features of a target pMHC's TCR-facing surface.

Setup:
    1. git clone https://github.com/dauparas/ProteinMPNN.git
    2. Weights are included in the repo (vanilla_model_weights/)
    3. Set PROTEINMPNN_DIR env variable

Public API:
    check_available()           -> bool
    design_peptide(pdb, ...)    -> list[MPNNDesign]
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

log = logging.getLogger(__name__)

# ── Configuration ───────────────────────────────────────────────────────
from decoy_a.config import PROTEINMPNN_DIR

import sys
# Add local external path to sys.path
_EXTERNAL_DIR = Path(__file__).resolve().parent.parent / "external"
if str(_EXTERNAL_DIR) not in sys.path:
    sys.path.insert(0, str(_EXTERNAL_DIR))

PROTEINMPNN_WEIGHTS = PROTEINMPNN_DIR / "vanilla_model_weights"
PROTEINMPNN_RUN = PROTEINMPNN_DIR / "protein_mpnn_run.py"
PROTEINMPNN_PARSE = PROTEINMPNN_DIR / "helper_scripts" / "parse_multiple_chains.py"
PROTEINMPNN_FIXED = PROTEINMPNN_DIR / "helper_scripts" / "make_fixed_positions_dict.py"
PROTEINMPNN_CHAINS = PROTEINMPNN_DIR / "helper_scripts" / "assign_fixed_chains.py"


@dataclass
class MPNNDesign:
    """A single sequence design from ProteinMPNN."""
    sequence: str              # Designed peptide sequence
    score: float               # Negative log-probability (lower = better)
    global_score: float        # Score across all residues
    seq_recovery: float        # Fraction matching input sequence
    sample_index: int = 0
    temperature: float = 0.1


def check_available() -> bool:
    """Check if ProteinMPNN is installed."""
    return PROTEINMPNN_RUN.exists() and PROTEINMPNN_WEIGHTS.exists()


def _parse_pdb(pdb_path: Path, output_dir: Path) -> Path:
    """Parse a PDB file into ProteinMPNN's JSONL format."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "python", str(PROTEINMPNN_PARSE),
        "--input_path", str(pdb_path),
        "--output_path", str(output_dir),
    ]

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    if proc.returncode != 0:
        raise RuntimeError(f"PDB parsing failed: {proc.stderr[:300]}")

    # The output is a JSONL file
    jsonl_files = list(output_dir.glob("*.jsonl"))
    if not jsonl_files:
        raise RuntimeError("No JSONL output from PDB parsing")

    return output_dir


def _write_fixed_positions(
    peptide_chain_id: str,
    peptide_length: int,
    anchor_positions: List[int],
    output_path: Path,
    pdb_name: str,
) -> Path:
    """
    Create fixed_positions JSONL for ProteinMPNN.

    Anchors (e.g., p2/p9 for HLA-A*02:01) are fixed; TCR contact
    positions are designable.

    Format: {pdb_name: {chain_id: [1-indexed positions to FIX]}}
    """
    fixed = list(anchor_positions)  # 1-indexed positions to fix

    data = {pdb_name: {peptide_chain_id: fixed}}
    output_path.write_text(json.dumps(data) + "\n", encoding="utf-8")
    return output_path


def _write_chain_design(
    peptide_chain_id: str,
    all_chain_ids: List[str],
    output_path: Path,
    pdb_name: str,
) -> Path:
    """
    Create chain_id JSONL specifying which chains to design.

    Only the peptide chain is designable; MHC chains are fixed.
    """
    design_chains = [peptide_chain_id]
    data = {pdb_name: design_chains}
    output_path.write_text(json.dumps(data) + "\n", encoding="utf-8")
    return output_path


def _parse_mpnn_output(fasta_path: Path) -> List[MPNNDesign]:
    """Parse ProteinMPNN FASTA output into MPNNDesign objects."""
    designs = []
    current_header = ""
    current_seq = ""

    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if current_seq and current_header:
                designs.append(_header_to_design(current_header, current_seq))
            current_header = line[1:]
            current_seq = ""
        else:
            current_seq += line.strip()

    if current_seq and current_header:
        designs.append(_header_to_design(current_header, current_seq))

    return designs


def _header_to_design(header: str, full_sequence: str) -> MPNNDesign:
    """Parse MPNN FASTA header and extract the peptide chain sequence."""
    score = 0.0
    global_score = 0.0
    seq_recovery = 0.0
    sample_idx = 0
    temperature = 0.1

    # Parse header fields: "name, score=X, global_score=X, ..."
    for part in header.split(","):
        part = part.strip()
        if "=" in part:
            key, val = part.split("=", 1)
            key = key.strip()
            val = val.strip()
            try:
                if key == "score":
                    score = float(val)
                elif key == "global_score":
                    global_score = float(val)
                elif key == "seq_recovery":
                    seq_recovery = float(val)
                elif key == "sample":
                    sample_idx = int(val)
                elif key == "T":
                    temperature = float(val)
            except ValueError:
                pass

    # full_sequence contains all chains separated by "/"
    # The peptide chain is typically the last (shortest) one
    chains = full_sequence.split("/")
    peptide_seq = min(chains, key=len) if chains else full_sequence

    return MPNNDesign(
        sequence=peptide_seq,
        score=score,
        global_score=global_score,
        seq_recovery=seq_recovery,
        sample_index=sample_idx,
        temperature=temperature,
    )


def design_peptide(
    pdb_path: str,
    peptide_chain_id: str = "C",
    anchor_positions: Optional[List[int]] = None,
    num_designs: int = 100,
    temperatures: Optional[List[float]] = None,
    model_name: str = "v_48_020",
) -> List[MPNNDesign]:
    """
    Design peptide sequences that preserve the pMHC structural features.

    Parameters
    ----------
    pdb_path : str
        Path to the pMHC complex PDB (from tFold/AF3).
    peptide_chain_id : str
        Chain ID of the peptide in the PDB (default: "C").
    anchor_positions : list[int], optional
        1-indexed positions to keep fixed (HLA anchors).
        Default for HLA-A*02:01 9-mer: [2, 9] (p2 and p9).
    num_designs : int
        Number of sequences to generate.
    temperatures : list[float], optional
        Sampling temperatures. Default: [0.1, 0.2, 0.3].
    model_name : str
        ProteinMPNN model variant.

    Returns
    -------
    list[MPNNDesign]
        Designed peptide sequences with scores.
    """
    if anchor_positions is None:
        anchor_positions = [2, 9]  # HLA-A*02:01 9-mer anchors

    if temperatures is None:
        temperatures = [0.1, 0.2, 0.3]

    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    work_dir = Path(tempfile.mkdtemp(prefix="mpnn_"))

    try:
        # Step 1: Parse PDB
        parsed_dir = work_dir / "parsed"
        _parse_pdb(pdb_path, parsed_dir)

        pdb_name = pdb_path.stem

        # Step 2: Write fixed positions (anchor residues)
        fixed_path = work_dir / "fixed_positions.jsonl"
        _write_fixed_positions(
            peptide_chain_id, 9, anchor_positions, fixed_path, pdb_name,
        )

        # Step 3: Write chain design spec (only peptide chain is designable)
        chain_path = work_dir / "chain_ids.jsonl"
        _write_chain_design(
            peptide_chain_id, ["A", "B", peptide_chain_id],
            chain_path, pdb_name,
        )

        # Step 4: Run ProteinMPNN
        output_dir = work_dir / "output"
        output_dir.mkdir()

        temp_str = " ".join(str(t) for t in temperatures)

        cmd = [
            "python", str(PROTEINMPNN_RUN),
            "--jsonl_path", str(parsed_dir),
            "--chain_id_jsonl", str(chain_path),
            "--fixed_positions_jsonl", str(fixed_path),
            "--out_folder", str(output_dir),
            "--num_seq_per_target", str(num_designs),
            "--sampling_temp", temp_str,
            "--seed", "42",
            "--batch_size", "1",
            "--path_to_model_weights", str(PROTEINMPNN_WEIGHTS),
            "--model_name", model_name,
        ]

        log.info("Running ProteinMPNN: %d designs at T=%s", num_designs, temp_str)

        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )

        if proc.returncode != 0:
            raise RuntimeError(f"ProteinMPNN failed: {proc.stderr[:500]}")

        # Step 5: Parse output
        fasta_dir = output_dir / "seqs"
        fasta_files = list(fasta_dir.glob("*.fa")) if fasta_dir.exists() else []

        all_designs = []
        for fasta_file in fasta_files:
            designs = _parse_mpnn_output(fasta_file)
            all_designs.extend(designs)

        # Deduplicate by sequence
        seen = set()
        unique_designs = []
        for d in all_designs:
            if d.sequence not in seen:
                seen.add(d.sequence)
                unique_designs.append(d)

        # Sort by score (lower = better)
        unique_designs.sort(key=lambda d: d.score)

        log.info(
            "ProteinMPNN: %d total designs, %d unique sequences",
            len(all_designs), len(unique_designs),
        )

        return unique_designs

    finally:
        # Cleanup temp directory
        import shutil
        shutil.rmtree(work_dir, ignore_errors=True)


def design_and_filter(
    pdb_path: str,
    target_hla: str = "HLA-A*02:01",
    peptide_chain_id: str = "C",
    anchor_positions: Optional[List[int]] = None,
    num_designs: int = 1000,
    temperatures: Optional[List[float]] = None,
    el_rank_threshold: float = 2.0,
) -> List[MPNNDesign]:
    """
    Design sequences, then filter for HLA binding using mhcflurry.

    This is the standalone version of the HLA qualification filter.
    The full two-tier pipeline (proteome match + HLA filter) is in
    ``decoy_b.scanner.run_mpnn_design()``.

    Parameters
    ----------
    pdb_path : str
        pMHC structure PDB.
    target_hla : str
        HLA allele for binding prediction.
    peptide_chain_id : str
        Peptide chain in the PDB.
    anchor_positions : list[int], optional
        Fixed anchor positions.
    num_designs : int
        Number of initial designs.
    temperatures : list[float], optional
        Sampling temperatures.
    el_rank_threshold : float
        Max EL %Rank to qualify as HLA-presentable (default 2.0).

    Returns
    -------
    list[MPNNDesign]
        Designs that pass HLA binding filter.
    """
    # Generate designs
    designs = design_peptide(
        pdb_path, peptide_chain_id, anchor_positions,
        num_designs, temperatures,
    )

    if not designs:
        return []

    peptides = [d.sequence for d in designs]

    # Try mhcflurry (preferred — always available in our deployment)
    try:
        from decoy_a.tools.mhcflurry import (
            check_available as mhcf_available,
            predict_binding as mhcf_predict,
        )

        if mhcf_available():
            results = mhcf_predict(peptides, target_hla)
            binder_seqs = {
                r.sequence for r in results
                if r.el_rank <= el_rank_threshold
            }
            filtered = [d for d in designs if d.sequence in binder_seqs]
            log.info(
                "MPNN + mhcflurry filter: %d designs -> %d HLA binders "
                "(EL%%Rank ≤ %.1f)",
                len(designs), len(filtered), el_rank_threshold,
            )
            return filtered
    except ImportError:
        pass

    # Fallback: try NetMHCpan
    try:
        from decoy_a.tools.netmhcpan import (
            check_available as nmhc_available,
            predict_binding,
        )

        if nmhc_available():
            results = predict_binding(peptides, target_hla)
            binder_seqs = {
                r.sequence for r in results
                if r.el_rank <= el_rank_threshold
            }
            filtered = [d for d in designs if d.sequence in binder_seqs]
            log.info(
                "MPNN + NetMHCpan filter: %d designs -> %d HLA binders "
                "(EL%%Rank ≤ %.1f)",
                len(designs), len(filtered), el_rank_threshold,
            )
            return filtered
    except ImportError:
        pass

    log.warning("No HLA predictor available; returning unfiltered MPNN designs")
    return designs
