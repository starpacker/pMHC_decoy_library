"""
Step 1: Human K-mer Dictionary & Tissue Expression Mapping
==========================================================
Builds the multi-omics base layer for the Decoy A/B funnel:

    1. Download Swiss-Prot (reviewed) human proteome from UniProt
    2. Sliding-window generation of all 8/9/10/11-mer peptides
    3. De-duplication with gene/protein provenance tracking
    4. Tissue expression mapping via Human Protein Atlas (HPA)

This step runs once per server; subsequent pipeline runs reuse the
generated Parquet databases.

Public API
----------
    build_kmer_database()       -> Path  (to Parquet)
    build_expression_database() -> Path  (to Parquet)
    load_kmer_db()              -> pd.DataFrame
    load_expression_db()        -> pd.DataFrame
"""

from __future__ import annotations

import io
import logging
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import requests

from .config import (
    AB_DATA_DIR,
    HPA_RNA_TISSUE_TSV,
    HPA_RNA_TISSUE_URL,
    KMER_DB_PATH,
    KMER_LENGTHS,
    EXPRESSION_DB_PATH,
    LOW_RELEVANCE_TISSUES,
    SWISSPROT_FASTA,
    TPM_EXPRESSED_THRESHOLD,
    TPM_HIGH_RISK_THRESHOLD,
    TPM_SILENT_THRESHOLD,
    UNIPROT_PROTEOME_URL,
    UNIPROT_PROTEOME_PARAMS,
    VITAL_ORGANS,
)

log = logging.getLogger(__name__)


# ── FASTA Parsing ────────────────────────────────────────────────────────

def _parse_fasta(fasta_text: str) -> List[Tuple[str, str, str]]:
    """
    Parse a UniProt FASTA file.

    Returns
    -------
    list of (uniprot_id, gene_symbol, sequence)
    """
    entries = []
    current_id = ""
    current_gene = ""
    current_seq_parts: List[str] = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            # Save previous entry
            if current_id and current_seq_parts:
                entries.append((current_id, current_gene, "".join(current_seq_parts)))
            current_seq_parts = []

            # Parse header: >sp|Q8WZ42|TTN_HUMAN Titin OS=Homo sapiens ...
            # Extract UniProt ID
            parts = line.split("|")
            current_id = parts[1] if len(parts) >= 2 else ""

            # Extract gene name from GN=XXX
            gn_match = re.search(r"GN=(\S+)", line)
            current_gene = gn_match.group(1) if gn_match else ""
        else:
            current_seq_parts.append(line.strip())

    # Don't forget the last entry
    if current_id and current_seq_parts:
        entries.append((current_id, current_gene, "".join(current_seq_parts)))

    return entries


# ── UniProt Download ─────────────────────────────────────────────────────

def download_swissprot_fasta(force: bool = False) -> Path:
    """
    Download the reviewed human proteome (Swiss-Prot) from UniProt REST API.

    Parameters
    ----------
    force : bool
        Re-download even if local cache exists.

    Returns
    -------
    Path to the local FASTA file.
    """
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if SWISSPROT_FASTA.exists() and not force:
        log.info("Swiss-Prot FASTA already cached at %s", SWISSPROT_FASTA)
        return SWISSPROT_FASTA

    log.info("Downloading Swiss-Prot human proteome from UniProt...")
    resp = requests.get(
        UNIPROT_PROTEOME_URL,
        params=UNIPROT_PROTEOME_PARAMS,
        timeout=600,
        stream=True,
    )
    resp.raise_for_status()

    with open(SWISSPROT_FASTA, "w", encoding="utf-8") as f:
        for chunk in resp.iter_content(chunk_size=1024 * 1024, decode_unicode=True):
            f.write(chunk)

    size_mb = SWISSPROT_FASTA.stat().st_size / (1024 * 1024)
    log.info("Downloaded Swiss-Prot FASTA: %.1f MB -> %s", size_mb, SWISSPROT_FASTA)
    return SWISSPROT_FASTA


# ── K-mer Generation ─────────────────────────────────────────────────────

_VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def generate_kmers(
    sequence: str,
    lengths: List[int],
) -> List[str]:
    """
    Generate all k-mers of specified lengths from a protein sequence.

    Skips k-mers containing non-standard amino acids (X, U, B, Z, etc.).
    """
    kmers = []
    for k in lengths:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i : i + k]
            if all(aa in _VALID_AA for aa in kmer):
                kmers.append(kmer)
    return kmers


def build_kmer_database(force: bool = False) -> Path:
    """
    Build the de-duplicated human k-mer dictionary from Swiss-Prot.

    Creates a Parquet file with columns:
        sequence, length, source_proteins (list), gene_symbols (list)

    Uses a SQLite-backed disk approach for minimal memory footprint:
    k-mers are inserted into a temporary SQLite database, then exported
    to Parquet. This avoids holding ~11M k-mers in Python memory.

    Returns
    -------
    Path to the generated Parquet file.
    """
    try:
        import pandas as pd
    except ImportError:
        raise RuntimeError("pandas is required: pip install pandas pyarrow")

    if KMER_DB_PATH.exists() and not force:
        log.info("K-mer database already exists at %s", KMER_DB_PATH)
        return KMER_DB_PATH

    import gc
    import sqlite3

    # Step 1: Get the proteome
    fasta_path = download_swissprot_fasta(force=force)
    fasta_text = fasta_path.read_text(encoding="utf-8")
    proteins = _parse_fasta(fasta_text)
    log.info("Parsed %d proteins from Swiss-Prot", len(proteins))
    del fasta_text
    gc.collect()

    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    db_path = AB_DATA_DIR / "_kmer_build.db"
    if db_path.exists():
        db_path.unlink()

    # Step 2: Insert k-mers into SQLite (disk-based, minimal RAM)
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute("""
        CREATE TABLE kmer_raw (
            sequence TEXT NOT NULL,
            source_protein TEXT NOT NULL,
            gene_symbol TEXT NOT NULL
        )
    """)

    batch = []
    batch_size = 100_000
    total_inserted = 0

    for i, (uniprot_id, gene_symbol, seq) in enumerate(proteins):
        for kmer in generate_kmers(seq, KMER_LENGTHS):
            batch.append((kmer, uniprot_id, gene_symbol or ""))
            if len(batch) >= batch_size:
                conn.executemany(
                    "INSERT INTO kmer_raw VALUES (?, ?, ?)", batch
                )
                total_inserted += len(batch)
                batch.clear()

        if (i + 1) % 2000 == 0:
            if batch:
                conn.executemany(
                    "INSERT INTO kmer_raw VALUES (?, ?, ?)", batch
                )
                total_inserted += len(batch)
                batch.clear()
            conn.commit()
            log.info(
                "K-mer insert: %d/%d proteins, %d rows",
                i + 1, len(proteins), total_inserted,
            )

    if batch:
        conn.executemany("INSERT INTO kmer_raw VALUES (?, ?, ?)", batch)
        total_inserted += len(batch)
        batch.clear()
    conn.commit()
    del proteins
    gc.collect()

    log.info("Total k-mer rows inserted: %d. Creating index...", total_inserted)
    conn.execute("CREATE INDEX idx_seq ON kmer_raw(sequence)")
    conn.commit()

    # Step 3: Group and export in chunked Parquet files (low memory)
    log.info("Counting unique k-mers...")
    n_unique = conn.execute(
        "SELECT COUNT(DISTINCT sequence) FROM kmer_raw"
    ).fetchone()[0]
    log.info("Unique k-mers: %d", n_unique)

    import pyarrow as pa
    import pyarrow.parquet as pq

    log.info("Exporting to Parquet (streaming)...")
    cursor = conn.execute("""
        SELECT
            sequence,
            GROUP_CONCAT(DISTINCT source_protein) as prots,
            GROUP_CONCAT(DISTINCT gene_symbol) as genes
        FROM kmer_raw
        GROUP BY sequence
    """)

    chunk_export_size = 2_000_000
    chunk_seqs, chunk_lens, chunk_prots, chunk_genes = [], [], [], []
    export_count = 0
    chunk_files: List[Path] = []

    for seq, prots_str, genes_str in cursor:
        prot_list = sorted(set(prots_str.split(","))) if prots_str else []
        gene_list = sorted(set(g for g in genes_str.split(",") if g)) if genes_str else []
        chunk_seqs.append(seq)
        chunk_lens.append(len(seq))
        chunk_prots.append(prot_list)
        chunk_genes.append(gene_list)
        export_count += 1

        if len(chunk_seqs) >= chunk_export_size:
            table = pa.table({
                "sequence": chunk_seqs,
                "length": chunk_lens,
                "source_proteins": chunk_prots,
                "gene_symbols": chunk_genes,
            })
            chunk_path = AB_DATA_DIR / f"_kmer_part_{len(chunk_files):03d}.parquet"
            pq.write_table(table, chunk_path)
            chunk_files.append(chunk_path)
            chunk_seqs.clear()
            chunk_lens.clear()
            chunk_prots.clear()
            chunk_genes.clear()
            del table
            gc.collect()
            log.info("  Exported %d/%d k-mers (%d chunks)", export_count, n_unique, len(chunk_files))

    # Write remaining
    if chunk_seqs:
        table = pa.table({
            "sequence": chunk_seqs,
            "length": chunk_lens,
            "source_proteins": chunk_prots,
            "gene_symbols": chunk_genes,
        })
        chunk_path = AB_DATA_DIR / f"_kmer_part_{len(chunk_files):03d}.parquet"
        pq.write_table(table, chunk_path)
        chunk_files.append(chunk_path)
        del table, chunk_seqs, chunk_lens, chunk_prots, chunk_genes
        gc.collect()

    conn.close()
    log.info("Exported %d k-mers in %d chunk files", export_count, len(chunk_files))

    # Merge chunk files into final Parquet
    if len(chunk_files) == 1:
        chunk_files[0].rename(KMER_DB_PATH)
    else:
        log.info("Merging %d chunk files...", len(chunk_files))
        tables = [pq.read_table(f) for f in chunk_files]
        merged = pa.concat_tables(tables)
        pq.write_table(merged, KMER_DB_PATH)
        del tables, merged
        gc.collect()
        for f in chunk_files:
            f.unlink(missing_ok=True)

    # Cleanup SQLite temp
    db_path.unlink(missing_ok=True)

    final_size = KMER_DB_PATH.stat().st_size / (1024 * 1024)
    log.info(
        "K-mer database saved: %d unique peptides -> %s (%.1f MB)",
        export_count,
        KMER_DB_PATH,
        final_size,
    )
    return KMER_DB_PATH


# ── Expression Database ──────────────────────────────────────────────────

def download_hpa_expression(force: bool = False) -> Path:
    """
    Download the HPA RNA tissue consensus expression data.

    Returns
    -------
    Path to the local TSV file.
    """
    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)

    if HPA_RNA_TISSUE_TSV.exists() and not force:
        log.info("HPA expression data already cached at %s", HPA_RNA_TISSUE_TSV)
        return HPA_RNA_TISSUE_TSV

    log.info("Downloading HPA RNA tissue consensus data...")
    resp = requests.get(HPA_RNA_TISSUE_URL, timeout=300)
    resp.raise_for_status()

    # The download is a zip file containing the TSV
    with zipfile.ZipFile(io.BytesIO(resp.content)) as zf:
        names = zf.namelist()
        tsv_name = [n for n in names if n.endswith(".tsv")][0]
        with zf.open(tsv_name) as f:
            content = f.read()
        HPA_RNA_TISSUE_TSV.write_bytes(content)

    log.info("HPA expression data saved: %s", HPA_RNA_TISSUE_TSV)
    return HPA_RNA_TISSUE_TSV


def build_expression_database(force: bool = False) -> Path:
    """
    Build the gene-tissue expression database from HPA data.

    Creates a Parquet file with columns:
        gene_symbol, tissue, tpm

    And a pivoted summary with max vital organ TPM and expression category.

    Returns
    -------
    Path to the generated Parquet file.
    """
    try:
        import pandas as pd
    except ImportError:
        raise RuntimeError("pandas is required: pip install pandas pyarrow")

    if EXPRESSION_DB_PATH.exists() and not force:
        log.info("Expression database already exists at %s", EXPRESSION_DB_PATH)
        return EXPRESSION_DB_PATH

    # Download HPA data
    tsv_path = download_hpa_expression(force=force)

    # Parse: columns are Gene, Gene name, Tissue, TPM, nTPM
    df_raw = pd.read_csv(tsv_path, sep="\t")
    log.info("Loaded HPA expression data: %d rows", len(df_raw))

    # Rename columns for consistency
    # HPA format: Gene (ENSG), Gene name (HGNC symbol), Tissue, nTPM
    col_map = {}
    for col in df_raw.columns:
        col_lower = col.strip().lower()
        if col_lower == "gene name":
            col_map[col] = "gene_symbol"
        elif col_lower == "gene":
            col_map[col] = "ensembl_id"
        elif col_lower == "tissue":
            col_map[col] = "tissue"
        elif col_lower in ("tpm", "ntpm"):
            col_map[col] = "tpm"
    if "tpm" not in col_map.values():
        for col in df_raw.columns:
            if col.strip().lower() in ("ntpm",):
                col_map[col] = "tpm"
                break

    df_raw = df_raw.rename(columns=col_map)

    # Ensure we have required columns
    required = {"gene_symbol", "tissue", "tpm"}
    available = set(df_raw.columns)
    if not required.issubset(available):
        # Fallback: use positional columns (Gene, Gene name, Tissue, nTPM)
        log.warning(
            "HPA columns don't match expected format. Available: %s. "
            "Attempting positional mapping.",
            list(df_raw.columns),
        )
        if len(df_raw.columns) >= 4:
            df_raw.columns = ["ensembl_id", "gene_symbol", "tissue", "tpm"] + list(
                df_raw.columns[4:]
            )

    # Keep only needed columns and ensure TPM is numeric
    df = df_raw[["gene_symbol", "tissue", "tpm"]].copy()
    df["tpm"] = pd.to_numeric(df["tpm"], errors="coerce").fillna(0.0)
    df["tissue"] = df["tissue"].str.strip().str.lower()

    # Compute per-gene summary
    vital_set = {t.lower() for t in VITAL_ORGANS}
    low_set = {t.lower() for t in LOW_RELEVANCE_TISSUES}

    def categorize_gene(group):
        vital_tpm = group[group["tissue"].isin(vital_set)]["tpm"]
        max_vital = vital_tpm.max() if len(vital_tpm) > 0 else 0.0
        all_tpm = group["tpm"]

        if max_vital >= TPM_HIGH_RISK_THRESHOLD:
            category = "high_risk"
        elif max_vital >= TPM_EXPRESSED_THRESHOLD:
            category = "expressed"
        elif all_tpm.max() >= TPM_EXPRESSED_THRESHOLD:
            # Expressed somewhere, but not in vital organs
            expressed_tissues = set(
                group[group["tpm"] >= TPM_EXPRESSED_THRESHOLD]["tissue"]
            )
            if expressed_tissues.issubset(low_set):
                category = "restricted"
            else:
                category = "expressed"
        else:
            category = "silent"

        return pd.Series({
            "max_vital_organ_tpm": max_vital,
            "expression_category": category,
        })

    gene_summary = df.groupby("gene_symbol").apply(
        categorize_gene, include_groups=False,
    ).reset_index()

    # Merge summary back
    df = df.merge(gene_summary, on="gene_symbol", how="left")

    AB_DATA_DIR.mkdir(parents=True, exist_ok=True)
    df.to_parquet(EXPRESSION_DB_PATH, index=False, engine="pyarrow")

    n_genes = df["gene_symbol"].nunique()
    log.info(
        "Expression database saved: %d genes, %d tissue records -> %s",
        n_genes, len(df), EXPRESSION_DB_PATH,
    )
    return EXPRESSION_DB_PATH


# ── Loaders ──────────────────────────────────────────────────────────────

def load_kmer_db():
    """Load the k-mer database as a pandas DataFrame."""
    import pandas as pd

    if not KMER_DB_PATH.exists():
        raise FileNotFoundError(
            f"K-mer database not found at {KMER_DB_PATH}. "
            "Run `build_kmer_database()` first."
        )
    return pd.read_parquet(KMER_DB_PATH)


def load_expression_db():
    """Load the expression database as a pandas DataFrame."""
    import pandas as pd

    if not EXPRESSION_DB_PATH.exists():
        raise FileNotFoundError(
            f"Expression database not found at {EXPRESSION_DB_PATH}. "
            "Run `build_expression_database()` first."
        )
    return pd.read_parquet(EXPRESSION_DB_PATH)


def get_gene_expression(gene_symbol: str, expr_df=None) -> Optional[Dict]:
    """
    Get tissue expression profile for a gene.

    Returns
    -------
    dict with keys: gene_symbol, tissue_tpm, max_vital_organ_tpm,
                    vital_organ_expressed, expression_category
    Or None if gene not found.
    """
    if expr_df is None:
        expr_df = load_expression_db()

    gene_data = expr_df[expr_df["gene_symbol"] == gene_symbol]
    if gene_data.empty:
        return None

    tissue_tpm = dict(zip(gene_data["tissue"], gene_data["tpm"]))
    max_vital = gene_data["max_vital_organ_tpm"].iloc[0]
    category = gene_data["expression_category"].iloc[0]

    vital_set = {t.lower() for t in VITAL_ORGANS}
    vital_expressed = any(
        tpm >= TPM_EXPRESSED_THRESHOLD
        for tissue, tpm in tissue_tpm.items()
        if tissue in vital_set
    )

    return {
        "gene_symbol": gene_symbol,
        "tissue_tpm": tissue_tpm,
        "max_vital_organ_tpm": float(max_vital),
        "vital_organ_expressed": vital_expressed,
        "expression_category": category,
    }
