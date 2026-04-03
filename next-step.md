# Decoy B 界面描述符升级 — Agent 执行指南

> **本文档是 agent 可直接执行的实施方案。**
> 每个 Task 包含：目标、涉及的精确文件和行号、要写/改的完整代码、验证命令。
> Agent 应按 Task 编号顺序执行，每完成一个 Task 运行对应的验证命令。

---

## 项目定位

代码库路径（服务器端，所有修改在此进行）:
```
/share/liuyutian/pMHC_decoy_library/
```

本地参考副本（只读，用于确认代码结构）:
```
C:\Users\30670\AppData\Local\Temp\pMHC_decoy_library\
```

## 当前状态

Decoy B 的结构评分只用了 3 个特征：
1. **Atchley 余弦相似度** — 序列层面 TCR 接触区理化性质
2. **肽段 Cα RMSD** — 骨架原子位置偏差（双超叠法：MHC→肽段 + 肽段→凹槽）
3. **B-factor 相关系数** — tFold 预测的温度因子

评分公式:
```
surface_correlation = max(avg(sim_a, sim_b), bf_corr)     # scanner.py:463-472
similarity_score    = 0.4 * cos_sim + 0.6 * surface_corr  # scanner.py:928-929
risk_score          = similarity × (1/EL_Rank) × TPM_Weight # risk_scorer.py:123
```

**问题**: 完全没有侧链非共价相互作用、界面埋藏面积、结合亲和力估计等关键特征。

## 升级目标

新增 3 个界面描述符（Phase 1）:
1. **PLIP 非共价相互作用指纹** — 氢键/疏水/盐桥/π-stacking 的 Tanimoto 相似度
2. **FreeSASA 界面埋藏面积 (BSA)** — peptide-MHC 界面的溶剂可及面积差
3. **PRODIGY 结合亲和力** — 基于界面接触类型的 ΔG 预测

---

## Task 0: 环境准备

### 0.1 安装 Python 依赖

```bash
cd /share/liuyutian/pMHC_decoy_library

# PLIP — 非共价相互作用分析
# PLIP 依赖 OpenBabel, 推荐 conda 安装
conda install -c conda-forge openbabel -y
pip install plip

# FreeSASA — 溶剂可及面积
pip install freesasa

# PRODIGY — 接触基亲和力预测
pip install prodigy-prot
```

### 0.2 验证安装

```bash
python -c "
from plip.structure.preparation import PDBComplex
print('[OK] PLIP')
"

python -c "
import freesasa
print('[OK] FreeSASA version:', freesasa.__version__ if hasattr(freesasa, '__version__') else 'installed')
"

python -c "
import subprocess
result = subprocess.run(['prodigy', '--help'], capture_output=True, text=True)
print('[OK] PRODIGY' if result.returncode == 0 else '[WARN] PRODIGY CLI not found, will use Python API')
"
```

### 0.3 验证现有 PDB 数据可用

```bash
# 确认有可用的 tFold 输出 PDB 用于测试
ls data/decoy_b/stepwise_demo/tfold_pdbs/*.pdb 2>/dev/null | head -5
# 或
ls data/decoy_b/*.pdb 2>/dev/null | head -5
```

如果没有 PDB 文件，需要先运行一次 Decoy B 生成结构:
```bash
python -m decoy_b scan-b --target GILGFVFTL --hla "HLA-A*02:01"
```

---

## Task 1: 创建界面描述符模块

### 目标
新建 `decoy_b/tools/interface_descriptors.py`，实现 PLIP / FreeSASA / PRODIGY 三个描述符的计算和相似度比较。

### 涉及文件
- **新建**: `decoy_b/tools/interface_descriptors.py`

### 实现

创建文件 `decoy_b/tools/interface_descriptors.py`，内容如下:

```python
"""
Interface descriptor computation for pMHC structures.

Computes three classes of descriptors from PDB structures:
1. PLIP: non-covalent interaction fingerprint (H-bonds, hydrophobic, salt bridges, π-stacking)
2. FreeSASA: buried surface area (BSA / ΔSASA) at the peptide-MHC interface
3. PRODIGY: contact-based binding affinity prediction (ΔG)

All pairwise similarity functions return normalized [0, 1] scores.

Usage:
    from decoy_b.tools.interface_descriptors import compute_interface_similarity
    sim = compute_interface_similarity(target_pdb, candidate_pdb)
    # sim.plip_tanimoto, sim.bsa_similarity, sim.prodigy_similarity, sim.combined
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

# tFold / AF3 output chain convention
DEFAULT_PEPTIDE_CHAIN = "P"
DEFAULT_MHC_CHAINS = ("M", "N")


# ═══════════════════════════════════════════════════════════════════════════
#  Data Structures
# ═══════════════════════════════════════════════════════════════════════════

class InteractionFingerprint:
    """Vectorized non-covalent interaction profile from PLIP.

    Stores both global interaction counts and a per-residue interaction map
    for computing Tanimoto similarity between two pMHC interfaces.
    """

    __slots__ = (
        "hbond_count", "hydrophobic_count", "salt_bridge_count",
        "pi_stacking_count", "pi_cation_count", "residue_interactions",
    )

    def __init__(
        self,
        hbond_count: int = 0,
        hydrophobic_count: int = 0,
        salt_bridge_count: int = 0,
        pi_stacking_count: int = 0,
        pi_cation_count: int = 0,
        residue_interactions: Optional[Dict[int, List[str]]] = None,
    ):
        self.hbond_count = hbond_count
        self.hydrophobic_count = hydrophobic_count
        self.salt_bridge_count = salt_bridge_count
        self.pi_stacking_count = pi_stacking_count
        self.pi_cation_count = pi_cation_count
        self.residue_interactions = residue_interactions or {}

    def to_vector(self) -> np.ndarray:
        """Convert to a 5-element count vector for Tanimoto computation."""
        return np.array([
            self.hbond_count,
            self.hydrophobic_count,
            self.salt_bridge_count,
            self.pi_stacking_count,
            self.pi_cation_count,
        ], dtype=float)

    def to_dict(self) -> dict:
        return {
            "hbond_count": self.hbond_count,
            "hydrophobic_count": self.hydrophobic_count,
            "salt_bridge_count": self.salt_bridge_count,
            "pi_stacking_count": self.pi_stacking_count,
            "pi_cation_count": self.pi_cation_count,
        }


class InterfaceDescriptors:
    """All interface descriptors computed for a single pMHC structure."""

    __slots__ = (
        "fingerprint", "bsa_total", "bsa_polar", "bsa_apolar",
        "prodigy_dg", "prodigy_kd",
    )

    def __init__(
        self,
        fingerprint: Optional[InteractionFingerprint] = None,
        bsa_total: Optional[float] = None,
        bsa_polar: Optional[float] = None,
        bsa_apolar: Optional[float] = None,
        prodigy_dg: Optional[float] = None,
        prodigy_kd: Optional[float] = None,
    ):
        self.fingerprint = fingerprint
        self.bsa_total = bsa_total
        self.bsa_polar = bsa_polar
        self.bsa_apolar = bsa_apolar
        self.prodigy_dg = prodigy_dg
        self.prodigy_kd = prodigy_kd


class InterfaceSimilarity:
    """Pairwise similarity between target and candidate interface descriptors."""

    __slots__ = (
        "plip_tanimoto", "bsa_similarity", "prodigy_similarity", "combined",
    )

    def __init__(
        self,
        plip_tanimoto: float = 0.0,
        bsa_similarity: float = 0.0,
        prodigy_similarity: float = 0.0,
        combined: float = 0.0,
    ):
        self.plip_tanimoto = plip_tanimoto
        self.bsa_similarity = bsa_similarity
        self.prodigy_similarity = prodigy_similarity
        self.combined = combined

    def to_dict(self) -> dict:
        return {
            "plip_tanimoto": self.plip_tanimoto,
            "bsa_similarity": self.bsa_similarity,
            "prodigy_similarity": self.prodigy_similarity,
            "combined": self.combined,
        }


# ═══════════════════════════════════════════════════════════════════════════
#  1. PLIP: Non-Covalent Interaction Fingerprint
# ═══════════════════════════════════════════════════════════════════════════

def compute_plip_fingerprint(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[InteractionFingerprint]:
    """
    Run PLIP on a pMHC PDB to extract non-covalent interactions
    at the peptide–MHC interface.

    PLIP identifies: hydrogen bonds, hydrophobic contacts, salt bridges,
    π-stacking, and cation-π interactions.

    Parameters
    ----------
    pdb_path : str
        Path to pMHC PDB file (tFold/AF3 output, chains M/N/P).
    peptide_chain : str
        Chain ID for peptide (default "P").
    mhc_chains : tuple of str
        Chain IDs for MHC heavy chain + β2m (default ("M", "N")).

    Returns
    -------
    InteractionFingerprint or None
        None if PLIP is not installed or analysis fails.

    Notes
    -----
    PLIP 的分析以 "binding site" 为单位。对于 pMHC，我们遍历所有 binding site
    中的相互作用，只保留跨越 peptide ↔ MHC 界面的那些。

    PLIP 需要 OpenBabel 支持，如果未安装会返回 None 并记录 warning。
    """
    try:
        from plip.structure.preparation import PDBComplex
    except ImportError:
        log.warning("PLIP not installed. Run: conda install -c conda-forge openbabel && pip install plip")
        return None

    try:
        mol = PDBComplex()
        mol.load_pdb(pdb_path)
        mol.analyze()

        hbond_count = 0
        hydrophobic_count = 0
        salt_bridge_count = 0
        pi_stacking_count = 0
        pi_cation_count = 0
        residue_interactions: Dict[int, List[str]] = {}

        mhc_set = set(mhc_chains)

        for bsid, interaction in mol.interaction_sets.items():
            # ── Hydrogen bonds ──
            for hb in interaction.hbonds_pdon + interaction.hbonds_ldon:
                if _crosses_interface(hb, peptide_chain, mhc_set):
                    hbond_count += 1
                    _record_peptide_residue(residue_interactions, hb, peptide_chain, "hbond")

            # ── Hydrophobic contacts ──
            for hc in interaction.hydrophobic_contacts:
                if _crosses_interface(hc, peptide_chain, mhc_set):
                    hydrophobic_count += 1
                    _record_peptide_residue(residue_interactions, hc, peptide_chain, "hydrophobic")

            # ── Salt bridges ──
            for sb in interaction.saltbridge_lneg + interaction.saltbridge_pneg:
                if _crosses_interface(sb, peptide_chain, mhc_set):
                    salt_bridge_count += 1
                    _record_peptide_residue(residue_interactions, sb, peptide_chain, "salt_bridge")

            # ── π-stacking ──
            for ps in interaction.pistacking:
                if _crosses_interface(ps, peptide_chain, mhc_set):
                    pi_stacking_count += 1
                    _record_peptide_residue(residue_interactions, ps, peptide_chain, "pi_stacking")

            # ── Cation-π ──
            for pc in interaction.pication_laro + interaction.pication_paro:
                if _crosses_interface(pc, peptide_chain, mhc_set):
                    pi_cation_count += 1
                    _record_peptide_residue(residue_interactions, pc, peptide_chain, "pi_cation")

        fp = InteractionFingerprint(
            hbond_count=hbond_count,
            hydrophobic_count=hydrophobic_count,
            salt_bridge_count=salt_bridge_count,
            pi_stacking_count=pi_stacking_count,
            pi_cation_count=pi_cation_count,
            residue_interactions=residue_interactions,
        )

        log.debug(
            "PLIP %s: %d hbonds, %d hydrophobic, %d salt bridges, %d π-stack, %d π-cation",
            Path(pdb_path).name,
            hbond_count, hydrophobic_count, salt_bridge_count,
            pi_stacking_count, pi_cation_count,
        )
        return fp

    except Exception as exc:
        log.warning("PLIP analysis failed for %s: %s", pdb_path, exc)
        return None


def compute_plip_tanimoto(
    fp_target: InteractionFingerprint,
    fp_candidate: InteractionFingerprint,
) -> float:
    """
    Generalized Tanimoto similarity between two interaction fingerprints.

    Uses count vectors [hbond, hydrophobic, salt_bridge, pi_stack, pi_cation]:
        T(A, B) = A·B / (|A|² + |B|² − A·B)

    Returns
    -------
    float in [0, 1]. Returns 1.0 if both vectors are zero (no interactions).
    """
    vec_t = fp_target.to_vector()
    vec_c = fp_candidate.to_vector()

    dot = float(np.dot(vec_t, vec_c))
    norm_t = float(np.dot(vec_t, vec_t))
    norm_c = float(np.dot(vec_c, vec_c))
    denom = norm_t + norm_c - dot

    if denom < 1e-12:
        # Both zero vectors → identical (both have no interactions)
        return 1.0 if (norm_t < 1e-12 and norm_c < 1e-12) else 0.0

    return max(0.0, min(1.0, dot / denom))


# ═══════════════════════════════════════════════════════════════════════════
#  2. FreeSASA: Buried Surface Area (BSA)
# ═══════════════════════════════════════════════════════════════════════════

def compute_bsa(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> Optional[float]:
    """
    Compute buried surface area (BSA) at the peptide–MHC interface.

    BSA = SASA(peptide alone) + SASA(MHC alone) − SASA(complex)

    This measures how much surface area is "buried" (hidden from solvent)
    when the peptide binds into the MHC groove. For 9-mer pMHC-I complexes,
    typical BSA values are 800–1500 Å².

    Implementation:
        FreeSASA 不直接支持 "计算单独链的 SASA"，所以我们通过写出
        临时 PDB 文件（只包含指定链）来分别计算。

    Parameters
    ----------
    pdb_path : str
        Path to pMHC PDB file.

    Returns
    -------
    float or None
        BSA in Å², or None if computation fails.
    """
    try:
        import freesasa
    except ImportError:
        log.warning("FreeSASA not installed. Run: pip install freesasa")
        return None

    try:
        # Parse PDB lines by chain
        pdb_lines = Path(pdb_path).read_text().splitlines()
        atom_lines = [l for l in pdb_lines if l.startswith(("ATOM", "HETATM"))]

        peptide_lines = [l for l in atom_lines if len(l) > 21 and l[21] == peptide_chain]
        mhc_lines = [l for l in atom_lines if len(l) > 21 and l[21] in set(mhc_chains)]

        if not peptide_lines or not mhc_lines:
            log.warning("Could not find peptide chain %s or MHC chains %s in %s",
                        peptide_chain, mhc_chains, pdb_path)
            return None

        def _sasa_from_lines(lines: List[str]) -> float:
            """Write lines to temp PDB, compute SASA, return total area."""
            with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
                f.write("\n".join(lines) + "\nEND\n")
                tmp_path = f.name
            try:
                structure = freesasa.Structure(tmp_path)
                result = freesasa.calc(structure)
                return result.totalArea()
            finally:
                Path(tmp_path).unlink(missing_ok=True)

        # SASA of the full complex
        complex_lines = peptide_lines + mhc_lines
        sasa_complex = _sasa_from_lines(complex_lines)

        # SASA of peptide alone
        sasa_peptide = _sasa_from_lines(peptide_lines)

        # SASA of MHC alone
        sasa_mhc = _sasa_from_lines(mhc_lines)

        # BSA = SASA(parts) − SASA(complex)
        bsa = sasa_peptide + sasa_mhc - sasa_complex

        log.debug("FreeSASA %s: peptide=%.1f, mhc=%.1f, complex=%.1f, BSA=%.1f Å²",
                  Path(pdb_path).name, sasa_peptide, sasa_mhc, sasa_complex, bsa)

        return max(0.0, bsa)

    except Exception as exc:
        log.warning("FreeSASA computation failed for %s: %s", pdb_path, exc)
        return None


def compute_bsa_similarity(
    bsa_target: float,
    bsa_candidate: float,
    max_bsa_diff: float = 800.0,
) -> float:
    """
    BSA similarity: 1.0 when identical, 0.0 when difference >= max_bsa_diff.

    For 9-mer pMHC, BSA typically ranges 800–1500 Å², so we set max_bsa_diff=800
    as the normalization range. Two pMHC complexes with identical BSA get
    score 1.0; those differing by 800+ Å² get 0.0.

    Parameters
    ----------
    bsa_target, bsa_candidate : float
        BSA values in Å².
    max_bsa_diff : float
        Maximum BSA difference for normalization (default 800 Å²).

    Returns
    -------
    float in [0, 1]
    """
    diff = abs(bsa_target - bsa_candidate)
    return max(0.0, 1.0 - diff / max_bsa_diff)


# ═══════════════════════════════════════════════════════════════════════════
#  3. PRODIGY: Contact-Based Binding Affinity
# ═══════════════════════════════════════════════════════════════════════════

def compute_prodigy_affinity(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chain: str = "M",
) -> Optional[float]:
    """
    Predict binding ΔG (kcal/mol) using PRODIGY.

    PRODIGY classifies interfacial contacts by residue polarity
    (charged/polar/apolar) and uses a linear model to predict ΔG.
    Runtime: ~3-5 seconds per structure.

    Parameters
    ----------
    pdb_path : str
        Path to pMHC PDB.
    peptide_chain : str
        Chain ID for peptide.
    mhc_chain : str
        Chain ID for MHC heavy chain (we compare peptide vs heavy chain).

    Returns
    -------
    float or None
        Predicted ΔG in kcal/mol (negative = favorable binding), or None.

    Notes
    -----
    PRODIGY 有两种调用方式:
    1. 命令行: `prodigy <pdb> --selection <chain1> <chain2>`
    2. Python API: `from prodigy.predict import predict_IC`

    我们先尝试 Python API，失败则回退到命令行。
    """
    # Attempt 1: Python API
    try:
        from prodigy.predict import predict_IC
        result = predict_IC(pdb_path, [peptide_chain], [mhc_chain])
        if result and hasattr(result, "dg_predicted"):
            dg = float(result.dg_predicted)
            log.debug("PRODIGY (API) %s: ΔG = %.2f kcal/mol", Path(pdb_path).name, dg)
            return dg
    except (ImportError, AttributeError):
        pass
    except Exception as exc:
        log.debug("PRODIGY Python API failed: %s, trying CLI", exc)

    # Attempt 2: CLI fallback
    try:
        result = subprocess.run(
            ["prodigy", pdb_path, "--selection", peptide_chain, mhc_chain],
            capture_output=True, text=True, timeout=30,
        )
        for line in result.stdout.splitlines():
            if "binding affinity" in line.lower() or "predicted" in line.lower():
                match = re.search(r"(-?\d+\.?\d*)\s*kcal", line)
                if match:
                    dg = float(match.group(1))
                    log.debug("PRODIGY (CLI) %s: ΔG = %.2f kcal/mol", Path(pdb_path).name, dg)
                    return dg
        # Try parsing the last numeric value as ΔG
        for line in reversed(result.stdout.splitlines()):
            match = re.search(r"(-?\d+\.?\d+)", line)
            if match:
                dg = float(match.group(1))
                if -30.0 < dg < 0.0:  # Sanity check: ΔG should be negative
                    log.debug("PRODIGY (CLI-fallback) %s: ΔG = %.2f kcal/mol", Path(pdb_path).name, dg)
                    return dg
    except FileNotFoundError:
        log.warning("PRODIGY CLI not found in PATH. Run: pip install prodigy-prot")
    except subprocess.TimeoutExpired:
        log.warning("PRODIGY timed out for %s", pdb_path)
    except Exception as exc:
        log.warning("PRODIGY failed for %s: %s", pdb_path, exc)

    return None


def compute_prodigy_similarity(
    dg_target: float,
    dg_candidate: float,
    max_dg_diff: float = 5.0,
) -> float:
    """
    ΔG similarity: candidates with similar predicted binding affinity
    to the target are more likely to compete for TCR binding.

    For pMHC-I, ΔG typically ranges −5 to −15 kcal/mol.
    We use max_dg_diff=5.0 kcal/mol as the normalization range.

    Returns
    -------
    float in [0, 1]
    """
    diff = abs(dg_target - dg_candidate)
    return max(0.0, 1.0 - diff / max_dg_diff)


# ═══════════════════════════════════════════════════════════════════════════
#  Composite: All Descriptors
# ═══════════════════════════════════════════════════════════════════════════

def compute_all_descriptors(
    pdb_path: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
) -> InterfaceDescriptors:
    """
    Compute all interface descriptors for a single pMHC structure.

    This is a convenience function that runs PLIP, FreeSASA, and PRODIGY
    in sequence and collects the results.

    Parameters
    ----------
    pdb_path : str
        Path to pMHC PDB file.

    Returns
    -------
    InterfaceDescriptors
    """
    fp = compute_plip_fingerprint(pdb_path, peptide_chain, mhc_chains)
    bsa = compute_bsa(pdb_path, peptide_chain, mhc_chains)
    dg = compute_prodigy_affinity(pdb_path, peptide_chain, mhc_chains[0])

    return InterfaceDescriptors(
        fingerprint=fp,
        bsa_total=bsa,
        prodigy_dg=dg,
    )


def compute_interface_similarity(
    target_pdb: str,
    candidate_pdb: str,
    peptide_chain: str = DEFAULT_PEPTIDE_CHAIN,
    mhc_chains: Tuple[str, ...] = DEFAULT_MHC_CHAINS,
    weights: Optional[Dict[str, float]] = None,
) -> InterfaceSimilarity:
    """
    Compute pairwise interface similarity between target and candidate pMHC.

    Runs all three descriptors on both structures, computes per-descriptor
    similarity, and returns a weighted combination.

    Parameters
    ----------
    target_pdb, candidate_pdb : str
        Paths to PDB files.
    peptide_chain : str
        Peptide chain ID (default "P").
    mhc_chains : tuple
        MHC chain IDs (default ("M", "N")).
    weights : dict, optional
        Weight for each descriptor in the combined score.
        Default: {"plip": 0.40, "bsa": 0.25, "prodigy": 0.35}

    Returns
    -------
    InterfaceSimilarity
        Per-descriptor and combined similarity scores, all in [0, 1].

    Notes
    -----
    权重设计考量:
    - PLIP (0.40): 非共价相互作用模式是 TCR 识别的最直接决定因素
    - PRODIGY (0.35): 结合亲和力反映热力学稳定性
    - BSA (0.25): 界面面积是必要条件但非充分条件

    如果某个描述符计算失败，其权重会被重新分配给成功的描述符。
    """
    if weights is None:
        weights = {"plip": 0.40, "bsa": 0.25, "prodigy": 0.35}

    desc_target = compute_all_descriptors(target_pdb, peptide_chain, mhc_chains)
    desc_candidate = compute_all_descriptors(candidate_pdb, peptide_chain, mhc_chains)

    # ── PLIP similarity ──
    plip_sim = 0.0
    plip_ok = False
    if desc_target.fingerprint is not None and desc_candidate.fingerprint is not None:
        plip_sim = compute_plip_tanimoto(desc_target.fingerprint, desc_candidate.fingerprint)
        plip_ok = True

    # ── BSA similarity ──
    bsa_sim = 0.0
    bsa_ok = False
    if desc_target.bsa_total is not None and desc_candidate.bsa_total is not None:
        bsa_sim = compute_bsa_similarity(desc_target.bsa_total, desc_candidate.bsa_total)
        bsa_ok = True

    # ── PRODIGY similarity ──
    prodigy_sim = 0.0
    prodigy_ok = False
    if desc_target.prodigy_dg is not None and desc_candidate.prodigy_dg is not None:
        prodigy_sim = compute_prodigy_similarity(desc_target.prodigy_dg, desc_candidate.prodigy_dg)
        prodigy_ok = True

    # ── Adaptive weighting: redistribute failed descriptor weights ──
    active = {}
    if plip_ok:
        active["plip"] = weights["plip"]
    if bsa_ok:
        active["bsa"] = weights["bsa"]
    if prodigy_ok:
        active["prodigy"] = weights["prodigy"]

    if not active:
        # All descriptors failed
        return InterfaceSimilarity(
            plip_tanimoto=plip_sim,
            bsa_similarity=bsa_sim,
            prodigy_similarity=prodigy_sim,
            combined=0.0,
        )

    # Normalize weights to sum to 1.0
    total_weight = sum(active.values())
    norm = {k: v / total_weight for k, v in active.items()}

    combined = (
        norm.get("plip", 0.0) * plip_sim
        + norm.get("bsa", 0.0) * bsa_sim
        + norm.get("prodigy", 0.0) * prodigy_sim
    )

    log.debug(
        "Interface similarity: PLIP=%.3f%s BSA=%.3f%s PRODIGY=%.3f%s → combined=%.3f",
        plip_sim, "" if plip_ok else "(N/A)",
        bsa_sim, "" if bsa_ok else "(N/A)",
        prodigy_sim, "" if prodigy_ok else "(N/A)",
        combined,
    )

    return InterfaceSimilarity(
        plip_tanimoto=plip_sim,
        bsa_similarity=bsa_sim,
        prodigy_similarity=prodigy_sim,
        combined=combined,
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Internal Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _crosses_interface(interaction, peptide_chain: str, mhc_chains: set) -> bool:
    """
    Check if a PLIP interaction crosses the peptide–MHC interface.

    A cross-interface interaction has one partner in the peptide chain
    and the other in an MHC chain. Intra-chain interactions are ignored.

    PLIP interaction objects vary by type, but generally have attributes like:
    - resnr / reschain (protein side)
    - resnr_l / reschain_l (ligand side)
    - or for symmetric types: reschain, reschain_l
    """
    try:
        chains = set()

        # Try various attribute patterns that PLIP uses
        for attr_chain in ("reschain", "reschain_l"):
            val = getattr(interaction, attr_chain, None)
            if val:
                chains.add(val)

        # Some interaction types use different naming
        for attr_chain in ("chain", "chain_l"):
            val = getattr(interaction, attr_chain, None)
            if val:
                chains.add(val)

        has_pep = peptide_chain in chains
        has_mhc = bool(chains & mhc_chains)
        return has_pep and has_mhc

    except Exception:
        return False


def _record_peptide_residue(
    residue_interactions: Dict[int, List[str]],
    interaction,
    peptide_chain: str,
    interaction_type: str,
) -> None:
    """Record which peptide residue is involved in an interaction."""
    try:
        # Find the residue number on the peptide chain side
        for attr_chain, attr_resnr in [
            ("reschain", "resnr"),
            ("reschain_l", "resnr_l"),
            ("chain", "resnr"),
        ]:
            chain_val = getattr(interaction, attr_chain, None)
            resnr_val = getattr(interaction, attr_resnr, None)
            if chain_val == peptide_chain and resnr_val is not None:
                resnum = int(resnr_val)
                if resnum not in residue_interactions:
                    residue_interactions[resnum] = []
                residue_interactions[resnum].append(interaction_type)
                return
    except Exception:
        pass
```

### 验证 Task 1

```bash
cd /share/liuyutian/pMHC_decoy_library

# 验证模块可导入
python -c "
from decoy_b.tools.interface_descriptors import (
    compute_plip_fingerprint,
    compute_bsa,
    compute_prodigy_affinity,
    compute_interface_similarity,
    InteractionFingerprint,
    InterfaceDescriptors,
    InterfaceSimilarity,
)
print('[OK] interface_descriptors module imports successfully')
"
```

---

## Task 2: 扩展 StructuralScore 数据模型

### 目标
在 `StructuralScore` 中新增字段存储界面描述符结果，保持向后兼容（所有新字段 Optional + 默认 None）。

### 涉及文件
- **修改**: `decoy_a/models.py` — `StructuralScore` 类 (L139-149)

### 精确修改

在 `decoy_a/models.py` 中，找到 `StructuralScore` 类（第 139-149 行），将其替换为:

```python
class StructuralScore(BaseModel):
    """3D structural and electrostatic similarity assessment."""
    modeling_tool: str = Field("", description="Tool used: biopython+dual_superposition | none | unavailable | error")
    pdb_path: Optional[str] = Field(None, description="Path to modeled pMHC PDB")
    surface_correlation: float = Field(
        0.0, description="Combined structural similarity score [0-1]"
    )
    rmsd: Optional[float] = Field(None, description="Peptide backbone RMSD to target pMHC (Å)")
    electrostatic_fingerprint_distance: Optional[float] = Field(
        None, description="Distance in electrostatic fingerprint space"
    )
    # ── Phase 1: Interface descriptors (new) ──
    plip_tanimoto: Optional[float] = Field(
        None, description="PLIP non-covalent interaction Tanimoto similarity [0-1]"
    )
    bsa_target: Optional[float] = Field(
        None, description="Target pMHC buried surface area (Å²)"
    )
    bsa_candidate: Optional[float] = Field(
        None, description="Candidate pMHC buried surface area (Å²)"
    )
    bsa_similarity: Optional[float] = Field(
        None, description="BSA similarity [0-1]"
    )
    prodigy_dg_target: Optional[float] = Field(
        None, description="Target PRODIGY ΔG prediction (kcal/mol)"
    )
    prodigy_dg_candidate: Optional[float] = Field(
        None, description="Candidate PRODIGY ΔG prediction (kcal/mol)"
    )
    prodigy_similarity: Optional[float] = Field(
        None, description="PRODIGY ΔG similarity [0-1]"
    )
    interface_combined: Optional[float] = Field(
        None, description="Weighted combination of all interface descriptors [0-1]"
    )
```

**定位方法**: 搜索 `class StructuralScore(BaseModel):` 替换整个类定义直到下一个 `class` 关键字。

**向后兼容**: 所有新字段都是 `Optional[float] = Field(None, ...)`, 不影响现有代码中 `StructuralScore(modeling_tool=..., surface_correlation=..., ...)` 的调用。

### 验证 Task 2

```bash
cd /share/liuyutian/pMHC_decoy_library

python -c "
from decoy_a.models import StructuralScore

# 旧式调用仍然工作
s1 = StructuralScore(modeling_tool='test', surface_correlation=0.5)
assert s1.plip_tanimoto is None  # 新字段默认 None
print('[OK] Backward compatible')

# 新式调用也工作
s2 = StructuralScore(
    modeling_tool='biopython+dual_superposition',
    surface_correlation=0.8,
    rmsd=1.2,
    plip_tanimoto=0.75,
    bsa_similarity=0.6,
    prodigy_similarity=0.8,
    interface_combined=0.72,
)
assert s2.plip_tanimoto == 0.75
print('[OK] New fields work')

# JSON 序列化正常
d = s2.model_dump(mode='json')
assert 'plip_tanimoto' in d
print('[OK] JSON serialization includes new fields')
print('[OK] StructuralScore model upgrade complete')
"
```

---

## Task 3: 集成界面描述符到结构比较函数

### 目标
在 `compute_structure_similarity()` 中，在现有 RMSD 计算之后调用界面描述符模块，并将结果融合到 `surface_correlation` 中。

### 涉及文件
- **修改**: `decoy_b/scanner.py` — `compute_structure_similarity()` 函数 (L336-498)

### 精确修改

在 `decoy_b/scanner.py` 的 `compute_structure_similarity()` 函数中，找到以下代码块（约第 448-486 行）:

```python
        # ── Combine scores ───────────────────────────────────────────────
        # Method A similarity: peptide RMSD → 0-1
        sim_a = 0.0
        rmsd_for_report = peptide_rmsd if peptide_rmsd is not None else mhc_fit_rmsd
        if peptide_rmsd is not None and peptide_rmsd < 3.0:
            sim_a = max(0.0, 1.0 - peptide_rmsd / 3.0)

        # Method B similarity: groove RMSD → 0-1
        sim_b = 0.0
        if groove_rmsd is not None and groove_rmsd < 3.0:
            sim_b = max(0.0, 1.0 - groove_rmsd / 3.0)

        # Final surface_correlation: average of both methods + B-factor floor
        # Both methods are highly correlated (r=0.95), but averaging captures
        # cases where one method flags a risk the other misses.
        if sim_a > 0 and sim_b > 0:
            combined = 0.5 * sim_a + 0.5 * sim_b
        elif sim_a > 0:
            combined = sim_a
        elif sim_b > 0:
            combined = sim_b
        else:
            combined = 0.0

        combined = max(combined, bf_corr)

        # Determine tool label
        n_tgt = len(target_pep_cas)
        n_cand = len(cand_pep_cas)
        if n_tgt == n_cand:
            tool = "biopython+dual_superposition"
        else:
            tool = f"biopython+dual_superposition({n_tgt}v{n_cand})"

        return StructuralScore(
            modeling_tool=tool,
            pdb_path=str(candidate_pdb),
            surface_correlation=max(0.0, combined),
            rmsd=rmsd_for_report,
        )
```

**替换为**:

```python
        # ── Combine RMSD scores ──────────────────────────────────────────
        # Method A similarity: peptide RMSD → 0-1
        sim_a = 0.0
        rmsd_for_report = peptide_rmsd if peptide_rmsd is not None else mhc_fit_rmsd
        if peptide_rmsd is not None and peptide_rmsd < 3.0:
            sim_a = max(0.0, 1.0 - peptide_rmsd / 3.0)

        # Method B similarity: groove RMSD → 0-1
        sim_b = 0.0
        if groove_rmsd is not None and groove_rmsd < 3.0:
            sim_b = max(0.0, 1.0 - groove_rmsd / 3.0)

        # RMSD geometric similarity (average of both methods)
        if sim_a > 0 and sim_b > 0:
            rmsd_geo = 0.5 * sim_a + 0.5 * sim_b
        elif sim_a > 0:
            rmsd_geo = sim_a
        elif sim_b > 0:
            rmsd_geo = sim_b
        else:
            rmsd_geo = 0.0
        rmsd_geo = max(rmsd_geo, bf_corr)

        # ── Interface descriptors (PLIP + FreeSASA + PRODIGY) ────────────
        interface_sim = None
        plip_tanimoto = None
        bsa_similarity = None
        bsa_target_val = None
        bsa_candidate_val = None
        prodigy_similarity = None
        prodigy_dg_target = None
        prodigy_dg_candidate = None
        interface_combined = None

        try:
            from decoy_b.tools.interface_descriptors import (
                compute_all_descriptors,
                compute_plip_tanimoto,
                compute_bsa_similarity,
                compute_prodigy_similarity as _prodigy_sim,
            )

            desc_tgt = compute_all_descriptors(
                str(target_pdb), PEPTIDE_CHAIN, tuple(MHC_CHAINS),
            )
            desc_cand = compute_all_descriptors(
                str(candidate_pdb), PEPTIDE_CHAIN, tuple(MHC_CHAINS),
            )

            # PLIP
            if desc_tgt.fingerprint is not None and desc_cand.fingerprint is not None:
                plip_tanimoto = compute_plip_tanimoto(
                    desc_tgt.fingerprint, desc_cand.fingerprint,
                )

            # BSA
            bsa_target_val = desc_tgt.bsa_total
            bsa_candidate_val = desc_cand.bsa_total
            if bsa_target_val is not None and bsa_candidate_val is not None:
                bsa_similarity = compute_bsa_similarity(bsa_target_val, bsa_candidate_val)

            # PRODIGY
            prodigy_dg_target = desc_tgt.prodigy_dg
            prodigy_dg_candidate = desc_cand.prodigy_dg
            if prodigy_dg_target is not None and prodigy_dg_candidate is not None:
                prodigy_similarity = _prodigy_sim(prodigy_dg_target, prodigy_dg_candidate)

            # Adaptive weighted combination of interface descriptors
            active_scores = {}
            if plip_tanimoto is not None:
                active_scores["plip"] = (plip_tanimoto, 0.40)
            if bsa_similarity is not None:
                active_scores["bsa"] = (bsa_similarity, 0.25)
            if prodigy_similarity is not None:
                active_scores["prodigy"] = (prodigy_similarity, 0.35)

            if active_scores:
                total_w = sum(w for _, w in active_scores.values())
                interface_combined = sum(
                    s * (w / total_w) for s, w in active_scores.values()
                )

        except ImportError:
            log.debug("Interface descriptors not available (missing PLIP/FreeSASA/PRODIGY)")
        except Exception as exc:
            log.warning("Interface descriptor computation failed: %s", exc)

        # ── Final combined score ─────────────────────────────────────────
        # If interface descriptors available: blend RMSD + interface
        # Weight rationale:
        #   RMSD geometric (0.50) — backbone conformational agreement
        #   Interface descriptors (0.50) — side-chain non-covalent interactions
        if interface_combined is not None and interface_combined > 0:
            combined = 0.50 * rmsd_geo + 0.50 * interface_combined
        else:
            combined = rmsd_geo

        # Determine tool label
        n_tgt = len(target_pep_cas)
        n_cand = len(cand_pep_cas)
        if n_tgt == n_cand:
            tool = "biopython+dual_superposition"
        else:
            tool = f"biopython+dual_superposition({n_tgt}v{n_cand})"
        if interface_combined is not None:
            tool += "+interface"

        return StructuralScore(
            modeling_tool=tool,
            pdb_path=str(candidate_pdb),
            surface_correlation=max(0.0, combined),
            rmsd=rmsd_for_report,
            plip_tanimoto=plip_tanimoto,
            bsa_target=bsa_target_val,
            bsa_candidate=bsa_candidate_val,
            bsa_similarity=bsa_similarity,
            prodigy_dg_target=prodigy_dg_target,
            prodigy_dg_candidate=prodigy_dg_candidate,
            prodigy_similarity=prodigy_similarity,
            interface_combined=interface_combined,
        )
```

### 不需要修改的部分

- `scan_decoy_b()` 中第 928-929 行的 `combined = 0.4 * cos_sim + 0.6 * structural.surface_correlation` **不需要改**，因为 `surface_correlation` 现在已经包含了界面描述符的贡献。
- `risk_scorer.py` **不需要改**，它读取 `similarity_score` 和 `surface_correlation`，两者都已经升级。

### 验证 Task 3

```bash
cd /share/liuyutian/pMHC_decoy_library

# 找两个已有的 PDB 测试
python -c "
from pathlib import Path
from decoy_b.scanner import compute_structure_similarity

# 寻找可用的 PDB 文件
pdb_dir = Path('data/decoy_b')
pdbs = sorted(pdb_dir.rglob('*.pdb'))
if len(pdbs) < 2:
    print('[SKIP] Need at least 2 PDB files. Run decoy_b first.')
else:
    target = str(pdbs[0])
    candidate = str(pdbs[1])
    print(f'Target: {target}')
    print(f'Candidate: {candidate}')

    result = compute_structure_similarity(target, candidate)
    print(f'surface_correlation: {result.surface_correlation:.4f}')
    print(f'rmsd: {result.rmsd}')
    print(f'plip_tanimoto: {result.plip_tanimoto}')
    print(f'bsa_similarity: {result.bsa_similarity}')
    print(f'prodigy_similarity: {result.prodigy_similarity}')
    print(f'interface_combined: {result.interface_combined}')
    print(f'modeling_tool: {result.modeling_tool}')
    print('[OK] Structure comparison with interface descriptors works')
"
```

---

## Task 4: 更新 DecoyABEntry 输出格式

### 目标
在最终输出 `DecoyABEntry` 中暴露界面描述符的明细数据，方便下游分析和可视化。

### 涉及文件
- **修改**: `decoy_a/models.py` — `DecoyABEntry` 类 (L188-233)

### 精确修改

在 `DecoyABEntry` 类中，找到以下行（约 L211）:

```python
    structural_similarity: float = Field(0.0)
```

在其后添加:

```python
    # Interface descriptor details (from Phase 1 upgrade)
    plip_tanimoto: Optional[float] = Field(None, description="PLIP interaction Tanimoto [0-1]")
    bsa_similarity: Optional[float] = Field(None, description="BSA similarity [0-1]")
    prodigy_similarity: Optional[float] = Field(None, description="PRODIGY ΔG similarity [0-1]")
    interface_combined: Optional[float] = Field(None, description="Interface descriptor combined [0-1]")
```

然后修改 `risk_scorer.py` 中 Decoy B hit 处理，使其传递新字段。

在 `risk_scorer.py` 第 211-231 行的 `DecoyABEntry` 构造中，找到:

```python
        entry = DecoyABEntry(
            decoy_ab_id=f"DAB-{counter:04d}",
            sequence=hit.sequence,
            target_sequence=hit.target_sequence,
            hla_allele=hit.hla_allele,
            source=DecoySource.DECOY_B,
            gene_symbols=hit.gene_symbols,
            source_proteins=hit.source_proteins,
            el_rank=hit.el_rank,
            hamming_distance=hit.hamming_distance,
            physicochemical_similarity=hit.physicochemical.cosine_similarity,
            structural_similarity=(
                hit.structural.surface_correlation if hit.structural else 0.0
            ),
            expression=hit.expression,
            vital_organ_tpm_weight=tpm_weight,
            critical_organs=critical_organs,
            total_risk_score=risk,
            structural=hit.structural,
        )
```

**替换为**:

```python
        entry = DecoyABEntry(
            decoy_ab_id=f"DAB-{counter:04d}",
            sequence=hit.sequence,
            target_sequence=hit.target_sequence,
            hla_allele=hit.hla_allele,
            source=DecoySource.DECOY_B,
            gene_symbols=hit.gene_symbols,
            source_proteins=hit.source_proteins,
            el_rank=hit.el_rank,
            hamming_distance=hit.hamming_distance,
            physicochemical_similarity=hit.physicochemical.cosine_similarity,
            structural_similarity=(
                hit.structural.surface_correlation if hit.structural else 0.0
            ),
            plip_tanimoto=(
                hit.structural.plip_tanimoto if hit.structural else None
            ),
            bsa_similarity=(
                hit.structural.bsa_similarity if hit.structural else None
            ),
            prodigy_similarity=(
                hit.structural.prodigy_similarity if hit.structural else None
            ),
            interface_combined=(
                hit.structural.interface_combined if hit.structural else None
            ),
            expression=hit.expression,
            vital_organ_tpm_weight=tpm_weight,
            critical_organs=critical_organs,
            total_risk_score=risk,
            structural=hit.structural,
        )
```

### 验证 Task 4

```bash
cd /share/liuyutian/pMHC_decoy_library

python -c "
from decoy_a.models import DecoyABEntry, StructuralScore

# 模拟构造一个带界面描述符的 entry
entry = DecoyABEntry(
    decoy_ab_id='DAB-0001',
    sequence='GILGFVFTL',
    target_sequence='GILGFVFTL',
    hla_allele='HLA-A*02:01',
    source='Decoy_B_Structural_Similarity',
    el_rank=0.5,
    hamming_distance=0,
    structural_similarity=0.85,
    plip_tanimoto=0.75,
    bsa_similarity=0.6,
    prodigy_similarity=0.8,
    interface_combined=0.72,
    total_risk_score=1.5,
)

d = entry.model_dump(mode='json')
assert d['plip_tanimoto'] == 0.75
assert d['interface_combined'] == 0.72
print('[OK] DecoyABEntry with interface descriptors works')
print('JSON keys:', [k for k in d if k.startswith(('plip', 'bsa_s', 'prodigy_s', 'interface'))])
"
```

---

## Task 5: 单元测试

### 目标
编写测试覆盖所有新增功能：PLIP 指纹提取、BSA 计算、PRODIGY 亲和力、Tanimoto 相似度、端到端集成。

### 涉及文件
- **新建**: `tests/test_interface_descriptors.py`

### 实现

创建文件 `tests/test_interface_descriptors.py`:

```python
"""
Tests for the interface descriptor module (Phase 1).

Uses existing tFold PDB outputs in data/decoy_b/ for integration tests.
Unit tests for similarity functions use synthetic data.
"""

import pytest
import numpy as np
from pathlib import Path
import sys

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from decoy_b.tools.interface_descriptors import (
    InteractionFingerprint,
    InterfaceSimilarity,
    compute_plip_fingerprint,
    compute_plip_tanimoto,
    compute_bsa,
    compute_bsa_similarity,
    compute_prodigy_affinity,
    compute_prodigy_similarity,
    compute_interface_similarity,
    compute_all_descriptors,
)


# ── Paths ──────────────────────────────────────────────────────────────

DATA_DIR = PROJECT_ROOT / "data" / "decoy_b"
STEPWISE_DIR = DATA_DIR / "stepwise_demo" / "tfold_pdbs"


def _find_pdbs() -> list:
    """Find available PDB files for testing."""
    pdbs = []
    for d in [STEPWISE_DIR, DATA_DIR]:
        pdbs.extend(sorted(d.rglob("*.pdb")))
    return pdbs


# ── Unit Tests (no PDB needed) ────────────────────────────────────────

class TestTanimoto:
    """Test Tanimoto similarity computation on synthetic fingerprints."""

    def test_identical_fingerprints(self):
        fp = InteractionFingerprint(
            hbond_count=5, hydrophobic_count=10, salt_bridge_count=2,
        )
        assert compute_plip_tanimoto(fp, fp) == pytest.approx(1.0)

    def test_zero_fingerprints(self):
        fp1 = InteractionFingerprint()  # all zeros
        fp2 = InteractionFingerprint()
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(1.0)

    def test_disjoint_fingerprints(self):
        fp1 = InteractionFingerprint(hbond_count=10)
        fp2 = InteractionFingerprint(hydrophobic_count=10)
        # Tanimoto of orthogonal vectors = 0
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(0.0)

    def test_similar_fingerprints(self):
        fp1 = InteractionFingerprint(hbond_count=5, hydrophobic_count=10)
        fp2 = InteractionFingerprint(hbond_count=4, hydrophobic_count=9)
        sim = compute_plip_tanimoto(fp1, fp2)
        assert 0.5 < sim < 1.0

    def test_symmetry(self):
        fp1 = InteractionFingerprint(hbond_count=3, salt_bridge_count=1)
        fp2 = InteractionFingerprint(hbond_count=5, hydrophobic_count=2)
        assert compute_plip_tanimoto(fp1, fp2) == pytest.approx(
            compute_plip_tanimoto(fp2, fp1)
        )


class TestBSASimilarity:
    """Test BSA similarity computation."""

    def test_identical(self):
        assert compute_bsa_similarity(1000.0, 1000.0) == pytest.approx(1.0)

    def test_max_diff(self):
        assert compute_bsa_similarity(1000.0, 1800.0) == pytest.approx(0.0)

    def test_mid_range(self):
        sim = compute_bsa_similarity(1000.0, 1400.0)
        assert sim == pytest.approx(0.5)


class TestProdigySimilarity:
    """Test PRODIGY ΔG similarity computation."""

    def test_identical(self):
        assert compute_prodigy_similarity(-10.0, -10.0) == pytest.approx(1.0)

    def test_max_diff(self):
        assert compute_prodigy_similarity(-10.0, -15.0) == pytest.approx(0.0)

    def test_mid_range(self):
        sim = compute_prodigy_similarity(-10.0, -12.5)
        assert sim == pytest.approx(0.5)


# ── Integration Tests (need PDB files) ────────────────────────────────

class TestPLIPIntegration:
    """Test PLIP fingerprint extraction on real tFold PDBs."""

    def test_fingerprint_extraction(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        fp = compute_plip_fingerprint(str(pdbs[0]))
        if fp is None:
            pytest.skip("PLIP not available (missing OpenBabel?)")

        assert fp.hbond_count >= 0
        assert fp.hydrophobic_count >= 0
        total = fp.hbond_count + fp.hydrophobic_count + fp.salt_bridge_count
        assert total > 0, "Expected at least some interactions in pMHC"

    def test_self_tanimoto_is_one(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        fp = compute_plip_fingerprint(str(pdbs[0]))
        if fp is None:
            pytest.skip("PLIP not available")

        assert compute_plip_tanimoto(fp, fp) == pytest.approx(1.0)


class TestFreeSASAIntegration:
    """Test FreeSASA BSA computation on real tFold PDBs."""

    def test_bsa_positive(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        bsa = compute_bsa(str(pdbs[0]))
        if bsa is None:
            pytest.skip("FreeSASA not available")

        # pMHC complex should have positive BSA (some surface is buried)
        assert bsa > 0, f"Expected positive BSA, got {bsa}"
        # Typical 9-mer pMHC BSA: 800-1500 Å²
        assert bsa < 3000, f"BSA unreasonably large: {bsa}"


class TestPRODIGYIntegration:
    """Test PRODIGY affinity prediction on real tFold PDBs."""

    def test_prodigy_returns_negative_dg(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        dg = compute_prodigy_affinity(str(pdbs[0]))
        if dg is None:
            pytest.skip("PRODIGY not available")

        # ΔG should be negative for a bound complex
        assert dg < 0, f"Expected negative ΔG, got {dg}"


class TestEndToEnd:
    """End-to-end test: full interface similarity between two pMHCs."""

    def test_interface_similarity(self):
        pdbs = _find_pdbs()
        if len(pdbs) < 2:
            pytest.skip("Need at least 2 PDB files")

        sim = compute_interface_similarity(str(pdbs[0]), str(pdbs[1]))

        assert 0.0 <= sim.combined <= 1.0
        assert 0.0 <= sim.plip_tanimoto <= 1.0
        assert 0.0 <= sim.bsa_similarity <= 1.0
        assert 0.0 <= sim.prodigy_similarity <= 1.0

    def test_self_similarity_is_one(self):
        pdbs = _find_pdbs()
        if not pdbs:
            pytest.skip("No PDB files available")

        sim = compute_interface_similarity(str(pdbs[0]), str(pdbs[0]))

        # Self-similarity should be 1.0 (or very close)
        if sim.plip_tanimoto is not None:
            assert sim.plip_tanimoto == pytest.approx(1.0)
        if sim.bsa_similarity is not None:
            assert sim.bsa_similarity == pytest.approx(1.0)
        if sim.prodigy_similarity is not None:
            assert sim.prodigy_similarity == pytest.approx(1.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
```

### 验证 Task 5

```bash
cd /share/liuyutian/pMHC_decoy_library

# 运行所有单元测试（不需要 PDB）
pytest tests/test_interface_descriptors.py -v -k "not Integration and not EndToEnd"

# 运行集成测试（需要 PDB 文件）
pytest tests/test_interface_descriptors.py -v -k "Integration or EndToEnd"

# 全部运行
pytest tests/test_interface_descriptors.py -v
```

---

## Task 6: 全管线回归测试

### 目标
在 GILGFVFTL 靶标上运行完整 Decoy B 管线，对比升级前后的排名变化，确认无回归。

### 步骤

```bash
cd /share/liuyutian/pMHC_decoy_library

# 1. 备份现有结果
cp data/decoy_b/final_ranked_decoys.json data/decoy_b/final_ranked_decoys_BACKUP.json 2>/dev/null

# 2. 运行完整管线
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"

# 3. 查看结果
python -m decoy_b show --format table

# 4. 对比新旧 Top 20
python -c "
import json
from pathlib import Path

backup = Path('data/decoy_b/final_ranked_decoys_BACKUP.json')
current = Path('data/decoy_b/final_ranked_decoys.json')

if backup.exists() and current.exists():
    old = json.loads(backup.read_text())[:20]
    new = json.loads(current.read_text())[:20]

    old_seqs = [e['sequence'] for e in old]
    new_seqs = [e['sequence'] for e in new]

    overlap = set(old_seqs) & set(new_seqs)
    print(f'Top 20 overlap: {len(overlap)}/20')
    print(f'Old Top 5: {old_seqs[:5]}')
    print(f'New Top 5: {new_seqs[:5]}')

    # Show new interface descriptor fields
    if new:
        e = new[0]
        print(f'\\nTop 1 details:')
        print(f'  sequence: {e[\"sequence\"]}')
        print(f'  total_risk_score: {e[\"total_risk_score\"]:.4f}')
        print(f'  structural_similarity: {e.get(\"structural_similarity\", \"N/A\")}')
        print(f'  plip_tanimoto: {e.get(\"plip_tanimoto\", \"N/A\")}')
        print(f'  bsa_similarity: {e.get(\"bsa_similarity\", \"N/A\")}')
        print(f'  prodigy_similarity: {e.get(\"prodigy_similarity\", \"N/A\")}')
        print(f'  interface_combined: {e.get(\"interface_combined\", \"N/A\")}')
else:
    print('No backup to compare against, or no current results.')
"
```

### 预期结果
- Top 20 中应有 >50% 与旧结果重叠（核心高危序列不应消失）
- 新字段 `plip_tanimoto` / `bsa_similarity` / `prodigy_similarity` / `interface_combined` 应出现在输出中
- `modeling_tool` 字段应包含 `+interface` 后缀

---

## Task 7: 可视化更新（可选）

### 目标
在现有可视化脚本中添加界面描述符的展示。

### 涉及文件
- **修改**: `scripts/visualize_GILGFVFTL_summary.py`（如果存在）
- **或新建**: `scripts/visualize_interface_descriptors.py`

### 实现思路

```python
"""
Visualize interface descriptor distributions for Decoy B results.

Generates:
1. Radar chart: per-candidate descriptor profile (PLIP, BSA, PRODIGY, RMSD, Atchley)
2. Scatter plot: interface_combined vs structural_similarity (with/without interface)
3. Heatmap: per-residue interaction types from PLIP
"""
# 具体实现根据现有可视化风格进行
# 加载 final_ranked_decoys.json，提取新字段，用 matplotlib/plotly 绑制
```

---

## 依赖总结

| 依赖 | 用途 | 安装命令 | 是否必须 |
|------|------|---------|---------|
| `plip` | 非共价相互作用分析 | `pip install plip` | 是 |
| `openbabel` | PLIP 的底层依赖 | `conda install -c conda-forge openbabel` | 是 |
| `freesasa` | 溶剂可及面积计算 | `pip install freesasa` | 是 |
| `prodigy-prot` | 结合亲和力预测 | `pip install prodigy-prot` | 推荐 |

---

## 评分公式变化总结

### 旧公式
```
surface_correlation = max(avg(sim_a, sim_b), bf_corr)
  其中 sim_a = 1 - peptide_rmsd / 3.0   (MHC 超叠后肽段 RMSD)
       sim_b = 1 - groove_rmsd / 3.0     (肽段超叠后凹槽 RMSD)

similarity_score = 0.4 × cos_sim + 0.6 × surface_correlation
risk = similarity_score × (1/EL_Rank) × TPM_Weight
```

### 新公式
```
rmsd_geo = max(avg(sim_a, sim_b), bf_corr)       # 与旧公式相同

interface_combined = weighted_avg(                  # 新增
    0.40 × plip_tanimoto,                          # 非共价相互作用
    0.25 × bsa_similarity,                         # 界面埋藏面积
    0.35 × prodigy_similarity,                     # 结合亲和力
)

surface_correlation = 0.50 × rmsd_geo + 0.50 × interface_combined  # 融合

similarity_score = 0.4 × cos_sim + 0.6 × surface_correlation      # 不变
risk = similarity_score × (1/EL_Rank) × TPM_Weight                # 不变
```

**核心变化**: `surface_correlation` 从纯 RMSD 扩展为 RMSD + 界面描述符的融合。上下游公式不变，最小化改动范围。

---

## 执行顺序

```
Task 0: 环境准备 (pip install)
  ↓
Task 1: 创建 decoy_b/tools/interface_descriptors.py
  ↓
Task 2: 扩展 StructuralScore 数据模型
  ↓
Task 3: 集成到 compute_structure_similarity()
  ↓
Task 4: 更新 DecoyABEntry 输出格式 + risk_scorer.py
  ↓
Task 5: 单元测试
  ↓
Task 6: 全管线回归测试
  ↓
Task 7: 可视化更新（可选）
```

每个 Task 完成后运行对应的验证命令，全部通过再进入下一个 Task。
