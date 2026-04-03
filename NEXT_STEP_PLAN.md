# Next-Step Plan: Decoy B 界面描述符升级

## 背景与目标

当前 Decoy B 的评分仅使用 Atchley 余弦相似度 (40%) + 骨架 RMSD 相似度 (60%)，缺乏对**侧链非共价相互作用、静电互补、界面埋藏面积、相互作用指纹**等关键特征的捕捉。Deep Research 报告明确指出：序列相似性和 Cα RMSD 几乎无法预测交叉反应性，真正的分子模拟发生在侧链理化景观层面。

**目标**: 分三个阶段，将界面描述符从 2 维扩展到 ~10 维，最终用 LTR 模型替代线性加权。

---

## 代码架构概览（集成点）

```
关键文件与行号（基于服务器端代码 /share/liuyutian 下的 pMHC_decoy_library/）：

decoy_a/models.py
  ├── StructuralScore        (L139-149)  ← 需要扩展字段
  ├── DecoyBHit              (L152-183)  ← similarity_score 计算逻辑
  └── DecoyABEntry           (L188-233)  ← 最终输出，需要增加新特征列

decoy_b/scanner.py
  ├── compute_structure_similarity()  (L336-499)  ← 核心集成点，在这里添加新描述符
  ├── scan_decoy_b()                  (L746-979)  ← Stage 4 hit 组装，combined score 在 L928
  └── run_mpnn_design()               (L505-739)  ← MPNN 分支，同样需要更新评分

decoy_b/risk_scorer.py
  ├── compute_risk_score()    (L96-123)   ← 最终风险公式
  └── score_and_rank()        (L128-262)  ← A+B 合并排序

decoy_a/config.py
  └── 阈值与常量             (L131-166)  ← 新描述符的阈值
```

**当前评分数据流**:
```
Atchley cosine (cos_sim) ─┐
                          ├─→ combined = 0.4 * cos_sim + 0.6 * surface_correlation
surface_correlation ──────┘       ↓
                          risk = combined × (1/EL_Rank) × TPM_Weight
```

---

## Phase 1: 快速收益（PLIP + FreeSASA + PRODIGY）

> 目标：从纯骨架 RMSD 升级到包含非共价相互作用和界面埋藏面积的多维评分
> 预计工期：1-2 周
> 计算开销：<5s/结构，不影响现有吞吐量

### Step 1.1: 安装依赖

```bash
# 在服务器上执行
pip install plip freesasa prodigy-prot

# 验证安装
python -c "from plip.structure.preparation import PDBComplex; print('PLIP OK')"
python -c "import freesasa; print('FreeSASA OK')"
python -c "from prodigy import predict_IC; print('PRODIGY OK')"
```

**注意**: 
- PLIP 依赖 OpenBabel，可能需要 `conda install -c conda-forge openbabel`
- FreeSASA 是纯 C 扩展，pip 直接安装
- PRODIGY 如果 pip 版本有问题，可以从 GitHub 安装: `pip install git+https://github.com/haddocking/prodigy`

### Step 1.2: 创建界面描述符模块

**新建文件**: `decoy_b/tools/interface_descriptors.py`

```python
"""
Interface descriptor computation for pMHC structures.

Computes:
1. PLIP interaction fingerprint (H-bonds, hydrophobic, salt bridges, π-stacking)
2. FreeSASA buried surface area (BSA/ΔSASA)
3. PRODIGY contact-based affinity prediction

All functions accept two PDB paths (target, candidate) and return
normalized [0, 1] similarity scores.
"""

from __future__ import annotations
import logging
import tempfile
from pathlib import Path
from typing import Optional, NamedTuple
import numpy as np

log = logging.getLogger(__name__)


# ─── Data Structures ──────────────────────────────────────────────────────

class InteractionFingerprint(NamedTuple):
    """Vectorized non-covalent interaction profile from PLIP."""
    hbond_count: int
    hydrophobic_count: int
    salt_bridge_count: int
    pi_stacking_count: int
    pi_cation_count: int
    # Per-residue interaction map: {residue_position: [interaction_types]}
    residue_interactions: dict

class InterfaceDescriptors(NamedTuple):
    """All interface descriptors for a single pMHC structure."""
    # PLIP
    interaction_fingerprint: Optional[InteractionFingerprint]
    # FreeSASA
    bsa_total: Optional[float]          # Total buried surface area (Å²)
    bsa_polar: Optional[float]          # Polar BSA (Å²)
    bsa_apolar: Optional[float]         # Apolar BSA (Å²)
    # PRODIGY
    prodigy_dg: Optional[float]         # Predicted ΔG (kcal/mol)
    prodigy_kd: Optional[float]         # Predicted Kd (M)

class InterfaceSimilarity(NamedTuple):
    """Pairwise similarity between target and candidate interface descriptors."""
    plip_tanimoto: float                # Tanimoto coefficient on interaction fingerprints [0-1]
    bsa_similarity: float               # BSA similarity [0-1]
    prodigy_similarity: float           # ΔG similarity [0-1]
    combined: float                     # Weighted combination [0-1]


# ─── PLIP: Non-Covalent Interaction Fingerprint ───────────────────────────

def compute_plip_fingerprint(
    pdb_path: str,
    peptide_chain: str = "P",
    mhc_chains: tuple = ("M", "N"),
) -> Optional[InteractionFingerprint]:
    """
    Run PLIP on a pMHC PDB to extract non-covalent interactions
    at the peptide-MHC interface.
    
    Args:
        pdb_path: Path to pMHC PDB file (tFold/AF3 output)
        peptide_chain: Chain ID for peptide (default "P" for tFold)
        mhc_chains: Chain IDs for MHC heavy chain + β2m
    
    Returns:
        InteractionFingerprint or None if PLIP fails
    """
    try:
        from plip.structure.preparation import PDBComplex
        
        mol = PDBComplex()
        mol.load_pdb(pdb_path)
        mol.analyze()
        
        hbond_count = 0
        hydrophobic_count = 0
        salt_bridge_count = 0
        pi_stacking_count = 0
        pi_cation_count = 0
        residue_interactions = {}
        
        for bsid, interaction in mol.interaction_sets.items():
            # Filter: only peptide ↔ MHC interactions
            # PLIP returns interactions indexed by binding site
            
            # Hydrogen bonds
            for hb in interaction.hbonds_pdon + interaction.hbonds_ldon:
                if _is_peptide_mhc_interaction(hb, peptide_chain, mhc_chains):
                    hbond_count += 1
                    _record_residue(residue_interactions, hb, "hbond")
            
            # Hydrophobic contacts
            for hc in interaction.hydrophobic_contacts:
                if _is_peptide_mhc_interaction(hc, peptide_chain, mhc_chains):
                    hydrophobic_count += 1
                    _record_residue(residue_interactions, hc, "hydrophobic")
            
            # Salt bridges
            for sb in interaction.saltbridge_lneg + interaction.saltbridge_pneg:
                if _is_peptide_mhc_interaction(sb, peptide_chain, mhc_chains):
                    salt_bridge_count += 1
                    _record_residue(residue_interactions, sb, "salt_bridge")
            
            # π-stacking
            for ps in interaction.pistacking:
                if _is_peptide_mhc_interaction(ps, peptide_chain, mhc_chains):
                    pi_stacking_count += 1
                    _record_residue(residue_interactions, ps, "pi_stacking")
            
            # π-cation
            for pc in interaction.pication_laro + interaction.pication_paro:
                if _is_peptide_mhc_interaction(pc, peptide_chain, mhc_chains):
                    pi_cation_count += 1
                    _record_residue(residue_interactions, pc, "pi_cation")
        
        return InteractionFingerprint(
            hbond_count=hbond_count,
            hydrophobic_count=hydrophobic_count,
            salt_bridge_count=salt_bridge_count,
            pi_stacking_count=pi_stacking_count,
            pi_cation_count=pi_cation_count,
            residue_interactions=residue_interactions,
        )
    except Exception as exc:
        log.warning("PLIP analysis failed for %s: %s", pdb_path, exc)
        return None


def compute_plip_tanimoto(
    fp_target: InteractionFingerprint,
    fp_candidate: InteractionFingerprint,
) -> float:
    """
    Tanimoto similarity between two PLIP interaction fingerprints.
    
    Uses both global counts (interaction type totals) and per-residue 
    interaction patterns for TCR-contact positions (p4-p8).
    
    Returns: float in [0, 1]
    """
    # Convert to count vectors: [hbond, hydrophobic, salt_bridge, pi_stack, pi_cation]
    vec_t = np.array([
        fp_target.hbond_count,
        fp_target.hydrophobic_count,
        fp_target.salt_bridge_count,
        fp_target.pi_stacking_count,
        fp_target.pi_cation_count,
    ], dtype=float)
    
    vec_c = np.array([
        fp_candidate.hbond_count,
        fp_candidate.hydrophobic_count,
        fp_candidate.salt_bridge_count,
        fp_candidate.pi_stacking_count,
        fp_candidate.pi_cation_count,
    ], dtype=float)
    
    # Generalized Tanimoto for count vectors:
    # T(A,B) = A·B / (|A|² + |B|² - A·B)
    dot = np.dot(vec_t, vec_c)
    denom = np.dot(vec_t, vec_t) + np.dot(vec_c, vec_c) - dot
    
    if denom == 0:
        return 1.0 if np.sum(vec_t) == 0 and np.sum(vec_c) == 0 else 0.0
    
    return float(dot / denom)


# ─── FreeSASA: Buried Surface Area ────────────────────────────────────────

def compute_bsa(
    pdb_path: str,
    peptide_chain: str = "P",
) -> tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Compute buried surface area (BSA) at the peptide-MHC interface.
    
    BSA = SASA(MHC alone) + SASA(peptide alone) - SASA(complex)
    
    Args:
        pdb_path: Path to pMHC PDB
        peptide_chain: Chain ID for peptide
    
    Returns:
        (bsa_total, bsa_polar, bsa_apolar) in Å², or (None, None, None) on failure
    """
    try:
        import freesasa
        
        # Compute SASA for the full complex
        structure_complex = freesasa.Structure(pdb_path)
        result_complex = freesasa.calc(structure_complex)
        
        # Compute SASA for peptide chain alone
        structure_pep = freesasa.Structure(pdb_path)
        # FreeSASA: select atoms by chain
        selections_pep = freesasa.selectArea(
            {"peptide": f"chain {peptide_chain}"},
            structure_complex,
            result_complex,
        )
        
        selections_mhc = freesasa.selectArea(
            {"mhc": f"chain M, chain N"},
            structure_complex,
            result_complex,
        )
        
        # BSA approximation from selections
        # For a more precise calculation, we'd need to compute SASA
        # on isolated chains, but this gives a reasonable estimate
        # TODO: 更精确的实现应该分别计算孤立链的 SASA
        
        total_sasa = result_complex.totalArea()
        pep_sasa = selections_pep.get("peptide", 0.0)
        mhc_sasa = selections_mhc.get("mhc", 0.0)
        
        # BSA = SASA(parts) - SASA(complex)
        # This is an approximation; precise BSA needs isolated-chain computation
        bsa_total = max(0.0, (pep_sasa + mhc_sasa) - total_sasa)
        
        return bsa_total, None, None  # polar/apolar 需要更细致的原子分类
        
    except Exception as exc:
        log.warning("FreeSASA computation failed for %s: %s", pdb_path, exc)
        return None, None, None


def compute_bsa_similarity(
    bsa_target: float,
    bsa_candidate: float,
    max_bsa: float = 1500.0,
) -> float:
    """
    Similarity of BSA values. 
    
    Normalized by max expected BSA for pMHC (~1500 Å² for 9-mer).
    Uses 1 - |diff|/max to convert to [0, 1] similarity.
    """
    diff = abs(bsa_target - bsa_candidate)
    return max(0.0, 1.0 - diff / max_bsa)


# ─── PRODIGY: Contact-Based Affinity ──────────────────────────────────────

def compute_prodigy_affinity(
    pdb_path: str,
    peptide_chain: str = "P",
    mhc_chain: str = "M",
) -> tuple[Optional[float], Optional[float]]:
    """
    Predict binding affinity (ΔG, Kd) using PRODIGY.
    
    PRODIGY uses interfacial contact counts (charged/polar/apolar)
    to estimate binding affinity. ~3-5s per structure.
    
    Args:
        pdb_path: Path to pMHC PDB
        peptide_chain: Chain ID for peptide
        mhc_chain: Chain ID for MHC heavy chain
    
    Returns:
        (delta_g_kcal, kd_molar) or (None, None) on failure
    
    Note:
        PRODIGY 的 Python API 可能需要适配。如果 prodigy-prot 包的 API
        与下面的调用方式不匹配，需要调整为命令行调用:
        `prodigy --contact_list <pdb> --chain1 P --chain2 M`
    """
    try:
        # PRODIGY Python API (may vary by version)
        # 如果 API 不可用，回退到命令行调用
        import subprocess
        result = subprocess.run(
            ["prodigy", pdb_path, "--selection", peptide_chain, mhc_chain],
            capture_output=True, text=True, timeout=30,
        )
        # Parse output for ΔG and Kd
        # PRODIGY output format: "Predicted binding affinity: -XX.X kcal/mol"
        dg = None
        kd = None
        for line in result.stdout.splitlines():
            if "binding affinity" in line.lower():
                # Extract number
                import re
                match = re.search(r"(-?\d+\.?\d*)\s*kcal", line)
                if match:
                    dg = float(match.group(1))
            if "dissociation constant" in line.lower():
                match = re.search(r"(\d+\.?\d*[eE][+-]?\d+)", line)
                if match:
                    kd = float(match.group(1))
        return dg, kd
        
    except Exception as exc:
        log.warning("PRODIGY failed for %s: %s", pdb_path, exc)
        return None, None


def compute_prodigy_similarity(
    dg_target: float,
    dg_candidate: float,
    max_dg_diff: float = 10.0,
) -> float:
    """
    ΔG similarity: candidates with similar predicted binding affinity
    to the target are more likely to compete for TCR binding.
    """
    diff = abs(dg_target - dg_candidate)
    return max(0.0, 1.0 - diff / max_dg_diff)


# ─── Composite Interface Similarity ──────────────────────────────────────

def compute_interface_similarity(
    target_pdb: str,
    candidate_pdb: str,
    peptide_chain: str = "P",
    mhc_chains: tuple = ("M", "N"),
    weights: dict = None,
) -> InterfaceSimilarity:
    """
    Compute all interface descriptors for both structures and return
    pairwise similarity scores.
    
    Args:
        target_pdb: Path to target pMHC PDB
        candidate_pdb: Path to candidate pMHC PDB
        peptide_chain: Peptide chain ID
        mhc_chains: MHC chain IDs
        weights: Optional weight dict, default {"plip": 0.4, "bsa": 0.3, "prodigy": 0.3}
    
    Returns:
        InterfaceSimilarity with per-descriptor and combined scores
    """
    if weights is None:
        weights = {"plip": 0.4, "bsa": 0.3, "prodigy": 0.3}
    
    # 1. PLIP fingerprints
    fp_target = compute_plip_fingerprint(target_pdb, peptide_chain, mhc_chains)
    fp_candidate = compute_plip_fingerprint(candidate_pdb, peptide_chain, mhc_chains)
    
    plip_sim = 0.0
    if fp_target is not None and fp_candidate is not None:
        plip_sim = compute_plip_tanimoto(fp_target, fp_candidate)
    
    # 2. BSA
    bsa_t, _, _ = compute_bsa(target_pdb, peptide_chain)
    bsa_c, _, _ = compute_bsa(candidate_pdb, peptide_chain)
    
    bsa_sim = 0.0
    if bsa_t is not None and bsa_c is not None:
        bsa_sim = compute_bsa_similarity(bsa_t, bsa_c)
    
    # 3. PRODIGY
    dg_t, _ = compute_prodigy_affinity(target_pdb, peptide_chain, mhc_chains[0])
    dg_c, _ = compute_prodigy_affinity(candidate_pdb, peptide_chain, mhc_chains[0])
    
    prodigy_sim = 0.0
    if dg_t is not None and dg_c is not None:
        prodigy_sim = compute_prodigy_similarity(dg_t, dg_c)
    
    # Combined
    combined = (
        weights["plip"] * plip_sim
        + weights["bsa"] * bsa_sim
        + weights["prodigy"] * prodigy_sim
    )
    
    return InterfaceSimilarity(
        plip_tanimoto=plip_sim,
        bsa_similarity=bsa_sim,
        prodigy_similarity=prodigy_sim,
        combined=combined,
    )


# ─── Internal Helpers ─────────────────────────────────────────────────────

def _is_peptide_mhc_interaction(interaction, peptide_chain, mhc_chains):
    """Check if an interaction spans the peptide-MHC interface."""
    # PLIP interaction objects have resnr, reschain attributes
    # Implementation depends on PLIP version; adapt as needed
    try:
        chains = set()
        if hasattr(interaction, 'resnr'):
            chains.add(getattr(interaction, 'reschain', ''))
        if hasattr(interaction, 'resnr_l'):
            chains.add(getattr(interaction, 'reschain_l', ''))
        # Check if one side is peptide and other is MHC
        has_pep = peptide_chain in chains
        has_mhc = any(c in chains for c in mhc_chains)
        return has_pep and has_mhc
    except Exception:
        return False


def _record_residue(residue_interactions, interaction, interaction_type):
    """Record a residue-level interaction for fingerprinting."""
    try:
        resnum = getattr(interaction, 'resnr', None)
        if resnum is not None:
            if resnum not in residue_interactions:
                residue_interactions[resnum] = []
            residue_interactions[resnum].append(interaction_type)
    except Exception:
        pass
```

### Step 1.3: 扩展 StructuralScore 数据模型

**文件**: `decoy_a/models.py`, 修改 `StructuralScore` (L139-149)

```python
class StructuralScore(BaseModel):
    """3D structural and electrostatic similarity assessment."""
    modeling_tool: str
    pdb_path: Optional[str] = None
    surface_correlation: float              # 现有：RMSD-based similarity
    rmsd: Optional[float] = None
    electrostatic_fingerprint_distance: Optional[float] = None
    
    # ═══ 新增界面描述符 ═══
    plip_tanimoto: Optional[float] = None           # PLIP 非共价相互作用 Tanimoto 相似度
    bsa_total: Optional[float] = None               # 界面埋藏面积 (Å²)
    bsa_similarity: Optional[float] = None           # BSA 相似度 [0-1]
    prodigy_dg: Optional[float] = None              # PRODIGY 预测 ΔG (kcal/mol)
    prodigy_similarity: Optional[float] = None       # ΔG 相似度 [0-1]
    interface_combined: Optional[float] = None       # 界面描述符综合分 [0-1]
    
    # Phase 2 预留
    esmif_score: Optional[float] = None             # ESM-IF 逆折叠兼容性分数
    hermes_score: Optional[float] = None            # HERMES 侧链倾向性分数
```

### Step 1.4: 集成到 compute_structure_similarity()

**文件**: `decoy_b/scanner.py`, 修改 `compute_structure_similarity()` (L336-499)

在现有的 RMSD 计算之后（约 L460），插入界面描述符计算：

```python
# ═══ 现有代码 L448-487: sim_a, sim_b, bf_corr 计算 ═══
# ... (保持不变) ...

# ═══ 新增: 界面描述符计算 ═══
interface_sim = None
try:
    from decoy_b.tools.interface_descriptors import compute_interface_similarity
    interface_sim = compute_interface_similarity(
        str(target_pdb), str(candidate_pdb),
        peptide_chain=PEPTIDE_CHAIN,
        mhc_chains=tuple(MHC_CHAINS),
    )
except ImportError:
    log.debug("Interface descriptors not available (missing PLIP/FreeSASA/PRODIGY)")
except Exception as exc:
    log.warning("Interface descriptor computation failed: %s", exc)

# ═══ 更新 combined score ═══
# 旧公式: combined = 0.5 * sim_a + 0.5 * sim_b  (纯 RMSD)
# 新公式: 加入界面描述符
if interface_sim is not None and interface_sim.combined > 0:
    # 三部分加权: RMSD 几何 (40%) + 界面描述符 (40%) + B-factor 柔性 (20%)
    rmsd_geo = 0.5 * sim_a + 0.5 * sim_b if (sim_a > 0 or sim_b > 0) else 0.0
    combined = (
        0.40 * rmsd_geo
        + 0.40 * interface_sim.combined
        + 0.20 * max(0.0, bf_corr)
    )
else:
    # 回退到旧逻辑
    if sim_a > 0 and sim_b > 0:
        combined = 0.5 * sim_a + 0.5 * sim_b
    elif sim_a > 0:
        combined = sim_a
    elif sim_b > 0:
        combined = sim_b
    else:
        combined = 0.0
    combined = max(combined, bf_corr)

return StructuralScore(
    modeling_tool=tool,
    pdb_path=str(candidate_pdb),
    surface_correlation=max(0.0, combined),
    rmsd=rmsd_for_report,
    # 新增字段
    plip_tanimoto=interface_sim.plip_tanimoto if interface_sim else None,
    bsa_similarity=interface_sim.bsa_similarity if interface_sim else None,
    prodigy_similarity=interface_sim.prodigy_similarity if interface_sim else None,
    interface_combined=interface_sim.combined if interface_sim else None,
)
```

### Step 1.5: 更新最终评分公式

**文件**: `decoy_b/scanner.py`, 修改 `scan_decoy_b()` 中 hit 组装 (L928)

```python
# 旧:
# combined = 0.4 * cos_sim + 0.6 * structural.surface_correlation

# 新: surface_correlation 已经包含了界面描述符（在 compute_structure_similarity 中融合）
# 所以这行不需要改，但权重可以调整：
if structural is not None:
    combined = 0.35 * cos_sim + 0.65 * structural.surface_correlation
else:
    combined = cos_sim
```

### Step 1.6: 单元测试

**新建文件**: `tests/test_interface_descriptors.py`

```python
"""
Test interface descriptor computation on existing tFold PDB outputs.

使用 data/decoy_b/ 下已有的 GILGFVFTL PDB 结构进行测试。
"""
import pytest
from pathlib import Path

DATA_DIR = Path("data/decoy_b/stepwise_demo/tfold_pdbs")

def test_plip_fingerprint():
    """验证 PLIP 能正确解析 tFold 输出的 PDB 格式 (Chain M/N/P)."""
    from decoy_b.tools.interface_descriptors import compute_plip_fingerprint
    
    pdbs = list(DATA_DIR.glob("*.pdb"))
    if not pdbs:
        pytest.skip("No PDB files in stepwise_demo")
    
    fp = compute_plip_fingerprint(str(pdbs[0]))
    assert fp is not None
    assert fp.hbond_count >= 0
    assert fp.hydrophobic_count >= 0
    print(f"PLIP: {fp.hbond_count} H-bonds, {fp.hydrophobic_count} hydrophobic, "
          f"{fp.salt_bridge_count} salt bridges")

def test_plip_tanimoto_self():
    """同一结构的 Tanimoto 应为 1.0."""
    from decoy_b.tools.interface_descriptors import (
        compute_plip_fingerprint, compute_plip_tanimoto,
    )
    pdbs = list(DATA_DIR.glob("*.pdb"))
    if not pdbs:
        pytest.skip("No PDB files")
    
    fp = compute_plip_fingerprint(str(pdbs[0]))
    if fp is None:
        pytest.skip("PLIP not available")
    
    sim = compute_plip_tanimoto(fp, fp)
    assert sim == pytest.approx(1.0)

def test_freesasa_bsa():
    """验证 FreeSASA 能计算 BSA."""
    from decoy_b.tools.interface_descriptors import compute_bsa
    
    pdbs = list(DATA_DIR.glob("*.pdb"))
    if not pdbs:
        pytest.skip("No PDB files")
    
    bsa, _, _ = compute_bsa(str(pdbs[0]))
    assert bsa is not None
    assert bsa > 0  # pMHC 复合物应有正的 BSA
    print(f"BSA = {bsa:.1f} Å²")

def test_full_interface_similarity():
    """端到端测试: 两个不同 pMHC 的界面相似度."""
    from decoy_b.tools.interface_descriptors import compute_interface_similarity
    
    pdbs = list(DATA_DIR.glob("*.pdb"))
    if len(pdbs) < 2:
        pytest.skip("Need at least 2 PDB files")
    
    sim = compute_interface_similarity(str(pdbs[0]), str(pdbs[1]))
    assert 0.0 <= sim.combined <= 1.0
    assert 0.0 <= sim.plip_tanimoto <= 1.0
    print(f"Interface similarity: PLIP={sim.plip_tanimoto:.3f}, "
          f"BSA={sim.bsa_similarity:.3f}, combined={sim.combined:.3f}")
```

### Step 1.7: 验证与回归测试

```bash
# 1. 运行新的单元测试
pytest tests/test_interface_descriptors.py -v

# 2. 在 GILGFVFTL 靶标上运行完整管线，确认无回归
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01" --skip-structural
# 对比新旧结果的 Top 20 排名变化

# 3. 带结构的完整测试（需要 tFold）
python -m decoy_b run --target GILGFVFTL --hla "HLA-A*02:01"
```

---

## Phase 2: 深层特征（ESM-IF + HERMES + 静电势）

> 目标：引入学习型描述符，捕捉侧链兼容性和静电互补性
> 预计工期：3-4 周
> 前置条件：Phase 1 完成且通过回归测试

### Step 2.1: ESM-IF 逆折叠兼容性分数

**原理**: ESM-IF (Inverse Folding) 给定一个 3D 骨架，预测每个位置最可能的氨基酸序列。对于一个候选肽段，它在靶标 pMHC 骨架上的 ESM-IF perplexity 越低，说明该肽段越"适合"这个结构环境——即越可能被同一个 TCR 识别。

**安装**:
```bash
pip install fair-esm
# ESM-IF1 模型权重会自动下载 (~700MB)
```

**新增函数** (添加到 `interface_descriptors.py`):
```python
def compute_esmif_compatibility(
    pdb_path: str,
    candidate_sequence: str,
    peptide_chain: str = "P",
) -> Optional[float]:
    """
    Compute ESM-IF inverse folding compatibility score.
    
    Given the target pMHC backbone structure, how well does the 
    candidate peptide sequence "fit" structurally?
    
    Lower perplexity = better fit = higher cross-reactivity risk.
    
    Returns: normalized score [0, 1] where 1 = perfect fit
    """
    try:
        import esm
        # Load ESM-IF1
        model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        model = model.eval()
        
        # Extract backbone coordinates for peptide chain
        # Compute log-likelihood of candidate sequence given structure
        # Normalize to [0, 1]
        
        # TODO: 具体实现需要参考 ESM-IF 的 API
        # 核心调用: model.forward(coords, tokens) → logits
        # perplexity = exp(-mean(log_probs))
        # score = 1.0 / (1.0 + perplexity)  # sigmoid-like normalization
        
        pass  # 实现细节
    except Exception as exc:
        log.warning("ESM-IF failed: %s", exc)
        return None
```

### Step 2.2: HERMES 侧链倾向性

**原理**: HERMES 是等变 GNN，预测每个残基位置在给定结构环境中的侧链旋转异构体偏好。如果靶标和候选 pMHC 在 TCR 接触面的侧链倾向性分布相似，说明它们呈现给 TCR 的"化学面孔"相似。

**安装**:
```bash
pip install hermes-protein  # 或从 GitHub 安装
```

**设计要点**:
- 提取 TCR 接触面（p4-p8）每个位置的侧链 χ 角分布
- 计算两个 pMHC 的侧链倾向性向量的 KL 散度
- 转换为相似度: `sim = 1 / (1 + KL_divergence)`

### Step 2.3: APBS 静电势

**原理**: TCR 识别强烈依赖于 pMHC 表面的静电分布模式。两个 pMHC 如果在 TCR 接触面的静电势分布相似，即使序列不同也可能被同一 TCR 识别。

**安装**:
```bash
conda install -c conda-forge apbs  # 或下载二进制
pip install pdb2pqr  # PDB → PQR 转换（添加电荷和半径）
```

**计算流程**:
```
PDB → pdb2pqr (添加电荷) → APBS (解 PB 方程) → 静电势网格 (.dx)
→ 在 TCR 接触面上采样静电势 → 计算两个 pMHC 的相关系数
```

**设计要点**:
- 只在 peptide 上表面（朝向 TCR 的一侧）采样静电势
- 使用球面采样点（距肽段 CA 4-6 Å 的半球壳）
- 相似度 = Pearson correlation of sampled potentials

### Step 2.4: 更新 StructuralScore 和评分公式

扩展 `compute_structure_similarity()` 以包含 Phase 2 描述符：

```python
# Phase 2 评分公式
# 6 维特征:
#   1. RMSD 几何 (sim_a + sim_b)
#   2. B-factor 柔性
#   3. PLIP 非共价相互作用
#   4. BSA 埋藏面积
#   5. ESM-IF 逆折叠兼容性
#   6. 静电势相关

# 权重建议 (基于 deep research 报告):
weights_phase2 = {
    "rmsd_geo": 0.20,       # 降权：粗糙但鲁棒
    "bfactor": 0.05,        # 降权：预测结构的 B-factor 不太可靠
    "plip": 0.20,           # 非共价相互作用模式
    "bsa": 0.10,            # 界面面积
    "esmif": 0.25,          # 学习型，对预测结构鲁棒
    "electrostatic": 0.20,  # 静电互补性
}
```

---

## Phase 3: LTR 模型替代线性加权

> 目标：用数据驱动的 Learning-to-Rank 模型替代手工权重
> 预计工期：1-2 月
> 前置条件：Phase 2 完成，有足够的特征维度

### Step 3.1: 训练数据收集

| 数据源 | 内容 | 用途 | 获取方式 |
|--------|------|------|---------|
| **BATCAVE** | 25 个表位、>100 TCR 的深度突变扫描激活数据 | 正负样本训练 | 公开下载 |
| **ATLAS** | TCR-pMHC 结构 + 亲和力 (Kd/ΔG) | 能量特征验证 | 公开下载 |
| **VDJdb** | 实验验证的 TCR-peptide 配对 | 正样本 | https://vdjdb.cdr3.net |
| **STCRDab** | TCR 结构数据库 | 结构验证 | 公开下载 |
| **IEDB** | 免疫表位数据 | 交叉验证 | https://www.iedb.org |
| **Decoy C 文献库** | 我们自己的 250 条临床/实验证据 | 真阳性验证 | 已有 |

**数据准备流程**:
```
1. 从 BATCAVE 提取 (TCR, peptide, activation_score) 三元组
2. 对每个 peptide，用 tFold 预测 pMHC 结构
3. 计算 Phase 1 + Phase 2 的所有界面描述符
4. 标签: activation_score > threshold → positive, else → negative
5. 生成训练矩阵: [N_samples × ~10 features + 1 label]
```

### Step 3.2: XGBoost LTR 模型

```python
"""
decoy_b/ltr_scorer.py — Learning-to-Rank safety scorer
"""
import xgboost as xgb
import numpy as np
from pathlib import Path

FEATURE_COLUMNS = [
    "atchley_cosine",           # Phase 0: 现有
    "rmsd_similarity",          # Phase 0: 现有
    "bfactor_correlation",      # Phase 0: 现有
    "plip_tanimoto",            # Phase 1: 新增
    "bsa_similarity",           # Phase 1: 新增
    "prodigy_similarity",       # Phase 1: 新增
    "esmif_compatibility",      # Phase 2: 新增
    "electrostatic_correlation",# Phase 2: 新增
    "el_rank",                  # HLA 呈递强度
    "tpm_weight",               # 表达谱风险
]

class LTRScorer:
    """XGBoost-based Learning-to-Rank scorer for TCR cross-reactivity."""
    
    def __init__(self, model_path: Optional[Path] = None):
        self.model = None
        if model_path and model_path.exists():
            self.model = xgb.Booster()
            self.model.load_model(str(model_path))
    
    def train(self, X: np.ndarray, y: np.ndarray, groups: list):
        """
        Train LTR model.
        
        Args:
            X: Feature matrix [N, len(FEATURE_COLUMNS)]
            y: Relevance labels (0=non-binder, 1=weak, 2=strong)
            groups: Group sizes for LTR (each group = one target peptide)
        """
        dtrain = xgb.DMatrix(X, label=y)
        dtrain.set_group(groups)
        
        params = {
            "objective": "rank:ndcg",
            "eval_metric": "ndcg@20",
            "eta": 0.1,
            "max_depth": 6,
            "min_child_weight": 10,
            "subsample": 0.8,
            "colsample_bytree": 0.8,
        }
        
        self.model = xgb.train(params, dtrain, num_boost_round=200)
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict risk scores for candidates."""
        if self.model is None:
            raise RuntimeError("Model not trained or loaded")
        dtest = xgb.DMatrix(X)
        return self.model.predict(dtest)
    
    def save(self, path: Path):
        self.model.save_model(str(path))
```

### Step 3.3: 集成到 risk_scorer.py

```python
# risk_scorer.py 中新增逻辑:

def compute_risk_score_v2(hit: DecoyBHit, ltr_model: Optional[LTRScorer] = None):
    """
    V2 风险评分: 如果有训练好的 LTR 模型则使用，否则回退到 V1 线性加权。
    """
    if ltr_model is not None and ltr_model.model is not None:
        features = _extract_feature_vector(hit)
        return float(ltr_model.predict(features.reshape(1, -1))[0])
    else:
        # V1 回退
        return compute_risk_score(
            hit.similarity_score, hit.el_rank, compute_tpm_weight(hit.expression)
        )
```

### Step 3.4: 评估与验证

```bash
# 交叉验证
python -m decoy_b evaluate --model ltr --cv 5 --data batcave

# 与 V1 线性加权对比
python -m decoy_b evaluate --compare v1 v2 --target GILGFVFTL

# 输出指标:
# - NDCG@20 (排序质量)
# - Recall@50 (前 50 中有多少真阳性)
# - Precision@10 (前 10 中有多少是真正危险的)
```

---

## 里程碑检查清单

| 里程碑 | 验收标准 | 优先级 |
|--------|---------|--------|
| **M1**: PLIP 集成 | tFold PDB 能正确提取 H-bond/hydrophobic/salt bridge 计数 | P0 |
| **M2**: FreeSASA 集成 | BSA 计算值在合理范围 (500-1500 Å² for 9-mer pMHC) | P0 |
| **M3**: PRODIGY 集成 | ΔG 预测值在合理范围 (-5 to -15 kcal/mol) | P1 |
| **M4**: 新评分公式上线 | GILGFVFTL Top20 排名合理 + 无回归 | P0 |
| **M5**: ESM-IF 集成 | 逆折叠 perplexity 能区分相似/不相似肽段 | P1 |
| **M6**: 静电势集成 | 两个已知交叉反应对的静电相似度 > 0.7 | P2 |
| **M7**: LTR 模型训练 | BATCAVE 上 NDCG@20 > 0.8 | P2 |
| **M8**: LTR 模型部署 | 替代线性加权，生产环境可用 | P2 |

---

## 风险与缓解

| 风险 | 影响 | 缓解措施 |
|------|------|---------|
| PLIP 无法解析 tFold PDB 格式 (Chain M/N/P) | Phase 1 阻塞 | tFold PDB 链命名标准，应兼容；备选: 用 BioPython 自行计算接触 |
| PRODIGY 不支持单链 peptide 作为 partner | Phase 1 降级 | 改用命令行模式 或 直接用 PRODIGY 的核心算法（接触计数）手动实现 |
| ESM-IF 推理太慢 (>10s/结构) | Phase 2 吞吐不足 | 仅对 Top 200 候选运行，或 batch 推理 + GPU 加速 |
| BATCAVE 数据量不足训练 LTR | Phase 3 过拟合 | 合并 VDJdb + IEDB 扩大训练集；使用特征选择减少维度 |
| 预测结构的侧链建模不准确 | 所有 Phase | 优先使用对侧链鲁棒的描述符（PLIP 接触计数、ESM-IF、图方法）|

---

## 附录: 关键参考文献

从 Deep Research 报告中提取的最重要参考：

1. **BATCAVE** — 平衡正负样本的 TCR 激活数据库，25 表位、>100 TCR
2. **STAG** — 结构感知 GNN，RBF 编码距离边，专为 TCR-pMHC 设计
3. **t2pmhc** — 解决 unseen antigen 泛化问题的结构 GNN
4. **PPIscreenML** — AF 置信度 + Rosetta 能量的 LTR 实现示例
5. **IMPRINT/MaSIF** — 表面网格学习，用于界面相似度（Phase 2 备选）
6. **ESM-IF1** — 逆折叠模型，对预测结构鲁棒
7. **PRODIGY** — 基于接触的快速亲和力预测
8. **PLIP** — 蛋白质相互作用分析器，提取非共价相互作用指纹
