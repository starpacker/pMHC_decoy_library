#!/usr/bin/env python3
"""
Build interactive visualizations for the Decoy Library pipeline.

Generates:
  1. decoy_3d_comparison.html  — Enhanced Decoy B/D 3D structure viewer (3Dmol.js)
  2. decoy_a_overview.html     — Decoy A bubble chart + sequence alignment

Usage:
    PYTHONUTF8=1 python build_viz.py [--target GILGFVFTL] [--dir data/GILGFVFTL_summary]
"""
from __future__ import annotations

import argparse
import json
import re
import csv
from pathlib import Path
from typing import Any


# ═══════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════

def _parse_pdb_remarks(pdb_text: str) -> dict:
    """Extract tFold confidence scores from REMARK lines."""
    info = {}
    for line in pdb_text.splitlines():
        if not line.startswith("REMARK 250"):
            continue
        if "lDDT" in line:
            m = re.search(r"([\d.]+)\s*$", line)
            if m:
                info["lddt"] = float(m.group(1))
        if "pTM" in line and "ipTM" not in line:
            m = re.search(r"([\d.]+)\s*$", line)
            if m:
                info["ptm"] = float(m.group(1))
        if "ipTM" in line:
            m = re.search(r"([\d.]+)\s*$", line)
            if m:
                info["iptm"] = float(m.group(1))
    return info


def _parse_peptide_bfactors(pdb_text: str, chain: str = "P") -> list[dict]:
    """Extract per-residue B-factors (pLDDT) for the peptide chain."""
    residues = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        resnum = int(line[22:26])
        resname = line[17:20].strip()
        bfactor = float(line[60:66])
        residues[resnum] = {"resnum": resnum, "resname": resname, "plddt": bfactor}
    return [residues[k] for k in sorted(residues)]


# ═══════════════════════════════════════════════════════════════════════
#  1. Enhanced Decoy B/D 3D Viewer
# ═══════════════════════════════════════════════════════════════════════

def build_3d_viewer(base: Path, target_seq: str) -> Path:
    struct_dir_b = base / "Decoy_B" / "3D_structures"
    struct_dir_d = base / "Decoy_D" / "3D_structures"
    target_pdb_path = struct_dir_b / f"pmhc_{target_seq}_HLA-A0201.pdb"

    if not target_pdb_path.exists():
        # try Decoy D target
        target_pdb_path = base / "Decoy_D" / "target_pdb" / f"pmhc_{target_seq}_HLA-A0201.pdb"
    if not target_pdb_path.exists():
        print(f"  [WARN] Target PDB not found for {target_seq}, skipping 3D viewer")
        return None

    target_pdb = target_pdb_path.read_text(encoding="utf-8")
    target_remarks = _parse_pdb_remarks(target_pdb)
    target_bfactors = _parse_peptide_bfactors(target_pdb)

    # Collect Decoy B
    ranked_path = base / "Decoy_B" / "final_ranked_decoys.json"
    decoy_b_ranked = []
    if ranked_path.exists():
        decoy_b_ranked = json.loads(ranked_path.read_text(encoding="utf-8"))

    entries = []

    for e in decoy_b_ranked:
        seq = e["sequence"]
        pdb_file = struct_dir_b / f"pmhc_{seq}_HLA-A0201.pdb"
        if pdb_file.exists() and seq != target_seq:
            pdb_text = pdb_file.read_text(encoding="utf-8")
            entries.append({
                "sequence": seq,
                "source": "Decoy B",
                "hamming_distance": e.get("hamming_distance", "?"),
                "physchem_sim": f"{e.get('physicochemical_similarity', 0):.3f}",
                "struct_sim": f"{e.get('structural_similarity', 0):.3f}",
                "risk_score": f"{e.get('total_risk_score', 0):.2f}",
                "genes": ", ".join(e.get("gene_symbols", [])[:3]),
                "pdb": pdb_text,
                "remarks": _parse_pdb_remarks(pdb_text),
                "bfactors": _parse_peptide_bfactors(pdb_text),
            })
        if len(entries) >= 10:
            break

    # Collect Decoy D
    if struct_dir_d.exists():
        for pdb_file in sorted(struct_dir_d.glob("pmhc_*.pdb")):
            seq = pdb_file.stem.replace("pmhc_", "").replace("_HLA-A0201", "")
            if seq == target_seq:
                continue
            pdb_text = pdb_file.read_text(encoding="utf-8")
            hd = sum(1 for a, b in zip(target_seq, seq) if a != b) if len(seq) == len(target_seq) else "N/A"
            entries.append({
                "sequence": seq,
                "source": "Decoy D (MPNN)",
                "hamming_distance": hd,
                "physchem_sim": "-",
                "struct_sim": "-",
                "risk_score": "-",
                "genes": "MPNN designed",
                "pdb": pdb_text,
                "remarks": _parse_pdb_remarks(pdb_text),
                "bfactors": _parse_peptide_bfactors(pdb_text),
            })
            if len(entries) >= 20:
                break

    entries_json = json.dumps(entries, ensure_ascii=False)
    target_pdb_json = json.dumps(target_pdb, ensure_ascii=False)
    target_remarks_json = json.dumps(target_remarks, ensure_ascii=False)
    target_bfactors_json = json.dumps(target_bfactors, ensure_ascii=False)

    html = _build_3d_html(target_seq, target_pdb_json, target_remarks_json,
                          target_bfactors_json, entries_json, len(entries))

    out_path = base / "decoy_3d_comparison.html"
    out_path.write_text(html, encoding="utf-8")
    print(f"  [OK] 3D viewer: {out_path} ({len(entries)} decoys, {len(html):,} bytes)")
    return out_path


def _build_3d_html(target_seq, target_pdb_json, target_remarks_json,
                   target_bfactors_json, entries_json, n_entries):
    """Generate the full HTML for the enhanced 3D viewer."""
    # Using a multi-line string approach that avoids shell escaping issues
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>pMHC Decoy 3D — {target_seq}</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
:root {{
  --bg-primary: #0d1117; --bg-secondary: #161b22; --bg-tertiary: #21262d;
  --border: #30363d; --text-primary: #f0f6fc; --text-secondary: #c9d1d9;
  --text-muted: #8b949e; --accent-blue: #58a6ff; --accent-green: #7ee787;
  --accent-orange: #ffa657; --accent-red: #ff7b72; --accent-purple: #d2a8ff;
  --accent-yellow: #e3b341;
}}
* {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:"Segoe UI",system-ui,-apple-system,sans-serif; background:var(--bg-primary); color:var(--text-secondary); overflow:hidden; }}

/* Header */
.header {{ background:var(--bg-secondary); padding:14px 24px; border-bottom:1px solid var(--border); display:flex; align-items:center; justify-content:space-between; }}
.header-left h1 {{ font-size:17px; color:var(--accent-blue); font-weight:600; }}
.header-left p {{ font-size:12px; color:var(--text-muted); margin-top:2px; }}
.header-right {{ display:flex; gap:8px; }}
.mode-btn {{ padding:5px 14px; border:1px solid var(--border); border-radius:20px; background:var(--bg-tertiary); color:var(--text-muted); cursor:pointer; font-size:11px; font-weight:500; transition:all .15s; }}
.mode-btn:hover {{ border-color:var(--accent-blue); color:var(--text-secondary); }}
.mode-btn.active {{ background:rgba(88,166,255,0.15); border-color:var(--accent-blue); color:var(--accent-blue); }}

/* Layout */
.container {{ display:flex; height:calc(100vh - 56px); }}
.sidebar {{ width:320px; background:var(--bg-secondary); border-right:1px solid var(--border); display:flex; flex-direction:column; flex-shrink:0; }}
.sidebar-header {{ padding:10px 12px; border-bottom:1px solid var(--border); }}
.filter-pills {{ display:flex; gap:4px; }}
.pill {{ padding:4px 12px; border-radius:12px; border:1px solid var(--border); background:transparent; color:var(--text-muted); cursor:pointer; font-size:11px; transition:all .12s; }}
.pill:hover {{ border-color:var(--accent-blue); }}
.pill.active {{ background:rgba(88,166,255,0.12); border-color:var(--accent-blue); color:var(--accent-blue); }}
.card-list {{ flex:1; overflow-y:auto; padding:8px; }}
.card-list::-webkit-scrollbar {{ width:6px; }}
.card-list::-webkit-scrollbar-thumb {{ background:var(--border); border-radius:3px; }}

/* Cards */
.card {{ padding:10px 12px; margin-bottom:5px; border-radius:8px; cursor:pointer; border:1px solid transparent; transition:all .12s; position:relative; }}
.card:hover {{ background:rgba(88,166,255,0.06); border-color:var(--border); }}
.card.active {{ background:rgba(255,166,87,0.08); border-color:var(--accent-orange); }}
.card .seq {{ font-family:"Cascadia Code",Consolas,monospace; font-size:14px; font-weight:700; letter-spacing:1.5px; }}
.card .seq .m {{ color:var(--accent-green); }} .card .seq .x {{ color:var(--accent-red); }}
.card .meta {{ font-size:10px; color:var(--text-muted); margin-top:3px; display:flex; gap:8px; }}
.badge {{ display:inline-block; padding:1px 7px; border-radius:10px; font-size:9px; font-weight:700; letter-spacing:.5px; vertical-align:middle; margin-right:4px; }}
.badge-b {{ background:rgba(88,166,255,0.15); color:var(--accent-blue); }}
.badge-d {{ background:rgba(255,166,87,0.15); color:var(--accent-orange); }}

/* Viewer area */
.viewer-area {{ flex:1; display:flex; flex-direction:column; min-width:0; }}
.viewers {{ flex:1; display:flex; position:relative; }}
.viewer-box {{ flex:1; position:relative; }}
.vlabel {{ position:absolute; top:10px; left:14px; z-index:10; background:rgba(0,0,0,0.8); backdrop-filter:blur(8px); padding:5px 12px; border-radius:6px; font-size:12px; font-weight:500; pointer-events:none; }}
.vlabel.tgt {{ color:var(--accent-green); }}
.vlabel.dec {{ color:var(--accent-orange); }}
.divider {{ width:1px; background:var(--border); }}

/* Info bar */
.info-bar {{ background:var(--bg-secondary); border-top:1px solid var(--border); padding:10px 20px; display:flex; gap:24px; font-size:12px; flex-wrap:wrap; align-items:center; }}
.info-item {{ display:flex; gap:4px; align-items:center; }}
.info-item .l {{ color:var(--text-muted); }}
.info-item .v {{ color:var(--text-primary); font-weight:600; }}

/* Residue heatmap bar */
.heatmap-bar {{ background:var(--bg-secondary); border-top:1px solid var(--border); padding:10px 20px; }}
.heatmap-title {{ font-size:11px; color:var(--text-muted); margin-bottom:6px; font-weight:500; }}
.heatmap-row {{ display:flex; gap:2px; align-items:end; }}
.heatmap-row.labels {{ margin-bottom:2px; }}
.res-cell {{ width:52px; text-align:center; border-radius:4px; font-family:"Cascadia Code",Consolas,monospace; font-size:11px; font-weight:600; padding:4px 0; transition:all .15s; }}
.res-cell.header {{ background:transparent; color:var(--text-muted); font-size:9px; font-weight:400; }}
.plddt-bar {{ height:4px; border-radius:2px; margin-top:2px; }}

/* Confidence badge */
.conf-badge {{ display:inline-flex; align-items:center; gap:3px; padding:2px 8px; border-radius:10px; font-size:10px; font-weight:600; }}
.conf-badge .dot {{ width:6px; height:6px; border-radius:50%; }}

/* Legend panel */
.sidebar-footer {{ padding:10px 12px; border-top:1px solid var(--border); font-size:10px; color:var(--text-muted); }}
.legend-row {{ display:flex; align-items:center; gap:6px; margin:3px 0; }}
.legend-swatch {{ width:12px; height:12px; border-radius:2px; flex-shrink:0; }}

/* Color scale overlay (pLDDT / Surface modes) */
.color-scale {{ position:absolute; bottom:14px; right:14px; z-index:10; background:rgba(0,0,0,0.82); backdrop-filter:blur(8px); border:1px solid var(--border); border-radius:8px; padding:10px 14px; pointer-events:none; display:none; }}
.color-scale .cs-title {{ font-size:10px; color:var(--text-muted); font-weight:500; margin-bottom:6px; letter-spacing:.3px; }}
.color-scale .cs-bar {{ width:160px; height:12px; border-radius:3px; }}
.color-scale .cs-labels {{ display:flex; justify-content:space-between; font-size:9px; color:var(--text-muted); margin-top:3px; }}
</style>
</head>
<body>
<div class="header">
  <div class="header-left">
    <h1>pMHC Decoy Structure Comparison</h1>
    <p>Target: <b style="color:var(--accent-green);font-family:monospace;letter-spacing:1px">{target_seq}</b> &nbsp;&bull;&nbsp; HLA-A*02:01</p>
  </div>
  <div class="header-right">
    <button class="mode-btn active" id="btn-sidebyside" onclick="setMode('side')">Side-by-Side</button>
    <button class="mode-btn" id="btn-super" onclick="setMode('super')">Superpose</button>
    <button class="mode-btn" id="btn-plddt" onclick="setMode('plddt')">pLDDT Confidence</button>
    <button class="mode-btn" id="btn-surface" onclick="setMode('surface')">Surface Compare</button>
  </div>
</div>
<div class="container">
  <div class="sidebar">
    <div class="sidebar-header">
      <div class="filter-pills">
        <button class="pill active" onclick="filterCards('all',this)">All ({n_entries})</button>
        <button class="pill" onclick="filterCards('Decoy B',this)">Decoy B</button>
        <button class="pill" onclick="filterCards('Decoy D',this)">Decoy D</button>
      </div>
    </div>
    <div class="card-list" id="card-list"></div>
    <div class="sidebar-footer">
      <div class="legend-row"><div class="legend-swatch" style="background:var(--accent-green)"></div>Target peptide</div>
      <div class="legend-row"><div class="legend-swatch" style="background:var(--accent-orange)"></div>Decoy peptide</div>
      <div class="legend-row"><div class="legend-swatch" style="background:var(--accent-red)"></div>Mismatch residues</div>
      <div class="legend-row"><div class="legend-swatch" style="background:#30363d"></div>MHC groove (translucent)</div>
      <div style="margin-top:6px;line-height:1.5">Drag = rotate &bull; Scroll = zoom<br>Right-drag = pan &bull; Middle = slab</div>
    </div>
  </div>
  <div class="viewer-area">
    <div class="viewers">
      <div class="viewer-box" id="box-left">
        <div class="vlabel tgt" id="label-left">Target: {target_seq}</div>
        <div class="color-scale" id="cs-left">
          <div class="cs-title" id="cs-left-title"></div>
          <div class="cs-bar" id="cs-left-bar"></div>
          <div class="cs-labels"><span id="cs-left-lo"></span><span id="cs-left-hi"></span></div>
        </div>
        <div id="viewer-left" style="width:100%;height:100%"></div>
      </div>
      <div class="divider" id="divider"></div>
      <div class="viewer-box" id="box-right">
        <div class="vlabel dec" id="label-right">Select a decoy</div>
        <div class="color-scale" id="cs-right">
          <div class="cs-title" id="cs-right-title"></div>
          <div class="cs-bar" id="cs-right-bar"></div>
          <div class="cs-labels"><span id="cs-right-lo"></span><span id="cs-right-hi"></span></div>
        </div>
        <div id="viewer-right" style="width:100%;height:100%"></div>
      </div>
    </div>
    <div id="heatmap-panel" class="heatmap-bar" style="display:none">
      <div class="heatmap-title">Per-Residue Comparison &nbsp;
        <span id="heatmap-subtitle" style="color:var(--text-primary)"></span>
      </div>
      <div class="heatmap-row labels" id="hm-labels"></div>
      <div class="heatmap-row" id="hm-target"></div>
      <div class="heatmap-row" id="hm-decoy"></div>
      <div class="heatmap-row" id="hm-plddt-tgt" style="margin-top:4px"></div>
      <div class="heatmap-row" id="hm-plddt-dec"></div>
    </div>
    <div class="info-bar" id="info-bar">
      <div class="info-item"><span class="l">Source:</span><span class="v" id="i-src">&mdash;</span></div>
      <div class="info-item"><span class="l">HD:</span><span class="v" id="i-hd">&mdash;</span></div>
      <div class="info-item"><span class="l">Physchem:</span><span class="v" id="i-phys">&mdash;</span></div>
      <div class="info-item"><span class="l">Structural:</span><span class="v" id="i-struct">&mdash;</span></div>
      <div class="info-item"><span class="l">Risk:</span><span class="v" id="i-risk">&mdash;</span></div>
      <div class="info-item"><span class="l">Gene:</span><span class="v" id="i-gene">&mdash;</span></div>
      <div class="info-item" id="conf-container"></div>
    </div>
  </div>
</div>

<script>
var TGT = "{target_seq}";
var TGT_PDB = ''' + target_pdb_json + ''';
var TGT_RMK = ''' + target_remarks_json + ''';
var TGT_BF = ''' + target_bfactors_json + ''';
var DECOYS = ''' + entries_json + ''';

var vL, vR, mode = "side", curIdx = -1;

/* ── AA property color (charge-based for ESP-like view) ── */
var AA_CHARGE = {D:-1,E:-1,K:1,R:1,H:0.5,S:0,T:0,N:0,Q:0,C:0,G:0,P:0,A:0,V:0,I:0,L:0,M:0,F:0,Y:0,W:0};
var AA_HYDRO = {I:4.5,V:4.2,L:3.8,F:2.8,C:2.5,M:1.9,A:1.8,G:-0.4,T:-0.7,S:-0.8,W:-0.9,Y:-1.3,P:-1.6,H:-3.2,D:-3.5,E:-3.5,N:-3.5,Q:-3.5,K:-3.9,R:-4.5};
var THREE2ONE = {ALA:"A",CYS:"C",ASP:"D",GLU:"E",PHE:"F",GLY:"G",HIS:"H",ILE:"I",LYS:"K",LEU:"L",MET:"M",ASN:"N",PRO:"P",GLN:"Q",ARG:"R",SER:"S",THR:"T",VAL:"V",TRP:"W",TYR:"Y"};

function plddtColor(v) {
  if (v >= 0.90) return "#4477ff";
  if (v >= 0.70) return "#65cbf3";
  if (v >= 0.50) return "#e3b341";
  return "#ff7b72";
}
function chargeColor(aa) {
  var c = AA_CHARGE[aa] || 0;
  if (c > 0.3) return "#6eaaff";
  if (c < -0.3) return "#ff6e6e";
  return "#e8e8e8";
}
function hydroColor(aa) {
  var h = AA_HYDRO[aa] || 0;
  if (h > 2) return "#ffa657";
  if (h > 0) return "#e3b341";
  if (h > -2) return "#65cbf3";
  return "#58a6ff";
}

function colorSeq(seq) {
  var h = "";
  for (var i = 0; i < seq.length; i++)
    h += '<span class="' + (i < TGT.length && seq[i] === TGT[i] ? "m" : "x") + '">' + seq[i] + "</span>";
  return h;
}

function confBadge(rmk) {
  if (!rmk || !rmk.iptm) return "";
  var v = rmk.iptm, c = v > 0.85 ? "var(--accent-green)" : v > 0.7 ? "var(--accent-yellow)" : "var(--accent-red)";
  return '<span class="conf-badge" style="background:rgba(255,255,255,0.06)"><span class="dot" style="background:'+c+'"></span>ipTM '+v.toFixed(3)+'</span>';
}

/* ── Card list ── */
function buildCards() {
  var list = document.getElementById("card-list");
  list.innerHTML = "";
  DECOYS.forEach(function(d, i) {
    var badge = d.source.indexOf("MPNN") >= 0 ? '<span class="badge badge-d">D</span>' : '<span class="badge badge-b">B</span>';
    var card = document.createElement("div");
    card.className = "card"; card.dataset.source = d.source; card.dataset.idx = i;
    card.innerHTML = badge + '<span class="seq">' + colorSeq(d.sequence) + '</span>'
      + '<div class="meta"><span>HD=' + d.hamming_distance + '</span><span>' + d.genes + '</span></div>';
    card.onclick = function() { selectDecoy(i); };
    list.appendChild(card);
  });
}
function filterCards(src, btn) {
  document.querySelectorAll(".pill").forEach(function(b){b.classList.remove("active")});
  btn.classList.add("active");
  document.querySelectorAll(".card").forEach(function(c) {
    c.style.display = (src === "all" || c.dataset.source.indexOf(src) >= 0) ? "block" : "none";
  });
}

/* ── Rendering ── */
function styleMHC(v, model) {
  var sel = model !== undefined ? {model:model} : {};
  v.setStyle(Object.assign({chain:"M"}, sel), {cartoon:{color:"#30363d",opacity:0.20}});
  v.setStyle(Object.assign({chain:"N"}, sel), {cartoon:{color:"#30363d",opacity:0.10}});
}
function stylePeptide(v, color, model) {
  var sel = model !== undefined ? {model:model, chain:"P"} : {chain:"P"};
  v.setStyle(sel, {stick:{radius:0.16,color:color}, cartoon:{color:color,opacity:0.65,tube:true}});
  v.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.28,color:color}, sel);
}
function stylePlddt(v, model) {
  /* Color each residue by B-factor (pLDDT) using AlphaFold-style scheme:
     blue(>0.9) -> cyan(>0.7) -> yellow(>0.5) -> red(<0.5) */
  var sel = model !== undefined ? {model:model, chain:"P"} : {chain:"P"};
  var atoms = v.getModel(model !== undefined ? model : 0).selectedAtoms({chain:"P"});
  var seen = {};
  for (var i = 0; i < atoms.length; i++) {
    var resi = atoms[i].resi;
    if (seen[resi]) { atoms[i].b = seen[resi]; continue; }
    /* CA b-factor is the pLDDT; copy it so all atoms in the residue share it */
    if (atoms[i].atom === "CA") seen[resi] = atoms[i].b;
  }
  v.setStyle(sel, {
    stick:{radius:0.22, colorscheme:{prop:"b",gradient:new $3Dmol.Gradient.ROYGB(0.5,1.0)}},
    cartoon:{style:"tube",opacity:0.65, colorscheme:{prop:"b",gradient:new $3Dmol.Gradient.ROYGB(0.5,1.0)}}
  });
  v.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.55,
    colorscheme:{prop:"b",gradient:new $3Dmol.Gradient.ROYGB(0.5,1.0)}}, sel);
}
function styleCharge(v, model) {
  /* Color each peptide residue by amino acid electrostatic character.
     Red = negative (D,E), Blue = positive (K,R,H), White = neutral/hydrophobic.
     Use per-residue coloring via setStyle with a custom colorMap. */
  var sel = model !== undefined ? {model:model, chain:"P"} : {chain:"P"};
  var chargeColor = {
    ASP:"#ff4444", GLU:"#ff4444",                    /* negative: red */
    LYS:"#4466ff", ARG:"#4466ff", HIS:"#7799ff",     /* positive: blue */
    SER:"#e8e8e8", THR:"#e8e8e8", ASN:"#dde0e0", GLN:"#dde0e0",  /* polar: light grey */
    CYS:"#cccc88", TYR:"#bbcc99", TRP:"#aacc88",     /* special */
    GLY:"#ffffff", ALA:"#f5f5f5", PRO:"#eeddcc",     /* small/neutral: white-ish */
    VAL:"#ffe8cc", ILE:"#ffddaa", LEU:"#ffddaa", MET:"#ffcc88",  /* hydrophobic: warm */
    PHE:"#ffcc77"                                     /* aromatic hydrophobic */
  };
  /* Apply per-residue color via individual setStyle calls */
  var mdl = v.getModel(model !== undefined ? model : 0);
  var atoms = mdl.selectedAtoms({chain:"P"});
  var residues = {};
  for (var i = 0; i < atoms.length; i++) {
    var key = atoms[i].resi;
    if (!residues[key]) residues[key] = atoms[i].resn;
  }
  /* First set all peptide to white stick as base */
  v.setStyle(sel, {stick:{radius:0.16, color:"#ffffff"}});
  /* Then color each residue */
  for (var resi in residues) {
    var resn = residues[resi];
    var col = chargeColor[resn] || "#ffffff";
    var rsel = Object.assign({}, sel, {resi: parseInt(resi)});
    v.setStyle(rsel, {stick:{radius:0.16, color:col}});
  }
  /* Surface colored per-atom by the residue color */
  v.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.88,
    colorfunc: function(atom) {
      var c = chargeColor[atom.resn] || "#ffffff";
      return $3Dmol.CC.color(c);
    }
  }, sel);
}
function markMismatches(v, seq, model) {
  for (var i = 0; i < seq.length && i < TGT.length; i++) {
    if (seq[i] !== TGT[i]) {
      var sel = {chain:"P", resi:i+1};
      if (model !== undefined) sel.model = model;
      v.setStyle(sel, {stick:{radius:0.28, color:"#ff7b72"}});
    }
  }
}

function renderSide(pdbT, pdbD, seq) {
  vL.removeAllModels(); vL.removeAllSurfaces();
  vL.addModel(pdbT,"pdb"); styleMHC(vL); stylePeptide(vL,"#7ee787");
  vL.zoomTo({chain:"P"}); vL.render();

  vR.removeAllModels(); vR.removeAllSurfaces();
  vR.addModel(pdbD,"pdb"); styleMHC(vR); stylePeptide(vR,"#ffa657");
  markMismatches(vR, seq); vR.zoomTo({chain:"P"}); vR.render();
}

function renderSuper(pdbT, pdbD, seq) {
  vL.removeAllModels(); vL.removeAllSurfaces();
  vL.addModel(pdbT,"pdb"); styleMHC(vL,0); stylePeptide(vL,"#7ee787",0);
  vL.addModel(pdbD,"pdb");
  vL.setStyle({model:1,chain:"M"},{}); vL.setStyle({model:1,chain:"N"},{});
  stylePeptide(vL,"#ffa657",1); markMismatches(vL, seq, 1);
  vL.zoomTo({chain:"P"}); vL.render();
}

function hideMHC(v, model) {
  /* Hide MHC chains entirely — used in pLDDT/Surface modes where MHC is distracting */
  var m = model !== undefined ? {model:model} : {};
  v.setStyle(Object.assign({chain:"M"}, m), {});
  v.setStyle(Object.assign({chain:"N"}, m), {});
}

function renderPlddt(pdbT, pdbD, seq) {
  vL.removeAllModels(); vL.removeAllSurfaces();
  vL.addModel(pdbT,"pdb"); hideMHC(vL,0); stylePlddt(vL,0);
  vL.zoomTo({chain:"P"}); vL.render();

  vR.removeAllModels(); vR.removeAllSurfaces();
  vR.addModel(pdbD,"pdb"); hideMHC(vR,0); stylePlddt(vR,0);
  vR.zoomTo({chain:"P"}); vR.render();
}

function renderSurface(pdbT, pdbD, seq) {
  vL.removeAllModels(); vL.removeAllSurfaces();
  vL.addModel(pdbT,"pdb"); hideMHC(vL,0); styleCharge(vL,0);
  vL.zoomTo({chain:"P"}); vL.render();

  vR.removeAllModels(); vR.removeAllSurfaces();
  vR.addModel(pdbD,"pdb"); hideMHC(vR,0); styleCharge(vR,0);
  vR.zoomTo({chain:"P"}); vR.render();
}

/* ── Heatmap ── */
function buildHeatmap(seq, decoyBF) {
  var panel = document.getElementById("heatmap-panel");
  panel.style.display = "block";
  document.getElementById("heatmap-subtitle").textContent = TGT + "  vs  " + seq;

  var lbls = document.getElementById("hm-labels"); lbls.innerHTML = "";
  var rowT = document.getElementById("hm-target"); rowT.innerHTML = "";
  var rowD = document.getElementById("hm-decoy"); rowD.innerHTML = "";
  var pT = document.getElementById("hm-plddt-tgt"); pT.innerHTML = "";
  var pD = document.getElementById("hm-plddt-dec"); pD.innerHTML = "";

  var maxLen = Math.max(TGT.length, seq.length);
  for (var i = 0; i < maxLen; i++) {
    // Position label
    var lbl = document.createElement("div");
    lbl.className = "res-cell header";
    lbl.textContent = "p" + (i+1);
    lbls.appendChild(lbl);

    // Target residue
    var tAA = i < TGT.length ? TGT[i] : "-";
    var dAA = i < seq.length ? seq[i] : "-";
    var match = tAA === dAA;

    var tCell = document.createElement("div");
    tCell.className = "res-cell";
    tCell.textContent = tAA;
    tCell.style.background = match ? "rgba(126,231,135,0.15)" : "rgba(126,231,135,0.08)";
    tCell.style.color = "var(--accent-green)";
    tCell.title = "Target: " + tAA + " | Hydro: " + (AA_HYDRO[tAA]||0).toFixed(1) + " | Charge: " + (AA_CHARGE[tAA]||0);
    rowT.appendChild(tCell);

    var dCell = document.createElement("div");
    dCell.className = "res-cell";
    dCell.textContent = dAA;
    dCell.style.color = match ? "var(--accent-green)" : "var(--accent-red)";
    dCell.style.background = match ? "rgba(126,231,135,0.15)" : "rgba(255,123,114,0.15)";
    dCell.title = "Decoy: " + dAA + " | Hydro: " + (AA_HYDRO[dAA]||0).toFixed(1) + " | Charge: " + (AA_CHARGE[dAA]||0);
    rowD.appendChild(dCell);

    // pLDDT bars
    var tPlddt = (TGT_BF[i] || {}).plddt || 0;
    var dPlddt = (decoyBF[i] || {}).plddt || 0;

    var tBar = document.createElement("div");
    tBar.className = "res-cell";
    tBar.innerHTML = '<div class="plddt-bar" style="background:'+plddtColor(tPlddt)+';width:'+Math.round(tPlddt*100)+'%"></div>';
    tBar.title = "Target pLDDT: " + tPlddt.toFixed(2);
    tBar.style.background = "transparent";
    pT.appendChild(tBar);

    var dBar = document.createElement("div");
    dBar.className = "res-cell";
    dBar.innerHTML = '<div class="plddt-bar" style="background:'+plddtColor(dPlddt)+';width:'+Math.round(dPlddt*100)+'%"></div>';
    dBar.title = "Decoy pLDDT: " + dPlddt.toFixed(2);
    dBar.style.background = "transparent";
    pD.appendChild(dBar);
  }
}

/* ── Color scale legend ── */
function showScale(side, title, gradientCSS, lo, hi) {
  var el = document.getElementById("cs-" + side);
  el.style.display = "block";
  document.getElementById("cs-" + side + "-title").textContent = title;
  document.getElementById("cs-" + side + "-bar").style.background = gradientCSS;
  document.getElementById("cs-" + side + "-lo").textContent = lo;
  document.getElementById("cs-" + side + "-hi").textContent = hi;
}
function hideScales() {
  document.getElementById("cs-left").style.display = "none";
  document.getElementById("cs-right").style.display = "none";
}
var PLDDT_GRAD = "linear-gradient(to right, #ff0000, #ffaa00, #ffff00, #00ff00, #0000ff)";
var CHARGE_GRAD = "linear-gradient(to right, #ff3333, #ff8888, #ffffff, #8888ff, #3333ff)";

/* ── Mode switching ── */
function setMode(m) {
  mode = m;
  document.querySelectorAll(".mode-btn").forEach(function(b){b.classList.remove("active")});
  if (m === "side") document.getElementById("btn-sidebyside").classList.add("active");
  else if (m === "super") document.getElementById("btn-super").classList.add("active");
  else if (m === "plddt") document.getElementById("btn-plddt").classList.add("active");
  else if (m === "surface") document.getElementById("btn-surface").classList.add("active");

  var showRight = (m === "side" || m === "plddt" || m === "surface");
  document.getElementById("box-right").style.display = showRight ? "block" : "none";
  document.getElementById("divider").style.display = showRight ? "block" : "none";

  /* Switch background: dark grey for surface mode so the white/red/blue surface pops */
  var bg = (m === "surface") ? "#2b2b2b" : "#0d1117";
  if (vL) vL.setBackgroundColor(bg);
  if (vR) vR.setBackgroundColor(bg);

  if (curIdx >= 0) selectDecoy(curIdx);
}

/* ── Selection ── */
function selectDecoy(idx) {
  curIdx = idx;
  var d = DECOYS[idx];
  document.querySelectorAll(".card").forEach(function(c){c.classList.remove("active")});
  var ac = document.querySelector('.card[data-idx="'+idx+'"]');
  if (ac) { ac.classList.add("active"); ac.scrollIntoView({block:"nearest"}); }

  document.getElementById("i-src").textContent = d.source;
  document.getElementById("i-hd").textContent = d.hamming_distance;
  document.getElementById("i-phys").textContent = d.physchem_sim;
  document.getElementById("i-struct").textContent = d.struct_sim;
  document.getElementById("i-risk").textContent = d.risk_score;
  document.getElementById("i-gene").textContent = d.genes;
  document.getElementById("conf-container").innerHTML = confBadge(d.remarks);

  hideScales();
  if (mode === "super") {
    document.getElementById("label-left").innerHTML = '<span style="color:var(--accent-green)">'+TGT+'</span> vs <span style="color:var(--accent-orange)">'+d.sequence+'</span>';
    renderSuper(TGT_PDB, d.pdb, d.sequence);
  } else if (mode === "plddt") {
    document.getElementById("label-left").textContent = "Target pLDDT";
    document.getElementById("label-right").textContent = "Decoy pLDDT: " + d.sequence;
    showScale("left", "pLDDT Confidence", PLDDT_GRAD, "0.5 (low)", "1.0 (high)");
    showScale("right", "pLDDT Confidence", PLDDT_GRAD, "0.5 (low)", "1.0 (high)");
    renderPlddt(TGT_PDB, d.pdb, d.sequence);
  } else if (mode === "surface") {
    document.getElementById("label-left").textContent = "Target Surface";
    document.getElementById("label-right").textContent = "Decoy Surface: " + d.sequence;
    showScale("left", "Electrostatic Character", CHARGE_GRAD, "Negative (acidic)", "Positive (basic)");
    showScale("right", "Electrostatic Character", CHARGE_GRAD, "Negative (acidic)", "Positive (basic)");
    renderSurface(TGT_PDB, d.pdb, d.sequence);
  } else {
    document.getElementById("label-left").textContent = "Target: " + TGT;
    document.getElementById("label-right").textContent = "Decoy: " + d.sequence;
    renderSide(TGT_PDB, d.pdb, d.sequence);
  }

  buildHeatmap(d.sequence, d.bfactors || []);
}

/* ── Init ── */
window.addEventListener("load", function() {
  vL = $3Dmol.createViewer("viewer-left", {backgroundColor:"#0d1117"});
  vR = $3Dmol.createViewer("viewer-right", {backgroundColor:"#0d1117"});
  vL.addModel(TGT_PDB,"pdb"); styleMHC(vL); stylePeptide(vL,"#7ee787");
  vL.zoomTo({chain:"P"}); vL.render();
  buildCards();
  if (DECOYS.length > 0) selectDecoy(0);
});

/* Keyboard nav */
document.addEventListener("keydown", function(e) {
  if (e.key === "ArrowDown" && curIdx < DECOYS.length - 1) { selectDecoy(curIdx+1); e.preventDefault(); }
  if (e.key === "ArrowUp" && curIdx > 0) { selectDecoy(curIdx-1); e.preventDefault(); }
  if (e.key === "1") setMode("side");
  if (e.key === "2") setMode("super");
  if (e.key === "3") setMode("plddt");
  if (e.key === "4") setMode("surface");
});
</script>
</body>
</html>'''


# ═══════════════════════════════════════════════════════════════════════
#  2. Decoy A Overview (Bubble Chart + Alignment)
# ═══════════════════════════════════════════════════════════════════════

def build_decoy_a_overview(base: Path, target_seq: str) -> Path:
    results_path = base / "Decoy_A" / "decoy_a_results.json"
    if not results_path.exists():
        print("  [WARN] Decoy A results not found, skipping")
        return None

    data = json.loads(results_path.read_text(encoding="utf-8"))
    if not data:
        print("  [WARN] Decoy A results empty, skipping")
        return None

    # Prepare entries for JS
    entries = []
    for e in data:
        expr = e.get("expression") or {}
        entries.append({
            "seq": e["sequence"],
            "hd": e["hamming_distance"],
            "el": round(e["el_rank"], 4),
            "risk": round(e.get("similarity_score", 0), 4),
            "genes": e.get("gene_symbols", [])[:3],
            "tpm": round(expr.get("max_vital_organ_tpm", 0), 1) if expr else 0,
            "vital": expr.get("vital_organ_expressed", False) if expr else False,
            "cat": expr.get("expression_category", "unknown") if expr else "unknown",
            "mm": e.get("mismatches", []),
            "tcr_mm": e.get("n_tcr_contact_mismatches", 0),
            "anc_mm": e.get("n_anchor_mismatches", 0),
        })

    entries_json = json.dumps(entries, ensure_ascii=False)

    html = _build_decoy_a_html(target_seq, entries_json, len(entries))
    out_path = base / "decoy_a_overview.html"
    out_path.write_text(html, encoding="utf-8")
    print(f"  [OK] Decoy A overview: {out_path} ({len(entries)} hits, {len(html):,} bytes)")
    return out_path


def _build_decoy_a_html(target_seq, entries_json, n_entries):
    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Decoy A Overview — {target_seq}</title>
<style>
:root {{
  --bg0:#0d1117; --bg1:#161b22; --bg2:#21262d; --bd:#30363d;
  --t1:#f0f6fc; --t2:#c9d1d9; --tm:#8b949e;
  --green:#7ee787; --orange:#ffa657; --red:#ff7b72; --blue:#58a6ff; --purple:#d2a8ff; --yellow:#e3b341;
}}
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:"Segoe UI",system-ui,sans-serif;background:var(--bg0);color:var(--t2)}}
.header{{background:var(--bg1);padding:16px 28px;border-bottom:1px solid var(--bd)}}
.header h1{{font-size:18px;color:var(--blue)}}
.header p{{font-size:12px;color:var(--tm);margin-top:3px}}
.main{{display:flex;height:calc(100vh - 64px)}}
.chart-area{{flex:1;padding:20px;position:relative;overflow:hidden}}
.detail-panel{{width:380px;background:var(--bg1);border-left:1px solid var(--bd);overflow-y:auto;display:flex;flex-direction:column}}
.detail-panel::-webkit-scrollbar{{width:5px}}
.detail-panel::-webkit-scrollbar-thumb{{background:var(--bd);border-radius:3px}}

/* Chart */
svg {{font-family:"Segoe UI",system-ui,sans-serif}}
.axis text {{fill:var(--tm);font-size:11px}}
.axis line,.axis path {{stroke:var(--bd)}}
.grid line {{stroke:var(--bd);stroke-opacity:0.3;stroke-dasharray:2,3}}
.bubble {{cursor:pointer;transition:opacity .15s}}
.bubble:hover {{opacity:1!important}}
.tooltip {{position:absolute;background:rgba(22,27,34,0.95);backdrop-filter:blur(8px);border:1px solid var(--bd);border-radius:8px;padding:10px 14px;font-size:12px;pointer-events:none;z-index:100;display:none;max-width:280px}}
.tooltip .seq{{font-family:Consolas,monospace;font-size:14px;font-weight:700;letter-spacing:1px}}
.tooltip .m{{color:var(--green)}} .tooltip .x{{color:var(--red)}}

/* Legend */
.chart-legend{{position:absolute;top:20px;right:20px;background:rgba(22,27,34,0.85);border:1px solid var(--bd);border-radius:8px;padding:12px 16px;font-size:11px;color:var(--tm)}}
.chart-legend .row{{display:flex;align-items:center;gap:6px;margin:3px 0}}
.chart-legend .circ{{width:10px;height:10px;border-radius:50%;flex-shrink:0}}

/* Detail panel */
.detail-header{{padding:14px 16px;border-bottom:1px solid var(--bd);background:var(--bg1)}}
.detail-header .seq{{font-family:Consolas,monospace;font-size:20px;font-weight:700;letter-spacing:2px}}
.detail-header .seq .m{{color:var(--green)}} .detail-header .seq .x{{color:var(--red)}}
.detail-body{{padding:16px;flex:1}}
.stat-grid{{display:grid;grid-template-columns:1fr 1fr;gap:8px;margin-bottom:16px}}
.stat{{background:var(--bg2);border-radius:6px;padding:10px 12px}}
.stat .label{{font-size:10px;color:var(--tm);text-transform:uppercase;letter-spacing:.5px}}
.stat .value{{font-size:18px;font-weight:700;color:var(--t1);margin-top:2px}}

/* Alignment viz */
.align-title{{font-size:11px;color:var(--tm);margin-bottom:8px;font-weight:500}}
.align-row{{display:flex;gap:2px;margin-bottom:3px}}
.align-cell{{width:36px;height:36px;display:flex;align-items:center;justify-content:center;border-radius:5px;font-family:Consolas,monospace;font-size:14px;font-weight:700;position:relative}}
.align-cell .pos{{position:absolute;top:-12px;font-size:8px;color:var(--tm);font-weight:400}}
.align-cell.match{{background:rgba(126,231,135,0.12);color:var(--green)}}
.align-cell.mismatch{{background:rgba(255,123,114,0.12);color:var(--red)}}
.align-cell.anchor{{border:1.5px solid var(--purple)}}
.align-cell.tcr{{border:1.5px solid var(--yellow)}}
.mm-tag{{display:inline-block;padding:1px 6px;border-radius:8px;font-size:9px;font-weight:600;margin:2px}}
.mm-tag.tcr{{background:rgba(227,179,65,0.15);color:var(--yellow)}}
.mm-tag.anc{{background:rgba(210,168,255,0.15);color:var(--purple)}}

/* Expression bar */
.expr-section{{margin-top:16px}}
.expr-bar-bg{{height:8px;background:var(--bg2);border-radius:4px;overflow:hidden;margin-top:4px}}
.expr-bar-fill{{height:100%;border-radius:4px;transition:width .3s}}
</style>
</head>
<body>
<div class="header">
  <h1>Decoy A — Sequence Homology Overview</h1>
  <p>Target: <b style="color:var(--green);font-family:monospace;letter-spacing:1px">{target_seq}</b> &bull; HLA-A*02:01 &bull; {n_entries} hits (Hamming &le; 4) &bull; Bubble size = risk score &bull; Color = vital organ expression</p>
</div>
<div class="main">
  <div class="chart-area" id="chart-area">
    <div class="chart-legend">
      <div style="font-weight:600;margin-bottom:4px">Expression Risk</div>
      <div class="row"><div class="circ" style="background:#ff7b72"></div>High (vital organ TPM &gt; 10)</div>
      <div class="row"><div class="circ" style="background:#e3b341"></div>Medium (TPM 1-10)</div>
      <div class="row"><div class="circ" style="background:#58a6ff"></div>Low / Unknown</div>
      <div style="margin-top:8px;font-weight:600;margin-bottom:4px">Border</div>
      <div class="row"><div class="circ" style="border:2px solid var(--yellow);background:transparent"></div>TCR contact mismatch</div>
      <div class="row"><div class="circ" style="border:2px solid var(--purple);background:transparent"></div>Anchor mismatch</div>
    </div>
    <div class="tooltip" id="tooltip"></div>
  </div>
  <div class="detail-panel" id="detail-panel">
    <div class="detail-header" id="detail-header">
      <div style="color:var(--tm);font-size:12px;margin-bottom:4px">Click a bubble to inspect</div>
      <div class="seq" id="detail-seq" style="color:var(--tm)">—</div>
    </div>
    <div class="detail-body" id="detail-body">
      <div class="stat-grid">
        <div class="stat"><div class="label">Hamming</div><div class="value" id="d-hd">—</div></div>
        <div class="stat"><div class="label">EL%Rank</div><div class="value" id="d-el">—</div></div>
        <div class="stat"><div class="label">TCR Mis</div><div class="value" id="d-tcr">—</div></div>
        <div class="stat"><div class="label">Anc Mis</div><div class="value" id="d-anc">—</div></div>
      </div>
      <div id="align-container"></div>
      <div class="expr-section" id="expr-section" style="display:none">
        <div class="align-title">Vital Organ Expression</div>
        <div style="display:flex;justify-content:space-between;font-size:12px">
          <span id="d-gene" style="color:var(--t1);font-weight:600"></span>
          <span id="d-tpm" style="color:var(--tm)"></span>
        </div>
        <div class="expr-bar-bg"><div class="expr-bar-fill" id="d-tpm-bar"></div></div>
      </div>
    </div>
  </div>
</div>

<script>
var TGT = "{target_seq}";
var DATA = {entries_json};

var W, H, svg;
var margin = {{top:40,right:40,bottom:50,left:60}};

function colorSeq(seq) {{
  return seq.split("").map(function(c,i){{
    return '<span class="'+(i<TGT.length&&c===TGT[i]?"m":"x")+'">'+c+'</span>';
  }}).join("");
}}

function exprColor(d) {{
  if (d.tpm > 10) return "#ff7b72";
  if (d.tpm > 1) return "#e3b341";
  return "#58a6ff";
}}

function drawChart() {{
  var area = document.getElementById("chart-area");
  var rect = area.getBoundingClientRect();
  W = rect.width; H = rect.height;
  var iw = W - margin.left - margin.right;
  var ih = H - margin.top - margin.bottom;

  // Clear
  area.querySelectorAll("svg").forEach(function(s){{s.remove()}});

  var ns = "http://www.w3.org/2000/svg";
  svg = document.createElementNS(ns, "svg");
  svg.setAttribute("width", W); svg.setAttribute("height", H);
  area.insertBefore(svg, area.firstChild);

  // Scales
  var hdMin = 0, hdMax = 5;
  var elMax = Math.min(2.0, Math.max.apply(null, DATA.map(function(d){{return d.el}})) * 1.1);
  var riskMax = Math.max.apply(null, DATA.map(function(d){{return d.risk}}));

  function sx(hd) {{ return margin.left + (hd - hdMin) / (hdMax - hdMin) * iw; }}
  function sy(el) {{ return margin.top + (1 - el / elMax) * ih; }}
  function sr(risk) {{ return 4 + (risk / (riskMax||1)) * 18; }}

  // Grid
  for (var i = 0; i <= 5; i++) {{
    var x = sx(i);
    var line = document.createElementNS(ns,"line");
    line.setAttribute("x1",x); line.setAttribute("x2",x);
    line.setAttribute("y1",margin.top); line.setAttribute("y2",margin.top+ih);
    line.setAttribute("class","grid"); svg.appendChild(line);
  }}

  // Axes
  // X axis
  for (var i = 0; i <= 5; i++) {{
    var t = document.createElementNS(ns,"text");
    t.setAttribute("x",sx(i)); t.setAttribute("y",margin.top+ih+28);
    t.setAttribute("text-anchor","middle"); t.setAttribute("class","axis");
    t.textContent = i; svg.appendChild(t);
  }}
  var xlabel = document.createElementNS(ns,"text");
  xlabel.setAttribute("x",margin.left+iw/2); xlabel.setAttribute("y",margin.top+ih+46);
  xlabel.setAttribute("text-anchor","middle"); xlabel.setAttribute("fill","#8b949e");
  xlabel.setAttribute("font-size","12"); xlabel.textContent = "Hamming Distance";
  svg.appendChild(xlabel);

  // Y axis
  var nTicks = 5;
  for (var i = 0; i <= nTicks; i++) {{
    var val = (elMax / nTicks * i);
    var y = sy(val);
    var t = document.createElementNS(ns,"text");
    t.setAttribute("x",margin.left-10); t.setAttribute("y",y+4);
    t.setAttribute("text-anchor","end"); t.setAttribute("class","axis");
    t.textContent = val.toFixed(2); svg.appendChild(t);
    var gl = document.createElementNS(ns,"line");
    gl.setAttribute("x1",margin.left); gl.setAttribute("x2",margin.left+iw);
    gl.setAttribute("y1",y); gl.setAttribute("y2",y);
    gl.setAttribute("class","grid"); svg.appendChild(gl);
  }}
  var ylabel = document.createElementNS(ns,"text");
  ylabel.setAttribute("transform","rotate(-90)");
  ylabel.setAttribute("x",-(margin.top+ih/2)); ylabel.setAttribute("y",16);
  ylabel.setAttribute("text-anchor","middle"); ylabel.setAttribute("fill","#8b949e");
  ylabel.setAttribute("font-size","12"); ylabel.textContent = "EL%Rank (lower = stronger binder)";
  svg.appendChild(ylabel);

  // Jitter for overlapping points
  var jitter = {{}};
  DATA.forEach(function(d,i) {{
    var key = d.hd + "_" + d.el.toFixed(3);
    if (!jitter[key]) jitter[key] = 0;
    jitter[key]++;
    d._jx = (jitter[key] - 1) * 3 - (jitter[key] > 1 ? jitter[key]*1.5 : 0);
    d._jy = (jitter[key] % 2 === 0 ? 1 : -1) * jitter[key] * 1.5;
    d._idx = i;
  }});

  // Bubbles
  DATA.forEach(function(d) {{
    var g = document.createElementNS(ns,"g");
    g.setAttribute("class","bubble");
    g.style.opacity = "0.75";

    var cx = sx(d.hd) + (d._jx || 0);
    var cy = sy(d.el) + (d._jy || 0);
    var r = sr(d.risk);

    var circle = document.createElementNS(ns,"circle");
    circle.setAttribute("cx", cx); circle.setAttribute("cy", cy);
    circle.setAttribute("r", r);
    circle.setAttribute("fill", exprColor(d));
    circle.setAttribute("fill-opacity", "0.6");
    circle.setAttribute("stroke", d.tcr_mm > 0 ? "#e3b341" : (d.anc_mm > 0 ? "#d2a8ff" : "rgba(255,255,255,0.15)"));
    circle.setAttribute("stroke-width", d.tcr_mm > 0 || d.anc_mm > 0 ? "2" : "1");

    g.appendChild(circle);

    g.onmouseenter = function(ev) {{
      g.style.opacity = "1";
      var tip = document.getElementById("tooltip");
      tip.style.display = "block";
      tip.innerHTML = '<div class="seq">' + colorSeq(d.seq) + '</div>'
        + '<div style="margin-top:4px;color:#8b949e">HD=' + d.hd + ' | EL=' + d.el.toFixed(4)
        + ' | ' + (d.genes.join(", ") || "—") + '</div>';
      tip.style.left = (ev.clientX + 14) + "px";
      tip.style.top = (ev.clientY - 10) + "px";
    }};
    g.onmouseleave = function() {{
      g.style.opacity = "0.75";
      document.getElementById("tooltip").style.display = "none";
    }};
    g.onclick = function() {{ selectHit(d._idx); }};

    svg.appendChild(g);
  }});
}}

function selectHit(idx) {{
  var d = DATA[idx];
  document.getElementById("detail-seq").innerHTML = colorSeq(d.seq);
  document.getElementById("d-hd").textContent = d.hd;
  document.getElementById("d-el").textContent = d.el.toFixed(4);
  document.getElementById("d-tcr").textContent = d.tcr_mm;
  document.getElementById("d-anc").textContent = d.anc_mm;

  // Alignment
  var ac = document.getElementById("align-container");
  ac.innerHTML = '<div class="align-title">Position-level Alignment</div>';
  var rowT = document.createElement("div"); rowT.className = "align-row";
  var rowD = document.createElement("div"); rowD.className = "align-row";

  for (var i = 0; i < TGT.length; i++) {{
    var tAA = TGT[i];
    var dAA = i < d.seq.length ? d.seq[i] : "-";
    var match = tAA === dAA;

    // Find mismatch info
    var mmInfo = null;
    if (d.mm) d.mm.forEach(function(m){{ if (m.position === i) mmInfo = m; }});

    var tCell = document.createElement("div");
    tCell.className = "align-cell match";
    tCell.innerHTML = '<span class="pos">p'+(i+1)+'</span>' + tAA;
    rowT.appendChild(tCell);

    var dCell = document.createElement("div");
    dCell.className = "align-cell " + (match ? "match" : "mismatch");
    if (mmInfo && mmInfo.is_anchor) dCell.classList.add("anchor");
    if (mmInfo && mmInfo.is_tcr_contact) dCell.classList.add("tcr");
    dCell.innerHTML = dAA;
    rowD.appendChild(dCell);
  }}

  ac.appendChild(rowT);
  ac.appendChild(rowD);

  // Mismatch tags
  if (d.tcr_mm > 0 || d.anc_mm > 0) {{
    var tags = document.createElement("div");
    tags.style.marginTop = "6px";
    if (d.tcr_mm > 0) tags.innerHTML += '<span class="mm-tag tcr">' + d.tcr_mm + ' TCR contact</span>';
    if (d.anc_mm > 0) tags.innerHTML += '<span class="mm-tag anc">' + d.anc_mm + ' anchor</span>';
    ac.appendChild(tags);
  }}

  // Expression
  var es = document.getElementById("expr-section");
  if (d.tpm > 0) {{
    es.style.display = "block";
    document.getElementById("d-gene").textContent = d.genes.join(", ") || "—";
    document.getElementById("d-tpm").textContent = d.tpm + " TPM";
    var bar = document.getElementById("d-tpm-bar");
    var pct = Math.min(100, d.tpm / 100 * 100);
    bar.style.width = pct + "%";
    bar.style.background = exprColor(d);
  }} else {{
    es.style.display = "none";
  }}
}}

window.addEventListener("load", function() {{
  drawChart();
  if (DATA.length > 0) selectHit(0);
}});
window.addEventListener("resize", drawChart);
</script>
</body>
</html>'''


# ═══════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Build Decoy Library visualizations")
    parser.add_argument("--target", default="GILGFVFTL", help="Target peptide sequence")
    parser.add_argument("--dir", default=None, help="Summary directory (default: data/<target>_summary)")
    args = parser.parse_args()

    target = args.target.strip().upper()
    base = Path(args.dir) if args.dir else Path(f"data/{target}_summary")

    if not base.exists():
        print(f"Error: directory {base} does not exist")
        return

    print(f"Building visualizations for {target} in {base}/")
    print()

    # Build 3D viewer (Decoy B/D)
    if (base / "Decoy_B").exists() or (base / "Decoy_D").exists():
        build_3d_viewer(base, target)
    else:
        print("  [SKIP] No Decoy B/D data found")

    # Build Decoy A overview
    if (base / "Decoy_A").exists():
        build_decoy_a_overview(base, target)
    else:
        print("  [SKIP] No Decoy A data found")

    print()
    print("Done. Open the HTML files in a browser.")


if __name__ == "__main__":
    main()
