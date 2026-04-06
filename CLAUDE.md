# Decoy Library - Project Guidelines

## Core Logic: Target Peptide Selection

**The decoy library pipeline takes a therapeutically promising target peptide as INPUT, then automatically screens for cross-reactive off-target risks around it.**

When proposing new peptide sequences for the pipeline:
- INPUT peptides should be **benign, therapeutically valuable** targets (e.g., viral epitopes, cancer-testis antigens with restricted normal tissue expression)
- The pipeline (Decoy A-D) **automatically discovers** dangerous off-target peptides nearby
- Do NOT select known fatal/toxic peptides as input targets — those belong in the Decoy C curated database as known risks, not as pipeline inputs
- Think: "this target looks safe, but what hidden cross-reactivity risks does the pipeline find around it?"

**Correct example**: GILGFVFTL (influenza viral epitope) as input -> pipeline finds similar human-proteome peptides that could cause cross-reactivity
**Wrong example**: FLWGPRALV (MAGE-A3, known fatal cardiac toxicity) as input target — this is itself a known danger, not a therapeutic candidate to screen

## Documentation Update Rule

**Every time the codebase is modified, the following documentation files MUST be updated:**

1. **`README.md`** (root) — Record only **important, high-level information**: architecture changes, new tools/descriptors, major pipeline updates. Keep it concise.
2. **`progress_and_report.md`** (root) — Record **detailed technical information**: implementation specifics, algorithm changes, deployment steps, test results, parameter tuning.
3. **Individual decoy folder READMEs** (`decoy_a/README.md`, `decoy_b/README.md`, `decoy_c/README.md`, `decoy_d/README.md`) — Update the README of any decoy module whose code was changed.

**Hierarchy**: README.md = executive summary | progress_and_report.md = full technical detail | decoy_*/README.md = module-specific documentation
