#!/usr/bin/env bash
#
# Batch run Decoy D on all candidate targets and generate detail figures.
#
# Usage:
#   bash scripts/batch_decoy_d.sh
#   bash scripts/batch_decoy_d.sh --designs 2000      # override MPNN design count
#   bash scripts/batch_decoy_d.sh --top-k 20          # override top-K structures
#
# Output:
#   data/decoy_d/{TARGET}/decoy_d_results.csv   (per target)
#   figures/{TARGET}_decoy_d_detail.png          (per target)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

export PYTHONUTF8=1

# ── Parse extra flags (forwarded to run_decoy.py) ────────────────────────
EXTRA_ARGS=("$@")

# ── Target list from candidate_targets.json ──────────────────────────────
# Each line: SEQUENCE HLA_ALLELE
TARGETS=$(python -c "
import json, pathlib
data = json.loads(pathlib.Path('data/candidate_targets.json').read_text(encoding='utf-8'))
for seq in data.get('existing_targets', []):
    print(seq, 'HLA-A*02:01')
for t in data.get('proposed_targets', []):
    print(t['sequence'], t.get('hla_allele', 'HLA-A*02:01'))
")

TOTAL=$(echo "$TARGETS" | wc -l | tr -d ' ')
IDX=0
FAILED=()

echo "=============================================="
echo " Batch Decoy D: $TOTAL targets"
echo "=============================================="

while IFS=' ' read -r SEQ HLA; do
    IDX=$((IDX + 1))
    echo ""
    echo "[$IDX/$TOTAL] $SEQ ($HLA)"
    echo "----------------------------------------------"

    if python run_decoy.py "$SEQ" d --hla "$HLA" "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"; then
        echo "[OK] $SEQ Decoy D completed"
    else
        echo "[FAIL] $SEQ Decoy D failed (exit $?)"
        FAILED+=("$SEQ")
    fi
done <<< "$TARGETS"

# ── Generate all figures ─────────────────────────────────────────────────
echo ""
echo "=============================================="
echo " Generating Decoy D detail figures..."
echo "=============================================="
python scripts/visualize_decoy_d_detail.py --all

# ── Summary ──────────────────────────────────────────────────────────────
echo ""
echo "=============================================="
echo " Batch Complete"
echo "  Total:  $TOTAL"
echo "  Failed: ${#FAILED[@]}"
if [ ${#FAILED[@]} -gt 0 ]; then
    echo "  Failed targets: ${FAILED[*]}"
fi
echo "  Figures: figures/*_decoy_d_detail.png"
echo "=============================================="
