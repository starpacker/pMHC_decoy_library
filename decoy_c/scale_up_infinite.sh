#!/usr/bin/env bash
# ═══════════════════════════════════════════════════════════════════════
# scale_up_infinite.sh — Launch scale_up.py in infinite mode
# ═══════════════════════════════════════════════════════════════════════
# Usage:
#   nohup bash scale_up_infinite.sh [TARGET] > scale_up_infinite.log 2>&1 &
#   nohup bash scale_up_infinite.sh 1000 > scale_up_infinite.log 2>&1 &
# ═══════════════════════════════════════════════════════════════════════

set -euo pipefail

TARGET="${1:-1000}"

cd /home/yjh/decoy_library

echo "════════════════════════════════════════════════════════════"
echo "  🚀 Launching scale_up.py in INFINITE mode"
echo "  Target: ${TARGET}"
echo "  Started: $(date)"
echo "════════════════════════════════════════════════════════════"

python scale_up.py \
    --required_num "$TARGET" \
    --strategy all \
    --multi-source \
    --batch-size 40 \
    --retmax 80 \
    --max-rounds 2000 \
    --save-every 3 \
    --no-quality-filter \
    --infinite \
    --resume

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  🏁 Process completed at $(date)"
echo "════════════════════════════════════════════════════════════"
