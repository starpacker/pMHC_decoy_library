#!/bin/bash
# Batch Boltz-2 predictions for all Decoy B candidates
# Skips already-completed predictions (non-empty result dirs)

set -e

BOLTZ_INPUT_DIR="/share/liuyutian/pMHC_decoy_library/data/decoy_b/pmhc_models/boltz/_boltz_inputs"
BOLTZ_OUTPUT_DIR="/share/liuyutian/pMHC_decoy_library/data/decoy_b/pmhc_models/boltz"
BOLTZ_CACHE="$HOME/.boltz"

TOTAL=$(ls "$BOLTZ_INPUT_DIR"/*.yaml 2>/dev/null | wc -l)
DONE=0
SKIP=0
FAIL=0

echo "=== Boltz-2 Batch Prediction ==="
echo "Total inputs: $TOTAL"
echo "Output dir: $BOLTZ_OUTPUT_DIR"
echo ""

for yaml_file in "$BOLTZ_INPUT_DIR"/*.yaml; do
    stem=$(basename "$yaml_file" .yaml)
    result_dir="$BOLTZ_OUTPUT_DIR/boltz_results_${stem}/predictions/${stem}"
    
    # Check if already completed (has PDB output)
    if [ -d "$result_dir" ] && ls "$result_dir"/*.pdb 1>/dev/null 2>&1; then
        SKIP=$((SKIP + 1))
        echo "[SKIP $SKIP] $stem — already completed"
        continue
    fi
    
    DONE=$((DONE + 1))
    echo ""
    echo "[RUN $DONE/$((TOTAL - SKIP))] $stem"
    echo "  Input: $yaml_file"
    
    # Clean up any partial results
    rm -rf "$BOLTZ_OUTPUT_DIR/boltz_results_${stem}"
    
    # Run Boltz-2
    if conda run -n boltz --no-capture-output boltz predict \
        "$yaml_file" \
        --out_dir "$BOLTZ_OUTPUT_DIR" \
        --cache "$BOLTZ_CACHE" \
        --accelerator gpu \
        --model boltz2 \
        --output_format pdb \
        --devices 1 \
        --recycling_steps 3 \
        --sampling_steps 200 \
        --diffusion_samples 1 \
        --no_kernels 2>&1 | tail -5; then
        echo "  ✓ Success"
    else
        FAIL=$((FAIL + 1))
        echo "  ✗ Failed"
    fi
done

echo ""
echo "=== Summary ==="
echo "Total: $TOTAL"
echo "Skipped (cached): $SKIP"
echo "Predicted: $DONE"
echo "Failed: $FAIL"
echo "Success: $((DONE - FAIL))"
