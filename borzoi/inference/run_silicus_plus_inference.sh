#!/usr/bin/env bash
# =============================================================================
# Inference on silicus+ variant genomes using mm10-trained flashzoi rf524k model
#
# These are "inference-only" experiments: the model was trained on mm10
# (configs/flashzoi_rf/flashzoi_rf524k.yaml) and is evaluated on synthetic
# genomes where specific genomic features have been spliced into the
# silicus background.
#
# Usage:
#   ./run_silicus_plus_inference.sh                    # Run all variants
#   ./run_silicus_plus_inference.sh silicusPlusCGD      # Run specific variant
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"
BORZOI_CODE="${BORZOI_CODE:-$(cd "$SCRIPT_DIR/../code" && pwd)}"

# Configuration — update these paths for your setup
CONFIG="${SCRIPT_DIR}/../configs/flashzoi_rf/flashzoi_rf524k.yaml"
PREDICTIONS_BASE="${PREDICTIONS_BASE:-predictions}"
MISHA_GROOT="${MISHA_GROOT:-misha_db/mm10}"
GENOMES_DIR="${GENOMES_DIR:-genomes}"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0,1,2,3,4,5,6,7}"
NUM_GPUS="${NUM_GPUS:-8}"

declare -A VARIANTS=(
    ["silicus"]="${GENOMES_DIR}/silicus/mus_silicus_cg_gc_lower.fa"
    ["silicusPlusCGD"]="${GENOMES_DIR}/silicus_plus/silicusPlusCGD/silicusPlusCGD.fa"
    ["silicusPlusCTCF"]="${GENOMES_DIR}/silicus_plus/silicusPlusCTCF/silicusPlusCTCF.fa"
    ["silicusPlusCRE"]="${GENOMES_DIR}/silicus_plus/silicusPlusCRE/silicusPlusCRE.fa"
    ["silicusPlusExon"]="${GENOMES_DIR}/silicus_plus/silicusPlusExon/silicusPlusExon.fa"
    ["silicusPlusTE"]="${GENOMES_DIR}/silicus_plus/silicusPlusTE/silicusPlusTE.fa"
    ["silicusPlusRandom"]="${GENOMES_DIR}/silicus_plus/silicusPlusRandom/silicusPlusRandom.fa"
)

run_inference() {
    local variant="$1"
    local genome_fasta="${VARIANTS[$variant]}"
    local output_dir="${PREDICTIONS_BASE}/flashzoi_rf524k_${variant}"
    local track_prefix="seq.IQ.pcg.flashzoi.${variant}.rf524k_"

    echo "Running inference: ${variant}"
    echo "  Genome: ${genome_fasta}"
    echo "  Output: ${output_dir}"

    cd "$BORZOI_CODE"
    torchrun --nproc_per_node="$NUM_GPUS" \
        infer_borzoi_pytorch.py \
        --config "$CONFIG" \
        --genome_fasta "$genome_fasta" \
        --output_dir "$output_dir" \
        --misha_groot "$MISHA_GROOT" \
        --misha_track_prefix "$track_prefix" \
        --misha_suffix "" \
        --species mouse
}

if [[ $# -gt 0 ]]; then
    for variant in "$@"; do
        run_inference "$variant"
    done
else
    for variant in "${!VARIANTS[@]}"; do
        run_inference "$variant"
    done
fi
