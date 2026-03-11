#!/usr/bin/env bash
# =============================================================================
# Inference on silicus telescopic variant genomes using mm10-trained flashzoi
#
# Telescopic series: progressively add features to silicus background.
# silicus → +CGD → +CGD+CRE → +CGD+CRE+CTCF → +CGD+CRE+CTCF+Exon → +CGD+CRE+CTCF+Exon+TE
#
# The base silicus and +CGD variants are covered by run_silicus_plus_inference.sh.
# This script handles the combined variants (CGD+CRE, CGD+CRE+CTCF, etc.).
#
# Usage:
#   ./run_silicus_telescopic_inference.sh
#   ./run_silicus_telescopic_inference.sh silicusPlusCGDCre
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BORZOI_CODE="${BORZOI_CODE:-$(cd "$SCRIPT_DIR/../code" && pwd)}"

CONFIG="${SCRIPT_DIR}/../configs/flashzoi_rf/flashzoi_rf524k.yaml"
PREDICTIONS_BASE="${PREDICTIONS_BASE:-predictions}"
MISHA_GROOT="${MISHA_GROOT:-misha_db/mm10}"
GENOMES_DIR="${GENOMES_DIR:-genomes}"

export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0,1,2,3,4,5,6,7}"
NUM_GPUS="${NUM_GPUS:-8}"

declare -A VARIANTS=(
    ["silicusPlusCGDCre"]="${GENOMES_DIR}/silicus_plus/silicusPlusCGDCre/silicusPlusCGDCre.fa"
    ["silicusPlusCGDCreCtcf"]="${GENOMES_DIR}/silicus_plus/silicusPlusCGDCreCtcf/silicusPlusCGDCreCtcf.fa"
    ["silicusPlusCGDCreCtcfExon"]="${GENOMES_DIR}/silicus_plus/silicusPlusCGDCreCtcfExon/silicusPlusCGDCreCtcfExon.fa"
    ["silicusPlusCGDCreCtcfExonTE"]="${GENOMES_DIR}/silicus_plus/silicusPlusCGDCreCtcfExonTE/silicusPlusCGDCreCtcfExonTE.fa"
)

run_inference() {
    local variant="$1"
    local genome_fasta="${VARIANTS[$variant]}"
    local output_dir="${PREDICTIONS_BASE}/flashzoi_rf524k_${variant}"
    local track_prefix="seq.IQ.pcg.flashzoi.${variant}.rf524k_"

    echo "Running inference: ${variant}"
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
