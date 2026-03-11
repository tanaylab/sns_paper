#!/usr/bin/env bash
# =============================================================================
# Inference on mm10-minus variant genomes using mm10-trained flashzoi rf524k
#
# mm10-minus: the mm10 genome with specific features replaced by silicus
# background sequence. Tests which features are necessary for prediction.
#
# Usage:
#   ./run_mm10_minus_inference.sh
#   ./run_mm10_minus_inference.sh mm10MinusCGD
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
    ["mm10MinusCGD"]="${GENOMES_DIR}/mm10_minus/mm10MinusCGD/mm10MinusCGD.fa"
    ["mm10MinusCGDpad1k"]="${GENOMES_DIR}/mm10_minus/mm10MinusCGDpad1k/mm10MinusCGDpad1k.fa"
    ["mm10MinusCGDpad2k"]="${GENOMES_DIR}/mm10_minus/mm10MinusCGDpad2k/mm10MinusCGDpad2k.fa"
    ["mm10MinusCTCF"]="${GENOMES_DIR}/mm10_minus/mm10MinusCTCF/mm10MinusCTCF.fa"
    ["mm10MinusCRE"]="${GENOMES_DIR}/mm10_minus/mm10MinusCRE/mm10MinusCRE.fa"
    ["mm10MinusExon"]="${GENOMES_DIR}/mm10_minus/mm10MinusExon/mm10MinusExon.fa"
    ["mm10MinusTE"]="${GENOMES_DIR}/mm10_minus/mm10MinusTE/mm10MinusTE.fa"
    ["mm10MinusRandom"]="${GENOMES_DIR}/mm10_minus/mm10MinusRandom/mm10MinusRandom.fa"
    ["mm10MinusNonCGDExon"]="${GENOMES_DIR}/mm10_minus/mm10MinusNonCGDExon/mm10MinusNonCGDExon.fa"
    ["mm10MinusNonCGDCre"]="${GENOMES_DIR}/mm10_minus/mm10MinusNonCGDCre/mm10MinusNonCGDCre.fa"
    ["mm10MinusNonCGDTE"]="${GENOMES_DIR}/mm10_minus/mm10MinusNonCGDTE/mm10MinusNonCGDTE.fa"
    ["mm10MinusNonCGDCTCF"]="${GENOMES_DIR}/mm10_minus/mm10MinusNonCGDCTCF/mm10MinusNonCGDCTCF.fa"
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
