#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 5: MultiQC Custom Quality Control Report
# ==============================================================================

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir>"
    echo "  -i  Input directory containing logs and reports to scan"
    echo "  -o  Output directory where the MultiQC report will be saved"
    exit 1
}

while getopts "i:o:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h|?) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "ERROR: Missing required arguments."
    usage
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Verifying MultiQC dependency..."
command -v multiqc >/dev/null || {
    echo "ERROR: MultiQC is not installed or not found in PATH."
    exit 1
}

log "Running Custom MultiQC on: $INPUT_DIR"

multiqc "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    --filename "viralriver_multiqc_report.html" \
    --config "multiqc_config.yaml" \
    --force

#sed -i 's/\.2f/\.0f/g' viralriver_multiqc_report_data/multiqc_data.json
#sed -i 's/,\.2f/,\.0f/g' viralriver_multiqc_report.html

log "=============================================================================="
log "MODULE 5 EXECUTION COMPLETED SUCCESSFULLY"
log "=============================================================================="
