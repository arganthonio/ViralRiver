#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 4: Sequencing Depth Summary
# ==============================================================================
# Extracts total RNA-seq reads from fastp JSON reports for downstream RPM
# normalization.
# ==============================================================================

set -euo pipefail

RESULTS_DIR=""
OUTPUT_FILE="viralriver_read_depth.tsv"

usage() {
    echo "Usage: $0 -i <results_dir> [-o <output_file>]"
    echo "  -i  ViralRiver results directory containing sample subdirectories"
    echo "  -o  Output TSV file (default: viralriver_read_depth.tsv)"
    exit 1
}

while getopts "i:o:h" opt; do
    case "$opt" in
        i) RESULTS_DIR="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        h|?) usage ;;
    esac
done

if [[ -z "$RESULTS_DIR" ]]; then
    echo "ERROR: Missing required argument: -i <results_dir>"
    usage
fi

[[ -d "$RESULTS_DIR" ]] || {
    echo "ERROR: results directory not found: $RESULTS_DIR"
    exit 1
}

echo -e "sample\ttotal_reads_before_filtering\ttotal_reads_after_filtering" > "$OUTPUT_FILE"

for sample_dir in "$RESULTS_DIR"/*/; do
    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")

    json_path="${sample_dir}/fastp_reports/${sample}_fastp.json"

    if [[ ! -f "$json_path" ]]; then
        echo -e "${sample}\tNA\tNA" >> "$OUTPUT_FILE"
        continue
    fi

    before_reads=$(grep -A 5 '"before_filtering"' "$json_path" | grep -m 1 -oP '"total_reads":\s*\K\d+' || echo "NA")
    after_reads=$(grep -A 5 '"after_filtering"' "$json_path" | grep -m 1 -oP '"total_reads":\s*\K\d+' || echo "NA")

    echo -e "${sample}\t${before_reads}\t${after_reads}" >> "$OUTPUT_FILE"
done

echo "Sequencing depth table generated: $OUTPUT_FILE"