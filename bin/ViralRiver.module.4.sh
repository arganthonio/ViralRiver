#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 4: Sequencing Depth Summary (Definitivo para Nextflow)
# ==============================================================================
# Extracts total RNA-seq reads from fastp JSON reports for downstream RPM
# normalization.
# ==============================================================================

set -euo pipefail

RESULTS_DIR=""
OUTPUT_FILE="viralriver_read_depth.tsv"

usage() {
    echo "Usage: $0 -i <results_dir> [-o <output_file>]"
    echo "  -i  ViralRiver results directory containing fastp json files"
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

# Buscamos directamente todos los archivos .json que Nextflow ha copiado en la carpeta
for json_path in "$RESULTS_DIR"/*_fastp.json; do
    # Si no hay archivos que coincidan, romper el bucle de forma segura
    [[ -f "$json_path" ]] || continue

    # Extraemos el nombre de la muestra quitando el sufijo '_fastp.json'
    filename=$(basename "$json_path")
    sample="${filename%_fastp.json}"

    # Extraemos las lecturas directamente del archivo JSON
    before_reads=$(grep -A 5 '"before_filtering"' "$json_path" | grep -m 1 -oP '"total_reads":\s*\K\d+' || echo "NA")
    after_reads=$(grep -A 5 '"after_filtering"' "$json_path" | grep -m 1 -oP '"total_reads":\s*\K\d+' || echo "NA")

    echo -e "${sample}\t${before_reads}\t${after_reads}" >> "$OUTPUT_FILE"
done

echo "Sequencing depth table generated: $OUTPUT_FILE"
