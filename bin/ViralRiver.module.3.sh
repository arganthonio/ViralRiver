#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 3: High-Sensitivity Viral Rescue with Bowtie2
# ==============================================================================

set -euo pipefail

THREADS=8
MAPQ=25
COMPLEXITY_THRESHOLD=30
VIRAL_DB_FASTA=""
INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -v <viral_db.fasta> [-t <threads>] [-q <mapq>] [-c <complexity_threshold>]"
    echo "  -i  Input base directory from Module 1"
    echo "  -o  Output base directory"
    echo "  -v  Curated viral reference FASTA"
    echo "  -t  Threads (default: 8)"
    echo "  -q  Minimum MAPQ for Bowtie2 viral alignments (default: 25)"
    echo "  -c  fastp low-complexity threshold (default: 30)"
    exit 1
}

while getopts "i:o:v:t:q:c:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        v) VIRAL_DB_FASTA="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MAPQ="$OPTARG" ;;
        c) COMPLEXITY_THRESHOLD="$OPTARG" ;;
        h|?) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$VIRAL_DB_FASTA" ]]; then
    echo "ERROR: Missing required arguments."
    usage
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Verifying software dependencies..."

for tool in fastp bowtie2 bowtie2-build samtools awk find sed sort; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Dependency '$tool' not found."
        exit 1
    }
done

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input directory not found: $INPUT_DIR"; exit 1; }
[[ -f "$VIRAL_DB_FASTA" ]] || { echo "ERROR: viral FASTA not found: $VIRAL_DB_FASTA"; exit 1; }

mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Bowtie2 index
# ------------------------------------------------------------------------------

REF_DIR=$(dirname "$VIRAL_DB_FASTA")
DB_NAME="human.virus.selected_idx"
BT2_PREFIX="${REF_DIR}/${DB_NAME}"

if [[ ! -f "${BT2_PREFIX}.1.bt2" && ! -f "${BT2_PREFIX}.1.bt2l" ]]; then
    log "Bowtie2 index not found. Building index: ${BT2_PREFIX}"
    bowtie2-build "$VIRAL_DB_FASTA" "$BT2_PREFIX" > /dev/null 2>&1
else
    log "Found existing Bowtie2 index: ${BT2_PREFIX}"
fi

# ------------------------------------------------------------------------------
# Sample discovery
# ------------------------------------------------------------------------------

log "Locating sample directories in ${INPUT_DIR}..."

mapfile -t SAMPLES < <(
    find "$INPUT_DIR" -mindepth 1 -maxdepth 1 -type d \
    | sed "s|${INPUT_DIR}/||" \
    | sort
)

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: No sample directories found in ${INPUT_DIR}"
    exit 1
fi

log "Found ${#SAMPLES[@]} sample directories."

# ------------------------------------------------------------------------------
# Core function
# ------------------------------------------------------------------------------

process_rescue() {
    local ID="$1"
    local IN_SUBDIR="${INPUT_DIR}/${ID}"
    local OUT_SUBDIR="${OUTPUT_DIR}/${ID}"

    mkdir -p "$OUT_SUBDIR"

    local R1=""
    local R2=""

    if [[ -f "${IN_SUBDIR}/${ID}_candidate_reads_1.fq.gz" ]]; then
        R1="${IN_SUBDIR}/${ID}_candidate_reads_1.fq.gz"
        R2="${IN_SUBDIR}/${ID}_candidate_reads_2.fq.gz"

    elif [[ -f "${IN_SUBDIR}/${ID}_extracted_viral_reads_1.fq.gz" ]]; then
        R1="${IN_SUBDIR}/${ID}_extracted_viral_reads_1.fq.gz"
        R2="${IN_SUBDIR}/${ID}_extracted_viral_reads_2.fq.gz"

    elif [[ -f "${IN_SUBDIR}/${ID}_KRAKEN_READS_1.fq.gz" ]]; then
        R1="${IN_SUBDIR}/${ID}_KRAKEN_READS_1.fq.gz"
        R2="${IN_SUBDIR}/${ID}_KRAKEN_READS_2.fq.gz"

    else
        log "[$ID] WARNING: candidate reads not found. Skipping."
        return 0
    fi

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        log "[$ID] WARNING: paired candidate reads incomplete. Skipping."
        return 0
    fi

    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        log "[$ID] WARNING: candidate reads are empty. Writing empty report."
        echo -e "VIRUS_TAXON\tRESCUED_READS_COUNT" > "${OUT_SUBDIR}/${ID}_bowtie2_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_high_qual_reads.fasta"
        return 0
    fi

    log "[$ID] Step 1/3: fastp low-complexity and poly-G/X filtering"

    local CLEAN_R1="${OUT_SUBDIR}/${ID}_temp_rescue_R1.fq.gz"
    local CLEAN_R2="${OUT_SUBDIR}/${ID}_temp_rescue_R2.fq.gz"
    local BAM_TEMP="${OUT_SUBDIR}/${ID}_bowtie2_rescue.bam"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$CLEAN_R1" \
        -O "$CLEAN_R2" \
        --low_complexity_filter \
        --complexity_threshold "$COMPLEXITY_THRESHOLD" \
        --trim_poly_g \
        --trim_poly_x \
        --qualified_quality_phred 20 \
        --thread "$THREADS" \
        --html "${OUT_SUBDIR}/${ID}_fastp_rescue_report.html" \
        --json "${OUT_SUBDIR}/${ID}_fastp_rescue_report.json" \
        > "${OUT_SUBDIR}/${ID}_fastp_rescue.log" 2>&1

    if [[ ! -s "$CLEAN_R1" || ! -s "$CLEAN_R2" ]]; then
        log "[$ID] No reads survived rescue fastp filtering."
        echo -e "VIRUS_TAXON\tRESCUED_READS_COUNT" > "${OUT_SUBDIR}/${ID}_bowtie2_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_high_qual_reads.fasta"
        rm -f "$CLEAN_R1" "$CLEAN_R2"
        return 0
    fi

    log "[$ID] Step 2/3: sensitive local Bowtie2 viral alignment"

    bowtie2 \
        -x "$BT2_PREFIX" \
        -1 "$CLEAN_R1" \
        -2 "$CLEAN_R2" \
        --local \
        --very-sensitive-local \
        --threads "$THREADS" \
        --no-unal \
        2> "${OUT_SUBDIR}/${ID}_bowtie2_rescue.log" | \
    samtools view \
        -@ "$THREADS" \
        -h \
        -F 4 \
        -q "$MAPQ" | \
    samtools sort \
        -@ "$THREADS" \
        -o "$BAM_TEMP" 

    if [[ -s "$BAM_TEMP" ]] && [[ "$(samtools view -c "$BAM_TEMP")" -gt 0 ]]; then

        log "[$ID] Step 3/3: exporting rescue outputs"

        samtools index "$BAM_TEMP"

        cp "$BAM_TEMP" "${OUT_SUBDIR}/${ID}_bowtie2_viral_aligned.bam"
        cp "${BAM_TEMP}.bai" "${OUT_SUBDIR}/${ID}_bowtie2_viral_aligned.bam.bai"

        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tRESCUED_READS_COUNT" > "${OUT_SUBDIR}/${ID}_bowtie2_viral_counts.tsv"

        samtools idxstats "$BAM_TEMP" | \
        awk '$3 > 0 {print $1"\t"$2"\t"$3}' >> "${OUT_SUBDIR}/${ID}_bowtie2_viral_counts.tsv"

        samtools view "$BAM_TEMP" | \
        awk 'BEGIN{OFS=""} {print ">"$1"|virus="$3"|pos="$4"|mapq="$5"|cigar="$6; print $10}' \
        > "${OUT_SUBDIR}/${ID}_rescued_high_qual_reads.fasta"

        log "[$ID] SUCCESS: sensitive rescue completed."

    else
        log "[$ID] 0 reads met Bowtie2 rescue criteria."
        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tRESCUED_READS_COUNT" > "${OUT_SUBDIR}/${ID}_bowtie2_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_high_qual_reads.fasta"
    fi

    rm -f "$BAM_TEMP" "${BAM_TEMP}.bai" "$CLEAN_R1" "$CLEAN_R2"
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

for ID in "${SAMPLES[@]}"; do
    log "=============================================================================="
    log ">>> RUNNING MODULE 3 ON SAMPLE: $ID <<<"
    log "=============================================================================="
    process_rescue "$ID"
done

log "=============================================================================="
log "MODULE 3 EXECUTION COMPLETED SUCCESSFULLY"
log "=============================================================================="
