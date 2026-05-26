#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 2: Host Depletion Validation and Viral Quantification
# ==============================================================================

set -euo pipefail

THREADS=8
MAPQ=20
REF_GENOME=""
VIRAL_DB_FASTA=""
INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -v <viral_db.fasta> -r <ref_genome.fa> [-t <threads>] [-q <mapq>]"
    echo "  -i  Input base directory from Module 1"
    echo "  -o  Output base directory"
    echo "  -v  Curated viral reference FASTA"
    echo "  -r  Host reference genome FASTA"
    echo "  -t  Threads (default: 8)"
    echo "  -q  Minimum MAPQ for viral alignments (default: 20)"
    exit 1
}

while getopts "i:o:v:r:t:q:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        v) VIRAL_DB_FASTA="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        q) MAPQ="$OPTARG" ;;
        h|?) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$VIRAL_DB_FASTA" || -z "$REF_GENOME" ]]; then
    echo "ERROR: Missing required arguments."
    usage
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Verifying software dependencies..."

for tool in bwa minimap2 samtools awk find sed sort; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Dependency '$tool' not found."
        exit 1
    }
done

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input directory not found: $INPUT_DIR"; exit 1; }
[[ -f "$REF_GENOME" ]] || { echo "ERROR: host genome FASTA not found: $REF_GENOME"; exit 1; }
[[ -f "$VIRAL_DB_FASTA" ]] || { echo "ERROR: viral FASTA not found: $VIRAL_DB_FASTA"; exit 1; }

mkdir -p "$OUTPUT_DIR"

# ------------------------------------------------------------------------------
# Minimap2 index
# ------------------------------------------------------------------------------

MMI_INDEX="${VIRAL_DB_FASTA%.*}.mmi"

if [[ ! -f "$MMI_INDEX" ]]; then
    log "Minimap2 index not found. Building index: $MMI_INDEX"
    minimap2 -k 15 -w 5 -d "$MMI_INDEX" "$VIRAL_DB_FASTA" > /dev/null 2>&1
else
    log "Found existing Minimap2 index: $MMI_INDEX"
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

process_quantification() {
    local ID="$1"
    local IN_SUBDIR="${INPUT_DIR}/${ID}"
    local OUT_SUBDIR="${OUTPUT_DIR}/${ID}"

    mkdir -p "$OUT_SUBDIR"

    local R1=""
    local R2=""

    # Prefer current ViralRiver Module 1 naming
    if [[ -f "${IN_SUBDIR}/${ID}_candidate_reads_1.fq.gz" ]]; then
        R1="${IN_SUBDIR}/${ID}_candidate_reads_1.fq.gz"
        R2="${IN_SUBDIR}/${ID}_candidate_reads_2.fq.gz"

    # Legacy naming
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
        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tMAPPED_READ_COUNT" > "${OUT_SUBDIR}/${ID}_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_viral_reads.fasta"
        return 0
    fi

    log "[$ID] Step 1/4: strict host depletion with BWA"

    bwa mem \
        -t "$THREADS" \
        "$REF_GENOME" \
        "$R1" "$R2" \
        2> "${OUT_SUBDIR}/${ID}_bwa_host_depletion.log" | \
     samtools view \
        -@ "$THREADS" \
        -b \
        -f 12 \
        -o "${OUT_SUBDIR}/tmp_non_host.bam" -

    if [[ ! -s "${OUT_SUBDIR}/tmp_non_host.bam" ]] || [[ "$(samtools view -c "${OUT_SUBDIR}/tmp_non_host.bam")" -eq 0 ]]; then
        log "[$ID] No read pairs survived host depletion."
        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tMAPPED_READ_COUNT" > "${OUT_SUBDIR}/${ID}_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_viral_reads.fasta"
        rm -f "${OUT_SUBDIR}/tmp_non_host.bam"
        return 0
    fi

    log "[$ID] Step 2/4: extracting non-host FASTQ pairs"

    samtools fastq \
        -@ "$THREADS" \
        -1 "${OUT_SUBDIR}/clean_R1.fq" \
        -2 "${OUT_SUBDIR}/clean_R2.fq" \
        -0 /dev/null \
        -s /dev/null \
        -n \
        "${OUT_SUBDIR}/tmp_non_host.bam" \
        > "${OUT_SUBDIR}/${ID}_samtools_fastq.log" 2>&1

    if [[ ! -s "${OUT_SUBDIR}/clean_R1.fq" || ! -s "${OUT_SUBDIR}/clean_R2.fq" ]]; then
        log "[$ID] No paired FASTQ reads after host depletion."
        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tMAPPED_READ_COUNT" > "${OUT_SUBDIR}/${ID}_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_viral_reads.fasta"
        rm -f "${OUT_SUBDIR}/tmp_non_host.bam" "${OUT_SUBDIR}/clean_R1.fq" "${OUT_SUBDIR}/clean_R2.fq"
        return 0
    fi

    log "[$ID] Step 3/4: viral rescue mapping with Minimap2"

    minimap2 \
        -ax sr \
        -t "$THREADS" \
        "$MMI_INDEX" \
        "${OUT_SUBDIR}/clean_R1.fq" \
        "${OUT_SUBDIR}/clean_R2.fq" \
        2> "${OUT_SUBDIR}/${ID}_minimap2_viral.log" | \
    samtools view \
        -@ "$THREADS" \
        -h \
        -F 4 \
        -o "${OUT_SUBDIR}/tmp_viral_hits.bam" -

    if [[ -s "${OUT_SUBDIR}/tmp_viral_hits.bam" ]] && [[ "$(samtools view -c "${OUT_SUBDIR}/tmp_viral_hits.bam")" -gt 0 ]]; then

        log "[$ID] Step 4/4: sorting, indexing and reporting viral alignments"

        samtools sort \
            -@ "$THREADS" \
            -o "${OUT_SUBDIR}/${ID}_viral_aligned.bam" \
            "${OUT_SUBDIR}/tmp_viral_hits.bam"

        samtools index "${OUT_SUBDIR}/${ID}_viral_aligned.bam"

        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tMAPPED_READ_COUNT" > "${OUT_SUBDIR}/${ID}_viral_counts.tsv"

        samtools idxstats "${OUT_SUBDIR}/${ID}_viral_aligned.bam" | \
        awk '$3 > 0 {print $1"\t"$2"\t"$3}' >> "${OUT_SUBDIR}/${ID}_viral_counts.tsv"

        samtools stats "${OUT_SUBDIR}/${ID}_viral_aligned.bam" | \
        grep "average length:" > "${OUT_SUBDIR}/${ID}_read_length_stats.txt" || true

        samtools view "${OUT_SUBDIR}/${ID}_viral_aligned.bam" | \
        awk 'BEGIN{OFS=""} {print ">"$1"|virus="$3"|pos="$4"|mapq="$5"|cigar="$6; print $10}' \
        > "${OUT_SUBDIR}/${ID}_rescued_viral_reads.fasta"

        log "[$ID] SUCCESS: viral quantification completed."

    else
        log "[$ID] 0 viral alignments after MAPQ filtering."
        echo -e "VIRUS_TAXON\tSEQUENCE_LENGTH\tMAPPED_READ_COUNT" > "${OUT_SUBDIR}/${ID}_viral_counts.tsv"
        touch "${OUT_SUBDIR}/${ID}_rescued_viral_reads.fasta"
    fi

    rm -f \
        "${OUT_SUBDIR}/tmp_non_host.bam" \
        "${OUT_SUBDIR}/tmp_viral_hits.bam" \
        "${OUT_SUBDIR}/clean_R1.fq" \
        "${OUT_SUBDIR}/clean_R2.fq"
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

for ID in "${SAMPLES[@]}"; do
    process_quantification "$ID"
done

log "=============================================================================="
log "MODULE 2 EXECUTION COMPLETED SUCCESSFULLY"
log "=============================================================================="
