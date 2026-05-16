#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 1: QC, Candidate Read Extraction, Assembly and Host Depletion
# ==============================================================================

set -euo pipefail

THREADS=8
KRAKEN_DB=""
REF_GENOME=""
INPUT_DIR=""
OUTPUT_DIR=""
KRAKEN_TAXID=1

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir> -d <kraken_db> -r <ref_genome.fa> [-t <threads>] [-x <kraken_taxid>]"
    echo "  -i  Directory containing paired-end fastq.gz files (*_1.fastq.gz/*_2.fastq.gz)"
    echo "  -o  Base directory for output results"
    echo "  -d  Path to Kraken2 database"
    echo "  -r  Path to host reference genome FASTA"
    echo "  -t  Number of threads (default: 8)"
    echo "  -x  Kraken taxid to extract candidate reads (default: 1; use 10239 for strict NCBI Viruses)"
    exit 1
}

while getopts "i:o:d:r:t:x:h" opt; do
    case ${opt} in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        d) KRAKEN_DB="$OPTARG" ;;
        r) REF_GENOME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        x) KRAKEN_TAXID="$OPTARG" ;;
        h|?) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$KRAKEN_DB" || -z "$REF_GENOME" ]]; then
    echo "ERROR: Missing required arguments."
    usage
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Verifying software dependencies..."

for tool in fastp kraken2 extract_kraken_reads.py megahit bwa samtools seqkit python awk gzip find sed sort df cut; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Dependency '$tool' not found."
        exit 1
    }
done

python -c "import Bio" &>/dev/null || {
    echo "ERROR: Biopython is not installed in the active Python environment."
    exit 1
}

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input directory not found: $INPUT_DIR"; exit 1; }
[[ -d "$KRAKEN_DB" ]] || { echo "ERROR: Kraken2 database not found: $KRAKEN_DB"; exit 1; }
[[ -f "$REF_GENOME" ]] || { echo "ERROR: host reference FASTA not found: $REF_GENOME"; exit 1; }

TMP_BASE="${OUTPUT_DIR}/tmp_work"
mkdir -p "$OUTPUT_DIR" "$TMP_BASE"

process_sample() {
    local ID="$1"
    local OUT_DIR="${OUTPUT_DIR}/${ID}"
    local WORKDIR="${TMP_BASE}/${ID}"
    local RAM_TMP="/dev/shm/viralriver_tmp_${ID}"

    local R1="${INPUT_DIR}/${ID}_1.fastq.gz"
    local R2="${INPUT_DIR}/${ID}_2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        log "[$ID] WARNING: Missing paired FASTQ files. Skipping."
        return 0
    fi

    mkdir -p "$OUT_DIR" "$WORKDIR" "$OUT_DIR/fastp_reports"

    if [[ -d /dev/shm ]]; then
        mkdir -p "$RAM_TMP" || RAM_TMP="${WORKDIR}/viral_tmp"
    else
        RAM_TMP="${WORKDIR}/viral_tmp"
    fi
    mkdir -p "$RAM_TMP"

    trap '[[ -d "${RAM_TMP:-}" ]] && rm -rf "$RAM_TMP"; log "[$ID] interrupted; temporary files removed."; exit 1' SIGINT SIGTERM

    log "[$ID] Processing sample..."
    cd "$WORKDIR"

    log "[$ID] Step 1/6: Quality filtering with fastp"
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o clean_1.fq.gz \
        -O clean_2.fq.gz \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --html "${OUT_DIR}/fastp_reports/${ID}_fastp.html" \
        --json "${OUT_DIR}/fastp_reports/${ID}_fastp.json" \
        > "${OUT_DIR}/${ID}_fastp.log" 2>&1

    log "[$ID] Step 2/6: Kraken2 taxonomic screening"
    kraken2 \
        --db "$KRAKEN_DB" \
        --threads "$THREADS" \
        --paired \
        --confidence 0 \
        --use-names \
        --gzip-compressed \
        --report "${OUT_DIR}/${ID}_kraken_report.txt" \
        --output kraken.out \
        clean_1.fq.gz clean_2.fq.gz \
        > "${OUT_DIR}/${ID}_kraken.log" 2>&1

    if [[ -d /dev/shm ]]; then
        local RAM_FREE
        RAM_FREE=$(df --output=avail /dev/shm | awk 'NR==2 {print $1}')

        if [[ "$RAM_FREE" -lt 1000000 ]]; then
            log "[$ID] WARNING: insufficient /dev/shm space. Using disk temporary directory."
            rm -rf "$RAM_TMP"
            RAM_TMP="${WORKDIR}/viral_tmp"
            mkdir -p "$RAM_TMP"
        fi
    fi

    log "[$ID] Step 3/6: Extracting candidate reads with Kraken taxid ${KRAKEN_TAXID}"
    extract_kraken_reads.py \
        -k kraken.out \
        -r "${OUT_DIR}/${ID}_kraken_report.txt" \
        -s1 clean_1.fq.gz \
        -s2 clean_2.fq.gz \
        -t "$KRAKEN_TAXID" \
        --include-children \
        --fastq-output \
        --output "$RAM_TMP/viral_1.fq" \
        --output2 "$RAM_TMP/viral_2.fq" \
        > "${OUT_DIR}/${ID}_extract_candidate_reads.log" 2>&1

    rm -f kraken.out

    if [[ ! -s "$RAM_TMP/viral_1.fq" ]]; then
        log "[$ID] No candidate reads detected. Skipping assembly."
        touch "${OUT_DIR}/${ID}_candidate_reads_1.fq.gz"
        touch "${OUT_DIR}/${ID}_candidate_reads_2.fq.gz"
        touch "${OUT_DIR}/${ID}_total_assembled_contigs.fasta"
        touch "${OUT_DIR}/${ID}_host_depleted_viral_contigs.fasta"
        rm -rf "$WORKDIR" "$RAM_TMP"
        trap - SIGINT SIGTERM
        return 0
    fi

    gzip -c "$RAM_TMP/viral_1.fq" > "${OUT_DIR}/${ID}_candidate_reads_1.fq.gz"
    gzip -c "$RAM_TMP/viral_2.fq" > "${OUT_DIR}/${ID}_candidate_reads_2.fq.gz"

    log "[$ID] Step 4/6: De novo assembly with MEGAHIT"
    megahit \
        -1 "$RAM_TMP/viral_1.fq" \
        -2 "$RAM_TMP/viral_2.fq" \
        -o megahit_out \
        -t "$THREADS" \
        --k-list 21,33,55,77 \
        --min-contig-len 150 \
        > "${OUT_DIR}/${ID}_megahit.log" 2>&1

    log "[$ID] Step 5/6: Host depletion of assembled contigs"

    local n_final=0

    if [[ -s "megahit_out/final.contigs.fa" ]]; then
        seqkit sort \
            --by-length \
            --reverse \
            megahit_out/final.contigs.fa \
            -o "${OUT_DIR}/${ID}_total_assembled_contigs.fasta"

        bwa mem \
            -t "$THREADS" \
            "$REF_GENOME" \
            megahit_out/final.contigs.fa \
            2> "${OUT_DIR}/${ID}_bwa_contigs_vs_host.log" | \
        samtools view \
            -@ "$THREADS" \
            -b \
            -o contigs_vs_host.bam -

        samtools view \
            -@ "$THREADS" \
            -f 4 \
            contigs_vs_host.bam | \
        cut -f1 | \
        sort -u -T "$WORKDIR" > pure_viral_ids.txt

        if [[ -s pure_viral_ids.txt ]]; then
            seqkit grep \
                -f pure_viral_ids.txt \
                megahit_out/final.contigs.fa \
                -o "${OUT_DIR}/${ID}_host_depleted_viral_contigs.fasta"

            n_final=$(grep -c "^>" "${OUT_DIR}/${ID}_host_depleted_viral_contigs.fasta" 2>/dev/null || echo 0)
            log "[$ID] SUCCESS: ${n_final} candidate non-host contigs recovered."
        else
            log "[$ID] No contigs survived host depletion."
            touch "${OUT_DIR}/${ID}_host_depleted_viral_contigs.fasta"
        fi
    else
        log "[$ID] Assembly generated 0 contigs."
        touch "${OUT_DIR}/${ID}_total_assembled_contigs.fasta"
        touch "${OUT_DIR}/${ID}_host_depleted_viral_contigs.fasta"
    fi

    log "[$ID] Step 6/6: Cleanup"
    rm -rf "$WORKDIR" "$RAM_TMP"
    trap - SIGINT SIGTERM
}

log "Discovering samples in ${INPUT_DIR}..."

mapfile -t SAMPLES < <(
    find "$INPUT_DIR" -maxdepth 1 -name "*_1.fastq.gz" \
    | sed 's|.*/||' \
    | sed 's/_1.fastq.gz//' \
    | sort
)

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
    echo "ERROR: No valid paired-end samples (*_1.fastq.gz) found in ${INPUT_DIR}"
    exit 1
fi

log "Found ${#SAMPLES[@]} samples to process."

for ID in "${SAMPLES[@]}"; do
    process_sample "$ID"
done

rm -rf "$TMP_BASE"

log "=============================================================================="
log "MODULE 1 EXECUTION COMPLETED SUCCESSFULLY"
log "=============================================================================="