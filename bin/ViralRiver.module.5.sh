#!/usr/bin/env bash
# ==============================================================================
# ViralRiver - Module 5: MultiQC Custom Quality Control Report
# ==============================================================================

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 -i <input_dir> -o <output_dir>"
    echo "  -i  Directory containing MultiQC input files"
    echo "  -o  Output directory"
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
    usage
fi

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

log "Verifying MultiQC dependency..."

command -v multiqc >/dev/null || {
    echo "ERROR: MultiQC not found."
    exit 1
}

###############################################################################
# Generate MultiQC matrices and ViralRiver summary
###############################################################################

python3 << EOF

import os
import glob
import csv
import shutil

INPUT_DIR = "${INPUT_DIR}"
OUTPUT_DIR = "${OUTPUT_DIR}"

###############################################################################
# Build matrices for MultiQC
###############################################################################

def process_files(files, out_filename):

    if not files:
        return

    viruses = set()
    data = {}

    for f in files:

        sample = os.path.basename(f).split("_")[0]
        data[sample] = {}

        with open(f) as fh:

            for line in fh:

                line = line.strip()

                if (
                    not line
                    or line.startswith("#")
                    or line.startswith("VIRUS_TAXON")
                ):
                    continue

                cols = line.split("\t")

                try:

                    virus = cols[0]
                    count = int(cols[-1])

                    if count > 0:
                        data[sample][virus] = count
                        viruses.add(virus)

                except:
                    pass

    if not viruses:
        return

    viruses = sorted(viruses)

    with open(out_filename, "w") as out:

        out.write("\t".join(["Sample"] + viruses) + "\n")

        for sample in sorted(data):

            row = [sample]

            for virus in viruses:
                row.append(str(data[sample].get(virus,0)))

            out.write("\t".join(row) + "\n")

###############################################################################
# Generate matrices
###############################################################################

bowtie_files = glob.glob(
    os.path.join(INPUT_DIR,"*_bowtie2_viral_counts.tsv")
)

process_files(
    bowtie_files,
    os.path.join(INPUT_DIR,"summary_bowtie2_mqc.tsv")
)

minimap_files = [

    f for f in glob.glob(
        os.path.join(INPUT_DIR,"*_viral_counts.tsv")
    )

    if "_bowtie2" not in f
    and "summary_" not in os.path.basename(f)

]

process_files(
    minimap_files,
    os.path.join(INPUT_DIR,"summary_minimap2_mqc.tsv")
)

###############################################################################
# ViralRiver summary table
###############################################################################

summary_file = os.path.join(
    OUTPUT_DIR,
    "viralriver_summary.tsv"
)

with open(summary_file,"w",newline="") as fout:

    writer = csv.writer(fout,delimiter="\t")

    writer.writerow([
        "Sample",
        "Minimap2_viruses",
        "Minimap2_reads",
        "Bowtie2_viruses",
        "Bowtie2_reads",
        "Total_viral_reads",
        "Status"
    ])

    for mm_file in sorted(minimap_files):

        sample = os.path.basename(mm_file).replace(
            "_viral_counts.tsv",
            ""
        )

        #######################################################################
        # Minimap2
        #######################################################################

        mm_reads = 0
        mm_viruses = 0

        with open(mm_file) as fh:

            next(fh)

            for line in fh:

                if not line.strip():
                    continue

                cols = line.rstrip().split("\t")

                mm_viruses += 1

                try:
                    mm_reads += int(cols[2])
                except:
                    pass

        #######################################################################
        # Bowtie2
        #######################################################################

        bt_reads = 0
        bt_viruses = 0

        bt_file = os.path.join(
            INPUT_DIR,
            sample + "_bowtie2_viral_counts.tsv"
        )

        if os.path.exists(bt_file):

            with open(bt_file) as fh:

                next(fh)

                for line in fh:

                    if not line.strip():
                        continue

                    cols = line.rstrip().split("\t")

                    bt_viruses += 1

                    try:
                        bt_reads += int(cols[2])
                    except:
                        pass

        total = mm_reads + bt_reads

        if total > 0:
            status = "Virus detected"
        else:
            status = "No viruses detected"

        writer.writerow([
            sample,
            mm_viruses,
            mm_reads,
            bt_viruses,
            bt_reads,
            total,
            status
        ])

###############################################################################
# Copy summary so MultiQC can detect it
###############################################################################

shutil.copy(
    summary_file,
    os.path.join(INPUT_DIR,"viralriver_summary.tsv")
)

print("Summary written:",summary_file)




EOF

###############################################################################
# Run MultiQC
###############################################################################

log "Running Custom MultiQC on: ${INPUT_DIR}"

multiqc "${INPUT_DIR}" \
    -o "${OUTPUT_DIR}" \
    --filename "viralriver_multiqc_report.html" \
    --config "multiqc_config.yaml" \
    --force

log "MultiQC report generated successfully."