# ViralRiver

**ViralRiver** is a modular RNA-seq virome analysis workflow designed for the detection, assembly, rescue, and quantification of low-abundance viral sequences from host-dominated transcriptomic datasets.

## Workflow Overview

```text
Raw paired-end RNA-seq
        ↓
Module 1: QC + Kraken2 + Assembly
        ↓
Module 2: Minimap2 viral quantification
        ↓
Module 3: Bowtie2 sensitive rescue
        ↓
Final virome abundance tables + rescued reads + viral contigs
```

## Installation

```bash
conda env create -f environment.yml
conda activate viralriver
```

## Input Format

FASTQ files should follow:

```text
SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz
```

Example `samples.csv`:

```csv
sample,fastq_1,fastq_2
SRR000001,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
SRR000002,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
```

## Run

```bash
nextflow run main.nf --samples samples.csv
```

## Outputs

```text
results/
├── module1/
│   ├── candidate reads
│   ├── assembled contigs
│   └── host-depleted viral contigs
├── module2/
│   ├── minimap2 viral counts
│   ├── rescued viral reads
│   └── BAM alignments
└── module3/
    ├── Bowtie2 rescue counts
    ├── high-quality rescued reads
    └── BAM alignments
```

## Main output files

| File | Description |
|---|---|
| `*_viral_counts.tsv` | Viral abundance table |
| `*_rescued_viral_reads.fasta` | Viral rescued reads |
| `*_host_depleted_viral_contigs.fasta` | Candidate assembled viral contigs |
| `*_viral_aligned.bam` | Viral alignments |
| `*_kraken_report.txt` | Kraken2 classification report |

## Citation

Caruz et al. ViralRiver: a modular workflow for RNA-seq virome characterization and viral rescue in host-dominated transcriptomic datasets.

## License

MIT License
