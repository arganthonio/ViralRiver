ViralRiver is a modular RNA-seq virome analysis workflow designed for the detection, assembly, rescue, and quantification of low-abundance viral sequences from host-dominated transcriptomic datasets.

The pipeline was developed to characterize the plasma and PBMC-associated virome in individuals at risk of HIV-1 infection, but it can be applied to any paired-end RNA-seq dataset.

ViralRiver combines:

taxonomic screening,
de novo assembly,
stringent host depletion,
sensitive viral rescue mapping,
and quantitative viral profiling

into a fully automated and reproducible workflow.

Features
Detection of low-abundance viral transcripts from RNA-seq
Compatible with plasma, PBMC, tissue, and bulk RNA-seq
Automated host depletion against hg38
Kraken2-based viral candidate extraction
MEGAHIT de novo assembly
Sensitive viral rescue using Minimap2 and Bowtie2
Quantification of viral abundance
Recovery of assembled viral contigs
Modular architecture
Nextflow-ready workflow
Includes curated human virome reference databases
Workflow Overview
Raw paired-end RNA-seq
        ↓
Module 1
QC → Kraken2 → Viral candidate extraction → Assembly → Host depletion
        ↓
Module 2
Strict viral rescue and quantification (Minimap2)
        ↓
Module 3
High-sensitivity viral rescue (Bowtie2)
        ↓
Final virome abundance tables + rescued reads + viral contigs
Installation
Requirements

The following software must be available in the environment:

fastp
Kraken2
MEGAHIT
BWA
Bowtie2
Minimap2
Samtools
SeqKit
Python ≥ 3.8
Biopython
Nextflow
Conda Environment

Create the environment:

conda env create -f environment.yml

Activate:

conda activate viralriver
Reference Databases

ViralRiver includes:

hg38 reference genome
Kraken2 human-virus database
curated human core virome FASTA
prebuilt BWA, Minimap2, and Bowtie2 indexes

Default structure:

refs/
├── hg38/
├── kraken_db/
└── core_virome/
Input Format

Input FASTQ files must follow:

SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz
Running ViralRiver
Step 1 — Prepare samplesheet

Example samples.csv:

sample,fastq_1,fastq_2
SRR000001,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
SRR000002,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
Step 2 — Run pipeline
nextflow run main.nf --samples samples.csv
Output Structure
results/
├── module1/
│   ├── candidate reads
│   ├── assembled contigs
│   └── host-depleted viral contigs
│
├── module2/
│   ├── minimap2 viral counts
│   ├── rescued viral reads
│   └── BAM alignments
│
└── module3/
    ├── Bowtie2 rescue counts
    ├── high-quality rescued reads
    └── BAM alignments
Main Output Files
File	Description
*_viral_counts.tsv	Viral abundance table
*_rescued_viral_reads.fasta	Viral rescued reads
*_host_depleted_viral_contigs.fasta	Candidate assembled viral contigs
*_viral_aligned.bam	Viral alignments
*_kraken_report.txt	Kraken2 classification report
Notes
ViralRiver is optimized for low-biomass virome detection.
Module 3 performs highly sensitive rescue mapping and may recover additional low-abundance viral reads.
The workflow prioritizes specificity by combining host depletion, assembly validation, and rescue-based quantification.
Empty output files are intentionally generated when no viral signal is detected to maintain workflow reproducibility.
Citation

If you use ViralRiver, please cite:

Caruz et al. ViralRiver: a modular workflow for RNA-seq virome characterization and viral rescue in host-dominated t
