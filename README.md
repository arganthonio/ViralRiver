# ViralRiver

Software DOI:
https://doi.org/10.5281/zenodo.20254960

hg38 reference bundle DOI:
https://doi.org/10.5281/zenodo.20255121

**ViralRiver** is a modular RNA-seq virome analysis workflow designed for the detection, assembly, rescue, and quantification of low-abundance viral sequences from host-dominated transcriptomic datasets.
The pipeline was developed to characterize the plasma and PBMC-associated virome in individuals at risk of HIV-1 infection, but it can be applied to any paired-end RNA-seq dataset.

**ViralRiver** combines:

* Taxonomic screening,
* *De novo* assembly,
* Stringent host depletion,
* Sensitive viral rescue mapping,
* And quantitative viral profiling

into a fully automated and reproducible workflow.

## Features

### Biological & Analytical Capabilities
* **Detection of low-abundance viral transcripts** from RNA-seq data.
* **Versatile sample compatibility:** Optimized for plasma, PBMC, tissue, and bulk RNA-seq.
* **Quantification of viral abundance** to build accurate taxonomic expression matrices.
* **Recovery of fully assembled viral contigs** for downstream evolutionary or variant analysis.

### Automated Computational Workflow
* **Automated host depletion** via strict cross-alignment against the `hg38` reference genome.
* **Kraken2-based viral candidate extraction** utilizing customizable taxonomic identifiers (TaxIDs).
* **Robust *de novo* assembly** powered by MEGAHIT to reconstruct non-host fragments.
* **Dual-engine sensitive viral rescue** integrating both Minimap2 and Bowtie2 mapping strategies.

### Architecture & Deployment
* **Modular architecture:** Clean separation of software logic, runtime environments, and data references.
* **Nextflow-ready workflow:** Out-of-the-box support for parallel processing, pipeline resuming (`-resume`), and containerization.
* **Turnkey references:** Includes a curated human virome reference database ready for deployment.

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
Module 4 extracts per-sample total RNA-seq read depth from fastp reports and generates a table required for RPM normalization.
        ↓
Final virome abundance tables + rescued reads + viral contigs
```

# Installation

## Requirements 1

ViralRiver requires:

- Conda or Miniconda
- Nextflow (>=22.10)

All other software dependencies are installed automatically through the provided Conda environment, including:

- fastp
- Kraken2
- MEGAHIT
- BWA
- Bowtie2
- Minimap2
- Samtools
- SeqKit
- Python ≥ 3.8
- Biopython

## Requirements 2

### hg38 reference bundle

A pre-indexed hg38 reference genome compatible with ViralRiver is available from Zenodo.

Download:

```bash
wget https://zenodo.org/records/20255121/files/ViralRiver_hg38_reference_bundle_v1.0.tar.gz
tar -xzf ViralRiver_hg38_reference_bundle_v1.0.tar.gz
```

## Included reference resources

ViralRiver includes:

- curated human core virome reference FASTA
- prebuilt Minimap2 and Bowtie2 indexes
- lightweight Kraken2 human-virus database

## Create environment

```bash
conda env create -f environment.yml
```

## Activate environment

```bash
conda activate viralriver
```

## Verify installation

```bash
nextflow -version
fastp --version
kraken2 --version
```

## Input Format

FASTQ files should follow:

```text
SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz
```

Example general `samples.csv`:

```csv
sample,fastq_1,fastq_2
SRR000001,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
SRR000002,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
```
Example specific `samples.csv`: 

````csv
sample,fastq_1,fastq_2
SRR15413671,/mnt/e/pipeline_test/HIV_test/vaccine/SRR15413671_1.fastq.gz,/mnt/e/pipeline_test/HIV_test/vaccine/SRR15413671_2.fastq.gz
SRR15413654,/mnt/e/pipeline_test/HIV_test/vaccine/SRR15413654_1.fastq.gz,/mnt/e/pipeline_test/HIV_test/vaccine/SRR15413654_2.fastq.gz
```



## Example Public RNA-seq Datasets

Example paired-end FASTQ files can be downloaded directly from the European Nucleotide Archive (ENA):

```bash
# Example 1
curl -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/024/SRR16948824/SRR16948824_1.fastq.gz"
curl -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/024/SRR16948824/SRR16948824_2.fastq.gz"

# Example 2
curl -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/071/SRR32014171/SRR32014171_1.fastq.gz"
curl -O "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/071/SRR32014171/SRR32014171_2.fastq.gz"
```

Example `samples.csv`:

```csv
sample,fastq_1,fastq_2
SRR16948824,SRR16948824_1.fastq.gz,SRR16948824_2.fastq.gz
SRR32014171,SRR32014171_1.fastq.gz,SRR32014171_2.fastq.gz
```

## Run ViralRiver:

```bash
nextflow run main.nf \
  --samples /path/to/samples.csv \
  --host_ref /path/to/ref_hg38/hg38_full.fa \
  --viral_fasta /path/to/viral_master/human.virus.selected.fasta \
  --kraken_db /path/to/kraken_humanvirus_db \
  -with-conda
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
