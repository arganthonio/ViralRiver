# 🧬 ViralRiver

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.20254960-blue)](https://doi.org/10.5281/zenodo.20254960)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A526.0-brightgreen)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)
[![Docker](https://img.shields.io/badge/Docker-supported-blue.svg)](https://www.docker.com/)

---

# ViralRiver

**ViralRiver** is an end-to-end RNA-seq virome analysis workflow that integrates taxonomic classification, host depletion, de novo assembly, dual viral alignment, viral rescue, abundance estimation and interactive visualization into a single reproducible pipeline.

The workflow was originally developed to characterize the plasma and PBMC-associated virome in individuals at risk of HIV-1 infection, although it can be applied to any paired-end RNA-seq dataset generated from host-dominated transcriptomic samples.

Unlike existing virome analysis pipelines, ViralRiver emphasizes **accessibility**, **reproducibility**, and **ease of use**. The software is distributed both as a standard Nextflow workflow for bioinformaticians and as a Docker-based graphical web application that allows researchers with limited computational experience to perform complete virome analyses directly from a web browser.

---

# Highlights

✔ End-to-end RNA-seq virome workflow

✔ Quality control, taxonomic classification and host depletion

✔ De novo assembly of candidate viral contigs

✔ Dual viral alignment using **Minimap2** and **Bowtie2**

✔ Sensitive viral read rescue

✔ Viral abundance estimation

✔ Interactive HTML report generated with **MultiQC**

✔ Multi-sample processing through **Nextflow**

✔ Browser-based graphical interface built with **Streamlit**

✔ Docker image for one-click deployment

✔ Fully reproducible workflow

---

# Workflow Overview

```text
                Raw paired-end RNA-seq
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│ Module 1                                                   │
│ Quality control, taxonomic classification,                 │
│ host depletion and de novo assembly                        │
└────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│ Module 2                                                   │
│ Viral alignment and abundance estimation                   │
│ using Minimap2                                             │
└────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│ Module 3                                                   │
│ Sensitive viral rescue using Bowtie2                       │
└────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│ Module 4                                                   │
│ Extraction of sequencing depth from fastp reports          │
│ and RPM normalization                                      │
└────────────────────────────────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│ Module 5                                                   │
│ Interactive MultiQC report generation                      │
└────────────────────────────────────────────────────────────┘
                         │
                         ▼
          Interactive HTML report + Output tables
```

> **Figure 1.** Overview of the ViralRiver workflow.  
> (Replace this diagram by the workflow figure used in the manuscript.)

---


The report is automatically generated as

results/
└── viralriver_multiqc_report.html


and can be opened locally using any modern web browser.

---

# Software Requirements

ViralRiver requires:

- Nextflow ≥ 26
- Conda or Miniconda
- Python 3.10

All third-party software dependencies are automatically installed from the provided Conda environment.

Main software components include:

- fastp
- Kraken2
- KrakenTools
- Minimap2
- Bowtie2
- Samtools
- MEGAHIT
- SeqKit
- BLAST+
- SRA Toolkit
- Biopython
- Pigz
- Wget


---

# Reference resources

Two Zenodo packages are available.

## ViralRiver software

Software DOI

https://doi.org/10.5281/zenodo.20254960

---

## hg38 reference bundle

A pre-indexed hg38 reference genome compatible with ViralRiver is available from Zenodo.

https://doi.org/10.5281/zenodo.20255121

Download:

```bash
wget https://zenodo.org/records/20255121/files/ViralRiver_hg38_reference_bundle_v1.0.tar.gz

tar -xzf ViralRiver_hg38_reference_bundle_v1.0.tar.gz
```

---

## Included reference resources

The ViralRiver package also includes:

- Curated human virome reference FASTA
- Pre-built Minimap2 indexes
- Pre-built Bowtie2 indexes
- Lightweight Kraken2 viral database

# Create the Conda Environment

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate ViralRiver
```

Verify the installation:

```bash
nextflow -version
fastp --version
kraken2 --version
minimap2 --version
bowtie2 --version
```

---

# Set execution permissions (Linux/macOS)

Before running ViralRiver for the first time, ensure that all shell scripts are executable.

```bash
chmod +x ViralRiver.module.*.sh
chmod +x main.nf
```

---

# Input Format

ViralRiver accepts paired-end RNA-seq FASTQ files.

Expected naming convention:

```text
sample_1.fastq.gz
sample_2.fastq.gz
```

Samples are specified through a CSV file.

Example:

```csv
sample,fastq_1,fastq_2
SRR000001,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
SRR000002,/path/sample_1.fastq.gz,/path/sample_2.fastq.gz
```

A single execution may contain one or many samples.

Nextflow automatically schedules the analyses in parallel whenever computational resources are available.

---

# Example Public RNA-seq datasets

Example paired-end RNA-seq datasets can be downloaded directly from the European Nucleotide Archive (ENA).

Example 1

```bash
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/024/SRR16948824/SRR16948824_1.fastq.gz

curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR169/024/SRR16948824/SRR16948824_2.fastq.gz
```

Example 2

```bash
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/071/SRR32014171/SRR32014171_1.fastq.gz

curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR320/071/SRR32014171/SRR32014171_2.fastq.gz
```

Example samples.csv

```csv
sample,fastq_1,fastq_2
SRR16948824,SRR16948824_1.fastq.gz,SRR16948824_2.fastq.gz
SRR32014171,SRR32014171_1.fastq.gz,SRR32014171_2.fastq.gz
```

It is recommended to place these files in the samples directory.

---

# Running ViralRiver

## Command-line version

Run ViralRiver using Nextflow:

```bash
nextflow run main.nf \
  --samples exmaples/samples.csv \
  --host_ref ref/hg38/hg38_full.fa \
  --viral_fasta ref/core_virome/human.virus.selected.fasta \
  --kraken_db ref/kraken_db
```


# Running ViralRiver using Docker (Graphical Interface)

The Streamlit graphical interface allows users with limited computational experience to execute complete virome analyses.

The interface provides:

- Input file validation.
- Parameter configuration.
- Live monitoring of the pipeline.
- Execution progress.
- Automatic report generation.
- Direct access to the MultiQC report.


![ViralRiver graphical interface](docs/images/gui.png)
**Figure 2. ViralRiver graphical interface.**

---

## Prerequisites

Make sure **Docker** is installed and running on your system before proceeding:
* **Windows / macOS:** [Docker Desktop](https://www.docker.com/products/docker-desktop/)
* **Linux:** [Docker Engine](https://docs.docker.com/engine/install/)

---

## Step 1: Set Up Your Working Directory

Create a workspace folder on your computer (e.g., `my_analysis/`) and structure your files and subdirectories as shown below:

```text
my_analysis/
├── viralriver.bat   <-- (or viralriver.sh for Linux/macOS)
├── samples/         <-- Place your input paired-end FASTQ files (.fastq.gz) here
└── ref/             <-- Downloaded from the ViralRiver repository
    └── hg38/        <-- ⚠️ Place the human host reference (hg38_full.fa) here
```

> 💡 **Note on the `ref/` folder:** You can copy the `ref/` directory directly from the root of the [ViralRiver Repository](https://github.com/arganthonio/ViralRiver). However, you must **manually download and place** the human reference genome file (`hg38_full.fa`) inside the `ref/hg38/` subdirectory.

---

## Step 2: Download the Launcher Script

Go to the [`docker/`](https://github.com/arganthonio/ViralRiver/tree/main/docker) directory in the GitHub repository and download the appropriate launcher script into your working directory (`my_analysis/`):

* **Windows:** Download `viralriver.bat`
* **Linux / macOS:** Download `viralriver.sh`

---

## Step 3: Launch the Application

### On Windows
1. Ensure **Docker Desktop** is active and running.
2. **Double-click** `viralriver.bat` (or execute it via `cmd`).
3. The script will automatically pull the required Docker image and open the web interface in your default browser at `http://localhost:8501`.

---

### On Linux / macOS
1. Open a terminal inside your working directory (`my_analysis/`).
2. Make the launcher script executable (only needed once):
   ```bash
   chmod +x viralriver.sh
   ```
3. Run the script:
   ```bash
   ./viralriver.sh
   ```
4. The application will pull the image and open the browser automatically at `http://localhost:8501`.

# Interactive MultiQC Report

After the pipeline finishes, ViralRiver automatically generates

```text
results/
    viralriver_multiqc_report.html
```

The report summarizes the results produced by all modules except Kraken2.

Unlike static reports, MultiQC allows users to interactively explore the data.

Available features include:

- Interactive bar plots.
- Viral abundance comparison.
- Minimap2 versus Bowtie2 comparison.
- Publication-quality figure export (PNG, SVG, PDF).
- Interactive tables.
- Download of processed data.
- Fast navigation through all samples.

No additional software is required.

Simply open

```text
viralriver_multiqc_report.html
```

using any modern web browser.

Insert here:

**Figure 3. MultiQC interactive report.**

---

# Output Directory Structure

The pipeline generates the following directory structure.

```text
results/

├── module1/
│   ├── quality_control/
│   ├── kraken2/
│   ├── assembly/
│   └── host_depletion/
│
├── module2/
│   ├── minimap2/
│   ├── viral_counts/
│   └── bam/
│
├── module3/
│   ├── bowtie2/
│   ├── rescued_reads/
│   └── bam/
│
├── module4/
│   ├── rpm_table.tsv
│   └── sequencing_depth.tsv
│
├── viralriver_multiqc_report.html
│
└── viralriver_multiqc_report_data/
```

---

# Main Output Files

| File | Description |
|------|-------------|
| `*_viral_counts.tsv` | Viral abundance estimates obtained with Minimap2 |
| `*_bowtie2_viral_counts.tsv` | Viral abundance after Bowtie2 rescue |
| `*_rescued_viral_reads.fasta` | Rescued viral reads |
| `*_host_depleted_viral_contigs.fasta` | Assembled candidate viral contigs |
| `*.bam` | Viral alignments |
| `rpm_table.tsv` | Reads per million normalization |
| `viralriver_multiqc_report.html` | Interactive MultiQC report |

---

# Why ViralRiver?

Compared with existing RNA-seq virome pipelines, ViralRiver provides:

- End-to-end automated analysis.
- Dual viral alignment using Minimap2 and Bowtie2.
- Sensitive viral read rescue.
- Multi-sample processing.
- Interactive visualization.
- Browser-based graphical interface.
- Docker deployment.
- Reproducible Nextflow workflow.
- Automatic organization of all results.

---

# Frequently Asked Questions

### Can ViralRiver analyse multiple samples?

Yes.

Simply include all paired-end datasets in the `samples.csv` file.

Nextflow automatically schedules the analyses.

---

### Can I use another viral reference database?

Yes.

Any FASTA database compatible with Minimap2 and Bowtie2 can be used.

---

### Does ViralRiver work on Windows?

Yes.

The recommended option is the Docker graphical interface.

---

### Does ViralRiver work on Linux?

Yes.

Both the command-line and Docker versions are supported.

---

### Does ViralRiver generate publication-quality figures?

Yes.

The MultiQC report allows exporting figures in PNG, SVG and PDF formats.

---

# Citation

If you use ViralRiver in your research, please cite:

```text
Caruz A. et al.

ViralRiver: an accessible end-to-end workflow for RNA-seq virome profiling using dual viral alignment and interactive visualization.

(Bioinformatics, under review)
```

Software DOI

https://doi.org/10.5281/zenodo.20254960

Reference bundle DOI

https://doi.org/10.5281/zenodo.20255121

---

# License

ViralRiver is distributed under the MIT License.

---

# Acknowledgements

The authors thank all contributors involved in the development and testing of ViralRiver.

The project has been developed to facilitate reproducible and accessible virome analyses from RNA-seq data for both bioinformaticians and experimental researchers.

---

# Roadmap

Future releases are expected to include:

- Additional viral reference databases.
- Long-read RNA sequencing support.
- Automatic taxonomic annotation.
- Extended MultiQC visualizations.
- Cloud deployment.
- Galaxy integration.
- Workflow benchmarking against additional virome pipelines.
