
# RNA-seq Analysis Pipeline Using Snakemake

An end-to-end RNA-seq analysis workflow built with **Snakemake** to process paired-end sequencing reads, quantify gene expression, and identify differentially expressed genes.

The pipeline integrates quality control, read trimming, alignment, gene quantification, and statistical analysis using reproducible Conda environments.

## 1. Workflow Overview

The pipeline performs the following steps:

1. **FastQC** — Quality assessment of raw sequencing reads.
2. **fastp** — Read trimming and filtering.
3. **FastQC** — Quality assessment of trimmed reads.
4. **STAR** — Alignment to the reference genome.
5. **SAMtools** — BAM indexing.
6. **featureCounts** — Gene-level read quantification.
7. **Python (pandas)** — Combining individual gene counts into a count matrix.
8. **DESeq2** — Differential gene expression analysis.
9. **ggplot2** — PCA and volcano plot visualization.

## 2. Dataset

The workflow was tested using publicly available *Saccharomyces cerevisiae* RNA-seq data from GEO accession **GSE110004**.

The experiment compares two conditions:

- **Control:** Rap1-AID cells treated with DMSO.
- **Treatment:** Rap1-AID cells treated with IAA for 2 hours.

Three biological replicates were included per condition.

| Sample | Condition | Replicate |
|--------|-----------|-----------|
| SRR6357073 | Control | 1 |
| SRR6357074 | Control | 2 |
| SRR6357075 | Control | 3 |
| SRR6357079 | Treatment | 1 |
| SRR6357080 | Treatment | 2 |
| SRR6357081 | Treatment | 3 |

**Note:** For workflow development and educational testing, only the first 50,000 read pairs from each sample were processed. These subsets are not randomly sampled and are insufficient for reliable biological conclusions.

## 3. Reference Genome

- Organism: *Saccharomyces cerevisiae*
- Genome assembly: R64-1-1
- Reference genome: Ensembl DNA FASTA
- Gene annotation: Ensembl GTF, release 115

Gene quantification uses reverse-stranded counting (`featureCounts -s 2`), consistent with the library orientation verified during testing.

## 4. Software and Dependencies

The pipeline uses:

- Snakemake
- Conda
- FastQC
- fastp
- STAR
- SAMtools
- Subread / featureCounts
- Python and pandas
- R
- DESeq2
- ggplot2

Individual software environments are specified in the `envs/` directory.

## 5. Running the Pipeline

Install Snakemake and a compatible Conda installation.

Download the paired-end FASTQ files and the reference genome and annotation files.

The workflow currently expects sequencing data in:

`../../datasets/yeast-rnaseq/`

and reference files in:

`../../datasets/reference/`

These paths are relative to the project directory and may need to be adapted for other systems.

Check the workflow using a dry run:

```bash
snakemake -n --use-conda
```

Execute the pipeline:

```bash
snakemake --cores 2 --use-conda
```

Snakemake automatically determines which processing steps need to be executed.

## 6. Main Outputs

| Output | Description |
|--------|-------------|
| `fastqc_results/` | Raw read quality reports |
| `trimmed/` | Trimmed FASTQ files |
| `fastqc_trimmed/` | Trimmed read quality reports |
| `aligned/` | STAR-aligned BAM files |
| `counts/` | Individual gene counts and combined count matrix |
| `deseq2_results/differential_expression.tsv` | Differential expression statistics |
| `deseq2_results/PCA_plot.pdf` | Principal component analysis |
| `deseq2_results/volcano_plot.pdf` | Volcano plot |

Generated sequencing and analysis outputs are excluded from Git tracking using `.gitignore`.

## 7. Differential Expression Analysis

Differential expression is performed using DESeq2 with the design:

```r
design = ~ condition
```

The comparison is:

**Treatment versus Control**

The analysis estimates log2 fold changes and adjusts p-values for multiple testing.

PCA and volcano plots are generated to visualize sample relationships and differential expression patterns.

## 8. Limitations

This project demonstrates a functional RNA-seq analysis workflow but is not intended to establish publication-ready biological findings.

Important limitations include:

- Only 50,000 read pairs were analyzed per sample.
- Read subsets were taken from the beginning of each sequencing run rather than randomly sampled.
- Sequencing depth is insufficient for comprehensive transcriptome analysis.
- Further quality assessment and validation would be needed for research-grade conclusions.

## 9. Future Improvements

Potential extensions include:

- MultiQC reporting
- Gene Ontology enrichment analysis
- Additional quality-control summaries
- Full-depth RNA-seq analysis
- Improved portability and workflow configuration
- Automated testing and validation

## 10. Project Purpose

This project was developed to practice reproducible bioinformatics workflow development, RNA-seq processing, differential expression analysis, Conda environment management, and Git/GitHub version control.

