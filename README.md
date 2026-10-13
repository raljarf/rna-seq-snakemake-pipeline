# RNA-seq Analysis Pipeline Using Snakemake

A reproducible RNA-seq processing pipeline built using Snakemake.

## Workflow

1. Raw read quality control - FastQC
2. Read trimming - fastp
3. Quality control after trimming - FastQC
4. Genome alignment - STAR
5. BAM indexing - SAMtools
6. Gene quantification - featureCounts
7. Gene count matrix generation - Python and pandas

## Technologies

- Snakemake
- Python
- Conda
- Linux (WSL2)
- Git and GitHub

## Project Status

The RNA-seq processing workflow is functional and has been tested
using a small practice dataset.

Future development will include analysis of real biological samples
and differential gene expression analysis.
