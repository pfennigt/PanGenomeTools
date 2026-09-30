# RNA-Seq Kallisto Mapping Workflows

This directory contains CWL workflows for processing RNA-Seq data by mapping reads to a reference transcriptome using kallisto.

## Features

- Processes RNA-Seq data from local FASTQ files using a CSV index
- Performs quality control with FastQC
- Maps reads to reference using kallisto
- Combines abundance estimates from all samples into a single TSV file
- Extracts CDS sequences from a genome in a pangenome
- Separate workflows for single-end and paired-end reads

## Available Workflows

- `rnaseq_single_workflow.cwl`: For single-end reads only
- `rnaseq_paired_workflow.cwl`: For paired-end reads only

## Requirements

- CWL runner (e.g., cwltool, toil, cwl-runner)
- Docker or Singularity/Apptainer with the `pangenometools-cwl` and `pangenometools-rnaseq-cwl` containers
- Pangenome folder and index files

## Quick Start

1. Prepare your CSV index file with sample information (see USAGE_GUIDE.md for format)
2. Adapt the `run.yml` file with the correct file paths and genotype to use
3. Run the appropriate workflow:

```bash
cwltool --parallel --log-dir logs --outdir dataset --singularity rnaseq_paired_workflow.cwl run.yml
```
