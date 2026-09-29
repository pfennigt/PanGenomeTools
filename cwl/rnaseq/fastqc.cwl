#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

# Metadata
label: FastQC - Quality Control Tool for Single Reads
doc: |
  Runs FastQC to assess the quality of single-end sequencing data.

# Requirements
requirements:
  - class: DockerRequirement
    dockerPull: pangenometools-rnaseq-cwl
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.ref_dir)

# Inputs
inputs:
  fastq_file:
    type: string
    doc: Single-end FASTQ file to analyze
  ref_dir:
    type: Directory
    doc: Reference directory for the fastq_file paths
  threads:
    type: int?
    doc: "Number of threads to use"
    default: 1
  output_dir:
    type: string?
    doc: "Directory for output files"
    default: "."

# Outputs
outputs:
  fastqc_report:
    type: File
    outputBinding:
      glob: "*_fastqc.zip"
    doc: FastQC report in zip format

# Standard output and error handling
stdout: fastqc.log
stderr: fastqc.error.log

# Base command
baseCommand: bash

arguments:
  - -c
  - |
    fastqc \
      -t $(inputs.threads) \
      -o $(inputs.output_dir) \
      $(inputs.fastq_file)

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: 1024