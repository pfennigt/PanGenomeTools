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
      - entry: $(inputs.fastqc_reports)
      - entry: $(inputs.fastqc_reports2)

# Inputs
inputs:
  fastqc_reports:
    type: File[]
    doc: FastQC reports of the RNASeq files
  fastqc_reports2:
    type: File[]?
    doc: Additional FastQC reports of the RNASeq files
  threads:
    type: int?
    doc: "Number of threads to use"
    default: 1

# Outputs
outputs:
  multiqc_report:
    type: File
    outputBinding:
      glob: "multiqc_report.html"
    doc: MultiQC report in zip format

# Standard output and error handling
stdout: multiqc.log
stderr: multiqc.error.log

# Base command
baseCommand: [multiqc, .]

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: 1024