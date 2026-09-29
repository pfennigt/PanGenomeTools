#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

# Metadata
label: Kallisto Quantification for Paired Reads
doc: |
  Runs kallisto quantification for paired-end sequencing data.

# Requirements
requirements:
  - class: DockerRequirement
    dockerPull: pangenometools-rnaseq-cwl
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.ref_dir)

# Inputs
inputs:
  sample_name:
    type: string
    doc: Name of the sample
  index:
    type: File
    doc: Kallisto index file
  ref_dir:
    type: Directory
    doc: Reference directory for the fastq_file paths
  fastq_file:
    type: string[]
    doc: First FASTQ file (R1)

# Outputs
outputs:
  abundance_tsv:
    type: File
    outputBinding:
      glob: .out/abundance.tsv
      outputEval: ${self[0].basename=inputs.sample_name+".tsv"; return self;}
    doc: Kallisto abundance TSV file

# Base command
baseCommand: bash

arguments:
  - -c
  - |
    kallisto quant \
      -i $(inputs.index.path) \
      --single \
      -l $(inputs.fragment_length_mean) -s $(inputs.fragment_length_sd) \
      -o .out <(cat $(inputs.fastq_file.join(' ')))

