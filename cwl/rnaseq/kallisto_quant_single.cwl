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
  fragment_length_mean:
    type: int
    doc: "Mean length of fragments in FASTQ file"
  fragment_length_sd:
    type: int
    doc: "Standard deviation of fragment lengths in FASTQ file"
  threads:
    type: int?
    default: 1
    doc: "Number of threads to use"

# Outputs
outputs:
  abundance_tsv:
    type: File
    outputBinding:
      glob: .out/abundance.tsv
      outputEval: ${self[0].basename=inputs.sample_name+".tsv"; return self;}
    doc: Kallisto abundance TSV file
  abundance_h5:
    type: File
    outputBinding:
      glob: .out/abundance.h5
      outputEval: ${self[0].basename=inputs.sample_name+".h5"; return self;}
    doc: Kallisto abundance HDF5 file

# Standard output and error handling
stdout: kallisto_quant.log
stderr: kallisto_quant.error.log

# Base command
baseCommand: bash

arguments:
  - -c
  - |
    kallisto quant \
      -i $(inputs.index.path) \
      -t $(inputs.threads) \
      --single \
      -l $(inputs.fragment_length_mean) -s $(inputs.fragment_length_sd) \
      -o .out <(cat $(inputs.fastq_file.join(' ')))

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: $(inputs.threads)
    ramMin: 1024
