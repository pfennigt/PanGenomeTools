#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Kallisto Index
doc: |
  Creates a kallisto index from a reference FASTA file.

# Requirements
requirements:
  - class: DockerRequirement
    dockerPull: pangenometools-rnaseq-cwl
  - class: ShellCommandRequirement

# Inputs
inputs:
  ref_fasta:
    type: File
    doc: Reference FASTA file to index

# Outputs
outputs:
  index:
    type: File
    outputBinding:
      glob: "ref_index"
    doc: Kallisto index file

# Standard output and error handling
stdout: kallisto_index.log
stderr: kallisto_index.error.log

# Base command
baseCommand: [kallisto, index]

# Arguments
arguments:
  - prefix: -i
    valueFrom: ref_index
  - valueFrom: $(inputs.ref_fasta.path)
    position: 1

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 1024