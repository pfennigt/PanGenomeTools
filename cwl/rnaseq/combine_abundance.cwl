#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

# Metadata
label: Combine Abundance Files
doc: |
  Combines all kallisto abundance TSV files into a single combined abundance file.
  Uses data.table if available, otherwise falls back to base R.

# Requirements
requirements:
  - class: DockerRequirement
    dockerPull: pangenometools-rnaseq-cwl
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.abundance_files)

# Inputs
inputs:
  abundance_files:
    type: File[]
    doc: List of kallisto abundance TSV files

# Outputs
outputs:
  combined_abundance:
    type: File
    outputBinding:
      glob: combined_abundance.tsv
    doc: Combined abundance TSV file

# Base command
baseCommand: [Rscript, /app/cwl/rnaseq/scripts/combine_abundance.R]
