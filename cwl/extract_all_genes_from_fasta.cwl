#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Extract all genes from FASTA using AGAT
doc: |
  Extracts all genes from FASTA files using AGAT for multiple genomes.
  This tool runs AGAT transcriptome extraction for each genome listed in the index.

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.output)
        writable: true

# Inputs
inputs:
  index:
    type: File
    inputBinding:
      position: 1
      prefix: --index
    doc: Path to pangenome_index.csv file

  output:
    type: string
    inputBinding:
      position: 2
      prefix: --output
    doc: Output directory for extracted FASTA files

  pangenome_dir:
    type: string?
    inputBinding:
      prefix: --pangenome-dir
    default: ""
    doc: Base directory for genome and annotation paths

  type:
    type: string?
    inputBinding:
      prefix: --type
    default: "gene"
    doc: Type of feature to extract

  extra:
    type: string?
    inputBinding:
      prefix: --extra
    default: ""
    doc: Additional arguments to pass to AGAT

# Outputs
outputs:
  extracted_fasta_files:
    type: File[]
    outputBinding:
      glob: $(inputs.output)/*.fa
    doc: Extracted FASTA files for each genome

  log_files:
    type: File[]
    outputBinding:
      glob: $(inputs.output)/logs/*.log
    doc: Log files from AGAT execution

# Base command - call the script directly
baseCommand: [/usr/local/bin/_entrypoint.sh, python]

# Arguments
arguments:
  - valueFrom: workflows/PanGenomeTools/scripts/extract_all_genes_from_fasta.py
  - prefix: --index
    valueFrom: $(inputs.index.path)
  - prefix: --output
    valueFrom: $(inputs.output)
  - prefix: --pangenome-dir
    valueFrom: $(inputs.pangenome_dir)
  - prefix: --type
    valueFrom: $(inputs.type)
  - prefix: --extra
    valueFrom: $(inputs.extra)