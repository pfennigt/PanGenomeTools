#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Create pseudo GFF annotations from FASTA
doc: |
  Creates pseudo GFF annotations from a FASTA file.
  This is useful for creating minimal annotations when none are available.

# Inputs
inputs:
  input_fasta:
    type: File
    inputBinding:
      position: 1
    doc: Path to input FASTA file

  output_gff:
    type: string
    inputBinding:
      position: 2
    doc: Path to output GFF file

  feature_type:
    type: string?
    inputBinding:
      prefix: --feature_type
    default: "gene"
    doc: Feature type for GFF entries

  upstream:
    type: int?
    inputBinding:
      prefix: --upstream
    default: 1500
    doc: Number of nucleotides extracted upstream of the feature

  downstream:
    type: int?
    inputBinding:
      prefix: --downstream
    default: 1500
    doc: Number of nucleotides extracted downstream of the feature

# Outputs
outputs:
  pseudo_annotation:
    type: File
    outputBinding:
      glob: $(inputs.output_gff)
    doc: Generated pseudo GFF annotation file

# Base command - call the script directly
baseCommand: python

# Arguments
arguments:
  - valueFrom: $(inputs.input_fasta.path)
  - valueFrom: $(inputs.output_gff)
  - prefix: --feature_type
    valueFrom: $(inputs.feature_type)
  - prefix: --upstream
    valueFrom: $(inputs.upstream)
  - prefix: --downstream
    valueFrom: $(inputs.downstream)
  - valueFrom: workflows/PanGenomeTools/scripts/make_pseudo_annotation.py