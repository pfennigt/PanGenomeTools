#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Extract genes from FASTA using PanGenomeTools
doc: |
  Extracts sequences from FASTA files using GFF coordinates.
  This tool uses the PanGenomeTools FastaHandler for efficient sequence extraction.

# Inputs
inputs:
  pangenome_folder:
    type: Directory
    inputBinding:
      position: 1
      prefix: --pangenome-folder
    doc: Path to pangenome folder containing genome assemblies

  pangenome_index:
    type: File
    inputBinding:
      position: 2
      prefix: --pangenome-index
    doc: Path to pangenome index file

  target_genes:
    type: File
    inputBinding:
      position: 3
      prefix: --target-genes
    doc: Path to target genes CSV file

  output:
    type: string
    inputBinding:
      position: 4
      prefix: --output
    doc: Output file path for extracted sequences

  # Optional parameters
  feature_type:
    type: string?
    inputBinding:
      prefix: --feature-type
    default: "gene"
    doc: Type of feature to use for coordinates

  upstream:
    type: int?
    inputBinding:
      prefix: --upstream
    default: 0
    doc: Nucleotides to include upstream (5' direction)

  downstream:
    type: int?
    inputBinding:
      prefix: --downstream
    default: 0
    doc: Nucleotides to include downstream (3' direction)

  inner_start:
    type: int?
    inputBinding:
      prefix: --inner-start
    default: 0
    doc: Nucleotides to include from start of feature

  inner_end:
    type: int?
    inputBinding:
      prefix: --inner-end
    default: 0
    doc: Nucleotides to include from end of feature

  pad:
    type: int?
    inputBinding:
      prefix: --pad
    default: 0
    doc: Number of Ns to pad between segments

  merge_strategy:
    type:
      - "null"
      - type: enum
        symbols: ["merge", "first", "all"]
    inputBinding:
      prefix: --merge-strategy
    default: "merge"
    doc: Strategy for merging multiple features

  use_five_prime_direction:
    type: boolean?
    inputBinding:
      prefix: --use-five-prime-direction
    default: false
    doc: Always interpret upstream/downstream in 5' direction regardless of strand

  include_coordinates:
    type: boolean?
    inputBinding:
      prefix: --include-coordinates
    default: false
    doc: Include coordinates in output header

  silent:
    type: boolean?
    inputBinding:
      prefix: --silent
    default: false
    doc: Suppress progress output

# Outputs
outputs:
  extracted_sequences:
    type: File
    outputBinding:
      glob: $(inputs.output)
    doc: Extracted sequences in FASTA format

# Base command - use Python to call the function directly
baseCommand: python

# Inline Python script that calls the function directly
stdin: |
  import sys
  from argparse import Namespace
  from pangenometools.cli.fasta_cli import extract_fasta_sequences

  # Create arguments namespace
  args = Namespace(
      pangenome_folder=$(inputs.pangenome_folder.path),
      pangenome_index=$(inputs.pangenome_index.path),
      target_genes=$(inputs.target_genes.path),
      output=$(inputs.output),
      feature_type=$(inputs.feature_type),
      upstream=$(inputs.upstream),
      downstream=$(inputs.downstream),
      inner_start=$(inputs.inner_start),
      inner_end=$(inputs.inner_end),
      pad=$(inputs.pad),
      merge_strategy=$(inputs.merge_strategy),
      use_five_prime_direction=$(inputs.use_five_prime_direction),
      include_coordinates=$(inputs.include_coordinates),
      silent=$(inputs.silent)
  )

  # Call the function directly
  extract_fasta_sequences(args)