#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Extract genes from GFF using PanGenomeTools
doc: |
  Extracts GFF features for target genes using PanGenomeTools.
  This tool uses the PanGenomeTools GFFHandler for efficient feature extraction.

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
    doc: Output file path for extracted GFF features

  # Optional parameters
  feature_type:
    type: string?
    inputBinding:
      prefix: --feature-type
    default: "gene"
    doc: Type of feature to extract

  search_mode:
    type:
      - "null"
      - type: enum
        symbols: ["strict", "children", "pattern"]
    inputBinding:
      prefix: --search-mode
    default: "children"
    doc: Search mode for finding features

  merge_strategy:
    type:
      - "null"
      - type: enum
        symbols: ["merge", "first", "all"]
    inputBinding:
      prefix: --merge-strategy
    default: "merge"
    doc: Strategy for merging multiple features

  silent:
    type: boolean?
    inputBinding:
      prefix: --silent
    default: false
    doc: Suppress progress output

# Outputs
outputs:
  extracted_features:
    type: File
    outputBinding:
      glob: $(inputs.output)
    doc: Extracted GFF features

# Base command - use Python to call the function directly
baseCommand: [/usr/local/bin/_entrypoint.sh, python]

# Inline Python script that calls the function directly
stdin: |
  import sys
  from argparse import Namespace
  from pangenometools.cli.gff_cli import extract_gff_features

  # Create arguments namespace
  args = Namespace(
      pangenome_folder=$(inputs.pangenome_folder.path),
      pangenome_index=$(inputs.pangenome_index.path),
      target_genes=$(inputs.target_genes.path),
      output=$(inputs.output),
      feature_type=$(inputs.feature_type),
      search_mode=$(inputs.search_mode),
      merge_strategy=$(inputs.merge_strategy),
      silent=$(inputs.silent)
  )

  # Call the function directly
  extract_gff_features(args)