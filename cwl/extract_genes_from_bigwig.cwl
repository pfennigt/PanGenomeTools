#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

# Metadata
label: Extract genes from BigWig using PanGenomeTools
doc: |
  Extracts signal values from BigWig files using GFF coordinates.
  This tool uses the PanGenomeTools BigWigHandler for efficient signal extraction.

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

  bigwig_folder:
    type: Directory
    inputBinding:
      position: 3
      prefix: --bigwig-folder
    doc: Path to BigWig folder

  bigwig_index:
    type: File
    inputBinding:
      position: 4
      prefix: --bigwig-index
    doc: Path to BigWig index CSV file

  target_genes:
    type: File
    inputBinding:
      position: 5
      prefix: --target-genes
    doc: Path to target genes CSV file

  output:
    type: string
    inputBinding:
      position: 6
      prefix: --output
    doc: Output JSON file path for extracted signals

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
    doc: Number of null values to pad between segments

  whole_seq:
    type: boolean?
    inputBinding:
      prefix: --whole-seq
    default: false
    doc: Extract entire feature sequence

  use_five_prime_direction:
    type: boolean?
    inputBinding:
      prefix: --use-five-prime-direction
    default: false
    doc: Always interpret upstream/downstream in 5' direction regardless of strand

  silent:
    type: boolean?
    inputBinding:
      prefix: --silent
    default: false
    doc: Suppress progress output

# Outputs
outputs:
  extracted_signals:
    type: File
    outputBinding:
      glob: $(inputs.output)
    doc: Extracted BigWig signals in JSON format

# Base command - use Python to call the function directly
baseCommand: [/usr/local/bin/_entrypoint.sh, python]

# Inline Python script that calls the function directly
stdin: |
  import sys
  from argparse import Namespace
  from pangenometools.cli.bigwig_cli import extract_bigwig_signals

  # Create arguments namespace
  args = Namespace(
      pangenome_folder=$(inputs.pangenome_folder.path),
      pangenome_index=$(inputs.pangenome_index.path),
      bigwig_folder=$(inputs.bigwig_folder.path),
      bigwig_index=$(inputs.bigwig_index.path),
      target_genes=$(inputs.target_genes.path),
      output=$(inputs.output),
      feature_type=$(inputs.feature_type),
      upstream=$(inputs.upstream),
      downstream=$(inputs.downstream),
      inner_start=$(inputs.inner_start),
      inner_end=$(inputs.inner_end),
      pad=$(inputs.pad),
      whole_seq=$(inputs.whole_seq),
      use_five_prime_direction=$(inputs.use_five_prime_direction),
      silent=$(inputs.silent)
  )

  # Call the function directly
  extract_bigwig_signals(args)