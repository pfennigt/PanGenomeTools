#!/usr/bin/env cwl-runner
cwlVersion: v1.2
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
    type: File?
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

  assays_dir:
    type: Directory
    doc: "Path to the assays directory containing input files"

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

  per_gene_group:
    type: boolean?
    inputBinding:
      prefix: --per-gene-group
    default: false
    doc: Write sequences into files per gene group, not per genotype

  genotypes:
    type:
    - "null"
    - string
    - type: array
      items: string 
    inputBinding:
      prefix: --genotypes
    doc: Genotype(s) to run the analysis for (all by default)

  silent:
    type: boolean?
    inputBinding:
      prefix: --silent
    default: false
    doc: Suppress progress output

# Outputs
outputs:
  extracted_sequences:
    type: 
    - File
    - type: array
      items: File
    outputBinding:
      glob: "$(inputs.output)/*.fa"
    doc: Extracted sequences in FASTA format

# Base command - use Python to call the function directly
baseCommand: [python, "-c", "from pangenometools.cli.fasta_cli import fasta_cli; fasta_cli()"]

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl:latest"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 10000

# Standard output and error handling
stdout: extract_genes_from_fasta.log
stderr: extract_genes_from_fasta.error.log
