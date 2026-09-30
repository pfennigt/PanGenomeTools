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

  # Optional parameters
  feature_type:
    type: string?
    inputBinding:
      prefix: --feature-type
    default: "gene"
    doc: Type of feature to use for coordinates

  genotypes:
    type:
    - "null"
    - string
    - type: array
      items: string 
    inputBinding:
      prefix: --genotypes
    doc: Genotype(s) to run the analysis for (all by default)

# Outputs
outputs:
  extracted_sequences:
    type: 
    - File
    - type: array
      items: File
    outputBinding:
      glob: "*.fa"
    doc: Extracted sequences in FASTA format

# Base command - use Python to call the function directly
baseCommand: [/usr/local/bin/_entrypoint.sh, python, "-c", "from pangenometools.cli.fasta_agat_cli import main; main()"]

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 10000

# Standard output and error handling
stdout: extract_from_fasta_agat.log
stderr: extract_from_fasta_agat.error.log
