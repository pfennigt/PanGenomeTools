#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "Gene Sequence Retriever"
doc: |
  Retrieves protein sequences from a genome using gene names in a CSV file.

  Subsets the genome GFF file to the requested genes and uses AGAT to get
  the protein sequences

# Base command - using Python module execution
baseCommand: [/usr/local/bin/_entrypoint.sh, python, -m, pangenometools.cli.retrieve_sequences_cli]

# Input parameters
inputs:
  gene_list:
    type: File
    inputBinding:
      prefix: --gene-list
      position: 1
    doc: "CSV file containing the names of genes to retrieve."
  
  # Reference genome
  reference_genome:
    type: File
    inputBinding:
      prefix: --reference-genome
    doc: "Genome fasta of the reference species to retrieve the gene sequences from."

  reference_gff:
    type: File
    inputBinding:
      prefix: --reference-gff
    doc: "GFF annotation of the reference genome to retrieve the gene sequences from."

# Output files
outputs:
  translated_sequences:
    type: File[]
    outputBinding:
      glob: "*.fa"
    doc: "FASTA file with the gene sequences translated to protein (one per family)"

# Standard output and error handling
stdout: retrieve_sequences.log
stderr: retrieve_sequences.error.log

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
    ramMin: 1024