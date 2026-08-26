#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "BLAST or DIAMOND Database Creator"
doc: |
  Creates a database for BLAST or DIAMOND searches.

  If a genome GFF file is provided, the FASTA is assumed to contain nucleotide sequences and
  AGAT is used to create translated protein sequences for the database creation

# Base command - using Python module execution
baseCommand: [/usr/local/bin/_entrypoint.sh, python, -m, pangenometools.cli.make_blast_diamond_db_cli]

# Input parameters
inputs: 
  # Target genome
  target_genome:
    type: File
    inputBinding:
      prefix: --target-genome
    doc: "Genome fasta of the target species to retrieve the gene sequences from."

  target_gff:
    type: File
    inputBinding:
      prefix: --target-gff
    doc: "GFF annotation of the target genome to retrieve the gene sequences from."

  # Method
  method:
    type: string
    inputBinding:
      prefix: --method
    default: "blast"
    doc: "Search method - 'blast' or 'diamond' (default: blast)"

# Output files
outputs:
  database:
    type: Directory
    outputBinding:
      glob: "db"
    doc: "Directory of the created database"

# Standard output and error handling
stdout: make_blast_diamond_db.log
stderr: make_blast_diamond_db.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl-blast"

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 1024