#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "Homolog Finder"
doc: |
  Uses BLAST or DIAMOND to find Homologs for a Query FASTA.

  Searches for matches for each protein sequence in the query file using BLAST or DIAMOND.

# Base command - using Python module execution
baseCommand: [/usr/local/bin/_entrypoint.sh, python, -m, pangenometools.cli.find_homologs_cli]

# Input parameters
inputs: 
  # Query FASTA
  query_fasta:
    type: 
    - File
    - {"type": "array", "items": "File"}
    inputBinding:
      prefix: --query-fasta
    doc: "Path to query FASTA file (TF protein sequences)"

  # Target genome
  target_genome:
    type: File
    inputBinding:
      prefix: --target-genome
      valueFrom: $(self.basename)
    doc: "Genome fasta of the reference species to retrieve the gene sequences from."

  reference_gff:
    type: File?
    inputBinding:
      prefix: --target-gff
    doc: "GFF annotation of the reference genome to retrieve the gene sequences from."

  # Pre-created database path
  db_path:
    type: Directory?
    inputBinding:
      prefix: --db-path
    doc: "Path to the BLAST/DIAMOND database"

  # Method
  method:
    type: string
    inputBinding:
      prefix: --method
    default: "blast"
    doc: "Search method - 'blast' or 'diamond' (default: blast)"

  # Search parameters
  evalue:
    type: float?
    inputBinding:
      prefix: --evalue
    default: 1e-5
    doc: "E-value threshold (default: 1e-5)"

  max_target_seqs:
    type: int?
    inputBinding:
      prefix: --max-target-seqs
    default: 10
    doc: "Maximum number of target sequences to report (default: 10)"
  
  num_threads:
    type: int?
    inputBinding:
      prefix: --num-threads
    default: 1
    doc: "Number of threads to use (default: 1)"

  # BLAST options
  blast_type:
    type: string?
    inputBinding:
      prefix: "--blast-type"
    default: "blastp"
    doc: "BLAST type (default: blastp)"

  # DIAMOND options
  diamond_sensitivity:
    type: string?
    inputBinding:
      prefix: --diamond-sensitivity
    default: "sensitive"
    doc: "DIAMOND sensitivity level (default: sensitive)"

  # General options
  keep_database:
    type: boolean
    inputBinding:
      prefix: --keep-database
    default: true
    doc: "Keep the database file after search"

# Output files
outputs:
  search_results:
    type: File
    outputBinding:
      glob: "*.tsv"
    doc: "Directory of the created database"

# Standard output and error handling
stdout: find_homologs.log
stderr: find_homologs.error.log

# Requirements
requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: "pangenometools-cwl-blast"
  - class: InitialWorkDirRequirement
    listing:
      - $(inputs.target_genome)

# Hints for better performance
hints:
  - class: ResourceRequirement
    coresMin: 1
    ramMin: 1024