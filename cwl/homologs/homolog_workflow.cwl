#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

# Metadata
label: "Homolog discovery Workflow"
doc: |
  Complete workflow for Homolog discovery from a list of gene names.
  This workflow prepares extracts the protein sequences of target genes from a reference genome and uses BLAST or DIAMOND to find homologs.

# Inputs
inputs:
  # Excel file (DeepCis)
  excel_file:
    type: File
    doc: "Excel file containing the list of TF genes"

  # Excel file specifications
  tf_family_column:
    type: string?
    default: "TF family"
    doc: "Name of column containing TF family names (default: 'TF family')"

  gene_name_column:
    type: string?
    default: "TF Gene"
    doc: "Name of column containing gene names (default: 'TF Gene')"

  sheet_name:
    type: string?
    doc: "Name of Excel sheet to read (default: first sheet)"

  families:
    type: string[]
    doc: "Specific TF families to extract (default: all families)"

  # Target genome
  target_genome:
    type: File
    doc: "Genome fasta of the target species to retrieve the gene sequences from."

  target_gff:
    type: File
    doc: "GFF annotation of the target genome to retrieve the gene sequences from."

  # Reference genome
  reference_genome:
    type: File
    doc: "Genome fasta of the reference species to retrieve the gene sequences from."

  reference_gff:
    type: File
    doc: "GFF annotation of the reference genome to retrieve the gene sequences from."

  # Method
  method:
    type: string
    default: "blast"
    doc: "Search method - 'blast' or 'diamond' (default: blast)"

  # Search parameters
  evalue:
    type: float?
    default: 1e-5
    doc: "E-value threshold (default: 1e-5)"

  max_target_seqs:
    type: int?
    default: 10
    doc: "Maximum number of target sequences to report (default: 10)"
  
  num_threads:
    type: int?
    default: 1
    doc: "Number of threads to use (default: 1)"

  # BLAST options
  blast_type:
    type: string?
    default: "blastp"
    doc: "BLAST type (default: blastp)"

  # DIAMOND options
  diamond_sensitivity:
    type: string?
    default: "sensitive"
    doc: "DIAMOND sensitivity level (default: sensitive)"

# Steps
steps:
  # Get the gene names for each TF family
  extract_tf_genes:
    run: extract_tf_genes.cwl
    in:
      excel_file: excel_file
      tf_family_column: tf_family_column
      gene_name_column: gene_name_column
      sheet_name: sheet_name
      families: families
    out: [gene_names]

  # Extract sequences
  extract_sequences:
    run: retrieve_sequences.cwl
    in:
      gene_list: extract_tf_genes/gene_names
      reference_genome: reference_genome
      reference_gff: reference_gff
    out: [translated_sequences]
    scatter: gene_list

  # Make the search database
  make_database:
    run: make_blast_diamond_db.cwl
    in:
      target_genome: target_genome
      target_gff: target_gff
      method: method
    out: [database]

  # Find the homologs
  find_homologs:
    run: find_homologs.cwl
    in:
      query_fasta: extract_sequences/translated_sequences
      target_genome: target_genome
      reference_gff: reference_gff
      db_path: make_database/database
      method: method
      evalue: evalue
      max_target_seqs: max_target_seqs
      num_threads: num_threads
      blast_type: blast_type
      diamond_sensitivity: diamond_sensitivity
    scatter: query_fasta
    out: [search_results]

# Outputs
outputs:
  search_results:
    type: File[]
    outputSource: find_homologs/search_results
    doc: "Found potential Homologs (TSV files)"

# Requirements
requirements:
  ScatterFeatureRequirement: {}
