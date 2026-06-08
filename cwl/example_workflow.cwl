#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

# Metadata
label: PanGenomeTools Example Workflow
doc: |
  Example workflow demonstrating how to combine PanGenomeTools CWL tools.
  This workflow extracts GFF features and FASTA sequences for target genes.

# Inputs
inputs:
  pangenome_folder:
    type: Directory
    doc: Path to pangenome folder containing genome assemblies

  pangenome_index:
    type: File
    doc: Path to pangenome index file

  target_genes:
    type: File
    doc: Path to target genes CSV file

  # Optional parameters
  feature_type:
    type: string?
    default: "gene"
    doc: Type of feature to extract

  upstream:
    type: int?
    default: 1000
    doc: Nucleotides to include upstream (5' direction)

  downstream:
    type: int?
    default: 500
    doc: Nucleotides to include downstream (3' direction)

  gff_output:
    type: string?
    default: "extracted_features.gff"
    doc: Output filename for GFF features

  fasta_output:
    type: string?
    default: "extracted_sequences.fasta"
    doc: Output filename for FASTA sequences

# Steps
steps:
  extract_gff:
    run: extract_genes_from_gff.cwl
    in:
      pangenome_folder: pangenome_folder
      pangenome_index: pangenome_index
      target_genes: target_genes
      feature_type: feature_type
      output: gff_output
    out: [extracted_features]
  extract_fasta:
    run: extract_genes_from_fasta.cwl
    in:
      pangenome_folder: pangenome_folder
      pangenome_index: pangenome_index
      target_genes: target_genes
      feature_type: feature_type
      upstream: upstream
      downstream: downstream
      output: fasta_output
    out: [extracted_sequences]
# Outputs
outputs:
  gff_features:
    type: File
    outputSource: extract_gff/extracted_features
    doc: Extracted GFF features for target genes

  fasta_sequences:
    type: File
    outputSource: extract_fasta/extracted_sequences
    doc: Extracted FASTA sequences for target genes