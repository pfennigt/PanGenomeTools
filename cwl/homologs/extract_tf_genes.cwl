#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# Metadata for the tool
label: "TF Gene Extractor"
doc: |
  Processes an Excel file to extract the names of TF genes in a given column.
  
  Allows to subset the table by a gene family column as well as through queries.

# Base command - using Python module execution
baseCommand: [/usr/local/bin/_entrypoint.sh, python, -m, pangenometools.cli.extract_tf_genes_cli]

# Input parameters
inputs:
  excel_file:
    type: File
    inputBinding:
      prefix: --excel-file
      position: 1
    doc: "Excel file containing the list of TF genes"
  
  output_dir:
    type: string?
    inputBinding:
      prefix: --output-dir
    default: "."
    doc: "Directory for output files (one per TF family)"

  # Excel file specifications
  tf_family_column:
    type: string?
    inputBinding:
      prefix: --tf-family-column
    default: "TF family"
    doc: "Name of column containing TF family names (default: 'TF family')"

  gene_name_column:
    type: string?
    inputBinding:
      prefix: --gene-name-column
    default: "TF Gene"
    doc: "Name of column containing gene names (default: 'TF Gene')"

  sheet_name:
    type: string?
    inputBinding:
      prefix: --sheet-name
    doc: "Name of Excel sheet to read (default: first sheet)"

  families:
    type: string[]
    inputBinding:
      prefix: --families
    doc: "Specific TF families to extract (default: all families)"

# Output files
outputs:
  gene_names:
    type: File[]
    outputBinding:
      glob: "*.csv"
    doc: "Files with the gene names in the regerence species (one per family)"

# Standard output and error handling
stdout: extract_tf_genes.log
stderr: extract_tf_genes.error.log

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