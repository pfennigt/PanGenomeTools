#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

# Metadata
label: RNA-Seq Paired-End Workflow
doc: |
  RNA-Seq workflow for paired-end reads that processes data from a CSV index,
  performs quality control, and maps to reference using kallisto.

# Inputs
inputs:
  # Reference inputs
  pangenome_index_path:
    type: File
    doc: "Path to pangenome index CSV file"
  pangenome_path:
    type: Directory
    doc: "Directory containing pangenome assemblies"
  rnaseq_folder:
    type: Directory
    doc: "Directory containing RNA-Seq FASTQ files"
  rnaseq_index:
    type: File
    doc: "CSV file with columns: sample, read1, read2"
  genotype:
    type: string
    doc: "Pan-genome genotype to map against"

  # Fragment length distribution
  fragment_length_mean:
    type: int
    doc: "Mean length of fragments in FASTQ file"
  fragment_length_sd:
    type: int
    doc: "Standard deviation of fragment lengths in FASTQ file"

  # Settings for runs
  threads_fastqc:
    type: int
    default: 1
    doc: "Number of threads to use for each FASTQC process"
  threads_kallisto:
    type: int
    default: 1
    doc: "Number of threads to use for each Kallisto process"
  batch_size_multiqc:
    type: int
    default: 100
    doc: "Number of FastQC files that are combined in one MultiQC report"

# Outputs
outputs:
  sample_abundances:
    type: Directory
    outputSource: collect_kallisto_files/out
    doc: Abundance HDF5 files with kallisto estimates for each sample
  combined_abundance:
    type: File
    outputSource: combine_step/combined_abundance
    doc: Combined abundance TSV file with kallisto estimates for all samples
  multiqc_report:
    type: File[]
    outputSource: multiqc/multiqc_report
    doc: MultiQC report of the RNA-seq files

# Steps
steps:
  ################################################################################
  #                                    FASTQC                                    #
  ################################################################################

  # Parse the rnaseq_index CSV file
  parse_csv:
    run: rnaseq_single_parse_csv.cwl
    in:
      csv: rnaseq_index
    out: [sample, read]

  fastqc:
    run: fastqc.cwl
    in: 
      ref_dir: rnaseq_folder
      fastq_file: parse_csv/read
      threads: threads_fastqc
    out: [fastqc_report]
    scatter: fastq_file

  batch_fastqc_results:
    run: ../batch_files.cwl
    in:
      files: fastqc/fastqc_report
      batch_size: batch_size_multiqc
    out: [batches]

  multiqc:
    run: multiqc.cwl
    in:
      fastqc_reports: batch_fastqc_results/batches
    scatter: fastqc_reports
    out: [multiqc_report]

  ################################################################################
  #                                   Kallisto                                   #
  ################################################################################

  parse_csv_grouped:
    run: rnaseq_single_parse_csv_grouped.cwl
    in:
      csv: rnaseq_index
      groupby: {default: "sample"}
    out: [sample, read]

  # Extract genes from pangenome
  extract_genes:
    run: ../extract_from_fasta_agat.cwl
    in:
      pangenome_index: pangenome_index_path
      feature_type: {default: "CDS"}
      pangenome_folder: pangenome_path
      genotypes: genotype
    out: [extracted_sequences]

  # Create kallisto index from extracted genes
  kallisto_index:
    run: kallisto_index.cwl
    in:
      ref_fasta: extract_genes/extracted_sequences
    out: [index]

  # Perform Kallisto mapping
  kallisto:
    run: kallisto_quant_single.cwl
    in: 
      sample_name: parse_csv_grouped/sample
      index: kallisto_index/index 
      ref_dir: rnaseq_folder
      fastq_file: parse_csv_grouped/read
      fragment_length_mean: fragment_length_mean
      fragment_length_sd: fragment_length_sd
      threads: threads_kallisto
    scatter: [sample_name, fastq_file]
    scatterMethod: dotproduct
    out: [abundance_tsv, abundance_h5]

  # Collect the kallisto files in a directory
  collect_kallisto_files:
    run: ../files_to_dir.cwl
    in:
      files: kallisto/abundance_h5
      dir_name: {default: "abundances"}
    out: [out]

  # Combine all abundance files
  combine_step:
    run: combine_abundance.cwl
    in:
      abundance_files: kallisto/abundance_tsv
    out: [combined_abundance]

requirements:
  - class: ScatterFeatureRequirement
