"""
GFF CLI module for PanGenomeTools.

This module provides command-line interface functionality for GFF operations.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List
from ..models.gff import GFFHandler
from ..models.base import PangenomeFileHandler
import pandas as pd

def setup_gff_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for GFF extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract GFF features for target genes."
    )

    # Required arguments
    parser.add_argument("--pangenome-folder", required=True,
                       help="Path to pangenome folder")
    parser.add_argument("--pangenome-index", required=True,
                       help="Path to pangenome index file")
    parser.add_argument("--target-genes", required=True,
                       help="Path to target genes CSV file")
    parser.add_argument("--output", required=True,
                       help="Output file path")

    # Feature extraction options
    parser.add_argument("--feature-type", default="gene",
                       help="Type of feature to extract (default: gene)")
    parser.add_argument("--search-mode", default="children",
                       choices=["strict", "children", "pattern"],
                       help="Search mode for finding features (default: children)")
    parser.add_argument("--merge-strategy", default="merge",
                       choices=["merge", "first", "all"],
                       help="Strategy for merging multiple features (default: merge)")

    # Additional options
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser

def read_target_genes(target_path: Path) -> tuple:
    """
    Read target genes from CSV file.

    Args:
        target_path: Path to target genes CSV file

    Returns:
        Tuple of (genotypes, rows)
    """
    rows = []
    with open(target_path, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        geno_cols = [h for h in headers if h.startswith("gene_ID_")]
        if "gene_name" not in headers:
            raise ValueError("target_genes.csv must contain a gene_name column")
        for r in reader:
            rows.append(r)

    genotypes = [h.replace("gene_ID_", "") for h in geno_cols]
    return genotypes, rows

def extract_gff_features(args: argparse.Namespace) -> None:
    """
    Extract GFF features for target genes.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_folder = Path(args.pangenome_folder)
    pangenome_index = Path(args.pangenome_index)

    gff_handler = GFFHandler(pangenome_folder, pangenome_index)

    # Read target genes
    genotypes, target_rows = read_target_genes(Path(args.target_genes))

    # Open output file
    with open(args.output, "w") as out_fh:
        # Process each target row
        for row in target_rows:
            gene_name = row.get("gene_name", "")

            # Process each genotype
            for genotype in genotypes:
                gene_id_col = f"gene_ID_{genotype}"
                gene_id_raw = row.get(gene_id_col, "")

                if not gene_id_raw or gene_id_raw in ["", "[]", "[""]", "['']"]:
                    continue

                # Handle list format
                if gene_id_raw.startswith("["):
                    gene_id_raw = gene_id_raw[1:-1].split(",")
                else:
                    gene_id_raw = [gene_id_raw]

                # Process each gene ID
                for gene_id_raw_item in gene_id_raw:
                    gene_id = gene_id_raw_item.strip().strip("'\"")

                    if not gene_id:
                        continue

                    try:
                        # Get features for the gene
                        features = gff_handler.get_features_for_gene(
                            genotype, gene_id, args.feature_type, args.merge_strategy
                        )
                        if features is None or len(features) == 0:
                            if not args.silent:
                                print(f"Warning: No features found for {gene_id} in {genotype}", file=sys.stderr)
                            continue
                        
                        # If a single feature is returned, package it in a list
                        if isinstance(features, pd.DataFrame):
                            features = [features]

                        # Write features to output
                        for feature in features:
                            if isinstance(features, pd.DataFrame):
                                # It's a PyRanges object
                                for _, row in feature.iterrows():
                                    line = f"{row['Chromosome']}\t.\t{row['Feature']}\t{row['Start'] + 1}\t{row['End']}\t.\t{row['Strand']}\t.\tID={gene_id};gene_name={gene_name}\n"
                                    out_fh.write(line)
                            else:
                                # It's a single feature tuple
                                chrom, start, end, strand, attr = feature
                                line = f"{chrom}\t.\t{args.feature_type}\t{start + 1}\t{end}\t.\t{strand}\t.\tID={gene_id};gene_name={gene_name}\n"
                                out_fh.write(line)

                    except Exception as e:
                        if not args.silent:
                            print(f"Error processing {gene_id} in {genotype}: {e}", file=sys.stderr)
                        continue

def gff_cli():
    """
    Main entry point for GFF CLI.

    Parses arguments and executes GFF extraction.
    """
    parser = setup_gff_parser()
    args = parser.parse_args()

    try:
        extract_gff_features(args)
        if not args.silent:
            print(f"GFF extraction completed. Output written to {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    gff_cli()