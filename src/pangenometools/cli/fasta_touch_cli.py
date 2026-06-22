"""
FASTA CLI module for PanGenomeTools.

This module provides command-line interface functionality for FASTA sequence extraction.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List
from ..models.fasta import FastaHandler
from ..models.gff import GFFHandler
from pangenometools.utils import replace_genotypes

def setup_fasta_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for FASTA sequence extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences from FASTA files using GFF coordinates."
    )

    # Required arguments
    parser.add_argument("--pangenome-folder", required=True,
                       help="Path to pangenome folder")
    parser.add_argument("--pangenome-index", required=True,
                       help="Path to pangenome index file")
    parser.add_argument("--target-genes", default=None,
                       help="Path to target genes CSV file")
    parser.add_argument("--genotypes", default=None, nargs="+",
                       help="Genotype(s) to run the analysis for (all by default) ")

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

# FASTA extraction using a target genes file
def touch_fasta(fasta_handler:FastaHandler, genotypes, args):
    # Process each genotype
    for g, genotype in enumerate(genotypes):

        try:
            print(f"Testing loading of {genotype}", file=sys.stdout)
            # Extract sequence
            fasta_handler.touch_fasta(
                genotype,
            )
        except Exception as e:
            if not args.silent:
                print(f"Error processing {genotype}: {e}", file=sys.stderr)
            continue

def create_fasta_fai(args: argparse.Namespace) -> None:
    """
    Extract FASTA sequences for target genes.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_folder = Path(args.pangenome_folder)
    pangenome_index = Path(args.pangenome_index)

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Read target genes or use all genotypes in the pangenome index
    if args.target_genes is not None:
        genotypes, target_rows = read_target_genes(Path(args.target_genes))

        # Replace the genotypes if provided
        genotypes = replace_genotypes(genotypes, args.genotypes)

    else:
        # Use all genotypes in the pangenome index
        genotypes = fasta_handler.pangenome_index.keys()

        # Replace the genotypes if provided
        genotypes = replace_genotypes(genotypes, args.genotypes)

    # Touch the fasta files
    touch_fasta(fasta_handler, genotypes,args)

def fasta_cli():
    """
    Main entry point for FASTA CLI.

    Parses arguments and executes FASTA sequence extraction.
    """
    parser = setup_fasta_parser()
    args = parser.parse_args()

    try:
        create_fasta_fai(args)
        if not args.silent:
            print("FASTA fai creation completed", file=sys.stdout)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    fasta_cli()