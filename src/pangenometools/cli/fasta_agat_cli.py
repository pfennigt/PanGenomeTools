"""
FASTA CLI module for PanGenomeTools.

This module provides command-line interface functionality for FASTA sequence extraction.
"""

import argparse
import sys
from pathlib import Path

from pangenometools.utils import replace_genotypes

from ..models.fasta import FastaHandler


def setup_parser() -> argparse.ArgumentParser:
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
    parser.add_argument("--output-dir", default=".",
                       help="Output file path")

    # Sequence extraction options
    parser.add_argument("--feature-type", default="gene",
                       help="Type of feature to use for coordinates (default: gene)")
    parser.add_argument("--genotypes", default=None, nargs="+",
                       help="Genotype(s) to run the analysis for (all by default) ")

    # Extra options to AGAT
    parser.add_argument(
        "--extra",
        default="",
        help="Additional arguments to pass to AGAT (e.g. '--upstream 1000 --downstream 1000')"
    )

    # Additional options
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser

def extract_agat(args: argparse.Namespace) -> None:
    """
    Extract FASTA sequences for target genes.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_folder = Path(args.pangenome_folder)
    pangenome_index = Path(args.pangenome_index)

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Use all genotypes in the pangenome index
    genotypes = fasta_handler.pangenome_index.keys()

    # Replace the genotypes if provided
    genotypes = replace_genotypes(genotypes, args.genotypes)

    # Extract using AGAT
    for genotype in genotypes:
        if not args.silent:
            print(f"Extracting for genotype {genotype}")
        
        fasta_handler.extract_agat(
            genotype=genotype,
            type=args.feature_type,
            extra_args=args.extra,
            output_dir=args.output_dir
        )

        if not args.silent:
            print(f"FASTA extraction completed. Output written to {args.output_dir}/{genotype}.fa")


def main():
    """
    Main entry point for FASTA CLI.

    Parses arguments and executes FASTA sequence extraction.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        extract_agat(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()