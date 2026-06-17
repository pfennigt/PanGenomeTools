"""
FASTA CLI module for PanGenomeTools.

This module provides command-line interface functionality for FASTA sequence extraction.
"""

import argparse
import json
from pathlib import Path
from ..models.base import PangenomeFileHandler

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
    parser.add_argument("--pangenome-index", required=True,
                       help="Path to pangenome index file")
    
    parser.add_argument("--output-file", default="genotypes.json",
                       help="Path to output JSON file")

    parser.add_argument("--genotypes", default=None, nargs="+",
                       help="# Genotype(s) to run the analysis for (all by default)")

    return parser

def get_genotypes_from_pangenome_index(args: argparse.Namespace) -> dict:
    """
    Read a pangenome index and return the genotypes as JSON.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_index = Path(args.pangenome_index)

    return PangenomeFileHandler(Path("."), pangenome_index).pangenome_index


def genotypes_cli():
    """
    Main entry point for FASTA CLI.

    Parses arguments and executes FASTA sequence extraction.
    """
    parser = setup_fasta_parser()
    args = parser.parse_args()

    # Use the genotypes from the input, otherwise read them from the pangenome index
    if args.genotypes is not None:
        genotypes = args.genotypes
    else:
        # Read the genotypes from the pangenome index
        genotypes = list(get_genotypes_from_pangenome_index(args).keys())

    with open(args.output_file, "w") as f:
        json.dump(genotypes, f)

if __name__ == "__main__":
    genotypes_cli()