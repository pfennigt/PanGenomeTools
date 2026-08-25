"""
Sequence retrieval CLI module for PanGenomeTools.

This module provides command-line interface functionality for retrieving
protein sequences from reference genomes using AGAT.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional
from ..models.sequence_retriever import SequenceRetriever


def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for sequence retrieval.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Retrieve protein sequences for TF genes from reference genome"
    )

    # Required arguments
    parser.add_argument("--gene-list", required=True,
                       help="Path to CSV file with gene names")
    parser.add_argument("--reference-genome", required=True,
                       help="Path to reference genome FASTA file")
    parser.add_argument("--reference-gff", required=True,
                       help="Path to reference GFF annotation file")
    parser.add_argument("--output", default=None,
                       help="Output FASTA file path")

    # Optional arguments
    parser.add_argument("--feature-type", default="CDS",
                       help="Type of feature to extract (default: CDS)")
    parser.add_argument("--gene-id-column", default="gene_name",
                       help="Column name in CSV containing gene IDs (default: gene_name)")
    parser.add_argument("--keep-temp-files", action="store_true",
                       help="Keep temporary GFF files")

    # Logging
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def retrieve_sequences(args: argparse.Namespace) -> None:
    """
    Retrieve protein sequences for TF genes.

    Args:
        args: Parsed command line arguments
    """
    # Setup logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    elif not args.silent:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    try:
        # Initialize retriever
        retriever = SequenceRetriever(
            reference_genome=Path(args.reference_genome),
            reference_gff=Path(args.reference_gff)
        )

        # Check if AGAT is available
        if not retriever.check_agat_available():
            raise RuntimeError("AGAT is not installed or not in PATH. "
                             "Please install AGAT: conda install -c bioconda agat")

        # If no output name is given, use the input name
        if args.output is None:
            output = Path(args.gene_list).stem + "fa"
        else:
            output = args.output

        # Retrieve sequences
        num_genes, num_sequences = retriever.retrieve_from_gene_list_csv(
            gene_list_csv=Path(args.gene_list),
            output_fasta=Path(output),
            feature_type=args.feature_type,
            gene_id_column=args.gene_id_column,
            keep_temp_files=args.keep_temp_files
        )

        # Print summary
        if not args.silent:
            print(f"\nSuccessfully retrieved sequences to {args.output}")
            print(f"Found {num_genes} genes")
            print(f"Retrieved {num_sequences} protein sequences")

            if num_sequences == 0:
                print("\nWarning: No sequences were retrieved. Please check:")
                print("  1. Gene names in CSV match those in GFF file")
                print("  2. GFF file contains CDS features")
                print("  3. Reference genome and GFF are compatible")

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def retrieve_sequences_cli():
    """
    Main entry point for sequence retrieval CLI.

    Parses arguments and executes sequence retrieval.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        retrieve_sequences(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    retrieve_sequences_cli()
