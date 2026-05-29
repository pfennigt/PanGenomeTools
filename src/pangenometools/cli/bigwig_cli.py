"""
BigWig CLI module for PanGenomeTools.

This module provides command-line interface functionality for BigWig signal extraction.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List
from ..models.bigwig import BigWigHandler


def setup_bigwig_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for BigWig signal extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract signal values from BigWig files using GFF coordinates."
    )

    # Required arguments
    parser.add_argument("--pangenome-folder", required=True,
                       help="Path to pangenome folder")
    parser.add_argument("--pangenome-index", required=True,
                       help="Path to pangenome index file")
    parser.add_argument("--bigwig-folder", required=True,
                       help="Path to BigWig folder")
    parser.add_argument("--bigwig-index", required=True,
                       help="Path to BigWig index CSV file")
    parser.add_argument("--target-genes", required=True,
                       help="Path to target genes CSV file")
    parser.add_argument("--output", required=True,
                       help="Output JSON file path")

    # Signal extraction options
    parser.add_argument("--feature-type", default="gene",
                       help="Type of feature to use for coordinates (default: gene)")
    parser.add_argument("--upstream", type=int, default=0,
                       help="Nucleotides to include upstream (5' direction)")
    parser.add_argument("--downstream", type=int, default=0,
                       help="Nucleotides to include downstream (3' direction)")
    parser.add_argument("--inner-start", type=int, default=0,
                       help="Nucleotides to include from start of feature")
    parser.add_argument("--inner-end", type=int, default=0,
                       help="Nucleotides to include from end of feature")
    parser.add_argument("--pad", type=int, default=0,
                       help="Number of null values to pad between segments")

    # Feature handling options
    parser.add_argument("--whole-seq", action="store_true",
                       help="Extract entire feature sequence")
    parser.add_argument("--use-five-prime-direction", action="store_true",
                       help="Always interpret upstream/downstream in 5' direction regardless of strand")

    # Additional options
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def extract_bigwig_signals(args: argparse.Namespace) -> None:
    """
    Extract BigWig signals for target genes.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_folder = Path(args.pangenome_folder)
    pangenome_index = Path(args.pangenome_index)

    bigwig_handler = BigWigHandler(pangenome_folder, pangenome_index)

    # Extract signals
    bigwig_handler.extract_and_save_bigwig_signals(
        bigwig_folder=Path(args.bigwig_folder),
        bigwig_index_path=Path(args.bigwig_index),
        target_genes_path=Path(args.target_genes),
        output_path=Path(args.output),
        feature_type=args.feature_type,
        upstream=args.upstream,
        downstream=args.downstream,
        whole_seq=args.whole_seq,
        inner_start=args.inner_start,
        inner_end=args.inner_end,
        pad=args.pad,
        use_five_prime_direction=args.use_five_prime_direction
    )


def bigwig_cli():
    """
    Main entry point for BigWig CLI.

    Parses arguments and executes BigWig signal extraction.
    """
    parser = setup_bigwig_parser()
    args = parser.parse_args()

    try:
        extract_bigwig_signals(args)
        if not args.silent:
            print(f"BigWig extraction completed. Output written to {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    bigwig_cli()