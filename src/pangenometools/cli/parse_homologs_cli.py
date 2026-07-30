"""
Homolog parsing CLI module for PanGenomeTools.

This module provides command-line interface functionality for parsing
and summarizing homolog search results.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional
from ..models.homolog_parser import HomologParser


def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for homolog parsing.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Parse and summarize homolog search results"
    )

    # Required arguments
    parser.add_argument("--blast-results", required=True,
                       help="Path to BLAST/DIAMOND results file")
    parser.add_argument("--query-fasta", required=True,
                       help="Path to original query FASTA file")
    parser.add_argument("--output-dir", required=True,
                       help="Directory for output files")

    # Filtering parameters
    parser.add_argument("--min-identity", type=float, default=30.0,
                       help="Minimum percent identity threshold (default: 30.0)")
    parser.add_argument("--min-coverage", type=float, default=50.0,
                       help="Minimum query coverage threshold (default: 50.0)")
    parser.add_argument("--top-hits-only", action="store_true", default=True,
                       help="Only keep top hit per query (default: True)")
    parser.add_argument("--keep-all-hits", action="store_true",
                       help="Keep all hits (overrides --top-hits-only)")

    # Logging
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def parse_homologs(args: argparse.Namespace) -> None:
    """
    Parse and summarize homolog search results.

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
        # Initialize parser
        parser = HomologParser(
            blast_results=Path(args.blast_results),
            query_fasta=Path(args.query_fasta)
        )

        # Determine if we should keep all hits
        top_hits_only = not args.keep_all_hits

        # Parse and save results
        summary_file, detailed_file, stats_file = parser.save_results(
            output_dir=Path(args.output_dir),
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            top_hits_only=top_hits_only
        )

        # Print summary
        if not args.silent:
            print(f"\nSuccessfully parsed homolog results")
            print(f"Output directory: {args.output_dir}")
            print(f"\nGenerated files:")
            print(f"  - Summary: {summary_file}")
            print(f"  - Detailed: {detailed_file}")
            print(f"  - Statistics: {stats_file}")

            # Print statistics
            with open(stats_file, 'r') as f:
                print("\n" + f.read())

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def parse_homologs_cli():
    """
    Main entry point for homolog parsing CLI.

    Parses arguments and executes homolog parsing.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        parse_homologs(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parse_homologs_cli()
