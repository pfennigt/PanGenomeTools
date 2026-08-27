"""
Homolog finding CLI module for PanGenomeTools.

This module provides command-line interface functionality for finding
homologs using BLAST or DIAMOND.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional
from ..models.homolog_finder import HomologFinder


def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for homolog finding.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Find homologs of query sequences in target genome using BLAST or DIAMOND"
    )

    # Required arguments
    parser.add_argument("--query-fasta", required=True,
                       help="Path to query FASTA file (TF protein sequences)")
    parser.add_argument("--target-genome", required=True,
                       help="Path to target genome FASTA file")
    parser.add_argument("--output", default=None,
                       help="Output file path for results")
    parser.add_argument("--db-path", default=None,
                       help="Output file path for results")

    # GFF for Gene translation
    parser.add_argument("--target-gff", default=None, 
                    help="Path to target genome GFF file")

    # Method selection
    parser.add_argument("--method", choices=["blast", "diamond"], default="blast",
                       help="Search method - 'blast' or 'diamond' (default: blast)")

    # Search parameters
    parser.add_argument("--evalue", type=float, default=1e-5,
                       help="E-value threshold (default: 1e-5)")
    parser.add_argument("--max-target-seqs", type=int, default=10,
                       help="Maximum number of target sequences to report (default: 10)")
    parser.add_argument("--output-format", default="6",
                       help="Output format (default: '6' for tabular)")
    parser.add_argument("--num-threads", type=int, default=1,
                       help="Number of threads to use (default: 1)")

    # BLAST-specific parameters
    parser.add_argument("--blast-type", choices=["blastp", "blastx", "tblastn", "tblastx"],
                       default="blastp",
                       help="BLAST type (default: blastp)")

    # DIAMOND-specific parameters
    parser.add_argument("--diamond-sensitivity",
                       choices=["fast", "mid-sensitive", "sensitive", "more-sensitive",
                                "very-sensitive", "ultra-sensitive"],
                       default="sensitive",
                       help="DIAMOND sensitivity level (default: sensitive)")

    # Optional arguments
    parser.add_argument("--keep-database", action="store_true",
                       help="Keep the database file after search")

    # Logging
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def find_homologs(args: argparse.Namespace) -> None:
    """
    Find homologs of query sequences in target genome.

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
        # Initialize finder
        if args.db_path is not None and Path(args.db_path).is_dir():
            db_path = args.db_path + "/database"
        else:
            db_path = args.db_path

        finder = HomologFinder(method=args.method, db_path=db_path)

        # Check if required tool is available
        if args.method == "blast":
            if not finder.check_blast_available():
                raise RuntimeError("BLAST+ is not installed or not in PATH. "
                                 "Please install BLAST+: conda install -c bioconda blast")
        else:  # diamond
            if not finder.check_diamond_available():
                raise RuntimeError("DIAMOND is not installed or not in PATH. "
                                 "Please install DIAMOND: conda install -c bioconda diamond")

        # Use the quesry name for the output by default
        if args.output is None:
            output = Path(args.query_fasta).stem + ".tsv"
        else:
            output = args.output          
        

        # Find homologs
        stats = finder.find_homologs(
            query_fasta=Path(args.query_fasta),
            target_fasta=Path(args.target_genome),
            target_gff=Path(args.target_gff) if args.target_gff is not None else None,
            output_file=Path(output),
            evalue=args.evalue,
            max_target_seqs=args.max_target_seqs,
            output_format=args.output_format,
            blast_type=args.blast_type,
            diamond_sensitivity=args.diamond_sensitivity,
            num_threads=args.num_threads,
            keep_database=args.keep_database
        )

        # Print summary
        if not args.silent:
            print(f"\nSuccessfully completed homolog search using {args.method.upper()}")
            print(f"Results saved to: {args.output}")
            print(f"\nSearch Statistics:")
            print(f"  Query sequences: {stats['num_queries']}")
            print(f"  Total hits: {stats['num_hits']}")
            print(f"  Queries with hits: {stats['num_queries_with_hits']}")

            if stats['num_queries_with_hits'] == 0:
                print("\nWarning: No homologs found. Consider:")
                print("  1. Relaxing e-value threshold (--evalue)")
                print("  2. Increasing max target sequences (--max-target-seqs)")
                print("  3. Checking that target genome contains protein sequences")

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def find_homologs_cli():
    """
    Main entry point for homolog finding CLI.

    Parses arguments and executes homolog search.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        find_homologs(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    find_homologs_cli()
