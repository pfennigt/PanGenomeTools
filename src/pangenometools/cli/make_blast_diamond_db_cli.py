"""
BLAST/DIAMOND database creation CLI module for PanGenomeTools.

This module provides command-line interface functionality for creating BLAST or DIAMOND databases for streamlined homolog finding.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional
from ..models.homolog_finder import HomologFinder


def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for db creation.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Create a database of a protein translated genome using BLAST or DIAMOND"
    )

    # Required arguments
    parser.add_argument("--target-genome", required=True,
                       help="Path to target genome FASTA file")

    # GFF for Gene translation
    parser.add_argument("--target-gff", required=True, 
                    help="Path to target genome GFF file")

    # Database path
    parser.add_argument("--db-path", default=None, 
                    help="Path to output database")

    # Method selection
    parser.add_argument("--method", choices=["blast", "diamond"], default="blast",
                       help="Search method - 'blast' or 'diamond' (default: blast)")

    # Logging
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def make_db(args: argparse.Namespace) -> None:
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
        finder = HomologFinder(method=args.method)

        # Check if required tool is available
        if args.method == "blast":
            if not finder.check_blast_available():
                raise RuntimeError("BLAST+ is not installed or not in PATH. "
                                 "Please install BLAST+: conda install -c bioconda blast")
        else:  # diamond
            if not finder.check_diamond_available():
                raise RuntimeError("DIAMOND is not installed or not in PATH. "
                                 "Please install DIAMOND: conda install -c bioconda diamond")

        # Set a default name for the database
        if args.db_path is None:
            # Create the db folder
            Path("db").mkdir(exist_ok=True)
            db_name="db/database"
        else:
            db_name=args.db_path

        # Create the database
        db_path = finder.create_database(
            target_fasta=Path(args.target_genome),
            db_name=db_name,
            target_gff=Path(args.target_gff))

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def make_blast_diamond_db():
    """
    Main entry point for DB creation CLI.

    Parses arguments and executes DB creation.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        make_db(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    make_blast_diamond_db()
