"""
FASTA CLI module for PanGenomeTools.

This module provides command-line interface functionality for FASTA sequence extraction.
"""

import argparse
import json
from pathlib import Path
from ..models.csv import parse_csv_columns

def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for CSV parser extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract columns from a CSV and return them as JSON."
    )

    # Required arguments
    parser.add_argument("--csv", required=True,
                       help="Path to pangenome index file")
    
    parser.add_argument("--output-dir", default=None,
                       help="Output directory for JSON files")

    parser.add_argument("--groupby", default=None,
                       help="Column to group all columns by")

    return parser

def main():
    """
    Main entry point for FASTA CLI.

    Parses arguments and executes FASTA sequence extraction.
    """
    parser = setup_parser()
    args = parser.parse_args()

    parsed_columns = parse_csv_columns(csv_path=args.csv, output_dir=args.output_dir, groupby=args.groupby)

    return parsed_columns

if __name__ == "__main__":
    main()