"""
TF gene extraction CLI module for PanGenomeTools.

This module provides command-line interface functionality for extracting
TF genes from Excel files by family.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional
from ..models.tf_genes import TFGeneExtractor


def setup_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for TF gene extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract TF genes from Excel file by family"
    )

    # Required arguments
    parser.add_argument("--excel-file", required=True,
                       help="Path to Excel file containing TF gene list")
    parser.add_argument("--output-dir", required=True,
                       help="Directory to output CSV files (one per TF family)")

    # Column names
    parser.add_argument("--tf-family-column", default="TF family",
                       help="Name of column containing TF family names (default: 'TF family')")
    parser.add_argument("--gene-name-column", default="TF Gene",
                       help="Name of column containing gene names (default: 'TF Gene')")

    # Optional arguments
    parser.add_argument("--sheet-name", default=None,
                       help="Name of Excel sheet to read (default: first sheet)")
    parser.add_argument("--families", nargs='+', default=None,
                       help="Specific TF families to extract (default: all families)")

    # Logging
    parser.add_argument("--verbose", action="store_true",
                       help="Enable verbose output")
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser


def extract_tf_genes(args: argparse.Namespace) -> None:
    """
    Extract TF genes from Excel file by family.

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
        # Initialize extractor
        extractor = TFGeneExtractor(
            excel_path=Path(args.excel_file),
            tf_family_column=args.tf_family_column,
            gene_name_column=args.gene_name_column,
            sheet_name=args.sheet_name
        )

        # Extract and save genes by family
        output_files = extractor.extract_and_save_all(
            output_dir=Path(args.output_dir),
            families=args.families
        )

        # Print summary
        if not args.silent:
            print(f"\nSuccessfully extracted TF genes to {args.output_dir}")
            print(f"Created {len(output_files)} family-specific CSV files:")
            for family, filepath in output_files.items():
                print(f"  - {family}: {filepath}")

    except Exception as e:
        logging.error(f"Error: {e}")
        sys.exit(1)


def extract_tf_genes_cli():
    """
    Main entry point for TF gene extraction CLI.

    Parses arguments and executes TF gene extraction.
    """
    parser = setup_parser()
    args = parser.parse_args()

    try:
        extract_tf_genes(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    extract_tf_genes_cli()
