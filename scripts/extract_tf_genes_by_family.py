#!/usr/bin/env python3
"""
Extracts TF genes from Excel file by family.

This script reads an Excel file containing transcription factor gene information
and outputs CSV files for each TF family.
"""

from pangenometools.cli.extract_tf_genes_cli import extract_tf_genes_cli


def main():
    """Main entry point that calls the TF gene extraction CLI."""
    extract_tf_genes_cli()


if __name__ == "__main__":
    main()
