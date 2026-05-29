#!/usr/bin/env python3

"""
Extracts sequences from FASTA files using the new refactored CLI.
"""

from pangenometools.cli.fasta_cli import fasta_cli

def main():
    """Main entry point that calls the new FASTA CLI."""
    fasta_cli()

if __name__ == "__main__":
    main()

