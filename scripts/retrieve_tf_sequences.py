#!/usr/bin/env python3
"""
Retrieves protein sequences for TF genes from reference genome.

This script uses AGAT to extract and translate protein sequences from
a reference genome based on a list of gene names.
"""

from pangenometools.cli.retrieve_sequences_cli import retrieve_sequences_cli


def main():
    """Main entry point that calls the sequence retrieval CLI."""
    retrieve_sequences_cli()


if __name__ == "__main__":
    main()
