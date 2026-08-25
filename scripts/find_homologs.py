#!/usr/bin/env python3
"""
Find Homologs of genes in a target genome using BLAST or DIAMOND.

This script uses BLAST or DIAMOND to find significant matches of from a protein fasta
in the translated CDSs from a reference genome based. If a genome is given, it is 
translated to protein using AGAT.
"""

from pangenometools.cli.find_homologs_cli import find_homologs_cli


def main():
    """Main entry point that calls the homolog finding CLI."""
    find_homologs_cli()


if __name__ == "__main__":
    main()
