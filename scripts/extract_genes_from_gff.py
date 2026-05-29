#!/usr/bin/env python3

"""
Extracts portions of GFF files related to a target gene using the new refactored CLI.
"""

from pangenometools.cli.gff_cli import gff_cli

def main():
    """Main entry point that calls the new GFF CLI."""
    gff_cli()

if __name__ == "__main__":
    main()

