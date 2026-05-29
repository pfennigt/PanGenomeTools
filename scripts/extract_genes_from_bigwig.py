#!/usr/bin/env python3

"""
extract_bigwig_signals.py

Reads BigWig files and extracts float signal values for target genes.
Matches the coordinate logic of the FASTA extractor:
    (upstream → inner-start) + PAD + (inner-end → downstream)

Features:
    - Reads multiple BigWig tracks per genotype (as listed in bigwig_index.csv)
    - Adds condition label to each output entry
    - Uses pyBigWig for low-RAM, streaming extraction
    - Uses pangenome GFF + target_genes.csv + pangenome_index.csv
    - Coordinates clipped to chrom sizes
"""

import sys
from pathlib import Path

# Add the src directory to the path so we can import the module
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pangenometools.cli.bigwig_cli import bigwig_cli
if __name__ == "__main__":
    bigwig_cli()

