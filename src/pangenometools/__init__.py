"""
PanGenomeTools - A package for interacting with pan-genome data.

This package provides tools for extracting and analyzing genomic data
from pan-genome datasets.
"""

from .models.base import PangenomeFileHandler
from .models.gff import GFFHandler
from .models.fasta import FastaHandler
from .cli.gff_cli import gff_cli
from .cli.fasta_cli import fasta_cli

__version__ = "0.0.1"
__all__ = [
    "PangenomeFileHandler",
    "GFFHandler", 
    "FastaHandler",
    "gff_cli",
    "fasta_cli"
]