"""
FASTA file handler for PanGenomeTools.

This module provides functionality for extracting sequences from FASTA files
using the GFF handler for coordinate information.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union, Literal, Tuple
import pyfaidx
import numpy as np
from Bio.Seq import Seq
from .base import PangenomeFileHandler
from .gff import GFFHandler
import pyranges as pr
from ..utils import calculate_coordinate_boundaries, clip_coordinates


class FastaHandler(PangenomeFileHandler):
    """Handler for FASTA files in pangenome."""

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize FASTA handler.

        Args:
            pangenome_folder: Path to pangenome folder
            pangenome_index: Path to pangenome index file
        """
        super().__init__(pangenome_folder, pangenome_index)
        self.gff_handler = GFFHandler(pangenome_folder, pangenome_index)

    def extract_sequence(
        self,
        genotype: str,
        gene_id: str,
        feature_type: str = "gene",
        upstream: int = 0,
        downstream: int = 0,
        inner_start: int = 0,
        inner_end: int = 0,
        merge_strategy: Literal["merge", "first", "all"] = "merge",
        pad: int = 0,
        use_five_prime_direction: bool = False,
        return_coordinates: bool = False
        ) -> Union[str, Tuple[str, str, int, int, str]]:
        """
        Extract sequence for a gene with given parameters.

        Args:
            genotype: Genotype identifier
            gene_id: Gene ID to extract sequence for
            feature_type: Type of feature to use for coordinates
            upstream: Nucleotides to include upstream (5' direction)
            downstream: Nucleotides to include downstream (3' direction)
            inner_start: Nucleotides to include from start of feature
            inner_end: Nucleotides to include from end of feature
            merge_strategy: How to handle multiple features
            pad: Number of Ns to pad between segments
            use_five_prime_direction: If True, always interpret upstream/downstream
                                     in 5' direction regardless of strand.
                                     If False, reverse upstream/downstream for negative strand.
            return_coordinates: If True, return tuple with coordinates. If False, return sequence only.

        Returns:
            If return_coordinates is False: Extracted sequence as string
            If return_coordinates is True: Tuple of (sequence, chromosome, start, end, strand)

        Raises:
            ValueError: If gene or features not found
            FileNotFoundError: If FASTA file not found
        """
        # Get the feature coordinates
        features = self.gff_handler.get_features_for_gene(
            genotype, gene_id, feature_type, merge_strategy
        )

        if len(features) == 0:
            raise ValueError(f"No features found for gene {gene_id} in {genotype}")

        # For merged features, we'll have a single feature spanning the entire range
        feature = features if isinstance(features, pr.PyRanges) else features[0]
        chrom = feature.df.loc[:,"Chromosome"].iloc[0]
        start = feature.df.loc[:,"Start"].iloc[0]
        end = feature.df.loc[:,"End"].iloc[0]
        strand = feature.df.loc[:,"Strand"].iloc[0]

        # Load FASTA file
        fasta_path = self.resolve_path(genotype, "assembly")
        fasta = pyfaidx.Fasta(str(fasta_path))

        # Adjust for pyranges' 0-based indexing (start included, end excluded)
        # Convert to 1-based inclusive coordinates for pyfaidx
        start_1based = start + 1  # pyranges start is 0-based inclusive
        end_1based = end          # pyranges end is 0-based exclusive, so it's already correct

        # Use shared coordinate calculation function
        left_a, left_b, right_a, right_b = calculate_coordinate_boundaries(
            start_1based, end_1based, strand, upstream, downstream,
            inner_start, inner_end, False, use_five_prime_direction
        )

        # Clip coordinates to chromosome length
        chrom_len = len(fasta[chrom])
        ll, lh, rl, rh = clip_coordinates(left_a, left_b, right_a, right_b, chrom_len)

        # Extract sequences
        left_seq = str(fasta[chrom][ll-1:lh]) if ll <= lh else ""
        right_seq = str(fasta[chrom][rl-1:rh]) if rl <= rh else ""

        # Apply special rules for zero-length segments
        if upstream == 0 and inner_start == 0:
            left_seq = ""
        if downstream == 0 and inner_end == 0:
            right_seq = ""

        # Compose final sequence
        combined = left_seq + (pad * "N") + right_seq

        # Apply strand correction
        if strand == "-":
            combined = str(Seq(combined).reverse_complement())

        if return_coordinates:
            return combined, chrom, start, end, strand
        else:
            return combined

    def extract_sequence_with_coordinates(
        self,
        genotype: str,
        gene_id: str,
        feature_type: str = "gene",
        upstream: int = 0,
        downstream: int = 0,
        inner_start: int = 0,
        inner_end: int = 0,
        merge_strategy: Literal["merge", "first", "all"] = "merge",
        pad: int = 0,
        use_five_prime_direction: bool = False
        ) -> tuple:
        """
        Extract sequence and return with coordinates.

        Args:
            genotype: Genotype identifier
            gene_id: Gene ID to extract sequence for
            feature_type: Type of feature to use for coordinates
            upstream: Nucleotides to include upstream
            downstream: Nucleotides to include downstream
            inner_start: Nucleotides to include from start of feature
            inner_end: Nucleotides to include from end of feature
            merge_strategy: How to handle multiple features
            pad: Number of Ns to pad between segments
            use_five_prime_direction: If True, always interpret upstream/downstream
                                     in 5' direction regardless of strand.
                                     If False, reverse upstream/downstream for negative strand.

        Returns:
            Tuple of (sequence, chromosome, start, end, strand)
        """
        # This method is kept for backward compatibility but now calls the main method
        return self.extract_sequence(
            genotype, gene_id, feature_type, upstream, downstream,
            inner_start, inner_end, merge_strategy, pad, use_five_prime_direction,
            return_coordinates=True
        )
