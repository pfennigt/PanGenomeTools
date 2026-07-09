"""
FASTA file handler for PanGenomeTools.

This module provides functionality for extracting sequences from FASTA files
using the GFF handler for coordinate information.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union, Literal, Tuple
import pyfaidx
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from .base import PangenomeFileHandler
from .gff import GFFHandler
import pyranges as pr
from ..utils import calculate_coordinate_boundaries, clip_coordinates
import sys


class FastaHandler(PangenomeFileHandler):
    """Handler for FASTA files in pangenome."""

    _fasta: Union[None,pyfaidx.Fasta] = None
    _genotype: Union[None,str] = None

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize FASTA handler.

        Args:
            pangenome_folder: Path to pangenome folder
            pangenome_index: Path to pangenome index file
        """
        super().__init__(pangenome_folder, pangenome_index)
        self.gff_handler = GFFHandler(pangenome_folder, pangenome_index)

    def load_fasta(self, genotype, rebuild=False) -> pyfaidx.Fasta:
        # See if the required fasta is still open
        if self._fasta is None or self._genotype != genotype:

            # Close the open fasta
            if self._genotype != genotype:
                self.close_fasta()

            # Open the fasta if not
            fasta_path = self.resolve_path(genotype, "assembly")

            try:
                self._fasta = pyfaidx.Fasta(str(fasta_path), rebuild=rebuild)
            except Exception:
                self._fasta = pyfaidx.Fasta(str(fasta_path))
            self._genotype = genotype
        
        return self._fasta

    def close_fasta(self):
        # Close the fasta if it is open
        if self._fasta is not None:
            self._fasta.close()

            self._fasta = self._genotype = None

    def touch_fasta(self, genotype, rebuild=True):
        self.load_fasta(genotype, rebuild=rebuild)
        self.close_fasta()

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
        whole_seq = False,
        skip_short_chromosomes=True,
        skip_short_genes=False,
        use_five_prime_direction: bool = False,
        return_info: bool = False,
        _use_cache: bool = False
        ) -> Union[str, Tuple[str, dict]]:
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
            pad: Number of Ns to pad between segments (-1 treated as 0)
            whole_seq: Extract the whole sequence between start and end
            use_five_prime_direction: If True, always interpret upstream/downstream
                                     in 5' direction regardless of strand.
                                     If False, reverse upstream/downstream for negative strand.
            return_info: If True, return dict with extraction info. If False, return sequence only.

        Returns:
            If return_info is False: Extracted sequence as string
            If return_info is True: Tuple of (sequence, info_dict)

        Raises:
            ValueError: If gene or features not found
            FileNotFoundError: If FASTA file not found
        """

        # If pad=-1 flag is given, use pad=0 instead
        if pad == -1:
            pad=0

        # Get the feature coordinates
        features = self.gff_handler.get_features_for_gene(
            genotype, gene_id, feature_type, merge_strategy
        )

        if len(features) == 0:
            raise ValueError(f"No features found for gene {gene_id} in {genotype}")

        # For merged features, we'll have a single feature spanning the entire range
        feature = features if isinstance(features, pd.DataFrame) else features[0]
        chrom = feature.loc[:,"Chromosome"].iloc[0]
        start = feature.loc[:,"Start"].iloc[0]
        end = feature.loc[:,"End"].iloc[0]
        strand = feature.loc[:,"Strand"].iloc[0]

        # Load FASTA file
        fasta = self.load_fasta(genotype)

        # Adjust for pyranges' 0-based indexing (start included, end excluded)
        # Convert to 1-based inclusive coordinates for pyfaidx
        start_1based = start + 1  # pyranges start is 0-based inclusive
        end_1based = end          # pyranges end is 0-based exclusive, so it's already correct

        # Use shared coordinate calculation function
        left_a, left_b, right_a, right_b, additional_padding = calculate_coordinate_boundaries(
            start_1based, end_1based, strand, upstream, downstream,
            inner_start, inner_end, whole_seq, use_five_prime_direction
        )

        # Clip coordinates to chromosome length
        chrom_len = len(fasta[chrom])
        ll, lh, rl, rh = clip_coordinates(left_a, left_b, right_a, right_b, chrom_len)

        # Check if the outer bounds of the flanking regions are outside of the chromosome
        if ll!=left_a or rh!=right_b:
            # If set, skip the gene
            if skip_short_chromosomes:
                
                # Close the FASTA
                if not _use_cache:
                    self.close_fasta()

                if return_info:
                    return "", {}
                else:
                    return ""

            else:
                left_pad = "N" * np.abs(left_a - ll)
                right_pad = "N" * np.abs(right_b - rh)
        else:
            left_pad = ""
            right_pad = ""
        # Extract sequences
        left_seq = (left_pad + str(fasta[chrom][ll-1:lh])) if ll <= lh else ""
        right_seq = (str(fasta[chrom][rl-1:rh]) + right_pad) if rl <= rh else ""

        # Get the target length os the extracted sequences
        front_len = upstream + inner_start
        back_len = upstream + inner_start
        left_target_len = front_len if strand == "+" else back_len
        right_target_len = back_len if strand == "+" else front_len
        
        # Check if the target lengths were met
        if additional_padding != 0:
            # If set, skip the gene
            if skip_short_genes:
                # Close the FASTA
                if not _use_cache:
                    self.close_fasta()

                if return_info:
                    return "", {}
                else:
                    return ""
            else:
                # Apply padding to reach target length
                left_seq = left_seq + "N" * (left_target_len - len(left_seq))
                right_seq = right_seq + "N" * (right_target_len - len(right_seq))

                # Check that the calculated length now match
                if not len(left_seq) == left_target_len:
                    raise ValueError(f"left seq wrong size: {len(left_seq)} != {left_target_len}")
                if not len(right_seq) == right_target_len:
                    raise ValueError(f"right seq wrong size: {len(right_seq)} != {right_target_len}")
                if not additional_padding == (left_target_len - len(left_seq)) + (right_target_len - len(right_seq)):
                    raise ValueError(f"padding wrong size: {additional_padding} != {(left_target_len - len(left_seq)) + (right_target_len - len(right_seq))}")

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

        # Check that the final length is correct
        if not len(combined) == upstream + inner_start + pad + inner_end + downstream:
            raise ValueError(f"Final length is not correct: {gene_id} ({start}-{end} ({strand})): {len(combined)} != {upstream + inner_start + pad + inner_end + downstream}\n ll:{ll}, lh:{lh}, rl:{rl}, rh:{rh} (sum: {lh - ll + rh - rl})\n left_target_len:{left_target_len}, right_target_len:{right_target_len}")

        # Close the FASTA
        if not _use_cache:
            self.close_fasta()

        if return_info:
            info = {
                "chrom":chrom,
                "strand":strand,
                "ll":ll,
                "lh":lh,
                "rl":rl,
                "rh":rh,
                "left_len": len(left_seq),
                "right_len": len(right_seq),
                "additional_padding": additional_padding,
            }
            return combined, info
        else:
            return combined
