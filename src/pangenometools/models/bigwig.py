"""
BigWig handler module for PanGenomeTools.

This module provides functionality for extracting signal values from BigWig files
using the same coordinate logic as the FASTA extractor.
"""

import csv
import gzip
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Literal
import pyBigWig
import pyranges as pr
from .base import PangenomeFileHandler
from .gff import GFFHandler
from ..utils import calculate_coordinate_boundaries, clip_coordinates, read_target_genes


class BigWigHandler(PangenomeFileHandler):
    """Handler for extracting signal values from BigWig files."""

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize BigWig handler.

        Args:
            pangenome_folder: Path to pangenome folder
            pangenome_index: Path to pangenome index file
        """
        super().__init__(pangenome_folder, pangenome_index)
        self.gff_handler = GFFHandler(pangenome_folder, pangenome_index)

    def read_bigwig_index(self, path: Path) -> List[Dict[str, str]]:
        """Read BigWig index file."""
        with open(path, newline="") as fh:
            return list(csv.DictReader(fh))

    def extract_bigwig_signals_for_gene(
        self,
        genotype: str,
        gene_id: str,
        bigwig_path: Path,
        condition: str,
        feature_type: str = "gene",
        upstream: int = 0,
        downstream: int = 0,
        whole_seq: bool = False,
        inner_start: int = 0,
        inner_end: int = 0,
        pad: int = 0,
        use_five_prime_direction: bool = False
    ) -> Optional[Dict[str, Any]]:
        """Extract BigWig signal values for a single gene."""
        # Get the feature coordinates
        features = self.gff_handler.get_features_for_gene(
            genotype, gene_id, feature_type, "merge"
        )

        if len(features) == 0:
            self.logger.warning(f"No features found for gene {gene_id} in {genotype}")
            return None

        # For merged features, we'll have a single feature spanning the entire range
        feature = features if isinstance(features, pr.PyRanges) else features[0]
        chrom = feature.df.loc[:,"Chromosome"].iloc[0]
        start = feature.df.loc[:,"Start"].iloc[0]
        end = feature.df.loc[:,"End"].iloc[0]
        strand = feature.df.loc[:,"Strand"].iloc[0]

        # Adjust for pyranges' 0-based indexing (start included, end excluded)
        # Convert to 1-based inclusive coordinates for BigWig
        start_1based = start + 1  # pyranges start is 0-based inclusive
        end_1based = end          # pyranges end is 0-based exclusive, so it's already correct

        # Use shared coordinate calculation function
        left_a, left_b, right_a, right_b, additional_padding = calculate_coordinate_boundaries(
            start_1based, end_1based, strand, upstream, downstream,
            inner_start, inner_end, whole_seq, use_five_prime_direction
        )

        # Open BigWig file
        try:
            bw = pyBigWig.open(str(bigwig_path))
        except Exception as e:
            self.logger.error(f"Could not open BigWig {bigwig_path}: {e}")
            return None

        if chrom not in bw.chroms():
            self.logger.warning(f"Chromosome {chrom} missing in {bigwig_path}")
            bw.close()
            return None

        chrom_len = bw.chroms()[chrom]

        # Use shared coordinate clipping function
        ll, lh, rl, rh = clip_coordinates(left_a, left_b, right_a, right_b, chrom_len)

        # Extract signal
        left_vals = bw.values(chrom, ll - 1, lh) if ll <= lh else []
        right_vals = bw.values(chrom, rl - 1, rh) if rl <= rh else []

        # Convert nan to None
        left_vals = [None if (v is None or (isinstance(v, float) and v != v)) else v for v in left_vals]
        right_vals = [None if (v is None or (isinstance(v, float) and v != v)) else v for v in right_vals]

        pad_list = [None] * (pad + additional_padding)

        combined = left_vals + pad_list + right_vals

        # Orientation
        if strand == "-":
            combined = combined[::-1]

        # Return result
        return {
            "gene_id": gene_id,
            "condition": condition,
            "location": f"{chrom}:{start}-{end}({strand})",
            "data": combined,
        }

    def extract_bigwig_signals(
        self,
        bigwig_folder: Path,
        bigwig_index_path: Path,
        target_genes_path: Path,
        feature_type: str = "gene",
        upstream: int = 0,
        downstream: int = 0,
        whole_seq: bool = False,
        inner_start: int = 0,
        inner_end: int = 0,
        pad: int = 0,
        use_five_prime_direction: bool = False
    ) -> List[Dict[str, Any]]:
        """Extract BigWig signal values for target genes."""
        # Read index files
        bigwig_index_rows = self.read_bigwig_index(bigwig_index_path)
        genotypes_in_targets, target_rows = read_target_genes(target_genes_path)

        # Build genotype → list of bigwig tracks
        bw_map: Dict[str, List[Dict[str, str]]] = {}
        for r in bigwig_index_rows:
            geno = r["genotype"]
            bw_map.setdefault(geno, []).append({
                "condition": r["condition"],
                "path": bigwig_folder / r["bigwig"],
            })

        # Collect outputs here (one JSON object per row+condition)
        outputs = []

        # Process by genotype (low RAM)
        for geno in genotypes_in_targets:
            if geno not in self.pangenome_index:
                self.logger.warning(f"Genotype {geno} missing in pangenome index")
                continue

            if geno not in bw_map:
                self.logger.warning(f"Genotype {geno} has no BigWigs")
                continue

            # Process each requested gene
            for row in target_rows:
                gene_name = row.get("gene_name", "")
                gene_id = row.get(f"gene_ID_{geno}", "")
                if not gene_id:
                    continue

                # For each BigWig track for this genotype
                for bwinfo in bw_map[geno]:
                    condition = bwinfo["condition"]
                    bwpath = bwinfo["path"]

                    result = self.extract_bigwig_signals_for_gene(
                        genotype=geno,
                        gene_id=gene_id,
                        bigwig_path=bwpath,
                        condition=condition,
                        feature_type=feature_type,
                        upstream=upstream,
                        downstream=downstream,
                        whole_seq=whole_seq,
                        inner_start=inner_start,
                        inner_end=inner_end,
                        pad=pad,
                        use_five_prime_direction=use_five_prime_direction
                    )

                    if result:
                        result["genotype"] = geno
                        result["gene_name"] = gene_name
                        outputs.append(result)

        return outputs

    def extract_and_save_bigwig_signals(
        self,
        bigwig_folder: Path,
        bigwig_index_path: Path,
        target_genes_path: Path,
        output_path: Path,
        feature_type: str = "gene",
        upstream: int = 0,
        downstream: int = 0,
        whole_seq: bool = False,
        inner_start: int = 0,
        inner_end: int = 0,
        pad: int = 0,
        use_five_prime_direction: bool = False
    ) -> None:
        """Extract BigWig signals and save to JSON file."""
        outputs = self.extract_bigwig_signals(
            bigwig_folder=bigwig_folder,
            bigwig_index_path=bigwig_index_path,
            target_genes_path=target_genes_path,
            feature_type=feature_type,
            upstream=upstream,
            downstream=downstream,
            whole_seq=whole_seq,
            inner_start=inner_start,
            inner_end=inner_end,
            pad=pad,
            use_five_prime_direction=use_five_prime_direction
        )

        # Write JSON
        with open(output_path, "w") as fh:
            json.dump(outputs, fh, indent=2)