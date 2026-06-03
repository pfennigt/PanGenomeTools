"""
GFF/GTF file handler for PanGenomeTools.

This module provides functionality for working with GFF/GTF files using pyranges,
including feature finding, merging, and transcript-to-gene mapping.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union, Literal
import pyranges as pr
import pandas as pd
import re
import numpy as np
from .base import PangenomeFileHandler


class GFFHandler(PangenomeFileHandler):
    """Handler for GFF/GTF files in pangenome."""

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize GFF handler.

        Args:
            pangenome_folder: Path to pangenome folder
            pangenome_index: Path to pangenome index file
        """
        super().__init__(pangenome_folder, pangenome_index)

    def load_gff(self, genotype: str, feature_type: Optional[str] = None) -> pr.PyRanges:
        """
        Load GFF file for a genotype, optionally filtered by feature type.

        Args:
            genotype: Genotype identifier
            feature_type: Optional feature type to filter by

        Returns:
            PyRanges object containing GFF features
        """
        gff_path = self.resolve_path(genotype, "annotation")
        features = pr.read_gff3(gff_path)
        return features if feature_type is None else features.features[feature_type]

    def get_features_for_gene(self,
                            genotype: str,
                            gene_id: str,
                            feature_type: str = "CDS",
                            merge_strategy: Literal["merge", "first", "all"] = "merge"
                            ) -> Union[pr.PyRanges, List[pr.PyRanges]]:
        """
        Get features of a specific type for a gene, with optional merging.

        This is particularly useful for features like CDS that may have multiple
        segments per gene.

        Args:
            genotype: Genotype identifier
            gene_id: Gene ID to search for
            feature_type: Type of feature to retrieve (CDS, exon, etc.)
            merge_strategy: How to handle multiple features ('merge', 'first', 'all')

        Returns:
            Merged or individual features

        Raises:
            ValueError: If merge_strategy is invalid
        """
        # Load all features for this genotype
        all_features = self.load_gff(genotype)

        # Find gene features
        gene_features = all_features[all_features.ID == gene_id]

        if len(gene_features) == 0:
            self.logger.warning(f"Gene {gene_id} not found in {genotype}")
            return gene_features

        # Get the IDs of children of this gene
        children_ids = set(all_features.ID[all_features.Parent == gene_id])

        # Find all features of the requested type that are children or grandchildren of this gene
        gene_children = all_features[(all_features.Parent == gene_id) | (all_features.ID == gene_id) | all_features.Parent.isin(children_ids)]

        if feature_type != "all":
            gene_children = gene_children[gene_children.Feature == feature_type]

        if len(gene_children) == 0:
            self.logger.warning(f"No {feature_type} features found for gene {gene_id} in {genotype}")
            return gene_children

        if merge_strategy == "merge":
            return self._merge_features(gene_children)
        elif merge_strategy == "all":
            return self._split_by_strand_and_chrom(gene_children)
        elif merge_strategy == "first":
            return gene_children.head(1)
        else:
            raise ValueError(f"Unknown merge strategy: {merge_strategy}")

    def _merge_features(self, features: pr.PyRanges) -> pr.PyRanges:
        """
        Merge multiple features into a single feature spanning the entire range.

        For features on the same strand and chromosome:
        - Uses the minimum start position
        - Uses the maximum end position
        - Preserves the strand
        - Combines attributes

        Args:
            features: PyRanges object containing features to merge

        Returns:
            PyRanges object with merged features
        """
        if len(features) == 1:
            return features

        # Group by chromosome and strand
        grouped = features.df.groupby(['Chromosome', 'Strand'])

        merged_features = []
        for (chrom, strand), group in grouped:
            # Get min start and max end
            start = group['Start'].min()
            end = group['End'].max()

            # Combine attributes
            attributes = self._combine_attributes(group)

            # Create new feature
            new_feature = pd.DataFrame({
                'Chromosome': [chrom],
                'Start': [start],
                'End': [end],
                'Strand': [strand],
                **attributes
            })

            merged_features.append(pr.PyRanges(new_feature))

        if not merged_features:
            return pr.PyRanges(pd.DataFrame())  # Empty result

        return pr.PyRanges(pd.concat([m.df for m in merged_features], ignore_index=True))

    def _split_by_strand_and_chrom(self, features: pr.PyRanges) -> List[pr.PyRanges]:
        """
        Split features by chromosome and strand, returning a list.

        Args:
            features: PyRanges object to split

        Returns:
            List of PyRanges objects, one per chromosome/strand combination
        """
        grouped = features.df.groupby(['Chromosome', 'Strand'])
        return [pr.PyRanges(group) for _, group in grouped]

    def _combine_attributes(self, features: pd.DataFrame) -> Dict[str, List]:
        """
        Combine attributes from multiple features.

        Args:
            features: DataFrame containing features to combine

        Returns:
            Dictionary of combined attributes
        """
        # For ID, use the first one with a _merged suffix
        first_id = features.iloc[0]['ID']
        combined_attrs = {
            'ID': [f"{first_id}_merged"],
            'Parent': [features.iloc[0].get('Parent', '')]
        }

        # For other attributes, collect all unique values
        for col in features.columns:
            if col not in ['Chromosome', 'Start', 'End', 'Strand', 'ID', 'Parent']:
                unique_vals = features[col].dropna().unique()
                if len(unique_vals) > 0:
                    combined_attrs[col] = [';'.join(map(str, unique_vals))]

        return combined_attrs