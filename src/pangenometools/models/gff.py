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
import sys
from typing import Union
from pandas import IndexSlice as idx

class GFFHandler(PangenomeFileHandler):
    """Handler for GFF/GTF files in pangenome."""

    _gff: Union[None, pr.PyRanges] = None
    _gff_feature: Union[None, pd.DataFrame] = None
    _gff_all_features: Union[None, pd.DataFrame] = None
    _genotype: Union[None, str] = None
    _feature: Union[None, str] = None

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize GFF handler.

        Args:
            pangenome_folder: Path to pangenome folder
            pangenome_index: Path to pangenome index file
        """
        super().__init__(pangenome_folder, pangenome_index)

    def _load_gff(self, genotype) -> pr.PyRanges:
        # See if the required fasta is still open
        if self._gff is None or self._genotype != genotype:

            # Close the open GFF
            if self._genotype != genotype:
                self.close_gff()

            # Open the GFF if not
            gff_path = self.resolve_path(genotype, "annotation")

            self._gff = pr.read_gff3(gff_path)
            self._genotype = genotype
        
        return self._gff

    def close_gff(self):
        # Reset the GFF information
        if self._gff is not None:

            self._gff = self._genotype = None

    def get_feature(self, genotype: str, feature_type: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load GFF file for a genotype, optionally filtered by feature type.

        Args:
            genotype: Genotype identifier
            feature_type: Optional feature type to filter by

        Returns:
            DataFrame containing GFF features
        """
        # Use the cached rage if it matches
        if (self._genotype == genotype and self._feature == feature_type) and not self._gff_feature is None and not self._gff_all_features is None:
            return self._gff_feature, self._gff_all_features
        else:
            # Other wise, load the GFF - use cached
            features = self._load_gff(genotype).df

            # Set the loaded genotype and feaute
            self._genotype = genotype
            self._feature = feature_type

            # Determine where the feature types is as requested
            feature_type_bool = (features.Feature == feature_type).to_numpy()

            # Set an index for easier searching
            _idx = pd.MultiIndex.from_frame(features.loc[:,["ID", "Parent"]])
            features.index = _idx

            # Save the full gff features
            self._gff_all_features = features

            # Subset the features if requested
            if feature_type == "all":
                feature = features
            else:
                feature = features[feature_type_bool]
            
            if feature.shape[0] == 0:
                raise ValueError(f"Feature type {feature_type} not found in annotation")

            self._gff_feature = feature

            return self._gff_feature, self._gff_all_features

    def get_features_for_gene(self,
                            genotype: str,
                            gene_id: str,
                            feature_type: str = "CDS",
                            merge_strategy: Literal["merge", "first", "all"] = "merge"
                            ) -> Union[pd.DataFrame, List[pd.DataFrame]]:
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
        all_selected, all_features = self.get_feature(genotype, feature_type)

        # Find all features of the requested type that are children or grandchildren of this gene
        if feature_type == "gene":
            gene_children = all_selected.loc[idx[gene_id, :]]
        else:
            # Get the IDs of children of this gene
            # children_ids = set(all_features.ID[all_features.Parent == gene_id])
            children_ids = set(all_features.loc[idx[:, gene_id], "ID"])

            gene_children = all_selected.loc[(all_selected.Parent == gene_id) | all_selected.Parent.isin(children_ids)]

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

    def _merge_features(self, features: pd.DataFrame) -> pd.DataFrame:
        """
        Merge multiple features into a single feature spanning the entire range.

        For features on the same strand and chromosome:
        - Uses the minimum start position
        - Uses the maximum end position
        - Preserves the strand
        - Combines attributes

        Args:
            features: DataFrame containing features to merge

        Returns:
            DataFrame with merged features
        """
        if len(features) == 1:
            return features

        # Group by chromosome and strand
        grouped = features.groupby(['Chromosome', 'Strand'])

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

            merged_features.append(new_feature)

        if not merged_features:
            return pd.DataFrame()  # Empty result

        return pd.concat([m for m in merged_features], ignore_index=True)

    def _split_by_strand_and_chrom(self, features: pd.DataFrame) -> List[pd.DataFrame]:
        """
        Split features by chromosome and strand, returning a list.

        Args:
            features: DataFrame to split

        Returns:
            List of DataFrames, one per chromosome/strand combination
        """
        grouped = features.groupby(['Chromosome', 'Strand'])
        return [group for _, group in grouped]

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