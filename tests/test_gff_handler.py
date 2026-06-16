"""
Tests for GFFHandler class.
"""

import pytest
import pandas as pd
import pyranges as pr
from pathlib import Path
from pangenometools.models.gff import GFFHandler

class TestGFFHandler:
    """Test cases for GFFHandler."""

    def test_load_gff(self, setup_test_files):
        """Test basic GFF loading."""
        handler = GFFHandler(setup_test_files, setup_test_files / "index.csv")
        features = handler.load_gff("test1")
        
        assert len(features) > 0
        assert "GENE001" in features.ID.values

    def test_get_features_for_gene_cds_merging(self, setup_test_files):
        """Test CDS feature merging for a gene."""
        handler = GFFHandler(setup_test_files, setup_test_files / "index.csv")
        merged_cds = handler.get_features_for_gene("test1", "GENE001", "CDS", "merge")
        
        assert len(merged_cds) == 1
        # Note: pyranges uses 0-based indexing (start included, end excluded)
        # So our original coordinates get adjusted
        assert merged_cds.Start[0] == 149  # Min start (150-1)
        assert merged_cds.End[0] == 450    # Max end (already exclusive)
        assert merged_cds.ID[0] == "CDS001_merged"

    def test_strand_specific_merging(self, setup_test_files):
        """Test that features on different strands aren't merged together."""
        handler = GFFHandler(setup_test_files, setup_test_files / "index.csv")

        # Get all CDS features for both genes
        all_cds = handler.get_features_for_gene("test1", "GENE001", "CDS", "all")
        assert len(all_cds) == 1  # Should be grouped by strand
        
        # Test merging for each gene separately
        merged_cds1 = handler.get_features_for_gene("test1", "GENE001", "CDS", "merge")
        merged_cds2 = handler.get_features_for_gene("test1", "GENE002", "CDS", "merge")
        
        assert len(merged_cds1) == 1
        assert len(merged_cds2) == 1
        assert merged_cds1.Strand[0] == "+"
        assert merged_cds2.Strand[0] == "-"

    def test_gene_not_found(self, setup_test_files):
        """Test handling of non-existent genes."""
        handler = GFFHandler(setup_test_files, setup_test_files / "index.csv")
        # Test non-existent gene
        result = handler.get_features_for_gene("test1", "NONEXISTENT", "CDS", "merge")
        assert len(result) == 0

    def test_first_strategy(self, setup_test_files):
        """Test 'first' merge strategy."""
        handler = GFFHandler(setup_test_files, setup_test_files / "index.csv")

        # Test first strategy
        first_cds = handler.get_features_for_gene("test1", "GENE001", "CDS", "first")
        assert len(first_cds) == 1
        # Access the ID correctly - pyranges might return a different structure
        assert first_cds.iloc[0]['ID'] == "CDS001"  # Should be the first CDS

