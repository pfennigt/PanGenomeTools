"""
Tests for FASTA handler functionality.
"""

import pytest
from pathlib import Path
import numpy as np
from src.pangenometools.models.fasta import FastaHandler
from src.pangenometools.models.gff import GFFHandler

def test_fasta_handler_initialization(setup_test_files):
    """Test FASTA handler initialization."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    # Test successful initialization
    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)
    assert fasta_handler is not None
    assert fasta_handler.gff_handler is not None
    assert isinstance(fasta_handler.gff_handler, GFFHandler)

def test_fasta_handler_invalid_index(setup_test_files):
    """Test FASTA handler with invalid index file."""
    pangenome_folder = setup_test_files
    invalid_index = pangenome_folder / "nonexistent.csv"

    with pytest.raises(FileNotFoundError):
        FastaHandler(pangenome_folder, invalid_index)

def test_fasta_handler_missing_genotype(setup_test_files):
    """Test FASTA handler with missing genotype."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    with pytest.raises(ValueError):
        fasta_handler.extract_sequence("nonexistent", "GENE001")

def test_fasta_handler_missing_gene(setup_test_files):
    """Test FASTA handler with missing gene."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    with pytest.raises(ValueError):
        fasta_handler.extract_sequence("test1", "NONEXISTENT")

def test_fasta_handler_basic_extraction(setup_test_files):
    """Test basic sequence extraction."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test gene extraction for gene 1
    sequence = fasta_handler.extract_sequence("test1", "GENE001", "gene", upstream=20, downstream=40, pad=0)
    assert len(sequence) > 0
    assert isinstance(sequence, str)
    assert sequence == "TTTTTAAAAAAACCCCCCCGCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTT"

    # Test gene extraction for gene 2
    sequence = fasta_handler.extract_sequence("test1", "GENE002", "gene", upstream=20, downstream=40, pad=0)
    assert len(sequence) > 0
    assert isinstance(sequence, str)
    assert sequence == "TTTTTTTTAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGTT"

def test_fasta_handler_cds_extraction(setup_test_files):
    """Test CDS sequence extraction."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test CDS extraction
    sequence = fasta_handler.extract_sequence("test1", "GENE001", "CDS", upstream=10, downstream=10)
    assert len(sequence) > 0
    assert isinstance(sequence, str)

def test_fasta_handler_negative_strand(setup_test_files):
    """Test sequence extraction for negative strand genes."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test negative strand gene
    sequence = fasta_handler.extract_sequence("test1", "GENE002", "gene", upstream=10, downstream=10)
    assert len(sequence) > 0

def test_fasta_handler_use_five_prime_direction(setup_test_files):
    """Test use_five_prime_direction parameter."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test with use_five_prime_direction=True
    sequence_5prime = fasta_handler.extract_sequence(
        "test1", "GENE002", "gene", upstream=20, downstream=10,
        use_five_prime_direction=True
    )

    # Test with use_five_prime_direction=False (default)
    sequence_default = fasta_handler.extract_sequence(
        "test1", "GENE002", "gene", upstream=20, downstream=10,
        use_five_prime_direction=False
    )

    # These should be different for negative strand genes
    assert sequence_5prime != sequence_default

def test_fasta_handler_inner_regions(setup_test_files):
    """Test sequence extraction with inner start and end regions."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test with inner regions
    sequence = fasta_handler.extract_sequence(
        "test1", "GENE001", "gene", upstream=10, downstream=10, inner_start=10, inner_end=10
    )
    assert len(sequence) > 0

def test_fasta_handler_merge_strategies(setup_test_files):
    """Test different merge strategies."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test merge strategy
    sequence_merge = fasta_handler.extract_sequence(
        "test1", "GENE001", "CDS", merge_strategy="merge", upstream=10, downstream=10
    )

    # Test first strategy
    sequence_first = fasta_handler.extract_sequence(
        "test1", "GENE001", "CDS", merge_strategy="first", upstream=10, downstream=10
    )

    # Test all strategy
    sequence_all = fasta_handler.extract_sequence(
        "test1", "GENE001", "CDS", merge_strategy="all", upstream=10, downstream=10
    )

    # All should return valid sequences
    assert len(sequence_merge) > 0
    assert len(sequence_first) > 0
    assert len(sequence_all) > 0

def test_fasta_handler_edge_cases(setup_test_files):
    """Test edge cases and boundary conditions."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test with zero upstream/downstream
    sequence = fasta_handler.extract_sequence(
        "test1", "GENE001", "gene", upstream=0, downstream=0
    )
    assert len(sequence) == 0

    # Test with large upstream/downstream (should be clipped)
    sequence = fasta_handler.extract_sequence(
        "test1", "GENE001", "gene", upstream=1000, downstream=1000
    )
    assert len(sequence) > 0

def test_fasta_handler_padding(setup_test_files):
    """Test sequence extraction with padding."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Test with padding
    sequence_with_pad = fasta_handler.extract_sequence(
        "test1", "GENE001", "gene", pad=5, upstream=10, downstream=10
    )

    sequence_without_pad = fasta_handler.extract_sequence(
        "test1", "GENE001", "gene", pad=0, upstream=10, downstream=10
    )

    # With padding should be longer
    assert len(sequence_with_pad) >= len(sequence_without_pad)

def test_fasta_handler_sequence_content(setup_test_files):
    """Test that extracted sequences contain expected content."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Extract sequence and check it contains expected nucleotides
    sequence = fasta_handler.extract_sequence("test1", "GENE001", "gene", upstream=10, downstream=10)
    assert any(nucleotide in sequence for nucleotide in ['A', 'T', 'C', 'G'])

def test_fasta_handler_strand_specific_extraction(setup_test_files):
    """Test strand-specific sequence extraction."""
    pangenome_folder = setup_test_files
    pangenome_index = pangenome_folder / "index.csv"

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Extract from positive strand gene
    pos_sequence = fasta_handler.extract_sequence("test1", "GENE001", "gene", upstream=10, downstream=10)

    # Extract from negative strand gene
    neg_sequence = fasta_handler.extract_sequence("test1", "GENE002", "gene", upstream=10, downstream=10)

    # Both should be valid sequences
    assert len(pos_sequence) > 0
    assert len(neg_sequence) > 0
