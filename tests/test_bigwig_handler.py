"""
Tests for BigWig handler functionality.
"""

import pytest
import tempfile
import json
from pathlib import Path
from pangenometools.models.bigwig import BigWigHandler


def test_bigwig_handler_initialization(setup_test_files):
    """Test BigWig handler initialization."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")
    assert handler is not None
    assert handler.pangenome_folder == setup_test_files
    assert len(handler.pangenome_index) > 0


def test_bigwig_handler_read_bigwig_index(setup_test_files):
    """Test reading BigWig index."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")
    
    # Create a temporary BigWig index file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("genotype,condition,bigwig\n")
        f.write("test1,condition1,test1.bw\n")
        f.write("test2,condition2,test2.bw\n")
        temp_path = Path(f.name)
    
    try:
        index = handler.read_bigwig_index(temp_path)
        assert len(index) == 2
        assert index[0]["genotype"] == "test1"
        assert index[0]["condition"] == "condition1"
        assert index[0]["bigwig"] == "test1.bw"
    finally:
        temp_path.unlink()

def test_bigwig_handler_extraction_with_test_data(setup_test_files):
    """Test BigWig extraction with actual test data."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")

    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test extraction
    outputs = handler.extract_bigwig_signals(
        bigwig_folder=bigwig_folder,
        bigwig_index_path=setup_test_files / "bigwig_index.csv",
        target_genes_path=setup_test_files / "target_genes.csv",
        feature_type="gene",
        upstream=50,
        downstream=50,
        inner_start=100,
        inner_end=100,
        pad=10
    )

    # The extraction should run without errors (though it may not find valid data)
    assert isinstance(outputs, list)

def test_bigwig_handler_extraction_and_save(setup_test_files):
    """Test BigWig extraction and save functionality."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")

    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test extraction and save
    output_file = setup_test_files / "bigwig_output.json"

    outputs = handler.extract_and_save_bigwig_signals(
        bigwig_folder=bigwig_folder,
        bigwig_index_path=setup_test_files / "bigwig_index.csv",
        target_genes_path=setup_test_files / "target_genes.csv",
        output_path=output_file,
        feature_type="gene",
        upstream=50,
        downstream=50,
        inner_start=100,
        inner_end=100,
        pad=10
    )
    
    # Check that output file was created
    assert output_file.exists()
    
    # The file should contain valid JSON (even if empty)
    with open(output_file, "r") as f:
        data = json.load(f)
        assert isinstance(data, list)

def test_bigwig_handler_different_feature_types(setup_test_files):
    """Test BigWig extraction with different feature types."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")
    
    # Create BigWig folder and index
    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test with different feature types
    feature_types = ["gene", "CDS"]
    
    for feature_type in feature_types:
        outputs = handler.extract_bigwig_signals(
            bigwig_folder=bigwig_folder,
            bigwig_index_path=setup_test_files / "bigwig_index.csv",
            target_genes_path=setup_test_files / "target_genes.csv",
            feature_type=feature_type,
            upstream=50,
            downstream=50,
            inner_start=100,
            inner_end=100,
            pad=10
        )
        
        # Should run without errors
        assert isinstance(outputs, list)

def test_bigwig_handler_with_five_prime_direction(setup_test_files):
    """Test BigWig extraction with five prime direction option."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")

    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test extraction
    outputs = handler.extract_bigwig_signals(
        bigwig_folder=bigwig_folder,
        bigwig_index_path=setup_test_files / "bigwig_index.csv",
        target_genes_path=setup_test_files / "target_genes.csv",
        feature_type="gene",
        upstream=50,
        downstream=50,
        inner_start=100,
        inner_end=100,
        pad=10,
        use_five_prime_direction=True
    )

    # The extraction should run without errors (though it may not find valid data)
    assert isinstance(outputs, list)


def test_bigwig_handler_with_whole_seq(setup_test_files):
    """Test BigWig extraction with whole sequence option."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")

    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test extraction
    outputs = handler.extract_bigwig_signals(
        bigwig_folder=bigwig_folder,
        bigwig_index_path=setup_test_files / "bigwig_index.csv",
        target_genes_path=setup_test_files / "target_genes.csv",
        feature_type="gene",
        upstream=50,
        downstream=50,
        inner_start=100,
        inner_end=100,
        pad=10,
        whole_seq=True
    )

    # The extraction should run without errors (though it may not find valid data)
    assert isinstance(outputs, list)

def test_bigwig_handler_extraction_with_test_data2(setup_test_files):
    """Test BigWig extraction with actual test data."""
    handler = BigWigHandler(setup_test_files, setup_test_files / "index.csv")

    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    # Test extraction
    outputs = handler.extract_bigwig_signals(
        bigwig_folder=bigwig_folder,
        bigwig_index_path=setup_test_files / "bigwig_index.csv",
        target_genes_path=setup_test_files / "target_genes.csv",
        feature_type="gene",
        upstream=50,
        downstream=50,
        inner_start=100,
        inner_end=100,
        pad=10
    )

    # The extraction should run without errors (though it may not find valid data)
    assert isinstance(outputs, list)