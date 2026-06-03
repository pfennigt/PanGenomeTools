"""
Tests specifically for BigWig CLI functionality.
"""

import pytest
import tempfile
import sys
from pathlib import Path
from unittest.mock import patch
from pangenometools.cli.bigwig_cli import setup_bigwig_parser, bigwig_cli


def test_bigwig_cli_parser_setup():
    """Test BigWig CLI parser setup."""
    parser = setup_bigwig_parser()
    assert parser is not None
    
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_bigwig_cli_parser_arguments():
    """Test BigWig CLI parser argument parsing."""
    parser = setup_bigwig_parser()
    
    # Test with all required arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--bigwig-folder", "/test/bigwigs",
        "--bigwig-index", "/test/bigwig_index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.json",
        "--feature-type", "CDS",
        "--upstream", "100",
        "--downstream", "50",
        "--inner-start", "25",
        "--inner-end", "25",
        "--pad", "5",
        "--whole-seq",
        "--use-five-prime-direction",
        "--silent"
    ])
    
    assert args.pangenome_folder == "/test/folder"
    assert args.pangenome_index == "/test/index.csv"
    assert args.bigwig_folder == "/test/bigwigs"
    assert args.bigwig_index == "/test/bigwig_index.csv"
    assert args.target_genes == "/test/genes.csv"
    assert args.output == "/test/output.json"
    assert args.feature_type == "CDS"
    assert args.upstream == 100
    assert args.downstream == 50
    assert args.inner_start == 25
    assert args.inner_end == 25
    assert args.pad == 5
    assert args.whole_seq == True
    assert args.use_five_prime_direction == True
    assert args.silent == True


def test_bigwig_cli_error_handling():
    """Test BigWig CLI error handling."""
    # Test with missing required arguments
    with patch.object(sys, 'argv', ['bigwig_cli']):
        with pytest.raises(SystemExit):
            bigwig_cli()


def test_bigwig_cli_help_output(capsys):
    """Test BigWig CLI help output."""
    # Test help output
    with patch.object(sys, 'argv', ['bigwig_cli', '--help']):
        with pytest.raises(SystemExit):
            bigwig_cli()


def test_bigwig_cli_with_test_data(setup_test_files):
    """Test BigWig CLI with actual test data."""
    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    output_file = setup_test_files / "bigwig_output.json"
    
    # Test CLI execution
    with patch.object(sys, 'argv', [
        'bigwig_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--bigwig-folder', str(bigwig_folder),
        '--bigwig-index', str(setup_test_files / "bigwig_index.csv"),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--inner-start', '50',
        '--inner-end', '50',
        '--silent'
    ]):
        bigwig_cli()
    
    # Verify output file was created
    assert output_file.exists()

def test_bigwig_cli_different_feature_types(setup_test_files):
    """Test BigWig CLI with different feature types."""
    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    output_file = setup_test_files / "bigwig_output.json"
    
    # Test with different feature types
    feature_types = ["gene", "CDS"]
    
    # Test CLI execution with different gene features
    for feature_type in feature_types:
        print(f"trying {feature_type}")
        output_file = setup_test_files / f"bigwig_output_{feature_type}.json"
        
        with patch.object(sys, 'argv', [
            'bigwig_cli',
            '--pangenome-folder', str(setup_test_files),
            '--pangenome-index', str(setup_test_files / 'index.csv'),
            '--bigwig-folder', str(bigwig_folder),
            '--bigwig-index', str(setup_test_files / "bigwig_index.csv"),
            '--target-genes', str(setup_test_files / 'target_genes.csv'),
            '--output', str(output_file),
            '--feature-type', 'gene',
            '--inner-start', '50',
            '--inner-end', '50',
            '--silent'
        ]):
            bigwig_cli()
        
        # Verify output file was created
        assert output_file.exists()

def test_bigwig_cli_with_five_prime_direction(setup_test_files):
    """Test BigWig CLI with five prime direction option."""
    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    output_file = setup_test_files / "bigwig_output.json"
    
    # Test CLI execution with forced five-prime direction
    with patch.object(sys, 'argv', [
        'bigwig_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--bigwig-folder', str(bigwig_folder),
        '--bigwig-index', str(setup_test_files / "bigwig_index.csv"),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--inner-start', '50',
        '--inner-end', '50',
        '--use-five-prime-direction',
        '--silent'
    ]):
        bigwig_cli()
    
    # Verify output file was created
    assert output_file.exists()

def test_bigwig_cli_with_whole_seq(setup_test_files):
    """Test BigWig CLI with whole sequence option."""
    bigwig_folder = setup_test_files / "bigwigs"
    bigwig_folder.mkdir(exist_ok=True)
    
    output_file = setup_test_files / "bigwig_output.json"
    
    # Test CLI execution with whole sequence output
    with patch.object(sys, 'argv', [
        'bigwig_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--bigwig-folder', str(bigwig_folder),
        '--bigwig-index', str(setup_test_files / "bigwig_index.csv"),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--inner-start', '50',
        '--inner-end', '50',
        '--use-five-prime-direction',
        '--silent'
    ]):
        bigwig_cli()
    
    # Verify output file was created
    assert output_file.exists()

def test_bigwig_cli_argument_types():
    """Test BigWig CLI argument type handling."""
    parser = setup_bigwig_parser()
    
    # Test with valid arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--bigwig-folder", "/test/bigwigs",
        "--bigwig-index", "/test/bigwig_index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.json",
        "--upstream", "100",
        "--downstream", "50",
        "--inner-start", "25",
        "--inner-end", "25",
        "--pad", "5"
    ])
    
    assert isinstance(args.pangenome_folder, str)
    assert isinstance(args.pangenome_index, str)
    assert isinstance(args.bigwig_folder, str)
    assert isinstance(args.bigwig_index, str)
    assert isinstance(args.target_genes, str)
    assert isinstance(args.output, str)
    assert isinstance(args.feature_type, str)
    assert isinstance(args.upstream, int)
    assert isinstance(args.downstream, int)
    assert isinstance(args.inner_start, int)
    assert isinstance(args.inner_end, int)
    assert isinstance(args.pad, int)
    assert isinstance(args.whole_seq, bool)
    assert isinstance(args.use_five_prime_direction, bool)
    assert isinstance(args.silent, bool)