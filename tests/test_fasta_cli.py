"""
Tests specifically for FASTA CLI functionality.
"""

import pytest
import tempfile
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
from pangenometools.cli.fasta_cli import setup_fasta_parser, fasta_cli

def test_fasta_cli_parser_setup():
    """Test FASTA CLI parser setup."""
    parser = setup_fasta_parser()
    assert parser is not None
    
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parser.parse_args([])

def test_fasta_cli_parser_arguments():
    """Test FASTA CLI parser argument parsing."""
    parser = setup_fasta_parser()
    
    # Test with all required arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.fa",
        "--feature-type", "gene",
        "--upstream", "100",
        "--downstream", "50",
        "--inner-start", "25",
        "--inner-end", "25",
        "--pad", "5",
        "--merge-strategy", "merge",
        "--use-five-prime-direction",
        "--silent",
        "--include-coordinates"
    ])
    
    assert args.pangenome_folder == "/test/folder"
    assert args.pangenome_index == "/test/index.csv"
    assert args.target_genes == "/test/genes.csv"
    assert args.output == "/test/output.fa"
    assert args.feature_type == "gene"
    assert args.upstream == 100
    assert args.downstream == 50
    assert args.inner_start == 25
    assert args.inner_end == 25
    assert args.pad == 5
    assert args.merge_strategy == "merge"
    assert args.use_five_prime_direction == True
    assert args.silent == True
    assert args.include_coordinates == True

def test_fasta_cli_error_handling():
    """Test FASTA CLI error handling."""
    # Test with missing required arguments
    with patch.object(sys, 'argv', ['fasta_cli']):
        with pytest.raises(SystemExit):
            fasta_cli()

def test_fasta_cli_help_output(capsys):
    """Test FASTA CLI help output."""
    # Test help output
    with patch.object(sys, 'argv', ['fasta_cli', '--help']):
        with pytest.raises(SystemExit):
            fasta_cli()

def test_fasta_cli_with_test_data(setup_test_files):
    """Test FASTA CLI with actual test data."""
    output_file = setup_test_files / "output.fa"
    
    # Test CLI execution
    with patch.object(sys, 'argv', [
        'fasta_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--inner-start', '50',
        '--inner-end', '50',
        '--silent'
    ]):
                fasta_cli()
    
    # Verify output file was created and has content
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    
    # Debug: print file size
    file_size = output_file.stat().st_size
    print(f"FASTA CLI output file size: {file_size}")

def test_fasta_cli_different_feature_types(setup_test_files):
    """Test FASTA CLI with different feature types."""
    
    # Test with different feature types
    feature_types = ["gene", "CDS"]
    
    for feature_type in feature_types:
        output_file = setup_test_files / f"output_{feature_type}.fa"
        
        with patch.object(sys, 'argv', [
            'fasta_cli',
            '--pangenome-folder', str(setup_test_files),
            '--pangenome-index', str(setup_test_files / 'index.csv'),
            '--target-genes', str(setup_test_files / 'target_genes.csv'),
            '--output', str(output_file),
            '--feature-type', feature_type,
            '--merge-strategy', 'merge',
            '--inner-start', '25',
            '--inner-end', '25',
            '--silent'
        ]):
                    fasta_cli()
        
        # Verify output file was created and has content
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        print(f"FASTA CLI {feature_type} feature type output file size: {output_file.stat().st_size}")

def test_fasta_cli_different_merge_strategies(setup_test_files):
    """Test FASTA CLI with different merge strategies."""
    
    # Test with different merge strategies
    merge_strategies = ["merge", "first", "all"]
    
    for merge_strategy in merge_strategies:
        output_file = setup_test_files / f"output_{merge_strategy}.fa"
        
        with patch.object(sys, 'argv', [
            'fasta_cli',
            '--pangenome-folder', str(setup_test_files),
            '--pangenome-index', str(setup_test_files / 'index.csv'),
            '--target-genes', str(setup_test_files / 'target_genes.csv'),
            '--output', str(output_file),
            '--feature-type', 'CDS',
            '--merge-strategy', merge_strategy,
            '--inner-start', '10',
            '--inner-end', '10',
            '--silent'
        ]):
                    fasta_cli()
        
        # Verify output file was created and has content
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        print(f"FASTA CLI {merge_strategy} strategy output file size: {output_file.stat().st_size}")

def test_fasta_cli_with_upstream_downstream(setup_test_files):
    """Test FASTA CLI with upstream and downstream options."""
    
    output_file = setup_test_files / "output_with_flanking.fa"
    
    with patch.object(sys, 'argv', [
        'fasta_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--upstream', '10',
        '--downstream', '10',
        '--inner-start', '25',
        '--inner-end', '25',
        '--silent'
    ]):
        fasta_cli()
        
    # Verify output file was created and has content
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    print(f"FASTA CLI with flanking regions output file size: {output_file.stat().st_size}")

def test_fasta_cli_with_coordinates_option(setup_test_files):
    """Test FASTA CLI with coordinates option."""
    
    output_file = setup_test_files / "output_with_coords.fa"
    
    with patch.object(sys, 'argv', [
        'fasta_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--inner-start', '50',
        '--inner-end', '50',
        '--include-coordinates',
        '--silent'
    ]):
        fasta_cli()
        
    # Verify output file was created and has content
    assert output_file.exists()
    assert output_file.stat().st_size > 0
        
    # Debug: print file size and content
    file_size = output_file.stat().st_size
    print(f"FASTA CLI with coordinates output file size: {file_size}")
        
    content = output_file.read_text()
    print(f"FASTA CLI with coordinates output content: {content[:200]}...")
    assert "chrom=" in content
    assert "start=" in content
    assert "end=" in content
    assert "strand=" in content

def test_fasta_cli_with_five_prime_direction(setup_test_files):
    """Test FASTA CLI with five prime direction option."""
    
    # Test with use-five-prime-direction
    output_file_5prime = setup_test_files / "output_5prime.fa"
    
    with patch.object(sys, 'argv', [
        'fasta_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file_5prime),
        '--feature-type', 'gene',
        '--upstream', '10',
        '--downstream', '5',
        '--inner-start', '25',
        '--inner-end', '25',
        '--use-five-prime-direction',
        '--silent'
    ]):
                fasta_cli()
    
    # Test without use-five-prime-direction
    output_file_default = setup_test_files / "output_default.fa"
    
    with patch.object(sys, 'argv', [
        'fasta_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file_default),
        '--feature-type', 'gene',
        '--upstream', '10',
        '--downstream', '5',
        '--inner-start', '25',
        '--inner-end', '25',
        '--silent'
    ]):
                fasta_cli()
    
    # Verify both output files were created and have content
    assert output_file_5prime.exists()
    assert output_file_default.exists()
    assert output_file_5prime.stat().st_size > 0
    assert output_file_default.stat().st_size > 0
    
    print(f"FASTA CLI 5prime direction output file size: {output_file_5prime.stat().st_size}")
    print(f"FASTA CLI default direction output file size: {output_file_default.stat().st_size}")

def test_fasta_cli_invalid_arguments():
    """Test FASTA CLI with invalid arguments."""
    parser = setup_fasta_parser()
    
    # Test with invalid merge strategy
    with pytest.raises(SystemExit):
        parser.parse_args([
            "--pangenome-folder", "/test/folder",
            "--pangenome-index", "/test/index.csv",
            "--target-genes", "/test/genes.csv",
            "--output", "/test/output.fa",
            "--merge-strategy", "invalid"
        ])

def test_fasta_cli_argument_types():
    """Test FASTA CLI argument type handling."""
    parser = setup_fasta_parser()
    
    # Test with valid arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.fa",
        "--upstream", "100",
        "--downstream", "50",
        "--inner-start", "25",
        "--inner-end", "25",
        "--pad", "5"
    ])
    
    assert isinstance(args.pangenome_folder, str)
    assert isinstance(args.pangenome_index, str)
    assert isinstance(args.target_genes, str)
    assert isinstance(args.output, str)
    assert isinstance(args.feature_type, str)
    assert isinstance(args.upstream, int)
    assert isinstance(args.downstream, int)
    assert isinstance(args.inner_start, int)
    assert isinstance(args.inner_end, int)
    assert isinstance(args.pad, int)
    assert isinstance(args.merge_strategy, str)
    assert isinstance(args.use_five_prime_direction, bool)
    assert isinstance(args.silent, bool)
    assert isinstance(args.include_coordinates, bool)