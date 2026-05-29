"""
Tests specifically for GFF CLI functionality.
"""

import pytest
import tempfile
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
from pangenometools.cli.gff_cli import setup_gff_parser, gff_cli

def test_gff_cli_parser_setup():
    """Test GFF CLI parser setup."""
    parser = setup_gff_parser()
    assert parser is not None
    
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parser.parse_args([])

def test_gff_cli_parser_arguments():
    """Test GFF CLI parser argument parsing."""
    parser = setup_gff_parser()
    
    # Test with all required arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.gff",
        "--feature-type", "CDS",
        "--search-mode", "children",
        "--merge-strategy", "merge",
        "--silent"
    ])
    
    assert args.pangenome_folder == "/test/folder"
    assert args.pangenome_index == "/test/index.csv"
    assert args.target_genes == "/test/genes.csv"
    assert args.output == "/test/output.gff"
    assert args.feature_type == "CDS"
    assert args.search_mode == "children"
    assert args.merge_strategy == "merge"
    assert args.silent == True

def test_gff_cli_error_handling():
    """Test GFF CLI error handling."""
    # Test with missing required arguments
    with patch.object(sys, 'argv', ['gff_cli']):
        with pytest.raises(SystemExit):
            gff_cli()

def test_gff_cli_help_output(capsys):
    """Test GFF CLI help output."""
    # Test help output
    with patch.object(sys, 'argv', ['gff_cli', '--help']):
        with pytest.raises(SystemExit):
            gff_cli()

def test_gff_cli_with_test_data(setup_test_files):
    """Test GFF CLI with actual test data."""
    output_file = setup_test_files / "output.gff"
    
    # Test CLI execution
    with patch.object(sys, 'argv', [
        'gff_cli',
        '--pangenome-folder', str(setup_test_files),
        '--pangenome-index', str(setup_test_files / 'index.csv'),
        '--target-genes', str(setup_test_files / 'target_genes.csv'),
        '--output', str(output_file),
        '--feature-type', 'gene',
        '--silent'
    ]):
                gff_cli()
    
    # Verify output file was created and has content
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    
    # Debug: print file size
    file_size = output_file.stat().st_size
    print(f"GFF CLI output file size: {file_size}")

def test_gff_cli_different_search_modes(setup_test_files):
    """Test GFF CLI with different search modes."""
    # Create test files
    target_content = """gene_name,gene_ID_test1
GENE001,GENE001"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test with different search modes
    search_modes = ["strict", "children", "pattern"]
    
    for search_mode in search_modes:
        output_file = setup_test_files / f"output_{search_mode}.gff"
        
        with patch.object(sys, 'argv', [
            'gff_cli',
            '--pangenome-folder', str(setup_test_files),
            '--pangenome-index', str(setup_test_files / 'index.csv'),
            '--target-genes', str(setup_test_files / 'target_genes.csv'),
            '--output', str(output_file),
            '--feature-type', 'gene',
            '--search-mode', search_mode,
            '--silent'
        ]):
                    gff_cli()
        
        # Verify output file was created and has content
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        print(f"GFF CLI {search_mode} mode output file size: {output_file.stat().st_size}")

def test_gff_cli_different_merge_strategies(setup_test_files):
    """Test GFF CLI with different merge strategies."""
    # Create test files
    target_content = """gene_name,gene_ID_test1
GENE001,GENE001"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test with different merge strategies
    merge_strategies = ["merge", "first", "all"]
    
    for merge_strategy in merge_strategies:
        print(f"Merge strategy: {merge_strategy}")
        output_file = setup_test_files / f"output_{merge_strategy}.gff"
        
        with patch.object(sys, 'argv', [
            'gff_cli',
            '--pangenome-folder', str(setup_test_files),
            '--pangenome-index', str(setup_test_files / 'index.csv'),
            '--target-genes', str(setup_test_files / 'target_genes.csv'),
            '--output', str(output_file),
            '--feature-type', 'CDS',
            '--merge-strategy', merge_strategy,
            '--silent'
        ]):
                    gff_cli()
        
        # Verify output file was created and has content
        assert output_file.exists()
        assert output_file.stat().st_size > 0
        print(f"GFF CLI {merge_strategy} strategy output file size: {output_file.stat().st_size}")

def test_gff_cli_invalid_arguments():
    """Test GFF CLI with invalid arguments."""
    parser = setup_gff_parser()
    
    # Test with invalid search mode
    with pytest.raises(SystemExit):
        parser.parse_args([
            "--pangenome-folder", "/test/folder",
            "--pangenome-index", "/test/index.csv",
            "--target-genes", "/test/genes.csv",
            "--output", "/test/output.gff",
            "--search-mode", "invalid"
        ])
    
    # Test with invalid merge strategy
    with pytest.raises(SystemExit):
        parser.parse_args([
            "--pangenome-folder", "/test/folder",
            "--pangenome-index", "/test/index.csv",
            "--target-genes", "/test/genes.csv",
            "--output", "/test/output.gff",
            "--merge-strategy", "invalid"
        ])

def test_gff_cli_argument_types():
    """Test GFF CLI argument type handling."""
    parser = setup_gff_parser()
    
    # Test with valid arguments
    args = parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.gff"
    ])
    
    assert isinstance(args.pangenome_folder, str)
    assert isinstance(args.pangenome_index, str)
    assert isinstance(args.target_genes, str)
    assert isinstance(args.output, str)
    assert isinstance(args.feature_type, str)
    assert isinstance(args.search_mode, str)
    assert isinstance(args.merge_strategy, str)
    assert isinstance(args.silent, bool)