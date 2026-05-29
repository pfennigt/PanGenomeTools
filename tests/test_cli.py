"""
Tests for CLI modules.
"""

import pytest
import tempfile
import csv
from pathlib import Path
from unittest.mock import patch, MagicMock
from src.pangenometools.cli.gff_cli import setup_gff_parser, read_target_genes as read_gff_target_genes
from src.pangenometools.cli.fasta_cli import setup_fasta_parser, read_target_genes as read_fasta_target_genes

def test_gff_parser_setup():
    """Test GFF parser setup."""
    parser = setup_gff_parser()
    assert parser is not None
    
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parser.parse_args([])

def test_fasta_parser_setup():
    """Test FASTA parser setup."""
    parser = setup_fasta_parser()
    assert parser is not None
    
    # Test that required arguments are present
    with pytest.raises(SystemExit):
        parser.parse_args([])

def test_gff_target_genes_parsing(setup_test_files):
    """Test GFF target genes parsing."""
    # Create a test target genes file
    target_content = """gene_name,gene_ID_test1
GENE001,GENE001
GENE002,GENE002"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test parsing
    genotypes, rows = read_gff_target_genes(target_file)
    
    assert len(genotypes) == 1
    assert "test1" in genotypes
    assert len(rows) == 2

def test_fasta_target_genes_parsing(setup_test_files):
    """Test FASTA target genes parsing."""
    # Create a test target genes file
    target_content = """gene_name,gene_ID_test1
GENE001,GENE001
GENE002,GENE002"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test parsing
    genotypes, rows = read_fasta_target_genes(target_file)
    
    assert len(genotypes) == 1
    assert "test1" in genotypes
    assert len(rows) == 2

def test_gff_parser_arguments():
    """Test GFF parser argument parsing."""
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

def test_fasta_parser_arguments():
    """Test FASTA parser argument parsing."""
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

def test_target_genes_with_multiple_genotypes(setup_test_files):
    """Test target genes parsing with multiple genotypes."""
    # Create a test target genes file with multiple genotypes
    target_content = """gene_name,gene_ID_genotype1,gene_ID_genotype2
GENE001,GENE001,GENE001
GENE002,GENE002,GENE002"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test parsing
    genotypes, rows = read_gff_target_genes(target_file)
    
    assert len(genotypes) == 2
    assert "genotype1" in genotypes
    assert "genotype2" in genotypes
    assert len(rows) == 2

def test_target_genes_missing_gene_name(setup_test_files):
    """Test target genes parsing with missing gene_name column."""
    # Create a test target genes file without gene_name
    target_content = """gene_ID_test1
GENE001
GENE002"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test that it raises ValueError
    with pytest.raises(ValueError) as exc_info:
        read_gff_target_genes(target_file)
    
    assert "gene_name column" in str(exc_info.value)

def test_cli_argument_validation():
    """Test CLI argument validation."""
    # Test GFF parser
    gff_parser = setup_gff_parser()
    
    # Test missing required arguments
    with pytest.raises(SystemExit):
        gff_parser.parse_args(["--pangenome-folder", "/test"])
    
    # Test FASTA parser
    fasta_parser = setup_fasta_parser()
    
    # Test missing required arguments
    with pytest.raises(SystemExit):
        fasta_parser.parse_args(["--pangenome-folder", "/test"])

def test_cli_help_output(capsys):
    """Test CLI help output."""
    # Test GFF help
    gff_parser = setup_gff_parser()
    with pytest.raises(SystemExit):
        gff_parser.parse_args(["--help"])
    
    # Test FASTA help
    fasta_parser = setup_fasta_parser()
    with pytest.raises(SystemExit):
        fasta_parser.parse_args(["--help"])

def test_target_genes_list_format(setup_test_files):
    """Test target genes parsing with list format."""
    # Create a test target genes file with list format
    target_content = """gene_name,gene_ID_test1
GENE001,['GENE001','GENE002']
GENE002,['GENE003']"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test parsing
    genotypes, rows = read_gff_target_genes(target_file)
    
    assert len(genotypes) == 1
    assert len(rows) == 2
    
    # Check that list format is handled correctly
    first_row = rows[0]
    gene_id_raw = first_row.get("gene_ID_test1", "")
    assert gene_id_raw.startswith("['")

def test_target_genes_empty_values(setup_test_files):
    """Test target genes parsing with empty values."""
    # Create a test target genes file with empty values
    target_content = """gene_name,gene_ID_test1,gene_ID_test2
GENE001,GENE001,
GENE002,,GENE002"""
    
    target_file = setup_test_files / "target_genes.csv"
    target_file.write_text(target_content)
    
    # Test parsing
    genotypes, rows = read_gff_target_genes(target_file)
    
    assert len(genotypes) == 2
    assert len(rows) == 2

def test_cli_argument_types():
    """Test CLI argument type handling."""
    # Test GFF parser
    gff_parser = setup_gff_parser()
    
    # Test with valid arguments
    args = gff_parser.parse_args([
        "--pangenome-folder", "/test/folder",
        "--pangenome-index", "/test/index.csv",
        "--target-genes", "/test/genes.csv",
        "--output", "/test/output.gff"
    ])
    
    assert isinstance(args.pangenome_folder, str)
    assert isinstance(args.pangenome_index, str)
    assert isinstance(args.target_genes, str)
    assert isinstance(args.output, str)
    
    # Test FASTA parser with numeric arguments
    fasta_parser = setup_fasta_parser()
    
    args = fasta_parser.parse_args([
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
    
    assert isinstance(args.upstream, int)
    assert isinstance(args.downstream, int)
    assert isinstance(args.inner_start, int)
    assert isinstance(args.inner_end, int)
    assert isinstance(args.pad, int)