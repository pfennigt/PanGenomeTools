"""
Sequence retriever model for PanGenomeTools.

This module provides functionality for retrieving protein sequences from reference genomes
using AGAT tools with GFF subsetting.
"""

from pathlib import Path
from typing import List, Optional, Tuple
import logging
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pyranges as pr


class SequenceRetriever:
    """Handler for retrieving protein sequences from reference genomes."""

    def __init__(self, reference_genome: Path, reference_gff: Path):
        """
        Initialize sequence retriever.

        Args:
            reference_genome: Path to reference genome FASTA file
            reference_gff: Path to reference GFF annotation file
        """
        self.reference_genome = Path(reference_genome)
        self.reference_gff = Path(reference_gff)
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Validate files exist
        if not self.reference_genome.exists():
            raise FileNotFoundError(f"Reference genome not found: {self.reference_genome}")
        if not self.reference_gff.exists():
            raise FileNotFoundError(f"Reference GFF not found: {self.reference_gff}")

    def check_agat_available(self) -> bool:
        """
        Check if AGAT is installed and available.

        Returns:
            True if AGAT is available, False otherwise
        """
        try:
            subprocess.run(["agat", "--version"], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def subset_gff_by_genes(self, gene_list: List[str], output_gff: Path) -> int:
        """
        Subset GFF file to include only specified genes and their features.

        Args:
            gene_list: List of gene names to extract
            output_gff: Path to output GFF file

        Returns:
            Number of features extracted
        """
        self.logger.info(f"Subsetting GFF for {len(gene_list)} genes")
        
        # Read GFF file using pyranges
        try:
            gff = pr.read_gff3(str(self.reference_gff))
        except Exception:
            # Try reading as GFF if GFF3 fails
            gff = pr.read_gff(str(self.reference_gff))
        
        # Convert to DataFrame for easier manipulation
        df = gff.df
        
        # Find all gene IDs that match the gene names
        # Gene names might be in ID, Name, or gene_name attributes
        gene_ids = set()
        
        for gene_name in gene_list:
            # Search in ID column
            matches = df[df['ID'].str.contains(gene_name, case=False, na=False, regex=False)]
            if len(matches) > 0:
                gene_ids.update(matches['ID'].tolist())
            
            # Search in Name attribute if present
            if 'Name' in df.columns:
                matches = df[df['Name'].str.contains(gene_name, case=False, na=False, regex=False)]
                if len(matches) > 0:
                    gene_ids.update(matches['ID'].tolist())
            
            # Search in gene_name attribute if present
            if 'gene_name' in df.columns:
                matches = df[df['gene_name'].str.contains(gene_name, case=False, na=False, regex=False)]
                if len(matches) > 0:
                    gene_ids.update(matches['ID'].tolist())
        
        self.logger.info(f"Found {len(gene_ids)} matching gene IDs")
        
        if len(gene_ids) == 0:
            self.logger.warning("No matching genes found in GFF file")
            return 0
        
        # Get all features that are children of these genes
        # This includes mRNA, exon, CDS, etc.
        all_feature_ids = set(gene_ids)
        
        # Iteratively find all children
        max_iterations = 10  # Prevent infinite loops
        iteration = 0
        
        while iteration < max_iterations:
            iteration += 1
            # Find features whose Parent is in our current set
            parent_matches = df[df['Parent'].isin(all_feature_ids)]
            new_ids = set(parent_matches['ID'].tolist())
            
            if new_ids.issubset(all_feature_ids):
                # No new IDs found, we're done
                break
            
            all_feature_ids.update(new_ids)
        
        self.logger.info(f"Found {len(all_feature_ids)} total features (including children)")
        
        # Subset the GFF
        subset_df = df[df['ID'].isin(all_feature_ids)]
        
        # Write to output GFF file
        # Convert back to PyRanges and write
        subset_gff = pr.PyRanges(subset_df)
        subset_gff.to_gff3(str(output_gff))
        
        self.logger.info(f"Wrote {len(subset_df)} features to {output_gff}")
        
        return len(subset_df)

    def retrieve_sequences_with_agat(self, subset_gff: Path, output_fasta: Path,
                                     feature_type: str = "CDS") -> Tuple[int, int]:
        """
        Retrieve protein sequences using AGAT tools on a subset GFF.

        Args:
            subset_gff: Path to subset GFF file
            output_fasta: Path to output FASTA file
            feature_type: Type of feature to extract (CDS, mRNA, etc.)

        Returns:
            Tuple of (number of genes, number of protein sequences retrieved)

        Raises:
            RuntimeError: If AGAT is not installed or fails
        """
        if not self.check_agat_available():
            raise RuntimeError("AGAT is not installed or not in PATH. "
                             "Please install AGAT: conda install -c bioconda agat")

        # Use AGAT to extract and translate sequences
        cmd = [
            "agat_sp_extract_sequences.pl",
            "--gff", str(subset_gff),
            "--fasta", str(self.reference_genome),
            "--output", str(output_fasta),
            "--type", feature_type,
            "--aa"  # Extract CDS and translate to protein
        ]

        self.logger.info(f"Running AGAT command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Count sequences retrieved
            sequences_retrieved = 0
            if output_fasta.exists():
                sequences_retrieved = len(list(SeqIO.parse(output_fasta, "fasta")))
            
            # Count unique genes in the subset GFF
            subset_gff_data = pr.read_gff3(str(subset_gff))
            gene_features = subset_gff_data.df[subset_gff_data.df['Feature'] == 'gene']
            num_genes = len(gene_features['ID'].unique())
            
            self.logger.info(f"Retrieved {sequences_retrieved} protein sequences for {num_genes} genes")
            
            return num_genes, sequences_retrieved
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AGAT failed: {e.stderr}")
            raise RuntimeError(f"AGAT extraction failed: {e.stderr}")

    def retrieve_from_gene_list_csv(self, gene_list_csv: Path, output_fasta: Path,
                                    feature_type: str = "CDS",
                                    gene_id_column: str = "gene_name",
                                    keep_temp_files: bool = False) -> Tuple[int, int]:
        """
        Retrieve protein sequences from a CSV file containing gene names.

        Args:
            gene_list_csv: Path to CSV file with gene names
            output_fasta: Path to output FASTA file
            feature_type: Type of feature to extract (CDS, mRNA, etc.)
            gene_id_column: Column name in CSV containing gene IDs
            keep_temp_files: Whether to keep temporary GFF files

        Returns:
            Tuple of (number of genes found, number of sequences retrieved)
        """
        # Read gene list from CSV
        df = pd.read_csv(gene_list_csv)
        
        if gene_id_column not in df.columns:
            raise ValueError(f"Column '{gene_id_column}' not found in CSV file")

        gene_list = df[gene_id_column].dropna().unique().tolist()
        self.logger.info(f"Read {len(gene_list)} unique gene names from {gene_list_csv}")

        # Create temporary GFF file
        temp_gff = output_fasta.parent / f"{output_fasta.stem}_temp.gff3"
        
        try:
            # Subset GFF to include only specified genes
            num_features = self.subset_gff_by_genes(gene_list, temp_gff)
            
            if num_features == 0:
                self.logger.warning("No features found in subset GFF, skipping sequence extraction")
                return 0, 0
            
            # Retrieve sequences using AGAT
            num_genes, num_sequences = self.retrieve_sequences_with_agat(temp_gff, output_fasta, feature_type)
            
            return num_genes, num_sequences
            
        finally:
            # Clean up temporary GFF file unless requested to keep
            if not keep_temp_files and temp_gff.exists():
                temp_gff.unlink()
                self.logger.debug(f"Removed temporary GFF file: {temp_gff}")
            elif keep_temp_files and temp_gff.exists():
                self.logger.info(f"Kept temporary GFF file: {temp_gff}")