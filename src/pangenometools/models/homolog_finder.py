"""
Homolog finder model for PanGenomeTools.

This module provides functionality for finding homologs using BLAST or DIAMOND.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Dict
import logging
import subprocess
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


class HomologFinder:
    """Handler for finding homologs using BLAST or DIAMOND."""

    def __init__(self, method: str = "blast"):
        """
        Initialize homolog finder.

        Args:
            method: Search method - "blast" or "diamond"
        """
        self.method = method.lower()
        self.logger = logging.getLogger(self.__class__.__name__)
        
        if self.method not in ["blast", "diamond"]:
            raise ValueError(f"Unknown method: {method}. Use 'blast' or 'diamond'")

    def check_blast_available(self) -> bool:
        """
        Check if BLAST+ is installed and available.

        Returns:
            True if BLAST+ is available, False otherwise
        """
        try:
            subprocess.run(["blastp", "-version"], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def check_diamond_available(self) -> bool:
        """
        Check if DIAMOND is installed and available.

        Returns:
            True if DIAMOND is available, False otherwise
        """
        try:
            subprocess.run(["diamond", "version"], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def create_blast_database(self, target_fasta: Path, db_name: str) -> Path:
        """
        Create BLAST database from target sequences.

        Args:
            target_fasta: Path to target FASTA file
            db_name: Name for the BLAST database

        Returns:
            Path to created BLAST database

        Raises:
            RuntimeError: If BLAST is not installed or database creation fails
        """
        if not self.check_blast_available():
            raise RuntimeError("BLAST+ is not installed or not in PATH. "
                             "Please install BLAST+: conda install -c bioconda blast")

        cmd = [
            "makeblastdb",
            "-in", str(target_fasta),
            "-dbtype", "prot",
            "-out", db_name
        ]

        self.logger.info(f"Creating BLAST database: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"BLAST database created: {db_name}")
            return Path(db_name)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST database creation failed: {e.stderr}")
            raise RuntimeError(f"BLAST database creation failed: {e.stderr}")

    def create_diamond_database(self, target_fasta: Path, db_name: str) -> Path:
        """
        Create DIAMOND database from target sequences.

        Args:
            target_fasta: Path to target FASTA file
            db_name: Name for the DIAMOND database

        Returns:
            Path to created DIAMOND database

        Raises:
            RuntimeError: If DIAMOND is not installed or database creation fails
        """
        if not self.check_diamond_available():
            raise RuntimeError("DIAMOND is not installed or not in PATH. "
                             "Please install DIAMOND: conda install -c bioconda diamond")

        cmd = [
            "diamond", "makedb",
            "--in", str(target_fasta),
            "--db", db_name
        ]

        self.logger.info(f"Creating DIAMOND database: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"DIAMOND database created: {db_name}")
            return Path(db_name)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"DIAMOND database creation failed: {e.stderr}")
            raise RuntimeError(f"DIAMOND database creation failed: {e.stderr}")

    def run_blast_search(self, query_fasta: Path, db_path: Path, output_file: Path,
                        evalue: float = 1e-5, max_target_seqs: int = 10,
                        output_format: str = "6", blast_type: str = "blastp",
                        num_threads: int = 1) -> None:
        """
        Run BLAST search.

        Args:
            query_fasta: Path to query FASTA file
            db_path: Path to BLAST database
            output_file: Path to output file
            evalue: E-value threshold
            max_target_seqs: Maximum number of target sequences to report
            output_format: Output format (default: "6" for tabular)
            blast_type: BLAST type - "blastp", "blastx", "tblastn", or "tblastx"
            num_threads: Number of threads to use

        Raises:
            RuntimeError: If BLAST search fails
        """
        if not self.check_blast_available():
            raise RuntimeError("BLAST+ is not installed or not in PATH")

        # Define output format fields for tabular format
        format_fields = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
        
        cmd = [
            blast_type,
            "-query", str(query_fasta),
            "-db", str(db_path),
            "-out", str(output_file),
            "-evalue", str(evalue),
            "-max_target_seqs", str(max_target_seqs),
            "-outfmt", f"{output_format} {format_fields}",
            "-num_threads", str(num_threads)
        ]

        self.logger.info(f"Running BLAST search: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"BLAST search completed: {output_file}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLAST search failed: {e.stderr}")
            raise RuntimeError(f"BLAST search failed: {e.stderr}")

    def run_diamond_search(self, query_fasta: Path, db_path: Path, output_file: Path,
                          evalue: float = 1e-5, max_target_seqs: int = 10,
                          output_format: int = 6, sensitivity: str = "sensitive",
                          num_threads: int = 1) -> None:
        """
        Run DIAMOND search.

        Args:
            query_fasta: Path to query FASTA file
            db_path: Path to DIAMOND database
            output_file: Path to output file
            evalue: E-value threshold
            max_target_seqs: Maximum number of target sequences to report
            output_format: Output format (default: 6 for tabular)
            sensitivity: DIAMOND sensitivity level
            num_threads: Number of threads to use

        Raises:
            RuntimeError: If DIAMOND search fails
        """
        if not self.check_diamond_available():
            raise RuntimeError("DIAMOND is not installed or not in PATH")

        cmd = [
            "diamond", "blastp",
            "--query", str(query_fasta),
            "--db", str(db_path),
            "--out", str(output_file),
            "--evalue", str(evalue),
            "--max-target-seqs", str(max_target_seqs),
            "--outfmt", str(output_format),
            "--sensitive",  # Default sensitivity
            "--threads", str(num_threads)
        ]

        # Adjust sensitivity if specified
        if sensitivity != "sensitive":
            # Remove default --sensitive flag
            cmd = [c for c in cmd if c != "--sensitive"]
            
            # Add appropriate sensitivity flag
            sensitivity_flags = {
                "fast": "--fast",
                "mid-sensitive": "--mid-sensitive",
                "more-sensitive": "--more-sensitive",
                "very-sensitive": "--very-sensitive",
                "ultra-sensitive": "--ultra-sensitive"
            }
            
            if sensitivity in sensitivity_flags:
                cmd.append(sensitivity_flags[sensitivity])

        self.logger.info(f"Running DIAMOND search: {' '.join(cmd)}")
        
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"DIAMOND search completed: {output_file}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"DIAMOND search failed: {e.stderr}")
            raise RuntimeError(f"DIAMOND search failed: {e.stderr}")

    def find_homologs(self, query_fasta: Path, target_fasta: Path, output_file: Path,
                     evalue: float = 1e-5, max_target_seqs: int = 10,
                     output_format: str = "6", blast_type: str = "blastp",
                     diamond_sensitivity: str = "sensitive", num_threads: int = 1,
                     keep_database: bool = False) -> Dict[str, any]:
        """
        Find homologs using BLAST or DIAMOND.

        Args:
            query_fasta: Path to query FASTA file
            target_fasta: Path to target FASTA file
            output_file: Path to output file
            evalue: E-value threshold
            max_target_seqs: Maximum number of target sequences to report
            output_format: Output format
            blast_type: BLAST type (for BLAST method)
            diamond_sensitivity: DIAMOND sensitivity level (for DIAMOND method)
            num_threads: Number of threads to use
            keep_database: Whether to keep the database file

        Returns:
            Dictionary with search statistics

        Raises:
            RuntimeError: If search fails
        """
        # Count query sequences
        num_queries = len(list(SeqIO.parse(query_fasta, "fasta")))
        self.logger.info(f"Searching for homologs of {num_queries} query sequences")

        # Create database
        db_name = str(output_file.parent / output_file.stem)
        
        if self.method == "blast":
            db_path = self.create_blast_database(target_fasta, db_name)
        else:  # diamond
            db_path = self.create_diamond_database(target_fasta, db_name)

        try:
            # Run search
            if self.method == "blast":
                self.run_blast_search(query_fasta, db_path, output_file, evalue,
                                     max_target_seqs, output_format, blast_type, num_threads)
            else:  # diamond
                self.run_diamond_search(query_fasta, db_path, output_file, evalue,
                                      max_target_seqs, output_format, diamond_sensitivity, num_threads)

            # Count results
            if output_file.exists():
                results_df = pd.read_csv(output_file, sep='\t', header=None,
                                       names=['qseqid', 'sseqid', 'pident', 'length',
                                             'mismatch', 'gapopen', 'qstart', 'qend',
                                             'sstart', 'send', 'evalue', 'bitscore'])
                num_hits = len(results_df)
                num_queries_with_hits = results_df['qseqid'].nunique()
            else:
                num_hits = 0
                num_queries_with_hits = 0

            stats = {
                'method': self.method,
                'num_queries': num_queries,
                'num_hits': num_hits,
                'num_queries_with_hits': num_queries_with_hits,
                'database_path': str(db_path)
            }

            self.logger.info(f"Search completed: {num_hits} hits for {num_queries_with_hits} queries")

            return stats

        finally:
            # Clean up database unless requested to keep
            if not keep_database:
                db_extensions = [".phr", ".pin", ".psq", ".pdb", ".ptf", ".pto", ".dmnd"]
                for ext in db_extensions:
                    db_file = Path(str(db_path) + ext)
                    if db_file.exists():
                        db_file.unlink()
                        self.logger.debug(f"Removed database file: {db_file}")
            elif keep_database:
                self.logger.info(f"Kept database file: {db_path}")