"""
Homolog parser model for PanGenomeTools.

This module provides functionality for parsing and summarizing homolog search results.
"""

from pathlib import Path
from typing import List, Optional, Tuple, Dict
import logging
import pandas as pd
from Bio import SeqIO


class HomologParser:
    """Handler for parsing and summarizing homolog search results."""

    def __init__(self, blast_results: Path, query_fasta: Path):
        """
        Initialize homolog parser.

        Args:
            blast_results: Path to BLAST/DIAMOND results file
            query_fasta: Path to original query FASTA file
        """
        self.blast_results = Path(blast_results)
        self.query_fasta = Path(query_fasta)
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Validate files exist
        if not self.blast_results.exists():
            raise FileNotFoundError(f"BLAST results file not found: {self.blast_results}")
        if not self.query_fasta.exists():
            raise FileNotFoundError(f"Query FASTA file not found: {self.query_fasta}")

    def read_blast_results(self) -> pd.DataFrame:
        """
        Read BLAST/DIAMOND results file.

        Returns:
            DataFrame containing BLAST results
        """
        # Define column names for BLAST tabular format
        columns = [
            'qseqid', 'sseqid', 'pident', 'length',
            'mismatch', 'gapopen', 'qstart', 'qend',
            'sstart', 'send', 'evalue', 'bitscore'
        ]
        
        df = pd.read_csv(self.blast_results, sep='\t', header=None, names=columns)
        self.logger.info(f"Read {len(df)} BLAST hits from {self.blast_results}")
        
        return df

    def get_query_lengths(self) -> Dict[str, int]:
        """
        Get lengths of query sequences.

        Returns:
            Dictionary mapping query IDs to sequence lengths
        """
        query_lengths = {}
        for record in SeqIO.parse(self.query_fasta, "fasta"):
            query_lengths[record.id] = len(record.seq)
        
        self.logger.info(f"Retrieved lengths for {len(query_lengths)} query sequences")
        return query_lengths

    def calculate_coverage(self, df: pd.DataFrame, query_lengths: Dict[str, int]) -> pd.DataFrame:
        """
        Calculate query coverage for each hit.

        Args:
            df: DataFrame containing BLAST results
            query_lengths: Dictionary mapping query IDs to sequence lengths

        Returns:
            DataFrame with added coverage column
        """
        def calc_coverage(row):
            query_id = row['qseqid']
            if query_id in query_lengths:
                query_len = query_lengths[query_id]
                alignment_length = abs(row['qend'] - row['qstart']) + 1
                return (alignment_length / query_len) * 100
            return 0.0
        
        df['qcovs'] = df.apply(calc_coverage, axis=1)
        return df

    def filter_results(self, df: pd.DataFrame, min_identity: float = 30.0,
                      min_coverage: float = 50.0) -> pd.DataFrame:
        """
        Filter results by identity and coverage thresholds.

        Args:
            df: DataFrame containing BLAST results
            min_identity: Minimum percent identity threshold
            min_coverage: Minimum query coverage threshold

        Returns:
            Filtered DataFrame
        """
        initial_count = len(df)
        
        # Filter by identity
        df_filtered = df[df['pident'] >= min_identity].copy()
        
        # Filter by coverage
        df_filtered = df_filtered[df_filtered['qcovs'] >= min_coverage].copy()
        
        filtered_count = len(df_filtered)
        self.logger.info(f"Filtered {initial_count} hits to {filtered_count} hits "
                          f"(identity >= {min_identity}%, coverage >= {min_coverage}%)")
        
        return df_filtered

    def get_best_hits(self, df: pd.DataFrame, top_hits_only: bool = True) -> pd.DataFrame:
        """
        Get best hits for each query sequence.

        Args:
            df: DataFrame containing BLAST results
            top_hits_only: If True, only keep the best hit per query

        Returns:
            DataFrame with best hits
        """
        if top_hits_only:
            # Sort by bitscore (descending) and evalue (ascending)
            df_sorted = df.sort_values(['qseqid', 'bitscore', 'evalue'],
                                      ascending=[True, False, True])
            
            # Keep only the first hit per query
            best_hits = df_sorted.groupby('qseqid').first().reset_index()
            
            self.logger.info(f"Selected {len(best_hits)} best hits (one per query)")
            return best_hits
        else:
            return df

    def create_summary_table(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Create a summary table with key information.

        Args:
            df: DataFrame containing BLAST results

        Returns:
            Summary DataFrame
        """
        summary = df[[
            'qseqid', 'sseqid', 'pident', 'qcovs',
            'length', 'evalue', 'bitscore'
        ]].copy()
        
        # Rename columns for clarity
        summary.columns = [
            'query_id', 'subject_id', 'percent_identity',
            'query_coverage', 'alignment_length', 'evalue', 'bitscore'
        ]
        
        return summary

    def create_detailed_table(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Create a detailed table with all information.

        Args:
            df: DataFrame containing BLAST results

        Returns:
            Detailed DataFrame
        """
        detailed = df.copy()
        
        # Rename columns for clarity
        detailed.columns = [
            'query_id', 'subject_id', 'percent_identity',
            'alignment_length', 'mismatches', 'gap_opens',
            'query_start', 'query_end', 'subject_start',
            'subject_end', 'evalue', 'bitscore', 'query_coverage'
        ]
        
        return detailed

    def calculate_statistics(self, df: pd.DataFrame) -> Dict[str, any]:
        """
        Calculate summary statistics.

        Args:
            df: DataFrame containing BLAST results

        Returns:
            Dictionary with statistics
        """
        stats = {
            'total_hits': len(df),
            'unique_queries': df['qseqid'].nunique(),
            'unique_subjects': df['sseqid'].nunique(),
            'mean_identity': df['pident'].mean(),
            'median_identity': df['pident'].median(),
            'mean_coverage': df['qcovs'].mean(),
            'median_coverage': df['qcovs'].median(),
            'mean_evalue': df['evalue'].mean(),
            'median_evalue': df['evalue'].median()
        }
        
        return stats

    def save_results(self, output_dir: Path, min_identity: float = 30.0,
                    min_coverage: float = 50.0, top_hits_only: bool = True) -> Tuple[Path, Path, Path]:
        """
        Parse and save homolog search results.

        Args:
            output_dir: Directory to save output files
            min_identity: Minimum percent identity threshold
            min_coverage: Minimum query coverage threshold
            top_hits_only: If True, only keep the best hit per query

        Returns:
            Tuple of (summary_file, detailed_file, stats_file) paths
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Read results
        df = self.read_blast_results()

        # Get query lengths and calculate coverage
        query_lengths = self.get_query_lengths()
        df = self.calculate_coverage(df, query_lengths)

        # Filter results
        df_filtered = self.filter_results(df, min_identity, min_coverage)

        # Get best hits
        df_best = self.get_best_hits(df_filtered, top_hits_only)

        # Create summary and detailed tables
        summary_df = self.create_summary_table(df_best)
        detailed_df = self.create_detailed_table(df_filtered)

        # Calculate statistics
        stats = self.calculate_statistics(df_filtered)

        # Save files
        summary_file = output_dir / "homologs_summary.csv"
        detailed_file = output_dir / "homologs_detailed.csv"
        stats_file = output_dir / "statistics.txt"

        summary_df.to_csv(summary_file, index=False)
        detailed_df.to_csv(detailed_file, index=False)

        # Write statistics file
        with open(stats_file, 'w') as f:
            f.write("Homolog Search Statistics\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total hits: {stats['total_hits']}\n")
            f.write(f"Unique queries: {stats['unique_queries']}\n")
            f.write(f"Unique subjects: {stats['unique_subjects']}\n\n")
            f.write("Identity Statistics:\n")
            f.write(f"  Mean: {stats['mean_identity']:.2f}%\n")
            f.write(f"  Median: {stats['median_identity']:.2f}%\n\n")
            f.write("Coverage Statistics:\n")
            f.write(f"  Mean: {stats['mean_coverage']:.2f}%\n")
            f.write(f"  Median: {stats['median_coverage']:.2f}%\n\n")
            f.write("E-value Statistics:\n")
            f.write(f"  Mean: {stats['mean_evalue']:.2e}\n")
            f.write(f"  Median: {stats['median_evalue']:.2e}\n")

        self.logger.info(f"Saved results to {output_dir}")
        self.logger.info(f"  Summary: {summary_file}")
        self.logger.info(f"  Detailed: {detailed_file}")
        self.logger.info(f"  Statistics: {stats_file}")

        return summary_file, detailed_file, stats_file