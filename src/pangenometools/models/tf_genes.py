"""
TF gene extraction model for PanGenomeTools.

This module provides functionality for extracting TF genes from Excel files
and organizing them by transcription factor family.
"""

from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import logging


class TFGeneExtractor:
    """Handler for extracting TF genes from Excel files."""

    def __init__(self, excel_path: Path, tf_family_column: str = "TF family", 
                 gene_name_column: str = "TF Gene", sheet_name: Optional[str] = None):
        """
        Initialize TF gene extractor.

        Args:
            excel_path: Path to Excel file containing TF gene list
            tf_family_column: Name of column containing TF family names
            gene_name_column: Name of column containing gene names
            sheet_name: Name of Excel sheet to read (None for first sheet)
        """
        self.excel_path = Path(excel_path)
        self.tf_family_column = tf_family_column
        self.gene_name_column = gene_name_column
        self.sheet_name = sheet_name
        self.logger = logging.getLogger(self.__class__.__name__)
        
        # Validate file exists
        if not self.excel_path.exists():
            raise FileNotFoundError(f"Excel file not found: {self.excel_path}")

    def read_excel_file(self) -> pd.DataFrame:
        """
        Read Excel file and return as DataFrame.

        Returns:
            DataFrame containing TF gene information

        Raises:
            ValueError: If required columns are missing
        """
        try:
            # Read Excel file
            if self.sheet_name:
                df = pd.read_excel(self.excel_path, sheet_name=self.sheet_name)
            else:
                df = pd.read_excel(self.excel_path)
            
            self.logger.info(f"Read {len(df)} rows from {self.excel_path}")
            
            # Validate required columns
            required_columns = [self.tf_family_column, self.gene_name_column]
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                raise ValueError(f"Missing required columns: {missing_columns}")
            
            return df
            
        except Exception as e:
            self.logger.error(f"Error reading Excel file: {e}")
            raise

    def get_tf_families(self, df: pd.DataFrame) -> List[str]:
        """
        Get list of unique TF families in the data.

        Args:
            df: DataFrame containing TF gene information

        Returns:
            List of unique TF family names
        """
        families = df[self.tf_family_column].dropna().unique().tolist()
        self.logger.info(f"Found {len(families)} unique TF families")
        return families

    def extract_genes_by_family(self, df: pd.DataFrame, 
                                families: Optional[List[str]] = None) -> Dict[str, List[str]]:
        """
        Extract gene names grouped by TF family.

        Args:
            df: DataFrame containing TF gene information
            families: Optional list of specific families to extract (None for all)

        Returns:
            Dictionary mapping family names to lists of gene names
        """
        # Filter by families if specified
        if families:
            df = df[df[self.tf_family_column].isin(families)]
        
        # Group by family and extract gene names
        family_genes = {}
        for family, group in df.groupby(self.tf_family_column):
            genes = group[self.gene_name_column].dropna().unique().tolist()
            family_genes[family] = genes
            self.logger.info(f"Family '{family}': {len(genes)} genes")
        
        return family_genes

    def save_family_csv(self, family: str, genes: List[str], output_dir: Path) -> Path:
        """
        Save gene list for a single family to CSV file.

        Args:
            family: TF family name
            genes: List of gene names
            output_dir: Directory to save CSV file

        Returns:
            Path to saved CSV file
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Sanitize filename
        safe_family_name = "".join(c for c in family if c.isalnum() or c in ('_', '-'))
        output_path = output_dir / f"{safe_family_name}.csv"
        
        # Create DataFrame and save
        df = pd.DataFrame({"gene_name": genes})
        df.to_csv(output_path, index=False)
        
        self.logger.info(f"Saved {len(genes)} genes for family '{family}' to {output_path}")
        return output_path

    def extract_and_save_all(self, output_dir: Path, 
                            families: Optional[List[str]] = None) -> Dict[str, Path]:
        """
        Extract genes by family and save to CSV files.

        Args:
            output_dir: Directory to save CSV files
            families: Optional list of specific families to extract (None for all)

        Returns:
            Dictionary mapping family names to output file paths
        """
        # Read Excel file
        df = self.read_excel_file()
        
        # Extract genes by family
        family_genes = self.extract_genes_by_family(df, families)
        
        # Save each family to CSV
        output_files = {}
        for family, genes in family_genes.items():
            output_path = self.save_family_csv(family, genes, output_dir)
            output_files[family] = output_path
        
        return output_files