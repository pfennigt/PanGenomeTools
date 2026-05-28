"""
Base class for pangenome file handlers.

This module provides the foundation for all file type handlers in PanGenomeTools.
"""

from pathlib import Path
from typing import Dict, Any
import logging
import csv


class PangenomeFileHandler:
    """Base class for handling pangenome files."""

    def __init__(self, pangenome_folder: Path, pangenome_index: Path):
        """
        Initialize the pangenome file handler.

        Args:
            pangenome_folder: Path to the pangenome folder
            pangenome_index: Path to the pangenome index CSV file
        """
        self.pangenome_folder = pangenome_folder
        self.pangenome_index = self._load_pangenome_index(pangenome_index)
        self.logger = self._setup_logger()

    def _load_pangenome_index(self, index_path: Path) -> Dict[str, Dict[str, str]]:
        """
        Load and validate pangenome index.

        Args:
            index_path: Path to the pangenome index CSV file

        Returns:
            Dictionary mapping genotype names to file paths

        Raises:
            ValueError: If index file is missing required columns
        """
        data = {}
        with open(index_path, newline="") as fh:
            reader = csv.DictReader(fh)
            required = {"genotype", "annotation", "assembly"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"pangenome_index.csv must contain columns {required}")
            for r in reader:
                data[r["genotype"]] = {
                    "annotation": r["annotation"],
                    "assembly": r["assembly"],
                }
        return data

    def _setup_logger(self) -> logging.Logger:
        """
        Setup logger with standard configuration.

        Returns:
            Configured logger instance
        """
        logger = logging.getLogger(self.__class__.__name__)
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger

    def resolve_path(self, genotype: str, file_type: str) -> Path:
        """
        Resolve file path for a given genotype and file type.

        Args:
            genotype: Genotype identifier
            file_type: Type of file ('annotation' or 'assembly')

        Returns:
            Resolved file path

        Raises:
            ValueError: If genotype not found in index
            FileNotFoundError: If resolved file doesn't exist
        """
        if genotype not in self.pangenome_index:
            raise ValueError(f"Genotype {genotype} not found in pangenome index")

        file_path = self.pangenome_index[genotype][file_type]
        full_path = self.pangenome_folder / file_path

        if not full_path.exists():
            raise FileNotFoundError(f"File not found: {full_path}")

        return full_path