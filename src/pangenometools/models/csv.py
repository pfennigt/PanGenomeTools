"""
CSV module for PanGenomeTools.

This module provides functionality for parsing CSV files and converting them to JSON format.
"""

import json
from pathlib import Path
import pandas as pd
from typing import Any


def parse_csv_columns(csv_path: Path, output_dir: Path|str|None = None, groupby: str|None=None) -> list[str]:
    """
    Parse a CSV file and return each column as a separate list.
    
    Args:
        csv_path: Path to the CSV file
        output_dir: Optional directory to save individual column JSON files
        
    Returns:
        Dictionary mapping column names to lists of values
    """
    # Read CSV using pandas
    df = pd.read_csv(csv_path)

    # Make the output directory
    if output_dir is None:
        output_dir = Path(".")
    else:
        output_dir=Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)


    # If a grouping column is given, group by that column
    if groupby is not None: 
        groups = {}
        for name, values in df.groupby(groupby):
            values=values.drop(groupby, axis=1)
            groups[name] = {k:v.to_numpy().flatten() for k,v in values.T.iterrows()}

        df = pd.DataFrame(groups).T
        df[groupby] = df.index

    # Save individual columns as files
    for column_name, values in df.T.iterrows():
            output_file = output_dir / f"{column_name}.json"
            values.to_json(output_file, mode="w", index=False, orient="records")
    
    return df.columns.to_list()