"""
Utility functions for PanGenomeTools.

This module provides shared functionality used across different handlers.
"""

import re
import pandas as pd
import numpy as np
import sys
import contextlib
from tqdm.contrib import DummyTqdmFile
import csv
from pathlib import Path

from functools import partial

from typing import Union, Tuple, Dict, List, Optional

import logging

# Set up a logger
_formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
def setup_logger(name, log_file, level=logging.INFO, formatter=_formatter):
    """Get a preexisting logger or create it if it doesn't exist yet"""

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if len(logger.handlers) == 0:
        handler = logging.FileHandler(log_file)        
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger

# Function to extract information from FASTA-like headers
def extract_fasta_header_like(header:Union[str,list,pd.Series], key_prefix:str="", auto_key_prefix="_", chunk_suffix:Union[bool,None]=None)->pd.Series:
    """Function to extract information in a FASTA-like header as produced by the "extract_genes" functions. If a list-like of headers is given, it assumes the same structure for ever header.

    Args:
        header (str): Header of a sequence
        key_prefix (str, optional): String to be added to the beginning of the output indices. Defaults to "".
        chunk_suffix (Union[bool,None], optional): Does the header contain a suffix indicating the chunk position produced by the DeepCistrome function 'split_sequence_chunks_to_250bp'. If None, the suffix is automatically detected. Defaults to None.

    Raises:
        ValueError: If the header structure dos not follow the convention '[>]<_id> [<key>=<value>] ... [<key>=<value>][_<chunk_suffix>]'.

    Returns:
        pd.Series: A Series object storing the header information as key:value pairs.
    """

    # Turn the header into a pd.Series
    if isinstance(header, str):
        header = pd.DataFrame([header])
    elif isinstance(header, list):
        header = pd.DataFrame(header)
    elif isinstance(header, pd.Series):
        header = pd.DataFrame(header)
    else:
        raise ValueError(f"Unknown input type for header: {type(header)}")

    header.columns = ["_hdr"]

    # Remove flanking whitespaces and the initial ">"
    header.loc[:,"_hdr"] = header.loc[:,"_hdr"].str.strip()
    header.loc[:,"_hdr"] = header.loc[:,"_hdr"].str.removeprefix(">")

    if chunk_suffix is None:
        chunk_suffix = re.search("_chunk[0-9_\\-]+$", header.iloc[0,0])

    # If there is a chunk suffix, extract it
    if chunk_suffix:
        # Split off the chunk information 
        _header_split = header.loc[:,"_hdr"].str.split("_chunk").apply(pd.Series)

        # Set the header without the chunk information as the new header
        header = _header_split.loc[:,[0]]
        header.columns = ["_hdr"]

        # Get the chunk range information
        chunk_range = _header_split.loc[:,1].str.split("_").apply(pd.Series).iloc[:,1].str.split("-").apply(pd.Series)
        chunk_range.columns = [f"{key_prefix}{auto_key_prefix}chunk_start", f"{key_prefix}{auto_key_prefix}chunk_end"]
        chunk_range = chunk_range.astype(int)

    # Split the header by certain spaces
    # Replace first occurrence of space
    sub_first = partial(re.sub, "^([^ ]+) ", "\\1___")
    header.loc[:,"_hdr"] = header.loc[:,"_hdr"].apply(sub_first)

    # Replace space before keys
    sub_pairs = partial(re.sub, " ([^ =]+=)", "___\\1")
    header.loc[:,"_hdr"] = header.loc[:,"_hdr"].apply(sub_pairs)

    # Split by the introduced spacer to extract the information fields
    info = header.loc[:,"_hdr"].str.split("___").apply(pd.Series)

    # Create a container for the column names
    col_names = np.full(info.shape[1], "", dtype="U64")

    # Set the first column name as id
    col_names[0] = f"{auto_key_prefix}id"

    # Determine the column names form the first row
    is_keypair = [False] + (info.iloc[0,1:].str.find("=") > 0).to_list()

    # Extract the key names from the keypair columns
    col_names[is_keypair] = info.iloc[0, is_keypair].str.split("=").apply(pd.Series).iloc[:,0].to_list()

    # Set general names for the remaining columns
    is_remaining = np.invert(is_keypair[1:])
    if np.any(is_remaining):
        col_names[1:][is_remaining] = [f"{auto_key_prefix}descriptor{i if i>0 else ''}" for i in range(np.sum(is_remaining))]
    col_names = [key_prefix + x for x in col_names]
    info.columns = col_names

    # Remove the prefixes
    for _name in col_names[1:]:
        info.loc[:, _name] = info.loc[:, _name].str.removeprefix(f"{_name.removeprefix(key_prefix)}=")

    # Add the chunk information
    if chunk_suffix:
        info = pd.concat([info, chunk_range], axis=1)

    # Get meta info about the extraction and chunking
    # Automatically detect the positions of the sequences
    extract_info = { 
        
    }
    _extract_type = info.loc[0,f"{key_prefix}type"]
    _extract_region = info.loc[0,f"{key_prefix}{auto_key_prefix}id"].split("_")[-1]
    _extract_options = info.loc[0,f"{key_prefix}extraction_options"]
    _extract_options = {a:(int(b) if a in ["upstream", "inner_start", "inner_end", "downstream"] else b) for a,b in [x.split(":") for x in _extract_options.split("&")]}

    # Collect the information
    extract_info = {
        "type": _extract_type.split("_")[0],
        "region": _extract_region,
        "options": _extract_options
    }

    # Determine if there is a gap
    if len(_extract_options["pad"]) > 0 and not "whole-seq" in _extract_options:
        # Get the gap padding and the position
        pad_len = len(_extract_options["pad"])
        gap_info = {
            "gap_start": _extract_options["upstream"] + _extract_options["inner_start"],
            "gap_pad": pad_len,
            "gap_pos": _extract_options["upstream"] + _extract_options["inner_start"] + (pad_len // 2)
        }
    else:
        gap_info = None
        pad_len = 0

    # Determine the positions
    if _extract_region == "flanking":
        ref_pos = [
            _extract_options["upstream"],
            _extract_options["upstream"] + _extract_options["inner_start"] + pad_len + _extract_options["inner_end"]
        ]
    elif _extract_region == "promoter":
        ref_pos = [
            _extract_options["upstream"]
        ]
    elif _extract_region == "terminator":
        ref_pos = [
            _extract_options["inner_end"]
        ]

    # Set the labels
    if _extract_type.startswith("CDS"):
        ref_label = np.array(["AUG", "Ter", "AUG/Ter"])
    elif _extract_type.startswith("gene"):
        ref_label = np.array(["TSS", "TTS", "TSS/TTS"])
    else:
        ref_label = np.array(["<start>", "<end>", "<start>/<end>"])

    if len(ref_pos) == 1:
        ref_label = ref_label[[0]]

    # Get the info about the chunk size and their distance
    if chunk_suffix:
        chunk_info = {
            "chunk_length": info.loc[0, "chunk_end"] - info.loc[0, "chunk_start"] + 1,
            "window_stride": info.loc[1, "chunk_start"] - info.loc[0, "chunk_start"]
        }
    else:
        chunk_info = None

    # Collect all meta information
    meta_info = {
        "extraction": extract_info,
        "reference_points":{
            "positions": ref_pos,
            "labels": ref_label
        },
        "gap": gap_info,
        "chunks": chunk_info,
    }

    return info, meta_info

def get_genotype_alias_map(file_path):
    # Load the alias map to use common identifiers for the genotypes
    _genotype_alias_df = pd.read_csv(file_path).fillna("")

    # Add the name of the genotype to the alias list
    _genotype_alias_df["alias"] = (
        _genotype_alias_df["id"] + ";" + _genotype_alias_df["alias"]
    )

    # Set all aliases to uppercase and split the various aliases
    _genotype_alias_df["alias"] = _genotype_alias_df["alias"].str.upper().str.split(";")

    # Create a mapping to the common names
    genotype_alias_map = {}

    for _, (name, aliases) in _genotype_alias_df.iterrows():
        for alias in aliases:
            if len(alias) > 0:
                genotype_alias_map[alias] = name

    return genotype_alias_map

# Define a rerouter for stdout for tqdm
@contextlib.contextmanager
def std_out_err_redirect_tqdm():
    orig_out_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = map(DummyTqdmFile, orig_out_err)
        yield orig_out_err[0]
    # Relay exceptions
    except Exception as exc:
        raise exc
    # Always restore sys.stdout/err if necessary
    finally:
        sys.stdout, sys.stderr = orig_out_err

# Reader for the pangenome index
def read_index(index_path: Path) -> Dict[str, Dict[str, str]]:
    d = {}
    with open(index_path, newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"genotype", "annotation", "assembly"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"pangenome_index.csv must contain columns {required}")
        for r in reader:
            d[r["genotype"]] = {
                "annotation": r["annotation"],
                "assembly": r["assembly"],
            }
    return d

# Shared functions for coordinate calculations
def calculate_coordinate_boundaries(
    start: int,
    end: int,
    strand: str,
    upstream: int = 0,
    downstream: int = 0,
    inner_start: int = 0,
    inner_end: int = 0,
    whole_seq: bool = False,
    use_five_prime_direction: bool = False
) -> Tuple[int, int, int, int]:
    """
    Calculate coordinate boundaries for sequence extraction.
    
    Args:
        start: Feature start coordinate
        end: Feature end coordinate  
        strand: Feature strand (+ or -)
        upstream: Nucleotides to include upstream (5' direction)
        downstream: Nucleotides to include downstream (3' direction)
        inner_start: Nucleotides to include from start of feature
        inner_end: Nucleotides to include from end of feature
        whole_seq: If True, extract entire feature sequence
        use_five_prime_direction: If True, always interpret upstream/downstream
                                   in 5' direction regardless of strand.
                                   If False, reverse upstream/downstream for negative strand.
    
    Returns:
        Tuple of (left_a, left_b, right_a, right_b) coordinates
    """
    # Interior coordinate logic
    bstart, bend = (end, start) if strand == "-" else (start, end)

    if whole_seq:
        left_a, left_b = min(start, end), max(start, end)
        right_a, right_b = 1, 0
    else:
        # Determine direction based on strand and use_five_prime_direction
        if use_five_prime_direction:
            # Always use 5' direction (upstream = 5', downstream = 3')
            (upstream, inner_start, inner_end, downstream) = (upstream, inner_start, inner_end, downstream)
        elif strand == "-":
            # For negative strand without use_five_prime_direction, reverse the direction
            (upstream, inner_start, inner_end, downstream) = (downstream, inner_end, inner_start, upstream)
        else:
            # For positive strand, use as-is
            (upstream, inner_start, inner_end, downstream) = (upstream, inner_start, inner_end, downstream)

        left_a = bstart - upstream
        left_b = bstart + inner_start - 1
        
        right_a = bend - inner_end + 1
        right_b = bend + downstream

    return left_a, left_b, right_a, right_b

def clip_coordinates(
    left_a: int,
    left_b: int,
    right_a: int,
    right_b: int,
    chrom_len: int
) -> Tuple[int, int, int, int]:
    """
    Clip coordinates to valid chromosome range.
    
    Args:
        left_a: Left segment start coordinate
        left_b: Left segment end coordinate
        right_a: Right segment start coordinate
        right_b: Right segment end coordinate
        chrom_len: Chromosome length
    
    Returns:
        Tuple of clipped coordinates (ll, lh, rl, rh)
    """
    def clip_coord(coord: int, chrom_len: int) -> int:
        """Clip coordinate to valid range."""
        return max(1, min(coord, chrom_len))

    ll, lh = clip_coord(left_a, chrom_len), clip_coord(left_b, chrom_len)
    rl, rh = clip_coord(right_a, chrom_len), clip_coord(right_b, chrom_len)

    return ll, lh, rl, rh

def read_target_genes(target_path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """
    Read target genes from CSV file.
    
    Args:
        target_path: Path to target genes CSV file
    
    Returns:
        Tuple of (genotypes, rows)
    
    Raises:
        ValueError: If target_genes.csv doesn't contain required columns
    """
    rows = []
    with open(target_path, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        geno_cols = [h for h in headers if h.startswith("gene_ID_")]
        if "gene_name" not in headers:
            raise ValueError("target_genes.csv must contain a gene_name column")
        for r in reader:
            rows.append(r)

    genotypes = [h.replace("gene_ID_", "") for h in geno_cols]
    return genotypes, rows