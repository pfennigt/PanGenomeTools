#!/usr/bin/env python3

import re
import pandas as pd
import numpy as np

from functools import partial

from typing import Union

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
        chunk_suffix = re.search("_chunk[0-9_\-]+$", header.iloc[0,0])

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
    col_names[0] = f"{key_prefix}{auto_key_prefix}id"

    # Determine the column names form the first row
    is_keypair = [False] + (info.iloc[0,1:].str.find("=") > 0).to_list()

    # Extract the key names from the keypair columns
    col_names[is_keypair] = info.iloc[0, is_keypair].str.split("=").apply(pd.Series).iloc[:,0].to_list()

    # Set general names for the remaining columns
    is_remaining = np.invert(is_keypair[1:])
    if np.any(is_remaining):
        col_names[1:][is_remaining] = [f"{key_prefix}{auto_key_prefix}descriptor{i if i>0 else ''}" for i in range(np.sum(is_remaining))]
    info.columns = [key_prefix + x for x in col_names]

    # Remove the prefixes
    for _name in col_names[1:]:
        info.loc[:, _name] = info.loc[:, _name].str.removeprefix(f"{_name}=")

    # Add the chunk information
    if chunk_suffix:
        info = pd.concat([info, chunk_range], axis=1)
    return info