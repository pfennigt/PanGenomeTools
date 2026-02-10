#!/usr/bin/env python3

import re
import pandas as pd
import numpy as np

from typing import Union

# Function to extract information from FASTA-like headers
def extract_fasta_header_like(header:str, key_prefix:str="", chunk_suffix:Union[bool,None]=None)->pd.Series:
    """Function to extract information in a FASTA-like header as produced by the "extract_genes" functions.

    Args:
        header (str): Header of a sequence
        key_prefix (str, optional): String to be added to the beginning of the output indices. Defaults to "".
        chunk_suffix (Union[bool,None], optional): Does the header contain a suffix indicating the chunk position produced by the DeepCistrome function 'split_sequence_chunks_to_250bp'. If None, the suffix is automatically detected. Defaults to None.

    Raises:
        ValueError: If the header structure dos not follow the convention '[>]<_id> [<key>=<value>] ... [<key>=<value>][_<chunk_suffix>]'.

    Returns:
        pd.Series: A Series object storing the header information as key:value pairs.
    """
    # Remove the initial ">"
    header= header.strip()
    if header.startswith(">"):
        header = header[1:]

    # Initialize an info Series to save the header information
    info = pd.Series(dtype=object)

    if chunk_suffix is None:
        chunk_suffix = re.search("_chunk[0-9_\-]+$", header)

    # If there is a chunk suffix, extract it
    if chunk_suffix:
        header, chunk_range = header.split("_chunk")
        chunk_start, chunk_end = np.array(chunk_range.split("_")[1].split("-")).astype(int)


    # Split the header by spaces
    header_spc = header.split(" ")

    # The ID is the first element
    id = header_spc.pop(0)
    info[key_prefix + "_id"] = id

    # Initialise a string to hold the current piece of information
    curr_info = []

    for i, header_field in enumerate(header_spc):
        # Add the new field to the current info
        curr_info.append(header_field)

        # If the next field includes a "=", the current info is done
        if len(header_spc) == i+1 or re.search("=", header_spc[i+1]):
            # Combine the current information to a string
            curr_info_str = " ".join(curr_info)

            if len(info) == 1 and not re.search("=", curr_info_str):
                info_key = "_names"
                info_val = curr_info_str
            elif re.search("=", curr_info_str):
                info_key, info_val = curr_info_str.split("=")
            else:
                raise ValueError(f"Unrecognized fasta header structure.\n<key>=<value> par expected in {curr_info_str}")

            # Add the current info to the dict
            info.loc[key_prefix + info_key] = info_val

            # Reset the current information list
            curr_info = []

    if chunk_suffix:
        info[key_prefix + "_chunk_start"] = chunk_start
        info[key_prefix + "_chunk_end"] = chunk_end

    return info

def extract_fasta_header_like_series(header:Union[str,list,pd.Series], key_prefix:str="", auto_key_prefix="_", chunk_suffix:Union[bool,None]=None)->pd.Series:
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

    # Split the header by spaces
    info = header.loc[:,"_hdr"].str.split(" ").apply(pd.Series)

    # Determine the column names form the first row
    col_names = [f"{auto_key_prefix}id"] + info.iloc[0,1:].str.split("=").apply(pd.Series).iloc[:,0].to_list()
    info.columns = [key_prefix + x for x in col_names]

    # Remove the prefixes
    for _name in col_names[1:]:
        info.loc[:, _name] = info.loc[:, _name].str.removeprefix(f"{_name}=")

    # Add the chunk information
    if chunk_suffix:
        info = pd.concat([info, chunk_range], axis=1)
    return info