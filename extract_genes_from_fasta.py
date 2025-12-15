#!/usr/bin/env python3

"""
extract_genes_pyfaidx.py

Memory-efficient extractor using pyfaidx for random-access FASTA slicing,
processing one genotype at a time.

The extraction logic is identical to your original script:
    (upstream → inner-start) + PAD + (inner-end → downstream)
with whole-seq optional, and proper strand handling.

pyfaidx ensures low RAM usage even for huge genomes.
"""

import argparse
import csv
import gzip
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import numpy as np

try:
    from pyfaidx import Fasta
except Exception as e:
    sys.exit("pyfaidx is required (pip install pyfaidx). Error: " + str(e))

from Bio.Seq import Seq  # only for reverse complement


# ---------------------- argument parsing ----------------------

def parse_args():
    p = argparse.ArgumentParser(description="Extract sequences around genes using pyfaidx.")
    p.add_argument("--pangenome_folder", required=True)
    p.add_argument("--pangenome_index", required=True)
    p.add_argument("--target_genes", required=True)
    p.add_argument("--output", required=True)

    p.add_argument("--type", default="gene")
    p.add_argument("--upstream", type=int, default=0)
    p.add_argument("--downstream", type=int, default=0)
    p.add_argument("--whole-seq", action="store_true")
    p.add_argument("--inner-start", type=int, default=0)
    p.add_argument("--inner-end", type=int, default=0)
    p.add_argument("--pad", default="")

    p.add_argument("--search-mode", default="children")

    return p.parse_args()


# ---------------------- pangenome index ----------------------

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


def read_target_genes(target_path: Path):
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


# ---------------------- GFF parsing ----------------------

def open_gff(gff_path: Path):
    return gzip.open(gff_path, "rt") if str(gff_path).endswith(".gz") else open(gff_path, "r")


def parse_gff_features(gff_path: Path, feature_type: str):
    with open_gff(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != feature_type:
                continue

            start_i = int(start)
            end_i = int(end)
            attr = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr[k] = v

            yield seqid, start_i, end_i, strand, attr

def find_feature_for_id(gff_path: Path, feature_type: str, target_id: str, search_mode:str="strict", return_all:bool=False):
    if not target_id:
        return None

    escaped = re.escape(target_id)

    if search_mode == "pattern":
        pattern = re.compile(escaped)
    elif search_mode == "strict":
        pattern = re.compile(r"(^|[;\t,])" + escaped + r"($|[;\t,])")
    elif search_mode == "children":
        pattern = re.compile(r"(^|[;\t,])" + escaped + r"($|[;\t,]|\.)")
    else:
        raise ValueError(f"undefined search_mode: '{search_mode}'")

    res = []

    for seqid, start, end, strand, attr in parse_gff_features(gff_path, feature_type):
        for v in attr.values():
            if v == target_id or pattern.search(v):
                if return_all:
                    res.append((seqid, start, end, strand, attr))
                else:
                    return [(seqid, start, end, strand, attr)]
    if return_all and len(res) > 0:
        return res
    else:
        return None


def merge_features(feats:Union[list,None], strategy="merge") -> Tuple:
    if feats is None:
        return None, None
    elif len(feats) == 1:
        # If there is only one entry, return it
        return feats[0], None
    elif strategy=="first":
        return feats[0], "_first"
    elif strategy == "merge":
        seqid, start, end, strand, attr = feats[0]

        # Iterate through the features
        for feat in feats[1:]:
            _seqid, _start, _end, _strand, _attr = feat

            # Assert that the ID ands strand must be the same for merging
            if _seqid != seqid:
                raise RuntimeError(f"differing seqids in selected features: {_seqid} != {seqid}")
            if _strand != strand:
                raise RuntimeError(f"differing strands in selected features")
            
            # Update start and end
            start = np.min([start, _start])
            end = np.max([end, _end])

            # Update the attributes
            attr.update(_attr)

        # Return the merged feature
        return (seqid, start, end, strand, attr), "_merged"
    else:
        raise RuntimeError("error in merging features")

# ---------------------- coordinate helpers ----------------------

def clip_coords(a: int, b: int, seq_len: int):
    """Clip a 1-based inclusive interval to valid FASTA length."""
    a2 = max(1, min(a, seq_len))
    b2 = max(1, min(b, seq_len))
    return (a2, b2) if a2 <= b2 else (1, 0)


# ---------------------- main extraction ----------------------

def main():
    args = parse_args()

    pangenome_folder = Path(args.pangenome_folder)
    index = read_index(Path(args.pangenome_index))
    genotypes_in_targets, target_rows = read_target_genes(Path(args.target_genes))

    # map genotype → files
    geno_files = {}
    for g in genotypes_in_targets:
        if g not in index:
            print(f"Warning: genotype {g} missing in index.", file=sys.stderr)
            continue
        geno_files[g] = {
            "gff": pangenome_folder / index[g]["annotation"],
            "fasta": pangenome_folder / index[g]["assembly"],
        }

    out_fh = open(args.output, "w")

    # --------- PROCESS ONE GENOTYPE AT A TIME (pyfaidx does streaming) ---------
    for g in genotypes_in_targets:

        if g not in geno_files:
            continue

        gff_path = geno_files[g]["gff"]
        fasta_path = geno_files[g]["fasta"]

        try:
            # pyfaidx loads only FASTA index; bases are streamed
            fa = Fasta(str(fasta_path), rebuild=False)
        except Exception:
            # rebuild if missing index
            fa = Fasta(str(fasta_path), rebuild=True)

        # Process all rows referring to this genotype
        for row in target_rows:

            gene_name = row.get("gene_name", "")
            gene_id = row.get(f"gene_ID_{g}", "")
            if not gene_id:
                continue

            # Get all the features
            feats = find_feature_for_id(gff_path, args.type, gene_id, search_mode=args.search_mode, return_all=True)

            # Merge the features
            feat, merge_text = merge_features(feats, strategy="merge")

            if feat is None:
                print(f"Warning: gene_id {gene_id} with feature {args.type} not found in {g}" + " consider sing search_mode 'children' or 'pattern'" if args.search_mode=="strict" else "", file=sys.stderr)
                continue

            seqid, start, end, strand, attr = feat

            # orientation
            bstart, bend = (end, start) if strand == "-" else (start, end)

            # determine inner offsets
            if args.whole_seq:
                inner_start, inner_end = 0, 0
            else:
                inner_start, inner_end = args.inner_start, args.inner_end

            # boundaries (exact original logic)
            left_a = bstart - args.upstream
            left_b = bstart + inner_start - 1

            right_a = bend - inner_end + 1
            right_b = bend + args.downstream

            if args.whole_seq:
                left_a = min(start, end)
                left_b = max(start, end)
                right_a, right_b = 1, 0

            # find the correct chromosome name
            if seqid in fa:
                chrom = seqid
            elif "chr" + seqid in fa:
                chrom = "chr" + seqid
            elif seqid.replace("chr", "") in fa:
                chrom = seqid.replace("chr", "")
            else:
                print(f"Warning: seqid {seqid} not found in FASTA for {g}", file=sys.stderr)
                continue

            seqlen = len(fa[chrom])

            # clip
            ll, lh = clip_coords(left_a, left_b, seqlen)
            rl, rh = clip_coords(right_a, right_b, seqlen)

            # extract (pyfaidx returns strings)
            left_seq = fa[chrom][ll - 1: lh].seq if ll <= lh else ""
            right_seq = fa[chrom][rl - 1: rh].seq if rl <= rh else ""

            print(f"ll:{ll}, lh:{lh}, rl:{rl}, rh:{rh}", file=sys.stdout)
            print(f"left_seq:{len(left_seq)}, right_seq:{len(right_seq)}", file=sys.stdout)

            # special zero rules
            if args.upstream == 0 and inner_start == 0 and not args.whole_seq:
                left_seq = ""
            if args.downstream == 0 and inner_end == 0 and not args.whole_seq:
                right_seq = ""

            # compose
            combined = left_seq + args.pad + right_seq

            # strand correction
            if strand == "-":
                combined = str(Seq(combined).reverse_complement())

            # write
            header_location = [
                f"{seqid}:{ll}-{lh}({strand})" if len(left_seq) >0 else None,
                f"{seqid}:{rl}-{rh}({strand})" if len(right_seq) >0 and not args.whole_seq else None,
                ]
            # Remove either location if it is irrelevant
            header_location = [x for x in header_location if x is not None]
            # Join the locations
            header_location = "&".join(header_location)

            # Set the label for the sequence ID
            # Determine if the sequence is a promoter and/or terminator
            if (len(left_seq) >0 and len(right_seq) >0) or args.whole_seq:
                label="flanking"
            elif (len(left_seq) >0 and strand == "+") or (len(right_seq) >0 and strand == "-"):
                label="promoter"
            elif (len(left_seq) >0 and strand == "-") or (len(right_seq) >0 and strand == "+"):
                label="terminator"
            else:
                raise RuntimeError(f"error in determining sequence type for {gene_id}")

            # Get the extraction options
            ex_options = [
                f"upstream:{args.upstream}" if args.upstream is not None else "-",
                f"inner_start:{args.inner_start}" if (args.inner_start is not None and not args.whole_seq) else "-",
                f"inner_end:{args.inner_end}" if (args.inner_end is not None and not args.whole_seq) else "-",
                f"downstream:{args.downstream}" if args.downstream is not None else "-",
            ]
            # Join the options
            ex_options = "&".join(ex_options)

            # Create the header
            header = (
                f"{gene_id}_{label} genotype={g} gene_name={gene_name} type={args.type}{merge_text} "
                f"location={header_location} extraction_options={ex_options}"
            )
            out_fh.write(f">{header}\n")
            for i in range(0, len(combined), 80):
                out_fh.write(combined[i:i+80] + "\n")

        # unload FASTA completely
        fa.close()
        del fa

    out_fh.close()


if __name__ == "__main__":
    main()
