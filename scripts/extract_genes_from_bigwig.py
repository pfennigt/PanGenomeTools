#!/usr/bin/env python3

"""
extract_bigwig_signals.py

Reads BigWig files and extracts float signal values for target genes.
Matches the coordinate logic of the FASTA extractor:
    (upstream → inner-start) + PAD + (inner-end → downstream)

Features:
    - Reads multiple BigWig tracks per genotype (as listed in bigwig_index.csv)
    - Adds condition label to each output entry
    - Uses pyBigWig for low-RAM, streaming extraction
    - Uses pangenome GFF + target_genes.csv + pangenome_index.csv
    - Coordinates clipped to chrom sizes
"""

import argparse
import csv
import gzip
import json
import sys
import re
from pathlib import Path
from typing import Dict, List
import pyBigWig


# ---------------------- argument parsing ----------------------

def parse_args():
    p = argparse.ArgumentParser(description="Extract BigWig signal values for target genes.")
    p.add_argument("--pangenome_folder", required=True)
    p.add_argument("--pangenome_index", required=True)
    p.add_argument("--bigwig_folder", required=True)
    p.add_argument("--bigwig_index", required=True)
    p.add_argument("--target_genes", required=True)
    p.add_argument("--output", required=True)

    p.add_argument("--type", default="gene")
    p.add_argument("--upstream", type=int, default=0)
    p.add_argument("--downstream", type=int, default=0)
    p.add_argument("--whole-seq", action="store_true")
    p.add_argument("--inner-start", type=int, default=0)
    p.add_argument("--inner-end", type=int, default=0)
    p.add_argument("--pad", type=int, default=0, help="Pad with this many null values between segments")

    return p.parse_args()


# ---------------------- helpers ----------------------

def read_csv_index(path: Path) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_pangenome_index(path: Path) -> Dict[str, Dict[str, str]]:
    data = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            data[r["genotype"]] = {
                "annotation": r["annotation"],
                "assembly": r["assembly"],
            }
    return data


def read_target_genes(path: Path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        geno_cols = [h for h in headers if h.startswith("gene_ID_")]
        for r in reader:
            rows.append(r)
    genotypes = [x.replace("gene_ID_", "") for x in geno_cols]
    return genotypes, rows


# ---------------------- GFF parsing ----------------------

def open_gff(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_gff(path: Path, feature_type: str):
    with open_gff(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seqid, src, ftype, start, end, score, strand, phase, attrs = cols
            if ftype != feature_type:
                continue
            start = int(start)
            end = int(end)
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v
            yield seqid, start, end, strand, attr_dict


def find_feature(path: Path, feature_type: str, target_id: str):
    if not target_id:
        return None

    esc = re.escape(target_id)
    pat = re.compile(r"(^|[;\t,])" + esc + r"($|[;\t,])")

    for seqid, s, e, strand, attr in parse_gff(path, feature_type):
        for v in attr.values():
            if v == target_id or pat.search(v):
                return seqid, s, e, strand, attr
    return None


# ---------------------- coord helpers ----------------------

def clip_interval(a: int, b: int, maxlen: int):
    a2 = max(1, min(a, maxlen))
    b2 = max(1, min(b, maxlen))
    if a2 <= b2:
        return a2, b2
    return 1, 0


# ---------------------- main ----------------------

def main():
    args = parse_args()

    pan_folder = Path(args.pangenome_folder)
    bigwig_folder = Path(args.bigwig_folder)

    pangenome_index = read_pangenome_index(Path(args.pangenome_index))
    bigwig_index_rows = read_csv_index(Path(args.bigwig_index))
    genotypes_in_targets, target_rows = read_target_genes(Path(args.target_genes))

    # build genotype → list of bigwig tracks
    bw_map: Dict[str, List[Dict[str, str]]] = {}
    for r in bigwig_index_rows:
        geno = r["genotype"]
        bw_map.setdefault(geno, []).append({
            "condition": r["condition"],
            "path": bigwig_folder / r["bigwig"],
        })

    # collect outputs here (one JSON object per row+condition)
    outputs = []

    # ---------------- PROCESS BY GENOTYPE (low RAM) ----------------
    for geno in genotypes_in_targets:

        if geno not in pangenome_index:
            print(f"Warning: genotype {geno} missing in pangenome index", file=sys.stderr)
            continue

        gff_path = pan_folder / pangenome_index[geno]["annotation"]

        if geno not in bw_map:
            print(f"Warning: genotype {geno} has no BigWigs", file=sys.stderr)
            continue

        # per BigWig file
        bw_files = bw_map[geno]

        # process each requested gene
        for row in target_rows:

            gene_name = row.get("gene_name", "")
            gene_id = row.get(f"gene_ID_{geno}", "")
            if not gene_id:
                continue

            # locate gene feature
            feat = find_feature(gff_path, args.type, gene_id)
            if feat is None:
                print(f"Warning: gene_id {gene_id} not found in {geno}", file=sys.stderr)
                continue

            seqid, start, end, strand, attr = feat

            # interior coordinate logic
            bstart, bend = (end, start) if strand == "-" else (start, end)

            if args.whole_seq:
                left_a, left_b = min(start, end), max(start, end)
                right_a, right_b = 1, 0
            else:
                left_a = bstart - args.upstream
                left_b = bstart + args.inner_start - 1

                right_a = bend - args.inner_end + 1
                right_b = bend + args.downstream

            # for each BigWig track for this genotype
            for bwinfo in bw_files:

                condition = bwinfo["condition"]
                bwpath = bwinfo["path"]

                try:
                    bw = pyBigWig.open(str(bwpath))
                except Exception as e:
                    print(f"Could not open BigWig {bwpath}: {e}", file=sys.stderr)
                    continue

                if seqid not in bw.chroms():
                    print(f"Warning: seqid {seqid} missing in {bwpath}", file=sys.stderr)
                    bw.close()
                    continue

                chrom_len = bw.chroms()[seqid]

                # clip intervals
                ll, lh = clip_interval(left_a, left_b, chrom_len)
                rl, rh = clip_interval(right_a, right_b, chrom_len)

                # extract signal
                left_vals = bw.values(seqid, ll - 1, lh) if ll <= lh else []
                right_vals = bw.values(seqid, rl - 1, rh) if rl <= rh else []

                # convert nan to None
                left_vals = [None if (v is None or (isinstance(v, float) and v != v)) else v for v in left_vals]
                right_vals = [None if (v is None or (isinstance(v, float) and v != v)) else v for v in right_vals]

                pad_list = [None] * args.pad

                combined = left_vals + pad_list + right_vals

                # orientation
                if strand == "-":
                    combined = combined[::-1]

                # store JSON object
                outputs.append({
                    "genotype": geno,
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "condition": condition,
                    "location": f"{seqid}:{start}-{end}({strand})",
                    "data": combined,
                })

                bw.close()

    # write JSON
    with open(args.output, "w") as fh:
        json.dump(outputs, fh, indent=2)


if __name__ == "__main__":
    main()
