#!/usr/bin/env python3

"""
Extracts portions of GFF files related to a target gene
"""

import argparse
import csv
import gzip
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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

def parse_gff_features(gff_path: Path, feature_type: "str|None"):
    with open_gff(gff_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if feature_type is not None and ftype != feature_type:
                continue

            start_i = int(start)
            end_i = int(end)
            attr = {"feature_type":ftype}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr[k] = v

            yield seqid, start_i, end_i, strand, attr, line

def find_feature_for_id(gff_path: Path, feature_type: str, target_id: str, *, return_first=True, search_mode="strict", return_line=False):
    if not target_id:
        return None

    escaped = re.escape(target_id)

    if search_mode == "pattern":
        pattern = re.compile(escaped)
    elif search_mode=="children":
        pattern = re.compile(r"(^|[;\t,])" + escaped + r"($|[;\t,]|\.)")
    elif search_mode == "strict":
        pattern = re.compile(r"(^|[;\t,])" + escaped + r"($|[;\t,])")

    res = []

    for seqid, start, end, strand, attr, line in parse_gff_features(gff_path, feature_type):
        for v in attr.values():
            if v == target_id or pattern.search(v):
                if return_first:
                    if return_line:
                        return line
                    else:
                        return seqid, start, end, strand, attr
                else:
                    if return_line:
                        res.append(line)
                    else:
                        res.append((seqid, start, end, strand, attr))
                    break

    if not return_first and len(res) > 0:
        return res
    else:
        return None
    

def extract_genes_from_gff(gff_path: Path, feature_type: str, target_id: str, custom_coords=True, custom_coords_offset=1500):
    pass


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

            feat = find_feature_for_id(gff_path, args.type, gene_id)
            if feat is None:
                print(f"Warning: gene_id {gene_id} not found in {g}", file=sys.stderr)
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
            header = (
                f"genotype={g}|gene_name={gene_name}|gene_id={gene_id}|"
                f"location={seqid}:{start}-{end}({strand})"
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
