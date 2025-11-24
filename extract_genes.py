#!/usr/bin/env python3

"""
extract_genes.py

Extracts specified subsequences around genes listed in a target CSV from multiple genotypes described
in a pangenome index. Writes a multi-FASTA with headers containing genotype and gene_name.

Usage (example):

python extract_genes.py \
  --pangenome_folder /path/to/pangenome \
  --pangenome_index pangenome_index.csv \
  --target_genes target_genes.csv \
  --output extracted_sequences.fasta \
  --type gene \
  --upstream 1000 --downstream 200 --inner-start 50 --inner-end 50 --pad NNN

Notes / behavior:
- The script reads assembly FASTA files from <pangenome_folder>/Assembly and GFF files from
  <pangenome_folder>/Annotation/GFF as specified in the pangenome_index.csv file.
- The pangenome_index.csv must contain columns: genotype,annotation,assembly
- The target_genes.csv must contain at least a column `gene_name` and one or more columns of the
  form `gene_ID_<genotype>` where `<genotype>` matches the "genotype" value in pangenome_index.csv.
- Coordinates are handled **relative to the gene's biological orientation**:
    * upstream is upstream of the gene's biological start (5' side of the transcript)
    * downstream is downstream of the gene's biological end (3' side of the transcript)
    * inner-start is measured from the biological start toward the biological end
    * inner-end is measured from the biological end toward the biological start
- If --whole-seq is provided, inner-start and inner-end are ignored and the entire feature of
  type `--type` is considered the internal region.

This script is pure Python and uses Biopython for FASTA handling.
"""

import argparse
import csv
import gzip
import sys
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except Exception as e:
    sys.exit('Biopython is required (pip install biopython). Error: ' + str(e))


def parse_args():
    p = argparse.ArgumentParser(description='Extract sequences around genes from a pangenome.')
    p.add_argument('--pangenome_folder', required=True, help='Path to pangenome folder')
    p.add_argument('--pangenome_index', required=True, help='Path to pangenome_index.csv')
    p.add_argument('--target_genes', required=True, help='Path to target_genes.csv')
    p.add_argument('--output', required=True, help='Output FASTA file path')
    p.add_argument('--type', default='gene', help='GFF feature type to extract (e.g., gene, mRNA, CDS)')
    p.add_argument('--upstream', type=int, default=0, help='N nucleotides upstream of the biological start')
    p.add_argument('--downstream', type=int, default=0, help='N nucleotides downstream of the biological end')
    p.add_argument('--whole-seq', action='store_true', help='Use whole feature as internal region (ignore inner-start/inner-end)')
    p.add_argument('--inner-start', type=int, default=0, help='N nucleotides downstream from biological start')
    p.add_argument('--inner-end', type=int, default=0, help='N nucleotides upstream from biological end')
    p.add_argument('--pad', default='', help='String to place between the left and right extracted segments')
    return p.parse_args()


def read_index(index_path: Path) -> Dict[str, Dict[str, str]]:
    """Return mapping genotype -> {'annotation': <gff_filename>, 'assembly': <fasta_filename>}"""
    d = {}
    with open(index_path, newline='') as fh:
        reader = csv.DictReader(fh)
        required = {'genotype', 'annotation', 'assembly'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'pangenome_index.csv must contain columns: {required}')
        for r in reader:
            geno = r['genotype']
            d[geno] = {'annotation': r['annotation'], 'assembly': r['assembly']}
    return d


def read_target_genes(target_path: Path) -> Tuple[List[str], List[Dict[str,str]]]:
    """Return tuple(column_genotypes, rows)
    column_genotypes: list of genotype names deduced from columns gene_ID_<genotype>
    rows: list of dicts for each row with keys as column headers
    """
    rows = []
    with open(target_path, newline='') as fh:
        reader = csv.DictReader(fh)
        headers = reader.fieldnames or []
        # detect gene_ID_<genotype> columns
        geno_cols = []
        for h in headers:
            if h.startswith('gene_ID_'):
                geno_cols.append(h)
        if 'gene_name' not in headers:
            raise ValueError('target_genes.csv must contain a "gene_name" column')
        for r in reader:
            rows.append(r)
    # deduce genotype names
    genotypes = [h.replace('gene_ID_', '') for h in geno_cols]
    return genotypes, rows


def open_gff(gff_path: Path):
    if str(gff_path).endswith('.gz'):
        return gzip.open(gff_path, 'rt')
    else:
        return open(gff_path, 'r')


def parse_gff_features(gff_path: Path, feature_type: str):
    """Yield tuples (seqid, start, end, strand, attributes_dict) for features of feature_type."""
    with open_gff(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_str = parts[:9]
            if ftype != feature_type:
                continue
            start_i = int(start)
            end_i = int(end)
            # parse attributes
            attr = {}
            for item in attr_str.split(';'):
                if not item:
                    continue
                if '=' in item:
                    k, v = item.split('=', 1)
                elif ' ' in item:
                    k, v = item.split(' ', 1)
                else:
                    k, v = item, ''
                attr[k] = v
            yield seqid, start_i, end_i, strand, attr


def find_feature_for_id(gff_path: Path, feature_type: str, target_id: str) -> Optional[Tuple[str,int,int,str,Dict[str,str]]]:
    """Search GFF for a feature where any attribute value equals or contains target_id. Return first match."""
    if target_id is None or target_id == '':
        return None
    # treat target_id literally but match either exact or as semicolon-separated value
    escaped = re.escape(target_id)
    pattern = re.compile(r'(^|[;\t,])' + escaped + r'($|[;\t,])')
    for seqid, start, end, strand, attr in parse_gff_features(gff_path, feature_type):
        for v in attr.values():
            if v == target_id or pattern.search(v):
                return seqid, start, end, strand, attr
    return None


def clip_coords(a: int, b: int, seq_len: int) -> Tuple[int,int]:
    # ensure 1-based inclusive and within [1, seq_len]
    aa = max(1, min(a, seq_len))
    bb = max(1, min(b, seq_len))
    if aa <= bb:
        return aa, bb
    else:
        # return empty interval
        return 1, 0


def extract_region(seq_record, start:int, end:int) -> Seq:
    # seq_record is a Biopython SeqRecord, positions are 1-based inclusive
    # python slicing is 0-based, end exclusive
    if start > end:
        return Seq('')
    return seq_record.seq[start-1:end]


def main():
    args = parse_args()
    pangenome_folder = Path(args.pangenome_folder)
    index = read_index(Path(args.pangenome_index))
    genotypes_in_targets, target_rows = read_target_genes(Path(args.target_genes))

    # Build map genotype -> files
    geno_files = {}
    for g in genotypes_in_targets:
        if g not in index:
            print(f'Warning: genotype {g} not found in pangenome_index.csv - skipping', file=sys.stderr)
            continue
        ann = Path(pangenome_folder) / index[g]['annotation']
        asm = Path(pangenome_folder) / index[g]['assembly']
        if not ann.exists():
            print(f'Warning: annotation file {ann} for genotype {g} not found', file=sys.stderr)
        if not asm.exists():
            print(f'Warning: assembly file {asm} for genotype {g} not found', file=sys.stderr)
        geno_files[g] = {'gff': ann, 'fasta': asm}

    # Pre-index fasta files with SeqIO.index for random access
    fasta_index_cache = {}

    out_fh = open(args.output, 'w')

    # For each row in target rows, iterate genotype columns
    for row in target_rows:
        gene_name = row.get('gene_name', '')
        for col in [h for h in row.keys() if h and h.startswith('gene_ID_')]:
            g = col.replace('gene_ID_', '')
            gene_id = row.get(col, '')
            if not gene_id:
                continue
            if g not in geno_files:
                print(f'Warning: genotype {g} not available - skipping {gene_id}', file=sys.stderr)
                continue
            gff_path = geno_files[g]['gff']
            fasta_path = geno_files[g]['fasta']
            # find feature
            feat = find_feature_for_id(gff_path, args.type, gene_id)
            if feat is None:
                print(f'Warning: feature of type {args.type} with id {gene_id} not found in {gff_path} (genotype {g})', file=sys.stderr)
                continue
            seqid, start, end, strand, attr = feat

            # biological coordinates depend on strand
            if strand == '-':
                bstart = end
                bend = start
            else:
                bstart = start
                bend = end

            # inner region handling
            if args.whole_seq:
                inner_start = 0
                inner_end = 0
            else:
                inner_start = args.inner_start
                inner_end = args.inner_end
                # compute inner genomic coords later per segments

            # left segment: biological start - upstream  ..  biological start + inner_start - 1
            left_a = bstart - args.upstream
            left_b = bstart + inner_start - 1

            # right segment: biological end - inner_end + 1  .. biological end + downstream
            right_a = bend - inner_end + 1
            right_b = bend + args.downstream

            # if whole_seq: produce left as beginning..bstart-1? Simpler: when whole_seq, take entire feature as inner region
            if args.whole_seq:
                # When whole_seq, we set left and right to cover the whole feature split by pad at the internal boundary
                # We'll produce left = full feature and right = empty (so pad will not be used)
                left_a = min(start, end)
                left_b = max(start, end)
                right_a = 1
                right_b = 0

            # index fasta if not already
            if fasta_path not in fasta_index_cache:
                try:
                    fasta_index_cache[fasta_path] = SeqIO.index(str(fasta_path), 'fasta')
                except Exception as e:
                    print(f'Error indexing fasta {fasta_path}: {e}', file=sys.stderr)
                    continue
            seq_index = fasta_index_cache[fasta_path]

            # find seq record id matching seqid
            if seqid in seq_index:
                seqrec = seq_index[seqid]
            else:
                # try some common fallbacks
                alt = seqid
                if seqid.startswith('chr'):
                    alt = seqid.replace('chr', '')
                else:
                    alt = 'chr' + seqid
                if alt in seq_index:
                    seqrec = seq_index[alt]
                else:
                    print(f'Warning: seqid {seqid} not found in fasta {fasta_path} (tried {seqid} and {alt})', file=sys.stderr)
                    continue

            seqlen = len(seqrec.seq)
            # clip and extract left and right
            left_low, left_high = clip_coords(left_a, left_b, seqlen)
            right_low, right_high = clip_coords(right_a, right_b, seqlen)

            left_seq = extract_region(seqrec, left_low, left_high) if left_low <= left_high else Seq('')
            right_seq = extract_region(seqrec, right_low, right_high) if right_low <= right_high else Seq('')

            # If both upstream AND inner-start are 0, then no start should be extracted (i.e., left_seq empty)
            if args.upstream == 0 and inner_start == 0 and not args.whole_seq:
                left_seq = Seq('')
            # If both downstream AND inner-end are 0, then no end should be extracted
            if args.downstream == 0 and inner_end == 0 and not args.whole_seq:
                right_seq = Seq('')

            # Compose result sequence: left + pad + right
            pad_seq = Seq(args.pad) if args.pad else Seq('')
            combined = left_seq + pad_seq + right_seq

            # If feature is on '-' strand, return reverse complement of the combined sequence so the result is in
            # transcript (5'->3') orientation.
            if strand == '-':
                combined = combined.reverse_complement()

            # Build header
            header_fields = [f'genotype={g}', f'gene_name={gene_name}']
            # include gene_id and coordinates for convenience
            header_fields.append(f'gene_id={gene_id}')
            header_fields.append(f'location={seqid}:{start}-{end}({strand})')
            header = '|' .join(header_fields)

            # write FASTA
            out_fh.write(f'>{header}\n')
            # wrap at 80 chars
            seq_str = str(combined)
            for i in range(0, len(seq_str), 80):
                out_fh.write(seq_str[i:i+80] + '\n')

    out_fh.close()
    # close fasta index objects
    for idx in fasta_index_cache.values():
        try:
            idx.close()
        except Exception:
            pass


if __name__ == '__main__':
    main()
