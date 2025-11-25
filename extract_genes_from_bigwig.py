#!/usr/bin/env python3
import argparse
import csv
import json
from pathlib import Path
import gzip
from pyBigWig import open as bwopen

# helper: load annotations

def load_annotations(gff_path):
    feats = {}
    opener = gzip.open if gff_path.suffix == '.gz' else open
    with opener(gff_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('	')
            if len(parts) < 9:
                continue
            seqid, src, ftype, start, end, score, strand, phase, attrs = parts
            start, end = int(start), int(end)
            attrd = {}
            for a in attrs.split(';'):
                if '=' in a:
                    k, v = a.split('=',1)
                    attrd[k] = v
            gid = attrd.get('ID') or attrd.get('Parent')
            if gid:
                feats.setdefault(gid, []).append({
                    'seqid': seqid, 'type': ftype, 'start': start,
                    'end': end, 'strand': strand
                })
    return feats

# main

def main():
    p = argparse.ArgumentParser()
    p.add_argument('pangenome_folder')
    p.add_argument('pangenome_index')
    p.add_argument('target_genes')
    p.add_argument('bigwig_folder')
    p.add_argument('bigwig_index')
    p.add_argument('output_json')
    p.add_argument('--type', required=True)
    p.add_argument('--upstream', type=int, default=0)
    p.add_argument('--downstream', type=int, default=0)
    p.add_argument('--whole-seq', action='store_true')
    p.add_argument('--inner-start', type=int, default=0)
    p.add_argument('--inner-end', type=int, default=0)
    p.add_argument('--pad', type=int, default=0)
    args = p.parse_args()

    pangenome_folder = Path(args.pangenome_folder)
    bigwig_folder = Path(args.bigwig_folder)

    # read pangenome index
    index = {}
    with open(args.pangenome_index) as f:
        for row in csv.DictReader(f):
            index[row['genotype']] = row

    # read bigwig index
    bwindex = {}
    with open(args.bigwig_index) as f:
        for row in csv.DictReader(f):
            g = row['genotype']
            bwindex.setdefault(g, []).append({
                'condition': row['condition'],
                'path': bigwig_folder / row['bigwig']
            })

    # read target genes
    tg = []
    genotypes_needed = set()
    with open(args.target_genes) as f:
        for row in csv.DictReader(f):
            tg.append(row)
            for k, v in row.items():
                if k.startswith('gene_ID_') and v.strip():
                    genotypes_needed.add(k.replace('gene_ID_', ''))

    output = []

    for g in sorted(genotypes_needed):
        if g not in index or g not in bwindex:
            continue

        gff_path = pangenome_folder / index[g]['annotation']
        feats = load_annotations(gff_path)

        for bwinfo in bwindex[g]:
            condition = bwinfo['condition']
            bw = bwopen(str(bwinfo['path']))

            for row in tg:
                gene_id = row.get(f'gene_ID_{g}', '').strip()
                if not gene_id:
                    continue
                gene_name = row.get('gene_name', '')

                if gene_id not in feats:
                    continue
                regions = [f for f in feats[gene_id] if f['type'] == args.type]
                if not regions:
                    continue

                for rgn in regions:
                    seqid = rgn['seqid']
                    start, end, strand = rgn['start'], rgn['end'], rgn['strand']

                    if args.whole_seq:
                        istart, iend = start, end
                    else:
                        istart = start + args.inner_start
                        iend = end - args.inner_end

                    up_s = max(1, istart - args.upstream) if args.upstream > 0 else istart
                    dn_e = iend + args.downstream if args.downstream > 0 else iend

                    up_vals = []
                    mid_vals = []
                    dn_vals = []

                    if not (args.upstream == 0 and args.inner_start == 0) and up_s < istart:
                        up_vals = bw.values(seqid, up_s-1, istart-1) or []
                    mid_vals = bw.values(seqid, istart-1, iend) or []
                    if not (args.downstream == 0 and args.inner_end == 0) and dn_e > iend:
                        dn_vals = bw.values(seqid, iend, dn_e) or []

                    pad_block = [None] * args.pad if args.pad > 0 else []
                    seq_vals = up_vals + pad_block + mid_vals + pad_block + dn_vals

                    if strand == '-':
                        seq_vals = list(reversed(seq_vals))

                    output.append({
                        'genotype': g,
                        'condition': condition,
                        'gene_name': gene_name,
                        'gene_id': gene_id,
                        'values': seq_vals
                    })

            bw.close()

        feats.clear()

    with open(args.output_json, 'w') as o:
        json.dump(output, o, indent=2)

if __name__ == '__main__':
    main()
