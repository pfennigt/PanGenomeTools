#!/usr/bin/env python3

import argparse
from pyfaidx import Fasta


def make_pseudo_annotation(input_fasta, output_gff, feature_type="gene", upstream=1500, downstream=1500):
    fasta = Fasta(input_fasta)

    with open(output_gff, "w") as out:
        for record_name in fasta.keys():
            seq = fasta[record_name]
            length = len(seq)

            gff_line = [
                record_name,   # seqname
                ".",           # source
                feature_type,  # feature
                str(1 + upstream),           # start
                str(length - downstream),   # end
                ".",           # score
                "+",           # strand
                ".",           # frame
                f"ID={'_'.join(record_name.split('_')[:-1])}"            # attribute
            ]

            out.write("\t".join(gff_line) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Create pseudo GFF annotations from a FASTA file."
    )
    parser.add_argument("input_fasta", help="Path to input FASTA file")
    parser.add_argument("output_gff", help="Path to output GFF file")
    parser.add_argument(
        "--feature_type",
        default="gene",
        help='Feature type for GFF entries (default: "gene")'
    )
    parser.add_argument(
        "--upstream",
        default="1500",
        help='Number of nucleotides extracted upstream of the feature (default: 500)'
    )
    parser.add_argument(
        "--downstream",
        default="1500",
        help='Number of nucleotides extracted downstream of the feature (default: 500)'
    )

    args = parser.parse_args()

    make_pseudo_annotation(
        args.input_fasta,
        args.output_gff,
        args.feature_type,
        int(args.upstream),
        int(args.downstream),
    )


if __name__ == "__main__":
    main()
