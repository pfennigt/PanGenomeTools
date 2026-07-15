"""
FASTA CLI module for PanGenomeTools.

This module provides command-line interface functionality for FASTA sequence extraction.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List
from ..models.fasta import FastaHandler
from ..models.gff import GFFHandler
from pangenometools.utils import replace_genotypes

def setup_fasta_parser() -> argparse.ArgumentParser:
    """
    Set up argument parser for FASTA sequence extraction.

    Returns:
        Configured argument parser
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences from FASTA files using GFF coordinates."
    )

    # Required arguments
    parser.add_argument("--pangenome-folder", required=True,
                       help="Path to pangenome folder")
    parser.add_argument("--pangenome-index", required=True,
                       help="Path to pangenome index file")
    parser.add_argument("--target-genes", default=None,
                       help="Path to target genes CSV file")
    parser.add_argument("--output", required=True,
                       help="Output file path")

    # Sequence extraction options
    parser.add_argument("--feature-type", default="gene",
                       help="Type of feature to use for coordinates (default: gene)")
    parser.add_argument("--upstream", type=int, default=0,
                       help="Nucleotides to include upstream (5' direction)")
    parser.add_argument("--downstream", type=int, default=0,
                       help="Nucleotides to include downstream (3' direction)")
    parser.add_argument("--inner-start", type=int, default=0,
                       help="Nucleotides to include from start of feature")
    parser.add_argument("--inner-end", type=int, default=0,
                       help="Nucleotides to include from end of feature")
    parser.add_argument("--pad", type=int, default=0,
                       help="Number of Ns to pad between segments (use -1 to separately extract flanking regions)")
    parser.add_argument("--whole-seq", action="store_true",
                       help="Extract the whole sequence between start and end")
    parser.add_argument("--genotypes", default=None, nargs="+",
                       help="Genotype(s) to run the analysis for (all by default) ")

    # Feature handling options
    parser.add_argument("--merge-strategy", default="merge",
                       choices=["merge", "first", "all"],
                       help="Strategy for merging multiple features (default: merge)")
    parser.add_argument("--use-five-prime-direction", action="store_true",
                       help="Always interpret upstream/downstream in 5' direction regardless of strand")
    parser.add_argument("--per-gene-group", action="store_true",
                       help="Write sequences into files per gene group, not per genotype")

    parser.add_argument("--skip-short-chromosomes", action="store_true",
                       help="Skip genes where the chromosome is too short to fit the flanking regions (alternative: padding with N)")
    parser.add_argument("--skip-short-genes", action="store_true",
                       help="Skip genes where the gene is too short to fit the inner extraction regions (alternative: padding with N)")

    # Additional options
    parser.add_argument("--silent", action="store_true",
                       help="Suppress progress output")

    return parser

def read_target_genes(target_path: Path) -> tuple:
    """
    Read target genes from CSV file.

    Args:
        target_path: Path to target genes CSV file

    Returns:
        Tuple of (genotypes, rows)
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

def write_sequence_to_file(sequence, gene_id, genotype, gene_name, args, info, outfile, write_mode): #FIXME: When Padding is used on either side, this needs to be handled better

    if not sequence:
        if not args.silent:
            print(f"Warning: No sequence extracted for {gene_id} in {genotype}", file=sys.stderr)
        return False

    # If pad == -1 is given, extract the flanking sequences separately
    if args.pad == -1 and info['left_len'] >0 and info['right_len'] >0:
        # Split the sequence into left and right parts
        sequences = {"left": sequence[:info['left_len']], "right": sequence[info['left_len']:]}

        # Use effective pad = 0
        pad = 0

    else:
        # Use the combined sequence
        sequences = {"combined": sequence}
        pad=args.pad
    # Use the full sequence

    # Open the out FASTA file
    with open(outfile, write_mode) as out_fh: #FIXME
        for sequence_loc, sequence in sequences.items():
            # Create header
            header_location = [
                f"{info['chrom']}:{info['ll']}-{info['lh']}({info['strand']})" if (info['left_len'] >0 and sequence_loc in ["combined", "left"]) else None,
                f"{info['chrom']}:{info['rl']}-{info['rh']}({info['strand']})" if (info['right_len'] >0 and sequence_loc in ["combined", "right"]) and not args.whole_seq else None,
                ]
            # Remove either location if it is irrelevant
            header_location = [x for x in header_location if x is not None]
            # Join the locations
            header_location = "&".join(header_location)

            # Set the label for the sequence ID
            # Determine if the sequence is a promoter and/or terminator
            if sequence_loc == "combined":
                if (info['left_len'] >0 and info['right_len'] >0) or args.whole_seq:
                    label=""
                elif (info['left_len'] >0 and info['strand'] == "+") or (info['right_len'] >0 and info['strand'] == "-"):
                    label="_promoter"
                elif (info['left_len'] >0 and info['strand'] == "-") or (info['right_len'] >0 and info['strand'] == "+"):
                    label="_terminator"
                else:
                    raise RuntimeError(f"error in determining sequence type for {gene_id}")
            elif (sequence_loc == "left" and info['strand'] == "+") or (sequence_loc == "right" and info['strand'] == "-"):
                label="_promoter"
            elif (sequence_loc == "right" and info['strand'] == "+") or (sequence_loc == "left" and info['strand'] == "-"):
                label="_terminator"
            else:
                raise RuntimeError(f"error in determining sequence type for {gene_id}")

            # Get the extraction options

            ex_options = [
                f"upstream:{args.upstream}" if args.upstream is not None else "",
                f"inner_start:{args.inner_start}" if not args.whole_seq else "",
                f"inner_end:{args.inner_end}" if not args.whole_seq else "",
                f"downstream:{args.downstream}" if args.downstream is not None else "",
                "whole-seq:True" if args.whole_seq else "",
                "use-five-prime-direction:True" if args.use_five_prime_direction else "",
                f"pad:{int(pad)+info['additional_padding']}" if not args.whole_seq else "",
            ]
            ex_options = [x for x in ex_options if len(x)>0]

            # Join the options
            ex_options = "&".join(ex_options)

            # Create the header
            header = (
                f"{gene_id}{label} genotype={genotype} gene_name={gene_name} type={args.feature_type} "
                f"location={header_location} extraction_options={ex_options}"
            )

            # Write to output
            out_fh.write(f">{header}\n")
            for i in range(0, len(sequence), 80): #FIXME
                out_fh.write(f"{sequence[i:i+80]}\n")
    return True

# FASTA extraction using a target genes file
def extract_with_target_genes(fasta_handler, genotypes, target_rows, args):
    # Process each genotype
    for g, genotype in enumerate(genotypes):
        # Track the genes that were already extracted
        extracted_genes=set()

        # Process each target row -> Gene groups
        for r, row in enumerate(target_rows):

            gene_name = row.get("gene_name", "")

            # Determine the output file to write to
            if args.per_gene_group:
                outfile = Path(args.output) / f"{gene_name}.fa"
                if g == 0:
                    write_mode = "w"
            else:
                outfile = Path(args.output) / f"{genotype}.fa"
                if r == 0 :
                    write_mode = "w"

            gene_id_col = f"gene_ID_{genotype}"
            gene_id_raw = row.get(gene_id_col, "")

            if not gene_id_raw or gene_id_raw in ["", "[]", "[""]", "['']"]:
                continue

            # Handle list format
            if gene_id_raw.startswith("["):
                gene_id_raw = gene_id_raw[1:-1].split(",")
            else:
                gene_id_raw = [gene_id_raw]

            # Process each gene ID
            for gene_id_raw_item in gene_id_raw:

                gene_id = gene_id_raw_item.strip().strip("'\"")

                if not gene_id or gene_id in extracted_genes:
                    continue

                try:
                    # Extract sequence
                    sequence, info = fasta_handler.extract_sequence(
                        genotype=genotype,
                        gene_id=gene_id,
                        feature_type=args.feature_type,
                        upstream=args.upstream,
                        downstream=args.downstream,
                        inner_start=args.inner_start,
                        inner_end=args.inner_end,
                        merge_strategy=args.merge_strategy,
                        pad=args.pad,
                        whole_seq=args.whole_seq,
                        use_five_prime_direction=args.use_five_prime_direction,
                        skip_short_chromosomes=args.skip_short_chromosomes,
                        skip_short_genes=args.skip_short_genes,
                        return_info=True,
                        _use_cache = True
                    )

                    # Skip the gene if no sequence was extracted
                    if len(sequence) == 0:
                        continue

                    # Write sequence to file
                    _written= write_sequence_to_file(
                        sequence, 
                        gene_id, 
                        genotype, 
                        gene_name, 
                        args, 
                        info, 
                        outfile, 
                        write_mode
                    )
                    if _written:
                        write_mode = "a"

                        # Add the gene to the lsit of extracted genes
                        extracted_genes.add(gene_id)
                except Exception as e:
                    if not args.silent:
                        print(f"Error processing {gene_id} in {genotype}: {e}", file=sys.stderr)
                    continue

# FASTA extraction using a target genes file
def extract_all_genes(fasta_handler:FastaHandler, genotypes, args):
    # Process each genotype
    for g, genotype in enumerate(genotypes):
        # Set the name and write mdoe for the out FASTA file
        outfile = Path(args.output) / f"{genotype}.fa"
        write_mode = "w"


        # Get the matching "gene" GFF features for that genotype
        features, _ = fasta_handler.gff_handler.get_feature(genotype, "gene")

        for gene_id in features.ID:

            try:
                # Extract sequence
                sequence, info = fasta_handler.extract_sequence(
                    genotype=genotype,
                    gene_id=gene_id,
                    feature_type=args.feature_type,
                    upstream=args.upstream,
                    downstream=args.downstream,
                    inner_start=args.inner_start,
                    inner_end=args.inner_end,
                    merge_strategy=args.merge_strategy,
                    pad=args.pad,
                    whole_seq=args.whole_seq,
                    use_five_prime_direction=args.use_five_prime_direction,
                    skip_short_chromosomes=args.skip_short_chromosomes,
                    skip_short_genes=args.skip_short_genes,
                    return_info=True,
                    _use_cache = True
                )

                # Write sequence to file
                _written= write_sequence_to_file(
                    sequence, 
                    gene_id, 
                    genotype, 
                    "None", 
                    args, 
                    info, 
                    outfile, 
                    write_mode
                )
                if _written:
                    write_mode = "a"
            except Exception as e:
                if not args.silent:
                    print(f"Error processing {gene_id} in {genotype}: {e}", file=sys.stderr)
                continue

def extract_fasta_sequences(args: argparse.Namespace) -> None:
    """
    Extract FASTA sequences for target genes.

    Args:
        args: Parsed command line arguments
    """
    # Initialize handlers
    pangenome_folder = Path(args.pangenome_folder)
    pangenome_index = Path(args.pangenome_index)

    fasta_handler = FastaHandler(pangenome_folder, pangenome_index)

    # Read target genes or use all genotypes in the pangenome index
    if args.target_genes is not None:
        genotypes, target_rows = read_target_genes(Path(args.target_genes))

        # Replace the genotypes if provided
        genotypes = replace_genotypes(genotypes, args.genotypes)

        # Extract using a target genes file
        extract_with_target_genes(fasta_handler, genotypes, target_rows, args)

    else:
        if args.per_gene_group:
            raise ValueError("A target_genes file must be provided to use per-gene-group")
        # Use all genotypes in the pangenome index
        genotypes = fasta_handler.pangenome_index.keys()

        # Replace the genotypes if provided
        genotypes = replace_genotypes(genotypes, args.genotypes)

        # Extract using a target genes file
        extract_all_genes(fasta_handler, genotypes, args)

    # Close the FASTA
    fasta_handler.close_fasta()


def fasta_cli():
    """
    Main entry point for FASTA CLI.

    Parses arguments and executes FASTA sequence extraction.
    """
    parser = setup_fasta_parser()
    args = parser.parse_args()

    try:
        extract_fasta_sequences(args)
        if not args.silent:
            print(f"FASTA extraction completed. Output written to {args.output}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    fasta_cli()