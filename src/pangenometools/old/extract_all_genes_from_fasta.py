#!/usr/bin/env python3
import sys
import subprocess
import pandas as pd
from pathlib import Path
from tqdm.auto import tqdm
from tqdm.contrib import DummyTqdmFile
import contextlib
import argparse

# --- AGAT base command ---
AGAT_CMD = (
    "agat_sp_extract_sequences.pl "
    "-g {annotation} -f {genome} -t {type} "
    "{extra} "
    "-o {output}"
)

@contextlib.contextmanager
def std_out_err_redirect_tqdm():
    orig_out_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = map(DummyTqdmFile, orig_out_err)
        yield orig_out_err[0]
    except Exception as exc:
        raise exc
    finally:
        sys.stdout, sys.stderr = orig_out_err


def extract_transcriptome(pangenome_index: pd.DataFrame,
                          pangenome_dir: Path,
                          output_dir: Path,
                          type:str,
                          extra_args: str,
                          tqdm_file=None):
    """Run AGAT for each genome listed in the index."""
    
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    for genotype, row in tqdm(pangenome_index.iterrows(),
                              total=pangenome_index.shape[0],
                              file=tqdm_file):

        annotation = pangenome_dir / row["annotation"]
        genome = pangenome_dir / row["assembly"]
        out_fa = output_dir / f"{genotype}.fa"

        if out_fa.exists() and out_fa.stat().st_size > 0:
            print(f"✅ {genotype}: already extracted ({out_fa})", file=sys.stdout)
            continue

        if not annotation.exists():
            print(f"⚠️ Missing annotation: {annotation}", file=sys.stdout)
            continue
        if not genome.exists():
            print(f"⚠️ Missing genome: {genome}", file=sys.stdout)
            continue

        cmd = AGAT_CMD.format(
            annotation=annotation.resolve(),
            genome=genome.resolve(),
            type=type,
            extra=extra_args,
            output=out_fa
        )

        print(f"\n🔹 Extracting transcriptome for {genotype} ...", file=sys.stdout)

        try:
            subprocess.run(cmd, shell=True, check=True)

            # Move AGAT log file if it exists
            log_file = Path(f"{annotation.stem}.agat.log")
            if log_file.exists():
                log_file.rename(logs_dir / log_file.name)

            # Move AGAT log folder if it exists
            log_folder = Path(f"agat_log_{annotation.stem}")
            if log_folder.exists():
                log_folder.rename(logs_dir / log_folder.name)

        except subprocess.CalledProcessError as e:
            print(f"❌ AGAT failed for {genotype}: {e}", file=sys.stdout)


def main():
    parser = argparse.ArgumentParser(description="Run AGAT transcriptome extraction for multiple genomes.")

    parser.add_argument(
        "--index",
        required=True,
        help="Path to pangenome_index.csv"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory for extracted FASTA files"
    )

    parser.add_argument(
        "--pangenome-dir",
        default="",
        help="Base directory for genome and annotation paths (default: directory of pangenome_index)"
    )

    parser.add_argument(
        "--type",
        default="gene",
        help="Type of feature to extract (default: 'gene')"
    )

    parser.add_argument(
        "--extra",
        default="",
        help="Additional arguments to pass to AGAT (e.g. '--upstream 1000 --downstream 1000')"
    )

    args = parser.parse_args()

    index_path = Path(args.index)

    pangenome_dir = Path(args.pangenome_dir if args.pangenome_dir != "" else index_path.parent)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    with std_out_err_redirect_tqdm() as orig_stdout:
        print("🔧 Loading pangenome index ...", file=sys.stdout)
        pangenome_index = pd.read_csv(index_path, index_col=0)

        print("\n🚀 Extracting transcriptomes ...", file=sys.stdout)
        extract_transcriptome(
            pangenome_index,
            pangenome_dir,
            output_dir,
            args.type,
            args.extra,
            tqdm_file=orig_stdout
        )

        print("\n🎉 Extraction complete", file=sys.stdout)


if __name__ == "__main__":
    main()