#!/usr/bin/env python3
import argparse
import multiprocessing
from itertools import repeat
from pathlib import Path
from typing import Optional

import hmmsearch
import prodigal
import summarize
import utils


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Operon finder")
    required_args = parser.add_argument_group("REQUIRED")
    required_args.add_argument(
        "-i", "--input", nargs="+", required=True, help="fasta-formatted genomes"
    )
    required_args.add_argument("-o", "--outdir", required=True, help="output directory")
    required_args.add_argument(
        "-d", "--databases", nargs="+", required=True, help="name of HMM databases"
    )
    required_args.add_argument(
        "-n", "--operon-name", required=True, help="name of operon"
    )
    required_args.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="tab-delimited metadata file that connects hmm names gene names and desired color in plots",
    )

    operon_args = parser.add_argument_group("OPERON ARGS")
    operon_args.add_argument(
        "-a",
        "--anchor-gene",
        help="anchor gene to orient operon plots to",
    )
    operon_args.add_argument(
        "-t",
        "--distance-threshold",
        type=int,
        default=20,
        help="minimum allowed distance in units of orfs to connect genes in an operon (default: %(default)s)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=10,
        help="number of parallel jobs to run (default: %(default)s)",
    )
    # parser.add_argument(
    #     "-s",
    #     "--summary-mode",
    #     choices=set(summarize.Mode),
    #     default=summarize.Mode.genome,
    #     help="DEPRECATED -- mode to choose the best annotation (default: %(default)s)",
    # )

    return parser.parse_args()


def main(
    fastafiles: list[Path],
    databases: list[Path],
    outdir: Path,
    operon_name: str,
    metadata: Path,
    distance_thresh: int,
    jobs: int,
    anchor_gene: Optional[str] = None,
):
    outdir.mkdir(exist_ok=True)
    # check number of input hmm databases
    if len(databases) == 1:
        database = databases[0]
    else:
        database = outdir.joinpath(f"{operon_name}.hmm")
        utils.concat_files(databases, database)

    with multiprocessing.Pool(jobs) as pool:
        prodigal_files = pool.starmap(prodigal.main, zip(fastafiles, repeat(outdir)))

        faa_files: list[Path] = list()
        orf_summaries: list[Path] = list()
        for faa_file, orf_summary in prodigal_files:
            faa_files.append(faa_file)
            orf_summaries.append(orf_summary)

        hmmtbls = pool.starmap(
            hmmsearch.main,
            zip(
                faa_files,
                repeat(database),
                repeat(outdir),
            ),
        )

    hmmtbl_summary = outdir.joinpath(f"{operon_name}_summary.tsv")
    utils.concat_files(hmmtbls, hmmtbl_summary, True)

    orf_table = outdir.joinpath("orf_summary.tsv")
    utils.concat_files(orf_summaries, orf_table, True)

    summary = outdir.joinpath(f"{operon_name}_filtered_summary.tsv")
    summarize.main(
        hmmtbl_summary,
        orf_table,
        summary,
        outdir,
        metadata,
        distance_thresh,
        operon_name,
        anchor_gene,
    )


if __name__ == "__main__":
    # check if running from cli
    args = parse_args()

    fastafiles = [Path(file) for file in args.input]
    databases = [Path(db) for db in args.databases]
    outdir = Path(args.outdir)
    metadata = Path(args.metadata)
    main(
        fastafiles,
        databases,
        outdir,
        args.operon_name,
        metadata,
        args.distance_threshold,
        args.jobs,
        args.anchor_gene,
    )
