#!/usr/bin/env python3
import argparse
import multiprocessing
from dataclasses import dataclass
from itertools import repeat
from pathlib import Path
from typing import Optional

from operon_finder import hmmsearch, prodigal, summarize, utils


@dataclass
class Args:
    input: list[Path]
    outdir: Path
    databases: list[Path]
    operon_name: str
    metadata: Path
    distance_threshold: int
    jobs: int
    anchor_gene: Optional[str] = None

    @classmethod
    def from_namespace(cls, namespace: argparse.Namespace):
        return cls(**vars(namespace))


def parse_args() -> Args:
    parser = argparse.ArgumentParser(description="Operon finder")
    required_args = parser.add_argument_group("REQUIRED")
    required_args.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="fasta-formatted genomes",
        type=Path,
        metavar="FILE(s)",
    )
    required_args.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="output directory",
        type=Path,
        metavar="DIR",
    )
    required_args.add_argument(
        "-d",
        "--databases",
        nargs="+",
        required=True,
        help="name of HMM databases",
        type=Path,
        metavar="FILE(s)",
    )
    required_args.add_argument(
        "-n",
        "--operon-name",
        required=True,
        help="name of operon",
        metavar="STR",
    )
    required_args.add_argument(
        "-m",
        "--metadata",
        required=True,
        help="tab-delimited metadata file that connects hmm names gene names and desired color in plots",
        type=Path,
        metavar="FILE",
    )

    operon_args = parser.add_argument_group("OPERON ARGS")
    operon_args.add_argument(
        "-a",
        "--anchor-gene",
        metavar="STR",
        help="anchor gene to orient operon plots to",
    )
    operon_args.add_argument(
        "-t",
        "--distance-threshold",
        type=int,
        metavar="INT",
        default=summarize.SENTINEL_DISTANCE_THRESHOLD,
        help="minimum allowed distance in units of orfs to connect genes in an operon (default: %(default)s = no distance filtering)",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=10,
        help="number of parallel jobs to run (default: %(default)s)",
        metavar="INT",
    )
    args = parser.parse_args()
    return Args.from_namespace(args)


def main():
    args = parse_args()

    fastafiles = args.input
    databases = args.databases
    outdir = args.outdir
    metadata = args.metadata
    operon_name = args.operon_name
    distance_thresh = args.distance_threshold
    jobs = args.jobs
    anchor_gene = args.anchor_gene

    outdir.mkdir(exist_ok=True)
    # check number of input hmm databases
    if len(databases) == 1:
        database = databases[0]
    else:
        database = outdir.joinpath(f"{operon_name}.hmm")
        utils.concat_files(databases, database)

    jobs = min(len(fastafiles), jobs)
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
    main()
