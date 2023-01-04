"""Run hmmsearch with a custom cleaning script for the output table"""
import re
import subprocess
import shlex

from pathlib import Path

HMMER_header = [
    "target_name",
    "target_accession",
    "query_name",
    "query_accession",
    "full_E-value",
    "full_score",
    "full_bias",
    "dom_E-value",
    "dom_score",
    "dom_bias",
    "exp",
    "reg",
    "clu",
    "ov",
    "env",
    "dom",
    "rep",
    "inc",
    "description",
]

SPACES_PATTERN = re.compile(r"\s+")


def run_hmmsearch(fastafile: Path, database: Path, outdir: Path) -> Path:
    """Run hmmsearch with the supplied protein fasta file and the hmm
    database

    Args:
        fastafile (Path): protein fasta file (.faa)
        database (Path): hmm database
        outdir (Path): output directory for entire tool

    Returns:
        Path: hmmsearch tabular output file path
    """
    outdir = outdir.joinpath("hmmsearch")
    outdir.mkdir(exist_ok=True)

    filename = f"{fastafile.stem}_{database.stem}.hmmtbl"
    hmmtbl = outdir.joinpath(filename)

    if not hmmtbl.exists():
        command = f"hmmsearch --noali --cut_tc --tblout {hmmtbl} -o /dev/null --cpu 0 {database} {fastafile}"

        subprocess.run(shlex.split(command))
    return hmmtbl


def clean_hmmtbl(hmmtbl: str | Path, output: Path, header: list[str]) -> None:
    """Clean hmmsearch table output convert it from space- to tab-delimited,
    fix the multiline headers, and remove the comments.

    Args:
        hmmtbl (str | Path): file path to hmmsearch tabular output
        output (Path): cleaned output path
        header (list[str]): contains desired names of the headers
    """
    if not output.exists():
        with open(hmmtbl) as infile, open(output, "w") as outfile:
            headerline = "\t".join(header)
            outfile.write(f"{headerline}\n")
            for line in infile:
                if line[0] != "#":
                    tabbed_line = SPACES_PATTERN.sub("\t", line.rstrip())
                    outfile.write(f"{tabbed_line}\n")


def main(fastafile: Path, database: Path, outdir: Path) -> Path:
    """Run hmmsearch on the input protein fasta file with the provided
    hmm databases. Then, cleans the tabular output.

    Args:
        fastafile (Path): protein fasta file (.faa)
        database (Path): hmm database
        outdir (Path): output directory for entire tool

    Returns:
        Path: cleaned hmmsearch tabular output file path
    """
    hmmtbl = run_hmmsearch(fastafile, database, outdir)
    cleaned_hmmtbl = hmmtbl.with_suffix(".cleaned.tsv")
    clean_hmmtbl(hmmtbl, cleaned_hmmtbl, HMMER_header)
    return cleaned_hmmtbl
