"""Run prodigal orf prediction with a custom cleaning step for the fasta headers"""
import shlex
import subprocess
from pathlib import Path

from operon_finder.fastaparser import fastaparser, write_fasta


def run_prodigal(fastafile: Path, outdir: Path) -> Path:
    """Run prodigal on input genome fasta file

    Args:
        fastafile (Path): genome fasta file (.fna)
        outdir (Path): output directory for entire tool

    Returns:
        Path: predicted protein orfs fasta file (.faa)
    """
    outdir = outdir.joinpath("prodigal")
    outdir.mkdir(exist_ok=True)

    # fastafile       .fna
    # proteinsfile    .faa
    proteinsfile = outdir.joinpath(fastafile.with_suffix(".faa").name)
    output = outdir.joinpath(fastafile.with_suffix(".txt").name)

    if not proteinsfile.exists():
        command = f"prodigal -i {fastafile} -a {proteinsfile} -o {output} -q -m"

        subprocess.run(shlex.split(command))

    return proteinsfile


def clean_header(fastafile: Path) -> tuple[Path, Path]:
    """Clean fasta headers of the predicted prodigal orfs to only keep
    the genome and orf number. Additionally, this writes a summary file
    that maps the original header to the new cleaned header.

    Args:
        fastafile (Path): predicted protein orfs fasta file (.faa)

    Returns:
        tuple[Path, Path]: cleaned fasta file path, summary file path that
            contains the mapping from original to new headers
    """
    outdir = fastafile.parent
    output = outdir.joinpath(f"{fastafile.stem}_cleaned.faa")
    summaryfile = outdir.joinpath(f"{fastafile.stem}_orf_summary.tsv")

    genome = fastafile.stem.rsplit("_genomic")[0]

    if not output.exists():
        with output.open("w") as outfile, summaryfile.open("w") as summary:
            summary.write("genome\tprotein\tstart\tend\tstrand\n")
            for name, seq in fastaparser(fastafile):
                header = name.split(" # ")
                cleaned_name = header[0]
                start = header[1]
                end = header[2]
                strand = header[3]

                write_fasta(outfile, cleaned_name, seq)

                summary.write(f"{genome}\t{cleaned_name}\t{start}\t{end}\t{strand}\n")

    return output, summaryfile


def concat_fasta():
    # ignore for now
    pass


def main(fastafile: Path, outdir: Path) -> tuple[Path, Path]:
    """Run prodigal and clean the headers.

    Args:
        fastafile (Path): genome fasta file (.fna)
        outdir (Path): output directory for entire tool

    Returns:
        tuple[Path, Path]: cleaned fasta file path, summary file path that
            contains the mapping from original to new headers
    """
    proteinsfile = run_prodigal(fastafile, outdir)
    cleaned_output, orf_summary = clean_header(proteinsfile)
    return cleaned_output, orf_summary
