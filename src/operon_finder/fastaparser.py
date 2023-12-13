#!/usr/bin/env python3
import textwrap
from pathlib import Path
from typing import Iterator, List, TextIO, Tuple, Union


def fastaparser(file: Union[str, Path]) -> Iterator[Tuple[str, str]]:
    """Parse a FASTA file by yielding (header, sequence) tuples.

    Args:
        file (str | Path): path to a valid FASTA file

    Yields:
        Tuple[str, str]: (header, sequence) tuple
    """
    with open(file) as f:
        header = f.readline().rstrip()[1:]
        seq: List[str] = list()
        for line in f:
            line = line.rstrip()
            if line[0] == ">":
                yield header, "".join(seq)
                header = line[1:]
                seq = list()
            else:
                seq.append(line)
        yield header, "".join(seq)


def wrap_fasta(sequence: str, width: int = 75) -> Iterator[str]:
    """Wrap a sequence when outputting to a FASTA file.

    Args:
        sequence (str): biological sequence
        width (int, optional): width of a sequence line. 
            Defaults to 75 characters.

    Yields:
        str: a single sequence line
    """
    yield from textwrap.wrap(sequence, width=width)


def write_fasta(fobj: TextIO, name: str, sequence: str, width: int = 75) -> None:
    """Write a fasta sequence to file with line wrapping for the sequence.

    Args:
        fobj (TextIO): open file object in text write mode
        name (str): name of fasta sequence
        sequence (str): fasta sequence
        width (int, optional): text wrapping width. Defaults to 75.
    """
    fobj.write(f">{name}\n")
    for seqline in wrap_fasta(sequence, width):
        fobj.write(f"{seqline}\n")
