from pathlib import Path
from shutil import copyfileobj


def concat_files(files: list[Path], output: Path, cut_header: bool = False):
    """Concatenates files from the `files` argument into the output file.

    Args:
        files (list[Path]): list of file paths to concatenate
        output (Path): output file name
        cut_header (bool): Use if cutting a header from tabular data to 
            keep a single header in the concatenated output. Defaults to False.
    """
    # cat *.hmm > output
    with output.open("w") as outfile:
        for i, file in enumerate(files):
            with file.open("r") as infile:
                if cut_header and i > 0:
                    infile.readline()
                copyfileobj(infile, outfile)
