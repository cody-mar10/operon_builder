"""Summarize hmmsearch results to locate operons"""
from enum import Enum
from functools import reduce
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord

KEEP_COLS = [
    "protein",
    "query_name",
    "full_E-value",
    "full_score",
    "dom_E-value",
    "dom_score",
]

COL_ORDER = [
    "genome",
    "scaffold",
    "protein_id",
    "protein",
    "start",
    "end",
    "strand",
    "query_name",
    "gene_name",
    "full_E-value",
    "full_score",
    "dom_E-value",
    "dom_score",
    "color",
]

SENTINEL_DISTANCE_THRESHOLD = -1


class Mode(str, Enum):
    genome = "genome"
    scaffold = "scaffold"
    orf = "orf"

    def __str__(self):
        return str(self.value)


def choose_best(hits: pd.DataFrame, output: Path, mode: Mode):
    if mode.value == Mode.genome:
        (
            hits.sort_values(
                by=["genome", "query_name", "full_score"], ascending=[True, True, False]
            )
            .drop_duplicates(subset=["genome", "query_name"], keep="first")
            .sort_values(by=["genome", "scaffold", "protein_id"])
            .to_csv(output, sep="\t", index=False)
        )
    elif mode.value == Mode.scaffold:
        pass
    else:
        # mode.value == Mode.orf
        pass


def filter_operons(hits: pd.DataFrame, distance_thresh: int = 20) -> pd.DataFrame:
    """Filter out genes that are greater than `distance_thresh` away
    from all other genes.

    Args:
        hits (pd.DataFrame): dataframe of hmmsearch hits for all genomes
        distance_thresh (int, optional): minimum allowed distance for
            genes belonging to an operon. Defaults to 20.

    Returns:
        pd.DataFrame: filtered dataframe of hits, removing the genes
            that are above the distance threshold away from all other
            genes
    """
    # distance to upstream gene
    upstream = (
        hits.groupby("scaffold")
        .rolling(window=2)["protein_id"]
        .apply(np.diff)
        .reset_index()
        .set_index("level_1")
        .rename_axis(None, axis=0)
        .rename({"protein_id": "upstream_dist"}, axis=1)
    )["upstream_dist"]

    # downstream gen
    downstream = (
        hits.groupby("scaffold", as_index=False)
        .rolling(window=2)["protein_id"]
        .apply(np.diff)
        .shift(-1)
        .reset_index()
        .set_index("level_1")
        .rename_axis(None, axis=0)
        .rename({"protein_id": "downstream_dist"}, axis=1)
    )["downstream_dist"]

    hits = reduce(
        lambda left, right: pd.concat([left, right], axis=1),
        [hits, upstream, downstream],
    ).fillna(distance_thresh + 1)

    if distance_thresh == SENTINEL_DISTANCE_THRESHOLD:
        return hits.query(
            "upstream_dist <= @distance_thresh | downstream_dist <= @distance_thresh"
        )

    return hits


def plot_operons(
    hits: pd.DataFrame, outdir: Path, operon: str, anchor_gene: Optional[str] = None
):
    """Make operon plots for all hits detected

    Args:
        hits (pd.DataFrame): filtered hits dataframe with faraway
            genes removed
        outdir (Path): image output directory
        operon (str): name of operon
        anchor_gene (Optional[str], optional): name of an anchor gene
            to orient all plots based on this gene. Defaults to None.
            If no gene is provided, do not re-orient.
    """
    genomes = hits.groupby("genome", as_index=False)
    for genome, df in genomes:
        scaffolds = df.groupby("scaffold", as_index=False)
        n_scaffolds = scaffolds.ngroups
        fig, axes = plt.subplots(
            nrows=n_scaffolds, ncols=1, figsize=(10, 2 * n_scaffolds)
        )
        fig.suptitle(f"Genome: {genome}")
        for (scaffold, scaf_df), ax in zip(
            scaffolds, np.array(axes, ndmin=1, copy=False)
        ):
            features = list()

            if anchor_gene is not None:
                invert = df.query("gene_name == @anchor_gene")["strand"].values[0] == -1
            else:
                invert = False

            for row in scaf_df.itertuples():
                feature = GraphicFeature(
                    start=row.start,
                    end=row.end,
                    strand=row.strand,
                    label=row.gene_name,
                    color=row.color,
                )
                features.append(feature)

            start = scaf_df.iloc[0]["start"]
            end = scaf_df.iloc[-1]["end"]
            record = GraphicRecord(
                first_index=start,
                sequence_length=end - start + 1,
                features=features,
            )
            record.plot(ax=ax, level_offset=0)
            if invert:
                ax.invert_xaxis()
            ax.set_title(f"Scaffold: {scaffold}")
        fig.tight_layout(pad=0.75)
        fig.savefig(outdir.joinpath(f"{genome}_{operon}.png").as_posix())


def main(
    hmmtbl_file: Path,
    orf_file: Path,
    output: Path,
    outdir: Path,
    metadata_file: Path,
    distance_thresh: int,
    operon: str,
    anchor_gene: Optional[str] = None,
):
    """Identify operons based on proximity in units of orf distance.
    Then make operon plots for all identified operons per scaffold per genome.

    Args:
        hmmtbl_file (Path): cleaned and concatenated hmmsearch table
        orf_file (Path): orf summary file
        output (Path): name of output for hits summary
        outdir (Path): entire results directory
        metadata_file (Path): file with metadata for the operon
        distance_thresh (int): minimum allowed distance in orf units
            to define an operon
        operon (str): name of operon
        anchor_gene (Optional[str], optional): name of an anchor gene
            to orient all plots based on this gene. Defaults to None.
            If no gene is provided, do not re-orient.
    """
    hmmtbl = pd.read_table(hmmtbl_file).rename({"target_name": "protein"}, axis=1)[
        KEEP_COLS
    ]

    if len(hmmtbl) > 0:
        hmmtbl[["scaffold", "protein_id"]] = hmmtbl["protein"].str.rsplit(
            "_", n=1, expand=True
        )
        hmmtbl = hmmtbl.astype({"protein_id": int})

        orf_summary = pd.read_table(orf_file)
        metadata = pd.read_table(metadata_file).rename(
            {"hmm_name": "query_name", "gene": "gene_name"}, axis=1
        )

        hits = (
            hmmtbl.merge(orf_summary, on="protein")
            .merge(metadata)
            .sort_values(by=["genome", "scaffold", "protein_id"])
        )[COL_ORDER]

        # choose_best(hits, output, mode)
        hits = filter_operons(hits, distance_thresh)
        hits.to_csv(output, sep="\t", index=False)

        imgdir = outdir.joinpath("operon_plots")
        imgdir.mkdir(exist_ok=True)
        plot_operons(hits, imgdir, operon, anchor_gene)
    else:
        print("No hmmsearch hits were detected, so no operons can be found.")
