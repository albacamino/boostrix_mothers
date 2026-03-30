#!/usr/bin/env python3

import os
import pickle as pkl

import pandas as pd

from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def _load_data(path_to_counts_csv, path_to_metadata_csv):
    # Load data
    counts = pd.read_csv(path_to_counts_csv, index_col=0)
    metadata = pd.read_csv(path_to_metadata_csv, index_col=0)

    counts = counts.T

    return counts, metadata

def _filter_data(counts, metadata):

    # Filter genes that have less than 10 reads counts in total
    genes_to_keep = counts.columns[counts.sum(axis=0) >= 10]
    counts = counts[genes_to_keep]

    conteo = metadata.groupby("Ext_ID").size()
    # Filter the Ext_ID that have 3 registers
    ext_ids_validos = conteo[conteo == 3].index
    metadata = metadata[metadata["Ext_ID"].isin(ext_ids_validos)]
    
    counts = counts.loc[metadata.index]
    
    return counts, metadata


if __name__ == "__main__":
    script_dir = Path(__file__).parent
    counts, metadata = _load_data(
        (script_dir / "data" / "counts.csv"),
        (script_dir / "data" / "metadata.csv")
    )

    counts, metadata = _filter_data(counts, metadata)


