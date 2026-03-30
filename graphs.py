#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def _load_data(path_to_counts_csv, path_to_metadata_csv):
    # Load data
    counts = pd.read_csv(path_to_counts_csv, index_col=0)
    metadata = pd.read_csv(path_to_metadata_csv, index_col=0)

    counts = counts.T

    return counts, metadata

def _filter_data(counts, metadata):

    common_samples = counts.index.intersection(metadata.index)

    counts = counts.loc[common_samples]
    metadata = metadata.loc[common_samples]
    
    return counts, metadata

def _plot_pca(counts, metadata, output_dir):
    # Asegurar alineación
    common_samples = counts.index.intersection(metadata.index)
    counts = counts.loc[common_samples]
    metadata = metadata.loc[common_samples]
    condition = metadata["Condition"].unique()[0]

    # Escalar datos (muy importante en PCA)
    X = StandardScaler().fit_transform(counts)

    # PCA
    pca = PCA(n_components=2)
    areas_pca = pca.fit_transform(X)

    # Colores por timepoint
    color_map = {
        "1": "#1f77b4",  # basal (embarazo)
        "2": "#ff7f0e",  # post-vacuna
        "3": "#2ca02c",  # post-parto
    }

    sample_colors = metadata["Timepoint"].astype(str).map(color_map)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(
        areas_pca[:, 0],
        areas_pca[:, 1],
        c=sample_colors,
        edgecolor="k",
        s=350,
        alpha=0.9,
    )

    # LaTeX-style axis labels (matching volcano)
    ax.set_xlabel(
        rf"$\mathrm{{PC1}}\ (\mathrm{{{pca.explained_variance_ratio_[0] * 100:.1f}\%}})$",
        fontsize=40,
    )
    ax.set_ylabel(
        rf"$\mathrm{{PC2}}\ (\mathrm{{{pca.explained_variance_ratio_[1] * 100:.1f}\%}})$",
        fontsize=40,
    )
            # Remove ticks for consistency with volcano
    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False,
    )

    # Leyenda manual
    for tp, color in color_map.items():
        ax.scatter([], [], c=color, label=f"TP {tp}", s=150)

    ax.legend(title=f"TP - {condition}", title_fontsize=15,fontsize=15)

    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"PCA_TP_{condition}.png"), dpi=600)
    plt.close()

if __name__ == "__main__":
    script_dir = Path(__file__).parent
    counts, metadata = _load_data(
        (script_dir / "data" / "vsd_plac_matrix.csv"),
        (script_dir / "data" / "metadata.csv")
    )

    output_dir = script_dir / "output"
    print("Counts samples:", counts.index[:5])
    print("Metadata samples:", metadata.index[:5])

    print("Common samples:", len(set(counts.index) & set(metadata.index)))  
    
    counts, metadata = _filter_data(counts, metadata)

    #samples_to_remove = metadata[metadata["Ext_ID"] == "A4316"].index
    #counts = counts.drop(index=samples_to_remove)
    #metadata = metadata.drop(index=samples_to_remove)
    _plot_pca(counts, metadata, output_dir)
    



