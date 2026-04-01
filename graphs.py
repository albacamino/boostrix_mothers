#!/usr/bin/env python3

import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def _load_data(path_to_counts_csv, path_to_metadata_csv, path_to_results_csv):
    # Load data
    counts = pd.read_csv(path_to_counts_csv, index_col=0)
    metadata = pd.read_csv(path_to_metadata_csv, index_col=0)
    results = pd.read_csv(path_to_results_csv, index_col=6)

    counts = counts.T

    return counts, metadata, results

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
    # X = StandardScaler().fit_transform(counts)

    # PCA
    pca = PCA(n_components=2)
    areas_pca = pca.fit_transform(counts)

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

def _plot_heatmap(counts, metadata, results, output_dir, timepoints_to_plot=["1","3"]):
    # Filtrar muestras comunes
    common_samples = counts.index.intersection(metadata.index)
    counts_filtered = counts.loc[common_samples]
    metadata_filtered = metadata.loc[common_samples]

    # Filtrar por timepoints que queremos comparar
    metadata_filtered = metadata_filtered[metadata_filtered["Timepoint"].astype(str).isin(timepoints_to_plot)]
    counts_filtered = counts_filtered.loc[metadata_filtered.index]

    condition = metadata_filtered["Condition"].unique()[0]

    # Seleccionar genes/proteínas diferencialmente expresados
    diff_proteins = results[
        (results["padj"] < 0.01) & (results["log2FoldChange"].abs() > 1.5)
    ].index

    # Filtrar la matriz de counts
    diff_counts = counts_filtered.loc[:, counts_filtered.columns.intersection(diff_proteins)]

    # Colores por timepoint
    color_map = {
        "1": "#1f77b4",  # basal (embarazo)
        "2": "#ff7f0e",  # post-vacuna
        "3": "#2ca02c",  # post-parto
    }
    sample_colors = metadata_filtered["Timepoint"].astype(str).map(color_map)

    # Generar heatmap
    g = sns.clustermap(
        diff_counts.transpose(),
        cmap="viridis",
        vmin=-2,
        vmax=2,
        center=0,
        linewidths=0.5,
        linecolor="black",
        yticklabels=False,
        xticklabels=False,
        col_colors=sample_colors,
        row_cluster=True,
        col_cluster=True,
        figsize=(8, 10),
        z_score=1,
        dendrogram_ratio=(0.1, 0.05),
        colors_ratio=(0.02, 0.02),
    )
    g.figure.subplots_adjust(top=1, bottom=0, left=0, right=1)

    # Ajustes visuales
    ax = g.ax_heatmap
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
    ax.set_xlabel("")
    ax.set_ylabel("")

    if g.cax is not None:
        g.cax.remove()
    for sub_ax in [g.ax_col_dendrogram, g.ax_row_dendrogram, g.ax_heatmap]:
        legend = sub_ax.get_legend()
        if legend is not None:
            legend.remove()
        sub_ax.set_title("")

    # Guardar
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(
        os.path.join(output_dir, f"heatmap_{condition}_{'_'.join(timepoints_to_plot)}.png"),
        dpi=600,
        bbox_inches="tight",
    )
    plt.close()

if __name__ == "__main__":
    script_dir = Path(__file__).parent
    counts, metadata, results = _load_data(
        (script_dir / "data" / "vsd_vac_matrix.csv"),
        (script_dir / "data" / "metadata.csv"),
        (script_dir / "data" / "results_Vaccinated.csv")
    )

    output_dir = script_dir / "output"
    counts, metadata = _filter_data(counts, metadata)

    #samples_to_remove = metadata[metadata["Ext_ID"] == "A4316"].index
    #counts = counts.drop(index=samples_to_remove)
    #metadata = metadata.drop(index=samples_to_remove)
    _plot_pca(counts, metadata, output_dir)
    _plot_heatmap(counts, metadata, results, output_dir)
    



