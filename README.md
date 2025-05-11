# scRNA-seq analysis in murine prostate cancer samples
This repository contains a structured and reproducible pipeline for analyzing single-cell RNA-seq data in the context of prostate cancer progression. The scripts are designed to handle multiple biological conditions and generate publication-ready outputs.

## Repository structure

### `01_scRNAseq_processing_Harmony.R`

Processes filtered gene-barcode matrices from Cell Ranger using Seurat and Harmony. Generates cluster annotations and saves intermediate objects.

### `02_GO_analysis.R`

Performs Gene Ontology (GO) enrichment analysis using the output of differential expression. Functions allow input from either Seurat objects or Excel files.

### `03_Pseudotime_TSCAN.R`

Subsamples balanced cell numbers per cluster and computes pseudotime trajectories using TSCAN. Includes gene expression overlay and trajectory visualizations.

### `04_CellChat_createObjects.R`

Generates individual CellChat objects for each sample and processes them using the mouse CellChatDB. Outputs are saved as `.rds`.

### `05_CellChat_comparison.Rmd`

Compares ligand-receptor interaction networks across conditions using chord diagrams, heatmaps, and differential expression analysis of signaling pathways.

### `06_Visualization_Seurat.R`

Generates common single-cell visualizations such as UMAPs, violin plots, dot plots, heatmaps, and stacked bar charts. Based on Seurat outputs.

> **Note**: This script is optional and intended for figure generation.

## Data organization

* All scripts assume a consistent metadata column named `sample` for biological condition.
* Cluster annotations are stored in the `sub_` or `clusters` columns.

## Reproducibility

* Random seeds are set in scripts involving subsampling (e.g., TSCAN pseudotime) to ensure reproducibility.
* All scripts are modular and independently executable.
* Figures are generated using `ggplot2`, `Seurat`, `plotly`, or `CellChat`.

## Requirements

All required packages are loaded at the beginning of each script. A list of versions for all packages used in the analysis can be found in package_versions.txt.

### Input examples

The `data/` folder contains example files used across the pipeline:

* `example_object.rds`: a small Seurat object subset for testing.
* `RMP_annotated.xlsx`: gene list used for dot plots.

### Output files

The pipeline generates the following types of outputs:

* `.rds`: Seurat and CellChat objects
* `.xlsx`: marker and dotplot tables
* `.pdf` / `.html`: visualizations from CellChat

Author: Ana Añazco
Institution: CIC (University of Salamanca - CSIC)
Year: 2025

