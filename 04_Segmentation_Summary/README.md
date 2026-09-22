# **04 - SEGMENTATION SUMMARY**

## **Overview and objectives**

The aim of this analysis step is to provide a global overview of the Proseg segmentation using the optimal set of hyperparameters selected in the previous step (03_Proseg_Hyperparameter_Search). No specific downstream analysis is performed at this stage; instead, the goal is to visualize the improvement in segmentation quality and obtain a first view of the cell clusters present in the spatial data.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load Xenium data with refined segmentation

Since the optimal set of Proseg hyperparameters was selected in the previous step (03_Proseg_Hyperparameter_Search), the spatial transcriptomic dataset with the refined segmentation corresponding to these parameters can now be loaded as a Seurat object for downstream analysis.

In parallel, the initial dataset with the segmentation provided by Xenium is also loaded as a Seurat object, enabling direct comparison with the Proseg-refined data.

### 2. Plot TMA view

In this case, the spatial transcriptomics data correspond to tissue sections arranged on a tissue microarray (TMA). Four distinct patients are represented, each in duplicate and distributed across two different TMAs, resulting in a total of eight tissue sections.

This section provides a global overview of these two TMAs.

## **02_quality_control.R**

This step consists of performing quality control on cells from both the initial Xenium segmentation and the Proseg-refined segmentation, in order to highlight differences between the two conditions.

### 1. Quality control

For now, quality control is based on two criteria: the number of features per cell and cell area.

  - Number of conserved cells (**nb_conserved**): Conserved cells are those that pass quality control (QC) criteria based on feature counts (CELL_MINIMUM_FEATURES) and cell size (CELL_MAXIMUM_AREA).
  - Number of low feature cells (**nb_low_feature**): Low feature cells are cells that fail the feature count filter but pass the cell size filter.
  - Number of high area cells (**nb_high_area**): High area cells are cells that pass the feature count filter but fail the cell size filter.
  - Number of low feature and high area cells (**nb_low_features_high_area**): Cells that fail both quality thresholds.

### 2. Plot quality control

This section generates all figures related to quality control. For each TMA and each segmentation method, the following visualizations are produced:

  - Violin plots: Distributions of the main cell-level statistics.
  - Scatter plots: Pairwise comparisons of key cell statistics.
  - Cell size distributions: Ridge plots showing the distribution of cell area across the different cell quality categories defined above.
  - QC images: TMA views of the tissue sections, with cells colored according to their quality category.
  - Top gene distribution: Visualization of the most frequently detected genes in the dataset and their relative proportions within the library.

A 'Histograms' section is also included, providing a direct comparison of the number of cells in each quality category across the TMAs and the segmentation methods under study.

### 3. Plot distributions of data statistics

This section generates ridge plots that compare the distributions of key cell-level statistics across each TMA and segmentation method under study.

## **03_compute_MECR.R**

This step consists of computing the MECR score for each gene pair across all TMAs and segmentation methods under study. This enables comparison of MECR values and highlights the improvement achieved by the Proseg-refined segmentation relative to the initial Xenium segmentation.

### 1. Retrieve MECR pairs

Retrieve mutually exclusive gene pairs computed in 01_Find_Mutually_Exclusive_Gene_Pairs.

### 2. Compute MECR

MECR scores for each gene pair are computed across all TMAs and segmentation methods under study. The following visualizations are then generated to enable comparison between conditions:

  - MECR score distributions: Distribution of MECR values for each mutually exclusive gene pair, stratified by segmentation method.
  - MECR scores by condition: A table reporting both the MECR score and the expression rate (i.e., the proportion of cells expressing at least one gene in the mutually exclusive pair) for each gene pair, across all TMAs and segmentation methods under study.

### 3. Compare MECR in segmentation methods

Display a scatter plot comparing the MECR scores obtained with each segmentation method for every mutually exclusive gene pair.

## **04_cluster_TMA_cells.R**

This step aims to visualize cell clusters in the data prior to applying any quality control filters. It allows for a quick visual comparison between the refined segmentation obtained with Proseg and the initial Xenium segmentation.

### 1. Cluster TMAs cells

For each TMA and segmentation method, cluster cells using predefined dimensionality and clustering resolution parameters.

### 2. Plot TMAs clusters

For each TMA, segmentation method, and clustering resolution:

  - DimPlot: UMAP visualization of clustered cells.
  - ElbowPlot: Visualization used to determine the optimal number of dimensions.
  - TMA: Visualization of clusters across the tissue (global TMA view).

### 3. Plot QC on clusters

For each TMA and segmentation method, display cell quality control results on the UMAP. This allows assessment of whether low-quality cells cluster together or are distributed across clusters.

## **05_marker_genes.R**

This step aims to identify which markers are associated with the clusters obtained in the previous step. It is performed prior to any quality control filtering; therefore, the results and cluster assignments may change in later steps after quality control is applied.

### 1. AddModuleScore on UMAP

Since a reference single-cell dataset is available for this spatial dataset, AddModuleScore can be used to evaluate the enrichment of gene signatures for each broad cell type found in the reference. These signatures are based on mutually exclusive gene sets that were computed in '01_Find_Mutually_Exclusive_Gene_Pairs'.

### 2. FindAllMarkers

For each TMA, segmentation method, and clustering resolution, FindAllMarkers is used to identify marker genes for each cluster.

  - DimPlot: UMAP visualization of clustered cells.
  - Cluster index:
    * Features: UMAP-based visualization of the top 9 marker genes per cluster, selected based on log fold-change expression.
    * Violin plots: Display the expression distribution of the top 9 marker genes across clusters.

## **90_finalizeSession.R**

Save the Seurat objects corresponding to the Proseg segmentation for each TMA under study. Note that these objects will not be reused in downstream analyses, as the clustering performed here was intended solely for data visualization and was conducted without applying any quality control filtering.
