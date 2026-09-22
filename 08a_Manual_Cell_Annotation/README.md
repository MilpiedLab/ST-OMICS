# **08a - MANUAL CELL ANNOTATION**

## **Overview and objectives**

The aim of this analysis step is to manually annotate the cell clusters obtained from the Harmony-integrated assay by assigning them to their corresponding cell types. This annotation process is iterative.

First, cells are clustered and annotated into broad cell types (Level 1). Once these major populations are identified, each can be subsetted and re-clustered to reveal finer heterogeneity. Based on this refined clustering, more specific cell types (Level 2) can be annotated. Note that for the re-clustering step, a new Harmony integration must be performed on this subset of cells prior to clustering.

An important aspect of this process is the choice of clustering resolution. The selected resolution strongly influences the granularity and biological interpretability of the resulting clusters. The optimal resolution is therefore chosen based on the one that yields the most biologically meaningful clusters, as assessed by marker gene expression patterns and their coherence with known cell type signatures.

This process can be repeated iteratively to define increasingly granular levels of annotation (e.g., Level 3 and beyond), until no additional biologically meaningful subpopulations can be identified.

Manual annotation is guided by the marker genes identified for each cluster using the FindAllMarkers function from Seurat. In this case, because a reference annotated single-cell dataset is available, correspondence between marker genes from the spatial dataset and clusters in the single-cell reference can be used to facilitate and refine the annotation process. This comparison helps to identify matching cell types and improves annotation accuracy.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load merged TMAs data

Load the merged TMA dataset containing the Seurat object with the Harmony-integrated assay obtained in step '07_Integrate_TMAs'.

### 2. Display cell type hierarchy

Display a hierarchical tree of cell types based on the annotations from the reference single-cell dataset.

## **02_cell_annotation_lvl1.R**

This step aims to annotate clusters obtained from the full Harmony-integrated dataset. The annotated clusters will correspond to the Level 1 cell types described previously.

### 1. Clustering of harmony-integrated dataset 

Cluster cells using the Harmony-integrated assay at each tested resolution. 

### 2. Plot integrated merged TMA data

For each clustering resolution tested, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - Clusters frequencies: Cell frequencies per cluster.
  - UMAP by sample ID: UMAP colored by sample ID.
  - Individual clusters: UMAP visualization highlighting each cluster at the current resolution.

### 3. Characterize lvl1 clusters

This step aims to identify marker genes and corresponding clusters in the reference single-cell dataset that can help annotate clusters in the merged Xenium data. To do so, marker genes are computed for each cluster. The top marker genes are then used as a gene module to assess similarity with cells in the reference single-cell dataset using the AddModuleScore function from Seurat.

For each clustering resolution tested on the Harmony-integrated assay, and for each cluster within these resolutions, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Feature plots: UMAP-based visualization of the top 9 marker genes per cluster, selected based on log fold-change expression.
  - Violin plots: Expression distribution of the top 9 marker genes across clusters.
  - Reference UMAP: Reference single-cell UMAP visualization, with cells colored according to the module score derived from the marker genes of the Harmony cluster under study.
  - Reference AddModuleScore table: Table reporting the mean module score per cell type in the reference single-cell dataset, computed using the marker-derived gene module of the Harmony cluster under study.

### 4. Display available manual annotation

In this case, two of the clustering resolutions tested earlier were used for manual annotation. This is done to illustrate how the choice of resolution influences the resulting annotation and the biological interpretation of the data.

For each resolution for which manual annotation is provided, the following sections are displayed:

  - DimPlot: UMAP visualization of cells colored by the corresponding manual annotation.
  - Cell type distribution: Cell frequencies for each manually identified cell type.

### 5. Final lvl1 annotation

Add a metadata column to the Seurat object corresponding to the Level 1 manual annotation selected as the most appropriate based on the chosen clustering resolution.

### 6. Display final lvl1 annotation

Display the Level 1 manual annotation selected using the clustering resolution deemed most appropriate, which will serve as the basis for proceeding to Level 2 annotation. The following sections are displayed:

  - DimPlot: UMAP visualization of cells colored by the corresponding manual annotation.
  - Cell type distribution: Cell frequencies for each manually identified cell type.

## **03_process_segmentation_errors.R**

In the previous step, some cells were annotated as segmentation errors because they exhibited marker genes from multiple distinct cell types, suggesting potential issues such as poor segmentation.

The aim of this analysis is to re-examine this subset of cells in order to recover as many valid cells as possible and reduce the number of cells classified as segmentation errors. To achieve this, cells annotated as segmentation errors at Level 1 are subsetted and re-clustered.

The goal of this finer clustering is to identify and recover genuine biological cell populations that were previously misclassified due to higher levels of segmentation uncertainty. Recovered cells can then be reassigned to the appropriate Level 1 manual annotation clusters whenever supported by their marker gene expression profiles.

### 1. Characterize segmentation errors

This step aims to characterize the segmentation error cluster identified during Level 1 manual annotation. The following sections are displayed:

  - Markers: Marker genes identified for the segmentation error cluster using the FindAllMarkers function from Seurat.
  - Feature plots: UMAP-based visualization of marker gene expression within the segmentation error cluster, highlighting the non-specific and mixed expression patterns characteristic of this population.

### 2. Segmentation error cells clustering

As explained earlier, this step aims to re-cluster cells initially identified as segmentation errors in order to recover as many valid cells as possible. This section provides a more detailed view of the re-clustering process, and the following outputs are displayed:

  - Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding.
  - Re-clustering resolutions tested:
    * DimPlot: UMAP visualization of re-clustered cells.
    * TMA ID batch effect: UMAP visualization colored by TMA ID to assess whether integration has been successful and batch effects have been mitigated.
    * Bar plot: Cell proportions per cluster.
  - Annotation resolution: UMAP visualization of re-clustered cells at the resolution selected for manual annotation.
  - Annotations: UMAP visualization of manual annotations assigned to the re-clustered segmentation error cells.
  - Cell type distribution: Frequencies of annotated cell types within the re-clustered segmentation error population.

### 3. Merge lvl1 and segmentation errors annotations

This step aims to update the full dataset by incorporating the revised manual annotations for a subset of cells previously classified as segmentation errors. Cells that could be confidently re-annotated are reassigned to the appropriate Level 1 cell types, thereby refining the overall annotation and reducing the number of cells labeled as segmentation errors.

The following sections are displayed to visualize the updated annotations:

  - Annotations before correction: Level 1 manual annotations prior to reassigning cells previously labeled as segmentation errors.
  - Annotations after correction: Level 1 manual annotations after incorporating the updated annotations for the corrected cells.
  - Sample IDs: Visualization of each level 1 cell type on the sample tissues.

## **04_cell_annotation_lvl2.R**

Now that the Level 1 manual annotation has been completed, this step focuses on subsetting each Level 1 cell type and re-clustering them individually, as described previously. This refined clustering approach enables the identification of finer heterogeneity within closely related cell populations.

The objective of this step is therefore to define and manually annotate Level 2 cell types within the dataset.

Note that some closely related Level 1 cell types were pooled prior to re-clustering.

### 1. Characterize lvl2 clusters

This step aims to visualize the re-clustering performed for each Level 1 cell type or group of related Level 1 cell types.

For each subset, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - For each cluster within the subset:
    * DimPlot: UMAP visualization highlighting the cells of the current cluster.
    * Feature plots: UMAP-based visualization of the top 9 marker genes for the current cluster, selected based on log fold-change expression.
    * Violin plots: Expression distribution of the top 9 marker genes across clusters.
    * Reference UMAP: UMAP visualization of the reference single-cell dataset, with cells colored according to the module score derived from the marker genes of the cluster under study. Only reference cell types relevant to the current Level 1 subset are included in the module score computation.
    * Reference AddModuleScore table: Table reporting the mean module score per cell type in the reference single-cell dataset, computed using the marker-derived gene module of the cluster under study (via the AddModuleScore function in Seurat).
  - Subset Level 2 final annotation: Final Level 2 manual annotation assigned to the current Level 1 subset.

### 2. Final lvl2 annotation

Add a metadata column to the Seurat object corresponding to the Level 2 manual annotation.

### 3. Display final lvl2 annotation

Display the Level 2 manual annotation. The following sections are provided:

  - DimPlot: UMAP visualization of cells colored according to the corresponding Level 2 manual annotation.
  - Cell type distribution: Cell proportions for each manually identified Level 2 cell type.
  - Sample IDs: Visualization of each Level 2 cell type on the sample tissues.
  
## **90_finalizeSession.R**

Save the Seurat object containing both the Level 1 and Level 2 manual annotations.
  