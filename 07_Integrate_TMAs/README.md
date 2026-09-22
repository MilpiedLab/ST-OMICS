# **07 - INTEGRATE TMAs**

## **Overview and objectives**

The aim of this analysis step is to integrate the data within the merged Seurat object in order to remove any batch effects that may be observed. Indeed, in Section 06_Merge_TMAs, a batch effect related to the TMA was identified when clustering cells from the merged dataset.

Therefore, Harmony will be used for the integration. The batch-corrected low-dimensional embedding obtained from Harmony will help reduce batch effects between the TMAs and enable better characterization of cell types during manual annotation in downstream analyses.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load merged data

Load the data after merging and after filtering out low-quality cells and samples, as described in Section '06_Merge_TMAs'.

### 2. Plot merged data before integration

Visualize the merged and filtered data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells for each clustering resolution tested.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by sample ID: UMAP colored by sample ID.
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect associated with the TMA.

## **02_harmony_integration.R**

This step performs data integration using Harmony. The integration is carried out using the TMA ID as a covariate representing the source of batch effects.

### 1. Integration by TMA ID

The following sections are provided to visualize and assess the effect of the integration:

  - Harmony convergence plot: Figure showing the evolution of the objective function across Harmony iterations, used to assess convergence of the algorithm.
  - Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - UMAP by TMA ID: UMAP colored by TMA ID. Visually, the batch effect appears to be substantially reduced; however, quantitative metrics should be included for a more rigorous assessment.
  - UMAP by sample ID: UMAP colored by sample ID.

### 2. Clustering

This step aims to cluster cells based on the Harmony-corrected embedding (i.e., the batch-corrected low-dimensional representation).

### 3. Plot clusters and marker genes

This step aims to visualize the clusters obtained after Harmony integration. For each clustering resolution tested, the following sections are provided:

  - DimPlot: UMAP visualization of clustered cells.
  - Bar plot: Cell frequencies per cluster.
  - Individual clusters: UMAP visualization highlighting each cluster at the current resolution.

For each cluster shown above, marker genes are identified using the FindAllMarkers function from Seurat. The results of the marker analysis will be saved in the next step as a .xlsx file.

### 4. Save marker genes

Save the marker genes identified for the clusters above as a .xlsx file

## **90_finalizeSession.R**

Save the Seurat object containing the Harmony-integrated assay.
