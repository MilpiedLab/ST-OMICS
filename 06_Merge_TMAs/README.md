# **06 - MERGE TMAs**

## **Overview and objectives**

The aim of this analysis step is to merge the Seurat objects generated in the '05_TMA_Processing' step where cells have been assigned to their corresponding samples and annotated with quality control metrics into a single Seurat object.

Once merged, both samples and cells are filtered based on the previously defined quality control criteria.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load Seurat objects

Load the Seurat objects for each TMA processed in the 05_TMA_Processing step, which now have metadata assinging each of their cells to their corresponding sample.

## **02_merge_data.R**

### 1. Merge data

Merge every available TMA into a single Seurat object.

### 2. Save merged data

Save the merged Seurat object before applying any filtering steps.

## **03_quality_control.R**

This step performs quality control on the merged dataset. Although low-quality cells were previously identified within individual TMAs, they have not yet been removed. Here, filtering is applied so that only samples and cells labeled as conserved are retained in the dataset.

### 1. Quality control

Assign a quality label to each cell based on the defined quality control criteria.

### 2. Filter low quality samples and cells

Apply quality filters.

## **04_cluster_unfiltered_merged_data.R**

This step performs cell clustering on the fully merged dataset prior to any quality filtering of samples or cells. It is intended to illustrate the impact of low-quality cells and samples, emphasizing the importance of applying quality control early in the analysis.

### 1. Clustering

Cluster cells using predefined dimensionality and clustering resolution parameters.

### 2. Visualize unfiltered merged data

Visualize the data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by QC: UMAP colored by cell quality, highlighting that low-feature cells tend to cluster together.
  - Clusters by sample quality: UMAP colored by sample quality, showing that low-quality samples cluster together. These samples did not pass the quality filters based on the transcript density analysis from the '05_TMA_Processing' step. They actually correspond to a different tissue type than the other samples, which explains the poor RNA capture by the gene panel designed for this Xenium experiment.
  - Clusters by sample ID: UMAP colored by sample ID, consistent with the clustering observed when colored by sample quality.
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect associated with the TMA.

### 3. Save clustered unfiltered merged data

Save data processed in current section.

## **05_cluster_filtered_merged_data.R**

This step performs cell clustering on the fully merged dataset after quality filtering. Both low-quality cells and low-quality samples, as previously identified, have been removed from the data.

### 1. Clustering

Cluster cells using predefined dimensionality and clustering resolution parameters.

### 2. Visualize filtered merged data

Visualize the data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells for every clustering resolution covered.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by sample ID: UMAP colored by sample ID
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect associated with the TMA.

### 3. Save clustered filtered merged data

Save data processed in current section.

## **06_compute_filtered_merged_data_markers.R**

This step computes gene markers for the clusters obtained earlier from the filtered, merged dataset. This allows identification of the broad cell types present in the data and evaluation of how well they can be distinguished from one another, providing insight into the quality of the segmentation.

For each clustering resolution and corresponding clusters, the following figures are available:
    
  - DimPlot: UMAP visualization of current cluster of interest.
  - Features: UMAP-based visualization of the top 9 marker genes per cluster, selected based on log fold-change expression.
  - Violins: Display the expression distribution of the top 9 marker genes across clusters.
  
Marker genes obtained in this section are saved as a '.xlsx' file.
