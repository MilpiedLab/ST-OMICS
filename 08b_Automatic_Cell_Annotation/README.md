# **08b - AUTOMATIC CELL ANNOTATION**

## **Overview and objectives**

The aim of this analysis step is to automatically annotate cell clusters and individual cells obtained from the Harmony-integrated assay using two distinct methods: scType and SingleR, which rely on different principles.

scType performs annotation at the cluster level by evaluating the enrichment of known marker genes from curated cell type–specific gene sets. In contrast, SingleR operates at the single-cell level and assigns labels by comparing each cell’s expression profile to annotated reference cells from a single-cell dataset using correlation-based approaches.

These two methods are applied independently to provide alternative automatic annotations based on different underlying assumptions.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load merged TMAs data

Load the merged TMA dataset containing the Seurat object with the Harmony-integrated assay obtained in step '07_Integrate_TMAs'.

### 2. Display cell type hierarchy

Display a hierarchical tree of cell types based on the annotations from the reference single-cell dataset.

## **02_cluster_merged_data.R**

This step aims to cluster cells from the full Harmony-integrated dataset. The resulting clusters will subsequently be used for annotation with scType.

### 1. Clustering

Cluster cells using the Harmony-integrated assay across the range of tested resolutions.

### 2. Plot integrated merged TMA data

For each clustering resolution tested, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - Cluster frequencies: Cell proportions per cluster.
  - UMAP by sample ID: UMAP colored by sample ID.
  - Individual clusters: UMAP visualization highlighting each cluster at the current resolution.

### 3. Harmony clusters markers

This step aims to identify marker genes associated with the Harmony-derived clusters obtained in the previous step.

For each clustering resolution tested on the Harmony-integrated assay, and for each cluster within these resolutions, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Feature plots: UMAP-based visualization of the top 9 marker genes per cluster, selected based on log fold-change expression.
  - Violin plots: Expression distribution of the top 9 marker genes across clusters.

## **03_scType_annotation.R**

This step aims to automatically annotate the Harmony-derived clusters obtained earlier using scType. For each clustering resolution tested, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - scType annotation: Cluster-level annotations assigned by scType.
  - scType results table: Annotation scores provided by scType for each cluster at the current resolution.
  - Sample IDs: Visualization of each cell type identified by scType on the sample tissues.

## **04_singleR_annotation.R**

This step aims to automatically annotate cells using SingleR. As SingleR operates at the single-cell level, cells within the same cluster may be assigned to different, and sometimes inconsistent, cell types.

To obtain consistent cluster-level annotations, each cluster is labeled based on the most frequent SingleR annotation among its constituent cells (Level 1). In addition, SingleR is run using a reference single-cell dataset restricted to Level 1 cell type annotations, ensuring that the predicted labels correspond to broad cell populations.

Although SingleR can directly provide fine-grained (e.g., Level 2) annotations at the single-cell level, this analysis adopts an alternative strategy. Because annotations are harmonized at the cluster level, clusters are subsequently subsetted and re-clustered to define new clusters capturing finer levels of heterogeneity, following an approach similar to the manual annotation described in Section 08a_Manual_Cell_Annotation.

Based on this refined clustering, more specific cell types (Level 2) can then be assigned. These annotations are obtained using scType, with cluster labels determined by the most frequent predicted cell type within each refined cluster.

### 1. SingleR lvl1 clustering

In this section, a DimPlot is displayed showing the clusters obtained at a high resolution, in order to generate a large number of clusters from the full Harmony-integrated dataset.

### 2. SingleR lvl1 annotation

This step aims to visualize the Level 1 SingleR annotations using the following sections:

  - Raw SingleR Level 1 annotation: UMAP visualization of the raw SingleR cell-level annotations.
  - Cluster annotation refinement: For each cluster identified in the previous section:
    * DimPlot: UMAP visualization highlighting the current cluster.
    * Cell type frequencies: Bar plot showing the distribution of SingleR cell annotations within the current cluster.
  - Refined SingleR Level 1 annotation:
    * Global view: UMAP visualization of the harmonized (cluster-level) SingleR annotations.
    * Sample IDs: Visualization of each singleR Level 1 annotations on the sample tissues.

### 3. SingleR lvl2 annotation

As explained earlier, clusters sharing the same refined Level 1 SingleR annotation are subsetted together and re-clustered to define new clusters capturing finer levels of heterogeneity. These finer clusters are then annotated using SingleR following the same procedure as for Level 1 clusters.

This step aims to visualize the Level 2 SingleR annotations using the following sections:

  - Level 1 annotation (per subset of refined SingleR Level 1 clusters):
    * DimPlot: UMAP visualization of clustered cells.
    * Marker genes: Table reporting the results from FindAllMarkers for all clusters.
    * Elbow plot: Visualization used to determine the number of dimensions to retain from the Harmony-corrected embedding.
    * Raw SingleR Level 2 annotation: UMAP visualization of raw cell-level SingleR annotations.
    * Cluster annotation refinement: For each cluster identified within the current subset:
      - DimPlot: UMAP visualization highlighting the current cluster.
      - Cell type frequencies: Bar plot showing the distribution of SingleR annotations within the cluster.
    * Refined SingleR Level 2 annotation: UMAP visualization of harmonized (cluster-level) SingleR annotations.
    * Sample IDs: Visualization of Level 2 SingleR annotations mapped onto tissue samples for the current subset.
  - Refined global Level 2 SingleR annotation: UMAP visualization of harmonized Level 2 SingleR annotations for the full merged dataset.
  
## **90_finalizeSession.R**

Save the Seurat object containing both scType and singleR annotations.
