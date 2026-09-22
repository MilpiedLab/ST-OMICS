# **01 - FIND MUTUALLY EXCLUSIVE GENE PAIRS**

## **Overview and objectives**

Mutually Exclusive Co-expression Rate (MECR) is a metric that can be used to assess the quality of cell segmentation in imaging spatial transcriptomics data.

MECR is computed by dividing the fraction of cells that express both markers by the fraction of cells that express at least one of the markers, for each pair of mutually exclusive genes (Satija et al., bioRxiv, 2024).

Computing this metric requires access to a list of mutually exclusive gene pairs for the dataset being analyzed. These correspond to genes that are not expected to be co-expressed within the same cell type (i.e., genes whose expression is typically restricted to distinct cell types).

The aim of this analysis step is therefore to retrieve a list of mutually exclusive gene pairs from a reference single-cell dataset. The cell types present in this reference are expected to also be represented in the spatial dataset. The gene pairs identified in this analysis will be used in subsequent steps to assess the quality of cell segmentation.

Docker use: "**STOmics_Seurat5**".

## **01_prepareData.R**

### 1. Load and prepare single-cell dataset

Load the reference single-cell dataset and manually complete the cell type annotations in the corresponding Seurat object metadata.

### 2. Load Xenium spatial dataset

Load the Xenium dataset to be analyzed in this experiment. From the Xenium Seurat object, retrieve the genes of interest (i.e., the gene panel) used to generate the data, in order to restrict the search for mutually exclusive gene pairs to this specific set of genes. Constraining the exclusive gene pairs search to the Xenium gene panel reduces computation time.

### 3. Cell type tree vizualisation

Build a hierarchical tree of cell types based on the annotations in the Seurat object metadata. This hierarchical tree will be displayed in the HTML report corresponding to this analysis step.

## **02_find_mutually_exclusive_genes.R**

### 1. Retrieve mutually exclusive genes

This script aims to compute mutually exclusive gene pairs from single-cell reference data. Here, the FindMarkers function is used alongside predefined thresholds, described below, to identify mutually exclusive genes for each cell type found in the reference dataset. Other methods than the one described below may also be considered to define these gene pairs.

Even though four levels of cell type annotation are available for this dataset, only the first level will be considered to retrieve mutually exclusive gene pairs. While it is feasible to identify exclusive pairs for broad cell types, doing so becomes increasingly difficult as the cell type annotations become more refined.

Once cell type–specific genes are identified using the FindMarkers function, markers are filtered based on log fold change (logFC), pct.1, and pct.2. This approach allows the identification of genes that should be specific to each cell type.

A marker gene is considered cell type–specific if it is expressed in more than 50% of the cells within that cell type and in less than 5% of the other cells. Mutually exclusive gene pairs are then derived in subsequent analysis steps by computing all heterotypic pairs of these exclusive genes.

## **90_finalizeSession.R**

This script aims to save the reference single-cell dataset with cell type annotations added to its metadata. It also saves the results from the FindMarkers step, as well as the mutually exclusive gene pairs obtained from the analysis. These gene pairs are exported to a .xlsx file, along with an indication of the cell types to which they correspond.
