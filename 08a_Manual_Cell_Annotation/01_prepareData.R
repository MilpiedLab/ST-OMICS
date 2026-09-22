
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to manually annotate the cell clusters 
obtained from the Harmony-integrated assay by assigning them to their corresponding 
cell types. This annotation process is iterative.

First, cells are clustered and annotated into broad cell types (Level 1). Once 
these major populations are identified, each can be subsetted and re-clustered 
to reveal finer heterogeneity. Based on this refined clustering, more specific 
cell types (Level 2) can be annotated.

This process can be repeated iteratively to define increasingly granular levels 
of annotation (e.g., Level 3 and beyond), until no additional biologically 
meaningful subpopulations can be identified.

An important aspect of this process is the choice of clustering resolution. The 
selected resolution strongly influences the granularity and biological interpretability 
of the resulting clusters. The optimal resolution is therefore chosen based on 
the one that yields the most biologically meaningful clusters, as assessed by marker 
gene expression patterns and their coherence with known cell type signatures.

Manual annotation is guided by the marker genes identified for each cluster 
using the FindAllMarkers function from Seurat. In this case, because a reference 
annotated single-cell dataset is available, correspondence between marker genes 
from the spatial dataset and clusters in the single-cell reference can be used 
to facilitate and refine the annotation process. This comparison helps to identify 
matching cell types and improves annotation accuracy.")
cat("\n\n")

################################################################################
#'                              1. LOAD MERGED DATA
#' 
#' Load the merged TMA dataset containing the Seurat object with the 
#' Harmony-integrated assay obtained in step '07_Integrate_TMAs'.
#' 
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Read merged TMA data
xenium_data = readRDS(PATH_MERGED_TMA_DATA)

# Load cell type hierarchy
reference_cell_type_hierarchy = read.csv(PATH_REFERENCE_CELL_TYPE_HIERARCHY)

# Load reference single cell data for annotation
sc_reference_data = readRDS(file = REFERENCE_SC_DATA_PATH)

################################################################################
#'                          2. DISPLAY CELL TYPE HIERARCHY
################################################################################

cat("\n\n")
cat(sprintf("%s **Cell types hierarchy**", strrep("#", 1)))
cat("\n\n")

cat("\n\n")
cat("Hierarchical tree of cell types based on the annotations from the reference 
single-cell dataset.")
cat("\n\n")

cat("\n\n")
collapsibleTree(reference_cell_type_hierarchy,
                hierarchy = c("LEVEL1", "LEVEL2", "LEVEL3", "LEVEL4"))
cat("\n\n")
