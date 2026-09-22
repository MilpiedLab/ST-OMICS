
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to characterize cell neighborhoods in each 
sample of the merged TMA dataset. In spatial transcriptomics, neighborhood analysis 
consists of clustering cells by jointly leveraging their molecular profiles and 
spatial context. In this framework, neighborhoods are defined as groups of cells 
sharing a similar local microenvironment.

Specifically, this is achieved using BANKSY, which constructs a neighbor-augmented 
expression matrix by integrating each cell’s own transcriptome with the transcriptomic 
information of its spatial neighbors. Standard clustering methods are then applied 
to this augmented representation, similarly to single-cell RNA-seq analysis. 
However, in this context, the resulting clusters do not represent cell types, 
but instead correspond to spatially defined cellular neighborhoods.

Once these neighborhoods are defined, their cellular composition can be characterized 
in terms of cell types, and correspondence between similar neighborhood structures 
across different samples can be assessed.")
cat("\n\n")

################################################################################
#'                            1. LOAD ANNOTATED DATA
#' 
#' Load the merged TMAs Seurat object containing manual cell type annotations.
#' 
################################################################################

# Load annotated Xenium data
xenium_data = readRDS(file = PATH_ANNOTATED_MERGED_TMA_DATA)
