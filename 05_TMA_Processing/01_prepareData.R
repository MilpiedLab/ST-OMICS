
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("This script splits the TMA into individual samples. Currently, the Seurat 
objects do not contain metadata linking each cell to its corresponding sample on 
the TMA.

The goal is therefore to identify the spatial location of each sample on the TMA 
and assign cells accordingly.

Sample detection is performed using a nearest-neighbor graph. Cells are first 
connected to their spatial nearest neighbors, and the resulting connected 
components are then classified as either individual samples or background based 
on their size (i.e., number of cells per component). See the split_xenium_tma_samples 
function for more details.")
cat("\n\n")

################################################################################
#'                              1. LOAD SEURAT OBJECTS
#' 
#' This step loads the spatial transcriptomics data with Proseg-refined 
#' segmentation for each TMA under study. For each TMA, it uses the optimal 
#' Proseg parameters identified in the '03_Proseg_Hyperparameter_Search' step.
#' 
################################################################################

# Initialise list to store the TMAs segmented data
xenium_data = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){

  # Load Xenium data
  xenium_data[[TMA]] = LoadXenium_new(TMA_METADATA[[TMA]][[PROSEG_SEGMENTATION_NAME]][["data_dir"]])
  
}
