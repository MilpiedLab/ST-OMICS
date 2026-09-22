
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to provide a global overview of the Proseg 
segmentation using the optimal set of hyperparameters selected in the previous 
step (03_Proseg_Hyperparameter_Search). No specific downstream analysis is 
performed at this stage; instead, the goal is to visualize the improvement in 
segmentation quality and obtain a first view of the cell clusters present in the 
spatial data.")
cat("\n\n")

################################################################################
#'                  1. LOAD XENIUM DATA WITH REFINED SEGMENTATION
#' 
#' Since the optimal set of Proseg hyperparameters was selected in the previous 
#' step (03_Proseg_Hyperparameter_Search), the spatial transcriptomic dataset 
#' with the refined segmentation corresponding to these parameters can now be 
#' loaded as a Seurat object for downstream analysis.
#' 
#' In parallel, the initial dataset with the segmentation provided by Xenium is 
#' also loaded as a Seurat object, enabling direct comparison with the 
#' Proseg-refined data.
#' 
################################################################################

# Load Xenium data in list
xenium_data = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Load Xenium data
    xenium_data[[TMA]][[segmentation_method]] = LoadXenium_new(TMA_METADATA[[TMA]][[segmentation_method]][["data_dir"]])
    
  }
}

################################################################################
#'                                2. PLOT TMA VIEW
################################################################################

cat("\n\n")
cat(sprintf("%s **TMA global view**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("In this case, the spatial transcriptomics data correspond to tissue sections 
arranged on a tissue microarray (TMA). Four distinct patients are represented, 
each in duplicate and distributed across two different TMAs, resulting in a total 
of eight tissue sections.

This section provides a global overview of these two TMAs.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA))
  cat("\n\n")
  
  # Display captured TMA section
  cat("\n\n")
  print(ImageDimPlot(xenium_data[[TMA]][["xenium_segmentation"]]))
  cat("\n\n")
}
