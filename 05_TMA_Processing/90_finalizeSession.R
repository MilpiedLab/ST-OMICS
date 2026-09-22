
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' Save the Seurat objects corresponding to the optimal Proseg segmentation for 
#' each TMA, with each cell now assigned to its corresponding sample on the TMA. 
#' These objects also include metadata from both sample-level quality assessment 
#' and cell-level quality control. At this stage, no filtering has been applied 
#' to either cells or samples. All filtering steps will be performed downstream 
#' on these Seurat objects.      
#'                     
################################################################################

########################################
## SAVE PROSEG SEGMENTATION AS RDS ##
########################################

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Create output directory
  dir.create(PATH_ANALYSIS_OUTPUT, 
             recursive = TRUE, 
             showWarnings = FALSE)
  
  # Save RDS
  SaveSeuratRds(xenium_data[[TMA]],
                file = file.path(PATH_ANALYSIS_OUTPUT, 
                                 paste0(TMA, "_splitted_", PROSEG_SEGMENTATION_NAME, ".Rds")))
}
