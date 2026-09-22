
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' Save the Seurat object containing both the Level 1 and Level 2 manual 
#' annotations.  
#'                     
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Save RDS
SaveSeuratRds(xenium_data,
              file = file.path(PATH_ANALYSIS_OUTPUT, 
                               "manually_annotated_merged_TMA.Rds"))
