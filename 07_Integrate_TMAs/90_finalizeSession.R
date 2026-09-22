
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' Save the Seurat object containing the Harmony-integrated assay.  
#'                     
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Save RDS
SaveSeuratRds(xenium_data,
              file = file.path(PATH_ANALYSIS_OUTPUT, 
                               "integrated_merged_TMA.Rds"))
