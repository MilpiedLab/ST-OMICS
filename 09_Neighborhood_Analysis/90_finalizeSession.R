
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' Save the merged Seurat object containing the BANKSY assay, the BANKSY 
#' neighborhood annotations across all tested resolutions, and the unified 
#' hierarchical clusters derived from the Jensen–Shannon pairwise comparison of 
#' BANKSY neighborhoods.
#'                     
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Save RDS object
SaveSeuratRds(banksy_xenium_data, file = file.path(PATH_ANALYSIS_OUTPUT, 
                                                   "banksy_merged_TMAs.Rds"))
