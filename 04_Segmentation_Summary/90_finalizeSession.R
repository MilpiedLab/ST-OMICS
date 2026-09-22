
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' Save the Seurat objects corresponding to the Proseg segmentation for each TMA 
#' under study. Note that these objects will not be reused in downstream analyses, 
#' as the clustering performed here was intended solely for data visualization 
#' and was conducted without applying any quality control filtering.
#'                           
################################################################################

# Create output directory
if (!(dir.exists(PATH_ANALYSIS_OUTPUT))){
  dir.create(PATH_ANALYSIS_OUTPUT)
}

#####################################
## SAVE PROSEG SEGMENTATION AS RDS ##
#####################################

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Create output directory
  out_dir = file.path(TMA_METADATA[[TMA]][["proseg_segmentation"]][["output_path"]])
  dir.create(out_dir, 
             recursive = TRUE, 
             showWarnings = FALSE)
  
  # Save RDS
  SaveSeuratRds(xenium_data[[TMA]][[names(TMA_METADATA[[TMA]])[2]]],
                file = file.path(out_dir, 
                                 paste0(names(TMA_METADATA[[TMA]])[2], ".Rds")))
  
}
