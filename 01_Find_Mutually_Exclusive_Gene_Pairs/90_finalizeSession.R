
## @knitr final_saveSessionImage

################################################################################
#'                              1. FINALIZE SESSION           
#'                
#' This script aims to save the reference single-cell dataset with cell type 
#' annotations added to its metadata. It also saves the results from the 
#' FindMarkers step, as well as the mutually exclusive gene pairs obtained from 
#' the analysis. These gene pairs are exported to a .xlsx file, along with an 
#' indication of the cell types to which they correspond.       
#'        
################################################################################

# Create output directory
if (!(dir.exists(PATH_ANALYSIS_OUTPUT))){
  dir.create(PATH_ANALYSIS_OUTPUT)
}

# Save reference single cell ATLAS data as RDS
SaveSeuratRds(reference_sc_data, file = file.path(PATH_ANALYSIS_OUTPUT, 
                                                  "reference_single_cell_data.rds"))

# Save results from FindMarkers
saveRDS(cell_type_lvl1_markers, file = file.path(PATH_ANALYSIS_OUTPUT,
                                                 "cell_type_lvl1_markers.rds"))

# Save mutually exclusive gene pairs as xlsx
write.xlsx(markers_as_strings, file.path(PATH_ANALYSIS_OUTPUT, 
                                         "mutually_exclusive_gene_pairs.xlsx"), 
           rowNames = TRUE)
