
## @knitr quality_control

################################################################################
#'                                1. QUALITY CONTROL
#'                                
#' This step performs quality control on the merged dataset. Although low-quality 
#' cells were previously identified within individual TMAs, they have not yet 
#' been removed. Here, filtering is applied so that only samples and cells labeled 
#' as conserved are retained in the dataset.
#' 
#' Therefore, assign a quality label to each cell based on the defined quality 
#' control criteria.        
#'                        
################################################################################

# Tag cells that should be removed based on filtering parameters for current condition
complete_xenium_data = pastis_quality_control(complete_xenium_data,
                                              minimum_features = CELL_MINIMUM_FEATURES,
                                              maximum_area = CELL_MAXIMUM_AREA)

################################################################################
#'                      2. FILTER LOW QUALITY SAMPLES AND CELLS
#'                      
#' Apply quality filters.
#'                     
################################################################################

# Keep cells tagged as "conserved"
filtered_complete_xenium_data = subset(complete_xenium_data,
                                       subset = QC_result == "conserved")

# Keep only high quality samples
filtered_complete_xenium_data = subset(filtered_complete_xenium_data,
                                       subset = sample_tissue_quality == "conserved")
