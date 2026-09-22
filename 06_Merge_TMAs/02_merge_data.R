
## @knitr merge_data

################################################################################
#'                                 1. MERGE DATA
#'                                 
#' Merge every available TMA into a single Seurat object.
#'                                  
################################################################################

# Merge all available TMA data sets
complete_xenium_data = merge(x = xenium_data[[1]], y = xenium_data[2:length(xenium_data)])
complete_xenium_data = JoinLayers(complete_xenium_data)

################################################################################
#'                               2. SAVE MERGED DATA
#'                               
#' Save the merged Seurat object before applying any filtering steps.
#'                               
################################################################################

## @knitr save_merged_data

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Save merged TMA data
saveRDS(complete_xenium_data,
        file = file.path(PATH_ANALYSIS_OUTPUT,
                         "unfiltered_merged_TMAs.Rds"))
