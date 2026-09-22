
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to merge the Seurat objects generated in 
the '05_TMA_Processing' step where cells have been assigned to their corresponding 
samples and annotated with quality control metrics into a single Seurat object.

Once merged, both samples and cells are filtered based on the previously defined 
quality control criteria.")
cat("\n\n")

################################################################################
#'                               1. LOAD SEURAT OBJECTS
#'                               
#' Load the Seurat objects for each TMA processed in the 05_TMA_Processing step, 
#' which now have metadata assinging each of their cells to their corresponding 
#' sample.            
#'                   
################################################################################


# Store data in list
xenium_data = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Read RDS for current condition
  xenium_data[[TMA]] = readRDS(file = PATH_SPLITTED_TMA_RDS[[TMA]])
  
  # Add meta data layer with TMA name
  xenium_data[[TMA]]$TMA_id = TMA
  
}
