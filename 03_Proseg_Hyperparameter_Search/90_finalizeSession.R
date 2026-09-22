
## @knitr final_saveSessionImage

################################################################################
#'                               90. FINALIZE SESSION
#'                               
#' This script aims to save the table containing all quality statistics computed 
#' for each Proseg run. It also saves the MECR scores calculated for each run.
#'                           
################################################################################

# Create output directory
if (!(dir.exists(PATH_ANALYSIS_OUTPUT))){
  dir.create(PATH_ANALYSIS_OUTPUT)
}

###################################################
# SAVE METRICS FROM HYPERPARAMETER SEARCH AS XLSX #
###################################################

# Create workbook
wb = createWorkbook()

# Create sheet for each TMA
for (TMA in TMA_COVERED) {
  addWorksheet(wb, TMA)
  writeData(wb, TMA, TMA_statistics[[TMA]])
}

# Save file
saveWorkbook(wb, 
             file.path(PATH_ANALYSIS_OUTPUT, "hyperparameter_search.xlsx"), 
             overwrite = TRUE)

###########################
# SAVE MECR SCORES AS RDS #
###########################

# Save MECR scores
saveRDS(MECR_scores, 
        file.path(PATH_ANALYSIS_OUTPUT, 
                  "MECR_scores_hyperparameter_search.rds"))
