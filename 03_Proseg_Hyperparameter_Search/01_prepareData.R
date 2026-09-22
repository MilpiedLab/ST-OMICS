
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to identify the set of Proseg hyperparameters, 
processed in the previous step (02_Cell_Segmentation), that yields the highest-quality 
cell segmentation. Segmentation quality is evaluated using several metrics, which 
will be detailed below.")
cat("\n\n")

################################################################################
#'                    1. RETRIEVE MUTUALLY EXCLUSIVE GENE PAIRS
################################################################################

cat("\n\n")
cat(sprintf("%s **Mutually exclusive gene pairs**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Retrieve the mutually exclusive genes identified in step 
01_Find_Mutually_Exclusive_Gene_Pairs. These genes are required to compute the 
MECR metric, which is the primary measure used here to assess cell segmentation 
quality.

From this set of genes, all heterotypic pairs (i.e., pairs between different 
cell types) are generated. These pairs correspond to the mutually exclusive 
gene pairs used to evaluate data quality.")
cat("\n\n")

# Load Xenium data as Seurat object
Xenium_data = LoadXenium_new(XENIUM_DATA_PATH[["TMA002"]])

# Retrieve gene panel
gene_panel_480 = rownames(Xenium_data)

# Clean env
rm(Xenium_data)
gc()

# Retrieve mutually exclusive genes obtained in "01_Find_Mutually_Exclusive_Gene_Pairs"
markers_as_strings = as.data.frame(read_excel(MUTUALLY_EXCLUSIVE_GENE_PAIRS))
rownames(markers_as_strings) = markers_as_strings$...1
markers_as_strings$...1 = NULL

# FACULTATIVE: Here, remove cell types that will not be considered for MECR computation
markers_as_strings = markers_as_strings[!rownames(markers_as_strings) %in% c("Granulo", "LAMP3cDC_9"), , drop = FALSE]

cat("\n\n")
cat(sprintf("%s Markers by cell type", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(htmltools::tagList(DT::datatable(markers_as_strings)))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Mutually exclusive gene pairs", strrep("#", 2)))
cat("\n\n")

# Retrieve heterotypic mutually exclusive pairs
MECR_pairs = retrieve_MECR_pairs(markers_as_strings, gene_panel_480)
colnames(MECR_pairs) = c("gene1", "gene2")

cat("\n\n")
print(htmltools::tagList(DT::datatable(MECR_pairs)))
cat("\n\n")

################################################################################
#'                         2. RETRIEVE PATHS TO PROSEG OUTPUTS
################################################################################

cat("\n\n")
cat(sprintf("%s **Proseg outputs paths**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Retrieve the paths to the segmentation outputs from each Proseg run 
generated during the hyperparameter grid search in the 02_Cell_Segmentation 
step. A table containing the run metadata is also initialized at this step.")
cat("\n\n")

# Retrieve all proseg hyperparameter combinations
hyperparameter_combinations = expand.grid(hyperparameters_values)

# Initialize statistics for each TMA
TMA_statistics = list()

# Retrieve paths for each TMA
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA))
  cat("\n\n")
  
  # Initialize data frame summarizing proseg runs
  TMA_statistics[[TMA]] = hyperparameter_combinations
  
  # Retrieve condition name and segmentation path of each hyperparameter combination 
  tmp = lapply(seq_len(nrow(TMA_statistics[[TMA]])), function(i) {
    
    # Retrieve condition name of current row
    row = TMA_statistics[[TMA]][i, ]
    current_params_combination = vapply(colnames(hyperparameter_combinations),
                                        function(hp) paste0(hp, row[[hp]]),
                                        character(1))
    
    # Return condition name and proseg output path for current params combination
    return(c(condition = paste(current_params_combination, collapse = "_"),
             data_dir = file.path(PATH_HYPERPARAMETER_SEARCH_XENIUM_OUTPUT[[TMA]],
                                  paste(current_params_combination, collapse = "/"),
                                  PROSEG_OUTPUTS_RELATIVE_PATHS)))
  })
  tmp = do.call(rbind, tmp)
  
  # Add condition name and proseg outputs paths to statistics
  TMA_statistics[[TMA]]$condition = tmp[, "condition"]
  TMA_statistics[[TMA]]$data_dir = tmp[, "data_dir"]
  
  # Add row corresponding to Xenium data before proseg segmentation
  xenium_condition_columns = c(colnames(hyperparameter_combinations), "condition", "data_dir")
  xenium_condition = as.data.frame(as.list(rep(NA, length(xenium_condition_columns))))
  colnames(xenium_condition) = xenium_condition_columns
  xenium_condition$data_dir = XENIUM_DATA_PATH[[TMA]]
  xenium_condition$condition = XENIUM_CONDITION_NAME
  TMA_statistics[[TMA]] = rbind(TMA_statistics[[TMA]], xenium_condition)
  
  # Display TMA_statistics[[TMA]] in report
  cat("\n\n")
  print(htmltools::tagList(DT::datatable(TMA_statistics[[TMA]])))
  cat("\n\n")
  
}

################################################################################
#'          3. RETRIEVE PROSEG AND INITIAL XENIUM SEGMENTATION STATISTICS
################################################################################

cat("\n\n")
cat(sprintf("%s **MECR by hyperparameter combination**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Retrieve statistics for both the Proseg-generated segmentations and the 
initial Xenium segmentation. These metrics will be used to evaluate the effect 
of hyperparameter choices on segmentation quality.

For each segmentation run, the following metrics are retrieved:

  - Number of cells in the data (**nb_cells**): Since Proseg works by refining 
existing cells rather than generating new ones, this value should remain close 
to the number of cells in the initial Xenium segmentation.
  - Number of conserved cells (**nb_conserved**): Conserved cells are those that 
pass quality control (QC) criteria based on feature counts (CELL_MINIMUM_FEATURES) 
and cell size (CELL_MAXIMUM_AREA). Because Proseg performs additional filtering 
on the transcripts, it typically reduces the number of conserved cells obtained 
during QC. Ideally, the selected hyperparameters should preserve a conserved 
cell count as close as possible to that of the initial Xenium segmentation. 
Since only conserved cells are retained for downstream analyses, excessive cell 
loss may create gaps in the tissue and negatively impact subsequent steps, such 
as neighborhood analysis.
  - Number of low feature cells (**nb_low_feature**): Low feature cells are cells 
that fail the feature count filter but pass the cell size filter. Proseg tends 
to increase the number of low feature cells during segmentation refinement, and 
the selected hyperparameters should minimize their occurrence.
  - Number of high area cells (**nb_high_area**): High area cells are cells that 
pass the feature count filter but fail the cell size filter.
  - Number of low feature and high area cells (**nb_low_features_high_area**): 
Cells that fail both quality thresholds.
  - Mean mutually exclusive co-expression rate (**mean_MECR**): The mutually 
exclusive co-expression rate is computed per cell, and the mean value is then 
calculated across the dataset. MECR is defined as the fraction of cells expressing 
both markers divided by the fraction of cells expressing at least one of the 
markers (Satija et al. Biorxiv 2024).
  - Median mutually exclusive co-expression rate (**median_MECR**): The median 
value of the mutually exclusive co-expression rate across all cells, computed 
in the same way as above.")
cat("\n\n")

# Store MECR in a list
MECR_scores = list()

# Initialize metrics columns in statistics
for (TMA in TMA_COVERED){
  TMA_statistics[[TMA]]$nb_cells = NA
  TMA_statistics[[TMA]]$nb_conserved = NA
  TMA_statistics[[TMA]]$nb_low_feature = NA
  TMA_statistics[[TMA]]$nb_high_area = NA
  TMA_statistics[[TMA]]$nb_low_features_high_area = NA
  TMA_statistics[[TMA]]$mean_MECR = NA
  TMA_statistics[[TMA]]$median_MECR = NA
}

cat("\n\n")
cat(sprintf("%s Gene pairs MECR by condition", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over all proseg hyperparameter combinations
  for (i in 1:nrow(TMA_statistics[[TMA]])){
    
    sprintf("Processing %s ...", TMA_statistics[[TMA]][i, ]$condition)
    
    # Make sure segmentation is available for current condition
    if (file.exists(TMA_statistics[[TMA]][i, ]$data_dir)){
      
      # Retrieve name for current condition
      condition_name = TMA_statistics[[TMA]][i, ]$condition
      
      ##################
      # 3.1. READ DATA #
      ##################
      
      # Read Xenium data for current condition
      temp_xenium = LoadXenium_new(TMA_statistics[[TMA]][i, ]$data_dir)
      
      ########################
      # 3.2. QUALITY CONTROL #
      ########################
      
      # Tag cells that should be removed based on filtering parameters for current condition
      temp_xenium = pastis_quality_control(temp_xenium,
                                           minimum_features = CELL_MINIMUM_FEATURES,
                                           maximum_area = CELL_MAXIMUM_AREA)
      
      #####################
      # 3.3. COMPUTE MECR #
      #####################
      
      # Compute MECR scores for each gene pair in current condition
      MECR_scores[[TMA]][[condition_name]] = calculate_MECR_score_parallelized(seurat_object = temp_xenium, 
                                                                               mutually_exclusive_gene_pairs = MECR_pairs, 
                                                                               sample_name = condition_name,
                                                                               batch_size = BATCH_SIZE,
                                                                               workers = CORES)
      
      cat("\n\n")
      cat(sprintf("%s %s ", strrep("#", 4), condition_name))
      cat("\n\n")
      
      # Display MECR scores for current condition
      cat("\n\n")
      print(htmltools::tagList(DT::datatable(MECR_scores[[TMA]][[condition_name]])))
      cat("\n\n")
      
      ##############################################
      # 3.4. SAVE STATISTICS FOR CURRENT CONDITION #
      ##############################################
      
      # Retrieve quality control related data
      quality_control_stats = as.data.frame(table(temp_xenium@meta.data[["QC_result"]]))
      
      # Fill data df with corresponding metrics
      TMA_statistics[[TMA]][i, ]$nb_cells = length(colnames(temp_xenium))
      TMA_statistics[[TMA]][i, ]$nb_conserved = quality_control_stats[quality_control_stats$Var1 == "conserved", ]$Freq
      TMA_statistics[[TMA]][i, ]$nb_low_feature = quality_control_stats[quality_control_stats$Var1 == "low_feature", ]$Freq
      TMA_statistics[[TMA]][i, ]$nb_high_area = quality_control_stats[quality_control_stats$Var1 == "high_area", ]$Freq
      TMA_statistics[[TMA]][i, ]$nb_low_features_high_area = quality_control_stats[quality_control_stats$Var1 == "low_features_high_area", ]$Freq
      TMA_statistics[[TMA]][i, ]$mean_MECR = mean(MECR_scores[[TMA]][[condition_name]]$MECR)
      TMA_statistics[[TMA]][i, ]$median_MECR = median(MECR_scores[[TMA]][[condition_name]]$MECR)
      
    }
  }
  
  cat("\n\n")
  cat(sprintf("%s %s global statistics", strrep("#", 3), TMA))
  cat("\n\n")
  
  cat("\n\n")
  print(htmltools::tagList(DT::datatable(TMA_statistics[[TMA]])))
  cat("\n\n")
  
}
