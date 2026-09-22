
## @knitr compute_mecr

cat("\n\n")
cat(sprintf("%s **MECR**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step consists of computing the MECR score for each gene pair across all 
TMAs and segmentation methods under study. This enables comparison of MECR values 
and highlights the improvement achieved by the Proseg-refined segmentation 
relative to the initial Xenium segmentation.")
cat("\n\n")

################################################################################
#'                              1. RETRIEVE MECR PAIRS 
#'                              
#' Retrieve mutually exclusive gene pairs computed in 01_Find_Mutually_Exclusive_Gene_Pairs.   
#'                        
################################################################################

# Retrieve MECR genes per cell type
markers_as_strings = as.data.frame(read_excel(MECR_GENES_PER_CELL_TYPE))
markers_as_strings = markers_as_strings[!(markers_as_strings$...1 %in% c("LAMP3cDC_9", "Granulo")), ] # Remove cell types that would lead to shared MECR genes between cell types
rownames(markers_as_strings) = markers_as_strings$...1
markers_as_strings$...1 = NULL

# Retrieve gene panel
gene_panel_480 = rownames(xenium_data[["TMA002"]][["xenium_segmentation"]])

# Retrieve heterotypic MECR pairs
MECR_pairs = retrieve_MECR_pairs(markers_as_strings, gene_panel_480)

################################################################################
#'                                 2. COMPUTE MECR      
#'                                 
#' MECR scores for each gene pair are computed across all TMAs and segmentation 
#' methods under study. The following visualizations are then generated to enable 
#' comparison between conditions:
#' 
#'     MECR score distributions: Distribution of MECR values for each mutually 
#'     exclusive gene pair, stratified by segmentation method.
#'     
#'     MECR scores by condition: A table reporting both the MECR score and the 
#'     expression rate (i.e., the proportion of cells expressing at least one 
#'     gene in the mutually exclusive pair) for each gene pair, across all TMAs 
#'     and segmentation methods under study.     
#'                        
################################################################################

###########################
## 2.1. COMPUTE ALL MECR ##
###########################

# Store MECR for each condition
MECR_scores = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Compute MECR scores for each gene pair in current condition
    MECR_scores[[TMA]][[segmentation_method]] = calculate_MECR_score_parallelized(seurat_object = xenium_data[[TMA]][[segmentation_method]], 
                                                                                  mutually_exclusive_gene_pairs = MECR_pairs, 
                                                                                  sample_name = paste0(TMA, "_", segmentation_method),
                                                                                  batch_size = BATCH_SIZE,
                                                                                  workers = CORES)
    
  }
}

####################################################
## 2.2. CONCATENATE ALL AVAILABLE MECR DATAFRAMES ##
####################################################

# Prepare data frame to merge all MECR scores
merged_MECR_scores = data.frame()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Add MECR scores from current condition
    merged_MECR_scores = rbind(merged_MECR_scores, MECR_scores[[TMA]][[segmentation_method]])
    
  }
}

############################
## 2.3. PLOT MECR RESULTS ##
############################

cat("\n\n")
cat("MECR scores for each gene pair are computed across all TMAs and segmentation 
methods under study. The following visualizations are then generated to enable 
comparison between conditions:

  - MECR score distributions: Distribution of MECR values for each mutually 
exclusive gene pair, stratified by segmentation method.
  - MECR scores by condition: A table reporting both the MECR score and the 
expression rate (i.e., the proportion of cells expressing at least one gene in 
the mutually exclusive pair) for each gene pair, across all TMAs and segmentation 
methods under study.")
cat("\n\n")

###############################
## 2.3.1. MERGED MECR SCORES ##
###############################

cat("\n\n")
cat(sprintf("%s MECR scores ", strrep("#", 2)))
cat("\n\n")

# Retrieve outliers gene pairs for each segmentation
MECR_outliers = merged_MECR_scores %>%
  group_by(sample_ID) %>%
  mutate(
    Q1 = quantile(MECR, 0.25),
    Q3 = quantile(MECR, 0.75),
    IQR = Q3 - Q1,
    is_outlier = MECR < Q1 - 1.5 * IQR | MECR > Q3 + 1.5 * IQR
  ) %>%
  filter(is_outlier)

# Plot MECR boxplot containing all studied conditions
cat("\n\n")
merged_MECR_plot = ggplot(merged_MECR_scores, aes(x = sample_ID, y = MECR, fill = sample_ID)) +
  geom_boxplot(outlier.shape = NA) +  # Hide default outliers
  geom_point(data = MECR_outliers, aes(x = sample_ID, y = MECR), color = "black", size = 2) +
  scale_fill_manual(values = COLOR_PALETTE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "MECR scores over segmentation methods", x = "Condition", y = "MECR score")
print(merged_MECR_plot)
cat("\n\n")

#####################################
## 2.3.2. MECR SCORES BY CONDITION ##
#####################################

cat("\n\n")
cat(sprintf("%s MECR scores by condition ", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s ", strrep("#", 4), segmentation_method))
    cat("\n\n")
    
    # Plot MECR scores by condition
    cat("\n\n")
    print(htmltools::tagList(DT::datatable(MECR_scores[[TMA]][[segmentation_method]])))
    cat("\n\n")
    
  }
}

################################################################################
#'                     3. COMPARE MECR IN SEGMENTATION METHODS
################################################################################

cat("\n\n")
cat(sprintf("%s **MECR Segmentation method comparison**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Display a scatter plot comparing the MECR scores obtained with each 
segmentation method for every mutually exclusive gene pair.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA))
  cat("\n\n")

  # Define methods to compare
  method1 = names(TMA_METADATA[[TMA]])[1]
  method2 = names(TMA_METADATA[[TMA]])[2]
  
  # Create data frame with MECR values to compare
  MECR_comparison = data.frame(setNames(list(MECR_scores[[TMA]][[method1]]$MECR,
                                             MECR_scores[[TMA]][[method2]]$MECR),
                                        c(method1,
                                          method2)))
  
  # Scatter plot of MECR vs MECR of current segmentation methods
  cat("\n\n")
  print(ggplot(MECR_comparison, aes(x=.data[[method1]], y=.data[[method2]])) + 
          geom_point() +
          ggtitle(sprintf("%s MECR vs %s MECR (TMA = %s)", method1, method2, TMA)) +
          geom_abline(color = "red") +
          xlab(sprintf("%s MECR", method1)) +
          ylab(sprintf("%s MECR", method2)))
  cat("\n\n")
  
}
