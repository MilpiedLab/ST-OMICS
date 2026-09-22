
## @knitr compare_inter_sample_neighborhoods

cat("\n\n")
cat(sprintf("%s **COMPARE INTER-SAMPLE NEIGHBORHOODS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Similarly to the intra-sample comparison of neighborhoods, inter-sample 
comparisons of average cell-neighbors cell-type composition in BANKSY neighborhoods 
can be performed using the Jensen–Shannon divergence.

This analysis enables the identification of similarities between neighborhoods 
across different samples and helps to highlight corresponding or shared tissue 
regions.

For each clustering resolution, a heatmap is displayed showing pairwise comparisons 
between neighborhoods across samples using the Jensen–Shannon divergence.

On the right side, the cell-type proportions corresponding to each neighborhood 
are shown. On the left and top margins, a color bar indicates the sample of origin 
for each neighborhood, along with the corresponding BANKSY neighborhood index 
within the sample.")
cat("\n\n")

#########################################
## 1. RETRIEVE NEIGHBORHOODS META DATA ##
#########################################

# Retrieve general cell meta data
global_neighborhoods_meta_data = data.frame(cell = Cells(banksy_xenium_data),
                                            cell_type = banksy_xenium_data[[NEIGHBORS_ANNOTATION]],
                                            sample = banksy_xenium_data$SAMPLE_ID,
                                            x = banksy_xenium_data$x_centroid,
                                            y = banksy_xenium_data$y_centroid)

# Retrieve neighborhoods assignment for each clustering resolution covered
for (res in RES_RANGE){
  temp_neighborhood_layer_name = sprintf("%s_snn_res.%s", 
                                         BANKSY_ASSAY_NAME, 
                                         res)
  global_neighborhoods_meta_data[[temp_neighborhood_layer_name]] = banksy_xenium_data@meta.data[[temp_neighborhood_layer_name]]
  
}

############################################
## 2. RETRIEVE NEIGHBORHOODS COMPOSITIONS ##
############################################

# Init list to store neighborhood compositions for each condition
global_neighborhoods_composition = list()

# Iterate over clustering resolutions covered
for (res in RES_RANGE){
  
  # Define the neighborhood column to consider
  temp_neighborhood_layer_name = sprintf("%s_snn_res.%s", 
                                         BANKSY_ASSAY_NAME, 
                                         res)
  
  # Iterate over samples covered
  for (sample in unique(global_neighborhoods_meta_data$sample)){
    
    ########################
    ## 2.1. SUBSET SAMPLE ##
    ########################
    
    # Subset current sample from global merged data
    sample_subset = subset(banksy_xenium_data,
                           subset = seurat_clusters == sample)
    
    # Activate current neighborhood layer in meta data
    sample_subset = activateSeuratAssay(sample_subset, 
                                        ident = temp_neighborhood_layer_name)
    
    #############################################
    ## 2.2. RETRIEVE NEIGHBORHOOD COMPOSITIONS ##
    #############################################
    
    # Iterate over neighborhoods found in this sample
    for (neighborhood in 0:(length(unique(sample_subset$seurat_clusters))-1)){
      
      # Retrieve cell type composition of cell neighborhood for current neighborhood assignment
      current_neighborhood_composition = get_cell_counts_in_cells_neighborhood(seurat_object = sample_subset,
                                                                               cluster_layer_name = temp_neighborhood_layer_name,
                                                                               active_assay = XENIUM_ASSAY_NAME,
                                                                               cluster_of_interest = as.character(neighborhood),
                                                                               k_nearest_neighbors = KNN_NEIGHBORHOOD_COMPOSITION,
                                                                               neighbors_annotation = NEIGHBORS_ANNOTATION)
      
      # Keep composition expressed as proportions
      cell_type_proportions = current_neighborhood_composition$proportions
      rownames(cell_type_proportions) = sprintf("%s_%s_%s", sample, res, neighborhood)
      
      # Store distribution of cell types for current neighborhood
      global_neighborhoods_composition[[as.character(res)]][[sprintf("%s_%s", sample, neighborhood)]] = cell_type_proportions
      
    }
  }
}

##########################################
## 3. ASSESS NEIGHBORHOODS SIMILARITIES ##
##########################################

# Init list to store JS results for each covered clustering resolution
global_JS_results = list()

# Iterate over clustering resolutions covered
for (res in RES_RANGE){
  
  # Retrieve neighborhood composition for current clustering resolution
  current_neighborhoods_compositions = global_neighborhoods_composition[[as.character(res)]]
  
  # Init JS result matrix
  current_JS_results = data.frame(matrix(nrow = length(names(current_neighborhoods_compositions)), 
                                         ncol = length(names(current_neighborhoods_compositions))))
  colnames(current_JS_results) = names(current_neighborhoods_compositions)
  rownames(current_JS_results) = names(current_neighborhoods_compositions)
  
  # Compute JS divergence for every combination
  for (i in seq_len(nrow(current_JS_results))) {
    for (j in seq_len(ncol(current_JS_results))) {
      
      # Prepare data to compute JS
      js_data = as.matrix(rbind(current_neighborhoods_compositions[[rownames(current_JS_results)[i]]],
                                current_neighborhoods_compositions[[colnames(current_JS_results)[j]]]))
      
      # Compute JS for current combination of sample and neighborhood
      current_JS_results[i, j] = JSD(js_data)
      
    }
  }
  
  # Store JS results for current resolution
  global_JS_results[[as.character(res)]] = current_JS_results
  
}

########################################
## 4. PLOT NEIGHBORHOODS SIMILARITIES ##
########################################

# Init list to store heatmap objects to retrieve clusters in later step
neighborhood_similarities_heatmaps = list()

# Iterate over clustering resolutions covered
for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  #############################
  ## 4.1. RETRIEVE PLOT DATA ##
  #############################
  
  # Retrieve JS divergence matrix for current resolution
  current_JS_results = as.matrix(global_JS_results[[as.character(res)]])
  
  ############################
  ## 4.2. JS COLOR GRADIENT ##
  ############################
  
  # Define color function
  col_fun = colorRamp2(c(0, quantile(current_JS_results, 0.5), max(current_JS_results)),
                       c("blue", "white", "red"))

  ###########################
  ## 4.3. COLUMNS LABELING ##
  ###########################
  
  # Define sample col names
  col_sample_names = strsplit(colnames(current_JS_results), split = "_")
  col_sample_names = unlist(lapply(col_sample_names, FUN = function(x){
    paste0(x[1:length(x)-1], collapse = "_")
  }))
  
  # Define neighborhood col index
  col_neighborhood_index = strsplit(colnames(current_JS_results), split = "_")
  col_neighborhood_index = unlist(lapply(col_neighborhood_index, FUN = function(x){
    x[length(x)]
  }))
  
  # Define sample annotation on the left of heatmap
  top_annot = HeatmapAnnotation(Sample = col_sample_names,
                                Neighborhood = anno_text(col_neighborhood_index,
                                                         rot = 90,
                                                         gp = gpar(fontsize = 8)),
                                col = list(Sample = SAMPLE_COLORS),
                                show_legend = TRUE,
                                show_annotation_name = FALSE)
  
  ########################
  ## 4.4. ROWS LABELING ##
  ########################
  
  # Define sample row names
  row_sample_names = strsplit(rownames(current_JS_results), split = "_")
  row_sample_names = unlist(lapply(row_sample_names, FUN = function(x){
    paste0(x[1:length(x)-1], collapse = "_")
  }))
  
  # Define neighborhood row index
  row_neighborhood_index = strsplit(rownames(current_JS_results), split = "_")
  row_neighborhood_index = unlist(lapply(row_neighborhood_index, FUN = function(x){
    x[length(x)]
  }))
  
  # Define sample annotation on the left of heatmap
  left_annot = rowAnnotation(Sample = row_sample_names,
                             Neighborhood = anno_text(row_neighborhood_index,
                                                      gp = gpar(fontsize = 8)),
                             col = list(Sample = SAMPLE_COLORS),
                             show_legend = FALSE,
                             show_annotation_name = FALSE)
  
  ##############################
  ## 4.5. ROWS CELL TYPE BARS ##
  ##############################
  
  # Define cell types proportions
  barplot_list = lapply(unique(global_neighborhoods_meta_data$cell_type), function(cell_type) {
    
    # Retrieve proportions of current cell types for each neighborhood
    proportions = unlist(lapply(rownames(current_JS_results), function(x) {
      global_neighborhoods_composition[[as.character(res)]][[x]][[cell_type]]
    }))
    
    # Create complexheatmap barplot for current cell type
    row_anno_barplot(proportions,
                     gp = gpar(fill = COMPLETE_CELL_TYPES_COLOR_PALETTE[[cell_type]]),
                     border = FALSE,
                     ylim = c(0, 1),
                     axis_param = list(gp = gpar(fontsize = 5)))
  })
  
  # Define cell type names to assign
  names(barplot_list) = unique(global_neighborhoods_meta_data$cell_type)
  
  # Create row annotation objects
  right_annot = do.call(rowAnnotation, c(barplot_list, 
                                         list(annotation_name_rot = 45,
                                              annotation_name_gp = gpar(fontsize = 10),
                                              gap = unit(2, "mm"))))
  
  #######################
  ## 4.6. PLOT HEATMAP ##
  #######################
  
  # Plot heatmap
  cat("\n\n")
  neighborhood_similarities_heatmaps[[as.character(res)]] = Heatmap(current_JS_results,
                                                                    name = "JS",
                                                                    col = col_fun,
                                                                    column_title = sprintf("Cell neighborhoods JS divergence (res = %s)", res),
                                                                    top_annotation = top_annot,
                                                                    left_annotation = left_annot,
                                                                    right_annotation = right_annot,
                                                                    cluster_rows = TRUE,
                                                                    cluster_columns = TRUE,
                                                                    show_row_names = FALSE,
                                                                    show_column_names = FALSE)
  print(neighborhood_similarities_heatmaps[[as.character(res)]])
  cat("\n\n")
  
}
