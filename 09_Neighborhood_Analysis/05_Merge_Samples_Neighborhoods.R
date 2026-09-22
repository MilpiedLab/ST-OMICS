
## @knitr merge_samples_neighborhoods_heatmaps

cat("\n\n")
cat(sprintf("%s **MERGE SAMPLES NEIGHBORHOODS - HEATMAPS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to harmonize BANKSY neighborhoods across samples by assigning 
a common neighborhood index to neighborhoods considered to be similar.

To achieve this, a simple hierarchical approach is used to identify similar BANKSY 
neighborhoods. Specifically, the dendrogram obtained from the heatmap based on 
Jensen–Shannon divergence is used to define groups of similar neighborhoods. 
A predefined number of hierarchical clusters is extracted from this dendrogram, 
and BANKSY neighborhoods belonging to the same hierarchical cluster are considered 
equivalent. Their identifiers are then harmonized by assigning a shared index 
across samples.")
cat("\n\n")

################################################################################
#'                    1. MERGE SAMPLES NEIGHBORHOODS - HEATMAPS
################################################################################

cat("\n\n")
cat("In this step, the following sections are provided to visualize the neighborhood 
merging process:

For each resolution and each predefined number of hierarchical clusters:

  - A heatmap similar to the one described in the previous section is displayed. 
In addition, an additional color bar is shown on the left side of the heatmap, 
indicating the hierarchical cluster assignment of each BANKSY neighborhood.")
cat("\n\n")

# Init list to store plots with tissue view
hierarchical_clusters_on_tissue_plots = list()

# Iterate over clustering resolutions covered
for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s res = %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over number of cluster to retrieve from dendrogram
  for (number_of_clusters in HIERARCHICAL_CLUSTERS_RANGE){
    
    cat("\n\n")
    cat(sprintf("%s hclust = %s", strrep("#", 3), number_of_clusters), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    #########################################
    ## 1.1. RETRIEVE HIERARCHICAL CLUSTERS ##
    #########################################
    
    # Retrieve dendrogram for current condition
    temp_dendrogram = row_dend(neighborhood_similarities_heatmaps[[as.character(res)]])
    temp_dendrogram = as.hclust(temp_dendrogram)
    
    # Retrieve clusters from dendrogram
    temp_clusters = cutree(temp_dendrogram,
                           k = number_of_clusters)
    temp_clusters = data.frame(complete_tag = names(temp_clusters),
                               hierarchical_cluster = as.numeric(temp_clusters))
    
    ###############################
    ## 1.2. UPDATE SEURAT OBJECT ##
    ###############################
    
    # Retrieve current neighborhood assignment for samples
    temp_neighborhood_layer_name = sprintf("%s_snn_res.%s", 
                                           BANKSY_ASSAY_NAME, 
                                           res)
    temp_neighborhood_assigment = data.frame(sample = banksy_xenium_data$SAMPLE_ID,
                                             neighborhood = banksy_xenium_data@meta.data[[temp_neighborhood_layer_name]])
    temp_neighborhood_assigment$complete_tag = paste0(temp_neighborhood_assigment$sample,
                                                      "_",
                                                      temp_neighborhood_assigment$neighborhood)
    
    # Assign to each cell its hierarchical cluster
    temp_neighborhood_assigment$hierarchical_cluster = NA
    temp_neighborhood_assigment$hierarchical_cluster = apply(X = temp_neighborhood_assigment,
                                                             MARGIN = 1,
                                                             FUN = function(x){
                                                               x[["hierarchical_cluster"]] = temp_clusters[temp_clusters[["complete_tag"]] == x[["complete_tag"]], ][["hierarchical_cluster"]]
                                                             })
    
    # Update Seurat object with new hierarchical cluster layer
    temp_hierarchical_layer_name = sprintf("%s_hierarchical_clusters_res_%s_%s", 
                                           BANKSY_ASSAY_NAME, 
                                           res,
                                           number_of_clusters)
    banksy_xenium_data@meta.data[[temp_hierarchical_layer_name]] = temp_neighborhood_assigment$hierarchical_cluster
    
    ########################################
    ## 1.3. DISPLAY HIERARCHICAL CLUSTERS ##
    ########################################
    
    ## 1.3.1. Retrieve plot data ##
    
    # Retrieve JS divergence matrix for current resolution
    current_JS_results = as.matrix(global_JS_results[[as.character(res)]])
    
    ## 1.3.2. JS color gradient ##
    
    # Define color function
    col_fun = colorRamp2(c(0, quantile(current_JS_results, 0.5), max(current_JS_results)),
                         c("blue", "white", "red"))
    
    ## 1.3.3. Columns labeling ##
    
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
    
    ## 1.3.4. Rows labeling ##
    
    # Assign a color to each cluster
    clusters_colors = COLOR_PALETTE[1:length(unique(temp_clusters$hierarchical_cluster))]
    names(clusters_colors) = unique(temp_clusters$hierarchical_cluster)
    
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
    left_annot = rowAnnotation(Cluster = temp_clusters$hierarchical_cluster,
                               Sample = row_sample_names,
                               Neighborhood = anno_text(row_neighborhood_index,
                                                        gp = gpar(fontsize = 8)),
                               col = list(Sample = SAMPLE_COLORS,
                                          Cluster = clusters_colors),
                               show_legend = FALSE,
                               show_annotation_name = FALSE)
    
    ## 1.3.5. Rows cell type bars ##
    
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
    
    ## 1.3.6. Plot heatmap ##
    
    cat("\n\n")
    print(Heatmap(current_JS_results,
                  name = "JS",
                  col = col_fun,
                  column_title = sprintf("Cell neighborhoods JS divergence (res = %s)", res),
                  top_annotation = top_annot,
                  left_annotation = left_annot,
                  right_annotation = right_annot,
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_row_names = FALSE,
                  show_column_names = FALSE))
    cat("\n\n")
    
    ## 1.3.7. Create plots with hierarchical clusters on tissue ##
    
    # Iterate over available samples in merged data
    for (sample in unique(banksy_xenium_data@meta.data[["SAMPLE_ID"]])){

      # Subset current sample from global merged data
      sample_subset = subset(banksy_xenium_data,
                             subset = seurat_clusters == sample)
      
      # Activate current hierarchical cluster layer in meta data
      sample_subset = activateSeuratAssay(sample_subset, 
                                          ident = temp_hierarchical_layer_name)
      
      # Create tissue view
      hierarchical_clusters_on_tissue_plots[[as.character(res)]][[as.character(number_of_clusters)]][[sample]] = ImageDimPlot(sample_subset, 
                                                                                                                              cols = clusters_colors, 
                                                                                                                              size = 0.75) + 
        ggtitle(sprintf("%s merged neighborhoods from JS divergence hierarchical clustering",
                        sample)) +
        labs(fill = sprintf("Hierarchical clusters"))
      
    }
  }
}

## @knitr merge_samples_neighborhoods_tissue_views

################################################################################
#'                 2. MERGE SAMPLES NEIGHBORHOODS - TISSUE VIEW
################################################################################

cat("\n\n")
cat(sprintf("%s **MERGE SAMPLES NEIGHBORHOODS - TISSUE VIEW**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("In this step, the following sections are provided to visualize the 
neighborhood merging process:

For each resolution and each predefined number of hierarchical clusters:

  - A visualization of the sample tissue, where cells are colored according to 
their unified BANKSY neighborhood index.")
cat("\n\n")

# Iterate over clustering resolutions covered
for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s res = %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over number of cluster to retrieve from dendrogram
  for (number_of_clusters in HIERARCHICAL_CLUSTERS_RANGE){
    
    cat("\n\n")
    cat(sprintf("%s hclust = %s", strrep("#", 3), number_of_clusters), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Iterate over available samples in merged data
    for (sample in unique(banksy_xenium_data@meta.data[["SAMPLE_ID"]])){
      
      cat("\n\n")
      cat(sprintf("%s %s", strrep("#", 4), sample))
      cat("\n\n")
      
      cat("\n\n")
      print(hierarchical_clusters_on_tissue_plots[[as.character(res)]][[as.character(number_of_clusters)]][[sample]])
      cat("\n\n")
      
    }
  }
}
