
## @knitr compare_intra_sample_neighborhoods

cat("\n\n")
cat(sprintf("%s **COMPARE INTRA-SAMPLE NEIGHBORHOODS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to visualize the intra-sample comparison of neighborhoods using 
the sliced Wasserstein distance and the Jensen–Shannon divergence. The following 
outputs are provided to illustrate this comparison:

For each sample and clustering resolution combination, two heatmaps are displayed:

  - A heatmap showing pairwise comparisons between neighborhoods within the sample 
using the Jensen–Shannon divergence.
  - A heatmap showing pairwise comparisons between neighborhoods within the sample 
using the sliced Wasserstein distance.

For both heatmaps, the relative cell-type proportions of each neighborhood are 
displayed alongside the matrix. The neighborhood indices are indicated on both the 
rows and columns of the heatmaps.
    
The Jensen–Shannon divergence and sliced Wasserstein distance values are highlighted 
in blue in the heatmap when both metrics meet the thresholds defined in analysisParams.R. 
When both criteria are satisfied, this indicates BANKSY neighborhoods that may be 
candidates for merging, as they likely represent the same underlying tissue region.")
cat("\n\n")

################################################################################
#'                         PLOT NEIGHBORHOODS DIVERGENCE
################################################################################

# Store plot to display later
wasserstein_vs_js_plots = list()

# Activate current assay (sample id layer)
banksy_xenium_data = activateSeuratAssay(banksy_xenium_data, 
                                         ident = "SAMPLE_ID", 
                                         active_assay = XENIUM_ASSAY_NAME)

# Iterate over available samples in merged data
for (sample in unique(banksy_xenium_data@meta.data[["SAMPLE_ID"]])){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), sample), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Subset current sample from global merged data
  sample_subset = subset(banksy_xenium_data,
                         subset = seurat_clusters == sample)
  
  # Iterate over resolutions
  for (res in RES_RANGE){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), res), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    ########################################
    ## 1. FIND SIMILAR NEIGHBORHOOD PAIRS ##
    ########################################
    
    # Retrieve distance matrices for current condition
    temp_js_divergences = js_divergence_matrices[[sample]][[as.character(res)]]
    temp_wasserstein_distances = wasserstein_matrices[[sample]][[as.character(res)]]
    
    # Tag neighborhood pairs that are below the similarity threshold
    wasserstein_similar_domains = as.data.frame(which(temp_wasserstein_distances < WASSERSTEIN_THRESHOLD, arr.ind = TRUE))
    js_divergences_similar_domains = as.data.frame(which(temp_js_divergences < JS_THRESHOLD, arr.ind = TRUE))
    
    ## JS being asymmetric, let's consider the criteria met if one of the two comparison is below the threshold ##
    # Add reversed pairs
    inverse_rows = js_divergences_similar_domains[, c("col", "row")]
    # Combine and remove duplicates
    js_divergences_similar_domains = unique(rbind(js_divergences_similar_domains, inverse_rows))
    
    # Retrieve similar pairs of neighborhoods between the two distance metrics
    similar_domains = intersect(wasserstein_similar_domains, js_divergences_similar_domains)
    
    # Remove pairs comparing a neighborhood to itself
    # similar_domains = similar_domains[which(similar_domains$row != similar_domains$col),]
    
    ###################
    ## 2. JS HEATMAP ##
    ###################
    
    ##############################
    ## 2.1. rows cell type bars ##
    ##############################
    
    # Retrieve general cell meta data
    global_neighborhoods_meta_data = data.frame(cell = Cells(sample_subset),
                                                cell_type = sample_subset[[NEIGHBORS_ANNOTATION]],
                                                sample = sample_subset$SAMPLE_ID,
                                                x = sample_subset$x_centroid,
                                                y = sample_subset$y_centroid)
    
    # Define cell types proportions
    barplot_list = lapply(unique(global_neighborhoods_meta_data$cell_type), function(cell_type) {
      
      # Retrieve proportions of current cell types for each neighborhood
      proportions = unlist(lapply(rownames(temp_js_divergences), function(x) {
        sample_neighborhood_compositions[[as.character(res)]][[sample]][[x]][[cell_type]]
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
    
    ########################
    ## 2.2. rows labeling ##
    ########################
    
    # Define neighborhood annotation on the left of heatmap
    left_annot = rowAnnotation(Neighborhood = anno_text(rownames(temp_js_divergences),
                                                        gp = gpar(fontsize = 8)))
    
    #########################
    ## 2.3. create heatmap ##
    #########################
    
    # Create js divergence heatmap
    js_heatmap = Heatmap(temp_js_divergences,
                         name = "JS divergence",
                         col = colorRamp2(c(0, quantile(as.matrix(temp_js_divergences), 0.5), max(temp_js_divergences)),
                                          c("antiquewhite2", "brown1", "black")),
                         column_title = sprintf("JS divergences between neighborhoods (res = %s)", res),
                         cell_fun = function(j, i, x, y, width, height, fill){
                           if (any(similar_domains$row == i & similar_domains$col == j)){
                             grid.text(sprintf("%.3f", temp_js_divergences[i, j]), x, y, gp = gpar(fontsize = 10, col = "blue4"))
                           } else {
                             grid.text(sprintf("%.3f", temp_js_divergences[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
                           }},
                         left_annotation = left_annot,
                         right_annotation = right_annot,
                         cluster_rows = TRUE,
                         cluster_columns = TRUE,
                         show_row_names = FALSE,
                         show_column_names = TRUE)
    
    # Retrieve row and col order to apply on wasserstein heatmap
    #row_ord = row_order(js_heatmap)
    #col_ord = column_order(js_heatmap)
    
    ############################
    ## 3. WASSERSTEIN HEATMAP ##
    ############################
    
    # Create sliced wasserstein distance heatmap
    wasserstein_heatmap = Heatmap(temp_wasserstein_distances,
                                  name = "Sliced Wasserstein",
                                  col = c("antiquewhite2", "brown1", "black"),
                                  column_title = sprintf("Sliced wasserstein distances between neighborhoods (res = %s)", res),
                                  cell_fun = function(j, i, x, y, width, height, fill){
                                    if (any(similar_domains$row == i & similar_domains$col == j)){
                                      grid.text(sprintf("%.3f", temp_wasserstein_distances[i, j]), x, y, gp = gpar(fontsize = 10, col = "blue4"))
                                    } else {
                                      grid.text(sprintf("%.3f", temp_wasserstein_distances[i, j]), x, y, gp = gpar(fontsize = 10, col = "white"))
                                    }},
                                  cluster_rows = TRUE,
                                  cluster_columns = TRUE,
                                  right_annotation = right_annot,
                                  left_annotation = left_annot,
                                  #row_order = row_ord,
                                  #column_order = col_ord,
                                  show_row_names = FALSE,
                                  show_column_names = TRUE)
    
    ##########################
    ## 4. DISPLAY IN REPORT ##
    ##########################
    
    cat("\n\n")
    print(js_heatmap)
    cat("\n\n")
    
    cat("\n\n")
    print(wasserstein_heatmap)
    cat("\n\n")
    
  }
}
