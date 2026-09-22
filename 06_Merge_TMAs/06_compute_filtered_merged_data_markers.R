
## @knitr markers_filtered_clustered_data

################################################################################
#'                COMPUTE MARKERS FOR FILTERED MERGED DATA CLUSTERS
################################################################################

cat("\n\n")
cat(sprintf("%s **FILTERED MERGED TMA MARKER GENES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step computes gene markers for the clusters obtained earlier from the 
filtered, merged dataset. This allows identification of the broad cell types 
present in the data and evaluation of how well they can be distinguished from one 
another, providing insight into the quality of the segmentation.
    
For each clustering resolution and corresponding clusters, the following figures 
are available:
    
  - DimPlot: UMAP visualization of current cluster of interest.
  - Features: UMAP-based visualization of the top 9 marker genes per cluster, 
selected based on log fold-change expression.
  - Violins: Display the expression distribution of the top 9 marker genes across 
clusters.")
cat("\n\n")

# Store markers in list to save them later
markers_across_resolutions = list()

# Iterate over each computed clustering resolution
for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = filtered_complete_xenium_data, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  # Store current markers
  markers_across_resolutions[[as.character(res)]] = markers
  
  # Activate current assay (required resolution)
  filtered_complete_xenium_data = activateSeuratAssay(filtered_complete_xenium_data, 
                                                      ident = sprintf("%s_%s", 
                                                                      CLUSTER_LAYER_PREFIX, 
                                                                      gsub("\\.", "", as.character(res))), 
                                                      active_assay = XENIUM_ASSAY_NAME, 
                                                      umap_layer_name = UMAP_LAYER_NAME)
  
  # Get current cluster identities
  clusters = Idents(filtered_complete_xenium_data)
  
  # Convert identities to a factor ordered numerically
  Idents(filtered_complete_xenium_data) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
  
  # Iterate over clusters to generate corresponding plots
  for (cluster in sort(as.numeric(unique(markers$cluster)))){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), as.numeric(cluster) - 1), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s DimPlot ", strrep("#", 4)))
    cat("\n\n")
    
    # Set up cluster colors
    temp_clusters = levels(Idents(filtered_complete_xenium_data))
    temp_clusters_colors = rep("grey80", length(temp_clusters))
    names(temp_clusters_colors) = temp_clusters
    temp_clusters_colors[cluster] = "red"
    
    cat("\n\n")
    print(DimPlot(filtered_complete_xenium_data,
                  label = TRUE,
                  cols = temp_clusters_colors) +
            ggtitle(sprintf("Cluster %s (res = %s)", as.numeric(cluster) - 1, res)))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Features ", strrep("#", 4)))
    cat("\n\n")
    
    # Compute features for feature plot
    temp_features = get_marker_genes_feature_plots(seurat_object = filtered_complete_xenium_data, 
                                                   all_marker_genes = markers,
                                                   cluster_layer_name = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                                   cluster_number = (as.numeric(cluster) - 1), 
                                                   top_n_markers = TOP_N_MARKERS)
    
    cat("\n\n")
    grid.arrange(grobs=temp_features, top = sprintf("%s marker genes expression", as.numeric(cluster) - 1))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Violins ", strrep("#", 4)))
    cat("\n\n")
    
    temp_violins = get_marker_genes_violins(seurat_object = filtered_complete_xenium_data, 
                                            all_marker_genes = markers,
                                            cluster_layer_name = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                            cluster_number = (as.numeric(cluster) - 1), 
                                            top_n_markers = TOP_N_MARKERS)
    
    cat("\n\n")
    grid.arrange(grobs=temp_violins, top = sprintf("Cluster %s marker genes expression", as.numeric(cluster) - 1))
    cat("\n\n")
    
  }
}

# Save markers in xlsx file
save_marker_genes_xlsx(marker_genes = markers_across_resolutions, 
                       output_directory = PATH_ANALYSIS_OUTPUT, 
                       output_file_name = "filtered_merged_data_markers.xlsx")
