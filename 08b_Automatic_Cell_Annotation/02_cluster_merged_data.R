
## @knitr cluster_merged_data

cat("\n\n")
cat(sprintf("%s **INTEGRATED TMAs CLUSTERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to cluster cells from the full Harmony-integrated dataset. 
The resulting clusters will subsequently be used for annotation with scType.")
cat("\n\n")

################################################################################
#'                                 1. CLUSTERING
#' 
#' Cluster cells using the Harmony-integrated assay across the range of tested 
#' resolutions.
#' 
################################################################################

# Activate Harmony integrated PCs
xenium_data = activateSeuratAssay(xenium_data,
                                  umap_layer_name = HARMONY_UMAP_LAYER_NAME)

## Find clusters at required resolutions from harmony PCs ##
# Constructs a K-nearest neighbors (KNN) graph based on the selected PCs
xenium_data = FindNeighbors(xenium_data, 
                            reduction = HARMONY_ASSAY_NAME, 
                            dims = 1:HARMONY_DIM_USE, 
                            verbose = FALSE)
# Create output cluster layers names
cluster_layers = c()
for (res in HARMONY_RESOLUTION_RANGE){
  cluster_layers = c(cluster_layers, sprintf("%s_%s", 
                                             HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                             gsub("\\.", "", as.character(res))))
}
# Applies graph-based clustering (Louvain algorithm) on the KNN graph
xenium_data = FindClusters(xenium_data,
                           resolution = HARMONY_RESOLUTION_RANGE,
                           verbose = TRUE,
                           cluster.name = cluster_layers)

################################################################################
#'                       2. PLOT INTEGRATED MERGED TMA DATA
################################################################################

cat("For each clustering resolution tested, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to 
retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - Cluster frequencies: Cell proportions per cluster.
  - UMAP by sample ID: UMAP colored by sample ID.
  - Individual clusters: UMAP visualization highlighting each cluster at the 
current resolution.")

# Display the different clustering available
for (res in HARMONY_RESOLUTION_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Activate current condition
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = sprintf("%s_%s", 
                                                    HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                                    gsub("\\.", "", as.character(res))), 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = HARMONY_UMAP_LAYER_NAME)
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(xenium_data))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Elbow plot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(get_elbow_plot(seurat_object = xenium_data, 
                       dimension_principal_component_analysis = DIM_PCA,
                       dimension_clustering = HARMONY_DIM_USE))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Clusters frequencies", strrep("#", 3)))
  cat("\n\n")
  
  # Retrieve number of cells by clusters
  cluster_freq = as.data.frame(table(xenium_data$seurat_clusters))
  
  cat("\n\n")
  print(ggplot(cluster_freq, aes(x=Var1, 
                                 y=Freq)) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq, vjust = -0.3)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                legend.position = "none") +
          ylab("Cell count") +
          xlab("Clusters"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s UMAP by sample ID", strrep("#", 3)), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Whole", strrep("#", 4)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(xenium_data,
                group.by = "SAMPLE_ID") +
          theme(legend.position = "none"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Splitted", strrep("#", 4)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(xenium_data,
                split.by = "SAMPLE_ID",
                group.by = "SAMPLE_ID") +
          theme(legend.position = "none") +
          facet_wrap(~SAMPLE_ID, nrow = round(length(unique(xenium_data$SAMPLE_ID))/2)) +
          theme(legend.position = "none"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Individual clusters", strrep("#", 3)), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Retrieve data to plot
  complete_plot_data = data.frame(cell = Cells(xenium_data),
                                  umap_1 = xenium_data@reductions[["umap"]]@cell.embeddings[, 1],
                                  umap_2 = xenium_data@reductions[["umap"]]@cell.embeddings[, 2],
                                  cluster = xenium_data$seurat_clusters,
                                  color = "grey")
  
  # Iterate over clusters in this resolution
  for (cluster in sort(unique(as.numeric(complete_plot_data$cluster))) - 1){
    
    # Update color for plot
    complete_plot_data$color = ifelse(complete_plot_data$cluster == cluster, "red", "grey")
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 4), cluster))
    cat("\n\n")
    
    cat("\n\n")
    print(ggplot(complete_plot_data, aes(x = umap_1, y = umap_2)) +
            geom_scattermore(
              data = subset(complete_plot_data, color == "grey"),
              aes(color = color),
              pointsize = 1
            ) +
            geom_scattermore(
              data = subset(complete_plot_data, color == "red"),
              aes(color = color),
              pointsize = 1
            ) +
            scale_color_identity() +
            ggtitle(sprintf("Cluster %s (res = %s)", cluster, res)))
    cat("\n\n")
    
  }
}

################################################################################
#'                          3. HARMONY CLUSTERS MARKERS
################################################################################

cat("\n\n")
cat(sprintf("%s **HARMONY CLUSTERS MARKERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to identify marker genes associated with the Harmony-derived 
clusters obtained in the previous step.

For each clustering resolution tested on the Harmony-integrated assay, and for 
each cluster within these resolutions, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Feature plots: UMAP-based visualization of the top 9 marker genes per cluster, 
selected based on log fold-change expression.
  - Violin plots: Expression distribution of the top 9 marker genes across clusters.")
cat("\n\n")

# Store markers in list to save them later
markers_across_resolutions = list()

# Iterate over each computed clustering resolution
for (res in HARMONY_RESOLUTION_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Define layers for current condition
  current_cluster_layer = sprintf("%s_%s", 
                                  HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                  gsub("\\.", "", as.character(res)))
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = xenium_data, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = current_cluster_layer,
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  # Store current markers
  markers_across_resolutions[[as.character(res)]] = markers
  
  # Activate current assay (required resolution)
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = current_cluster_layer, 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = HARMONY_UMAP_LAYER_NAME)
  
  # Get current cluster identities
  clusters = Idents(xenium_data)
  
  # Convert identities to a factor ordered numerically
  Idents(xenium_data) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
  
  # Iterate over clusters to generate corresponding plots
  for (cluster in sort(as.numeric(unique(markers$cluster)))){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), as.numeric(cluster) - 1), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s DimPlot ", strrep("#", 4)))
    cat("\n\n")
    
    # Set up cluster colors
    temp_clusters = levels(Idents(xenium_data))
    temp_clusters_colors = rep("grey80", length(temp_clusters))
    names(temp_clusters_colors) = temp_clusters
    temp_clusters_colors[cluster] = "red"
    
    cat("\n\n")
    print(DimPlot(xenium_data,
                  label = TRUE,
                  cols = temp_clusters_colors) +
            ggtitle(sprintf("Cluster %s (res = %s)", as.numeric(cluster) - 1, res)))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Features ", strrep("#", 4)))
    cat("\n\n")
    
    # Compute features for feature plot
    temp_features = get_marker_genes_feature_plots(seurat_object = xenium_data, 
                                                   all_marker_genes = markers,
                                                   cluster_layer_name = current_cluster_layer,
                                                   cluster_number = (as.numeric(cluster) - 1), 
                                                   top_n_markers = TOP_N_MARKERS)
    
    cat("\n\n")
    grid.arrange(grobs=temp_features, top = sprintf("%s marker genes expression", as.numeric(cluster) - 1))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Violins ", strrep("#", 4)))
    cat("\n\n")
    
    temp_violins = get_marker_genes_violins(seurat_object = xenium_data, 
                                            all_marker_genes = markers,
                                            cluster_layer_name = current_cluster_layer,
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
                       output_file_name = "merged_data_markers.xlsx")
