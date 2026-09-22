
## @knitr harmony_integration

cat("\n\n")
cat(sprintf("%s **HARMONY INTEGRATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step performs data integration using Harmony. The integration is carried 
out using the TMA ID as a covariate representing the source of batch effects.")
cat("\n\n")

################################################################################
#'                            1. INTEGRATION BY TMA ID
################################################################################

cat("\n\n")
cat("The following sections are provided to visualize and assess the effect of 
the integration:

  - Harmony convergence plot: Figure showing the evolution of the objective 
function across Harmony iterations, used to assess convergence of the algorithm.
  - Elbow plot: Visualization used to determine the number of dimensions to 
retain from the Harmony-corrected embedding for downstream analyses such as clustering.
  - UMAP by TMA ID: UMAP colored by TMA ID. Visually, the batch effect appears 
to be substantially reduced; however, quantitative metrics should be included 
for a more rigorous assessment.
  - UMAP by sample ID: UMAP colored by sample ID.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Harmony plot convergence", strrep("#", 2)))
cat("\n\n")

# Run Harmony integration
xenium_data = RunHarmony(xenium_data, 
                         HARMONY_INTEGRATION_COVARIABLE, 
                         early_stop = FALSE, 
                         plot_convergence = TRUE)

# Run UMAP from Harmony embedding
xenium_data = RunUMAP(xenium_data,
                      reduction = HARMONY_ASSAY_NAME,
                      dims = 1:HARMONY_DIM_USE)

cat("\n\n")
cat(sprintf("%s Elbow plot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(get_elbow_plot(xenium_data, 
                     dimension_principal_component_analysis = DIM_PCA,
                     dimension_clustering = HARMONY_DIM_USE,
                     reduction = HARMONY_ASSAY_NAME))
cat("\n\n")

# Save UMAP embeddings in meta.data
xenium_data = storeUMAPembeddings(seurat_object = xenium_data, 
                                  umap_layer_name = HARMONY_UMAP_LAYER_NAME)

cat("\n\n")
cat(sprintf("%s UMAP by TMA ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data, 
              group.by = "TMA_id") +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              split.by = "TMA_id",
              group.by = "TMA_id") +
        theme(legend.position = "none") +
        facet_wrap(~TMA_id, nrow = round(length(unique(xenium_data$TMA_id))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s UMAP by sample ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data, 
              group.by = "SAMPLE_ID"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              split.by = "SAMPLE_ID",
              group.by = "SAMPLE_ID") +
        theme(legend.position = "none") +
        facet_wrap(~SAMPLE_ID, nrow = round(length(unique(xenium_data$SAMPLE_ID))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

################################################################################
#'                                2. CLUSTERING
#' 
#' This step aims to cluster cells based on the Harmony-corrected embedding 
#' (i.e., the batch-corrected low-dimensional representation).
#' 
################################################################################

## Clustering ##
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
#'                       3. PLOT CLUSTERS AND MARKER GENES
################################################################################

cat("\n\n")
cat(sprintf("%s **HARMONY CLUSTERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to visualize the clusters obtained after Harmony integration. 
For each clustering resolution tested, the following sections are provided:

  - DimPlot: UMAP visualization of clustered cells.
  - Bar plot: Cell frequencies per cluster.
  - Individual clusters: UMAP visualization highlighting each cluster at the current resolution.

For each cluster shown above, marker genes are identified using the FindAllMarkers 
function from Seurat. The results of the marker analysis will be saved in the next 
step as a .xlsx file.")
cat("\n\n")

# Init storage of marker genes
markers_across_resolutions = list()

# Iterate over clustering resolutions covered
for (res in HARMONY_RESOLUTION_RANGE){
  
  # Define current cluster layer name
  temp_cluster_layer_name = sprintf("%s_%s", 
                                    HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                    gsub("\\.", "", as.character(res)))
  
  ########################
  ## 3.1. PLOT CLUSTERS ##
  ########################
  
  # Activate current assay (required resolution)
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = temp_cluster_layer_name, 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = HARMONY_UMAP_LAYER_NAME)
  
  cat("\n\n")
  cat(sprintf("%s Res %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(xenium_data) +
          ggtitle(sprintf("Clusters at res = %s", res)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Barplot", strrep("#", 3)))
  cat("\n\n")
  
  # Retrieve number of cells by clusters
  cluster_freq = as.data.frame(table(xenium_data$seurat_clusters))
  
  cat("\n\n")
  print(ggplot(cluster_freq, aes(x=Var1, 
                                 y=Freq)) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq, vjust = -0.3)) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                legend.position = "none") +
          ylab("Cell number") +
          ggtitle(sprintf("Cluster cell frequencies at res = %s", res)))
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
  
  ##################
  ## 3.2. MARKERS ##
  ##################
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = xenium_data, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = temp_cluster_layer_name,
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  # Store current markers
  markers_across_resolutions[[as.character(res)]] = markers
  
}

################################################################################
#'                              4. SAVE MARKER GENES
#'                              
#' Save the marker genes identified for the clusters above as a .xlsx file
#'                               
################################################################################

# Save markers in xlsx file
save_marker_genes_xlsx(marker_genes = markers_across_resolutions, 
                       output_directory = PATH_ANALYSIS_OUTPUT, 
                       output_file_name = "integrated_merged_TMA_markers.xlsx")
