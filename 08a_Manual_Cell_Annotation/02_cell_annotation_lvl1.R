
## @knitr cell_annotation_lvl1

cat("\n\n")
cat(sprintf("%s **INTEGRATED MERGED TMAs DATA CLUSTERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to annotate clusters obtained from the full Harmony-integrated 
dataset. The annotated clusters will correspond to the Level 1 cell types 
described previously.")
cat("\n\n")

################################################################################
#'                 1. CLUSTERING OF HARMONY-INTEGRATED DATASET
#' 
#' Cluster cells using the Harmony-integrated assay at each tested resolution.
#' 
################################################################################

# Activate harmony UMAP layer
xenium_data = activateSeuratAssay(xenium_data, 
                                  umap_layer_name = HARMONY_UMAP_LAYER_NAME)

## Clustering ##
# Constructs a K-nearest neighbors (KNN) graph based on the selected PCs
xenium_data = FindNeighbors(xenium_data, 
                            reduction = HARMONY_ASSAY_NAME, 
                            dims = 1:HARMONY_DIM_USE, 
                            verbose = FALSE)
# Create output cluster layers names
cluster_layers = c()
for (res in LEVEL_1_RES_RANGE){
  cluster_layers = c(cluster_layers, sprintf("%s_%s", 
                                             LEVEL_1_LAYER_PREFIX, 
                                             gsub("\\.", "", as.character(res))))
}
# Applies graph-based clustering (Louvain algorithm) on the KNN graph
xenium_data = FindClusters(xenium_data,
                           resolution = LEVEL_1_RES_RANGE,
                           verbose = TRUE,
                           cluster.name = cluster_layers)

################################################################################
#'                      2. PLOT INTEGRATED MERGED TMA DATA
################################################################################

cat("\n\n")
cat("For each clustering resolution tested, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to retain 
from the Harmony-corrected embedding for downstream analyses such as clustering.
  - Clusters frequencies: Cell frequencies per cluster.
  - UMAP by sample ID: UMAP colored by sample ID.
  - Individual clusters: UMAP visualization highlighting each cluster at the 
current resolution.")
cat("\n\n")

# Display the different clustering available
for (res in LEVEL_1_RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Activate current condition
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = sprintf("%s_%s", 
                                                    LEVEL_1_LAYER_PREFIX, 
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
                       dimension_principal_component_analysis = LEVEL_1_DIM_PCA,
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
#'                          3. CHARACTERISE LVL1 CLUSTERS
################################################################################

cat("\n\n")
cat(sprintf("%s **CHARACTERISE LVL1 CLUSTERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to identify marker genes and corresponding clusters in the 
reference single-cell dataset that can help annotate clusters in the merged 
Xenium data. To do so, marker genes are computed for each cluster. The top marker 
genes are then used as a gene module to assess similarity with cells in the 
reference single-cell dataset using the AddModuleScore function from Seurat.

For each clustering resolution tested on the Harmony-integrated assay, and for 
each cluster within these resolutions, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Feature plots: UMAP-based visualization of the top 9 marker genes per cluster, 
selected based on log fold-change expression.
  - Violin plots: Expression distribution of the top 9 marker genes across clusters.
  - Reference UMAP: Reference single-cell UMAP visualization, with cells colored 
according to the module score derived from the marker genes of the Harmony cluster 
under study.
  - Reference AddModuleScore table: Table reporting the mean module score per 
cell type in the reference single-cell dataset, computed using the marker-derived 
gene module of the Harmony cluster under study.")
cat("\n\n")

# Initialize list to store markers at the different resolutions
markers_across_resolutions = list()

# Iterate over resolution covered for lvl1 cell annotation
for (res in LEVEL_1_RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  ##########################
  ## 3.1. COMPUTE MARKERS ##
  ##########################
  
  # Activate current assay (required resolution)
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = sprintf("%s_%s", 
                                                    LEVEL_1_LAYER_PREFIX, 
                                                    gsub("\\.", "", as.character(res))), 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = HARMONY_UMAP_LAYER_NAME)
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = xenium_data, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = sprintf("%s_%s", 
                                                     LEVEL_1_LAYER_PREFIX,
                                                     gsub("\\.", "", as.character(res))),
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  # Store current resolution clusters markers
  markers_across_resolutions[[as.character(res)]] = markers
  
  ##########################
  ## 3.2. DISPLAY MARKERS ##
  ##########################
  
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
                                                   cluster_layer_name = paste0(LEVEL_1_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
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
                                            cluster_layer_name = paste0(LEVEL_1_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                            cluster_number = (as.numeric(cluster) - 1), 
                                            top_n_markers = TOP_N_MARKERS)
    
    cat("\n\n")
    grid.arrange(grobs=temp_violins, top = sprintf("Cluster %s marker genes expression", as.numeric(cluster) - 1))
    cat("\n\n")
    
    ##################################################
    ## 3.3. FIND CORRESPONDING REFERENCE SC CLUSTER ##
    ##################################################
    
    # Retrieve markers from current cluster
    current_cluster_markers_genes = markers[markers$cluster == (cluster-1), ]
    current_cluster_markers_genes = top_n(current_cluster_markers_genes, 
                                          TOP_N_ANNOTATION_MARKERS_LVL1, 
                                          avg_log2FC)
    current_cluster_markers_genes = current_cluster_markers_genes$gene
    module = list(current_cluster_markers_genes)
    
    # Compute module score for current cluster
    similarity_results = compute_reference_clusters_modulescores(query_seurat = xenium_data,
                                                                 reference_seurat = sc_reference_data,
                                                                 query_markers_list = module)
    
    cat("\n\n")
    cat(sprintf("%s Reference UMAP", strrep("#", 4)))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("Module: %s", module))
    print(similarity_results$plot)
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Reference AddModuleScore", strrep("#", 4)))
    cat("\n\n")
    
    cat("\n\n")
    print(htmltools::tagList(DT::datatable(similarity_results$similarity_scores)))
    cat("\n\n")
    
  }
}

# Save markers from all resolutions in xlsx file
save_marker_genes_xlsx(marker_genes = markers_across_resolutions, 
                       output_directory = PATH_ANALYSIS_OUTPUT, 
                       output_file_name = "level_1_clusters_markers.xlsx")

################################################################################
#'                   4. DISPLAY AVAILABLE MANUAL ANNOTATION
################################################################################

cat("\n\n")
cat(sprintf("%s **AVAILABLE LVL1 MANUAL ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("In this case, two of the clustering resolutions tested earlier were used for 
manual annotation. This is done to illustrate how the choice of resolution influences 
the resulting annotation and the biological interpretation of the data.

For each resolution for which manual annotation is provided, the following sections 
are displayed:

  - DimPlot: UMAP visualization of cells colored by the corresponding manual annotation.
  - Cell type distribution: Cell frequencies for each manually identified cell type.")
cat("\n\n")

# Iterate over resolutions with annotation available
for (res in names(LEVEL_1_TYPING_METADATA)){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Activate current assay (required resolution)
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = sprintf("%s_%s", 
                                                    LEVEL_1_LAYER_PREFIX, 
                                                    gsub("\\.", "", as.character(res))), 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = HARMONY_UMAP_LAYER_NAME)
  
  # Initialize cell annotation meta data
  cell_type_metadata = data.frame(cluster = xenium_data$seurat_clusters,
                                  cell_type = "Unknown")
  
  # Iterate over clusters to generate corresponding plots
  for (cluster in sort(as.numeric(unique(xenium_data$seurat_clusters))-1)){
    
    # Check if current cluster corresponds to a cell type
    if (cluster %in% names(LEVEL_1_TYPING_METADATA[[as.character(res)]])){
      
      cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = LEVEL_1_TYPING_METADATA[[as.character(res)]][[as.character(cluster)]]
      
    } else {
      
      cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = "Unknown"
      
    }
  }
  
  # Add cell typing to meta data
  xenium_data@meta.data[[LEVEL1_ANNOTATION_LAYER_NAME]] = cell_type_metadata$cell_type
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(xenium_data,
                group.by = LEVEL1_ANNOTATION_LAYER_NAME,
                cols = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
          ggtitle(sprintf("Level 1 cell annotation at res = %s", res)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Cell types distribution", strrep("#", 3)))
  cat("\n\n")
  
  # Prepare data
  cell_type_distribution_data = as.data.frame(table(xenium_data[[LEVEL1_ANNOTATION_LAYER_NAME]]))
  
  cat("\n\n")
  print(ggplot(cell_type_distribution_data, aes(x=.data[[LEVEL1_ANNOTATION_LAYER_NAME]], 
                                                y=Freq,
                                                fill = .data[[LEVEL1_ANNOTATION_LAYER_NAME]])) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq, vjust = -0.3)) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                legend.position = "none") +
          scale_fill_manual(values = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
          ylab("Cell number") +
          ggtitle(sprintf("Level 1 cell types counts (res = %s)", res)))
  cat("\n\n")
  
}

################################################################################
#'                            5. FINAL LVL1 ANNOTATION
#'                            
#' Add a metadata column to the Seurat object corresponding to the Level 1 manual 
#' annotation selected as the most appropriate based on the chosen clustering 
#' resolution.                           
#' 
################################################################################

# Activate current assay (required resolution)
xenium_data = activateSeuratAssay(xenium_data, 
                                  ident = sprintf("%s_%s", 
                                                  LEVEL_1_LAYER_PREFIX, 
                                                  gsub("\\.", "", as.character(LEVEL_1_ANNOTATION_RES))), 
                                  active_assay = XENIUM_ASSAY_NAME, 
                                  umap_layer_name = HARMONY_UMAP_LAYER_NAME)

# Initialize cell annotation meta data
cell_type_metadata = data.frame(cluster = xenium_data$seurat_clusters,
                                cell_type = "Unknown")

# Iterate over clusters to generate corresponding plots
for (cluster in sort(as.numeric(unique(xenium_data$seurat_clusters))-1)){
  
  # Check if current cluster corresponds to a cell type
  if (cluster %in% names(LEVEL_1_TYPING_METADATA[[as.character(LEVEL_1_ANNOTATION_RES)]])){
    
    cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = LEVEL_1_TYPING_METADATA[[as.character(LEVEL_1_ANNOTATION_RES)]][[as.character(cluster)]]
    
  } else {
    
    cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = "Unknown"
    
  }
}

# Final annotation step for level 1 cell types
xenium_data@meta.data[[LEVEL1_ANNOTATION_LAYER_NAME]] = cell_type_metadata$cell_type

################################################################################
#'                        6. DISPLAY FINAL LVL1 ANNOTATION
################################################################################

cat("\n\n")
cat(sprintf("%s **FINAL LVL1 ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Display the Level 1 manual annotation selected using the clustering resolution 
deemed most appropriate, which will serve as the basis for proceeding to Level 2 
annotation. The following sections are displayed:

  - DimPlot: UMAP visualization of cells colored by the corresponding manual annotation.
  - Cell type distribution: Cell frequencies for each manually identified cell type.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s DimPlot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = LEVEL1_ANNOTATION_LAYER_NAME,
              cols = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
        ggtitle(sprintf("Level 1 cell annotation at res = %s", LEVEL_1_ANNOTATION_RES)))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Cell types distribution", strrep("#", 2)))
cat("\n\n")

# Prepare data
cell_type_distribution_data = as.data.frame(table(xenium_data[[LEVEL1_ANNOTATION_LAYER_NAME]]))

cat("\n\n")
print(ggplot(cell_type_distribution_data, aes(x=.data[[LEVEL1_ANNOTATION_LAYER_NAME]], 
                                              y=Freq,
                                              fill = .data[[LEVEL1_ANNOTATION_LAYER_NAME]])) + 
        geom_bar(stat = "identity") +
        geom_text(aes(label = Freq, vjust = -0.3)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none") +
        scale_fill_manual(values = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
        ylab("Cell number") +
        ggtitle(sprintf("Level 1 cell types counts (res = %s)", LEVEL_1_ANNOTATION_RES)))
cat("\n\n")
