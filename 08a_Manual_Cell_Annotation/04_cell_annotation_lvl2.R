
## @knitr cell_annotation_lvl2

cat("\n\n")
cat(sprintf("%s **CHARACTERISE LVL2 CLUSTERS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Now that the Level 1 manual annotation has been completed, this step focuses 
on subsetting each Level 1 cell type and re-clustering them individually, as 
described previously. This refined clustering approach enables the identification 
of finer heterogeneity within closely related cell populations.

The objective of this step is therefore to define and manually annotate Level 2 
cell types within the dataset.
    
Note that some closely related Level 1 cell types were pooled prior to re-clustering.")
cat("\n\n")

################################################################################
#'                          1. CHARACTERISE LVL2 CLUSTERS
################################################################################

cat("\n\n")
cat("This step aims to visualize the re-clustering performed for each Level 1 cell 
type or group of related Level 1 cell types.

For each subset, the following sections are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the number of dimensions to retain 
from the Harmony-corrected embedding for downstream analyses such as clustering.
  - For each cluster within the subset:
    * DimPlot: UMAP visualization highlighting the cells of the current cluster.
    * Feature plots: UMAP-based visualization of the top 9 marker genes for the 
current cluster, selected based on log fold-change expression.
    * Violin plots: Expression distribution of the top 9 marker genes across clusters.
    * Reference UMAP: UMAP visualization of the reference single-cell dataset, 
with cells colored according to the module score derived from the marker genes of 
the cluster under study. Only reference cell types relevant to the current Level 1 
subset are included in the module score computation.
    * Reference AddModuleScore table: Table reporting the mean module score per 
cell type in the reference single-cell dataset, computed using the marker-derived 
gene module of the cluster under study (via the AddModuleScore function in Seurat).
  - Subset Level 2 final annotation: Final Level 2 manual annotation assigned to 
the current Level 1 subset.")
cat("\n\n")

# Create output directory
dir.create(file.path(PATH_ANALYSIS_OUTPUT,
                     CELL_TYPE_LVL2_SUBCLUSTER_OUTPUT_DIR), 
           recursive = TRUE, 
           showWarnings = FALSE)

# Initialize annotation results
global_lvl2_annotation = list()

# Iterate over lvl1 cell types for further annotation
for (cell_type in names(LVL_1_CELL_TYPES_SUBSETS)){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), cell_type), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("**Current lvl1 cell types concerned: %s**", paste(LVL_1_REFERENCE_CELL_TYPES_SUBSETS[[cell_type]], collapse = ", ")))
  cat("\n\n")
  
  ################################
  ## 1.1. SUBSET LVL1 CELL TYPE ##
  ################################
  
  # Subset current cell type
  cell_type_subset = subset(xenium_data, 
                            subset = cell_type_lvl1 %in% LVL_1_CELL_TYPES_SUBSETS[[cell_type]])
 
  # Subset reference single cell data
  reference_subset = subset(sc_reference_data, 
                            subset = Celltype_Level_1 %in% LVL_1_REFERENCE_CELL_TYPES_SUBSETS[[cell_type]])
  
  #########################
  ## 1.2. CLUSTER SUBSET ##
  #########################
  
  # Log normalize data
  cell_type_subset = log_normalisation(cell_type_subset,
                                       active_assay = XENIUM_ASSAY_NAME)
  
  # Integrate and cluster cells from current condition (TMA)
  cell_type_subset = wrapped_harmony_integration(seurat_object = cell_type_subset, 
                                                 active_assay = XENIUM_ASSAY_NAME, 
                                                 cluster_layer_prefix = LEVEL_2_LAYER_PREFIX, 
                                                 umap_layer_name = LEVEL_2_UMAP_LAYER_NAME,
                                                 dimension_principal_component_analysis = LEVEL_2_DIM_PCA[[cell_type]],
                                                 dimension_clustering = LEVEL_2_DIM_USE[[cell_type]],
                                                 resolution_range = LEVEL_2_ANNOTATION_RES_RANGE[[cell_type]],
                                                 harmony_integration_covariable = HARMONY_INTEGRATION_COVARIABLE,
                                                 reduction = HARMONY_ASSAY_NAME)
  
  # Initialize list to store markers at the different resolutions
  markers_across_resolutions = list()
  
  # Iterate over resolution covered for lvl1 cell annotation
  for (res in LEVEL_2_ANNOTATION_RES_RANGE[[cell_type]]){
    
    cat("\n\n")
    cat(sprintf("%s Res = %s", strrep("#", 3), res), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    ###########################
    ## 1.3. DISPLAY CLUSTERS ##
    ###########################
    
    # Activate current assay (required resolution)
    cell_type_subset = activateSeuratAssay(cell_type_subset, 
                                           ident = sprintf("%s_%s", 
                                                           LEVEL_2_LAYER_PREFIX, 
                                                           gsub("\\.", "", as.character(res))), 
                                           active_assay = XENIUM_ASSAY_NAME, 
                                           umap_layer_name = LEVEL_2_UMAP_LAYER_NAME)
    
    cat("\n\n")
    cat(sprintf("%s DimPlot", strrep("#", 4)))
    cat("\n\n")
    
    cat("\n\n")
    print(DimPlot(cell_type_subset))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Elbow plot", strrep("#", 4)))
    cat("\n\n")
    
    cat("\n\n")
    print(get_elbow_plot(seurat_object = cell_type_subset, 
                         dimension_principal_component_analysis = LEVEL_2_DIM_PCA[[cell_type]],
                         dimension_clustering = LEVEL_2_DIM_USE[[cell_type]],
                         reduction = HARMONY_ASSAY_NAME))
    cat("\n\n")
    
    ##########################
    ## 1.4. COMPUTE MARKERS ##
    ##########################
    
    # Compute markers for current condition
    markers = find_all_markers(seurat_object = cell_type_subset, 
                               active_assay = XENIUM_ASSAY_NAME, 
                               cluster_layer = sprintf("%s_%s", 
                                                       LEVEL_2_LAYER_PREFIX,
                                                       gsub("\\.", "", as.character(res))),
                               logfc_threshold = LOGFC_THRESHOLD, 
                               min_pct = MIN_PCT, 
                               only_pos = TRUE,
                               n_top_markers = NULL)
    
    # Store current resolution clusters markers
    markers_across_resolutions[[as.character(res)]] = markers
    
    ##########################
    ## 1.5. DISPLAY MARKERS ##
    ##########################
    
    # Get current cluster identities
    clusters = Idents(cell_type_subset)
    
    # Convert identities to a factor ordered numerically
    Idents(cell_type_subset) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
    
    # Iterate over clusters to generate corresponding plots
    for (cluster in sort(as.numeric(unique(markers$cluster)))){
      
      cat("\n\n")
      cat(sprintf("%s Cluster %s", strrep("#", 4), as.numeric(cluster) - 1), " {.tabset .tabset-fade}\n\n")
      cat("\n\n")
      
      cat("\n\n")
      cat(sprintf("%s DimPlot ", strrep("#", 5)))
      cat("\n\n")
      
      # Set up cluster colors
      temp_clusters = levels(Idents(cell_type_subset))
      temp_clusters_colors = rep("grey80", length(temp_clusters))
      names(temp_clusters_colors) = temp_clusters
      temp_clusters_colors[cluster] = "red"
      
      cat("\n\n")
      print(DimPlot(cell_type_subset,
                    label = TRUE,
                    cols = temp_clusters_colors) +
              ggtitle(sprintf("Cluster %s (res = %s)", as.numeric(cluster) - 1, res)))
      cat("\n\n")
      
      cat("\n\n")
      cat(sprintf("%s Features ", strrep("#", 5)))
      cat("\n\n")
      
      # Compute features for feature plot
      temp_features = get_marker_genes_feature_plots(seurat_object = cell_type_subset, 
                                                     all_marker_genes = markers,
                                                     cluster_layer_name = paste0(LEVEL_2_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                                     cluster_number = (as.numeric(cluster) - 1), 
                                                     top_n_markers = TOP_N_MARKERS)
      
      cat("\n\n")
      grid.arrange(grobs=temp_features, top = sprintf("%s marker genes expression", as.numeric(cluster) - 1))
      cat("\n\n")
      
      cat("\n\n")
      cat(sprintf("%s Violins ", strrep("#", 5)))
      cat("\n\n")
      
      temp_violins = get_marker_genes_violins(seurat_object = cell_type_subset, 
                                              all_marker_genes = markers,
                                              cluster_layer_name = paste0(LEVEL_2_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                              cluster_number = (as.numeric(cluster) - 1), 
                                              top_n_markers = TOP_N_MARKERS)
      
      cat("\n\n")
      grid.arrange(grobs=temp_violins, top = sprintf("Cluster %s marker genes expression", as.numeric(cluster) - 1))
      cat("\n\n")
      
      ##################################################
      ## 1.6. FIND CORRESPONDING REFERENCE SC CLUSTER ##
      ##################################################
      
      # Retrieve markers from current cluster
      current_cluster_markers_genes = markers[markers$cluster == (cluster-1), ]
      current_cluster_markers_genes = top_n(current_cluster_markers_genes, 
                                            TOP_N_ANNOTATION_MARKERS_LVL2, 
                                            avg_log2FC)
      current_cluster_markers_genes = current_cluster_markers_genes$gene
      module = list(current_cluster_markers_genes)
      
      # Compute module score for current cluster
      similarity_results = compute_reference_clusters_modulescores(query_seurat = cell_type_subset,
                                                                   reference_seurat = reference_subset,
                                                                   query_markers_list = module)
      
      cat("\n\n")
      cat(sprintf("%s Reference UMAP", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      cat(sprintf("Module: %s", module))
      print(similarity_results$plot)
      cat("\n\n")
      
      cat("\n\n")
      cat(sprintf("%s Reference AddModuleScore table", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      print(htmltools::tagList(DT::datatable(similarity_results$similarity_scores)))
      cat("\n\n")
      
    }
  }
  
  # Save markers from all resolutions for current cell type in xlsx file
  save_marker_genes_xlsx(marker_genes = markers_across_resolutions, 
                         output_directory = file.path(PATH_ANALYSIS_OUTPUT,
                                                      CELL_TYPE_LVL2_SUBCLUSTER_OUTPUT_DIR), 
                         output_file_name = sprintf("%s_markers.xlsx", 
                                                    paste0(LVL_1_CELL_TYPES_SUBSETS[[cell_type]], collapse = "+")))
  
  ##################################
  ## 1.7. ANNOTATE CURRENT SUBSET ##
  ##################################
  
  # Activate current assay (required resolution)
  cell_type_subset = activateSeuratAssay(cell_type_subset, 
                                    ident = sprintf("%s_%s", 
                                                    LEVEL_2_LAYER_PREFIX, 
                                                    gsub("\\.", "", as.character(LEVEL_2_ANNOTATION_RES[[cell_type]]))), 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = LEVEL_2_UMAP_LAYER_NAME)
  
  # Initialize cell annotation meta data
  cell_type_metadata = data.frame(cell_id = Cells(cell_type_subset),
                                  cluster = cell_type_subset$seurat_clusters,
                                  cell_type = "Unknown")
  
  # Iterate over clusters to generate corresponding plots
  for (cluster in sort(as.numeric(unique(cell_type_subset$seurat_clusters))-1)){
    
    # Check if current cluster corresponds to a cell type
    if (cluster %in% names(LEVEL_2_TYPING_METADATA[[cell_type]][[as.character(res)]])){
      
      cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = LEVEL_2_TYPING_METADATA[[cell_type]][[as.character(res)]][[as.character(cluster)]]
      
    } else {
      
      cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = "Unknown"
      
    }
  }
  
  # Save lvl2 cell types annotation for later
  global_lvl2_annotation[[cell_type]] = cell_type_metadata
  
  ####################################
  ## 1.8. DISPLAY SUBSET ANNOTATION ##
  ####################################
  
  # Add lvl2 annotation to current subset
  cell_type_subset@meta.data[[LEVEL2_ANNOTATION_LAYER_NAME]] = cell_type_metadata$cell_type
  
  cat("\n\n")
  cat(sprintf("%s Subset lvl2 final annotation", strrep("#", 3)), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(cell_type_subset,
                group.by = LEVEL2_ANNOTATION_LAYER_NAME,
                cols = COMPLETE_CELL_TYPES_COLOR_PALETTE))
  cat("\n\n")
  
}

################################################################################
#'                            2. FINAL LVL2 ANNOTATION
#'                            
#' Add a metadata column to the Seurat object corresponding to the Level 2 
#' manual annotation.             
#'               
################################################################################

# Merge all available lvl2 annotations
final_lvl2_annotations = do.call(rbind, global_lvl2_annotation)

# Retrieve cells from global data
global_xenium_metadata = data.frame(cell_id = Cells(xenium_data),
                                    cell_type_lvl1 = xenium_data$cell_type_lvl1,
                                    cell_type_lvl2 = "Unknown")

# Retrieve lvl2 annotation if available, else retrieve lvl1 annotation
global_xenium_metadata$cell_type_lvl2 = sapply(global_xenium_metadata$cell_id, FUN = function(x){
  if (x %in% final_lvl2_annotations$cell_id){
    return(final_lvl2_annotations[final_lvl2_annotations$cell_id == x, ]$cell_type)
  } else {
    return(global_xenium_metadata[global_xenium_metadata$cell_id == x, ]$cell_type_lvl1)
  }
})

# Add lvl2 annotation to merged xenium TMAs data
xenium_data@meta.data[[LEVEL2_ANNOTATION_LAYER_NAME]] = global_xenium_metadata$cell_type_lvl2

################################################################################
#'                        3. DISPLAY FINAL LVL2 ANNOTATION
################################################################################

cat("\n\n")
cat(sprintf("%s **FINAL LVL2 ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Display the Level 2 manual annotation. The following sections are provided:

  - DimPlot: UMAP visualization of cells colored according to the corresponding 
Level 2 manual annotation.
  - Cell type distribution: Cell proportions for each manually identified Level 2 
cell type.
  - Sample IDs: Visualization of each level 2 cell type on the sample tissues.")
cat("\n\n")

######################
## 3.1. GLOBAL VIEW ##
######################

cat("\n\n")
cat(sprintf("%s DimPlot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = LEVEL2_ANNOTATION_LAYER_NAME,
              cols = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
        guides(color = guide_legend(override.aes = list(size=2),
                                    ncol=1)))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Cell types distribution", strrep("#", 2)))
cat("\n\n")

# Prepare data
cell_type_distribution_data = as.data.frame(table(xenium_data[[LEVEL2_ANNOTATION_LAYER_NAME]]))

cat("\n\n")
print(ggplot(cell_type_distribution_data, aes(x=.data[[LEVEL2_ANNOTATION_LAYER_NAME]], 
                                              y=Freq,
                                              fill = .data[[LEVEL2_ANNOTATION_LAYER_NAME]])) + 
        geom_bar(stat = "identity") +
        geom_text(aes(label = Freq, vjust = -0.3)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none") +
        scale_fill_manual(values = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
        ylab("Cell number") +
        ggtitle("Level2 cell types counts"))
cat("\n\n")

#######################################
## 3.2. DISPLAY CELL TYPES ON TISSUE ##
#######################################

# Activate current assay (sample id layer)
xenium_data = activateSeuratAssay(xenium_data, 
                                  ident = "SAMPLE_ID", 
                                  active_assay = XENIUM_ASSAY_NAME)

# Iterate over samples in the data
for (sample in unique(xenium_data@meta.data[["SAMPLE_ID"]])){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), sample), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Subset current sample from global merged data
  sample_subset = subset(xenium_data,
                         subset = seurat_clusters == sample)
  
  # Activate lvl1 cell type annotation layer
  sample_subset = activateSeuratAssay(sample_subset, 
                                      ident = LEVEL2_ANNOTATION_LAYER_NAME, 
                                      active_assay = XENIUM_ASSAY_NAME)
  
  # Iterate over lvl1 cell types in the data
  for (cell_type in unique(xenium_data@meta.data[[LEVEL2_ANNOTATION_LAYER_NAME]])){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), cell_type))
    cat("\n\n")
    
    # Subset current cell type from sample subset
    cell_type_subset = subset(sample_subset,
                              subset = seurat_clusters == cell_type)
    
    # Display cell type on tissue
    cat("\n\n")
    print(ImageDimPlot(cell_type_subset, 
                       size = 0.75) + 
            ggtitle(sprintf("%s on tissue", cell_type)) +
            labs(fill = sprintf("cell type = %s", cell_type)))
    cat("\n\n")
    
  }
}
