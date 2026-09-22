
## @knitr cluster_filtered_merged_data

cat("\n\n")
cat(sprintf("%s **FILTERED MERGED TMA DATA**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step performs cell clustering on the fully merged dataset after quality 
filtering. Both low-quality cells and low-quality samples, as previously identified, 
have been removed from the data.")
cat("\n\n")

################################################################################
#'                                 1. CLUSTERING
#'                                 
#' Cluster cells using predefined dimensionality and clustering resolution 
#' parameters.   
#'                              
################################################################################

# Log normalize data
filtered_complete_xenium_data = log_normalisation(filtered_complete_xenium_data,
                                                  active_assay = XENIUM_ASSAY_NAME)

# Cluster cells from current condition (TMA)
filtered_complete_xenium_data = wrapped_clustering_stomics_imaging(seurat_object = filtered_complete_xenium_data, 
                                                                   active_assay = XENIUM_ASSAY_NAME, 
                                                                   cluster_layer_prefix = CLUSTER_LAYER_PREFIX, 
                                                                   umap_layer_name = UMAP_LAYER_NAME,
                                                                   dimension_principal_component_analysis = DIM_PCA,
                                                                   dimension_clustering = DIM_USE,
                                                                   resolution_range = RES_RANGE)

################################################################################
#'                         2. VISUALIZE FILTERED MERGED DATA
################################################################################

cat("\n\n")
cat("Visualize the data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells for every clustering resolution covered.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by sample ID: UMAP colored by sample ID
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect associated with the TMA.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Raw Clusters", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), res))
  cat("\n\n")
  
  # Activate current assay (required resolution)
  filtered_complete_xenium_data = activateSeuratAssay(filtered_complete_xenium_data, 
                                                      ident = sprintf("%s_%s", 
                                                                      CLUSTER_LAYER_PREFIX, 
                                                                      gsub("\\.", "", as.character(res))), 
                                                      active_assay = XENIUM_ASSAY_NAME, 
                                                      umap_layer_name = UMAP_LAYER_NAME)
  
  cat("\n\n")
  print(DimPlot(filtered_complete_xenium_data,
          label = TRUE) +
    ggtitle(sprintf("Merged TMAs clusters (res = %s)", res)))
  cat("\n\n")
}

cat("\n\n")
cat(sprintf("%s Elbow plot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(get_elbow_plot(filtered_complete_xenium_data, 
                     dimension_principal_component_analysis = DIM_PCA,
                     dimension_clustering = DIM_USE))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by sample ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(filtered_complete_xenium_data, 
              group.by = "SAMPLE_ID"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(filtered_complete_xenium_data,
              split.by = "SAMPLE_ID",
              group.by = "SAMPLE_ID") +
        theme(legend.position = "none") +
        facet_wrap(~SAMPLE_ID, nrow = round(length(unique(filtered_complete_xenium_data$SAMPLE_ID))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by TMA ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(filtered_complete_xenium_data, 
              group.by = "TMA_id") +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(filtered_complete_xenium_data,
              split.by = "TMA_id",
              group.by = "TMA_id") +
        theme(legend.position = "none") +
        facet_wrap(~TMA_id, nrow = round(length(unique(filtered_complete_xenium_data$TMA_id))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

## @knitr save_filtered_clustered_data

################################################################################
#'                      3. SAVE CLUSTERED FILTERED MERGED DATA
#'                     
#' Save data processed in current section.
#'                      
################################################################################

# Save merged TMA data
saveRDS(filtered_complete_xenium_data,
        file = file.path(PATH_ANALYSIS_OUTPUT,
                         "filtered_clustered_TMAs.Rds"))
