
## @knitr cluster_unfiltered_merged_data

cat("\n\n")
cat(sprintf("%s **UNFILTERED MERGED TMA DATA**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step performs cell clustering on the fully merged dataset prior to any 
quality filtering of samples or cells. It is intended to illustrate the impact of 
low-quality cells and samples, emphasizing the importance of applying quality 
control early in the analysis.")
cat("\n\n")

################################################################################
#'                                 1. CLUSTERING
#'                                 
#' Cluster cells using predefined dimensionality and clustering resolution 
#' parameters.   
#'                              
################################################################################

# Log normalize data
complete_xenium_data = log_normalisation(complete_xenium_data,
                                         active_assay = XENIUM_ASSAY_NAME)

# Cluster cells from current condition (TMA)
complete_xenium_data = wrapped_clustering_stomics_imaging(seurat_object = complete_xenium_data, 
                                                          active_assay = XENIUM_ASSAY_NAME, 
                                                          cluster_layer_prefix = MERGED_UNFILTERED_CLUSTER_LAYER_PREFIX,
                                                          umap_layer_name = MERGED_UNFILTERED_UMAP_LAYER_PREFIX,
                                                          dimension_principal_component_analysis = DIM_PCA,
                                                          dimension_clustering = DIM_USE,
                                                          resolution_range = MERGED_UNFILTERED_RES_RANGE)

################################################################################
#'                        2. VISUALIZE UNFILTERED MERGED DATA
################################################################################

cat("\n\n")
cat("Visualize the data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by QC: UMAP colored by cell quality, highlighting that low-feature 
cells tend to cluster together.
  - Clusters by sample quality: UMAP colored by sample quality, showing that 
low-quality samples cluster together. These samples did not pass the quality 
filters based on the transcript density analysis from the '05_TMA_Processing' 
step. They actually correspond to a different tissue type than the other samples, 
which explains the poor RNA capture by the gene panel designed for this Xenium 
experiment.
  - Clusters by sample ID: UMAP colored by sample ID, consistent with the 
clustering observed when colored by sample quality.
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect 
associated with the TMA.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Raw Clusters", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

for (res in MERGED_UNFILTERED_RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), res))
  cat("\n\n")
  
  # Activate current assay (required resolution)
  complete_xenium_data = activateSeuratAssay(complete_xenium_data, 
                                             ident = sprintf("%s_%s", 
                                                             MERGED_UNFILTERED_CLUSTER_LAYER_PREFIX, 
                                                             gsub("\\.", "", as.character(res))), 
                                             active_assay = XENIUM_ASSAY_NAME, 
                                             umap_layer_name = MERGED_UNFILTERED_UMAP_LAYER_PREFIX)
  
  cat("\n\n")
  print(DimPlot(complete_xenium_data,
                label = TRUE) +
          ggtitle(sprintf("Merged TMAs clusters (res = %s)", res)))
  cat("\n\n")
}

cat("\n\n")
cat(sprintf("%s Elbow plot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(get_elbow_plot(complete_xenium_data, 
                     dimension_principal_component_analysis = DIM_PCA,
                     dimension_clustering = DIM_USE))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by QC", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data, 
              group.by = "QC_result"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by sample quality", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data, 
              group.by = "sample_tissue_quality"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by sample ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data, 
              group.by = "SAMPLE_ID"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data,
              split.by = "SAMPLE_ID",
              group.by = "SAMPLE_ID") +
        theme(legend.position = "none") +
        facet_wrap(~SAMPLE_ID, nrow = round(length(unique(complete_xenium_data$SAMPLE_ID))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Clusters by TMA ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Whole", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data, 
              group.by = "TMA_id") +
        theme(legend.position = "none"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Splitted", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(complete_xenium_data,
              split.by = "TMA_id",
              group.by = "TMA_id") +
        theme(legend.position = "none") +
        facet_wrap(~TMA_id, nrow = round(length(unique(complete_xenium_data$TMA_id))/2)) +
        theme(legend.position = "none"))
cat("\n\n")

## @knitr save_unfiltered_clustered_data

################################################################################
#'                     3. SAVE CLUSTERED UNFILTERED MERGED DATA
#'                     
#' Save data processed in current section.
#'                      
################################################################################

# Save merged TMA data
saveRDS(complete_xenium_data,
        file = file.path(PATH_ANALYSIS_OUTPUT,
                         "unfiltered_clustered_TMAs.Rds"))
