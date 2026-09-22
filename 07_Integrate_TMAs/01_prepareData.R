
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to integrate the data within the merged 
Seurat object in order to remove any batch effects that may be observed. Indeed, 
in Section 06_Merge_TMAs, a batch effect related to the TMA was identified when 
clustering cells from the merged dataset.

Therefore, Harmony will be used for the integration. The batch-corrected 
low-dimensional embedding obtained from Harmony will help reduce batch effects 
between the TMAs and enable better characterization of cell types during manual 
annotation in downstream analyses.")
cat("\n\n")

################################################################################
#'                              1. LOAD MERGED DATA
#' 
#' Load the data after merging and after filtering out low-quality cells and 
#' samples, as described in Section '06_Merge_TMAs'.
#' 
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Read merged TMA data
xenium_data = readRDS(PATH_MERGED_TMA_DATA)

################################################################################
#'                    2. PLOT MERGED DATA BEFORE INTEGRATION
################################################################################

cat("\n\n")
cat(sprintf("%s **MERGED TMA DATA BEFORE INTEGRATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Visualize the merged and filtered data using the following sections:

  - Raw clusters: UMAP visualization of clustered cells for each clustering 
resolution tested.
  - Elbow plot: Visualization used to determine the optimal number of dimensions.
  - Clusters by sample ID: UMAP colored by sample ID.
  - Clusters by TMA ID: UMAP colored by TMA ID, highlighting a batch effect 
associated with the TMA.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Raw Clusters", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

for (res in RES_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), res))
  cat("\n\n")
  
  # Activate current assay (required resolution)
  xenium_data = activateSeuratAssay(xenium_data, 
                                    ident = sprintf("%s_%s", 
                                                    CLUSTER_LAYER_PREFIX, 
                                                    gsub("\\.", "", as.character(res))), 
                                    active_assay = XENIUM_ASSAY_NAME, 
                                    umap_layer_name = UMAP_LAYER_NAME)
  
  cat("\n\n")
  print(DimPlot(xenium_data,
                label = TRUE) +
          ggtitle(sprintf("Merged TMAs clusters (res = %s)", res)))
  cat("\n\n")
}

cat("\n\n")
cat(sprintf("%s Elbow plot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(get_elbow_plot(xenium_data, 
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

cat("\n\n")
cat(sprintf("%s Clusters by TMA ID", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
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
