
## @knitr process_segmentation_errors

cat("\n\n")
cat(sprintf("%s **CHARACTERIZE SEGMENTATION ERRORS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("In the previous step, some cells were annotated as segmentation errors because 
they exhibited marker genes from multiple distinct cell types, suggesting potential 
issues such as poor segmentation.

The aim of this analysis is to re-examine this subset of cells in order to recover 
as many valid cells as possible and reduce the number of cells classified as 
segmentation errors. To achieve this, cells annotated as segmentation errors at 
Level 1 are subsetted and re-clustered.

The goal of this finer clustering is to identify and recover genuine biological 
cell populations that were previously misclassified due to higher levels of 
segmentation uncertainty. Recovered cells can then be reassigned to the appropriate 
Level 1 manual annotation clusters whenever supported by their marker gene expression 
profiles.")
cat("\n\n")

################################################################################
##                      1. CHARACTERIZE SEGMENTATION ERRORS                   ##
################################################################################

cat("\n\n")
cat("This step aims to characterize the segmentation error cluster identified 
during Level 1 manual annotation. The following sections are displayed:

  - Markers: Marker genes identified for the segmentation error cluster using the 
FindAllMarkers function from Seurat.
  - Feature plots: UMAP-based visualization of marker gene expression within the 
segmentation error cluster, highlighting the non-specific and mixed expression 
patterns characteristic of this population.")
cat("\n\n")

##########################
## 1.1. COMPUTE MARKERS ##
##########################

# Retrieve cell meta data
cell_metadata = data.frame(cell = Cells(xenium_data),
                           umap_1 = xenium_data@reductions[["umap"]]@cell.embeddings[, 1],
                           umap_2 = xenium_data@reductions[["umap"]]@cell.embeddings[, 2],
                           cell_type = xenium_data$cell_type_lvl1,
                           color = "grey")
cell_metadata[cell_metadata$cell_type == "segmentation.error", ]$color = "red"

# Retrieve markers from cells tagged as segmentation errors
segmentation_error_markers = FindMarkers(xenium_data,
                                         ident.1 = cell_metadata[cell_metadata$cell_type == "segmentation.error", ]$cell,
                                         ident.2 = cell_metadata[cell_metadata$cell_type != "segmentation.error", ]$cell,
                                         logfc.threshold = LOGFC_THRESHOLD,
                                         min.pct = MIN_PCT, 
                                         only.pos = TRUE)

cat("\n\n")
cat(sprintf("%s Markers", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(htmltools::tagList(DT::datatable(segmentation_error_markers)))
cat("\n\n")

#######################
## 1.2. PLOT MARKERS ##
#######################

cat("\n\n")
cat(sprintf("%s FeaturePlots", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Plot clusters with highlighted error segmented cells
error_segmented_cells_plot = ggplot(cell_metadata, aes(x = umap_1, y = umap_2)) +
  geom_scattermore(
    data = subset(cell_metadata, color == "grey"),
    aes(color = color),
    pointsize = 1
  ) +
  geom_scattermore(
    data = subset(cell_metadata, color == "red"),
    aes(color = color),
    pointsize = 1
  ) +
  scale_color_identity() +
  ggtitle("Segmentation errors")

# Iterate over segmentation error cluster markers
for (feature in rownames(segmentation_error_markers)){
  
  # Retrieve current gene expression
  cell_metadata$current_gene_expression = FetchData(xenium_data, vars = feature)[, 1]
  
  # Plot expression of current feature
  feature_plot = ggplot(cell_metadata, aes(x = umap_1, y = umap_2)) +
    geom_scattermore(
      aes(color = current_gene_expression),
      pointsize = 1
    ) +
    scale_color_viridis_c(option = "magma") +
    ggtitle(paste("Expression of", feature)) +
    theme_classic() + 
    theme(legend.position = "none")
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), feature))
  cat("\n\n")
  
  # Merge plots
  cat("\n\n")
  print(grid.arrange(error_segmented_cells_plot, feature_plot, ncol = 2))
  cat("\n\n")
  
}

##########################################
## 1.3. SUBSET SEGMENTATION ERROR CELLS ##
##########################################

# Subset segmentation error cells
segmentation_error_subset = subset(xenium_data, 
                                   subset = cell_type_lvl1 %in% c("segmentation.error"))

################################################################################
#'                     2. SEGMENTATION ERROR CELLS CLUSTERING
################################################################################

cat("\n\n")
cat(sprintf("%s **SEGMENTATION ERRORS SUBCLUSTERING**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("As explained earlier, this step aims to re-cluster cells initially identified 
as segmentation errors in order to recover as many valid cells as possible. This 
section provides a more detailed view of the re-clustering process, and the following 
outputs are displayed:

  - Elbow plot: Visualization used to determine the number of dimensions to retain 
from the Harmony-corrected embedding.
  - Re-clustering resolutions tested:
    * DimPlot: UMAP visualization of re-clustered cells.
    * TMA ID batch effect: UMAP visualization colored by TMA ID to assess whether 
integration has been successful and batch effects have been mitigated.
    * Bar plot: Cell proportions per cluster.
  - Annotation resolution: UMAP visualization of re-clustered cells at the resolution 
selected for manual annotation.
  - Annotations: UMAP visualization of manual annotations assigned to the re-clustered 
segmentation error cells.
  - Cell type distribution: Frequencies of annotated cell types within the re-clustered 
segmentation error population.")
cat("\n\n")

# Integrate and cluster cells from current condition (TMA)
segmentation_error_subset = wrapped_harmony_integration(seurat_object = segmentation_error_subset, 
                                                        active_assay = XENIUM_ASSAY_NAME, 
                                                        cluster_layer_prefix = SEGMENTATION_ERROR_SUBSET_LAYER_PREFIX, 
                                                        umap_layer_name = SEGMENTATION_ERROR_SUBSET_UMAP_LAYER_NAME,
                                                        dimension_principal_component_analysis = SEGMENTATION_ERROR_SUBSET_DIM_PCA,
                                                        dimension_clustering = SEGMENTATION_ERROR_SUBSET_DIM_USE,
                                                        resolution_range = SEGMENTATION_ERROR_SUBSET_RES_RANGE,
                                                        harmony_integration_covariable = HARMONY_INTEGRATION_COVARIABLE,
                                                        reduction = HARMONY_ASSAY_NAME)

cat("\n\n")
cat(sprintf("%s Elbow plot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(get_elbow_plot(seurat_object = segmentation_error_subset, 
                     dimension_principal_component_analysis = SEGMENTATION_ERROR_SUBSET_DIM_PCA,
                     dimension_clustering = SEGMENTATION_ERROR_SUBSET_DIM_USE,
                     reduction = HARMONY_ASSAY_NAME))
cat("\n\n")

# Init markers list
markers_across_resolutions = list()

# Iterate over display resolutions
for (res in SEGMENTATION_ERROR_SUBSET_RES_RANGE){
  
  # Init cluster layer
  current_cluster_layer = sprintf("%s_%s", 
                                  SEGMENTATION_ERROR_SUBSET_LAYER_PREFIX, 
                                  gsub("\\.", "", as.character(res)))
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  ######################################
  ## 2.1. DISPLAY SEGMENTATION ERRORS ##
  ######################################
  
  # Activate current resolution
  segmentation_error_subset = activateSeuratAssay(segmentation_error_subset, 
                                                  ident = current_cluster_layer, 
                                                  active_assay = XENIUM_ASSAY_NAME, 
                                                  umap_layer_name = SEGMENTATION_ERROR_SUBSET_UMAP_LAYER_NAME)
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(segmentation_error_subset) +
          ggtitle(sprintf("Segmentation error clusters at resolution %s", res)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s TMA_ID batch effect", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(segmentation_error_subset,
                split.by = "TMA_id",
                group.by = "TMA_id") +
          theme(legend.position = "none") +
          facet_wrap(~TMA_id, nrow = round(length(unique(xenium_data$TMA_id))/2)) +
          theme(legend.position = "none"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Barplot", strrep("#", 3)))
  cat("\n\n")
  
  # Retrieve number of cells by clusters
  cluster_freq = as.data.frame(table(segmentation_error_subset$seurat_clusters))
  
  cat("\n\n")
  print(ggplot(cluster_freq, aes(x=Var1, 
                                 y=Freq)) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq, vjust = -0.3)) +
          theme(axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                legend.position = "none") +
          ylab("Cell number"))
  cat("\n\n")
  
  ##################
  ## 2.2. MARKERS ##
  ##################
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = segmentation_error_subset, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = current_cluster_layer,
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  # Store current markers
  markers_across_resolutions[[as.character(res)]] = markers
  
}

# Save markers in xlsx file
save_marker_genes_xlsx(marker_genes = markers_across_resolutions, 
                       output_directory = PATH_ANALYSIS_OUTPUT, 
                       output_file_name = "segmentation_error_subclustering_markers.xlsx")

##########################################################
## 2.3. MANUAL ANNOTATION OF SEGMENTATION ERRORS SUBSET ##
##########################################################

# Activate resolution at which segmentation errors are annotated
segmentation_error_subset = activateSeuratAssay(segmentation_error_subset, 
                                                ident = sprintf("%s_%s", 
                                                                SEGMENTATION_ERROR_SUBSET_LAYER_PREFIX, 
                                                                gsub("\\.", "", as.character(SEGMENTATION_ERROR_SUBSET_ANNOTATION_RES))), 
                                                active_assay = XENIUM_ASSAY_NAME, 
                                                umap_layer_name = SEGMENTATION_ERROR_SUBSET_UMAP_LAYER_NAME)

cat("\n\n")
cat(sprintf("%s Annotation res", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(segmentation_error_subset) +
        ggtitle(sprintf("Segmentation error clusters at resolution %s", 
                        SEGMENTATION_ERROR_SUBSET_ANNOTATION_RES)))
cat("\n\n")

# Initialize cell annotation meta data
cell_type_metadata = data.frame(cluster = segmentation_error_subset$seurat_clusters,
                                cell_type = "segmentation.error")

# Annotate cell type metadata
for (cluster in unique(cell_type_metadata$cluster)){
  
  # Check if current cluster corresponds to a cell type
  if (cluster %in% names(SEGMENTATION_ERROR_SUBSET_TYPING_METADATA)){
    
    cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = SEGMENTATION_ERROR_SUBSET_TYPING_METADATA[[as.character(cluster)]]
    
  } else {
    
    cell_type_metadata[cell_type_metadata$cluster == cluster, ]$cell_type = "segmentation.error"
    
  }
}

# Add cell typing to meta data
segmentation_error_subset[[LEVEL1_ANNOTATION_LAYER_NAME]] = cell_type_metadata$cell_type

cat("\n\n")
cat(sprintf("%s Annotations", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(segmentation_error_subset,
              group.by = LEVEL1_ANNOTATION_LAYER_NAME,
              cols = LEVEL_1_CELL_TYPE_COLORS))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Cell types distribution", strrep("#", 2)))
cat("\n\n")

# Prepare data
cell_type_distribution_data = as.data.frame(table(segmentation_error_subset[[LEVEL1_ANNOTATION_LAYER_NAME]]))

cat("\n\n")
print(ggplot(cell_type_distribution_data, aes(x=.data[[LEVEL1_ANNOTATION_LAYER_NAME]], 
                                              y=Freq,
                                              fill = .data[[LEVEL1_ANNOTATION_LAYER_NAME]])) + 
        geom_bar(stat = "identity") +
        geom_text(aes(label = Freq, vjust = -0.3)) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              legend.position = "none") +
        scale_fill_manual(values = LEVEL_1_CELL_TYPE_COLORS) +
        ylab("Cell number"))
cat("\n\n")

#############################################
## 2.4. SAVE ANNOTATED SEGMENTATION ERRORS ##
#############################################

# Save subset
saveRDS(segmentation_error_subset,
        file = file.path(PATH_ANALYSIS_OUTPUT,
                         "annotated_segmentation_error.Rds"))

################################################################################
##               3. MERGE LVL1 AND SEGMENTATION ERRORS ANNOTATIONS            ##
################################################################################

cat("\n\n")
cat(sprintf("%s **MERGE LVL1 AND SEGMENTATION ERRORS ANNOTATIONS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to update the full dataset by incorporating the revised manual 
annotations for a subset of cells previously classified as segmentation errors. 
Cells that could be confidently re-annotated are reassigned to the appropriate 
Level 1 cell types, thereby refining the overall annotation and reducing the number 
of cells labeled as segmentation errors.
    
The following sections are displayed to visualize the updated annotations:

  - Annotations before correction: Level 1 manual annotations prior to reassigning 
cells previously labeled as segmentation errors.
  - Annotations after correction: Level 1 manual annotations after incorporating 
the updated annotations for the corrected cells.
  - Sample IDs: Visualization of each level 1 cell type on the sample tissues.")
cat("\n\n")

#############################
## 3.1. RETRIEVE META DATA ##
#############################

# Retrieve level 1 cell annotations from global GAINED and segmentation cells
global_cell_type_lvl1_annotation = data.frame(cell = Cells(xenium_data),
                                              initial_cell_type = xenium_data$cell_type_lvl1)
segmentation_error_cell_type_lvl1_annotation = data.frame(cell = Cells(segmentation_error_subset),
                                                          updated_cell_type = segmentation_error_subset$cell_type_lvl1)

################################
## 3.2. DISPLAY INITIAL STATE ##
################################

cat("\n\n")
cat(sprintf("%s Annotations before merging", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = LEVEL1_ANNOTATION_LAYER_NAME,
              cols = LEVEL_1_CELL_TYPE_COLORS))
cat("\n\n")

############################
## 3.3. MERGE ANNOTATIONS ##
############################

# Merge annotations and update annotations
merged_annotations = left_join(global_cell_type_lvl1_annotation, segmentation_error_cell_type_lvl1_annotation, by = "cell")
merged_annotations[is.na(merged_annotations$updated_cell_type), ]$updated_cell_type = merged_annotations[is.na(merged_annotations$updated_cell_type), ]$initial_cell_type

# Update global cell type annotations
xenium_data[[LEVEL1_ANNOTATION_LAYER_NAME]] = merged_annotations$updated_cell_type

cat("\n\n")
cat(sprintf("%s Annotations after merging", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = LEVEL1_ANNOTATION_LAYER_NAME,
              cols = LEVEL_1_CELL_TYPE_COLORS))
cat("\n\n")

#######################################
## 3.4. DISPLAY CELL TYPES ON TISSUE ##
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
                                      ident = LEVEL1_ANNOTATION_LAYER_NAME, 
                                      active_assay = XENIUM_ASSAY_NAME)
  
  # Iterate over lvl1 cell types in the data
  for (cell_type in unique(xenium_data@meta.data[[LEVEL1_ANNOTATION_LAYER_NAME]])){
    
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
