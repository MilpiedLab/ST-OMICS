
## @knitr singler_lvl1_annotation

cat("\n\n")
cat(sprintf("%s **SINGLER LVL1 CLUSTERING**", strrep("#", 1)))
cat("\n\n")

cat("\n\n")
cat("This step aims to automatically annotate cells using SingleR. As SingleR 
operates at the single-cell level, cells within the same cluster may be assigned 
to different, and sometimes inconsistent, cell types.

To obtain consistent cluster-level annotations, each cluster is labeled based on 
the most frequent SingleR annotation among its constituent cells (Level 1). In 
addition, SingleR is run using a reference single-cell dataset restricted to 
Level 1 cell type annotations, ensuring that the predicted labels correspond to 
broad cell populations.

Although SingleR can directly provide fine-grained (e.g., Level 2) annotations 
at the single-cell level, this analysis adopts an alternative strategy. Because 
annotations are harmonized at the cluster level, clusters are subsequently subsetted 
and re-clustered to define new clusters capturing finer levels of heterogeneity, 
following an approach similar to the manual annotation described in Section 
08a_Manual_Cell_Annotation.

Based on this refined clustering, more specific cell types (Level 2) can then be 
assigned. These annotations are obtained using scType, with cluster labels determined 
by the most frequent predicted cell type within each refined cluster.")
cat("\n\n")

################################################################################
#'                           1. SINGLER LVL1 CLUSTERING
################################################################################

cat("\n\n")
cat("In this section, a DimPlot is displayed showing the clusters obtained at a 
high resolution, in order to generate a large number of clusters from the full 
Harmony-integrated dataset.")
cat("\n\n")

# Find clusters at high resolutions to obtain many clusters to annotate
xenium_data = FindClusters(xenium_data,
                           resolution = SINGLER_LEVEL1_RESOLUTION, 
                           verbose = TRUE,
                           cluster.name = sprintf("%s_%s", 
                                                  HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                                  gsub("\\.", "", as.character(SINGLER_LEVEL1_RESOLUTION))))

cat("\n\n")
cat(sprintf("%s DimPlot", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data))
cat("\n\n")

################################################################################
#'                           2. SINGLER LVL1 ANNOTATION
################################################################################

cat("\n\n")
cat(sprintf("%s **SINGLER LVL1 ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to visualize the Level 1 SingleR annotations using the 
following sections:

  - Raw SingleR Level 1 annotation: UMAP visualization of the raw SingleR cell-level 
annotations.
  - Cluster annotation refinement: For each cluster identified in the previous section:
    * DimPlot: UMAP visualization highlighting the current cluster.
    * Cell type frequencies: Bar plot showing the distribution of SingleR cell 
annotations within the current cluster.
  - Refined SingleR Level 1 annotation:
    * Global view: UMAP visualization of the harmonized (cluster-level) SingleR 
annotations.
    * Sample IDs: Visualization of each singleR Level 1 annotations on the sample 
tissues.")
cat("\n\n")

######################################
## 2.1. SINGLER RAW LVL1 ANNOTATION ##
######################################

# SingleR level1 cell annotation
singler_lvl1_annotations = singleRWrapper(seurat_object = xenium_data,
                                          singler_reference = sc_reference_data,
                                          active_assay = XENIUM_ASSAY_NAME,
                                          cluster_layer_name = sprintf("%s_%s", 
                                                                       HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                                                       gsub("\\.", "", as.character(SINGLER_LEVEL1_RESOLUTION))),
                                          annotation_layer_name = "Celltype_Level_1")

# Update xenium data with singleR annotations
xenium_data@meta.data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]] = singler_lvl1_annotations[["singler_results"]]@listData[["labels"]]

cat("\n\n")
cat(sprintf("%s Raw singleR lvl1 annotation", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = SINGLER_LVL1_ANNOTATIONS_LAYER_NAME))
cat("\n\n")

#############################################
## 2.2. SINGLER LVL1 ANNOTATION REFINEMENT ##
#############################################

# Activate current condition
xenium_data = activateSeuratAssay(xenium_data, 
                                  ident = sprintf("%s_%s", 
                                                  HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                                  gsub("\\.", "", as.character(SINGLER_LEVEL1_RESOLUTION))), 
                                  active_assay = XENIUM_ASSAY_NAME, 
                                  umap_layer_name = HARMONY_UMAP_LAYER_NAME)

cat("\n\n")
cat(sprintf("%s Cluster annotation refinement", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Retrieve raw lvl1 singleR annotation
singleR_lvl1_raw_annotation = data.frame("cell_id" = Cells(xenium_data),
                                         "cluster" = xenium_data$seurat_clusters,
                                         "cell_type" = xenium_data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]])
singleR_lvl1_refined_annotation = singleR_lvl1_raw_annotation

# Retrieve UMAP data to plot
complete_plot_data = data.frame(cell = Cells(xenium_data),
                                umap_1 = xenium_data@reductions[["umap"]]@cell.embeddings[, 1],
                                umap_2 = xenium_data@reductions[["umap"]]@cell.embeddings[, 2],
                                cluster = xenium_data$seurat_clusters,
                                color = "grey")

# Iterate over clusters to make singleR annotation
for (cluster in sort(as.numeric(unique(xenium_data@meta.data[["seurat_clusters"]]))-1)){
  
  # Retrieve meta data relative to current cluster
  cluster_raw_annotations = singleR_lvl1_raw_annotation[singleR_lvl1_raw_annotation$cluster == cluster, ]
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), cluster), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 4)))
  cat("\n\n")
  
  # Update color for plot
  complete_plot_data$color = ifelse(complete_plot_data$cluster == cluster, "red", "grey")
  
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
          ggtitle(sprintf("Cluster %s", cluster)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Cell type frequencies", strrep("#", 4)))
  cat("\n\n")
  
  # Retrieve number of cells by clusters
  cluster_freq = as.data.frame(table(cluster_raw_annotations[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]]))
  
  cat("\n\n")
  print(ggplot(cluster_freq, aes(x=Var1, 
                                 y=Freq)) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq, vjust = -0.3)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                legend.position = "none") +
          ylab("Cell count") +
          xlab("Cell types") +
          ggtitle(sprintf("Cluster %s cell type frequencies", cluster)))
  cat("\n\n")
  
  # Retrieve most frequent annotation in the cluster
  cluster_most_frequent_annotation = names(which.max(table(cluster_raw_annotations[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]])))
  
  # Update annotations for current cluster
  singleR_lvl1_refined_annotation[singleR_lvl1_refined_annotation$cluster == cluster, ]$singleR_lvl1_cell_type = cluster_most_frequent_annotation
  
}

# Update xenium data with final refined annotations
xenium_data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]] = singleR_lvl1_refined_annotation$singleR_lvl1_cell_type

cat("\n\n")
cat(sprintf("%s Refined singleR lvl1 annotation", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

#############################################
## 2.3. PLOT FINAL LVL1 SINGLER ANNOTATION ##
#############################################

cat("\n\n")
cat(sprintf("%s Global view", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = SINGLER_LVL1_ANNOTATIONS_LAYER_NAME))
cat("\n\n")

#################################################
## 2.4. PLOT SINGLER LVL1 ANNOTATION ON TISSUE ##
#################################################

# Iterate over samples in data
for (sample in unique(xenium_data@meta.data[["SAMPLE_ID"]])){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 3), sample), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Subset current sample from global merged data
  sample_subset = subset(xenium_data,
                         subset = SAMPLE_ID == sample)
  
  # Retrieve correspondance between cells and singleR annotation
  all_cells_singler_annotation = data.frame(cell_id = Cells(sample_subset),
                                            cell_type = sample_subset@meta.data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]])
  
  # Iterate over anootated cell types for current condition
  for (singler_cell_type in unique(all_cells_singler_annotation$cell_type)){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 4), singler_cell_type))
    cat("\n\n")
    
    # Retrieve cells with current annotation
    cells_of_interest = all_cells_singler_annotation[all_cells_singler_annotation$cell_type == singler_cell_type, ]$cell_id
    
    # Subset current cell type
    cell_type_subset = subset(sample_subset, cells = cells_of_interest)
    
    # Display cell type on tissue
    cat("\n\n")
    print(ImageDimPlot(cell_type_subset,
                       group.by = SINGLER_LVL1_ANNOTATIONS_LAYER_NAME,
                       size = 0.75) + 
            ggtitle(sprintf("%s on tissue", singler_cell_type)) +
            labs(fill = sprintf("cell type = %s", singler_cell_type)))
    cat("\n\n")
    
  }
}

## @knitr singler_lvl2_annotation

################################################################################
#'                           3. SINGLER LVL2 ANNOTATION
################################################################################

cat("\n\n")
cat(sprintf("%s **SINGLER LVL2 ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("As explained earlier, clusters sharing the same refined Level 1 SingleR 
annotation are subsetted together and re-clustered to define new clusters capturing 
finer levels of heterogeneity. These finer clusters are then annotated using SingleR 
following the same procedure as for Level 1 clusters.

This step aims to visualize the Level 2 SingleR annotations using the following 
sections:

  - Level 1 annotation (per subset of refined SingleR Level 1 clusters):
    * DimPlot: UMAP visualization of clustered cells.
    * Marker genes: Table reporting the results from FindAllMarkers for all 
clusters.
    * Elbow plot: Visualization used to determine the number of dimensions to 
retain from the Harmony-corrected embedding.
    * Raw SingleR Level 2 annotation: UMAP visualization of raw cell-level SingleR annotations.
    * Cluster annotation refinement: For each cluster identified within the current subset:
      - DimPlot: UMAP visualization highlighting the current cluster.
      - Cell type frequencies: Bar plot showing the distribution of SingleR 
annotations within the cluster.
    * Refined SingleR Level 2 annotation: UMAP visualization of harmonized 
(cluster-level) SingleR annotations.
    * Sample IDs: Visualization of Level 2 SingleR annotations mapped onto 
tissue samples for the current subset.
  - Refined global Level 2 SingleR annotation: UMAP visualization of harmonized 
Level 2 SingleR annotations for the full merged dataset.")
cat("\n\n")

# Init annotation for lvl2 cell types
singleR_lvl2_refined_annotation = data.frame("cell_id" = Cells(xenium_data),
                                             "singleR_lvl1_cell_type" = xenium_data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]],
                                             "singleR_lvl2_cell_type" = "Unknown")

# Iterate over level1 refined cell obtained with singleR
for (lvl1_cell_type in unique(xenium_data@meta.data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]])){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), lvl1_cell_type), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  ################################
  ## 3.1. SUBSET AND CLUSTERING ##
  ################################
  
  # Subset data specific to current cell type
  subset_xenium = xenium_data[, xenium_data[[SINGLER_LVL1_ANNOTATIONS_LAYER_NAME]] == lvl1_cell_type]
  
  # Integrate and cluster cells from current condition (TMA)
  subset_xenium = wrapped_harmony_integration(seurat_object = subset_xenium, 
                                              active_assay = XENIUM_ASSAY_NAME, 
                                              cluster_layer_prefix = SINGLER_LVL2_CLUSTER_LAYER_PREFIX, 
                                              umap_layer_name = SINGLER_LVL2_UMAP_LAYER_NAME,
                                              dimension_principal_component_analysis = SINGLER_LVL2_DIM_PCA[[lvl1_cell_type]],
                                              dimension_clustering = SINGLER_LVL2_DIM_USE[[lvl1_cell_type]],
                                              resolution_range = SINGLER_LVL2_RESOLUTION[[lvl1_cell_type]],
                                              harmony_integration_covariable = HARMONY_INTEGRATION_COVARIABLE,
                                              reduction = HARMONY_ASSAY_NAME)
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(subset_xenium) +
          ggtitle(sprintf("Cluster at res = %s", SINGLER_LVL2_RESOLUTION[[lvl1_cell_type]])))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Marker genes", strrep("#", 3)))
  cat("\n\n")
  
  # Compute markers for current condition
  markers = find_all_markers(seurat_object = subset_xenium, 
                             active_assay = XENIUM_ASSAY_NAME, 
                             cluster_layer = sprintf("%s_%s", 
                                                     SINGLER_LVL2_CLUSTER_LAYER_PREFIX, 
                                                     gsub("\\.", "", as.character(SINGLER_LVL2_RESOLUTION[[lvl1_cell_type]]))),
                             logfc_threshold = LOGFC_THRESHOLD, 
                             min_pct = MIN_PCT, 
                             only_pos = TRUE,
                             n_top_markers = NULL)
  
  cat("\n\n")
  print(htmltools::tagList(DT::datatable(markers)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Elbowplot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(get_elbow_plot(seurat_object = subset_xenium, 
                       dimension_principal_component_analysis = SINGLER_LVL2_DIM_PCA[[lvl1_cell_type]],
                       dimension_clustering = SINGLER_LVL2_DIM_USE[[lvl1_cell_type]]))
  cat("\n\n")
  
  ######################################
  ## 3.2. SINGLER RAW LVL2 ANNOTATION ##
  ######################################
  
  # Subset reference data to focus on current lvl1 cell type
  subset_reference = subset(sc_reference_data, subset = Celltype_Level_1 == lvl1_cell_type)
  
  # SingleR level2 cell annotation
  singler_lvl2_annotations = singleRWrapper(seurat_object = subset_xenium,
                                            singler_reference = subset_reference,
                                            active_assay = XENIUM_ASSAY_NAME,
                                            cluster_layer_name = sprintf("%s_%s", 
                                                                         SINGLER_LVL2_CLUSTER_LAYER_PREFIX, 
                                                                         gsub("\\.", "", as.character(SINGLER_LVL2_RESOLUTION[[lvl1_cell_type]]))),
                                            annotation_layer_name = "Celltype_Level_2")
  
  # Update xenium data with singleR annotations
  subset_xenium@meta.data[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]] = singler_lvl2_annotations[["singler_results"]]@listData[["labels"]]
  
  cat("\n\n")
  cat(sprintf("%s Raw singleR lvl2 annotation", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(subset_xenium,
                group.by = SINGLER_LVL2_ANNOTATIONS_LAYER_NAME))
  cat("\n\n")
  
  #############################################
  ## 3.3. SINGLER LVL2 ANNOTATION REFINEMENT ##
  #############################################
  
  # Activate current condition
  subset_xenium = activateSeuratAssay(subset_xenium, 
                                      ident = sprintf("%s_%s", 
                                                      SINGLER_LVL2_CLUSTER_LAYER_PREFIX, 
                                                      gsub("\\.", "", as.character(SINGLER_LVL2_RESOLUTION[[lvl1_cell_type]]))), 
                                      active_assay = XENIUM_ASSAY_NAME, 
                                      umap_layer_name = SINGLER_LVL2_UMAP_LAYER_NAME)
  
  cat("\n\n")
  cat(sprintf("%s Cluster annotation refinement", strrep("#", 3)), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Retrieve raw lvl2 singleR annotation
  current_cell_type_raw_annotation = data.frame("cell_id" = Cells(subset_xenium),
                                                "cluster" = subset_xenium$seurat_clusters,
                                                "cell_type" = subset_xenium[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]])
  current_cell_type_refined_annotation = current_cell_type_raw_annotation
  
  # Retrieve UMAP data to plot
  complete_plot_data = data.frame(cell = Cells(subset_xenium),
                                  umap_1 = subset_xenium@reductions[["umap"]]@cell.embeddings[, 1],
                                  umap_2 = subset_xenium@reductions[["umap"]]@cell.embeddings[, 2],
                                  cluster = subset_xenium$seurat_clusters,
                                  color = "grey")
  
  # Iterate over clusters to make singleR annotation
  for (cluster in sort(unique(subset_xenium@meta.data[["seurat_clusters"]]))){
    
    # Retrieve meta data relative to current cluster
    cluster_raw_annotations = current_cell_type_raw_annotation[current_cell_type_raw_annotation$cluster == cluster, ]
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 4), cluster), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s DimPlot", strrep("#", 5)))
    cat("\n\n")
    
    # Update color for plot
    complete_plot_data$color = ifelse(complete_plot_data$cluster == cluster, "red", "grey")
    
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
            ggtitle(sprintf("Cluster %s", cluster)))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Cell type frequencies", strrep("#", 5)))
    cat("\n\n")
    
    # Retrieve number of cells by clusters
    cluster_freq = as.data.frame(table(cluster_raw_annotations[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]]))
    
    cat("\n\n")
    print(ggplot(cluster_freq, aes(x=Var1, 
                                   y=Freq)) + 
            geom_bar(stat = "identity") +
            geom_text(aes(label = Freq, vjust = -0.3)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  legend.position = "none") +
            ylab("Cell count") +
            xlab("Cell types") +
            ggtitle(sprintf("Cluster %s cell type frequencies", cluster)))
    cat("\n\n")
    
    # Retrieve most frequent annotation in the cluster
    cluster_most_frequent_annotation = names(which.max(table(cluster_raw_annotations[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]])))
    
    # Update annotations for current cluster
    current_cell_type_refined_annotation[current_cell_type_refined_annotation$cluster == cluster, ]$singleR_lvl2_cell_type = cluster_most_frequent_annotation
    
  }
  
  # Update xenium data with final refined annotations
  subset_xenium[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]] = current_cell_type_refined_annotation$singleR_lvl2_cell_type
  
  cat("\n\n")
  cat(sprintf("%s Refined singleR lvl2 annotation", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(subset_xenium,
                group.by = SINGLER_LVL2_ANNOTATIONS_LAYER_NAME))
  cat("\n\n")
  
  # Iterate over samples in data
  for (sample in unique(subset_xenium@meta.data[["SAMPLE_ID"]])){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), sample), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Subset current sample from global merged data
    sample_subset = subset(subset_xenium,
                           subset = SAMPLE_ID == sample)
    
    # Retrieve correspondance between cells and singleR annotation
    all_cells_singler_annotation = data.frame(cell_id = Cells(sample_subset),
                                              cell_type = sample_subset@meta.data[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]])
    
    # Iterate over anootated cell types for current condition
    for (singler_cell_type in unique(all_cells_singler_annotation$cell_type)){
      
      cat("\n\n")
      cat(sprintf("%s %s", strrep("#", 4), singler_cell_type))
      cat("\n\n")
      
      # Retrieve cells with current annotation
      cells_of_interest = all_cells_singler_annotation[all_cells_singler_annotation$cell_type == singler_cell_type, ]$cell_id
      
      # Subset current cell type
      cell_type_subset = subset(sample_subset, cells = cells_of_interest)
      
      # Display cell type on tissue
      cat("\n\n")
      print(ImageDimPlot(cell_type_subset,
                         group.by = SINGLER_LVL2_ANNOTATIONS_LAYER_NAME,
                         size = 0.75) + 
              ggtitle(sprintf("%s on tissue", singler_cell_type)) +
              labs(fill = sprintf("cell type = %s", singler_cell_type)))
      cat("\n\n")
      
    }
  }
  
  # Add annotations from current cell type to global dataset annotation metadata
  idx = match(singleR_lvl2_refined_annotation$cell_id, current_cell_type_refined_annotation$cell_id)
  singleR_lvl2_refined_annotation$singleR_lvl2_cell_type[!is.na(idx)] = current_cell_type_refined_annotation$singleR_lvl2_cell_type[idx[!is.na(idx)]]
  
}

# Update xenium data with final refined annotations
xenium_data[[SINGLER_LVL2_ANNOTATIONS_LAYER_NAME]] = singleR_lvl2_refined_annotation$singleR_lvl2_cell_type

cat("\n\n")
cat(sprintf("%s Refined global singleR lvl2 annotation", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(DimPlot(xenium_data,
              group.by = SINGLER_LVL2_ANNOTATIONS_LAYER_NAME))
cat("\n\n")
