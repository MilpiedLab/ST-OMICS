
## @knitr sctype_annotation

cat("\n\n")
cat(sprintf("%s **SCTYPE ANNOTATION**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to automatically annotate the Harmony-derived clusters obtained 
earlier using scType. For each clustering resolution tested, the following sections 
are displayed:

  - DimPlot: UMAP visualization of clustered cells.
  - scType annotation: Cluster-level annotations assigned by scType.
  - scType results table: Annotation scores provided by scType for each cluster at 
the current resolution.
  - Sample IDs: Visualization of each cell type identified by scTpe on the sample 
tissues.")
cat("\n\n")

################################################################################
#'                                SCTYPE ANNOTATION
################################################################################

# Iterate over resolutions to annotate
for (res in HARMONY_RESOLUTION_RANGE){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), res), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  ##########################
  ## 1. SCTYPE ANNOTATION ##
  ##########################
  
  # Define layers for current condition
  current_cluster_layer = sprintf("%s_%s", 
                                  HARMONY_CLUSTER_LAYER_NAME_PREFIX, 
                                  gsub("\\.", "", as.character(res)))
  current_sctype_annotation_layer = sprintf("%s_%s", 
                                            SCTYPE_LAYER_PREFIX, 
                                            gsub("\\.", "", as.character(res)))
  
  # scType analysis
  sctype_results = scTypeWrapper(seurat_object = xenium_data, 
                                 marker_genes_database = CELL_TYPE_MARKERS_DB_PATH, 
                                 tissue_type = TISSUE_TYPE, 
                                 active_assay = XENIUM_ASSAY_NAME,
                                 cluster_layer_name = current_cluster_layer,
                                 annotation_layer_name = current_sctype_annotation_layer)
  
  ###########################
  ## 2. DISPLAY ANNOTATION ##
  ###########################
  
  cat("\n\n")
  cat(sprintf("%s DimPlot", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(sctype_results[["seurat_object"]],
                group.by = current_cluster_layer) +
          ggtitle(sprintf("Merged Xenium data at res = %s", res)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s scType annotation", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(DimPlot(sctype_results[["seurat_object"]],
                group.by = current_sctype_annotation_layer) +
          ggtitle(sprintf("scType annotation at res = %s", res)))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s scType result table", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(htmltools::tagList(DT::datatable(sctype_results[["sctype_all_scores"]])))
  cat("\n\n")
  
  #####################################
  ## 3. DISPLAY CELL TYPES ON TISSUE ##
  #####################################
  
  # Iterate over samples in data
  for (sample in unique(sctype_results[["seurat_object"]]@meta.data[["SAMPLE_ID"]])){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), sample), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Subset current sample from global merged data
    sample_subset = subset(sctype_results[["seurat_object"]],
                           subset = SAMPLE_ID == sample)
    
    # Retrieve correspondance between cells and sctype annotation
    all_cells_sctype_annotation = data.frame(cell_id = Cells(sample_subset),
                                             cell_type = sample_subset@meta.data[[current_sctype_annotation_layer]])
    
    # Iterate over anootated cell types for current condition
    for (sctype_cell_type in unique(all_cells_sctype_annotation$cell_type)){
      
      cat("\n\n")
      cat(sprintf("%s %s", strrep("#", 4), sctype_cell_type))
      cat("\n\n")
      
      # Retrieve cells with current annotation
      cells_of_interest = all_cells_sctype_annotation[all_cells_sctype_annotation$cell_type == sctype_cell_type, ]$cell_id
      
      # Subset current cell type
      cell_type_subset = subset(sample_subset, cells = cells_of_interest)
      
      # Display cell type on tissue
      cat("\n\n")
      print(ImageDimPlot(cell_type_subset,
                         group.by = current_sctype_annotation_layer,
                         size = 0.75) + 
              ggtitle(sprintf("%s on tissue", sctype_cell_type)) +
              labs(fill = sprintf("cell type = %s", sctype_cell_type)))
      cat("\n\n")
      
    }
  }
  
  ############################
  ## 4. FINALIZE ANNOTATION ##
  ############################
  
  # Update xenium data with annotation
  xenium_data = sctype_results[["seurat_object"]]
  
}
