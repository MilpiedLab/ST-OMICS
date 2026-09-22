
## @knitr addmodulescore

################################################################################
#'                                  MARKER GENES
#'                                   
#' This step aims to identify which markers are associated with the clusters 
#' obtained in the previous step. It is performed prior to any quality control 
#' filtering; therefore, the results and cluster assignments may change in later 
#' steps after quality control is applied.           
#'                        
################################################################################

################################################################################
#'                            1. ADDMODULESCORE ON UMAP
################################################################################

cat("\n\n")
cat(sprintf("%s **Cell types AddModuleScore**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Since a reference single-cell dataset is available for this spatial dataset, 
AddModuleScore can be used to evaluate the enrichment of gene signatures for each 
broad cell type found in the reference. These signatures are based on mutually 
exclusive gene sets that were computed in '01_Find_Mutually_Exclusive_Gene_Pairs'.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), segmentation_method), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Iterate over all cell types for which markers are available
    for (cell_type in rownames(markers_as_strings)){
      
      ############################################
      ## 1.1. RETRIEVE CELL TYPE GENE SIGNATURE ##
      ############################################
      
      # Keep markers of current broad cell type
      cell_type_genes = markers_as_strings[cell_type, ]
      cell_type_genes = strsplit(cell_type_genes, split = ",")
      cell_type_genes = cell_type_genes[[1]]
      
      #########################
      ## 1.2. ADDMODULESCORE ##
      #########################
      
      # AddModuleScore for gene signature of current cell type
      temp_seurat = AddModuleScore(xenium_data[[TMA]][[segmentation_method]], 
                                   features = list(cell_type_genes), 
                                   ctrl = ADD_MODULE_SCORE_CTRL)
      
      ######################
      ## 1.2. PLOT SCORES ##
      ######################
      
      cat("\n\n")
      cat(sprintf("%s %s ", strrep("#", 4), cell_type))
      cat("\n\n")
      
      cat("\n\n")
      print(FeaturePlot(temp_seurat, 
                        features = "Cluster1", 
                        cols = c("white", "blue4")) +
              ggtitle(sprintf("%s AddModuleScore", cell_type)))
      cat("\n\n")
      
    }
  }
}

## @knitr findallmarkers

################################################################################
#'                               2. FIND ALL MARKERS
################################################################################

cat("\n\n")
cat(sprintf("%s **Marker genes (find_all_markers)**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("For each TMA, segmentation method, and clustering resolution, FindAllMarkers 
is used to identify marker genes for each cluster.

  - DimPlot: UMAP visualization of clustered cells.
  - Cluster index:
    * Features: UMAP-based visualization of the top 9 marker genes per cluster, selected based on log fold-change expression.
    * Violin plots: Display the expression distribution of the top 9 marker genes across clusters.")
cat("\n\n")

# Store markers in list to save them later
markers_across_resolutions = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), segmentation_method), " {.tabset .tabset-fade}\n\n")
    cat("\n\n") 
    
    # Iterate over resolution that were used
    for (res in RES_RANGE){
      
      cat("\n\n")
      cat(sprintf("%s %s", strrep("#", 4), res), " {.tabset .tabset-fade}\n\n")
      cat("\n\n")
      
      # Compute markers for current condition
      markers = find_all_markers(seurat_object = xenium_data[[TMA]][[segmentation_method]], 
                                 active_assay = XENIUM_ASSAY_NAME, 
                                 cluster_layer = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                 logfc_threshold = LOGFC_THRESHOLD, 
                                 min_pct = MIN_PCT, 
                                 only_pos = TRUE,
                                 n_top_markers = NULL)
      
      # Store current markers
      markers_across_resolutions[[TMA]][[segmentation_method]][[as.character(res)]] = markers
      
      # Activate current assay (required resolution)
      xenium_data[[TMA]][[segmentation_method]] = activateSeuratAssay(xenium_data[[TMA]][[segmentation_method]], 
                                                                      ident = sprintf("%s_%s", 
                                                                                      CLUSTER_LAYER_PREFIX, 
                                                                                      gsub("\\.", "", as.character(res))), 
                                                                      active_assay = XENIUM_ASSAY_NAME, 
                                                                      umap_layer_name = UMAP_LAYER_NAME)
      
      # Get current cluster identities
      clusters = Idents(xenium_data[[TMA]][[segmentation_method]])
      
      # Convert identities to a factor ordered numerically
      Idents(xenium_data[[TMA]][[segmentation_method]]) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
      
      cat("\n\n")
      cat(sprintf("%s DimPlot ", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      print(DimPlot(xenium_data[[TMA]][[segmentation_method]],
                    label = TRUE) +
              ggtitle(sprintf("TMA clusters with %s segmentation (res = %s)", segmentation_method, res)))
      cat("\n\n")
      
      # Iterate over clusters to generate corresponding plots
      for (cluster in sort(as.numeric(unique(markers$cluster)))){
        
        cat("\n\n")
        cat(sprintf("%s %s", strrep("#", 5), as.numeric(cluster) - 1), " {.tabset .tabset-fade}\n\n")
        cat("\n\n")
        
        cat("\n\n")
        cat(sprintf("%s Features ", strrep("#", 6)))
        cat("\n\n")
        
        temp_features = get_marker_genes_feature_plots(seurat_object = xenium_data[[TMA]][[segmentation_method]], 
                                                       all_marker_genes = markers,
                                                       cluster_layer_name = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                                       cluster_number = (as.numeric(cluster) - 1),
                                                       top_n_markers = TOP_N_MARKERS)
        
        cat("\n\n")
        grid.arrange(grobs=temp_features, top = sprintf("%s marker genes expression", as.numeric(cluster) - 1))
        cat("\n\n")
        
        cat("\n\n")
        cat(sprintf("%s Violins ", strrep("#", 6)))
        cat("\n\n")
        
        temp_violins = get_marker_genes_violins(seurat_object = xenium_data[[TMA]][[segmentation_method]], 
                                                all_marker_genes = markers,
                                                cluster_layer_name = paste0(CLUSTER_LAYER_PREFIX, "_", gsub("\\.", "", as.character(res))),
                                                cluster_number = (as.numeric(cluster) - 1),
                                                top_n_markers = TOP_N_MARKERS)
        
        cat("\n\n")
        grid.arrange(grobs=temp_violins, top = sprintf("Cluster %s marker genes expression", as.numeric(cluster) - 1))
        cat("\n\n")
        
      }
    }
    
    # Create output directory
    if (!(dir.exists(PATH_ANALYSIS_OUTPUT))){
      dir.create(PATH_ANALYSIS_OUTPUT)
    }
    
    # Save markers in xlsx file
    save_marker_genes_xlsx(marker_genes = markers_across_resolutions[[TMA]][[segmentation_method]], 
                           output_directory = PATH_ANALYSIS_OUTPUT, 
                           output_file_name = paste0(TMA, "_", segmentation_method, "_", "markers.xlsx"))
    
  }
}
