
## @knitr cluster_TMA_cells

cat("\n\n")
cat(sprintf("%s **TMA clusters**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to visualize cell clusters in the data prior to applying any 
quality control filters. It allows for a quick visual comparison between the 
refined segmentation obtained with Proseg and the initial Xenium segmentation.")
cat("\n\n")

################################################################################
#'                              1. CLUSTER TMAs CELLS
#'                              
#' For each TMA and segmentation method, cluster cells using predefined 
#' dimensionality and clustering resolution parameters. 
#'                            
################################################################################

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Log normalize data
    xenium_data[[TMA]][[segmentation_method]] = log_normalisation(xenium_data[[TMA]][[segmentation_method]], 
                                                                  active_assay = XENIUM_ASSAY_NAME)
    
    # Cluster cells from current segmentation method
    xenium_data[[TMA]][[segmentation_method]] = wrapped_clustering_stomics_imaging(seurat_object = xenium_data[[TMA]][[segmentation_method]], 
                                                                                   active_assay = XENIUM_ASSAY_NAME, 
                                                                                   cluster_layer_prefix = CLUSTER_LAYER_PREFIX, 
                                                                                   umap_layer_name = UMAP_LAYER_NAME,
                                                                                   dimension_principal_component_analysis = DIM_PCA,
                                                                                   dimension_clustering = DIM_USE,
                                                                                   resolution_range = RES_RANGE)
    
  }
}

## @knitr plot_TMA_clusters

################################################################################
#'                               2. PLOT TMAs CLUSTERS
################################################################################

cat("\n\n")
cat("For each TMA, segmentation method, and clustering resolution:

  - DimPlot: UMAP visualization of clustered cells.
  - ElbowPlot: Visualization used to determine the optimal number of dimensions.
  - TMA: Visualization of clusters across the tissue (global TMA view).")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s ", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s ", strrep("#", 3), segmentation_method), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Iterate over resolution that were used
    for (res in RES_RANGE){
      
      cat("\n\n")
      cat(sprintf("%s %s ", strrep("#", 4), res), " {.tabset .tabset-fade}\n\n")
      cat("\n\n")
      
      ############################################
      ## 2.1. ACTIVATE CURRENT ASSAY CONDITIONS ##
      ############################################
      
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
      
      ###################
      ## 2.2. CLUSTERS ##
      ###################
      
      cat("\n\n")
      cat(sprintf("%s DimPlot ", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      print(DimPlot(xenium_data[[TMA]][[segmentation_method]],
                    label = TRUE) +
              ggtitle(sprintf("TMA clusters with %s segmentation (TMA = %s, res = %s)", segmentation_method, TMA, res)))
      cat("\n\n")
      
      #####################
      ## 2.3. ELBOW PLOT ##
      #####################
      
      cat("\n\n")
      cat(sprintf("%s Elbowplot ", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      print(get_elbow_plot(xenium_data[[TMA]][[segmentation_method]], 
                           dimension_principal_component_analysis = DIM_PCA,
                           dimension_clustering = DIM_USE))
      cat("\n\n")
      
      ##############
      ## 2.4. TMA ##
      ##############
      
      cat("\n\n")
      cat(sprintf("%s TMA ", strrep("#", 5)))
      cat("\n\n")
      
      cat("\n\n")
      print(ImageDimPlot(xenium_data[[TMA]][[segmentation_method]], 
                         size = 0.75) + 
              ggtitle(sprintf("TMA clusters with %s segmentation (TMA = %s, res = %s)", segmentation_method, TMA, res)) +
              labs(fill = sprintf("Resolution = %s", res)))
      cat("\n\n")
      
    }
  }
}

################################################################################
#'                              3. PLOT QC ON CLUSTERS
################################################################################

cat("\n\n")
cat(sprintf("%s **Quality control over TMA clusters**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("For each TMA and segmentation method, display cell quality control results 
on the UMAP. This allows assessment of whether low-quality cells cluster together 
or are distributed across clusters.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s ", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s ", strrep("#", 3), segmentation_method))
    cat("\n\n")
    
    cat("\n\n")
    print(DimPlot(xenium_data[[TMA]][[segmentation_method]],
                  cols = QUALITY_CONTROL_COLOR_PALETTE,
                  size = 0.75,
                  group.by = "QC_result") +
            ggtitle(sprintf("UMAP-projected QC results of %s (TMA = %s)", segmentation_method, TMA)))
    cat("\n\n")
    
  }
}
