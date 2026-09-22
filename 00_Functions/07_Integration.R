################################################################################
##                          INTEGRATION RELATED FUNCTIONS                     ##
################################################################################

########################################
## HARMONY INTEGRATION AND CLUSTERING ##
########################################

#' This function wraps the steps necessary for the harmony integration,
#' dimension reduction and clustering analysis of imaging-based spatial 
#' transcriptomics data.
#'
#' Libraries: Seurat
#'
#' @param seurat_object Seurat object from which cells will be clustered.
#' @param active_assay Name of the assay to set as default assay.
#' @param cluster_layer_prefix Name prefix of the layer(s) that will contain 
#' clusters assignment and store in meta.data.
#' @param umap_layer_name Name of the UMAP coordinates layer that will be stored
#' in meta.data.
#' @param dimension_principal_component_analysis Number of computed PCs during
#' the "RunPCA" step.
#' @param dimension_clustering Number of PCs kept for the subsequent clustering
#' step.
#' @param resolution_range Range of resolutions (can be a single value) to
#' compute in the "findClusters" step of the analysis.
#'
#' @return Seurat object with data resulting from clustering step.
#'
wrapped_harmony_integration = function(seurat_object, 
                                       active_assay, 
                                       cluster_layer_prefix, 
                                       umap_layer_name,
                                       dimension_principal_component_analysis,
                                       dimension_clustering,
                                       resolution_range,
                                       harmony_integration_covariable,
                                       reduction = "harmony")
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(umap_layer_name))) stop("Provided umap_layer_name is not a character.")
  if (!(is.character(cluster_layer_prefix))) stop("Provided cluster_layer_prefix is not a character.")
  if (!(is.numeric(dimension_principal_component_analysis))) stop("Provided dimension_principal_component_analysis is not a numeric")
  if (!(is.numeric(dimension_clustering))) stop("Provided dimension_clustering is not a numeric")
  if (!(is.numeric(resolution_range))) stop("Provided resolution_range is not a numeric")
  
  # Activate required assay from which clusters will be computed
  seurat_object = activateSeuratAssay(seurat_object, active_assay = active_assay)
  
  ## PCA ##
  # Extract the list of genes that will be used for the PCA
  all.genes = rownames(seurat_object)
  # Scale and center the expression values for all genes
  seurat_object = ScaleData(seurat_object, features = all.genes)
  # Run PCA on scaled expression values for specified genes
  seurat_object = RunPCA(seurat_object, 
                         npcs = dimension_principal_component_analysis,
                         features = all.genes)
  
  # Run Harmony integration
  seurat_object = RunHarmony(seurat_object, 
                             harmony_integration_covariable, 
                             early_stop = FALSE)
  
  ## UMAP ##
  # Run UMAP to project selected PCs in 2D space
  seurat_object = RunUMAP(seurat_object,
                          reduction = reduction,
                          dims = 1:dimension_clustering, 
                          verbose = FALSE)
  # Save UMAP embeddings in meta.data
  seurat_object = storeUMAPembeddings(seurat_object = seurat_object, 
                                      umap_layer_name = umap_layer_name)
  
  ## Clustering ##
  # Constructs a K-nearest neighbors (KNN) graph based on the selected PCs
  seurat_object = FindNeighbors(seurat_object, 
                                reduction = reduction, 
                                dims = 1:dimension_clustering, 
                                verbose = FALSE)
  # Create output cluster layers names
  cluster_layers = c()
  for (res in resolution_range){
    cluster_layers = c(cluster_layers, sprintf("%s_%s", 
                                               cluster_layer_prefix, 
                                               gsub("\\.", "", as.character(res))))
  }
  # Applies graph-based clustering (Louvain algorithm) on the KNN graph
  seurat_object = FindClusters(seurat_object,
                               resolution = resolution_range, 
                               verbose = TRUE,
                               cluster.name = cluster_layers)
  
  return(seurat_object)
}
