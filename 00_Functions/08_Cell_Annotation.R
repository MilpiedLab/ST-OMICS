################################################################################
##                       CELL-ANNOTATION RELATED FUNCTIONS                    ##
################################################################################

##############################################
## RETRIEVE CORRESPONDING REFERENCE CLUSTER ##
##############################################

#' From a cluster and its markers in the Xenium data, find corresponding
#' cluster in a reference annotated single cell data using AddModuleScore.
compute_reference_clusters_modulescores = function(query_seurat,
                                                   reference_seurat,
                                                   query_markers_list,
                                                   CTRL = 25)
{
  
  # Init result list
  results = list()
  
  # Filter genes from module that are not found in reference single cell data
  query_markers_list = lapply(query_markers_list, function(signature) {
    signature[signature %in% rownames(reference_seurat)]
  })
  
  # Retrieve data to plot
  complete_plot_data = data.frame(cell = Cells(reference_seurat),
                                  umap_1 = reference_seurat@reductions[["umap"]]@cell.embeddings[, 1],
                                  umap_2 = reference_seurat@reductions[["umap"]]@cell.embeddings[, 2],
                                  celltype = reference_seurat$celltype,
                                  modulescore = AddModuleScore(reference_seurat, 
                                                               features = query_markers_list,
                                                               ctrl = CTRL)$Cluster1)
  
  # Display module score over reference cells
  results[["plot"]] = ggplot(complete_plot_data, aes(x = umap_1, y = umap_2, color = modulescore)) +
    geom_scattermore(pointsize = 1) +
    scale_color_gradientn(colors = c("blue", "grey", "red")) +
    ggtitle("module score")
  
  # Retrieve statistics by cell type
  celltype_frequency = as.data.frame(table(complete_plot_data$celltype))
  names(celltype_frequency) = c("Cell_type", "Frequency")
  rownames(celltype_frequency) = celltype_frequency$celltype
  
  # Compute mean module score by cell type
  mean_module_score_by_celltype = as.data.frame(tapply(complete_plot_data$modulescore, 
                                                       complete_plot_data$celltype, 
                                                       mean))
  names(mean_module_score_by_celltype) = c("AddModuleScore")
  results[["similarity_scores"]] = cbind(celltype_frequency,
                                         mean_module_score_by_celltype)
  
  return(results)
  
}

#####################
## SCTYPE ANALYSIS ##
#####################

#' This function wraps the steps necessary for ScType automatic cell annotation.
#' This is an adaptation of the code found in ScType GitHub.
#' ScType requires to have a cluster layer to be computed. It will then annotate
#' each of the cluser found in this cluster layer.
#' ScType is a two step process for its annotation:
#'    - For each cell in the dataset, assign an annotation score for each
#'    available cell type based on the list of positive and negative marker
#'    genes.
#'    - For each cluster, it will sum the score for each cell type across the
#'    cells from this cluster. The cell type with the highest score will be the
#'    cell type assigned to the cluster.
#' 
#' @param seurat_object Seurat object containing cells to annotate.
#' @param marker_genes_database xlsx file containing cell types as rows and 
#' their corresponding marker genes (positive and negative). This file is read 
#' as a dataframe by ScType and should contain the following columns:
#'    - "tissueType": Column with the name of the tissue of interest. This
#'    column is there in case there are multiple tissue type within the excel 
#'    file. ScType will only keep the rows (cell types) that have the required
#'    tissue type in this column.
#'    - "cellName": Column with the cell types of interest. Each cell type
#'    corresponds to a row in the file.
#'    - "geneSymbolmore1": List of positive marker genes for the corresponding
#'    cell type. The gene names are all written in the same cell with a comma 
#'    separation.
#'    - "geneSymbolmore2": List of negative marker genes for the corresponding
#'    cell type. The gene names are all written in the same cell with a comma 
#'    separation.
#'    - "shortName": Abbreviation of cell type names in the "cellName" column.
#' @param tissue_type Tissue type name corresponding to the first column of
#' the excel file database.
#' @param active_assay Name of the assay to set as default assay.
#' @param cluster_layer_name Name of the cluster layer store in meta.data and
#' containing the clusters to annotate.
#' @param annotation_layer_name Name of the meta.data layer that will be created
#' upon annotation with ScType. It will contain the resulting annotations for
#' each cells in the seurat object.
#' @param color_palette Color palette to use for the bubble plot. 
#' 
#' @return A list with the following elements:
#'    - seurat_object: Updated Seurat object that was provided which now contains
#'    a layer in meta data with ScType annotations for each cell.
#'    - sctype_all_scores: A dataframe with top 10 cell types and corresponding
#'    scores obtained by ScType for each cluster.
#'    
scTypeWrapper = function(seurat_object, 
                         marker_genes_database, 
                         tissue_type, 
                         active_assay,
                         cluster_layer_name,
                         annotation_layer_name)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!file.exists(marker_genes_database)) stop("marker_genes_database file does not exist.")
  if (!(is.character(tissue_type))) stop("Provided tissue_type is not a character.")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.character(cluster_layer_name))) stop("Provided cluster_layer_name is not a character.")
  if (!(is.character(annotation_layer_name))) stop("Provided annotation_layer_name is not a character.")
  if (is.null(seurat_object[[active_assay]]$scale.data)) stop("Provided seurat object does not have a scaled data layer.")
  
  # Activate current assay (get back assay name and required initial clusters)
  seurat_object = activateSeuratAssay(seurat_object,
                                      ident = cluster_layer_name,
                                      active_assay = active_assay)
  
  # Retrieve list of positive and negative marker genes for each cell type in the required tissue from provided database
  gene_set = gene_sets_prepare(marker_genes_database, tissue_type)
  
  # Extract scaled scRNA-seq matrix
  scaled_data = as.matrix(seurat_object[[active_assay]]$scale.data)
  
  # Calculate ScType scores and assign cell types - Main matrix of ScType scores for all cell types and cells
  # It is a complete score matrix of size (cell types × cells)
  enrichment_scores = sctype_score(scRNAseqData = scaled_data, 
                                   scaled = TRUE, 
                                   gs = gene_set$gs_positive, 
                                   gs2 = gene_set$gs_negative)
  
  # Merge score by cluster (for each cluster, get an aggregated score for each required cell type)
  cluster_label_results = do.call("rbind", 
                                  lapply(unique(seurat_object@meta.data$seurat_clusters), 
                                         function(cluster_number){
                                           enrichment_scores_cluster_number = sort(rowSums(enrichment_scores[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cluster_number, ])]), decreasing = !0)
                                           head(data.frame(cluster = cluster_number, type = names(enrichment_scores_cluster_number), scores = enrichment_scores_cluster_number, ncells = sum(seurat_object@meta.data$seurat_clusters==cluster_number)), 10)
                                         }))
  
  # Get final ScType results - For each cluster, keep the top annotation
  sctype_scores = cluster_label_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # Set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  
  # Create annotation layer in meta.data
  seurat_object@meta.data[annotation_layer_name] = ""
  
  # Add annotation to meta data
  for(j in unique(sctype_scores$cluster)){
    # Extract current cluster
    cluster_label_type = sctype_scores[sctype_scores$cluster == j,]
    # Add annotations to meta data for current cluster (first cell type is kept in case of equal scores between types)
    seurat_object@meta.data[seurat_object@meta.data$seurat_clusters == j, annotation_layer_name] = as.character(cluster_label_type$type[1])
  }
  
  return(list(seurat_object = seurat_object,
              sctype_all_scores = cluster_label_results))
  
}

######################
## SINGLER ANALYSIS ##
######################

#' This function wraps the steps necessary for SingleR automatic cell annotation.
#' This will assign a label from the provided reference to each cell of in
#' the Seurat object. This is a correlation-based annotation method.
#'
#' @param seurat_object Seurat object containing cells to annotate.
#' @param singler_reference A numeric matrix of reference expression values, 
#' usually containing log-expression values. Alternatively, a SummarizedExperiment 
#' object containing such a matrix.
#' @param active_assay Name of the assay to set as default assay.
#' @param cluster_layer_name Name of the cluster layer store in meta.data and
#' containing the clusters to annotate.
#' @param annotation_layer_name Name of the meta.data layer that will be created
#' upon annotation with SingleR. It will contain the resulting annotations for
#' each cells in the seurat object.
#' 
#' @return A list with the following elements:
#'    - seurat_object: Updated Seurat object that was provided which now contains
#'    a layer in meta data with SingleR annotations for each cell.
#'    - singler_results: Detailed results provided by singleR (scores, ...).
#'    
singleRWrapper = function(seurat_object,
                          singler_reference,
                          active_assay,
                          cluster_layer_name,
                          annotation_layer_name)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.character(cluster_layer_name))) stop("Provided cluster_layer_name is not a character.")
  if (!(is.character(annotation_layer_name))) stop("Provided annotation_layer_name is not a character.")
  
  # Activate current assay (get back assay name and required initial clusters)
  seurat_object = activateSeuratAssay(seurat_object,
                                      ident = cluster_layer_name,
                                      active_assay = active_assay)
  
  # Convert Spatial to Single-Cell object
  sce = as.SingleCellExperiment(seurat_object)
  
  # Convert reference Seurat to Single-Cell object
  singler_reference = as.SingleCellExperiment(singler_reference)
  
  # SingleR automatic annotation
  singler_predictions = SingleR(test = sce, 
                                ref = singler_reference,
                                labels = singler_reference[[annotation_layer_name]])
  
  # Add the predictions to the Spatial Seurat object
  seurat_object@meta.data[[annotation_layer_name]] = singler_predictions@listData[["labels"]]
  
  return(list(seurat_object = seurat_object,
              singler_results = singler_predictions))
}
