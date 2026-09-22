################################################################################
##                             MARKER GENES FUNCTIONS                         ##
################################################################################

#######################
## FIND MARKER GENES ##
#######################

#' This function aims to compute marker genes corresponding to each cluster in
#' the provided cluster layer of a Seurat object.
#' 
#' Libraries: Seurat, immunogenomics/presto (optimized "FindAllMarkers")
#' 
#' @param seurat_object Seurat object on which marker genes will be assigned to
#' the active cluster layer.
#' @param active_assay Assay to set as default working assay in Seurat object.
#' @param cluster_layer Cluster layer to set as default layer and use for
#' finding marker genes.
#' @param logfc_threshold "FindAllMarkers": Limit testing to genes which show, 
#' on average, at least X-fold difference (log-scale) between the two groups 
#' of cells. 
#' @param min_pct "FindAllMarkers": only test genes that are detected in a 
#' minimum fraction of min.pct cells in either of the two populations.
#' @param only_pos "FindAllMarkers": Only return positive markers.
#' @param n_top_markers N top markers to retain for each cluster in the output 
#' results. This filtering is based on the average log fold 2 change. By 
#' default, all marker genes are outputted.
#' @param output_directory Output directory in which to save the resulting
#' marker genes (.xlsx file). By default, these results are not saved.
#' 
#' @return Top markers (all markers by default) of each cluster obtained after
#' the analysis. The output is a data.frame with the following columns:
#' p_val / avg_log2FC / pct.1 / pct.2 / p_val_adj / cluster / gene
#' 
find_all_markers = function(seurat_object, 
                            active_assay, 
                            cluster_layer,
                            logfc_threshold=0.25, 
                            min_pct=0.1, 
                            only_pos=TRUE,
                            n_top_markers=NULL)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.character(cluster_layer))) stop("Provided cluster_layer is not a character.")
  if (!(is.numeric(logfc_threshold))) stop("Provided logfc_threshold is not a numeric")
  if (!(is.numeric(min_pct))) stop("Provided min_pct is not a numeric")
  if (!(is.logical(only_pos))) stop("Provided only_pos is not logical")
  if (!(is.null(n_top_markers)) 
      && !(is.numeric(n_top_markers))) stop("Provided n_top_markers is not a numeric")
  
  # Activate corresponding assay and cluster layer
  seurat_object = activateSeuratAssay(seurat_object, active_assay = active_assay, ident = cluster_layer)
  
  # Get current cluster identities
  clusters = Idents(seurat_object)
  
  # Convert identities to a factor ordered numerically
  Idents(seurat_object) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
  
  # Compute all marker genes for each cluster
  marker_genes = FindAllMarkers(object = seurat_object, 
                                assay = active_assay, 
                                logfc.threshold = logfc_threshold,
                                min.pct = min_pct,
                                only.pos = only_pos)
  
  # Filter top N marker genes if necessary
  if (!(is.null(n_top_markers))){
    # Filter top N marker genes based on average log2 fold change
    top_markers = do.call(rbind, 
                          lapply(split(marker_genes, 
                                       marker_genes$cluster), 
                                 function(x) head(x[order(-x$avg_log2FC), ], n_top_markers)))
  } else {
    # Keep all marker genes to output
    top_markers = marker_genes
  }
  
  return(top_markers)
}

###############################
## SAVE MARKER GENES IN XLSX ##
###############################

#' This functions aims to save a list of data frames in a .xlsx file with a sheet
#' for each provided data frame. This function is used mainly to save output
#' from marker genes discovery.
#' 
#' Libraries: openxlsx
#' 
#' @param marker_genes List where keys are sheet names in the output .xlsx and
#' values are the data frames to write in the corresponding sheet.
#' @param output_directory Output directory in which to save the .xlsx file.
#' @param output_file_name Name of the output .xlsx file.
#' 
#' @return Returns NULL

save_marker_genes_xlsx = function(marker_genes, 
                                  output_directory, 
                                  output_file_name)
{
  # Verify inputs
  if (!(is.list(marker_genes) && all(sapply(marker_genes, is.data.frame)))) stop("Provided marker_genes is not a list of data frames.")
  if (!(is.character(output_directory)) 
      && !(dir.exists(output_directory))) stop("Provided output_directory is invalid/does not exists.")
  if (!(is.character(output_file_name))) stop("Provided output_file_name is not a character.")
  
  # Create output object
  wb = createWorkbook()
  
  # Create .xlsx sheets
  for (sheet in names(marker_genes)){
    # Create sheet
    addWorksheet(wb, 
                 sheet)
    # Add corresponding data to sheet
    writeData(wb, 
              sheet, 
              marker_genes[[sheet]])
  }
  
  # Save .xlsx
  saveWorkbook(wb, file = file.path(output_directory, output_file_name), overwrite = TRUE)
}

#################################
## PLOT MARKER GENES - VIOLINS ##
#################################

#' This function aims to plot the top N markers of a given cluster in a Seurat 
#' object as individual violin plots (gene expression distribution of marker 
#' genes from a given cluster across the clusters displayed as violin plots).
#' 
#' Libraries: Seurat
#' 
#' @param seurat_object Seurat object from which clusters and data are used to 
#' create corresponding plots
#' @param all_marker_genes Marker genes obtained from the provided cluster layer
#' after "FindAllMarkers" step. The dataframe is expected to have the following
#' format:
#' p_val / avg_log2FC / pct.1 / pct.2 / p_val_adj / cluster / gene
#' @param cluster_layer Name of the cluster layer in Seurat object meta data
#' from which marker genes were computed.
#' @param cluster_number Cluster of interest. The marker genes specific to this
#' cluster are the marker genes plotted by this function.
#' @param color_palette Color palette for plotting.
#' @param top_n_markers Plot only the top N marker genes based on the average 
#' log fold 2 change. Set as FALSE by default.
#' 
#' @return A list of violin plots. The keys in the list correspond to feature
#' (gene) names while the values are violin plot objects corresponding to the
#' feature.
get_marker_genes_violins = function(seurat_object, 
                                    all_marker_genes,
                                    cluster_layer_name,
                                    cluster_number, 
                                    top_n_markers=FALSE)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!is.data.frame(all_marker_genes)) stop("Input all_marker_genes must be a data frame.")
  if (nrow(all_marker_genes) == 0) stop("Data frame all_marker_genes is empty.")
  required_cols = c("avg_log2FC", "cluster", "gene")
  if (!all(required_cols %in% colnames(all_marker_genes))) stop(paste("Data frame must contain the following columns:", paste(required_cols, collapse = ", ")))
  if (is.null(seurat_object@meta.data[[cluster_layer_name]])) stop("Provided cluster_layer_name does not exist in object meta.data.")
  if (!(is.numeric(cluster_number))) stop("Provided cluster_number is not numeric")
  if (!(cluster_number %in% seurat_object@meta.data[[cluster_layer_name]])) stop("Provided cluster_number is not an existing cluster in the required cluster layer.")
  clusters_in_layer = unique(seurat_object@meta.data[[cluster_layer_name]])
  clusters_in_markers = unique(all_marker_genes$cluster)
  if (length(setdiff(clusters_in_layer, clusters_in_markers)) > 0) warning("Different number of clusters in cluster_layer and in all_marker_genes.")
  
  # Activate required cluster layer
  seurat_object = activateSeuratAssay(seurat_object,
                                      ident = cluster_layer_name)
  
  # Get current cluster identities
  clusters = Idents(seurat_object)
  
  # Convert identities to a factor ordered numerically
  Idents(seurat_object) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
  
  # Keep top N marker genes if required
  if (!(isFALSE(top_n_markers))){ # Filter top marker genes based on average log2 fold change
    top_markers = do.call(rbind, 
                          lapply(split(all_marker_genes, 
                                       all_marker_genes$cluster), 
                                 function(x) head(x[order(-x$avg_log2FC), ], top_n_markers)))
  } else { # No filtering
    top_markers = all_marker_genes
  }
  
  # Keep markers for current cluster
  cluster_top_markers = top_markers[top_markers$cluster == cluster_number,]
  
  # Cluster marker genes specific plots
  violins = list()
  for (feature in cluster_top_markers$gene){
    temp_violin_plot = VlnPlot(seurat_object, features = feature, pt.size = 0) + 
      theme(legend.position = 'none')
    violins[[feature]] = temp_violin_plot
  }
  
  return(violins)
}

#######################################
## PLOT MARKER GENES - FEATURE PLOTS ##
#######################################

#' This function aims to plot the top N markers of a given cluster in a Seurat 
#' object as individual feature plots (gene expression distribution of marker 
#' genes from a given cluster across the clusters displayed for each cell
#' in the UMAP dimensions).
#' 
#' Libraries: Seurat
#' 
#' @param seurat_object Seurat object from which clusters and data are used to 
#' create corresponding plots
#' @param all_marker_genes Marker genes obtained from the provided cluster layer
#' after "FindAllMarkers" step. The dataframe is expected to have the following
#' format:
#' p_val / avg_log2FC / pct.1 / pct.2 / p_val_adj / cluster / gene
#' @param cluster_layer Name of the cluster layer in Seurat object meta data
#' from which marker genes were computed.
#' @param cluster_number Cluster of interest. The marker genes specific to this
#' cluster are the marker genes plotted by this function.
#' @param color_palette Color palette for plotting.
#' @param top_n_markers Plot only the top N marker genes based on the average 
#' log fold 2 change. Set as FALSE by default.
#' 
#' @return A list of feature plots. The keys in the list correspond to feature
#' (gene) names while the values are feature plot objects corresponding to the
#' feature.
#' 
get_marker_genes_feature_plots = function(seurat_object,
                                          all_marker_genes,
                                          cluster_layer_name,
                                          cluster_number,
                                          top_n_markers=FALSE)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!is.data.frame(all_marker_genes)) stop("Input all_marker_genes must be a data frame.")
  required_cols <- c("avg_log2FC", "cluster", "gene")
  if (!all(required_cols %in% colnames(all_marker_genes))) stop(paste("Data frame must contain the following columns:", paste(required_cols, collapse = ", ")))
  if (is.null(seurat_object@meta.data[[cluster_layer_name]])) stop("Provided cluster_layer_name does not exist in object meta.data.")
  if (!(is.numeric(cluster_number))) stop("Provided cluster_number is not a numeric")
  if (!(cluster_number %in% seurat_object@meta.data[[cluster_layer_name]])) stop("Provided cluster_number is not an existing cluster in the required cluster layer.")
  if (nrow(all_marker_genes) == 0) stop("Data frame all_marker_genes is empty.")
  clusters_in_layer = unique(seurat_object@meta.data[[cluster_layer_name]])
  clusters_in_markers = unique(all_marker_genes$cluster)
  if (length(setdiff(clusters_in_layer, clusters_in_markers)) > 0) warning("Different number of clusters in cluster_layer and in all_marker_genes.")
  
  # Activate required cluster layer
  seurat_object = activateSeuratAssay(seurat_object,
                                      ident = cluster_layer_name)
  
  # Get current cluster identities
  clusters = Idents(seurat_object)
  
  # Convert identities to a factor ordered numerically
  Idents(seurat_object) = factor(clusters, levels = sort(unique(as.numeric(as.character(clusters)))))
  
  # Keep top N marker genes if required
  if (!(isFALSE(top_n_markers))){ # Filter top marker genes based on average log2 fold change
    top_markers = do.call(rbind, 
                          lapply(split(all_marker_genes, 
                                       all_marker_genes$cluster), 
                                 function(x) head(x[order(-x$avg_log2FC), ], top_n_markers)))
  } else { # No filtering
    top_markers = all_marker_genes
  }
  
  # Keep markers for current cluster
  cluster_top_markers = top_markers[top_markers$cluster == cluster_number,]
  
  # Cluster marker genes specific plots
  feature_plots = list()
  for (feature in cluster_top_markers$gene){
    temp_feature_plot = FeaturePlot(seurat_object, features = feature)
    feature_plots[[feature]] = temp_feature_plot
  }
  
  return(feature_plots)
}
