################################################################################
##                           SPLIT TMA RELATED FUNCTIONS                      ##
################################################################################

############################################
## SPLIT SAMPLES FROM GLOBAL XENIUM PLATE ##
############################################

#' This function aims to extract individual samples from a seurat object that
#' contains multiple samples. Here, a TMA is considered as a cluster of cell 
#' that are all interconnected (have a connection with at least one other cell 
#' of the cluster) and have no connections with other clusters. The connection
#' between cells are made based on their spatial coordinates.
#' 
#' Libraries: Seurat, igraph, RANN
#' 
#' @param seurat_object Seurat object from which samples will be extracted.
#' @param radius_threshold Maximum radius distance to which a cell can be
#' considered as a potential neighbor.
#' @param neighbors_number Maximum number of neighbors that are computed for
#' each cell.
#' @param sample_minimum_size Minimum number of cells the sample needs to have
#' not to be removed.
#' @param number_sample_rows Number of rows (on the TMA) on which the samples
#' were dispatched.
#' @param number_sample_columns Number of columns (on the TMA) on which the 
#' samples were dispatched.
#' 
#' @return A list containing the following elements:
#'        - "grid": A data frame that reproduce the plate structure with the
#'        same number of rows and columns, and samples assigned to the correct
#'        place.
#'        - "splitted_seurat": The initial seurat object with the newly assigned
#'        clusters ident.
#'
split_xenium_tma_samples = function(seurat_object,
                                    radius_threshold = 100,
                                    neighbors_number = 50,
                                    sample_minimum_size = 100,
                                    number_sample_rows,
                                    number_sample_columns)
{
  
  ######################
  ## BUILD CELL GRAPH ##
  ######################
  
  # Retrieve cell coordinates
  cell_coordinates = data.frame(x = seurat_object$x_centroid, 
                                y = seurat_object$y_centroid)
  
  # Radius neighbor search for each cell in the data set
  cell_neighbors = nn2(cell_coordinates, 
                       cell_coordinates, 
                       searchtype = "radius", 
                       radius = radius_threshold, 
                       k = neighbors_number)
  
  # Build edge list from neighbors
  graph_edges = do.call(rbind, lapply(1:nrow(cell_neighbors$nn.idx), function(i){
    ni = cell_neighbors$nn.idx[i, ]
    ni = ni[ni != 0] # & ni != i] Remove 0s (and self-links)
    if (length(ni) > 0){
      cbind(rep(i, length(ni)), ni)
    }
  }))
  
  # Build graph from edge list
  graph = graph_from_edgelist(graph_edges, directed = FALSE)
  
  #########################################
  ## RETRIEVE SAMPLES AND FILL META DATA ##
  #########################################
  
  ## TODO: ici dans le cas ou un cell n'a pas de voisins --> pose problème car
  ## alors il compl-ète avec 1 comme valeur --> faut trouver solution pour ça --> peut être en empêchant pas les selfs links justement (pcq au final ça n'a pas bcp d'importance de les enlever sur le résultat final il me semble)
  
  # Retrieve connected components
  # (graphs where each cell is connected to at least another cell from the corresponding graph)
  connected_components = components(graph)
  
  # Add the graph identity to meta data
  # (this identity corresponds to the sample the cell is corresponding to)
  seurat_object$sample_identity = connected_components$membership
  
  # Set the sample identity as the default ident
  Idents(seurat_object) = "sample_identity"
  
  ####################
  ## FILTER SAMPLES ##
  ####################
  
  # Retrieve cell counts by sample identity
  sample_sizes = table(seurat_object$sample_identity)
  
  # Label small samples as "excluded"
  filtered_samples = ifelse(sample_sizes[as.character(seurat_object$sample_identity)] < sample_minimum_size, "excluded", seurat_object$sample_identity)
  
  # Add final samples to meta data
  seurat_object$filtered_sample_identity = filtered_samples
  
  # Set the sample identity as the default ident
  Idents(seurat_object) = "filtered_sample_identity"
  
  ##########################
  ## RETRIEVE SAMPLE GRID ##
  ##########################
  
  # Remove samples to be excluded to retrieve the grid
  filtered_seurat_object = subset(seurat_object, subset = filtered_sample_identity != "excluded")
  
  # Get coordinates of remaining cells in the object
  filtered_coordinates = data.frame(x = filtered_seurat_object$x_centroid,
                                    y = filtered_seurat_object$y_centroid,
                                    sample_identity = filtered_seurat_object$filtered_sample_identity)
  
  # Add sample identity to cell coordinates
  sample_centroids = filtered_coordinates %>%
    group_by(sample_identity) %>%
    summarize(centroid_x = mean(x),
              centroid_y = mean(y),
              .groups = "drop")
  
  # Get range of x and y values in CELL centroids 
  # (using sample centroids is a bit too short and might assign uncorrectly some samples to their rows/columns)
  x_range = range(filtered_seurat_object$x_centroid)
  y_range = range(filtered_seurat_object$y_centroid)
  
  # Partition samples to their corresponding row and column
  # (In seurat, x corresponds to row while y corresponds to columns)
  sample_centroids = sample_centroids %>%
    mutate(
      # X is vertical (top to bottom)
      row_label = cut(centroid_x, 
                      breaks = seq(x_range[1], x_range[2], length.out = number_sample_rows + 1),
                      labels = paste0("row", number_sample_rows:1),
                      include.lowest = TRUE),
      # Y is horizontal (left to right)
      col_label = cut(centroid_y, 
                      breaks = seq(y_range[1], y_range[2], length.out = number_sample_columns + 1),
                      labels = paste0("col", 1:number_sample_columns),
                      include.lowest = TRUE))
  
  # Convert to factor to ensure correct order in the grid data frame
  sample_centroids = sample_centroids %>%
    mutate(row_label = factor(row_label, levels = paste0("row", 1:number_sample_rows)),
           col_label = factor(col_label, levels = paste0("col", 1:number_sample_columns)))
  
  # Build grid of sample identities
  sample_grid = sample_centroids %>%
    group_by(row_label, col_label) %>%
    summarise(
      sample_ids = paste(sample_identity, collapse = ", "),
      .groups = "drop"
    ) %>%
    complete(row_label, col_label) %>%
    pivot_wider(
      names_from = col_label,
      values_from = sample_ids
    )
  
  # Clean sample grid
  sample_grid = as.matrix(sample_grid[ , -1])
  colnames(sample_grid) = NULL
  rownames(sample_grid) = NULL 
  
  ####################
  ## OUTPUT RESULTS ##
  ####################
  
  # Add sample grid and complete image to output
  results = list()
  results[["grid"]] = sample_grid
  Idents(seurat_object) = "filtered_sample_identity"
  results[["splitted_seurat"]] = seurat_object
  
  return(results)
}
