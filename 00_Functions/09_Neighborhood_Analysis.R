################################################################################
##                    NEIGHBORHOOD ANALYSIS RELATED FUNCTIONS                 ##
################################################################################

#####################
## BANKSY ANALYSIS ##
#####################

#' This function wraps the steps necessary for BANKSY domain discovery analysis.
#' First, new embeddings are computed for cells in the data set based on 
#' neighboring cells (RunBanksy step). Then, classic clustering steps are made
#' in order to retrieve what will be considered as the domains (a.k.a. niches) 
#' of our sample.
#' 
#' Libraries: Banksy
#' 
#' @param seurat_object Seurat object from which domains will be inferred.
#' @param lambda Spatial weight parameter: Banksy parameter that defines how 
#' much neighboring cells will influence a given cell embeddings. If lambda = 1, 
#' there are no contributions of the cell of interest for its own embeddings.
#' @param active_assay Name of the assay to set as default assay.
#' @param active_slot Slot name within the "active_assay" to use for the
#' analysis (e.g., "counts" or "data").
#' @param k_geom kNN parameter - number of neighbors to use
#' @param banksy_assay_name Banksy analysis creates a new assay in the seurat 
#' object. Assay_name allows you to define what will be the name of this assay.
#' @param dimension_principal_component_analysis Number of computed PCs during
#' the "RunPCA" step.
#' @param dimension_clustering Number of PCs kept for the subsequent clustering
#' step.
#' @param umap_layer_name Name of the UMAP coordinates layer that will be stored
#' in meta.data.
#' @param resolution_range Range of resolutions (can be a single value) to
#' compute in the "findClusters" step of the analysis.
#' 
#' @return Seurat object with new clusters layer(s) corresponding to domains
#' found in the data. UMAP coordinates resulting from BANKSY embeddings are also
#' saved in the meta.data of the seurat object.
#' 
banksy_neighborhood_analysis = function(seurat_object,
                                        lambda,
                                        active_assay,
                                        active_slot,
                                        k_geom,
                                        banksy_assay_name,
                                        dimension_principal_component_analysis,
                                        dimension_clustering,
                                        umap_layer_name,
                                        resolution_range)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.numeric(lambda))) stop("Provided lambda is not a numeric")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.character(active_slot))) stop("Provided active_slot is not a character.")
  if (!(is.numeric(k_geom))) stop("Provided k_geom is not a numeric")
  if (!(is.character(banksy_assay_name))) stop("Provided banksy_assay_name is not a character.")
  if (!(is.numeric(dimension_principal_component_analysis))) stop("Provided dimension_principal_component_analysis is not a numeric")
  if (!(is.numeric(dimension_clustering))) stop("Provided dimension_clustering is not a numeric")
  if (!(is.character(umap_layer_name))) stop("Provided umap_layer_name is not a character.")
  if (!(is.numeric(resolution_range))) stop("Provided resolution_range is not a numeric")
  if (!exists("storeUMAPembeddings", mode = "function", inherits = TRUE)) stop("The function 'storeUMAPembeddings' does not exist or is not a function.")
  
  # Activate assay
  seurat_object = activateSeuratAssay(seurat_object, 
                                      active_assay = active_assay)
  
  # Find variable features in the data
  seurat_object = FindVariableFeatures(seurat_object)
  
  # Banksy analysis
  seurat_object = RunBanksy(seurat_object, 
                            lambda = lambda, 
                            verbose=TRUE, 
                            assay = active_assay, 
                            slot = active_slot, 
                            features = 'variable', 
                            k_geom = k_geom,
                            assay_name = banksy_assay_name)
  
  # PCA from Banksy embeddings
  seurat_object = RunPCA(seurat_object, 
                          assay = banksy_assay_name, 
                          features = rownames(seurat_object), 
                          npcs = dimension_principal_component_analysis)
  
  # UMAP
  seurat_object = RunUMAP(seurat_object, 
                           dims = 1:dimension_clustering)
  # Save UMAP embeddings in meta.data
  seurat_object = storeUMAPembeddings(seurat_object = seurat_object, 
                                      umap_layer_name = umap_layer_name)
  
  # Clustering from Banksy embeddings
  seurat_object = FindNeighbors(seurat_object, 
                                 dims = 1:dimension_clustering)
  seurat_object = FindClusters(seurat_object, 
                                resolution = resolution_range)
  
  return(seurat_object)
}

########################
## CELL NEIGHBORHOODS ##
########################

#' This function aims to retrieve the counts of k nearest neighboring cells to 
#' a given cell class in a seurat object.
#' 
#' @param seurat_object Seurat object from which clusters and other cells data
#' will be retrieved.
#' @param active_assay Name of the assay to set as default assay.
#' @param cluster_layer_name Seurat object's meta data layer to which the
#' cell classes currently under study are found.
#' @param cluster_of_interest Cell class of interest. The cell neighborhood of
#' these cells will be retrieved.
#' @param k_nearest_neighbors Number of k nearest neighbors considered in the
#' cells neighborhoods.
#' @param neighbors_annotation The layer in meta data that must be used to
#' annotate the neighbors that will be found during the process. By default, it
#' is set to NULL and the neighbors will be annotated based on the
#' "cluster_layer_name" meta data layer.
#' 
#' @return A list with two entries:
#'        - "counts": Total counts per cell classes in cells neighborhood of
#'        cluster of interest.
#'        - "proportions": Proportions of cells classes in cells neighborhood
#'        of cluster of interest.
#'        
get_cell_counts_in_cells_neighborhood = function(seurat_object,
                                                 cluster_layer_name,
                                                 active_assay,
                                                 cluster_of_interest,
                                                 k_nearest_neighbors,
                                                 neighbors_annotation = NULL)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.character(cluster_layer_name))) stop("Provided cluster_layer_name is not a character.")
  if (!(is.character(cluster_of_interest))) stop("Provided cluster_of_interest is not a character.")
  if (!(is.numeric(k_nearest_neighbors))) stop("Provided k_nearest_neighbors is not a numeric")
  
  #####################
  ## 1. PREPARE DATA ##
  #####################
  
  # Define meta data layer to use to annotate neighbors
  if (is.null(neighbors_annotation)){
    neighbors_annotation = cluster_layer_name
  }
  
  # Activate current assay (get back assay name and required initial clusters)
  seurat_object = activateSeuratAssay(seurat_object, 
                                      ident = cluster_layer_name,
                                      active_assay = active_assay)
  
  # Retrieve all data from seurat object
  complete_data = data.frame(cell_name = rownames(seurat_object@meta.data),
                             cluster_layer = seurat_object@meta.data[[cluster_layer_name]],
                             x = seurat_object@meta.data$x_centroid,
                             y = seurat_object@meta.data$y_centroid)
  
  # Retrieve cell coordinates from complete seurat object
  cell_coordinates = data.frame(x = seurat_object@meta.data[["x_centroid"]],
                                y = seurat_object@meta.data[["y_centroid"]])
  
  # Get cell coordinates of cluster of interest
  cells_of_interest = complete_data[complete_data$cluster_layer == cluster_of_interest, ]
  cells_of_interest_coordinates = data.frame(x = cells_of_interest$x,
                                             y = cells_of_interest$y)
  
  # Retrieve all possible cell classes (e.g., clusters or cell types) in data set
  all_cell_classes = unique(seurat_object@meta.data[[neighbors_annotation]])
  
  ################################
  ## 2. RETRIEVE CELL NEIGHBORS ##
  ################################
  
  # Retrieve K+1 nearest neighbors of cells of interest
  cells_of_interest_neighbors = nn2(data = cell_coordinates,
                                    query = cells_of_interest_coordinates, 
                                    k = k_nearest_neighbors + 1)
  cells_of_interest_neighbors = cells_of_interest_neighbors[["nn.idx"]]
  cells_of_interest_neighbors = as.data.frame(cells_of_interest_neighbors)
  
  # Remove first column as it corresponds to the cell itself (nearest cell)
  cells_of_interest_neighbors$V1 = NULL
  
  # Convert cell indices to equivalent in "neighbors_annotation" meta data layer
  cells_of_interest_neighbors = as.data.frame(apply(cells_of_interest_neighbors, 2, function(col) seurat_object@meta.data[[neighbors_annotation]][as.numeric(col)]))
  
  #############################
  ## 3. RETRIEVE CELL COUNTS ##
  #############################
  
  # Get total cell counts in neighborhoods
  cells_of_interest_neighbors_counts = as.list(table(unlist(cells_of_interest_neighbors)))
  
  # Store counts in temp results
  counts = setNames(rep(0, length(all_cell_classes)), all_cell_classes)
  for (class in names(cells_of_interest_neighbors_counts)) {
    if (class %in% names(counts)) {
      counts[class] = cells_of_interest_neighbors_counts[[class]]
    }
  }
  
  # Convert temp results to a data frame
  counts = as.data.frame(t(counts[ all_cell_classes ]))
  
  # Convert counts to proportions of cell classes in cell neighborhoods
  proportions = counts / (k_nearest_neighbors * nrow(cells_of_interest))
  
  return(list("counts" = counts,
              "proportions" = proportions))
}

###############################
## RANDOM CELL NEIGHBORHOODS ##
###############################

#' This function aims to retrieve the counts of k nearest neighboring cells to 
#' a given amount a randomly selected cells within the seurat object over a
#' specified amount of samplings.
#' 
#' @param seurat_object Seurat object from which clusters and other cells data
#' will be retrieved.
#' @param active_assay Name of the assay to set as default assay.
#' @param cluster_layer_name Seurat object's meta data layer to which the
#' cell classes currently under study are found.
#' @param k_nearest_neighbors Number of k nearest neighbors considered in the
#' cells neighborhoods.
#' @param number_of_sampling Number of random cell selections and neighborhood
#' assessments.
#' @param sampling_size Number of cells randomly selected by random sampling.
#' @param neighbors_annotation The layer in meta data that must be used to
#' annotate the neighbors that will be found during the process.
#' 
#' @return A list with two entries:
#'        - "counts": Total counts per cell classes in random cells 
#'        neighborhoods.
#'        - "proportions": Proportions of cells classes in random cells 
#'        neighborhoods.
#'        
get_cell_counts_in_random_cells_neighborhood = function(seurat_object,
                                                        active_assay,
                                                        k_nearest_neighbors,
                                                        number_of_sampling,
                                                        sampling_size,
                                                        neighbors_annotation,
                                                        verbose = FALSE)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
  if (!(is.numeric(k_nearest_neighbors))) stop("Provided k_nearest_neighbors is not a numeric")
  if (!(is.numeric(number_of_sampling))) stop("Provided number_of_sampling is not a numeric")
  if (!(is.numeric(sampling_size))) stop("Provided sampling_size is not a numeric")
  
  #####################
  ## 1. PREPARE DATA ##
  #####################
  
  # Activate current assay (get back assay name and required initial clusters)
  seurat_object = activateSeuratAssay(seurat_object, 
                                      ident = neighbors_annotation,
                                      active_assay = active_assay)
  
  # Retrieve all data from seurat object
  complete_data = data.frame(cell_name = rownames(seurat_object@meta.data),
                             x = seurat_object@meta.data$x_centroid,
                             y = seurat_object@meta.data$y_centroid)
  
  # Retrieve cell coordinates from complete seurat object
  cell_coordinates = data.frame(x = seurat_object@meta.data[["x_centroid"]],
                                y = seurat_object@meta.data[["y_centroid"]])
  
  # Get all possible cell types from metadata
  all_cell_classes = unique(seurat_object@meta.data[[neighbors_annotation]])
  
  ## Initialize an empty count data frame with those cell types as columns
  # Counts
  random_counts = data.frame(matrix(0, nrow = 0, ncol = length(all_cell_classes)))
  colnames(random_counts) = all_cell_classes
  # Proportions
  random_proportions = data.frame(matrix(0, nrow = 0, ncol = length(all_cell_classes)))
  colnames(random_proportions) = all_cell_classes
  
  #######################################
  ## 2. RETRIEVE RANDOM CELL NEIGHBORS ##
  #######################################
  
  # Process random sampling
  for (i in 1:number_of_sampling){
    
    # Print ongoing process
    if (verbose){
      cat(sprintf("Processing sampling %s ... \r", i))
    }
    
    # Sample random cells coordinates
    random_cells_coordinates = cell_coordinates[sample(nrow(cell_coordinates), sampling_size), ]
    
    # Retrieve K+1 nearest neighbors of current random cells
    random_cells_neighbors = nn2(data = cell_coordinates,
                                 query = random_cells_coordinates, 
                                 k = k_nearest_neighbors + 1)
    random_cells_neighbors = random_cells_neighbors[["nn.idx"]]
    random_cells_neighbors = as.data.frame(random_cells_neighbors)
    
    # Remove first column as it corresponds to the cell itself (nearest cell)
    random_cells_neighbors$V1 = NULL
    
    # Convert cell indices to equivalent in "neighbors_annotation" meta data layer
    random_cells_neighbors = as.data.frame(apply(random_cells_neighbors, 2, function(col) seurat_object@meta.data[[neighbors_annotation]][as.numeric(col)]))
    
    #############################
    ## 3. RETRIEVE CELL COUNTS ##
    #############################
    
    # Get total cell counts in neighborhoods
    random_cells_neighbors_counts = as.list(table(unlist(random_cells_neighbors)))
    
    # Store counts in temp results
    temp_counts = setNames(rep(0, length(all_cell_classes)), all_cell_classes)
    for (class in names(random_cells_neighbors_counts)) {
      if (class %in% names(temp_counts)) {
        temp_counts[class] = random_cells_neighbors_counts[[class]]
      }
    }
    
    # Convert temp results to a data frame
    temp_counts = as.data.frame(t(temp_counts[ all_cell_classes ]))
    
    # Convert counts to proportions of cell classes in cell neighborhoods
    temp_proportions = temp_counts / (k_nearest_neighbors * sampling_size)
    
    # Append to results data frames
    random_counts = bind_rows(random_counts, temp_counts)
    random_proportions = bind_rows(random_proportions, temp_proportions)
    
  }
  
  # Return distribution of cell classes counts and proportions in random cells
  return(list("counts" = random_counts,
              "proportions" = random_proportions))
  
}
