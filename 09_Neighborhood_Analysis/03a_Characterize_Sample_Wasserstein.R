
## @knitr characterize_sample_wasserstein

cat("\n\n")
cat(sprintf("%s **SLICED WASSERSTEIN DISTANCE**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Now that neighborhoods have been defined for each sample, they can be further 
characterized in greater detail.

This step aims to compare neighborhoods within each sample using the sliced 
Wasserstein distance. The Wasserstein distance can be interpreted as the amount 
of “mass” that must be transported to transform one neighborhood’s spatial 
distribution into another. In this context, a smaller distance indicates that 
two neighborhoods are more similar in their spatial organization within the 
tissue.

Here, the sliced Wasserstein distance is used, which is an approximation of the 
full Wasserstein distance that enables faster computation while preserving 
meaningful comparisons between distributions.")
cat("\n\n")

################################################################################
#'                           SLICED WASSERSTEIN DISTANCE
################################################################################

# Store wsd results to plot later
wasserstein_matrices = list()

# Activate current assay (sample id layer)
banksy_xenium_data = activateSeuratAssay(banksy_xenium_data, 
                                         ident = "SAMPLE_ID", 
                                         active_assay = XENIUM_ASSAY_NAME)

# Iterate over available samples in merged data
for (sample in unique(banksy_xenium_data@meta.data[["SAMPLE_ID"]])){
  
  ######################
  ## 1. SUBSET SAMPLE ##
  ######################
  
  # Subset current sample from global merged data
  sample_subset = subset(banksy_xenium_data,
                         subset = seurat_clusters == sample)
  
  # Iterate over resolutions
  for (res in RES_RANGE){
    
    # Define current neighborhood cluster layer
    temp_neighborhood_cluster_layer = sprintf("%s_snn_res.%s", 
                                              BANKSY_ASSAY_NAME, 
                                              res)
    
    # Activate current neighborhood layer in meta data
    sample_subset = activateSeuratAssay(sample_subset, 
                                        ident = temp_neighborhood_cluster_layer)
    
    #####################################
    ## 2. COMPUTE WASSERSTEIN DISTANCE ##
    #####################################
    
    # Retrieve required data for Wasserstein distance computation (cell position and neighborhood assignment)
    wasserstein_data = data.frame(x = sample_subset@meta.data[["x_centroid"]],
                                  y = sample_subset@meta.data[["y_centroid"]],
                                  neighborhood = as.character(sample_subset@meta.data[[temp_neighborhood_cluster_layer]]))
    
    # Init empty matrix to store Wasserstein distances for current condition
    nb_neighborhoods = length(unique(sample_subset$seurat_clusters))
    wasserstein_distances = matrix(, nb_neighborhoods, nb_neighborhoods)
    rownames(wasserstein_distances) = sort(as.numeric(as.character(unique(sample_subset$seurat_clusters))))
    colnames(wasserstein_distances) = sort(as.numeric(as.character(unique(sample_subset$seurat_clusters))))
    
    # Compute sliced Wasserstein distance for each pair of neighborhoods
    for (neighborhood_1 in 0:(length(unique(sample_subset$seurat_clusters))-1)){
      for (neighborhood_2 in neighborhood_1:(length(unique(sample_subset$seurat_clusters))-1)){
        
        # Retrieve cell coordinates for current neighborhoods
        coords_neighborhood_1 = wasserstein_data[wasserstein_data$neighborhood == neighborhood_1, c("x", "y")]
        coords_neighborhood_2 = wasserstein_data[wasserstein_data$neighborhood == neighborhood_2, c("x", "y")]
        
        ## Min-Max Scaling of the coordinates ##
        
        # Find min and max value in MERGED systems
        min_x = min(c(coords_neighborhood_1$x, coords_neighborhood_2$x))
        max_x = max(c(coords_neighborhood_1$x, coords_neighborhood_2$x))
        min_y = min(c(coords_neighborhood_1$y, coords_neighborhood_2$y))
        max_y = max(c(coords_neighborhood_1$y, coords_neighborhood_2$y))
        
        # Normalize coordinates in both domains
        coords_neighborhood_1$x = (coords_neighborhood_1$x - min_x) / (max_x - min_x)
        coords_neighborhood_2$x = (coords_neighborhood_2$x - min_x) / (max_x - min_x)
        coords_neighborhood_1$y = (coords_neighborhood_1$y - min_y) / (max_y - min_y)
        coords_neighborhood_2$y = (coords_neighborhood_2$y - min_y) / (max_y - min_y)
        
        ## Compute sliced Wasserstein ##
        
        # Convert cell coordinates to numpy arrays
        system1 = np$array(as.matrix(coords_neighborhood_1))
        system2 = np$array(as.matrix(coords_neighborhood_2))
        
        # Compute sliced Wasserstein distance between current neighborhoods
        temp_wasserstein_distance = ot$sliced_wasserstein_distance(system1, 
                                                                   system2, 
                                                                   n_projections = 1000L, 
                                                                   p = 1L)
        
        # Update wasserstein distances matrix
        wasserstein_distances[as.character(neighborhood_1), as.character(neighborhood_2)] = temp_wasserstein_distance
        wasserstein_distances[as.character(neighborhood_2), as.character(neighborhood_1)] = temp_wasserstein_distance
        
      }
    }
    
    # Store distance matrices to plot later
    wasserstein_matrices[[sample]][[as.character(res)]] = wasserstein_distances
    
  }
}
