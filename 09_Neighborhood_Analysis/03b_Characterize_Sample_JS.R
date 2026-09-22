
## @knitr characterize_sample_JS

cat("\n\n")
cat(sprintf("%s **JENSEN–SHANNON DIVERGENCE**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Similarly, BANKSY neighborhoods within a sample can also be compared using 
the Jensen–Shannon divergence, a metric used to quantify the similarity between 
two probability distributions. In this context, the cell-type composition of each 
neighborhood is considered as its probability distribution.

To construct these probability distributions, cells belonging to a given BANKSY 
neighborhood are first subsetted. For each cell, its k-nearest neighbors are 
identified, and the cell-type composition within this local neighborhood is 
retrieved. The probability distribution for a given BANKSY neighborhood is then 
obtained by averaging these local cell-type compositions across all cells in that 
BANKSY neighborhood.

A smaller Jensen–Shannon divergence between two BANKSY neighborhoods indicates 
that their average cell-neighbors cell-type compositions are more similar.

Using both the sliced Wasserstein distance and the Jensen–Shannon divergence 
allows comparison of neighborhoods within a sample by quantifying their similarity 
in terms of spatial organization and cell-type composition of local microenvironments.

One potential application of these similarity measures is to assess whether the 
clustering resolution used to define neighborhoods is too high. In such cases, 
neighborhoods may be over-segmented, and highly similar clusters could correspond 
to the same underlying tissue region. These neighborhoods can then be considered 
for merging to obtain a more biologically meaningful representation of tissue 
structure.")
cat("\n\n")

################################################################################
#'                           JENSEN–SHANNON DIVERGENCE
################################################################################

# Store JS results to plot later
js_divergence_matrices = list()

# Activate current assay (sample id layer)
banksy_xenium_data = activateSeuratAssay(banksy_xenium_data, 
                                         ident = "SAMPLE_ID", 
                                         active_assay = XENIUM_ASSAY_NAME)

# Store neighborhood compositions for each sample
sample_neighborhood_compositions = list()

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
    
    ############################################
    ## 2. RETRIEVE NEIGHBORHOODS COMPOSITIONS ##
    ############################################
    
    # Compute cell type composition in neighborhoods of each domain of current meta data layer
    for (neighborhood in 0:(length(unique(sample_subset$seurat_clusters))-1)){
      
      # Retrieve cell type composition of cell neighborhood for current neighborhood assignment
      current_neighborhood_composition = get_cell_counts_in_cells_neighborhood(seurat_object = sample_subset,
                                                                               cluster_layer_name = temp_neighborhood_cluster_layer,
                                                                               active_assay = XENIUM_ASSAY_NAME,
                                                                               cluster_of_interest = as.character(neighborhood),
                                                                               k_nearest_neighbors = KNN_NEIGHBORHOOD_COMPOSITION,
                                                                               neighbors_annotation = NEIGHBORS_ANNOTATION)
      
      # Store cell type proportions in current neighborhood
      sample_neighborhood_compositions[[as.character(res)]][[sample]][[as.character(neighborhood)]] = current_neighborhood_composition$proportions
      
    }
    
    ###############################
    ## 3. COMPUTE JS DIVERGENCES ##
    ###############################
    
    # Create empty matrix to store JS divergence for current condition
    nb_neighborhoods = length(unique(sample_subset$seurat_clusters))
    js_divergences = matrix(, nb_neighborhoods, nb_neighborhoods)
    rownames(js_divergences) = sort(as.numeric(as.character(unique(sample_subset$seurat_clusters))))
    colnames(js_divergences) = sort(as.numeric(as.character(unique(sample_subset$seurat_clusters))))
    
    # Compute JS divergence for each pair of domain
    for (neighborhood_1 in 0:(length(unique(sample_subset$seurat_clusters))-1)){
      for (neighborhood_2 in neighborhood_1:(length(unique(sample_subset$seurat_clusters))-1)){
        
        # Retrieve cell type probabilities for current domains
        neighborhood_1_probabilities = as.numeric(sample_neighborhood_compositions[[as.character(res)]][[sample]][[as.character(neighborhood_1)]][1, ])
        neighborhood_2_probabilities = as.numeric(sample_neighborhood_compositions[[as.character(res)]][[sample]][[as.character(neighborhood_2)]][1, ])
        
        # Compute JS divergence
        merged_probabilities = rbind(neighborhood_1_probabilities,
                                     neighborhood_2_probabilities)
        JS_divergence = JSD(merged_probabilities)
        
        # Compute JS divergence
        js_divergences[as.character(neighborhood_1), as.character(neighborhood_2)] = JS_divergence
        js_divergences[as.character(neighborhood_2), as.character(neighborhood_1)] = JS_divergence
        
        # KL divergence alternative
        # kl_divergences[as.character(neighborhood_1), as.character(neighborhood_2)] = distance(rbind(neighborhood_1_probabilities, neighborhood_2_probabilities), method = "kullback-leibler") # D_KL(p || q)
        # kl_divergences[as.character(neighborhood_2), as.character(neighborhood_1)] = distance(rbind(neighborhood_2_probabilities, neighborhood_1_probabilities), method = "kullback-leibler") # D_KL(q || p)
        
      }
    }
    
    # Store divergence matrices to plot later
    js_divergence_matrices[[sample]][[as.character(res)]] = js_divergences
    
  }
}
