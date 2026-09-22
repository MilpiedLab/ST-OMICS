
## @knitr random_neighborhood_comparison

cat("\n\n")
cat(sprintf("%s **COMPARISON TO RANDOM NEIGHBORHOODS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to compare the average cell-neighbors cell-type compositions 
of BANKSY neighborhoods to a random distribution of cell-neighbors cell-type 
compositions.

This random distribution is obtained by performing N random samplings of X cells, 
and computing the average cell-neighbors cell-type composition for each sampling. 
For each cell in a given sampling, its K nearest neighbors are identified, and their 
average cell-neighbors cell-type composition is computed as described previously.

The following sections are then provided to compare BANKSY neighborhood compositions 
with those of randomly sampled neighborhoods:

For each sample:

  - Cell type distribution in random neighborhoods: Boxplots showing the distribution 
of cell-type proportions obtained from randomly sampled cells within the current sample.

For each clustering resolution used in the BANKSY analysis:

  - A red dot is overlaid on each boxplot, representing the average cell-type 
proportion in the cell-neighborhoods of the corresponding BANKSY clusters. 
This allows direct comparison between observed neighborhood compositions and those 
expected under random sampling.")
cat("\n\n")

################################################################################
#'                        COMPARE TO RANDOM NEIGHBORHOODS
################################################################################



# Activate current assay (sample id layer)
banksy_xenium_data = activateSeuratAssay(banksy_xenium_data, 
                                         ident = "SAMPLE_ID", 
                                         active_assay = XENIUM_ASSAY_NAME)

# Iterate over available samples in merged data
for (sample in unique(banksy_xenium_data@meta.data[["SAMPLE_ID"]])){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), sample), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Subset current sample from global merged data
  sample_subset = subset(banksy_xenium_data,
                         subset = seurat_clusters == sample)
  
  #################################################
  ## 1. RETRIEVE RANDOM NEIGHBORHOOD COMPOSITION ##
  #################################################
  
  # Assess cell types distributions in neighborhood of randomly selected cells
  random_neighborhoods = get_cell_counts_in_random_cells_neighborhood(seurat_object = sample_subset,
                                                                      active_assay = XENIUM_ASSAY_NAME,
                                                                      k_nearest_neighbors = K_GEOM,
                                                                      number_of_sampling = NUMBER_RANDOM_SAMPLING,
                                                                      sampling_size = NUMBER_RANDOM_CELLS,
                                                                      neighbors_annotation = NEIGHBORS_ANNOTATION)
  
  # Retrieve cell type proportions in random neighborhoods
  boxplot_random_neighborhoods_data = as.data.frame(do.call(cbind, random_neighborhoods$proportions))
  
  # Format random neighborhoods composition in long format
  boxplot_random_neighborhoods_data = pivot_longer(boxplot_random_neighborhoods_data,
                                                   cols = colnames(boxplot_random_neighborhoods_data),
                                                   names_to = "cell_type",
                                                   values_to = "random_proportion")
  
  # Plot cell type distributions in random neighborhoods
  cat("\n\n")
  cat(sprintf("%s Cell types distribution in random neighborhoods", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(ggplot(boxplot_random_neighborhoods_data, aes(x = cell_type, y = random_proportion, fill = cell_type)) +
          geom_boxplot() +
          scale_fill_manual(values = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          labs(title = "Cell types distribution in random neighborhoods", x = "Cell type", y = "Proportion") +
          theme(legend.position="none"))
  cat("\n\n")
  
  ###########################################
  ## 2. RETRIEVE NEIGHBORHOODS COMPOSITION ##
  ###########################################
  
  # Iterate over resolutions
  for (res in RES_RANGE){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), res), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Define current neighborhood cluster layer
    temp_neighborhood_cluster_layer = sprintf("%s_snn_res.%s", 
                                              BANKSY_ASSAY_NAME, 
                                              res)
    
    # Activate current neighborhood layer in meta data
    sample_subset = activateSeuratAssay(sample_subset, 
                                        ident = temp_neighborhood_cluster_layer)
    
    # Compute cell type composition in neighborhoods of each domain of current meta data layer
    for (neighborhood in 0:(length(unique(sample_subset$seurat_clusters))-1)){
      
      cat("\n\n")
      cat(sprintf("%s %s", strrep("#", 4), neighborhood))
      cat("\n\n")
      
      # Retrieve cell type composition of cell neighborhood for current neighborhood assignment
      current_neighborhood_composition = get_cell_counts_in_cells_neighborhood(seurat_object = sample_subset,
                                                                               cluster_layer_name = temp_neighborhood_cluster_layer,
                                                                               active_assay = XENIUM_ASSAY_NAME,
                                                                               cluster_of_interest = as.character(neighborhood),
                                                                               k_nearest_neighbors = KNN_NEIGHBORHOOD_COMPOSITION,
                                                                               neighbors_annotation = NEIGHBORS_ANNOTATION)
      
      # Extract observed neighborhood cell type proportions for current condition
      observed_cell_type_proportions = data.frame(cell_type = names(current_neighborhood_composition$proportions),
                                                  proportions = as.numeric(current_neighborhood_composition$proportions))
      
      # Display comparison with random neighborhoods
      cat("\n\n")
      print(ggplot(boxplot_random_neighborhoods_data, aes(x = cell_type, y = random_proportion, fill = cell_type)) +
              geom_boxplot() +
              geom_point(data = observed_cell_type_proportions, aes(x = cell_type, 
                                                                    y = proportions), 
                         color = "red", size = 2) +
              scale_fill_manual(values = COMPLETE_CELL_TYPES_COLOR_PALETTE) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(title = sprintf("Cell types distribution in cell neighborhoods (res = %s, neighborhood = %s) compared to random neighborhoods",
                                   res,
                                   neighborhood), 
                   x = "Cell type", 
                   y = "Proportion") +
              theme(plot.title = element_text(size = 10),
                    legend.position="none"))
      cat("\n\n")
      
    }
  }
}
