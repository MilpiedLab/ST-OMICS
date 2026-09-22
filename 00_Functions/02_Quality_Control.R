################################################################################
##                            QUALITY CONTROL FUNCTIONS                       ##
################################################################################

#################################
## PASTIS QUALITY CONTROL STEP ##
#################################

#' This function aims to make the quality control of cells based number of
#' features and cell area. This is the quality control step that was implemented
#' in the PASTIS pipeline.
#'
#' @param seurat_object Seurat object on which QC step will be applied.
#' @param minimum_features Minimum number of features a cell should have to be
#' kept.
#' @param maximum_area Maximum area a cell should have to be kept.
#' 
#' @return Seurat object with a new meta.data layer corresponding to tagged 
#' cells and whether they fit or not the filtering parameters that were defined.
#'
#' Note: This function does not filter the data. It just tags the cells that
#' should be filtered based on a predefined criteria.
#' 
pastis_quality_control = function(seurat_object,
                                  minimum_features=10,
                                  maximum_area=500)
{
  # Tag cells that should be removed based on filtering parameters
  seurat_object$lowQCell = TRUE
  seurat_object@meta.data[seurat_object$nFeature_Xenium >= minimum_features & seurat_object$cell_area < maximum_area, "lowQCell"] = FALSE
  
  # Split in factors
  seurat_object@meta.data$QC_result = ifelse(seurat_object$cell_area > maximum_area & seurat_object$nFeature_Xenium < minimum_features, "low_features_high_area",
                                             ifelse(seurat_object$cell_area > maximum_area, "high_area",
                                                    ifelse(seurat_object$nFeature_Xenium < minimum_features, "low_feature", "conserved")))
  seurat_object$QC_result = factor(seurat_object$QC_result, 
                                   levels = c("conserved", 
                                              "low_feature", 
                                              "high_area", 
                                              "low_features_high_area"))
  
  return(seurat_object)
}

##################################
## RETRIEVE MECR PAIRS FROM CSV ##
##################################

#' This function aims to retrieve all possible mutually exclusive gene pairs 
#' from a provided list of exclusive gene markers per cell type.
#' 
#' @param mutally_exclusive_genes A data frame with classes (e.g., cell types) 
#' as row names. The data frame should contain the following columns:
#'        - "Gene_list": For each row (cell type), a string containing the
#'        marker genes exclusive to the corresponding cell type and all
#'        separated by a comma (",").
#' @param gene_panel A character vector containing the gene panel (i.e., the
#' genes available in the data set of interest).
#' 
#' @return A data frame with two columns, each containing exclusive marker
#' genes from provided classes. Each row of the data frame corresponds to a
#' mutually exclusive marker gene pair available in the data set of interest.
#' 
retrieve_MECR_pairs = function(mutally_exclusive_genes, 
                               gene_panel)
{
  # Verify inputs
  mutally_exclusive_genes_required_cols = c("Gene_list")
  if (!all(mutally_exclusive_genes_required_cols %in% colnames(mutally_exclusive_genes))) stop(paste("Data frame must contain the following columns:", paste(mutally_exclusive_genes_required_cols, collapse = ", ")))
  if (!(is.character(gene_panel))) stop("Provided gene_panel is not a vector of characters.")
  
  # Store list of mutually exclusive gene pairs in the data set
  mutally_exclusive_genes_panel = list()
  
  ## Find the intersection between the known MECR and your gene panel
  # Iterate over classes (e.g., cell types) of known mutually exclusive genes
  for (class in rownames(mutally_exclusive_genes)){
    
    # Retrieve exclusive marker genes for current class
    gene_list = strsplit(mutally_exclusive_genes[class,], split = ",")[[1]]
    
    # Store exclusive marker genes that are available in the data set
    if (length(intersect(gene_list, gene_panel)) >= 1){
      mutally_exclusive_genes_panel[[class]] = intersect(gene_list, gene_panel)
    }
    
  }
  
  # Get all possible pairs of classes
  pairs_of_classes = combn(names(mutally_exclusive_genes_panel), 2)
  
  # Initialize pairs of exclusive marker genes
  exclusive_pairs = NULL
  pairs_names = NULL
  
  # Iterate over all classes pairs
  for (i in 1:dim(pairs_of_classes)[2]){
    
    # Retrieve current classes of interest
    class1 = pairs_of_classes[1,i]
    class2 = pairs_of_classes[2,i]
    
    # Define classes pair name
    pairs_names = c(pairs_names, paste0(class1,"_",class2))
    
    # Add all possible pairs of exclusive gene markers for current classes pair to already computed pairs
    exclusive_pairs = rbind(exclusive_pairs, 
                            tidyr::crossing(mutally_exclusive_genes_panel[[pairs_of_classes[1,i]]],
                                            mutally_exclusive_genes_panel[[pairs_of_classes[2,i]]]))
  }
  
  return(exclusive_pairs)
}

#####################################
## COMPUTE MECR SCORE PARALLELIZED ##
#####################################

#' This function aims to compute the MECR score of a provided Seurat object with
#' implemented parallelization.
#' 
#' @param seurat_object Seurat object from which MECR scores will be derived.
#' @param mutually_exclusive_gene_pairs A data frame with two columns, each 
#' containing exclusive marker genes from provided classes. Each row of the data 
#' frame corresponds to a mutually exclusive marker gene pair available in the 
#' data set of interest.
#' @param sample_name Sample ID that will be added in the result data frame.
#' @param batch_size Number of gene pairs per block being processed in each
#' core.
#' @param workers Number of cores used for MECR computation.
#' 
#' @return A data frame where rownames are the gene pairs and with the following 
#' columns:
#'        - "sample_ID": The user-provided sample name/ID.
#'        - "Gene_pair": The gene pair corresponding to the row and for which
#'        the MECR was computed.
#'        - "MECR": Mutually exclusive co expression ratio for the corresponding
#'        gene pair.
#'        - "Expression_rate": The proportion of cells expressing the 
#'        corresponding pair of mutual exclusive genes.
#'        
calculate_MECR_score_parallelized = function(seurat_object, 
                                             mutually_exclusive_gene_pairs, 
                                             sample_name,
                                             batch_size,
                                             workers)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  
  # Extract expression matrix once (smaller than Seurat object)
  expression_matrix = GetAssayData(seurat_object)
  all_genes = rownames(expression_matrix)
  
  # Parallelize by block of genes of size batch_size
  plan("multisession", workers = workers)
  
  # Total number of gene pairs
  total_pairs = nrow(mutually_exclusive_gene_pairs)
  
  # Split gene pairs into batches
  batches = split(1:total_pairs, ceiling(seq_along(1:total_pairs) / batch_size))
  
  # Function to process one batch
  process_batch = function(pair_indices) {
    results = lapply(pair_indices, function(pair_idx) {
      gene1 = as.character(mutually_exclusive_gene_pairs[pair_idx, 1])
      gene2 = as.character(mutually_exclusive_gene_pairs[pair_idx, 2])
      gene_pair = paste0(gene1, "_", gene2)
      
      # Skip if genes not in matrix
      if (!(gene1 %in% all_genes) || !(gene2 %in% all_genes)) {
        return(NULL)
      }
      
      # Get expression data for gene pair
      gene_expressions = expression_matrix[c(gene1, gene2), , drop = FALSE]
      gene_expressions = as.data.frame(t(as.matrix(gene_expressions)))
      
      expressed_counts = rowSums(gene_expressions > 0)
      MECR = sum(expressed_counts == 2) / sum(expressed_counts >= 1)
      expression_rate = sum(expressed_counts > 0) / nrow(gene_expressions)
      
      data.frame(
        sample_ID = sample_name,
        Gene_pair = gene_pair,
        MECR = as.numeric(MECR),
        Expression_rate = as.numeric(expression_rate),
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, results)
  }
  
  # Parallel map over batches
  MECR_results_list = future_map(batches, process_batch, .progress = TRUE)
  
  # Combine all results
  MECR_results_df = bind_rows(MECR_results_list)
  
  # Cleanup: switch back to sequential mode
  plan("sequential")
  
  return(MECR_results_df)
}
