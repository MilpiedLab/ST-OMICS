################################################################################
##                             READ XENIUM FUNCTIONS                          ##
################################################################################

######################
## READ XENIUM DATA ##
######################

#' This function aims to read output data from Xenium.
#' 
#' @param data_directory Path to directory containing Xenium output files.
#' @param outs Features-related data to retrieve from Xenium outputs.
#' @param type Cell-related data to retrieve from Xenium outputs.
#' @param mols_qv_threshold Quality threshold to apply on molecules.
#' 
#' @return Required outputs as a list of matrices and dataframes.
#' 
ReadXenium_new = function(data_directory, 
                          outs = c("matrix", "microns"), 
                          type = c("centroids", "segmentations", "cell_metadata"),
                          mols_qv_threshold = 20)
{
  ##############################################
  ## DEFINE DATA TO RECOVER FROM XENIUM FILES ##
  ##############################################
  
  # Cell-related data
  type = match.arg(arg = type,
                   choices = c("centroids", "segmentations", "cell_metadata"),
                   several.ok = TRUE)
  
  # Features-related data
  outs = match.arg(arg = outs,
                   choices = c("matrix", "microns"),
                   several.ok = TRUE)
  
  # Define outputs
  outs = c(outs, type)
  
  ######################
  ## LOAD XENIUM DATA ##
  ######################
  
  # Retrieve each required data type
  data = setNames(lapply(X = outs, FUN = function(otype){
    
    # Depending on the otype value, retrieve specific data
    switch(EXPR = otype, matrix = {
      
      ###########################
      ## RETRIEVE COUNT MATRIX ##
      ###########################
      
      # Retrieve feature counts
      matrix = suppressWarnings(Read10X(data.dir = file.path(data_directory, 
                                                             "cell_feature_matrix/")))
      
      return(matrix)
      
    }, centroids = {
      
      #############################
      ## RETRIEVE CELL CENTROIDS ##
      #############################
      
      # Retrieve cell centroids
      cell_centroids_df = as.data.frame(read_parquet(file.path(data_directory, 
                                                               "cells.parquet")),
                                        stringsAsFactors = FALSE)
      
      # Keep useful columns
      cell_centroids_df = data.frame(x = cell_centroids_df$x_centroid, 
                                     y = cell_centroids_df$y_centroid, 
                                     cell = cell_centroids_df$cell_id, 
                                     stringsAsFactors = FALSE)
      
      return(cell_centroids_df)
      
    }, cell_metadata = {
      
      #############################
      ## RETRIEVE CELL CENTROIDS ##
      #############################
      
      # Retrieve cell centroids
      cell_metadata_df = as.data.frame(read_parquet(file.path(data_directory, 
                                                              "cells.parquet")),
                                       stringsAsFactors = FALSE)
      
      # Keep useful columns
      cell_metadata_df = data.frame(cell = cell_metadata_df$cell_id,
                                    transcript_count = cell_metadata_df$transcript_count,
                                    control_probe_counts = cell_metadata_df$control_probe_counts,
                                    genomic_control_counts = cell_metadata_df$genomic_control_counts,
                                    control_codeword_counts = cell_metadata_df$control_codeword_counts,
                                    unassigned_codeword_counts = cell_metadata_df$unassigned_codeword_counts,
                                    deprecated_codeword_counts = cell_metadata_df$deprecated_codeword_counts,
                                    total_counts = cell_metadata_df$total_counts,
                                    cell_area = cell_metadata_df$cell_area,
                                    nucleus_area = cell_metadata_df$nucleus_area,
                                    nucleus_count = cell_metadata_df$nucleus_count,
                                    segmentation_method = cell_metadata_df$segmentation_method,
                                    x_centroid = cell_metadata_df$x_centroid, 
                                    y_centroid = cell_metadata_df$y_centroid,
                                    stringsAsFactors = FALSE)
      
      return(cell_metadata_df)
      
    }, segmentations = {
      
      ##############################
      ## RETRIEVE CELL BOUNDARIES ##
      ##############################
      
      # Retrieve cell boundaries
      cell_boundaries_df = as.data.frame(read_parquet(file.path(data_directory, 
                                                                "cell_boundaries.parquet")),
                                         stringsAsFactors = FALSE)
      
      # Keep useful columns
      cell_boundaries_df = data.frame(cell = cell_boundaries_df$cell_id,
                                      x = cell_boundaries_df$vertex_x, 
                                      y = cell_boundaries_df$vertex_y, 
                                      stringsAsFactors = FALSE)
      
      return(cell_boundaries_df)
      
    }, microns = {
      
      ######################################
      ## RETRIEVE TRANSCRIPTS COORDINATES ##
      ######################################
      
      # Retrieve transcripts coordinates
      transcripts = as.data.frame(read_parquet(file.path(data_directory, 
                                                         "transcripts.parquet")),
                                  stringsAsFactors = FALSE)
      
      # Filter transcripts based on molecule quality
      transcripts = subset(transcripts, qv >= mols_qv_threshold)
      
      # Keep useful columns
      transcripts = data.frame(x = transcripts$x_location, 
                               y = transcripts$y_location,
                               gene = transcripts$feature_name,
                               stringsAsFactors = FALSE)
      
      return(transcripts)
      
    }, stop("Unknown Xenium input type: ", otype))
    
  }), outs)
  
  return(data)
}

######################
## LOAD XENIUM DATA ##
######################

#' This function aims to load Xenium data from the output directory.
#' 
#' @param data_directory Path to directory containing Xenium output files.
#' @param fov Name of the field of view (FOV) that will be stored in the 
#' Seurat object
#' @param assay Name of the Xenium assay that will be created within the
#' Seurat object.
#' 
#' @return A Seurat object with the specified Xenium data loaded in the object.
#' 
LoadXenium_new = function (data_directory, 
                           fov = "fov", 
                           assay = "Xenium")
{
  # Retrieve Xenium data in a list with the required data
  data = ReadXenium_new(data_directory = data_directory, type = c("centroids", 
                                                                  "segmentations",
                                                                  "cell_metadata"), )
  
  # Create cell-related data
  segmentations_data = list(centroids = CreateCentroids(data$centroids), 
                            segmentation = CreateSegmentation(data$segmentations))
  
  # Create FOV associated to the object
  coords = CreateFOV(coords = segmentations_data, type = c("segmentation", 
                                                           "centroids"), molecules = data$microns, assay = assay)
  
  # Create Seurat object
  xenium_obj = CreateSeuratObject(counts = data$matrix[["Gene Expression"]], 
                                  assay = assay)
  
  # Create additional assays
  if ("Blank Codeword" %in% names(data$matrix)) {
    xenium_obj[["BlankCodeword"]] = CreateAssayObject(counts = data$matrix[["Blank Codeword"]])
  } else if ("Unassigned Codeword" %in% names(data$matrix)) {
    xenium_obj[["BlankCodeword"]] = CreateAssayObject(counts = data$matrix[["Unassigned Codeword"]])
  }
  xenium_obj[["ControlCodeword"]] = CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium_obj[["ControlProbe"]] = CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  
  # Add cell metadata to Seurat object
  xenium_obj = AddMetaData(xenium_obj, metadata = data$cell_metadata)
  
  # Store FOV
  xenium_obj[[fov]] = coords
  
  return(xenium_obj)
}
