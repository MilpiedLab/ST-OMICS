################################################################################
##                               GENERAL FUNCTIONS                            ##
################################################################################

############################
## ACTIVATE SEURAT OBJECT ##
############################

#' This function aims to update the default assay, clusters layer and UMAP
#' coordinates of a Seurat object in function of provided inputs.
#' 
#' Libraries: Seurat
#' 
#' @param seurat_object Seurat object where default assay/layer/reductions must
#' be changed.
#' @param ident Name of the clusters layer to set as default ("seurat_clusters"
#' layer in metadata).
#' @param active_assay Name of the assay to set as default assay.
#' @param umap_layer_name Name of the UMAP coordinates to set as default
#' reductions. The UMAP coordinates must be stored in metadata.
#' 
#' @return Seurat object with updated default assay/layer/reductions.
#' 
activateSeuratAssay = function(seurat_object, 
                               ident = NULL, 
                               active_assay = NULL, 
                               umap_layer_name = NULL)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  
  # Set ident and active cluster layer
  if (!(is.null(ident))){
    # Verify inputs
    if (!(is.character(ident))) stop("Provided ident is not a character.")
    if (is.null(seurat_object@meta.data[[ident]])) stop("Provided ident does not exist in object meta.data.")
    # Update object
    seurat_object = SetIdent(seurat_object, value = ident)
    seurat_object@meta.data[["seurat_clusters"]] = seurat_object@meta.data[[ident]]
  }
  
  # Set active assay
  if (!(is.null(active_assay))){
    # Verify inputs
    if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
    if (is.null(seurat_object@assays[[active_assay]])) stop("Provided active assay does not exist in object assays.")
    # Update object
    seurat_object@active.assay = active_assay
  }
  
  # Set active UMAP coordinates (create new DimReduc object)
  if (!(is.null(umap_layer_name)) & !(is.null(active_assay))){
    # Verify inputs
    if (!(is.character(active_assay))) stop("Provided active_assay is not a character.")
    if (is.null(seurat_object@assays[[active_assay]])) stop("Provided active assay does not exist in object assays.")
    if (!(is.character(umap_layer_name))) stop("Provided umap_layer_name is not a character.")
    if (is.null(seurat_object@meta.data[[umap_layer_name]])) stop("Provided umap_layer_name does not exist in object meta.data.")
    if (is.null(ncol(seurat_object@meta.data[[umap_layer_name]])) || ncol(seurat_object@meta.data[[umap_layer_name]]) != 2) stop("Provided umap_layer_name does not have 2 coordinates dimensions (x/y)")
    
    # Update object
    UMAP_coord = seurat_object@meta.data[[umap_layer_name]]
    UMAP_coord_mat = as(UMAP_coord, "matrix")
    seurat_object[["umap"]] = CreateDimReducObject(embeddings = UMAP_coord_mat, key = "umap_", assay = active_assay)
  }
  return(seurat_object)
}

###########################
## SAVE UMAP COORDINATES ##
###########################

#' This functions aims to store the default UMAP coordinates in meta.data of the
#' provided Seurat object.
#' 
#' Libraries: Seurat
#' 
#' @param seurat_object Seurat object from which the current default UMAP coords
#' must be saved in meta.data.
#' @param umap_layer_name Meta.data layer name that will contain the saved
#' UMAP coordinates.
#' 
#' @return Updated Seurat object with saved UMAP coordinates.
#' 
storeUMAPembeddings = function(seurat_object, 
                               umap_layer_name)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  if (!(is.character(umap_layer_name))) stop("Provided umap_layer_name is not a character.")
  
  # Store UMAP coordinates in meta.data
  seurat_object@meta.data[[umap_layer_name]] = seurat_object@reductions[["umap"]]@cell.embeddings
  
  return(seurat_object)
}
