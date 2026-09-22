################################################################################
##                             NORMALISATION FUNCTIONS                        ##
################################################################################

##############################
## LOG NORMALISATION (LOGN) ##
##############################

#' This function aims to normalize raw counts from a Seurat object using log
#' normalization.
#' 
#' Libraries: Seurat
#' 
#' @param seurat_object Seurat object from which raw counts will be normalized.
#' @param active_assay Assay on which normalization will be applied.
#' 
#' @return Seurat with normalized assay (normalized counts added to "data").
#' 
log_normalisation = function(seurat_object, active_assay)
{
  # Verify inputs
  if (!inherits(seurat_object, "Seurat")) stop("This function requires a Seurat object.")
  
  # Activate required assay
  seurat_object = activateSeuratAssay(seurat_object, active_assay = active_assay)
  
  # Normalize data
  seurat_object_normalised = NormalizeData(seurat_object, normalization.method = "LogNormalize")
  
  return(seurat_object_normalised)
}
