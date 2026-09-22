
## @knitr retrieve_mutually_exclusive_genes

cat(sprintf("%s **Mutually exclusive genes**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("This script aims to compute mutually exclusive gene pairs from single-cell 
reference data. Here, the FindMarkers function is used alongside predefined 
thresholds, described below, to identify mutually exclusive genes for each cell 
type found in the reference dataset. Other methods than the one described below 
may also be considered to define these gene pairs.

Even though four levels of cell type annotation are available for this dataset, 
only the first level will be considered to retrieve mutually exclusive gene pairs. 
While it is feasible to identify exclusive pairs for broad cell types, doing so 
becomes increasingly difficult as the cell type annotations become more refined.

Once cell type–specific genes are identified using the FindMarkers function, 
markers are filtered based on log fold change (logFC), pct.1, and pct.2. This 
approach allows the identification of genes that should be specific to each cell 
type.

A marker gene is considered cell type–specific if it is expressed in more than 
50% of the cells within that cell type and in less than 5% of the other cells. 
Mutually exclusive gene pairs are then derived in subsequent analysis steps by 
computing all heterotypic pairs of these exclusive genes.")
cat("\n\n")

################################################################################
##                     1. RETRIEVE MUTUALLY EXCLUSIVE GENES                   ##
################################################################################

# Activate first level of cell types
reference_sc_data = activateSeuratAssay(reference_sc_data, 
                                        ident = "Celltype_Level_1")

# Create list to store broad cell type markers
cell_type_lvl1_markers = data.frame()
markers_as_strings = data.frame()

# Iterate over broad cell types
for (cell_type in unique(reference_sc_data@meta.data[["Celltype_Level_1"]])){
  
  print(sprintf("Computing %s markers", cell_type))
  
  # Compute markers
  temp_markers = FindMarkers(reference_sc_data,
                             ident.1 = cell_type,
                             features = gene_panel_480,
                             thresh.use = LOGFC_THRESHOLD,
                             min.pct = MIN_PCT,
                             only.pos = TRUE)
  
  # Filter markers based on provided thresholds
  temp_markers = temp_markers[temp_markers$avg_log2FC >= LOGFC_MECR_THRESHOLD, ]
  temp_markers = temp_markers[temp_markers$pct.1 >= MIN_EXPRESSION_SPECIFIC_GENE, ]
  temp_markers = temp_markers[temp_markers$pct.2 <= MAX_EXPRESSION_NONSPECIFIC_GENE, ]
  temp_markers$cell_type = cell_type
  
  # Save markers
  cell_type_lvl1_markers = rbind(cell_type_lvl1_markers, temp_markers)
  markers_as_strings = rbind(markers_as_strings, data.frame(Gene_list = paste(rownames(temp_markers), collapse = ","),
                                                            row.names = cell_type))
  
}

################################################################################
##                      2. PLOT MUTUALLY EXCLUSIVE GENES                      ##
################################################################################

cat("\n\n")
print(htmltools::tagList(DT::datatable(markers_as_strings)))
cat("\n\n")
