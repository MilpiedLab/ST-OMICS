
## @knitr prepare_data

cat(sprintf("%s **OVERVIEW AND OBJECTIVES**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")

cat("\n\n")
cat("The aim of this analysis step is to automatically annotate cell clusters and 
individual cells obtained from the Harmony-integrated assay using two distinct 
methods: scType and SingleR, which rely on different principles.

scType performs annotation at the cluster level by evaluating the enrichment of 
known marker genes from curated cell type–specific gene sets. In contrast, SingleR 
operates at the single-cell level and assigns labels by comparing each cell’s 
expression profile to annotated reference cells from a single-cell dataset using 
correlation-based approaches.

These two methods are applied independently to provide alternative automatic 
annotations based on different underlying assumptions.")
cat("\n\n")

################################################################################
#'                              1. LOAD MERGED DATA
#' 
#' Load the merged TMA dataset containing the Seurat object with the 
#' Harmony-integrated assay obtained in step '07_Integrate_TMAs'.
#' 
################################################################################

# Create output directory
dir.create(PATH_ANALYSIS_OUTPUT, 
           recursive = TRUE, 
           showWarnings = FALSE)

# Read merged TMA data
xenium_data = readRDS(PATH_MERGED_TMA_DATA)

# Load cell type hierarchy
reference_cell_type_hierarchy = read.csv(PATH_REFERENCE_CELL_TYPE_HIERARCHY)

# Load reference single cell data for annotation
sc_reference_data = readRDS(file = REFERENCE_SC_DATA_PATH)

# Remove malignant cells as these are not found in current data
sc_reference_data = subset(sc_reference_data, subset = Celltype_Level_2 != "Malignant_B.cell")

################################################################################
#'                          2. DISPLAY CELL TYPE HIERARCHY
################################################################################

cat("\n\n")
cat(sprintf("%s **Cell type hierarchy**", strrep("#", 1)))
cat("\n\n")

cat("\n\n")
collapsibleTree(reference_cell_type_hierarchy,
                hierarchy = c("LEVEL1", "LEVEL2", "LEVEL3", "LEVEL4"))
cat("\n\n")
