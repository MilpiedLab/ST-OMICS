###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#

######## CONSTANTS ON THE ANALYSIS STEP IDENTITY
######## Modify the ANALYSIS_STEP_NAME and LITERAL_TITLE variable

# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "Segmentation summary - PASTIS"

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"
ANALYSIS_STEP_NAME = basename( WORKING_DIR)


# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add here all the constants required for the analysis

#####################
## 1. PREPARE DATA ##
#####################

# Path to preprocessed merged TMA data
PATH_MERGED_TMA_DATA = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/07_Integrate_TMAs/integrated_merged_TMA.Rds"

# Define path to reference single cell data
REFERENCE_SC_DATA_PATH = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/01_Find_Mutually_Exclusive_Gene_Pairs/reference_single_cell_data.rds"

# Define path to cell type hierarchy
PATH_REFERENCE_CELL_TYPE_HIERARCHY = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/00_RawData/Cell_Type_Hierarchy/Reference_Cell_Types.csv"

# Define available cores for multi-threading
WORKERS = 4

###################
## 2. CLUSTERING ##
###################

# Harmony assay name
HARMONY_ASSAY_NAME = "harmony"

# UMAP
DIM_PCA = 50
HARMONY_DIM_USE = 20
HARMONY_UMAP_LAYER_NAME = "HARMONY_DATA_UMAP"

# Leiden clustering meta data and resolution
HARMONY_CLUSTER_LAYER_NAME_PREFIX = "HARMONY_CLUSTER"
HARMONY_RESOLUTION_RANGE = c(0.1, 0.3, 0.6)

# Gene signatures per cell type
PERCENT_EXPRESSION_CELL_1ST_THRESHOLD = 90
ADD_MODULE_SCORE_CTRL = 20

# Find all markers step
LOGFC_THRESHOLD = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 9

##########################
## 3. SCTYPE ANNOTATION ##
##########################

# Define path to cell type markers used for sctype annotation
CELL_TYPE_MARKERS_DB_PATH = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/01_Reference/scType_db/cell_type_markers_db.xlsx"

# Define tissue type used for sctype annotation
TISSUE_TYPE = "Immune system"

# Define layer prefix to store sctype annotation in seurat object
SCTYPE_LAYER_PREFIX = "SCTYPE_ANNOTATION"

###########################
## 4. SINGLER ANNOTATION ##
###########################

# Define resolution at which level1 singleR annotation will be made
SINGLER_LEVEL1_RESOLUTION = 1.5

# Covariable that needs to be considered for the integration
HARMONY_INTEGRATION_COVARIABLE = "TMA_id"

# Define layer name with singleR lvl1 annotation
SINGLER_LVL1_ANNOTATIONS_LAYER_NAME = "singleR_lvl1_cell_type"

# Clustering dimensions
SINGLER_LVL2_DIM_PCA = list("B.cell" = 50,
                            "LAMP3cDC_9" = 50,
                            "NK_ILC.cells" = 50,
                            "NESC.cell" = 50,
                            "T.cells" = 50,
                            "pDC" = 50,
                            "BEC.cell" = 50,
                            "LEC.cell" = 50,
                            "Myeloid.cell" = 50,
                            "Granulo" = 50)
SINGLER_LVL2_DIM_USE = list("B.cell" = 20,
                            "LAMP3cDC_9" = 20,
                            "NK_ILC.cells" = 20,
                            "NESC.cell" = 20,
                            "T.cells" = 20,
                            "pDC" = 20,
                            "BEC.cell" = 20,
                            "LEC.cell" = 20,
                            "Myeloid.cell" = 20,
                            "Granulo" = 20)

# General names
SINGLER_LVL2_CLUSTER_LAYER_PREFIX = "SINGLER_LVL2_DATA_CLUSTER" 
SINGLER_LVL2_UMAP_LAYER_NAME = "SINGLER_LVL2_DATA_UMAP"

# Define common resolution at which cluster cell subsets for lvl2 annotation
SINGLER_LVL2_RESOLUTION = list("B.cell" = 0.3,
                               "LAMP3cDC_9" = 0.3,
                               "NK_ILC.cells" = 0.3,
                               "NESC.cell" = 0.3,
                               "T.cells" = 0.3,
                               "pDC" = 0.3,
                               "BEC.cell" = 0.3,
                               "LEC.cell" = 0.3,
                               "Myeloid.cell" = 0.3,
                               "Granulo" = 0.3)

# Define layer name with singleR lvl2 annotation
SINGLER_LVL2_ANNOTATIONS_LAYER_NAME = "singleR_lvl2_cell_type"
