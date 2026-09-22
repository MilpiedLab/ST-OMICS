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
PATH_MERGED_TMA_DATA = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/06_Merge_TMAs/filtered_clustered_TMAs.Rds"

# Clustering dimensions
DIM_PCA = 50
DIM_USE = 30

# General names
CLUSTER_LAYER_PREFIX = "MERGED_DATA_CLUSTER" 
UMAP_LAYER_NAME = "MERGED_DATA_UMAP"

# Clustering resolution
RES_RANGE = c(0.1, 0.3, 0.6)

############################
## 2. HARMONY INTEGRATION ##
############################

# Covariable that needs to be considered for the integration
HARMONY_INTEGRATION_COVARIABLE = "TMA_id"

# Harmony assay name
HARMONY_ASSAY_NAME = "harmony"

# UMAP
HARMONY_DIM_USE = 20
HARMONY_UMAP_LAYER_NAME = "HARMONY_DATA_UMAP"

# Leiden clustering meta data and resolution
HARMONY_CLUSTER_LAYER_NAME_PREFIX = "HARMONY_CLUSTER"
HARMONY_RESOLUTION_RANGE = c(0.3)

# Find all markers step
LOGFC_THRESHOLD = 0.1
MIN_PCT = 0.1
TOP_N_MARKERS = 9
