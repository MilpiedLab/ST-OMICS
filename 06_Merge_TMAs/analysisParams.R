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

#################################
## 1. DEFINE SAMPLES META DATA ##
#################################

# Define available TMAs
TMA_COVERED = c("TMA002", "TMA004")

# Define paths to splitted TMAs Seurat objects
PATH_SPLITTED_TMA_RDS = list("TMA002" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/05_TMA_Processing/TMA002_splitted_proseg_segmentation.Rds",
                             "TMA004" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/05_TMA_Processing/TMA004_splitted_proseg_segmentation.Rds")

###############################
## 03 QUALITY CONTROL PARAMS ##
###############################

# Quality thresholds to tag low quality cells
CELL_MINIMUM_FEATURES = 10
CELL_MAXIMUM_AREA = 500
QUALITY_CONTROL_COLOR_PALETTE = c("conserved" = "#b03a2e", 
                                  "low_feature" = "gold1", 
                                  "high_area" = "green4", 
                                  "low_features_high_area" = "darkviolet")

############################
## 04 CLUSTER MERGED DATA ##
############################

# Cluster on merged non filtered data
MERGED_UNFILTERED_CLUSTER_LAYER_PREFIX = "MERGED_UNFILTERED_CLUSTER"
MERGED_UNFILTERED_UMAP_LAYER_PREFIX = "MERGED_UNFILTERED_UMAP"

# Clustering dimensions
DIM_PCA = 50
DIM_USE = 30

# Clustering resolution
MERGED_UNFILTERED_RES_RANGE = c(0.1)

############################
## 05 CLUSTER MERGED DATA ##
############################

# General names
CLUSTER_LAYER_PREFIX = "MERGED_DATA_CLUSTER" 
UMAP_LAYER_NAME = "MERGED_DATA_UMAP"

# Clustering resolution
RES_RANGE = c(0.1, 0.3, 0.6)

###########################################
## 06 CLUSTERED MERGED DATA GENE MARKERS ##
###########################################

# Gene signatures per cell type
PERCENT_EXPRESSION_CELL_1ST_THRESHOLD = 90
ADD_MODULE_SCORE_CTRL = 20

# Find all markers step
LOGFC_THRESHOLD = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 9
