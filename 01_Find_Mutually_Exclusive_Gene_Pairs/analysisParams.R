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
LITERAL_TITLE = "Find mutually exclusive gene pairs"

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
## 01 PREPARE DATA ##
#####################

# Define path to reference single cell data
REFERENCE_SC_SEURAT_INPUT = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/00_RawData/Reference_Single_Cell/complete.atlas.integrated.seu.rds"

######################################
## 02 FIND MUTUALLY EXCLUSIVE GENES ##
######################################

# Find markers step
LOGFC_THRESHOLD = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 9
LOGFC_MECR_THRESHOLD = 3

# Minimum % of cells expressing a gene to be considered expressed in this cell type
MIN_EXPRESSION_SPECIFIC_GENE = 0.5

# Maximum % of cells expressing a gene to be considered unexpressed in this cell type
MAX_EXPRESSION_NONSPECIFIC_GENE = 0.05
