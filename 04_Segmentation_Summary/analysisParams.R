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
LITERAL_TITLE = "Segmentation summary"

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

# Proseg hyperparameters selected as the optimal settings
OPTIMAL_PROSEG_PARAMS = list("TMA002" = list(VOXEL_LAYERS = 20,
                                             NUC_ASSIGNMENT = 0,
                                             DIFFUSION_PROBABILITY = 0.2),
                             "TMA004" = list(VOXEL_LAYERS = 15,
                                             NUC_ASSIGNMENT = 0.1,
                                             DIFFUSION_PROBABILITY = 0.2))

# Proseg segmentation output directory
PROSEG_SEGMENTATION_OUTPUT_DIRECTORY = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/02_Cell_Segmentation/"

# Add path to proseg segmentation for each TMA in metadata
TMA_METADATA[["TMA002"]][["proseg_segmentation"]][["data_dir"]] = file.path(PROSEG_SEGMENTATION_OUTPUT_DIRECTORY, 
                                                                            "01_TMA002", 
                                                                            "output-XETG00445__0055723__TMA002_TEMOIN__20250708__153012",
                                                                            "proseg3_hyperparameter_search",
                                                                            paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["VOXEL_LAYERS"]]),
                                                                            paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["NUC_ASSIGNMENT"]]),
                                                                            paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["DIFFUSION_PROBABILITY"]]),
                                                                            "04_Proseg_Xenium",
                                                                            "outs")
TMA_METADATA[["TMA004"]][["proseg_segmentation"]][["data_dir"]] = file.path(PROSEG_SEGMENTATION_OUTPUT_DIRECTORY, 
                                                                            "02_TMA004", 
                                                                            "output-XETG00445__0055738__TMA004_TEMOIN__20250708__153012",
                                                                            "proseg3_hyperparameter_search",
                                                                            paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["VOXEL_LAYERS"]]),
                                                                            paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["NUC_ASSIGNMENT"]]),
                                                                            paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["DIFFUSION_PROBABILITY"]]),
                                                                            "04_Proseg_Xenium",
                                                                            "outs")
                                                                            
# Prepare Seurat output path
TMA_METADATA[["TMA002"]][["proseg_segmentation"]][["output_path"]] = file.path(PATH_ANALYSIS_OUTPUT,
                                                                               "proseg_segmentation_rds",
                                                                               "TMA002",
                                                                               paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["VOXEL_LAYERS"]]),
                                                                               paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["NUC_ASSIGNMENT"]]),
                                                                               paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["DIFFUSION_PROBABILITY"]]))
TMA_METADATA[["TMA004"]][["proseg_segmentation"]][["output_path"]] = file.path(PATH_ANALYSIS_OUTPUT,
                                                                               "proseg_segmentation_rds",
                                                                               "TMA004",
                                                                               paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["VOXEL_LAYERS"]]),
                                                                               paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["NUC_ASSIGNMENT"]]),
                                                                               paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["DIFFUSION_PROBABILITY"]]))

# Define TMA to process
TMA_COVERED = c("TMA002", "TMA004")

##############################
## 2 QUALITY CONTROL PARAMS ##
##############################

# Quality thresholds to tag low quality cells
CELL_MINIMUM_FEATURES = 10
CELL_MAXIMUM_AREA = 500
QUALITY_CONTROL_COLOR_PALETTE = c("conserved" = "#b03a2e", 
                                  "low_feature" = "gold1", 
                                  "high_area" = "green4", 
                                  "low_features_high_area" = "darkviolet")

#############
## 03 MECR ##
#############

# Path to pre-computed MECR genes
MECR_GENES_PER_CELL_TYPE = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/01_Find_Mutually_Exclusive_Gene_Pairs/mutually_exclusive_gene_pairs.xlsx"

# MECR computation params
BATCH_SIZE = 500
CORES = 15

##########################
## 04 CLUSTER TMA CELLS ##
##########################

# General names
CLUSTER_LAYER_PREFIX = "TMA_CLUSTER" 
UMAP_LAYER_NAME = "TMA_UMAP"

# Clustering dimensions
DIM_PCA = 50
DIM_USE = 30

# CLustering resolution
RES_RANGE = c(0.1, 0.4)

############################################
## 05 CONDITION SPECIFIC UMAP PROJECTIONS ##
############################################

# Gene signatures per cell type
PERCENT_EXPRESSION_CELL_1ST_THRESHOLD = 90
ADD_MODULE_SCORE_CTRL = 20

# Find all markers step
LOGFC_THRESHOLD = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 9
