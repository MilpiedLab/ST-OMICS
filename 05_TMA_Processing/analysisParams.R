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
LITERAL_TITLE = "TMA processing"

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

# Define default name for proseg segmentation
PROSEG_SEGMENTATION_NAME = "proseg_segmentation"

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
TMA_METADATA[["TMA002"]][[PROSEG_SEGMENTATION_NAME]][["data_dir"]] = file.path(PROSEG_SEGMENTATION_OUTPUT_DIRECTORY, 
                                                                               "01_TMA002", 
                                                                               "output-XETG00445__0055723__TMA002_TEMOIN__20250708__153012",
                                                                               "proseg3_hyperparameter_search",
                                                                               paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["VOXEL_LAYERS"]]),
                                                                               paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["NUC_ASSIGNMENT"]]),
                                                                               paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA002"]][["DIFFUSION_PROBABILITY"]]),
                                                                               "04_Proseg_Xenium",
                                                                               "outs")
TMA_METADATA[["TMA004"]][[PROSEG_SEGMENTATION_NAME]][["data_dir"]] = file.path(PROSEG_SEGMENTATION_OUTPUT_DIRECTORY, 
                                                                               "02_TMA004", 
                                                                               "output-XETG00445__0055738__TMA004_TEMOIN__20250708__153012",
                                                                               "proseg3_hyperparameter_search",
                                                                               paste0("voxel_layers", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["VOXEL_LAYERS"]]),
                                                                               paste0("nuc_assignment", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["NUC_ASSIGNMENT"]]),
                                                                               paste0("diffusion_probability", OPTIMAL_PROSEG_PARAMS[["TMA004"]][["DIFFUSION_PROBABILITY"]]),
                                                                               "04_Proseg_Xenium",
                                                                               "outs")

# Define TMA to process
TMA_COVERED = c("TMA002", "TMA004")

#########################
## 02_SPLIT_TMA PARAMS ##
#########################

# Parameters to optimize the sample splitting on the TMA
SPLIT_TMA_PARAMS = list("TMA002" = list(RADIUS_THRESHOLD = 200,
                                        NEIGHBORS_NUMBER = 100,
                                        SAMPLE_MINIMUM_SIZE = 200,
                                        NUMBER_SAMPLE_ROWS = 4,
                                        NUMBER_SAMPLE_COLUMNS = 4),
                        "TMA004" = list(RADIUS_THRESHOLD = 200,
                                        NEIGHBORS_NUMBER = 100,
                                        SAMPLE_MINIMUM_SIZE = 200,
                                        NUMBER_SAMPLE_ROWS = 4,
                                        NUMBER_SAMPLE_COLUMNS = 4))

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

# Sample quality control thresholds
SAMPLE_MINIMUM_AREA = 750000
SAMPLE_MINIMUM_CONSERVED_CELLS_PERCENTAGE = 80
SAMPLE_MINIMUM_TRANSCRIPT_DENSITY = 2
