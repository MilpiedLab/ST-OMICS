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
LITERAL_TITLE = "Proseg hyperparameter search"

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"
ANALYSIS_STEP_NAME = basename( WORKING_DIR)


# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add here all the constants required for the analysis

################################################################################
##                               01 PREPARE DATA                              ##
################################################################################

#################################################
## 1. RETRIEVE MUTUALLY EXCLUSIVE GENE PAIRS ##
#################################################

# Mutually exclusive gene pairs
MUTUALLY_EXCLUSIVE_GENE_PAIRS = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/01_Find_Mutually_Exclusive_Gene_Pairs/mutually_exclusive_gene_pairs.xlsx"

#########################################
## 2. RETRIEVE PATHS TO PROSEG OUTPUTS ##
#########################################

# Define Xenium condition name
XENIUM_CONDITION_NAME = "Xenium_segmentation"

# Define hyperparameters values (ordered in folder hierarchy order)
hyperparameters_values = list("voxel_layers" = c(10, 15, 20),
                              "nuc_assignment" = c(0, 0.1),
                              "diffusion_probability" = c(0.2, 0.3, 0.4))

# Define path to hyperparameter search Xenium outputs
PATH_HYPERPARAMETER_SEARCH_XENIUM_OUTPUT = list("TMA002" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/02_Cell_Segmentation/01_TMA002/output-XETG00445__0055723__TMA002_TEMOIN__20250708__153012/proseg3_hyperparameter_search/",
                                                "TMA004" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/02_Cell_Segmentation/02_TMA004/output-XETG00445__0055738__TMA004_TEMOIN__20250708__153012/proseg3_hyperparameter_search/")

########################################################
## 3. RETRIEVE PROSEG RUNS AND BASE XENIUM STATISTICS ##
########################################################

# Define TMA to process
TMA_COVERED = c("TMA002", "TMA004")

# Set relative path to proseg outputs
PROSEG_OUTPUTS_RELATIVE_PATHS = "04_Proseg_Xenium/outs"

# Quality thresholds to tag low quality cells
CELL_MINIMUM_FEATURES = 10
CELL_MAXIMUM_AREA = 500
QUALITY_CONTROL_COLOR_PALETTE = c("conserved" = "#b03a2e", 
                                  "low_feature" = "gold1", 
                                  "high_area" = "green4", 
                                  "low_features_high_area" = "darkviolet")

# MECR computation
BATCH_SIZE = 20
CORES = 7

################################################################################
##                              03 PLOT CORRELATIONS                          ##
################################################################################

# Define statistics on which effect will be observed
RESPONSE_VARIABLES = c("nb_cells", 
                       "nb_conserved", 
                       "nb_low_feature", 
                       "nb_high_area", 
                       "nb_low_features_high_area",
                       "mean_MECR")
