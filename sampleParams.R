################################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 
################################################################################

#################################
## 1. DEFINE SAMPLES META DATA ##
#################################

##############
## 1.1. TMA ##
##############

# Define path to raw Xenium data
XENIUM_DATA_PATH = list("TMA002" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/00_RawData/250708_Xenium_152942_20250708_PM_TMA002_TMA004/output-XETG00445__0055723__TMA002_TEMOIN__20250708__153012/",
                        "TMA004" = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/00_RawData/250708_Xenium_152942_20250708_PM_TMA002_TMA004/output-XETG00445__0055738__TMA004_TEMOIN__20250708__153012/")

# Define TMA meta data
TMA_METADATA = list("TMA002" = list("xenium_segmentation" = list("data_dir" = XENIUM_DATA_PATH[["TMA002"]])),
                    "TMA004" = list("xenium_segmentation" = list("data_dir" = XENIUM_DATA_PATH[["TMA004"]])))

########################
## 1.2. PATIENT NAMES ##
########################

# Define all patients found in all TMAs (each patient was done in duplicate)
PATIENTS = c("25T006268-12_TEMOIN",
             "14T017106_TEMOIN",
             "25T011646-03_TEMOIN",
             "14T017106_TEMOIN")

######################
## 1.3. SAMPLE GRID ##
######################

# Define position (grid) of the samples (patient replicates) on TMA
TMA_GRID = list("TMA002" = data.frame("1" = c("25T006268-12_TEMOIN", NA, NA, NA),
                                      "2" = c(NA, "14T017106_TEMOIN.1", NA, NA),
                                      "3" = c(NA, NA, "25T011646-03_TEMOIN", NA),
                                      "4" = c(NA, NA, NA, "14T017106_TEMOIN.2")),
                "TMA004" = data.frame("1" = c("25T006268-12_TEMOIN", NA, NA, NA),
                                      "2" = c(NA, "14T017106_TEMOIN.1", NA, NA),
                                      "3" = c(NA, NA, "25T011646-03_TEMOIN", NA),
                                      "4" = c(NA, NA, NA, "14T017106_TEMOIN.2")))

# Rows / columns on the TMA grid (TMA splitting step)
GRID_DIMENSIONS = list("TMA002" = list(nb_cols = 4,
                                       nb_rows = 4),
                       "TMA004" = list(nb_cols = 4,
                                       nb_rows = 4))
