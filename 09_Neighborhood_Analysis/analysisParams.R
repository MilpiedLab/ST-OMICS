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
LITERAL_TITLE = "Spatial domains analysis with BANKSY"

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

# Define path to annotated data
PATH_ANNOTATED_MERGED_TMA_DATA = "/mnt/DOSI/PMLAB/BIOINFO/Project/STOmics/Xenium_STOmics_Final_Workflow/05_Output/08a_Manual_Cell_Annotation/manually_annotated_merged_TMA.Rds"

# Define colors for level 1 cell types
LEVEL_1_CELL_TYPE_COLORS = setNames(c("#2E5E00",
                                      "#8DDB35",
                                      "#731073",
                                      "#FF7AFF",
                                      "#F5C9FF",
                                      "#FF8E00",
                                      "#FAE625",
                                      "#3570F0",
                                      "#A1DEF0",
                                      "#D9F7FF",
                                      "#B03A2E",
                                      "#FF977D",
                                      "#D9D9D9",
                                      "#3F4045"),
                                    c("T.cells",
                                      "B.cell",
                                      "Myeloid.cell",
                                      "T.cells+NK",
                                      "NK_ILC.cells",
                                      "PC.PB",
                                      "BEC.cell",
                                      "pDC",
                                      "CLEC9AcDC1_8+LAMP3cDC_9.cell",
                                      "LAMP3cDC_9.cell",
                                      "LEC.cell",
                                      "NESC.cell",
                                      "segmentation.error",
                                      "Granulo"))

# Define colors for level 2 cell types
LEVEL_2_CELL_TYPE_COLORS = setNames(c("#B9DAAF",
                                      "#68AF50",
                                      "#305025",
                                      "#FBF48D",
                                      "#FAD850",
                                      "#F6C0F6",
                                      "#EB6FEB",
                                      "#E01FE0",
                                      "#731073",
                                      "#DE877D",
                                      "#B03A2E",
                                      "#B8FFFB",
                                      "#2EFFF5",
                                      "#00D1C7",
                                      "#7DA4FF",
                                      "#2E6DFF",
                                      "#003FD1"),
                                    c("GC.B",
                                      "Non_malignant_Mem_B.cells",
                                      "Non_malignant_activated_B.cells",
                                      "arterialBEC.cell",
                                      "VenousBEC.cell",
                                      "moMac.cell",
                                      "MacC1Q_7",
                                      "CD1EcDC2_3",
                                      "SIGLEC6cDC_4",
                                      "otherLEC.cell",
                                      "collectLEC",
                                      "MRC",
                                      "advSC",
                                      "FDC",
                                      "CD8_T.CD4_T.cell",
                                      "NK",
                                      "CD8_T.cell"))

# Define complete color palette
COMPLETE_CELL_TYPES_COLOR_PALETTE = c(LEVEL_1_CELL_TYPE_COLORS,
                                      LEVEL_2_CELL_TYPE_COLORS)

########################
## 2. BANKSY ANALYSIS ##
########################

## Banksy params
LAMBDA = 0.8
K_GEOM = 15
BANKSY_ASSAY_NAME = sprintf("banksy_lambda_%s_kgeom_%s",
                            gsub("\\.", "", as.character(LAMBDA)),
                            as.character(K_GEOM))

## Clustering parameters
DIM_PCA = 50
DIM_USE = 15
RES_RANGE = c(0.4, 0.6, 0.8, 1)

############################
## 3. CHARACTERIZE SAMPLE ##
############################

# Define number of neatest neighbors to consider when retrieving neighborhood compositions
KNN_NEIGHBORHOOD_COMPOSITION = 6

# Define which layer to use to annotate the neighboring cells
NEIGHBORS_ANNOTATION = "cell_type_lvl2"

# Define sample colors
SAMPLE_COLORS = c("25T011646-03_TEMOIN_TMA002" = "#b03a2e", 
                  "25T011646-03_TEMOIN_TMA004" = "#002FA7", 
                  "25T006268-12_TEMOIN_TMA004" = "#FFE842")

## Divergence thresholds
# Neighborhoods are considered similar at KL < 0.1 and wasserstein < 0.05
JS_THRESHOLD = 0.05
WASSERSTEIN_THRESHOLD = 0.05

## Random sampling
NUMBER_RANDOM_SAMPLING = 999
NUMBER_RANDOM_CELLS = 10000

###################################################
## 4. COMPARE SAMPLES NEIGHBORHOODS COMPOSITIONS ##
###################################################

# Clusters to retrieve from hierarchical clustering
HIERARCHICAL_CLUSTERS_RANGE = c(3, 6, 12, 16)
