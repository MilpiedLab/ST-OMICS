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

#############################
## 2. CELL ANNOTATION LVL1 ##
#############################

# Define layer name in which cell type lvl1 annotations will be stored
LEVEL1_ANNOTATION_LAYER_NAME = "cell_type_lvl1"

# Covariable that needs to be considered for the integration
HARMONY_INTEGRATION_COVARIABLE = "TMA_id"

# Harmony assay name
HARMONY_ASSAY_NAME = "harmony"

# Harmony UMAP
HARMONY_UMAP_LAYER_NAME = "HARMONY_DATA_UMAP"

# Define resolutions range from which clusters will be computed
LEVEL_1_RES_RANGE = c(0.1, 0.3, 0.6)

# Define resolution at which lvl1 cell annotation will be made
LEVEL_1_ANNOTATION_RES = 0.3

# Clustering
LEVEL_1_DIM_PCA = 50
HARMONY_DIM_USE = 20
LEVEL_1_LAYER_PREFIX = "LEVEL_1_CLUSTER_LAYER"

# Find all markers step
LOGFC_THRESHOLD = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 9

# Define how many top markers are used to find corresponding reference cluster
TOP_N_ANNOTATION_MARKERS_LVL1 = 10

# Define how to annotate clusters (depending on resolution chosen): clustering resolution --> cluster index --> cell type annotation
LEVEL_1_TYPING_METADATA = list("0.1" = list("0" = "T.cells",
                                            "1" = "B.cell",
                                            "2" = "segmentation.error",
                                            "3" = "Myeloid.cell",
                                            "4" = "PC.PB",
                                            "5" = "BEC.cell",
                                            "6" = "pDC",
                                            "7" = "B.cell"),
                               "0.3" = list("0" = "T.cells",
                                            "1" = "B.cell",
                                            "2" = "segmentation.error",
                                            "3" = "Myeloid.cell",
                                            "4" = "PC.PB",
                                            "5" = "BEC.cell",
                                            "6" = "T.cells+NK",
                                            "7" = "pDC",
                                            "8" = "B.cell",
                                            "9" = "LEC.cell",
                                            "10" = "CLEC9AcDC1_8+LAMP3cDC_9.cell",
                                            "11" = "T.cells"))

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

####################################
## 3. PROCESS SEGMENTATION ERRORS ##
####################################

# Segmentation error subset clustering params
SEGMENTATION_ERROR_SUBSET_RES_RANGE = c(0.3)
SEGMENTATION_ERROR_SUBSET_ANNOTATION_RES = 0.3
SEGMENTATION_ERROR_SUBSET_DIM_PCA = 50
SEGMENTATION_ERROR_SUBSET_DIM_USE = 30
SEGMENTATION_ERROR_SUBSET_UMAP_LAYER_NAME = "SEGMENTATION_ERROR_SUBSET_UMAP"
SEGMENTATION_ERROR_SUBSET_LAYER_PREFIX = "SEGMENTATION_ERROR_SUBSET_CLUSTER_LAYER"

# Define how to annotate clusters obtained from segmentation error subset
SEGMENTATION_ERROR_SUBSET_TYPING_METADATA = list("0" = "segmentation.error",
                                                 "1" = "T.cells",
                                                 "2" = "B.cell",
                                                 "3" = "NESC.cell",
                                                 "4" = "LAMP3cDC_9.cell",
                                                 "5" = "NK_ILC.cells",
                                                 "6" = "pDC",
                                                 "7" = "LEC.cell",
                                                 "8" = "PC.PB")

#############################
## 4. CELL ANNOTATION LVL2 ##
#############################

# Define layer name in which cell type lvl2 annotations will be stored
LEVEL2_ANNOTATION_LAYER_NAME = "cell_type_lvl2"

# Define how to subset data for lvl2 annotation
LVL_1_CELL_TYPES_SUBSETS = list("1" = c("B.cell"),
                                "2" = c("BEC.cell"),
                                "3" = c("Myeloid.cell", "CLEC9AcDC1_8", "LAMP3cDC_9.cell"),
                                "4" = c("LEC.cell"),
                                "5" = c("NESC.cell"),
                                "6" = c("T.cells", "T.cells+NK", "NK_ILC.cells"))

# Define cell types in reference to consider for comparision
LVL_1_REFERENCE_CELL_TYPES_SUBSETS = list("1" = c("B.cell"),
                                          "2" = c("BEC.cell"),
                                          "3" = c("Myeloid.cell", "LAMP3cDC_9.cell"),
                                          "4" = c("LEC.cell"),
                                          "5" = c("NESC.cell"),
                                          "6" = c("T.cells", "NK_ILC.cells"))

# Clustering
LEVEL_2_UMAP_LAYER_NAME = "LEVEL_2_UMAP"
LEVEL_2_LAYER_PREFIX = "LEVEL_2_CLUSTER_LAYER"
LEVEL_2_DIM_PCA = list("1" = 50,
                       "2" = 50,
                       "3" = 50,
                       "4" = 50,
                       "5" = 50,
                       "6" = 50)
LEVEL_2_DIM_USE = list("1" = 20,
                       "2" = 20,
                       "3" = 30,
                       "4" = 20,
                       "5" = 50,
                       "6" = 20)

# Define resolutions to cover to characterize lvl2 cell types
LEVEL_2_ANNOTATION_RES_RANGE = list("1" = c(0.1),
                                    "2" = c(0.1),
                                    "3" = c(0.1),
                                    "4" = c(0.1),
                                    "5" = c(0.3),
                                    "6" = c(0.1))

# Define output of lvl2 cell types related data
CELL_TYPE_LVL2_SUBCLUSTER_OUTPUT_DIR = "subclustering_lvl2_cell_types"

# Define resolutions selected for the annotation of level2 cell types
LEVEL_2_ANNOTATION_RES = list("1" = 0.1,
                              "2" = 0.1,
                              "3" = 0.1,
                              "4" = 0.1,
                              "5" = 0.3,
                              "6" = 0.1)

# Define cell types associated to clusters obtained with resolutions set right above
LEVEL_2_TYPING_METADATA = list("1" = list("0.1" = list("0" = "Non_malignant_Mem_B.cells",
                                                       "1" = "GC.B",
                                                       "2" = "Non_malignant_activated_B.cells",
                                                       "3" = "PC.PB",
                                                       "4" = "Non_malignant_Mem_B.cells",
                                                       "5" = "Non_malignant_activated_B.cells")),
                               "2" = list("0.1" = list("0" = "arterialBEC.cell",
                                                       "1" = "VenousBEC.cell",
                                                       "2" = "arterialBEC.cell")),
                               "3" = list("0.1" = list("0" = "moMac.cell",
                                                       "1" = "MacC1Q_7",
                                                       "2" = "CD1EcDC2_3",
                                                       "3" = "SIGLEC6cDC_4"),
                                          "0.3" = list("0" = "moMac.cell",
                                                       "1" = "MacC1Q_7",
                                                       "2" = "Unknown",
                                                       "3" = "SIGLEC6cDC_4",
                                                       "4" = "SIGLEC6cDC_4",
                                                       "5" = "SIGLEC6cDC_4",
                                                       "6" = "moMac.cell",
                                                       "7" = "MacC1Q_7",
                                                       "8" = "Unknown",
                                                       "9" = "moMac.cell",
                                                       "10" = "moMac.cell")),
                               "4" = list("0.1" = list("0" = "otherLEC.cell",
                                                       "1" = "collectLEC")),
                               "5" = list("0.3" = list("0" = "MRC",
                                                       "1" = "advSC",
                                                       "2" = "FDC")),
                               "6" = list("0.1" = list("0" = "CD8_T.CD4_T.cell",
                                                       "1" = "NK",
                                                       "2" = "CD8_T.cell",
                                                       "3" = "CD8_T.CD4_T.cell",
                                                       "4" = "CD8_T.cell")))

# Define how many top markers are used to find corresponding reference cluster
TOP_N_ANNOTATION_MARKERS_LVL2 = 5
