################################################################################
##            GLOBAL DECLARATIONS AND LIBRARIES FOR THE ANALYSIS              ##
################################################################################

################################################################################
##                               1. LIBRARIES                                 ##
################################################################################

############
## 1.1. R ##
############

# General libraries
reticulate::use_virtualenv("/opt/venv/",required = TRUE)
library(reticulate)
library(Seurat)

# Banksy libraries
library(Banksy)
library(SeuratWrappers)

# Plot libraries
library(ggplot2)
library(gridExtra)
library(pander)
library(digest)
library(clustree)

# Characterize neighborhoods libraries
library(RANN)
library(tidyr)
library(dplyr)
library(philentropy)
library(ComplexHeatmap)
library(plotly)
library(circlize)
library(grid)

#################
## 1.1. PYTHON ##
#################

# Import numpy and POT for sliced wasserstein computation
np = import("numpy")
ot = import("ot.sliced")
