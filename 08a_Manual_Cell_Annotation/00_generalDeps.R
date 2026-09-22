################################################################################
##            GLOBAL DECLARATIONS AND LIBRARIES FOR THE ANALYSIS              ##
################################################################################

################################################################################
##                               1. LIBRARIES                                 ##
################################################################################

# General libraries
library(digest)
library(pander)
library(dplyr)
library(tidyr)
library(ggplot2)

# Seurat-related
library(Seurat)

# Clustering multi-threading
library(future)
plan("multisession", workers = WORKERS)
options(future.globals.maxSize = 100 * 1024^3)

# Cluster stability
library(clustree)

# Plot
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(ggridges)
library(gridExtra)
library(scattermore)
library(collapsibleTree)

# Markers
library(openxlsx)

# Integration
library(harmony)
