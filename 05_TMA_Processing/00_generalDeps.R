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
options(future.globals.maxSize = 20 * 1024^3)
library(arrow)

# Split TMA
library(RANN)
library(igraph)

# Plot
library(cowplot)
library(reshape2)
library(scales)
library(ggridges)
library(gridExtra)
library(data.table)
