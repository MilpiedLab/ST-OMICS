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
options(future.globals.maxSize = 10 * 1024^3)
library(arrow)

# Markers
library(openxlsx)

# Plot
library(ggplot2)
library(cowplot)
library(reshape2)
library(scales)
library(ggridges)
library(gridExtra)
library(pheatmap)
