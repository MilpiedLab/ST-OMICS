###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples and several analysis steps
###############################################################################


#### General
#### Variables that contains general information on the experiment

# The experiment name this file is in
# It must be the same name as the name of the experiment folder you defined
# Example : "10x_190712_m_moFluMemB"
EXPERIMENT_NAME = basename( dirname( dirname( WORKING_DIR)))


#### Automatic definition of path
#### You may not modify the following section where standard path are defined in
#### constants so to be used in your code. Have a look at the declared variable
#### and use it in your code.

#### Input / Output

# This is the path of the current experiment this file is in
PATH_EXPERIMENT = file.path( PATH_PROJECT, EXPERIMENT_NAME)

# Those are the path to the main folder : RAW data, REFERENCE data and Output of analysis
PATH_EXPERIMENT_RAWDATA   = file.path( PATH_EXPERIMENT, "00_RawData")
PATH_EXPERIMENT_REFERENCE = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_OUTPUT    = file.path( PATH_EXPERIMENT, "05_Output")


# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_",     
                            EXPERIMENT_NAME
                            )

# Define colors for plotting
COLOR_PALETTE = c("#b03a2e", "#002FA7", "green4", "#6A3D9A", "#FF7F00", "gold1", "palegreen2", "orchid1", "skyblue2", "#FB9A99", 
                  "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "deeppink1", "dodgerblue2", "steelblue4", "darkturquoise", 
                  "green1", "yellow4", "yellow3", "darkorange4", "#E31A1C", "cyan", "darkgoldenrod", "deeppink2", "darkviolet")

# Assay names
SCTRANSFORM_ASSAY_NAME = "SCT"
XENIUM_ASSAY_NAME = "Xenium"

# Set random seed
SEED = 42

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



