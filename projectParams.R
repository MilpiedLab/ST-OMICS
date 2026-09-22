###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples and several analysis steps
###############################################################################


#### General
#### Variables that contains general information on the project and experiment

# A global description of the projet
# It will appear at the beginning of all analysis reports
# Example : "Study of Mouse Memory B Cell infected by Flu"
GLOBAL_DESCRIPTION = "<PUT A GLOBAL DESCRIPTION HERE>"

# The name of you scientific group as it appear in the DOSI folder
# Example : "MGLAB"
SCIENTIFIC_GROUP = "PMLAB"

# The name of the scientific projet. 
# It must be the same name as the name of the project folder you defined
# Example : "moFluMemB"
SCIENTIFIC_PROJECT_NAME = "STOmics"


#### Automatic definition of path
#### You may not modify the following section where standard path are defined in
#### constants so to be used in your code. Have a look at the declared variable
#### and use it in your code.

#### Input / Output

# This is the path of the global project folder
PATH_PROJECT = file.path( "/mnt/DOSI", 
                          SCIENTIFIC_GROUP,
                          "BIOINFO", 
                          "Project",
                          SCIENTIFIC_PROJECT_NAME)

#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



