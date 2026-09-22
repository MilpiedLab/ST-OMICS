# The aim of this code is to wrap the PASTIS analysis of the provided Xenium 
# TMAs

################################################################################
##                          1. INIT VIRTUAL ENV                               ##
################################################################################

# Python virtual Environment install 
apt-get install -y python3-pip 
pip install --upgrade pip 
pip install virtualenv 

# Create the .snakemake folder with the virtual env using python3 (to do once per virtualenv) 
virtualenv -p /usr/bin/python3 ~/.snakemake 

# Activate the virtualenv (everytime you want to use the virtual environment): 
source ~/.snakemake/bin/activate 

# Install downgraded pulp version
pip install pulp==2.7.0

# Install snakemake
pip install snakemake

################################################################################
##                           2. PASTIS PIPELINE                               ##
################################################################################

# Retrieve params related to launching PASTIS analysis
source 91_launch_PASTIS/params_proseg3_hyperparameter_search.env

# Launch PASTIS analysis for TMA002
snakemake -ps $TMA002_SNAKEFILE_PATH --use-singularity --singularity-args "-B $TMA002_VOLUMES" --cores $TMA002_NUMBER_CORES

# Launch PASTIS analysis for TMA004
snakemake -ps $TMA004_SNAKEFILE_PATH --use-singularity --singularity-args "-B $TMA004_VOLUMES" --cores $TMA004_NUMBER_CORES

################################################################################
##                             3. END ANALYSIS                                ##
################################################################################

# Deactivate the virtualenv once finished to use it: 
deactivate
