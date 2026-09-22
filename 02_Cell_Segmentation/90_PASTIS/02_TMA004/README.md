# PASTIS : Pipeline for Analyis of Spatial Transcriptomic with Improved Segmentation

## A pipeline that will help you perform resegmentation using Proseg coupled with basic analysis

### Proseg 

[Proseg](https://github.com/dcjones/proseg) is a cell segmentation method based on cellular Potts model

---

### Snakemake

You need to install Snakemake to run the complete analysis workflow. Use your preferred solution : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

exemple using virtual env

```
pip install --user virtualenv

cd <VIRTUALENV_HOME>
virtualenv -p /bin/python3 .snakemakeenv

source <VIRTUALENV_HOME>/.snakemakeenv/bin/activate

pip3 install snakemake
```

---

### Singularity  

You need to install Singularity v3.x on your system to run the complete analysis. Follow the instructions here : https://docs.sylabs.io/guides/3.0/user-guide/installation.html
Several docker environments stored in a .sif format will be called by snakemake. 

---

### Pipeline steps :
For each sample (or FOV) : 
- Run Xenium ranger nucleus
- Run proseg segmentation
- Run proseg to baysor 
- Run Xenium ranger on the proseg output
- Run a QC step
- Run a cell heterogeneity analysis of the dataset

(Facultative) For all samples (or FOV), if many :
- Perform merge plus integration

---

### PASTIS folder architecture

```
PASTIS
    ├── 00_RawData
    ├── 01_Reference
    ├── 02_Docs
    │       ├── dag_integrate.png
    │       └── docs.md
    ├── 03_Script
    │       ├── 01_Proseg_QC
    │       │       ├── 00_generalDeps.R
    │       │       ├── 01_Proseg_QC.R
    │       │       ├── analysisParams.R
    │       │       ├── launch_reports_compilation.R
    │       │       └── Report_proseg_QC.Rmd
    │       ├── 02_Proseg_analysis
    │       │       ├── 00_generalDeps.R
    │       │       ├── 01_Proseg_analysis.R
    │       │       ├── analysisParams.R
    │       │       ├── launch_reports_compilation.R
    │       │       └── Report_proseg_analysis.Rmd
    │       └── 03_Proseg_integration
    │           ├── 00_generalDeps.R
    │           ├── 01_Proseg_summary_and_merge.R
    │           ├── 02_Proseg_integration.R
    │           ├── analysisParams.R
    │           ├── launch_reports_compilation.R
    │           └── Report_proseg_integration.Rmd
    ├── 04_Workflow
    │       ├── config
    │       │       ├── basic_config.yaml
    │       │       ├── expert_config.yaml
    │       │       └── mandatory_config.yaml
    │       └── snakefile_proseg_pipeline.py
    ├── 05_Output
    └── README.md
```

### Tutorial

An in depth tutorial can be found [here](02_Docs/dos.md)

1. The first step of the pipeline is to do a symbolic link of your raw data (Default Xenium Ranger output) in the 00_RawData folder.

```
cd 00_RawData
ln -s <YOUR_RAW_DATA_PATH> # With feature matrix...
```
2. The second step is to put your cell type associated genes lists file in the 01_Reference folder (facultative).

3. The third step 

This snakemake pipeline is customizable by modifying the different config.yaml files. There is three config files for different degree of customization.

- mandatory_config.yaml

This config file must not be modified. It contain the seed, dockers and the folder names that will be used.

- basic_config.yaml

This config file contain the basic customisation option.
You need to fill the following parameters for the pipeline to work : 
- PATH_RAW_DATA_XENIUM
- PATH_ANALYSIS_OUTPUT
- samples
- MECR_PATH
- proseg: threads
The remaining parameters are help to customize the analysis, description of the effect and potential values are in the file.

- expert_config.yaml

This file contain the more complex list of parameters like the normalisation method, integration method, number of PC to use.

Once the config the are updated you can launch the snakemake dry run by running :

```
snakemake -ps 04_Workflow/snakefile_proseg_pipeline.py \
    --use-singularity --singularity-args "-B /<volume>:/<volume>" \
    --cores 40 -n
```

If everything is correct remove the `-n` to start the pipeline.


