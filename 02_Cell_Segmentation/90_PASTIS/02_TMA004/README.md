# PASTIS : Pipeline for Analysis of Spatial Transcriptomic with Improved Segmentation

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
- Run proseg segmentation
- Run proseg to baysor

---

### PASTIS folder architecture

```
PASTIS
    ├── 02_Docs
    │       └── docs.md
    ├── 04_Workflow
    │       ├── config_proseg3_hyperparameter_search
    │       │       ├── basic_config.yaml
    │       │       ├── expert_config.yaml
    │       │       └── mandatory_config.yaml
    │       └── snakefile_proseg3_hyperparameter_search_pipeline.py
    └── README.md
```

### Tutorial

This snakemake pipeline is customizable by modifying the different config.yaml files. There is three config files for different degree of customization.

- mandatory_config.yaml

This config file must not be modified. It contain the seed, dockers and the folder names that will be used.

- basic_config.yaml

This config file contain the basic customisation option.
You need to fill the following parameters for the pipeline to work : 
- PATH_RAW_DATA_XENIUM
- PATH_ANALYSIS_OUTPUT
- samples
- proseg: threads
The remaining parameters are help to customize the analysis, description of the effect and potential values are in the file.

- expert_config.yaml

This file contain the more complex list of parameters such as Proseg parameters.

Once the config the are updated you can launch the snakemake dry run by running :

```
snakemake -ps 04_Workflow/snakefile_proseg3_hyperparameter_search_pipeline.py \
    --use-singularity --singularity-args "-B /<volume>:/<volume>" \
    --cores 40 -n
```

If everything is correct remove the `-n` to start the pipeline.


