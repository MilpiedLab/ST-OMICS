# PASTIS in depth documentation

In this tutorial we will describe how to properly run the PASTIS pipeline, understand the input file and the different output files.

## Table of Contents

[TOCM]

[TOC]

## Folder architecture

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

## Input file 

PASTIS input file is the standard Xenium output directory that look like this (in the new Xenium version the transcripts.csv.gz file is no longer present). The cell_feature_matrix.tar.gz might need to be untarred.

```
<YOUR_RAW_XENIUM_DIRECTORY>
    ├── analysis_summary.html
    ├── analysis.tar.gz
    ├── analysis.zarr.zip
    ├── cell_boundaries.csv.gz
    ├── cell_boundaries.parquet
    ├── cell_feature_matrix
    │    ├── barcodes.tsv.gz
    │    ├── features.tsv.gz
    │    └── matrix.mtx.gz
    ├── cell_feature_matrix.h5
    ├── cell_feature_matrix.tar.gz
    ├── cell_feature_matrix.zarr.zip
    ├── cells.csv.gz
    ├── cells.parquet
    ├── cells.zarr.zip
    ├── experiment.xenium
    ├── gene_panel.json
    ├── metrics_summary.csv
    ├── morphology_focus.ome.tif
    ├── morphology_mip.ome.tif
    ├── morphology.ome.tif
    ├── nucleus_boundaries.csv.gz
    ├── nucleus_boundaries.parquet
    ├── transcripts.csv.gz
    ├── transcripts.parquet
    └── transcripts.zarr.zip
```

## Config file

In PASTIS we use 3 different config files in order to customize the analysis. 
They can be found in :
```   
04_Workflow/
└── config
        ├── basic_config.yaml
        ├── expert_config.yaml
        └── mandatory_config.yaml
```
And are named : 
1. mandatory_config.yaml
2. basic_config.yaml
3. expert_config.yaml

### Mandatory config

This config file contain the seed an the singularity images to use for the different analysis.

### Basic config

In this file we define the basic customization option. Modifying the following parameters is mandatory if you want PASTIS to run properly.

- PATH_RAW_DATA_XENIUM
- PATH_ANALYSIS_OUTPUT
- samples
- proseg: threads

#### PATH_RAW_DATA_XENIUM

Here you will set your **PATH_RAW_DATA_XENIUM** as `00_RawData/<YOUR_RAW_XENIUM_DIRECTORY>`. If you only have one sample, your  &lt;YOUR_RAW_XENIUM_DIRECTORY&gt; will be considered as a sample and **PATH_RAW_DATA_XENIUM** will be equal to `00_RawData/`

#### PATH_ANALYSIS_OUTPUT

**PATH_ANALYSIS_OUTPUT** must be set as `05_Output/<YOUR_ANALYSIS_NAME>` with &lt;YOUR_ANALYSIS_NAME&gt; being the name you want to set for your analysis directory (it is mandatory to set a name)

#### samples

The **sample** field is a coma separated list of your sample name. If you have multiple sample, &lt;SAMPLE_X&gt; must correspond to the subdirectory name in &lt;YOUR_RAW_XENIUM_DIRECTORY&gt;. If you only have one sample, **sample** must be set as &lt;YOUR_RAW_XENIUM_DIRECTORY&gt;

#### proseg threads

Proseg is allowing parallelisation, the **threads** parameter allows you to chose the amount of threads used by the proseg process.

#### Other parameters

The remaining parameters are standard transcriptomic analysis parameters 

### Expert config

This config file contain more complex parameters that should be modified only if you have some experience with single cell RNA sequencing. 

## Snakefile 

:warning: :warning: :warning:
The snakefile should not be modified
:warning: :warning: :warning:

To launch the dry-run of the pipeline launch the following command.

```
snakemake -ps 04_Workflow/snakefile_proseg3_hyperparameter_search_pipeline.py \
    --use-singularity --singularity-args "-B /<volume>:/<volume>" \
    --cores 40 -n
```

The `"-B <volume>:<volume>"` correspond to singularity option and &lt;volume&gt; correspond to the host port and container port e.g. `"-B /home:/home"` this port must be connected to the PASTIS directory and the original &lt;YOUR_RAW_DATA_PATH&gt; 

If you want to start the pipeline and not the dry run you need to remove the `--dry-run` option

## Output

All the output will be in 
```
05_Output/
└── <YOUR_ANALYSIS_NAME>
    └── <SAMPLE_1>
```

### Proseg

The *Proseg* step consist in running the *Proseg* tool on the transcript file, by default the RNA diffusion option is activated. The Proseg output folder contain polygons coordinate in geojson format, cell/transcript association and cell metadata.

```
02_Proseg_run/
├── cell-metadata.csv.gz
├── cell-polygons.geojson.gz
├── cell-polygons-layers.geojson.gz
├── expected-counts.csv.gz
├── transcript-metadata.csv.gz
└── union-cell-polygons.geojson.gz
```

### Proseg to baysor

The *Proseg* to *Baysor* step consist in generating *Baysor* like output that will be interpreted by Xeniumr Ranger

```
03_Proseg_to_Baysor/
├── baysor-cell-polygons.geojson
└── baysor-transcript-metadata.csv
```
