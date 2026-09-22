# PASTIS in depth documentation

In this tutorial we will describe how to properly run the PASTIS pipeline, undertand the input file and the different output files.

## Table of Contents

[TOCM]

[TOC]

## Folder architecture

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

You need to generate a symbolic link of your raw data folder in the 00_RawData folder.

```
cd 00_RawData
ln -s <YOUR_RAW_DATA_PATH>
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

:warning: :warning: :warning:
This config file must not be modified.
:warning: :warning: :warning:

This config file contain the seed, the dockers to pull from dockerhub and finally the name of the folder for the different analysis.

### Basic config

In this file we define the basic customization option. Modifying the following parameters is mandatory if you want PASTIS to run properly.

- PATH_RAW_DATA_XENIUM
- PATH_ANALYSIS_OUTPUT
- samples
- FINAL_STEP
- proseg: threads
- r_analysis_param : MECR_PATH

#### PATH_RAW_DATA_XENIUM

Here you will set your **PATH_RAW_DATA_XENIUM** as `00_RawData/<YOUR_RAW_XENIUM_DIRECTORY>`. If you only have one sample, your  &lt;YOUR_RAW_XENIUM_DIRECTORY&gt; will be considered as a sample and **PATH_RAW_DATA_XENIUM** will be equal to `00_RawData/`

#### PATH_ANALYSIS_OUTPUT

**PATH_ANALYSIS_OUTPUT** must be set as `05_Output/<YOUR_ANALYSIS_NAME>` with &lt;YOUR_ANALYSIS_NAME&gt; being the name you want to set for your analysis directory (it is mandatory to set a name)

#### samples

The **sample** field is a coma separated list of your sample name. If you have multiple sample, &lt;SAMPLE_X&gt; must correspond to the subdirectory name in &lt;YOUR_RAW_XENIUM_DIRECTORY&gt;. If you only have one sample, **sample** must be set as &lt;YOUR_RAW_XENIUM_DIRECTORY&gt;

#### FINAL_STEP

The **FINAL_STEP** is used to define the final step of the pipeline it can be set to "QC" if you want to stop after running the QC step, "ANALYSIS" to stop after the sample analysis step and "INTEGRATION" to stop after the integration (be careful if you have only one sample you can't run the integration part).

#### proseg threads

Proseg is allowing parallelisation, the **threads** parameter allows you to chose the amount of threads used by the proseg process.

#### MECR_PATH

Within the r_analysis_param section you must fill the **MECR_PATH** parameter. If you have a list of exclusive gene markers csv, first you need to copy it in the `01_Reference` directory then you must set **MECR_PATH** to `01_Reference/<MECR.tsv>` with &lt;MECR.csv&gt; to your exclusive gene marker file. If you don't have a exclusive pairs set **MECR_PATH** to `NULL` (without "").

##### Exclusive gene file

If you have an &lt;MECR.tsv&gt; file it need to follow this architecture :

Cell_type | Gene_list
----------|-----------
cell_type_1 | gene_1,gene_2,gene_x
cell_type_2 | gene_1,gene_2,gene_x

#### Other parameters

The remaining parameters are standard transcriptomic analysis parameters 

### Expert config

This config file contain more complex parameters that should be modified only if you have some experience with single cell RNA sequencing. 
The parameters that you might want to modify the most are the **DIM_USE** and **DIM_USE_INTEGRATED**, their default is 15 since spatial transcriptomic is using gene panels (~500 genes), you might want to increase it if you work on bigger gene panels 

## Snakefile 

:warning: :warning: :warning:
The snakefile should not be modified
:warning: :warning: :warning:

To launch the dry-run of the pipeline launch the following command.

```
snakemake -ps 04_Workflow/snakefile_proseg_pipeline.py \
    --use-singularity --singularity-args "-B /<volume>:/<volume>" \
    --cores 40 --dry-run
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

### Flowchart

![](02_Docs/dag_integrate.png)

### Xenium nuclear segmentation

This step take the original Xenium segmentation and resegment the data using only the nucleus. The output is a standard Xenium ranger output and is in the `01_Nucleus_segmentation` folder.

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

### Proseg Xenium Ranger

The Proseg Xenium Ranger step takes the Proseg output and create a Xenium analysis output, this is convenient for Xenium Explorer analysis and for later analysis.

```
04_Proseg_Xenium/
├── 04_Proseg_Xenium.mri.tgz
├── _cmdline
├── _filelist
├── _finalstate
├── _invocation
├── _jobmode
├── _log
├── _mrosource
├── outs
│    ├── analysis
│    ├── analysis_summary.html
│    ├── analysis.zarr.zip
│    ├── cell_boundaries.csv.gz
│    ├── cell_boundaries.parquet
│    ├── cell_feature_matrix
│    ├── cell_feature_matrix.h5
│    ├── cell_feature_matrix.zarr.zip
│    ├── cells.csv.gz
│    ├── cells.parquet
│    ├── cells.zarr.zip
│    ├── experiment.xenium
│    ├── gene_panel.json
│    ├── metrics_summary.csv
│    ├── morphology_focus
│    ├── morphology.ome.tif
│    ├── nucleus_boundaries.csv.gz
│    ├── nucleus_boundaries.parquet
│    ├── transcripts.parquet
│    └── transcripts.zarr.zip
├── _perf
├── _sitecheck
├── _tags
├── _timestamp
├── _uuid
├── _vdrkill
├── _versions
└── XENIUM_RANGER_CS
```
### R QC

This step performs the QC based on the minimum number of genes detected and the maximum size of a cell. The ouput generated consist in a seurat V5 object pre and post QC, an html report and the metadatas.

```
05_Proseg_QC/
├── metaData_no_QC.csv
├── proseg_seurat_no_QC.rds
├── proseg_seurat_QC.rds
└── Report_proseg_QC.html
```

### R analysis

This step performs the basic analysis (Normalizing, Scaling, PCA, UMAP, Neighbors, Clusters). The output generated are a seurat V5 object, an html report, a metadata csv, some metrics as csv and finally one Xenium Explorer compatible file per cluster resolution. This file allows you to visualize the seurat cluster in Xenium Explorer.

```
06_Proseg_analysis/
├── metaData.csv
├── metrics_summary.csv
├── proseg_seurat_analysis.rds
├── Report_proseg_analysis.html
├── XeniumExplorerClusterRes0.3.csv
├── XeniumExplorerClusterRes0.4.csv
└── XeniumExplorerClusterRes0.5.csv
```

### R merge and integration

This step merge of all the samples, then integrate them together using the chosen integration method. The output generated are a merged seurat V5 object, an integrated seurat V5 object, an html report, a metadata csv, all samples metrics concatenated and finally one Xenium Explorer compatible file per sample and cluster resolution. This file allows you to visualize the seurat cluster in Xenium Explorer.

```
07_Merge_and_integrate/
├── all_run_metrics_summary.csv
├── metaData.csv
├── proseg_rpca_integrated.rds
├── proseg_seurat_merged.rds
├── Report_proseg_integration.html
└── Xenium_Explorer_cluster
    ├── SAMPLE_1
    │    ├── ClusterRes0.2.csv
    │    └── ClusterRes0.3.csv
    ├── SAMPLE_2
    │    ├── ClusterRes0.2.csv
    │    └── ClusterRes0.3.csv
    └── SAMPLE_X
         ├── ClusterRes0.2.csv
         └── ClusterRes0.3.csv
```
