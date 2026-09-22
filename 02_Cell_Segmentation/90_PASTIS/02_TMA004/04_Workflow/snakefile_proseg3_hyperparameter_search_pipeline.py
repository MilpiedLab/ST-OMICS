import os
import re

configfile: "90_PASTIS/02_TMA004/04_Workflow/config_proseg3_hyperparameter_search/mandatory_config.yaml"
configfile: "90_PASTIS/02_TMA004/04_Workflow/config_proseg3_hyperparameter_search/basic_config.yaml"
configfile: "90_PASTIS/02_TMA004/04_Workflow/config_proseg3_hyperparameter_search/expert_config.yaml"

print(config)

# Import config parameters

PATH_RAW_DATA_XENIUM = config["PATH_RAW_DATA_XENIUM"]
PATH_ANALYSIS_OUTPUT = config["PATH_ANALYSIS_OUTPUT"]

print(PATH_ANALYSIS_OUTPUT)

samples_raw = config["samples"]
samples = samples_raw.split(",")

# Import ranges of parameters to explore
nuc_assignment = config["proseg"]["nuc_assignment"]
diffusion_probability = config["proseg"]["diffusion_probability"]
voxel_layers = config["proseg"]["voxel_layers"]

###################################################
#           Define default rule : all             #
###################################################

rule all:
  input:
    expand(f"{PATH_ANALYSIS_OUTPUT}/{{smp}}/proseg3_hyperparameter_search/voxel_layers{{voxel_layers}}/nuc_assignment{{nuc_assignment}}/diffusion_probability{{diffusion_probability}}/XR_Proseg_done.txt",
    PATH_ANALYSIS_OUTPUT = PATH_ANALYSIS_OUTPUT,
    smp=samples,
    nuc_assignment=nuc_assignment,
    diffusion_probability=diffusion_probability,
    voxel_layers=voxel_layers)

#####################################################
#           Run proseg on xenium dataset            #
#####################################################

if config["proseg"]["datatype"] == "xenium":
    # Using xenium dataset
    rule proseg:
        input: transcripts = lambda wildcards: f"{PATH_RAW_DATA_XENIUM}/{wildcards.smp}/transcripts.parquet"
        output: "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/transcript-metadata.csv.gz",
                "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/cell-polygons.geojson.gz",
                directory("{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/proseg-output.zarr")
        container: config["docker"]["proseg"]
        params : ncomponents = config["proseg"]["ncomponents"],
                 nuc_assignment = lambda wc: wc.nuc_assignment,
                 diffusion_probability = lambda wc: wc.diffusion_probability,
                 voxel_layers = lambda wc: wc.voxel_layers,
                 rna_diffusion = config["proseg"]["rna_diffusion"],
                 min_qv = config["proseg"]["min_qv"]
        threads: config["proseg"]["threads"]
        shell:
             """
             proseg --xenium {input.transcripts} \
             --nthreads {threads} \
             --ncomponents {params.ncomponents} \
             --voxel-layers {params.voxel_layers} \
             --nuclear-reassignment-prob {params.nuc_assignment} \
             --min-qv {params.min_qv} \
             --diffusion-probability {params.diffusion_probability} \
             --output-spatialdata {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/proseg-output.zarr \
             --output-expected-counts {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/expected-counts.csv.gz \
             --output-cell-metadata {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-metadata.csv.gz \
             --output-transcript-metadata {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/transcript-metadata.csv.gz \
             --output-cell-polygons {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-polygons.geojson.gz \
             --output-cell-polygon-layers {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-polygons-layers.geojson.gz \
             --output-union-cell-polygons {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/union-cell-polygons.geojson.gz
             """
elif config["proseg"]["datatype"] == "custom":
    # Using custom dataset
    rule proseg:
        input: transcripts = lambda wildcards: f"{PATH_RAW_DATA_XENIUM}/{wildcards.smp}/transcripts.parquet"
        output: "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/transcript-metadata.csv.gz",
                "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/cell-polygons.geojson.gz",
                directory("{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/proseg-output.zarr")
        container: config["docker"]["proseg"]
        params : x = config["proseg"]["x"],
                 y = config["proseg"]["y"],
                 z = config["proseg"]["z"],
                 transcript_id = config["proseg"]["transcript_id"],
                 gene_col = config["proseg"]["gene_col"],
                 cell_col = config["proseg"]["cell_col"],
                 nuc_name = config["proseg"]["nuc_name"],
                 fov_col = config["proseg"]["fov_col"],
                 cell_id = config["proseg"]["cell_id"],
                 cell_id_unassigned = config["proseg"]["cell_id_unassigned"],
                 qv = config["proseg"]["qv"],
                 ncomponents = config["proseg"]["ncomponents"],
                 voxel_layers = lambda wc: wc.voxel_layers,
                 nuc_assignment = lambda wc: wc.nuc_assignment,
                 rna_diffusion = config["proseg"]["rna_diffusion"],
                 min_qv = config["proseg"]["min_qv"],
                 diffusion_probability = lambda wc: wc.diffusion_probability
        threads: config["proseg"]["threads"]
        shell:
             """
             proseg --xenium {input.transcripts} \
             --nthreads {threads} \
             -x {params.x} -y {params.y} -z {params.z}\
             --transcript-id-column {params.transcript_id} \
             --gene-column {params.gene_col} \
             --compartment-column {params.cell_col} \
             --compartment-nuclear {params.nuc_name} \
             --fov-column {params.fov_col} \
             --cell-id-column {params.cell_id} \
             --cell-id-unassigned {params.cell_id_unassigned} \
             --qv-column {params.qv} \
             --ncomponents {params.ncomponents} \
             --voxel-layers {params.voxel_layers} \
             --nuclear-reassignment-prob {params.nuc_assignment} \
             --min-qv {params.min_qv} \
             --diffusion-probability {params.diffusion_probability} \
             --output-spatialdata {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/proseg-output.zarr \
             --output-expected-counts {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/expected-counts.csv.gz \
             --output-cell-metadata {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-metadata.csv.gz \
             --output-transcript-metadata {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/transcript-metadata.csv.gz \
             --output-cell-polygons {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-polygons.geojson.gz \
             --output-cell-polygon-layers {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/cell-polygons-layers.geojson.gz \
             --output-union-cell-polygons {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/02_Proseg_run/union-cell-polygons.geojson.gz
             """
             

#############################################
#           Run proseg to baysor            #
#############################################

rule proseg_to_baysor:
    input: zarr_output = "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/02_Proseg_run/proseg-output.zarr" 
    output: "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/03_Proseg_to_Baysor/baysor-transcript-metadata.csv",
            "{PATH_ANALYSIS_OUTPUT}/{smp}/proseg3_hyperparameter_search/voxel_layers{voxel_layers}/nuc_assignment{nuc_assignment}/diffusion_probability{diffusion_probability}/03_Proseg_to_Baysor/baysor-cell-polygons.geojson"
    container: config["docker"]["proseg"]
    threads: 10
    shell:
         """
         proseg-to-baysor {input.zarr_output} \
         --output-transcript-metadata {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/03_Proseg_to_Baysor/baysor-transcript-metadata.csv \
         --output-cell-polygons {wildcards.PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}/03_Proseg_to_Baysor/baysor-cell-polygons.geojson
         """


##########################################
#           Run Xenium ranger            #
##########################################

rule proseg_xenium_ranger:
    input: transcripts = f"{PATH_ANALYSIS_OUTPUT}/{{smp}}/proseg3_hyperparameter_search/voxel_layers{{voxel_layers}}/nuc_assignment{{nuc_assignment}}/diffusion_probability{{diffusion_probability}}/03_Proseg_to_Baysor/baysor-transcript-metadata.csv",
           polygons = f"{PATH_ANALYSIS_OUTPUT}/{{smp}}/proseg3_hyperparameter_search/voxel_layers{{voxel_layers}}/nuc_assignment{{nuc_assignment}}/diffusion_probability{{diffusion_probability}}/03_Proseg_to_Baysor/baysor-cell-polygons.geojson"
    output: f"{PATH_ANALYSIS_OUTPUT}/{{smp}}/proseg3_hyperparameter_search/voxel_layers{{voxel_layers}}/nuc_assignment{{nuc_assignment}}/diffusion_probability{{diffusion_probability}}/XR_Proseg_done.txt"
    container: config["docker"]["xenium_ranger"]
    params: xenium_ori = lambda wildcards: f"{PATH_RAW_DATA_XENIUM}/{wildcards.smp}/",
            transcript_path = "03_Proseg_to_Baysor/baysor-transcript-metadata.csv",
            polygons_path = "03_Proseg_to_Baysor/baysor-cell-polygons.geojson"
    threads: 10
    shell:
         """
         cd {PATH_ANALYSIS_OUTPUT}/{wildcards.smp}/proseg3_hyperparameter_search/voxel_layers{wildcards.voxel_layers}/nuc_assignment{wildcards.nuc_assignment}/diffusion_probability{wildcards.diffusion_probability}
         xeniumranger import-segmentation \
         --id 04_Proseg_Xenium \
         --xenium-bundle={params.xenium_ori} \
         --viz-polygons={params.polygons_path} \
         --transcript-assignment={params.transcript_path} \
         --units=microns
         touch XR_Proseg_done.txt
         """
