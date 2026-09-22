# **03 - PROSEG HYPERPARAMETER SEARCH**

## **Overview and objectives**

The aim of this analysis step is to identify the set of Proseg hyperparameters, processed in the previous step (02_Cell_Segmentation), that yields the highest-quality cell segmentation. Segmentation quality is evaluated using several metrics, which will be detailed below.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Retrieve mutually exclusive gene pairs

Retrieve the mutually exclusive genes identified in step 01_Find_Mutually_Exclusive_Gene_Pairs. These genes are required to compute the MECR metric, which is the primary measure used here to assess cell segmentation quality.

From this set of genes, all heterotypic pairs (i.e., pairs between different cell types) are generated. These pairs correspond to the mutually exclusive gene pairs used to evaluate data quality.

### 2. Retrieve paths to Proseg outputs

Retrieve the paths to the segmentation outputs from each Proseg run generated during the hyperparameter grid search in the 02_Cell_Segmentation step. A table containing the run metadata is also initialized at this step.

### 3. Retrieve Proseg and initial Xenium segmentations statistics

Retrieve statistics for both the Proseg-generated segmentations and the initial Xenium segmentation. These metrics will be used to evaluate the effect of hyperparameter choices on segmentation quality.

For each segmentation run, the following metrics are retrieved:

  - Number of cells in the data (**nb_cells**): Since Proseg works by refining existing cells rather than generating new ones, this value should remain close to the number of cells in the initial Xenium segmentation.
  - Number of conserved cells (**nb_conserved**): Conserved cells are those that pass quality control (QC) criteria based on feature counts (CELL_MINIMUM_FEATURES) and cell size (CELL_MAXIMUM_AREA). Because Proseg performs additional filtering on the transcripts, it typically reduces the number of conserved cells obtained during QC. Ideally, the selected hyperparameters should preserve a conserved cell count as close as possible to that of the initial Xenium segmentation. Since only conserved cells are retained for downstream analyses, excessive cell loss may create gaps in the tissue and negatively impact subsequent steps, such as neighborhood analysis.
  - Number of low feature cells (**nb_low_feature**): Low feature cells are cells that fail the feature count filter but pass the cell size filter. Proseg tends to increase the number of low feature cells during segmentation refinement, and the selected hyperparameters should minimize their occurrence.
  - Number of high area cells (**nb_high_area**): High area cells are cells that pass the feature count filter but fail the cell size filter.
  - Number of low feature and high area cells (**nb_low_features_high_area**): Cells that fail both quality thresholds.
  - Mean mutually exclusive co-expression rate (**mean_MECR**): The mutually exclusive co-expression rate is computed per gene pair, and the mean value is then calculated across the dataset. MECR is defined as the fraction of cells expressing both markers divided by the fraction of cells expressing at least one of the markers (Satija et al. Biorxiv 2024).
  - Median mutually exclusive co-expression rate (**median_MECR**): The median value of the mutually exclusive co-expression rate across all gene pairs, computed in the same way as above.

## **02_plot_MECR_scores.R**

Display the MECR score distribution for each hyperparameter set used during the 02_Cell_Segmentation step. This enables visualization of the improvements provided by Proseg compared to the initial Xenium segmentation. It also provides an overview of how changing Proseg hyperparameters influence MECR values across the dataset.

## **03_plot_correlations.R**

The goal is to plot the statistics retrieved in the previous step against each other. This kind of visualization enables rapid identification of the Proseg hyperparameter set that yields the highest-quality segmentation. For each figure, observations are colored according to one of the evaluated Proseg parameters.

The selected hyperparameter set is defined as the one that produces a segmentation with:

  - Low MECR values
  - Few low-feature cells
  - A number of conserved cells similar to that of the initial Xenium segmentation

In these Proseg runs, it is commonly observed that as MECR decreases, the number of low-feature cells increases. This trade-off makes Proseg parameter selection challenging. The objective is therefore to identify a parameter set that achieves sufficiently low MECR while maintaining a low number of low-feature cells (or, equivalently, preserving a relatively high number of conserved cells).

## **90_finalizeSession.R**

This script aims to save the table containing all quality statistics computed for each Proseg run. It also saves the MECR scores calculated for each run.
