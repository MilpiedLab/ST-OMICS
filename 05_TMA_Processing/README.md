# **05 - TMA PROCESSING**

## **Overview and objectives**

This script splits the TMA into individual samples. Currently, the Seurat objects do not contain metadata linking each cell to its corresponding sample on the TMA.

The goal is therefore to identify the spatial location of each sample on the TMA and assign cells accordingly.

Sample detection is performed using a nearest-neighbor graph. Cells are first connected to their spatial nearest neighbors, and the resulting connected components are then classified as either individual samples or background based on their size (i.e., number of cells per component). See the split_xenium_tma_samples function for more details.

Docker use: "STOmics_Seurat5".

## **01_prepareData.R**

### 1. Load Seurat objects

This step loads the spatial transcriptomics data with Proseg-refined segmentation for each TMA under study. For each TMA, it uses the optimal Proseg parameters identified in the '03_Proseg_Hyperparameter_Search' step.

## **02_split_tma.R**

### 1. Annotate samples on TMA

This step aims to identify and separate individual samples on each TMA.

For each TMA, the following visualizations are generated:

  - Global TMA view: Displays the identified tissue sections on the TMA. Cells labeled as excluded are considered not to belong to any sample, based on the parameters used for TMA splitting.
  - Raw TMA grid: Tabular representation of sample positions on the TMA.
  - Annotated TMA: Displays the segmented samples on the TMA along with their corresponding sample annotations.

### 2. Rename cells based on sample ID

Now that all samples are labeled, rename all cells accordingly to ensure unique identifiers when merging samples into a single Seurat object in downstream analyses.

### 3. Display sample statistics

This step aims to visualize sample-specific statistics and identify potentially low-quality samples.

  - Feature count: Ridge plot showing the distribution of feature counts per cell across samples.
  - Mean gene count per sample: For each gene, the mean expression across cells is computed within each sample:
    * Boxplot: Distribution of mean gene counts across samples.
    * Ridge plot: Distribution of mean gene counts across samples.
    * Mean vs. standard deviation: Relationship between mean gene count and its variability.
    * Complete statistics: Table containing all values used to generate the above plots.
  - Sample similarity based on mean gene expression: Pearson correlation is computed between samples using their vectors of mean gene expression per gene. This provides a measure of transcriptional similarity between samples.

## **03_sample_quality_control.R**

This step evaluates the overall quality of samples on the TMA. Low-quality cells (e.g., low feature counts) will be filtered out later, but an initial quality filter is first applied at the sample level.

For example, if most cells within a sample are classified as low quality, the entire sample is removed from the analysis. This prevents biased downstream neighborhood analyses that could result from retaining samples with too few high-quality cells and excessive gaps in the tissue.

### 1. Retrieve sample quality

Iterate over each sample to compute used for sample-level quality assessment.

### 2. Plot sample QC

The following sections are used for sample quality assessment:

  - QC table: Table summarizing sample-level statistics used for quality evaluation.
  - Total area: Distribution of sample areas; the red line indicates the filtering threshold.
  - Conserved cell percentage: Distribution of conserved cell percentages per sample; the red line indicates the filtering threshold.
  - Transcript density: Distribution of transcript density per sample; the red line indicates the filtering threshold.

### 3. Update sample quality metadata

This step adds metadata to each TMA Seurat object indicating whether samples passed or failed the sample-level quality control.

After updating the Seurat objects, the sample-level quality control is visualized for each TMA:

  - Sample quality: Final quality label for each sample on the TMA.
  - Low-quality cells: Visualization of low-quality cells within each sample, often illustrating why a sample was classified as low quality.

## **90_finalizeSession.R**

Save the Seurat objects corresponding to the optimal Proseg segmentation for each TMA, with each cell now assigned to its corresponding sample on the TMA. These objects also include metadata from both sample-level quality assessment and cell-level quality control. At this stage, no filtering has been applied to either cells or samples. All filtering steps will be performed downstream on these Seurat objects.
