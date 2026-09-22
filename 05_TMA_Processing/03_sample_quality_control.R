
## @knitr sample_quality_control

cat("\n\n")
cat(sprintf("%s **SAMPLE QUALITY CONTROL**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This script evaluates the overall quality of samples on the TMA. Low-quality 
cells (e.g., low feature counts) will be filtered out later, but an initial quality 
filter is first applied at the sample level.

For example, if most cells within a sample are classified as low quality, the 
entire sample is removed from the analysis. This prevents biased downstream 
neighborhood analyses that could result from retaining samples with too few 
high-quality cells and excessive gaps in the tissue.")
cat("\n\n")

################################################################################
#'                             1. RETRIEVE SAMPLE QUALITY
#'                             
#' Iterate over each sample to compute used for sample-level quality assessment.
#'                             
################################################################################

# Initialize QC data for each sample
QC_data = data.frame(sample_id = c(),
                     total_area = c(),
                     conserved_cells_percentage = c())

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Tag cells with their corresponding quality
  xenium_data[[TMA]] = pastis_quality_control(xenium_data[[TMA]],
                                              minimum_features = CELL_MINIMUM_FEATURES,
                                              maximum_area = CELL_MAXIMUM_AREA)
  
  # Retrieve QC data
  temp_qc_data = data.frame("sample_id" = xenium_data[[TMA]]$SAMPLE_ID,
                            "cell_area" = xenium_data[[TMA]]$cell_area,
                            "quality_tag" = xenium_data[[TMA]]$QC_result,
                            "transcript_count" = xenium_data[[TMA]]$transcript_count)
  
  # Iterate over samples in the data
  for (sample in unique(temp_qc_data$sample_id)){
    
    # Retrieve current sample data
    sample_data = temp_qc_data[temp_qc_data$sample_id == sample, ]
    
    ###############
    ## CELL AREA ##
    ###############
    
    # Retrieve total area covered by cells of current sample
    temp_area = sum(sample_data$cell_area)
    
    #################
    ## QUALITY TAG ##
    #################
    
    # Retrieve proportion of conserved cells
    temp_conserved_cells_percentage = (nrow(sample_data[sample_data$quality_tag == "conserved", ]) / nrow(sample_data)) * 100
    
    #####################################
    ## MEAN CELL TRANSCRIPT PER SAMPLE ##
    #####################################
    
    # Retrieve mean cell transcript count in current sample
    temp_transcript_density = sum(sample_data$transcript_count)/temp_area
    
    ##################
    ## FILL QC DATA ##
    ##################
    
    # Initialize new row
    temp_row = data.frame(sample_id = sample,
                          total_area = temp_area,
                          conserved_cells_percentage = temp_conserved_cells_percentage,
                          transcript_density = temp_transcript_density)
    
    # Update global QC data
    QC_data = rbind(QC_data, temp_row)
  }
}

################################################################################
#'                                2. PLOT SAMPLE QC
################################################################################

cat("\n\n")
cat("The following sections are used for sample quality assessment:

  - QC table: Table summarizing sample-level statistics used for quality evaluation.
  - Total area: Distribution of sample areas; the red line indicates the filtering threshold.
  - Conserved cell percentage: Distribution of conserved cell percentages per sample; the red line indicates the filtering threshold.
  - Transcript density: Distribution of transcript density per sample; the red line indicates the filtering threshold.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s QC table", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(htmltools::tagList(DT::datatable(QC_data)))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Total area", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(QC_data, aes(x = total_area)) +
        geom_density(alpha = .2, fill="#FF6655") +
        geom_vline(xintercept = SAMPLE_MINIMUM_AREA, 
                   color = "red") +
        labs(x = "Total Area",
             y = "Density") +
        theme_minimal())
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Conserved cell percentage", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(QC_data, aes(x = conserved_cells_percentage)) +
        geom_density(alpha = .2, fill="#FF6655") +
        geom_vline(xintercept = SAMPLE_MINIMUM_CONSERVED_CELLS_PERCENTAGE,
                   color = "red") +
        labs(x = "Conserved cell percentage",
             y = "Density") +
        theme_minimal())
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Transcript density", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(QC_data, aes(x = transcript_density)) +
        geom_density(alpha = .2, fill="#FF6655") +
        geom_vline(xintercept = SAMPLE_MINIMUM_TRANSCRIPT_DENSITY, 
                   color = "red") +
        labs(x = "Transcript density",
             y = "Density") +
        theme_minimal())
cat("\n\n")

################################################################################
#'                        3. UPDATE SAMPLE QUALITY METADATA
################################################################################

cat("\n\n")
cat(sprintf("%s **TMA sample quality view**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step adds metadata to each TMA Seurat object indicating whether samples 
passed or failed the sample-level quality control.

After updating the Seurat objects, the sample-level quality control is visualized for each TMA:

  - Sample quality: Final quality label for each sample on the TMA.
  - Low-quality cells: Visualization of low-quality cells within each sample, often illustrating why a sample was classified as low quality.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Initialize sample quality meta data
  qc_metadata = data.frame(sample_id = xenium_data[[TMA]]$SAMPLE_ID,
                           sample_tissue_quality = "excluded")
  
  # Iterate over samples in the data
  for (sample in unique(qc_metadata$sample_id)){
    
    # Get current sample QC values
    sample_area = QC_data[QC_data$sample_id == sample, ]$total_area
    conserved_cells_percentage = QC_data[QC_data$sample_id == sample, ]$conserved_cells_percentage
    transcript_density = QC_data[QC_data$sample_id == sample, ]$transcript_density
    
    # Check sample quality thresholds
    if ((sample_area > SAMPLE_MINIMUM_AREA) & 
        (conserved_cells_percentage > SAMPLE_MINIMUM_CONSERVED_CELLS_PERCENTAGE) &
        (transcript_density > SAMPLE_MINIMUM_TRANSCRIPT_DENSITY)){
      
      # Set current sample as conserved (high enough quality)
      qc_metadata[qc_metadata$sample_id == sample, ]$sample_tissue_quality = "conserved"
      
    }
  }
  
  # Add new meta data layer relative to sample QC in the seurat object
  xenium_data[[TMA]]$sample_tissue_quality = qc_metadata$sample_tissue_quality
  
  cat("\n\n")
  cat(sprintf("%s Sample quality", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(ImageDimPlot(xenium_data[[TMA]], 
                     group.by = "sample_tissue_quality"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Low quality cells", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(ImageDimPlot(subset(xenium_data[[TMA]], subset = QC_result != "conserved"),
                     group.by = "QC_result", 
                     boundaries = "segmentation", 
                     border.size = 0.05, 
                     cols = QUALITY_CONTROL_COLOR_PALETTE, 
                     crop = T))
  cat("\n\n")
  
}
