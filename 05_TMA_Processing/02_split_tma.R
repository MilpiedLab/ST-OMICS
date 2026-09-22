
## @knitr split_TMA

################################################################################
#'                           1. ANNOTATE SAMPLES ON TMA
################################################################################

cat("\n\n")
cat(sprintf("%s **TMA splitting**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to identify and separate individual samples on each TMA.

For each TMA, the following visualizations are generated:

  - Global TMA view: Displays the identified tissue sections on the TMA. Cells labeled as excluded are considered not to belong to any sample, based on the parameters used for TMA splitting.
  - Raw TMA grid: Tabular representation of sample positions on the TMA.
  - Annotated TMA: Displays the segmented samples on the TMA along with their corresponding sample annotations.")
cat("\n\n")

# Store TMA grids for later
xenium_grids = list()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s ", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  ####################
  ## 1.1. SPLIT TMA ##
  ####################
  
  # Split current TMA
  split_results = split_xenium_tma_samples(seurat_object = xenium_data[[TMA]],
                                           radius_threshold = SPLIT_TMA_PARAMS[[TMA]]$RADIUS_THRESHOLD,
                                           neighbors_number = SPLIT_TMA_PARAMS[[TMA]]$NEIGHBORS_NUMBER,
                                           sample_minimum_size = SPLIT_TMA_PARAMS[[TMA]]$SAMPLE_MINIMUM_SIZE,
                                           number_sample_rows = SPLIT_TMA_PARAMS[[TMA]]$NUMBER_SAMPLE_ROWS,
                                           number_sample_columns = SPLIT_TMA_PARAMS[[TMA]]$NUMBER_SAMPLE_COLUMNS)
  
  cat("\n\n")
  cat(sprintf("%s Global TMA view ", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(ImageDimPlot(split_results$splitted_seurat, 
                     group.by = "filtered_sample_identity"))
  cat("\n\n")
  
  cat("\n\n")
  cat(sprintf("%s Raw TMA grid ", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(htmltools::tagList(DT::datatable(split_results$grid)))
  cat("\n\n")
  
  # Update Seurat object with the new meta data corresponding to punches assignment
  xenium_data[[TMA]] = split_results$splitted_seurat
  
  # Remove cells that are tagged as excluded
  xenium_data[[TMA]] = subset(xenium_data[[TMA]], 
                              subset = filtered_sample_identity != "excluded")
  
  ###############################
  ## 1.2. RETRIEVE SAMPLES IDS ##
  ###############################
  
  # FACULTATIVE: Adjust samples orientation on TMA
  split_results$grid = split_results$grid[, ncol(split_results$grid):1]
  split_results$grid = split_results$grid[nrow(split_results$grid):1, ]
  
  # Retrieve temporary sample id association
  temp_punches_ids = xenium_data[[TMA]]$filtered_sample_identity
  temp_patients_ids = xenium_data[[TMA]]$filtered_sample_identity
  
  # Replace temporary id by the actual sample ID
  for (i in 1:nrow(TMA_GRID[[TMA]])){
    for (j in 1:ncol(TMA_GRID[[TMA]])){
      if (!(is.na(TMA_GRID[[TMA]][i, j]))){
        
        # Retrieve punch (tissue) id to assign
        temp_SAMPLE_ID = paste0(TMA_GRID[[TMA]][i, j], "_", TMA)
        temp_PATIENT_ID = paste0(TMA_GRID[[TMA]][i, j])
        
        # Replace temporary punch id by real SAMPLE ID and PATIENT ID
        temp_punches_ids[temp_punches_ids == split_results$grid[i, j]] = temp_SAMPLE_ID
        temp_patients_ids[temp_patients_ids == split_results$grid[i, j]] = temp_PATIENT_ID
        
      }
    }
  }
  
  # Store grid to plot later
  xenium_grids[[TMA]] = split_results$grid
  
  # Update Seurat object with the new meta data corresponding to punches assignment
  xenium_data[[TMA]]$SAMPLE_ID = temp_punches_ids
  
  # Add meta data corresponding to patient ID
  xenium_data[[TMA]]$PATIENT_ID = temp_patients_ids
  
  # Activate new meta data layer
  xenium_data[[TMA]] = activateSeuratAssay(xenium_data[[TMA]], 
                                           ident = "SAMPLE_ID")
  
  cat("\n\n")
  cat(sprintf("%s Annotated TMA ", strrep("#", 3)))
  cat("\n\n")
  
  cat("\n\n")
  print(ImageDimPlot(xenium_data[[TMA]], 
                     group.by = "SAMPLE_ID"))
  cat("\n\n")
}

################################################################################
#'                       2. RENAME CELLS BASED ON SAMPLE ID
#'                       
#' Now that all samples are labeled, rename all cells accordingly to ensure 
#' unique identifiers when merging samples into a single Seurat object in 
#' downstream analyses.  
#'                      
################################################################################

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Retrieve current cell names and tissue (punch) id
  temp_TMA_cell_data = data.frame(cell = Cells(xenium_data[[TMA]]),
                                  SAMPLE_ID = xenium_data[[TMA]]$SAMPLE_ID)
  
  # Rename cells in temp df
  temp_TMA_cell_data$cell = apply(temp_TMA_cell_data, 1, FUN = function(x) {
    x[1] = paste0(x[2], "_", x[1])
  })
  
  # Update Seurat object
  xenium_data[[TMA]] = RenameCells(xenium_data[[TMA]], new.names = temp_TMA_cell_data$cell)
  
}

################################################################################
#'                           3. DISPLAY SAMPLE STATISTICS
################################################################################

cat("\n\n")
cat(sprintf("%s **SAMPLE STATISTICS**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This step aims to visualize sample-specific statistics and identify potentially low-quality samples.

  - Feature count: Ridge plot showing the distribution of feature counts per cell across samples.
  - Mean gene count per sample: For each gene, the mean expression across cells is computed within each sample:
    * Boxplot: Distribution of mean gene counts across samples.
    * Ridge plot: Distribution of mean gene counts across samples.
    * Mean vs. standard deviation: Relationship between mean gene count and its variability.
    * Complete statistics: Table containing all values used to generate the above plots.
  - Sample similarity based on mean gene expression: Pearson correlation is computed between samples using their vectors of mean gene expression per gene. This provides a measure of transcriptional similarity between samples.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Feature count", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Init QC data
QC_data = data.frame()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Retrieve qc data for current TMA
  temp_qc_data = data.frame(sample_ID = xenium_data[[TMA]]$SAMPLE_ID,
                            n_feature = xenium_data[[TMA]]$nFeature_Xenium)
  
  # Add to global QC data
  QC_data = rbind(QC_data, temp_qc_data)
  
}

# Plot distribution of feature number per cell by sample
cat("\n\n")
print(ggplot(QC_data, aes(x = n_feature, y = sample_ID, fill = sample_ID)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(0, 150)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
        geom_text(data=QC_data %>% group_by(sample_ID) %>% 
                    summarise(n_feature=median(n_feature)),
                  aes(label=sprintf("%1.1f", n_feature)), 
                  position=position_nudge(x = 10, y=0.8), colour="black", size=3.5))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Mean gene count per sample", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Init gene expression data for all TMAs
samples_gene_expressions = data.frame()

# Fetch and prepare count data for all cells and all genes
for (TMA in TMA_COVERED){
  
  # Retrieve gene panel
  gene_panel = rownames(xenium_data[[TMA]])
  
  # Retrieve gene expressions
  temp_samples_gene_expressions = FetchData(xenium_data[[TMA]], 
                                       vars = gene_panel,
                                       layer = 'counts')
  
  # Add sample ID
  temp_samples_gene_expressions$SAMPLE_ID = xenium_data[[TMA]]$SAMPLE_ID
  
  # Add to global gene expression df
  samples_gene_expressions = rbind(samples_gene_expressions, temp_samples_gene_expressions)
  
}

# For each sample, compute mean count value and stdev of each gene in the panel
samples_statistics_list = list()
for (sample in unique(samples_gene_expressions$SAMPLE_ID)){
  
  # Subset current sample gene expressions
  current_sample = samples_gene_expressions[samples_gene_expressions$SAMPLE_ID == sample, ]
  current_sample$SAMPLE_ID = NULL
  
  # Retrieve mean counts and std dev for each gene of current sample 
  sample_stats = apply(X = current_sample, MARGIN = 2, FUN = function(x){
    
    return(c("mean" = mean(x), "stdev" = sd(x)))
    
  })
  
  # Update statistics list
  sample_stats = as.data.frame(t(sample_stats))
  sample_stats$SAMPLE_ID = sample
  samples_statistics_list[[sample]] = sample_stats
}

# Merge statistics of all samples
complete_sample_statistics = do.call(rbind, samples_statistics_list)

cat("\n\n")
cat(sprintf("%s Boxplot", strrep("#", 3)))
cat("\n\n")

# Retrieve outliers for boxplot
outliers = complete_sample_statistics %>%
  group_by(SAMPLE_ID) %>%
  mutate(
    Q1 = quantile(mean, 0.25),
    Q3 = quantile(mean, 0.75),
    IQR = Q3 - Q1,
    is_outlier = mean < Q1 - 1.5 * IQR | mean > Q3 + 1.5 * IQR
  ) %>%
  filter(is_outlier)

# Display boxplots of gene mean expression by sample
cat("\n\n")
ggplot(complete_sample_statistics, aes(x = SAMPLE_ID, y = mean, fill = SAMPLE_ID)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(data = outliers, aes(x = SAMPLE_ID, y = mean), color = "black", size = 2) +
  scale_fill_manual(values = COLOR_PALETTE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "", x = "Segmentation", y = "Mean gene count per cell")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Ridge plot", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(ggplot(complete_sample_statistics, aes(x = mean, y = SAMPLE_ID, fill = SAMPLE_ID)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        scale_fill_manual(values = COLOR_PALETTE) +
        scale_x_continuous(limits = c(0, 0.25)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
        geom_text(data=complete_sample_statistics %>% group_by(SAMPLE_ID) %>% 
                    summarise(mean=median(mean)),
                  aes(label=sprintf("%1.1f", mean)), 
                  position=position_nudge(x = 10, y=0.8), colour="black", size=3.5))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Mean vs. standard deviation", strrep("#", 3)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Complete", strrep("#", 4)))
cat("\n\n")

ggplot(complete_sample_statistics, aes(x = mean, y = stdev, color = SAMPLE_ID)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values=COLOR_PALETTE) +
  labs(
    title = "Gene counts (mean vs sd)",
    x = "Mean",
    y = "Standard deviation"
  ) +
  theme_minimal()

cat("\n\n")
cat(sprintf("%s Zoomed", strrep("#", 4)))
cat("\n\n")

ggplot(complete_sample_statistics, aes(x = mean, y = stdev, color = SAMPLE_ID)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values=COLOR_PALETTE) +
  labs(
    title = "Gene counts (mean vs sd)",
    x = "Mean",
    y = "Standard deviation"
  ) +
  xlim(c(0, 1)) +
  ylim(c(0, 1)) +
  theme_minimal()

cat("\n\n")
cat(sprintf("%s Complete statistics", strrep("#", 3)))
cat("\n\n")

cat("\n\n")
print(htmltools::tagList(DT::datatable(complete_sample_statistics)))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Sample similarity based on mean gene expression", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Initialize correlation matrix
sample_number = length(unique(complete_sample_statistics$SAMPLE_ID))
correlation_matrix = matrix(NA, nrow = sample_number, ncol = sample_number)
rownames(correlation_matrix) = unique(complete_sample_statistics$SAMPLE_ID)
colnames(correlation_matrix) = unique(complete_sample_statistics$SAMPLE_ID)
correlation_matrix = as.data.frame(correlation_matrix)

# Compute correlation values for each sample pair one by one
for (row in rownames(correlation_matrix)){
  for (col in colnames(correlation_matrix)){
    
    correlation_matrix[row, col] = cor(complete_sample_statistics[complete_sample_statistics$SAMPLE_ID == row, ]$mean,
                                       complete_sample_statistics[complete_sample_statistics$SAMPLE_ID == col, ]$mean,
                                       method = "pearson")
    
  }
}

# Plot heatmap
cat("\n\n")
print(pheatmap::pheatmap(correlation_matrix))
cat("\n\n")
