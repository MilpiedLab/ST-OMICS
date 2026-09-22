
## @knitr quality_control

################################################################################
#'                                1. QUALITY CONTROL
#' 
#' This step consists of performing quality control on cells from both the 
#' initial Xenium segmentation and the Proseg-refined segmentation, in order to 
#' highlight differences between the two conditions.
#' 
#' For now, quality control is based on two criteria: the number of features 
#' per cell and cell area.
#' 
#' Number of conserved cells (**nb_conserved**): Conserved cells are those that 
#' pass quality control (QC) criteria based on feature counts 
#' (CELL_MINIMUM_FEATURES) and cell size (CELL_MAXIMUM_AREA).
#' 
#' Number of low feature cells (**nb_low_feature**): Low feature cells are cells 
#' that fail the feature count filter but pass the cell size filter.
#' 
#' Number of high area cells (**nb_high_area**): High area cells are cells that 
#' pass the feature count filter but fail the cell size filter.
#' 
#' Number of low feature and high area cells (**nb_low_features_high_area**): 
#' Cells that fail both quality thresholds.
#' 
################################################################################

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Tag cells that should be removed based on filtering parameters for current condition
    xenium_data[[TMA]][[segmentation_method]] = pastis_quality_control(xenium_data[[TMA]][[segmentation_method]],
                                                                       minimum_features = CELL_MINIMUM_FEATURES,
                                                                       maximum_area = CELL_MAXIMUM_AREA)
    
  }
}

## @knitr plot_quality_control

################################################################################
#'                              2. PLOT QUALITY CONTROL
################################################################################

cat("\n\n")
cat(sprintf("%s **Quality control**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This section generates all figures related to quality control. For each TMA 
and each segmentation method, the following visualizations are produced:

  - Violin plots: Distributions of the main cell-level statistics.
  - Scatter plots: Pairwise comparisons of key cell statistics.
  - Cell size distributions: Ridge plots showing the distribution of cell area 
across the different cell quality categories defined above.
  - QC images: TMA views of the tissue sections, with cells colored according to 
their quality category.
  - Top gene distribution: Visualization of the most frequently detected genes in 
the dataset and their relative proportions within the library.

A 'Histograms' section is also included, providing a direct comparison of the 
number of cells in each quality category across the TMAs and the segmentation 
methods under study.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s ", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    cat("\n\n")
    cat(sprintf("%s %s ", strrep("#", 3), segmentation_method), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    #####################
    # 2.1. VIOLIN PLOTS #
    #####################
    
    nCount_vln = VlnPlot(xenium_data[[TMA]][[segmentation_method]], features = c("nCount_Xenium"), pt.size = 0, cols = QUALITY_CONTROL_COLOR_PALETTE[1]) +
      theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    nFeature_vln = VlnPlot(xenium_data[[TMA]][[segmentation_method]], features = c("nFeature_Xenium"), pt.size = 0, cols = QUALITY_CONTROL_COLOR_PALETTE[1]) +
      theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      geom_hline(yintercept = CELL_MINIMUM_FEATURES, linewidth = 1, linetype = "dashed")
    
    area_vln = VlnPlot(xenium_data[[TMA]][[segmentation_method]], features = c("cell_area"), pt.size = 0, cols = QUALITY_CONTROL_COLOR_PALETTE[1]) +
      theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      geom_hline(yintercept = CELL_MAXIMUM_AREA, linewidth = 1, linetype = "dashed")
    
    grouped_violin = plot_grid(nCount_vln, nFeature_vln, area_vln, ncol = 3)
    
    cat("\n\n")
    cat(sprintf("%s Violin plot ", strrep("#", 4)))
    cat("\n\n")
    
    cat("\n\n")
    print(grouped_violin)
    cat("\n\n")
    
    ######################
    # 2.2. SCATTER PLOTS #
    ######################
    
    cat("\n\n")
    cat(sprintf("%s Scatter plot ", strrep("#", 4)), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Count per Features ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    nCount_nFeature_scatter = ggplot(xenium_data[[TMA]][[segmentation_method]]@meta.data, 
                                     aes(x = nCount_Xenium, y = nFeature_Xenium, col = QC_result)) +
      geom_point() +
      theme_classic() +
      geom_hline(yintercept = CELL_MINIMUM_FEATURES, linewidth=1, linetype = "dashed") +
      scale_color_manual(values = QUALITY_CONTROL_COLOR_PALETTE) + 
      xlim(0, 8000) +
      ylim(0, 300)
    print(nCount_nFeature_scatter)
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Count per cell area ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    nCount_cell_area_scatter = ggplot(xenium_data[[TMA]][[segmentation_method]]@meta.data, 
                                      aes(x = nCount_Xenium, y = cell_area, col = QC_result)) +
      geom_point() +
      theme_classic() +
      geom_hline(yintercept = CELL_MAXIMUM_AREA, linewidth = 1, linetype = "dashed") +
      scale_color_manual(values = QUALITY_CONTROL_COLOR_PALETTE)
    print(nCount_cell_area_scatter)
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Features per cell area ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    nFeature_cell_area_scatter = ggplot(xenium_data[[TMA]][[segmentation_method]]@meta.data, 
                                        aes(x = nFeature_Xenium, y = cell_area, col = QC_result)) +
      geom_point() +
      theme_classic() +
      geom_hline(yintercept = CELL_MAXIMUM_AREA, linewidth = 1, linetype = "dashed") +
      geom_vline(xintercept = CELL_MINIMUM_FEATURES, linewidth = 1, linetype = "dashed") +
      scale_color_manual(values = QUALITY_CONTROL_COLOR_PALETTE)
    print(nFeature_cell_area_scatter)
    cat("\n\n")
    
    ###############################
    # 2.3. CELL SIZE DISTRIBUTION #
    ###############################
    
    cat("\n\n")
    cat(sprintf("%s Cell size distribution ", strrep("#", 4)), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    # Retrieve cell distributions by QC results
    QC_tag_data = data.frame(cell_area = xenium_data[[TMA]][[segmentation_method]]@meta.data$cell_area,
                             QC_tag = xenium_data[[TMA]][[segmentation_method]]@meta.data$QC_result)
    
    cat("\n\n")
    cat(sprintf("%s Total ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    print(ggplot(QC_tag_data, aes(x = cell_area, y = QC_tag, fill = QC_tag)) +
            geom_density_ridges() +
            theme_ridges() + 
            theme(legend.position = "none") +
            scale_x_continuous(limits = c(0, 2000)) +
            scale_fill_manual(values = QUALITY_CONTROL_COLOR_PALETTE) +
            stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
            geom_text(data=QC_tag_data %>% group_by(QC_tag) %>% 
                        summarise(cell_area=median(cell_area)),
                      aes(label=sprintf("%1.1f", cell_area)), 
                      position=position_nudge(x = 50, y=0.2), colour="black", size=3.5) +
            labs(title = sprintf("Cell area distribution by QC results (TMA = %s, segmentation = %s)", TMA, segmentation_method)))
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Zoomed ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    print(ggplot(QC_tag_data, aes(x = cell_area, y = QC_tag, fill = QC_tag)) +
            geom_density_ridges() +
            theme_ridges() + 
            theme(legend.position = "none") +
            scale_x_continuous(limits = c(0, 200)) +
            scale_fill_manual(values = QUALITY_CONTROL_COLOR_PALETTE) +
            stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
            geom_text(data=QC_tag_data %>% group_by(QC_tag) %>% 
                        summarise(cell_area=median(cell_area)),
                      aes(label=sprintf("%1.1f", cell_area)), 
                      position=position_nudge(x = 10, y=0.3), colour="black", size=3.5) +
            labs(title = sprintf("Cell area distribution by QC results (TMA = %s, segmentation = %s)", TMA, segmentation_method)))
    cat("\n\n")
    
    #################
    # 2.4. QC IMAGE #
    #################
    
    cat("\n\n")
    cat(sprintf("%s QC image ", strrep("#", 4)), " {.tabset .tabset-fade}\n\n")
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s All cells ", strrep("#", 5)))
    cat("\n\n")
    
    cat("\n\n")
    tissue_qc = ImageDimPlot(xenium_data[[TMA]][[segmentation_method]], 
                             group.by = "QC_result", 
                             boundaries = "segmentation", 
                             border.size = 0.05, 
                             cols = QUALITY_CONTROL_COLOR_PALETTE, 
                             crop = T)
    print(tissue_qc)
    cat("\n\n")
    
    cat("\n\n")
    cat(sprintf("%s Filtered cells ", strrep("#", 5)))
    cat("\n\n")
    
    # Retrieve cell tags that are present in the data
    Idents(xenium_data[[TMA]][[segmentation_method]]) = "QC_result"
    available_cell_tags = intersect(c("low_feature", 
                                      "high_area", 
                                      "low_features_high_area"),
                                    levels(Idents(xenium_data[[TMA]][[segmentation_method]])))
    
    # Subset cells that should be filtered
    if (length(available_cell_tags) > 0){
      temp_subset = subset(xenium_data[[TMA]][[segmentation_method]], idents = available_cell_tags)
    } else {
      temp_subset = xenium_data[[TMA]][[segmentation_method]]
    }
    
    # Set back original ident
    Idents(xenium_data[[TMA]][[segmentation_method]]) = "orig.ident"
    
    cat("\n\n")
    tissue_not_conserved_qc = ImageDimPlot(temp_subset, 
                                           group.by = "QC_result", 
                                           boundaries = "segmentation", 
                                           border.size = 0.05, 
                                           cols = QUALITY_CONTROL_COLOR_PALETTE, 
                                           crop = T)
    print(tissue_not_conserved_qc)
    cat("\n\n")
    
    ################################
    # 2.5. GENE COUNT DISTRIBUTION #
    ################################
    
    cat("\n\n")
    cat(sprintf("%s Top gene distribution ", strrep("#", 4)))
    cat("\n\n")
    
    sample_counts = LayerData(object = xenium_data[[TMA]][[segmentation_method]], layer = "counts")
    
    # % of library per gene
    sample_counts_geneSum = rowSums(sample_counts)
    sample_counts_geneSum = sample_counts_geneSum[order(sample_counts_geneSum, decreasing = T)]
    sample_counts_geneSumPercent = sample_counts_geneSum/sum(sample_counts)*100
    sample_counts_geneSumPercent_df = as.data.frame(sample_counts_geneSumPercent)
    colnames(sample_counts_geneSumPercent_df)[1] = "percent"
    sample_counts_geneSumPercent_df$Gene = rownames(sample_counts_geneSumPercent_df)
    sample_counts_geneSumPercent_df$Gene = factor(sample_counts_geneSumPercent_df$Gene, levels = sample_counts_geneSumPercent_df$Gene)
    
    # Count per gene all cells
    sample_counts_ordered = sample_counts[names(sample_counts_geneSum),]
    sample_counts_ordered_top25 = sample_counts_ordered[1:25,]
    sample_counts_ordered_top25_melt = melt(as.data.frame(t(sample_counts_ordered_top25)))
    
    # Top gene distribution
    percentLibTop25 = ggplot(sample_counts_geneSumPercent_df[1:25,], aes(x = forcats::fct_rev(Gene), y = percent)) + 
      geom_bar(stat = "identity", fill = COLOR_PALETTE[1], col = "black") + 
      coord_flip() +
      theme_bw() + 
      ylab("% of library") +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
    
    nCountGeneTop25 = ggplot(sample_counts_ordered_top25_melt, aes(x=forcats::fct_rev(variable), y = value, fill = variable)) +
      geom_boxplot() +
      scale_fill_manual(values = hue_pal()(nrow(sample_counts_ordered_top25))) +
      coord_flip() +
      theme_bw() +
      ylab("nCount") +
      theme(legend.position = "none", axis.title.y = element_blank())
    
    cat("\n\n")
    top_genes = plot_grid(nCountGeneTop25, percentLibTop25, rel_widths = c(1, 0.5))
    print(top_genes)
    cat("\n\n")
    
  }
}

##################################
# 2.6. DATA STATISTICS HISTOGRAM #
##################################

# Store data statistics to display later
data_qc = data.frame()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Create data frame with data statistics of current condition
    temp_data_qc = data.frame(sample_ID = paste0(TMA, "_", segmentation_method), 
                              qc_tag = xenium_data[[TMA]][[segmentation_method]]@meta.data$QC_result)
    
    # Concatenate to global data statistics
    data_qc = rbind(data_qc, temp_data_qc)
    
  }
}

cat("\n\n")
cat(sprintf("%s Histograms", strrep("#", 2)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

# Retrieve data to plot
plot_qc_data = as.data.frame(table(data_qc))

# Plot corresponding barplot for each qc tag
for (tag in c("conserved", "low_feature", "high_area", "low_features_high_area")){
  
  cat("\n\n")
  cat(sprintf("%s %s ", strrep("#", 3), tag))
  cat("\n\n")
  
  # Subset df to keep current qc tag
  temp_data_qc = plot_qc_data[plot_qc_data$qc_tag == tag, ]
  
  cat("\n\n")
  print(ggplot(temp_data_qc, aes(x = sample_ID, y = Freq, fill = sample_ID)) +
          geom_bar(stat = "identity") +
          geom_text(aes(label = Freq), vjust = -0.5, size = 4) +
          theme_minimal() +
          scale_fill_manual(values = COLOR_PALETTE) +
          ggtitle(sprintf("%s cells over segmentation methods", tag)) +
          labs(x = "Sample ID", y = "Frequency") + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  cat("\n\n")
  
}

## @knitr compare_statistics_distributions

################################################################################
#'                      3. PLOT DISTRIBUTIONS OF DATA STATISTICS
################################################################################

###################################
## 3.1. RETRIEVE DATA STATISTICS ##
###################################

# Store data stastitics to display later
data_statistics = data.frame()

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  # Iterate over segmentation methods covered
  for (segmentation_method in names(TMA_METADATA[[TMA]])){
    
    # Create data frame with data stastistics of current condition
    temp_data_statistics = data.frame(sample_ID = paste0(TMA, "_", segmentation_method), 
                                      cell_area = xenium_data[[TMA]][[segmentation_method]]@meta.data$cell_area,
                                      n_feature = xenium_data[[TMA]][[segmentation_method]]@meta.data$nFeature_Xenium,
                                      transcript_count = xenium_data[[TMA]][[segmentation_method]]@meta.data$transcript_count)
    
    # Concatenate to global data statistics
    data_statistics = rbind(data_statistics, temp_data_statistics)
    
  }
}

#################################
## 3.2. DISPLAY TMA STATISTICS ##
#################################

cat("\n\n")
cat(sprintf("%s **TMA statistics**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("This section generates ridge plots that compare the distributions of key 
cell-level statistics across each TMA and segmentation method under study.")
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Cell area ", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(data_statistics, aes(x = cell_area, y = sample_ID, fill = sample_ID)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(0, 300)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
        geom_text(data=data_statistics %>% group_by(sample_ID) %>% 
                    summarise(cell_area=median(cell_area)),
                  aes(label=sprintf("%1.1f", cell_area)), 
                  position=position_nudge(x = 10, y=0.8), colour="black", size=3.5) +
        labs(title = "Cell area distribution across segmentation methods"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Features per cell ", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(data_statistics, aes(x = n_feature, y = sample_ID, fill = sample_ID)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(0, 300)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
        geom_text(data=data_statistics %>% group_by(sample_ID) %>% 
                    summarise(n_feature=median(n_feature)),
                  aes(label=sprintf("%1.1f", n_feature)), 
                  position=position_nudge(x = 10, y=0.8), colour="black", size=3.5) +
        labs(title = "Features per cell distribution across segmentation methods"))
cat("\n\n")

cat("\n\n")
cat(sprintf("%s Transcripts per cell ", strrep("#", 2)))
cat("\n\n")

cat("\n\n")
print(ggplot(data_statistics, aes(x = transcript_count, y = sample_ID, fill = sample_ID)) +
        geom_density_ridges() +
        theme_ridges() + 
        theme(legend.position = "none") +
        scale_x_continuous(limits = c(0, 2000)) +
        stat_density_ridges(quantile_lines = TRUE, quantiles = 2, vline_color = "black") +
        geom_text(data=data_statistics %>% group_by(sample_ID) %>% 
                    summarise(transcript_count=median(transcript_count)),
                  aes(label=sprintf("%1.1f", transcript_count)), 
                  position=position_nudge(x = 55, y=0.6), colour="black", size=3.5) +
        labs(title = "Transcript count distribution across segmentation methods"))
cat("\n\n")
