
## @knitr plot_mecr_scores

################################################################################
#'                              2. PLOT MECR SCORES
################################################################################

cat("\n\n")
cat(sprintf("%s **MECR Boxplots across hyperparameter sets**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("Display the MECR score distribution for each hyperparameter set used during 
the 02_Cell_Segmentation step. This enables visualization of the improvements 
provided by Proseg compared to the initial Xenium segmentation. It also provides 
an overview of how changing Proseg hyperparameters influence MECR values across 
the dataset.")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA))
  cat("\n\n")
  
  # Prepare data to plot box
  boxplot_data = bind_rows(lapply(names(MECR_scores[[TMA]]), function(n){
    data.frame(sample_ID = n, MECR = MECR_scores[[TMA]][[n]]$MECR)
  }))
  
  # Chose color for labeling
  label_colors = ifelse(levels(factor(boxplot_data$sample_ID)) == XENIUM_CONDITION_NAME, "red", "black")
  
  # Make boxplot
  boxplot_figure = ggplot(boxplot_data, aes(x = sample_ID, y = MECR)) +
    geom_boxplot() +
    theme_bw() +
    theme(
      axis.text.x = element_text(color = label_colors, angle = 90, hjust = 1)
    )
  
  cat("\n\n")
  print(boxplot_figure)
  cat("\n\n")
  
}
