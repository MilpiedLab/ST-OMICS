
## @knitr plot_correlations

################################################################################
#'                                3. PLOT CORRELATIONS
################################################################################

cat("\n\n")
cat(sprintf("%s **Observe correlation between response variables**", strrep("#", 1)), " {.tabset .tabset-fade}\n\n")
cat("\n\n")

cat("\n\n")
cat("The goal is to plot the statistics retrieved in the previous step against 
each other. This kind of visualization enables rapid identification of the Proseg 
hyperparameter set that yields the highest-quality segmentation. For each figure, 
observations are colored according to one of the evaluated Proseg parameters.

The selected hyperparameter set is defined as the one that produces a segmentation with:

  - Low MECR values
  - Few low-feature cells
  - A number of conserved cells similar to that of the initial Xenium segmentation

In these Proseg runs, it is commonly observed that as MECR decreases, the number 
of low-feature cells increases. This trade-off makes Proseg parameter selection 
challenging. The objective is therefore to identify a parameter set that achieves 
sufficiently low MECR while maintaining a low number of low-feature cells (or, 
equivalently, preserving a relatively high number of conserved cells).")
cat("\n\n")

# Iterate over available TMAs
for (TMA in TMA_COVERED){
  
  cat("\n\n")
  cat(sprintf("%s %s", strrep("#", 2), TMA), " {.tabset .tabset-fade}\n\n")
  cat("\n\n")
  
  # Iterate over explanatory variable to color observations accordingly
  for (explanatory_variable in names(hyperparameters_values)){
    
    cat("\n\n")
    cat(sprintf("%s %s", strrep("#", 3), explanatory_variable))
    cat("\n\n")
    
    # Remove NA values for current explanatory variable
    temp_statistics_df = TMA_statistics[[TMA]][!(is.na(TMA_statistics[[TMA]][[explanatory_variable]])), ]
    
    # Convert explanatory variable in factor
    dp = factor(temp_statistics_df[[explanatory_variable]])
    
    # Colors per level (legend)
    levs = levels(dp)
    lev_cols = COLOR_PALETTE[seq_along(levs)]
    
    # Correlation panel
    panel.cor = function(x, y){
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = round(cor(x, y), digits=2)
      txt = paste0("R = ", r)
      cex.cor = 0.8/strwidth(txt)
      text(0.5, 0.5, txt)
    }
    # Customize upper panel
    upper.panel = function(x, y){
      points(x,
             y, 
             pch = 19, 
             col = COLOR_PALETTE[as.factor(temp_statistics_df[[explanatory_variable]])])
    }
    
    # Plot graph
    cat("\n\n")
    pairs(temp_statistics_df[, RESPONSE_VARIABLES], 
          lower.panel = panel.cor,
          upper.panel = upper.panel,
          oma=c(3,3,3,15))
    
    # Add legend to pairs graph
    par(xpd = NA)
    legend("bottomright",
           inset = c(-0.01, 0),
           fill = lev_cols,
           legend = levels(dp),
           title = sprintf("%s", explanatory_variable),
           cex = 0.8,
           bg = "white")
    cat("\n\n")
    
  }
}
