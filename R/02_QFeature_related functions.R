visualize_imputation_density <- function(data, 
                                         assayName_ = "NAfiltered_proteins", 
                                         imputation_methods = c("knn", "zero", "MinDet", "bpca"),
                                         colors = c("black", "red", "blue", "steelblue", "orange"),
                                         legend_position = "topright") {
  # Extract protein data
  protein_data <- assay(data[[assayName_]])
  
  # Remove NAs from original data for comparison
  orig_density <- density(na.omit(protein_data))
  
  # Plot original data density
  plot(orig_density, col = colors[1], main = "Imputation Density Comparison", 
       xlab = "Density", ylab = "Frequency")
  
  # Loop through imputation methods and add density lines
  for (i in seq_along(imputation_methods)) {
    method <- imputation_methods[i]
    imputed_data <- assay(impute(data[[assayName_]], method = method))
    lines(density(imputed_data), col = colors[i + 1])
  }
  
  # Add legend
  legend(legend_position, 
         legend = c("Original", imputation_methods), 
         col = colors[1:(length(imputation_methods) + 1)], 
         lwd = 2, bty = "n")
}
