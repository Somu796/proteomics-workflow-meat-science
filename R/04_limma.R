# library(limma)

# Get DAP List
getDAPlist <- function(fit_contrasts_bayes, contrasts_group, 
                       information_DE = c("logFC", "AveExpr", "adj.P.Val"), 
                       multiple_testing_adjusting_method = "none") {
  
  # Initialize the DAP_table as NULL
  DAP_table <- NULL
  not_in_topTable <- NULL
  DAP_table_contrast <- NULL
  
  # Iterate over each contrast in contrasts_group
  for (i in colnames(contrasts_group)) {
    # Get the topTable for the current contrast
    DAP_table_contrast <- topTable(
      fit_contrasts_bayes,
      coef = i,
      n = Inf,
      adjust.method = multiple_testing_adjusting_method,
      sort.by = "none"
    )
    
    # If DAP_table is NULL (first iteration), initialize it
    if (is.null(DAP_table)) {
      DAP_table <- data.frame(row.names = rownames(DAP_table_contrast))  # Set row names
    }
    
    # Add the specified columns (logFC, AveExpr, adj.P.Val) for the current contrast
    for (j in information_DE) {
      # checking user input is correct
      if (!j %in% colnames(DAP_table_contrast)){
        not_in_topTable <- append(not_in_topTable, j)
      }
      
      DAP_table[, paste0(j, "_", i)] <- DAP_table_contrast[rownames(DAP_table), j]
    }
  }
  
  if (length(not_in_topTable) > 0) {
    message(paste0(
      "\nWarning: The following columns were not found in topTable output and were skipped: ",
      paste(unique(not_in_topTable), collapse = ", "),
      "\nAvailable columns in topTable output are: ",
      paste(colnames(DAP_table_contrast), collapse = ", ")
    ))
  }
  
  return(DAP_table)
}

# Volcano Plot

getVolcanoPlot = function(DAP_table_factor, FactorLevel.1, FactorLevel.2, mrlimit = 1.2, pValue_Limit= 0.05, legend_title = NULL){
  #on-off proteins(only present in one group
  #Factor eg. Horns, RC, SC
  
  differential_data <- data.frame(gene_name = DAP_table_factor$Gene_primary , 
                                  meanRatio = 2^DAP_table_factor[[colnames(DAP_table_factor)[grepl("logFC",colnames(DAP_table_factor))]]], 
                                  pValue = DAP_table_factor[[colnames(DAP_table_factor)[grepl("adj.P.Val",colnames(DAP_table_factor))]]], 
                                  differential_group = "No") # forming a new differential data for plotting
  
  # Adding a row to inform Up regulated or Down regulated
  differential_data$differential_group[(differential_data$meanRatio> mrlimit) & differential_data$pValue< pValue_Limit] = "Up"
  differential_data$differential_group[(differential_data$meanRatio< 1/mrlimit) & differential_data$pValue< pValue_Limit] = "Down"
  
  # Making a factor of the differential group to assign color to it
  differential_data = differential_data %>%
    mutate(differential_group = factor(differential_group, levels = c("Up", "No", "Down")))
  
  # Assigning the color
  regulation_color <- setNames(c("#FF7F0E", "grey", "#1F77B4"), # #bb0c00
                               levels(differential_data$differential_group))
  
  # getting the plots for plotting
  differential_data$delabel <- ifelse(differential_data$differential_group == "Up" | differential_data$differential_group == "Down", 
                                      differential_data$gene_name, NA)
  # Total up, down and not sig entries 
  total_Up = sum(differential_data$differential_group == "Up")
  print(total_Up)
  total_Down = sum(differential_data$differential_group == "Down")
  print(total_Down)
  total_Not_Sig = sum(is.na(differential_data$delabel))
  print(total_Not_Sig)
  
  # # p-Adj Values Calculation
  # if (pAdjust == TRUE){
  #   differential_data$fdr.p.value <- p.adjust(differential_data$pValue, method = pAdjustMethod)
  # }
  
  # Plotting the result
  
  
  volcano_plot = ggplot(data = differential_data, 
                        mapping = aes(log2(meanRatio), 
                                      -log10(pValue), 
                                      col= differential_group #, 
                                      # label = delabel
                        )) +
    geom_point(alpha = 5/10, 
               size = 2) +
    geom_vline(xintercept = c(-log2(mrlimit), 
                              log2(mrlimit)), 
               col = "black", 
               linetype = 'dashed') +
    geom_hline(yintercept = -log10(pValue_Limit), 
               col = "black", 
               linetype = 'dashed') + 
    # to overcome the text overlap
    # geom_text_repel(aes(family = "serif"),
    # max.overlaps = Inf
    # position =
    #   position_nudge_to(x = 2.3),
    # min.segment.length = 0,
    # segment.color = "black",
    # arrow = arrow(length = unit(0.015, "npc")),
    # direction = "y"
    # ) +
    labs(color = legend_title, 
         x = substitute(log[2] ~ "Fold Change (" * F/C * ")", list(F = FactorLevel.1, C = FactorLevel.2)) ,
         y = expression("-log"[10]*"(Pvalue)"),
         caption = glue("Pvalue < {pValue_Limit}; Fold Change > {mrlimit} \n Up = {total_Up}; Down = {total_Down}; Not sig. = {total_Not_Sig}")) +
    scale_color_manual(values = regulation_color, # to set the colours of our variable  
                       labels = c("Up", "Not sig.", "Down")) +
    theme_minimal() + 
    theme(text = element_text(family = "serif", size = 12),
          axis.line = element_line(),
          axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
          axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
          plot.title = element_text(hjust = 0.5)
    )
  
  return(volcano_plot)
  
}

# Upste Plot

getUpSetPlot <- function(upSet_data_stress, fill_color = "grey", title = ""){
  # Convert the table into a data frame
  upSet_data_stress_df <- as.data.frame(upSet_data_stress)
  
  # Rename columns for clarity
  colnames(upSet_data_stress_df) <- c("Combination", "Count")
  
  # Split the "Combination" column into individual sets
  upSet_data_stress_df <- upSet_data_stress_df %>%
    mutate(Sets = strsplit(as.character(Combination), "&")) # Convert to a list column
  
  # Create the ggupset plot
  ggplot(upSet_data_stress_df, aes(x = Sets, y = Count)) +
    geom_bar(stat = "identity", fill = fill_color) +
    geom_text(aes(label = Count), vjust = -0.5, size = 4) + # Add counts on top of bars
    scale_x_upset() +
    labs(
      title = title,
      y = "Intersection Size"
    ) +
    theme_minimal()+
    theme(title = element_text(family = "serif", size = 14, colour = "black"),
          plot.title.position = "plot",
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(colour = "black"),
          # axis.text.x = element_text(colour = "black", size = 14),
          axis.title.y = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black", size = 10))
  
}

# 
# 
# # Generalized function for factorial design analysis
# factorial_limm <- function(data, outcome_var, factor_vars) {
#   # Check if outcome_var and factor_vars are in the data
#   if (!outcome_var %in% colnames(data)) {
#     stop("The specified outcome variable is not in the data frame.")
#   }
#   if (!all(factor_vars %in% colnames(data))) {
#     stop("One or more specified factor variables are not in the data frame.")
#   }
#   
#   # Combine factors into a single variable
#   data$Combined <- factor(do.call(paste, c(data[factor_vars], sep = "_")))
#   
#   # Create the design matrix without intercept
#   design <- model.matrix(~ 0 + Combined, data = data)
#   colnames(design) <- levels(data$Combined)
#   
#   # Fit the linear model
#   fit <- lmFit(data[[outcome_var]], design)
#   
#   # Generate all pairwise contrasts
#   contrast_names <- combn(levels(data$Combined), 2, FUN = function(x) paste(x[1], "vs", x[2], sep = "_"))
#   contrast_matrix <- makeContrasts(contrasts = paste(levels(data$Combined), collapse = " - "), 
#                                    levels = design)
#   
#   # Compute all pairwise contrasts
#   contrast_matrix <- combn(levels(data$Combined), 2, FUN = function(x) {
#     contrast <- rep(0, length(levels(data$Combined)))
#     contrast[which(levels(data$Combined) == x[1])] <- 1
#     contrast[which(levels(data$Combined) == x[2])] <- -1
#     contrast
#   }, simplify = TRUE)
#   colnames(contrast_matrix) <- contrast_names
#   
#   # Apply contrasts to the model
#   fit2 <- contrasts.fit(fit, contrast_matrix)
#   
#   # Perform empirical Bayes moderation
#   fit2 <- eBayes(fit2)
#   
#   # Return results
#   list(
#     fit = fit2,
#     top_table = topTable(fit2, number = Inf),
#     contrasts = contrast_names
#   )
# }
# 
# set.seed(123)  # For reproducibility
# 
# # Generate dummy data
# dummy_data <- data.frame(
#   Factor1 = factor(rep(c("A", "B"), each = 6)),         # Two levels in Factor1
#   Factor2 = factor(rep(c("X", "Y", "Z"), times = 4)),   # Three levels in Factor2
#   Outcome = rnorm(12, mean = 5, sd = 1)                # Random outcome variable
# )
# 
# # View the dummy data
# print(dummy_data)
# 
