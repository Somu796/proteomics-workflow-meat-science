
#' General Function to create

generate_dimension_reduction_plot <- function(model, x, y, model_output_data, label, labels, explained_variance = NULL) {
  
  # pca data and corresponding "quantCols" name
  model_output_data <- model_output_data  %>%
    rownames_to_column("quantCols")
  
  # group label and corresponding "quantCols" name
  group_names <- labels %>% 
    dplyr::select(all_of(label))%>% 
    rownames_to_column("quantCols")
  
  colnames(group_names)[2] <- "group" 
  
  # joining pca data and group label by "quantCols"
  model_output_data <- model_output_data %>% 
    full_join(group_names, by = "quantCols") %>% 
    column_to_rownames("quantCols")
  
  
  p <-ggplot(model_output_data,
           aes(x = model_output_data[, x], y = model_output_data[, y], color = group)) +
    geom_point(size = 4, alpha = 5 / 10) +  # Marker size
    geom_hline(yintercept = 0, linetype="dashed", alpha = 3/10) +
    geom_vline(xintercept = 0, linetype="dashed", alpha = 3/10) +
    geom_polygon(stat = "ellipse",
                 aes(fill = group),
                 alpha = 0.001) +
    theme_bw(base_family = "serif")
  
  if (tolower(model) == "pca") {
    p <- p +
      labs(
        x = glue("Principal Component {x} ({explained_variance[x]}%)"),
        y = glue("Principal Component {y} ({explained_variance[y]}%)")
      )
    
  } else if (tolower(model) == "tsne"){
    p <- p  +
      labs(
        x = glue("tSNE {x}"),
        y = glue("tSNE {y}")
      )
  } else if (tolower(model) == "umap"){
    p <- p  +
      labs(
        x = glue("UMAP {x}"),
        y = glue("UMAP {y}")
      )
  } else{
    p <- p  +
      labs(
        x = glue("Component {x}"),
        y = glue("Component {y}")
      )
  }
  
  # show(p)
  # fig_scorePlot <- ggplotly(p)
  # fig_scorePlot
  
  return(p)
}




generate_pca_scree_plot <- function(scree_plot_data, n){
  scree_plot <- ggplot(scree_plot_data[1:n,], mapping = aes(x = reorder(PC_component, explained_variance), y = explained_variance)) +
    geom_bar(stat = "identity", fill="#56B4E9", colour="black") +
    geom_line(aes(x = 1:n, y = explained_variance), linewidth = 1) +
    geom_text(aes(label = round(explained_variance, 2)), vjust = -0.5, size = 3)+
    labs(x = "Principal Componenets",
         y = "Explained Variance (%)") +
    theme_bw(base_family = "serif")
  return(scree_plot)
}
