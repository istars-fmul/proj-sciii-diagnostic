library(umap) # perform UMAP
library(dplyr)
library(plotly)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)

# plot_wilcox_adjusted_pheatmap function: Plot a heatmap of adjusted p-values from Wilcoxon tests
# Input: 
# - table_cluster_stats_wilcox: Data frame containing the results of the Wilcoxon tests
# - color_scheme: Color scheme for the heatmap (default: c("lightyellow", "red"))
# - cutoff: Cutoff value for the p-values (default: NULL)
# - output_path: Path to save the heatmap (default: NULL)
# Output: Heatmap of adjusted p-values
plot_wilcox_adjusted_pheatmap <- function(table_cluster_stats_wilcox, color_scheme = c("lightyellow", "red"), cutoff = NULL, output_path = NULL) {
  # Pivot the data to have features as rows and comparisons as columns
  p_adjusted_pivot <- table_cluster_stats_wilcox %>%
    pivot_wider(
      values_from = log_p_adjusted,
      id_cols = feature,
      names_from = comparison
    )

  # Set feature as row names
  p_adjusted_pivot <- p_adjusted_pivot %>%
    column_to_rownames(var = "feature")

  # Apply statistical test (-log10 values) with a lower cutoff (if specified)
  p_adjusted_pivot <- apply(p_adjusted_pivot, c(1, 2), function(x) ifelse(x < -log10(0.05), NA, x))

  # Plot heatmap with adjusted p-values
  heatmap_p_adjusted <- pheatmap(t(as.matrix(p_adjusted_pivot)),
            color = colorRampPalette(color_scheme)(50),
            na_col = "grey",
            main = "Wilcoxon test -log10(adjusted p-value)",
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            legend = TRUE,
            angle_col = 90,
            fontsize = 18
  )
  # Save the plot if an output path is specified
  if (!is.null(output_path)) {
       save_pheatmap_pdf(heatmap_p_adjusted, output_path, width = 20, height = 10)
  }

  return(heatmap_p_adjusted)
}

# feature_cluster_stats_t_test function: Perform a Wilcoxon rank sum test for each feature between clusters
# Input: data - a data frame containing the feature values for each sample
#        partition - a vector containing the cluster assignments for each sample
# Output: stats_test_df - a data frame containing the results of the statistical tests
cluster_stats_wilcox <- function(data, partition) {
  # Initialize an empty dataframe to store the results
  stats_test_df <- data.frame(comparison = character(),
                              p_value = numeric(),
                              cluster_n1 = numeric(),
                              cluster_n2 = numeric(),
                              statistic_value = numeric(),
                              feature = character(),
                              stringsAsFactors = FALSE)
  
  # get number of clusters
  k <- length(unique(partition))
  
  # Loop through each feature (column) in the data frame
  for (feature in colnames(data)) {
    # Loop through each unique cluster
    for (i in 1:(k - 1)) {
      feature_i <- data[which(partition == i), feature]
      for (j in (i+1):k) {
        feature_j <- data[which(partition == j), feature]
        comparison <- paste0("C", i, "-vs-C", j)

        # Perform Wilcoxon rank-sum test
        test_result <- wilcox.test(feature_i, feature_j)
        p_value <- test_result$p.value
        stats_value <- test_result$statistic


        # Append the result to the stats_feature_test data frame
        stats_test_df <- rbind(stats_test_df, data.frame(comparison = comparison, 
                                                         cluster_n1 = i,
                                                         cluster_n2 = j,
                                                         feature = feature,
                                                         statistic_value = stats_value,
                                                         p_value = p_value))
      }
    }
  }
  
  # Adjust p-values using the Bonferroni method and calculate log adjusted p-values
  stats_test_df$p_adjusted <- p.adjust(stats_test_df$p_value, method = "bonferroni")
  stats_test_df$log_p_adjusted <- -log10(stats_test_df$p_adjusted)
  stats_test_df$mean_is_different <- stats_test_df$p_adjusted < 0.05
  
  return(stats_test_df %>% arrange(comparison, feature))
}

# heatmap_aov_cluster function: generate a heatmap of AOV results per cluster
# Input:
#       - df: data frame which contains the data to perform AOV
#       - cluster: cluster information
#       - annotation_col: data frame with annotation colors 
#       - annotation_colors: vector of colors for the annotation
#       - color_scheme: color scheme for the heatmap
#       - title: title of the heatmap
# Output:
#       - pheatmap object of the heatmap

heatmap_aov_cluster <- function(df, cluster, annotation_col, annotation_colors, color_scheme = c("blue", "white", "red"), title = "AOV - Cluster", output_path) {
  # Perform AOV analysis - response variable; grouping variable is the cluster; -1 removes the intercept
  Cluster <- as.factor(cluster)  # Use the provided 'cluster' input
  aov_result <- aov(as.matrix(df) ~ Cluster - 1)
  
  # Scale the coefficients - represent the differences in means between each group and the reference group
  aov_scaled_coefficients <- scale(aov_result$coefficients)
  
  # Plot heatmap
  heatmap <- pheatmap(
    aov_scaled_coefficients,
    color = colorRampPalette(color_scheme)(300),
    main = title,
    annotation_colors = annotation_colors,  # Colors for the anatomical regions
    annotation_col = annotation_col,  # Mapping column names to anatomical regions
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    legend = TRUE,
    fontsize = 18,
    angle_col = 90
  )

  # Save the heatmap if output path is provided
  if (!is.null(output_path)) {
     save_pheatmap_pdf(heatmap, output_path, width = 20, height = 10)
  }
  return(heatmap)
}
# generate_umap_plot function: generate a plot of UMAP with color by group
# Input:
#       - df: data frame with UMAP layout (UMAP1, UMAP2)
#       - df_group: data frame with group information
#       - group_label: column name of the group information
#       - color_mapping: list of colors for the groups
# Output:
#       - plotly plot of UMAP with color by groups
generate_umap_plot <- function(df, df_group, group_label = "Group", color_mapping = NULL, path = NULL) {
  # Perform UMAP
  umap_result <- umap(df, n_components = 2, random_state = 15)
  
  # Transform the UMAP result to a data frame
  data <- as.data.frame(umap_result$layout)
  colnames(data) <- c("UMAP1", "UMAP2")
  rownames(data) <- rownames(df)

  # Bind the grouping variable to the data frame
  data[[group_label]] <- df_group[[group_label]]  # Use the correct grouping column

  # Convert the grouping variable to a factor if not already
  data[[group_label]] <- as.factor(data[[group_label]])
  
  # Generate the UMAP plot using plotly
  p <- plot_ly(
    data = data,
    x = ~UMAP1,
    y = ~UMAP2,
    color = ~data[[group_label]],
    colors = setNames(unlist(color_mapping), levels(data[[group_label]])),
    type = "scatter",
    mode = "markers",
    text = ~rownames(data),
    marker = list(size = 7),
    showlegend = TRUE,
    width = 600,
    height = 500
  ) %>%
    layout(
      title = "",
      legend = list(
        orientation = "v",  # Set legend orientation to vertical
        x = 1, y = 1, title = list(text = group_label),  # Position the legend dynamically
        xanchor = "center", yanchor = "top",
        font = list(family = "Helvetica", size = 15)  # Set font for the legend
      )
    )
  return(p)

  if (!is.null(path)) {
    save_image_svg(p, output_path = path, width = 600, height = 500)
    }
}

# generate_mean_sd_cluster function: generate a data frame with mean and standard deviation for each cluster
# Input:
#       - df: data frame with coordinates
#       - cluster: cluster information
#       - path: path to save the data frame as a .CSV file
# Output:
#       - data frame with mean and standard deviation for each cluster
generate_mean_sd_cluster <- function(df, cluster, path) {
  # Merge the data frame with the cluster information
  mean_sd_cluster <- merge(df, cluster, by = "row.names", all.x = TRUE) %>%
    # Remove the row names column
    select(-Row.names) %>%
    # Group by the cluster and calculate the mean and standard deviation
    group_by(Cluster) %>%
    summarize(across(everything(), ~paste0(round(mean(.), 2), " ± ", round(sd(.), 2))), .groups = "drop") %>%
    as.data.frame() %>%
    column_to_rownames("Cluster") %>%  # Set 'Cluster' column as row names
    t()

  
  # Save data frame to .CSV file
  if (!is.null(path)) {
    write.csv(mean_sd_cluster, path, row.names = TRUE)
  }

  return(mean_sd_cluster)
}

# violin_plot_cluster function: generate a violin plot of df per cluster
# Input:
#       - df: data frame with measurements and cluster information it should be in a format of a long data frame 
#       - color_mapping: vector of colors for the clusters
#       - y_var: y-axis variable name
# Output:
#       - ggplot object of the violin plot
violin_plot_function <- function(df, color_mapping, group_var = "Measurement", y_label ="Ceohalometric Variable", output_path = NULL, width = 13, height = 20) {
  # Ensure Cluster is a factor
  df$Cluster <- as.factor(df$Cluster)

  # Create violin plot
  violin_plot <- ggplot(df, aes(x = Cluster, y = Value, fill = Cluster)) +
    geom_violin() +
    stat_summary(fun.data = mean_sdl, geom = "crossbar", width = 0.1) +
    scale_fill_manual(values = color_mapping) +
    facet_wrap(as.formula(paste("~", group_var)), scales = "free", ncol = 4) +
    theme_minimal() +
    labs(x = "Cluster", y = y_label) +
  theme(
    text = element_text(size = 12, family = "Helvetica"),
    axis.title = element_text(size = 13, family = "Helvetica"),
    axis.text = element_text(size = 11, family = "Helvetica"),
    strip.text = element_text(size = 13, family = "Helvetica"),
    legend.title = element_text(size = 11, family = "Helvetica"),
    legend.text = element_text(size = 11, family = "Helvetica")
  )
  
  # Save the plot if output path is provided
  if (!is.null(output_path)) {
    ggsave(filename = output_path, plot = violin_plot, device = "svg", width = width, height = height)
  }

  return(violin_plot)
}