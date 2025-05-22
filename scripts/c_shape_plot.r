# Load required libraries
library(plotly)

# Load scripts
source("scripts/utils.r")

# cluster_shapes function: Plots the mean shape and individual shapes of clusters
# Input:
#   - data: procGPA object containing the rotated shapes and mean shape
#   - joinlines: List of lists containing indices for connecting landmarks
#   - clusters: Vector of cluster assignments for each shape
#   - color_mapping: Matrix of colors for each cluster
#   - show_landmarks_dots: Boolean indicating whether to show individual landmark dots
#   - output_path: Path to save the combined plot
#   - output_individual_path: Path to save individual cluster plots
#   - annotation_plot_labels: Labels for the clusters
#   - line_width: Width of the lines in the plot
#   - combined_plot_mode: Mode for combining plots ("h" for horizontal, "grid" for grid layout)
# Output:
#   - combined_plot: Combined plot of all clusters
#   - plot_list: List of individual cluster plots
cluster_shapes <- function(data, joinlines, clusters, color_mapping, show_landmarks_dots = FALSE,
                           output_path = NULL, output_individual_path = NULL,
                           annotation_plot_labels = "Cluster", line_width = 2, combined_plot_mode = "h") {

  # Validate inputs
  if (is.null(data) || is.null(data$rotated) || is.null(data$mshape)) {
    stop("Invalid data input. Ensure 'data' is a valid procGPA object with 'rotated' and 'mshape'.")
  }

  # Remove invalid cluster values (0, -1, NaN)
  valid_clusters <- unique(clusters[!clusters %in% c(0, -1, NaN)])
  k <- length(valid_clusters)  # Number of valid clusters

  # Create folder for individual plots if needed
  if (!is.null(output_individual_path)) {
    individual_path <- file.path(output_individual_path, paste0(annotation_plot_labels, "_shapes"))
    dir.create(individual_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Compute overall shape variation (dispersion)
  rotated_sd <- apply(data$rotated, c(1, 2), sd)
  polygon_x <- cbind(data$mshape[, 1] - rotated_sd[, 1], data$mshape[, 1] + rotated_sd[, 1],
                     data$mshape[, 1] + rotated_sd[, 1], data$mshape[, 1] - rotated_sd[, 1])
  polygon_y <- cbind(data$mshape[, 2] - rotated_sd[, 2], data$mshape[, 2] - rotated_sd[, 2],
                     data$mshape[, 2] + rotated_sd[, 2], data$mshape[, 2] + rotated_sd[, 2])
  
  # Initialize lists
  plot_list <- vector("list", k)
  annotations <- list()
  
  # Define legend positions
  legend_x_positions <- seq(0.1, 0.88, length.out = k)
  
  # Loop through each cluster
  for (cluster_idx in 1:k) {
    # Compute mean shape and standard deviation for the current cluster
    cluster_members <- which(clusters == cluster_idx)
    cluster_mean_shape <- apply(data$rotated[, , cluster_members], c(1, 2), mean)
    cluster_std <- apply(data$rotated[, , cluster_members], c(1, 2), sd)
    
    # Generate shape plot using plot_mshapes()
    p <- plot_mshapes(
      shapes = list(data$mshape, cluster_mean_shape),
      joinlines = joinlines,
      shape_labels = c("mean", paste(annotation_plot_labels, cluster_idx)),
      color = c("#948c8c", color_mapping[cluster_idx, ]),
      line_width = line_width
    ) %>%
      layout(showlegend = FALSE)
    
    # Add cluster annotation
    annotations <- c(annotations, list(
      list(
        x = legend_x_positions[cluster_idx], y = 1.05,
        text = paste0("<b>", annotation_plot_labels, " ", cluster_idx, "</b>"),
        showarrow = FALSE, xref = "paper", yref = "paper", align = "left"
      )
    ))
    
    # Add landmark dots if enabled
    if (show_landmarks_dots) {
      for (i in cluster_members) {
        p <- add_trace(p,
                       type = "scatter",
                       mode = "markers",
                       x = data$rotated[, 1, i],
                       y = data$rotated[, 2, i],
                       marker = list(color = color_mapping[cluster_idx, ], symbol = "circle", size = 3, opacity = 0.3),
                       showlegend = FALSE)
      }
    }
    
    # Draw landmark dispersion
    for (i in 1:nrow(cluster_mean_shape)) {
      p <- p %>%
        add_trace(
          x = c(cluster_mean_shape[i, 1], cluster_mean_shape[i, 1]), 
          y = c(cluster_mean_shape[i, 2] - cluster_std[i, 2], cluster_mean_shape[i, 2] + cluster_std[i, 2]),
          type = "scatter", mode = "lines",
          line = list(color = color_mapping[cluster_idx, ], width = line_width),
          showlegend = FALSE
        ) %>%
        add_trace(
          x = c(cluster_mean_shape[i, 1] - cluster_std[i, 1], cluster_mean_shape[i, 1] + cluster_std[i, 1]), 
          y = c(cluster_mean_shape[i, 2], cluster_mean_shape[i, 2]),
          type = "scatter", mode = "lines",
          line = list(color = color_mapping[cluster_idx, ], width = line_width),
          showlegend = FALSE
        ) %>%
        add_polygons(
          x = polygon_x[i, ], y = polygon_y[i, ],
          type = "scatter", mode = "lines",
          fill = "toself", fillcolor = "rgba(128, 128, 128, 0.3)",
          line = list(color = "transparent"), showlegend = FALSE
        )
    }
    
    # Store the plot
    plot_list[[cluster_idx]] <- p
    
    # Save individual cluster plot if required
    if (!is.null(output_individual_path)) {
      save_image_svg(
        p = p %>% layout(
          font = list(size = 19, family = "Helvetica"),
          title = list(text = paste("Cluster", cluster_idx), font = list(size = 40, family = "Helvetica"), x = 0.5, xanchor = "center"),
          margin = list(t = 50)
        ),
        output_path = file.path(individual_path, paste0(annotation_plot_labels, "_", cluster_idx, ".svg")),
        width = 700, height = 700
      )
    }
  }
  
  # Combine all cluster plots into one subplot
  if (combined_plot_mode == "h") {
  combined_plot <- subplot(plot_list, shareX = FALSE, shareY = FALSE, margin = 0.05) %>%
    layout(annotations = annotations, legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.1))} 
    else if (combined_plot_mode == "grid") {
      combined_plot <- subplot(plot_list, nrows = 2, margin = 0.05)
    }

  # Save combined plot if required
  if (!is.null(output_path)) {
    save_image_svg(p = combined_plot, output_path = paste0(output_path, "cluster_shape.svg"), width = 1400, height = 700)
  }
  return(list(combined_plot = combined_plot, plot_list = plot_list))

}
