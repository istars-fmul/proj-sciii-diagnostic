library(plotly)
library(grid)
library(clue) # for cluster matching (solve lsap)
library(reticulate)

py_install("kaleido==0.2.1")
if (!py_module_available("plotly")) py_install("plotly")


# euclidean_distance function: compute the euclidean distance between two points
# Input:
#       - p1: first point - a vector of coordinates
#       - p2: second point - a vector of coordinates
# Output:
#       - euclidean distance
euclidean_distance <- function(p1, p2) {
  sqrt(sum((p1 - p2)^2, na.rm = TRUE))
}

# gini function: compute the Gini coefficient
# Input:
#       - x: vector of values
# Output:
#       - Gini coefficient
# index near 0 = perfect equality, index near 1 = perfect inequality
gini <- function(x) {
  n <- length(x)
  x <- sort(x)
  G <- sum((2 * (1:n) - n - 1) * x)
  return(G / (n * sum(x)))
}

# save_pheatmap_pdf function: save pheatmap object to pdf
# Input:
#       - x: pheatmap object
#       - filename: name of the file to save
#       - width: width of the pdf
#       - height: height of the pdf
# Output:
#       - pdf file
# example: save_pheatmap_pdf(heat_map_scaled_coefficients, "test.pdf")
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# save_image_svg function: save image plotly to svg
# Input:
#       - p: plot
#       - output_path: path to save the plot
#       - width
#       - heigth
save_image_svg <- function(p, output_path = NULL, width = 1400, height = 700) {
  plotly::save_image(p, file = output_path, width = width, height = height, scale = 1)
  # Return the output path for reference
  return(output_path)
}

# find_positions function: find the positions of the names in the all_names vector
# Input:
#       - names: names to find
#       - all_names: all names
# Output:
#       - positions: positions of the names in the all_names vector
# Get lines to represent anatomical structures
find_positions <- function(names, all_names) {
  positions <- sapply(names, function(name) {
    # Find the position of the name in the all_names vector
    which(all_names == name)
  })
  return(positions)
}

# get_anatomical_lines function: get the anatomical lines to represent the anatomical structures
# Input:
#       - config: configuration object
# Output:
#       - get the position where in the matrix (eg. [[X-Sella, Y-Sella], [X-, Y- ]]) the coordinates of a anatomical line are
get_anatomical_lines <- function(config) {
  # Get the positions of the selected landmarks
  anatomical_lines <- lapply(config$anatomical_lines, find_positions, config$selected_landmark)
  return(anatomical_lines)
}

# get_anatomy_part_of_landmarks function: get the anatomical line of anatomical landmark 
# Input:
#       - config: configuration object
#       - coordinates: coordinates of the landmarks (list of names)
# Output:
#       - get the anatomical line (eg. cranial base, maxilla, mandible) of the coordinates
get_anatomy_part_of_landmarks <- function(config, coordinates) {
  # Get Landmarks and corresponding anatomical lines
  land_anatomy <- reverse_list(config$anatomical_lines)

  coord_anatomy <- list()
  # Get the anatomical lines of the coordinates
  for (coordinate in coordinates) {
    # Remove Prefix from the coordinates
    land <- gsub("^(X_|Y_)", "", coordinate)
    # Get the anatomical line of the landmark
    coord_anatomy[[coordinate]] <- land_anatomy[[land]]
  }
  return(coord_anatomy)
}

# reverse_list function: takes a named list of groups with correspondence to points and returns a named list of points with correspondence to groups
# This is usefull to get the anatomical line of coordinates and the category of cephalometric variables
reverse_list <- function(points_list) {
  reversed_correspondence <- list()
  # Loop through the list to reverse the correspondence
  for (group in names(points_list)) {
    for (point in points_list[[group]]) {
      reversed_correspondence[[point]] <- group
    }
  }
  return(reversed_correspondence)
}



# plot_mshapes function: plot multiple shapes
# Input:
#     - shapes: list of shapes to plot
#     - color: list of colors for each mean shape
#     - mshape_labels: list of labels for each mean shape
#     - joinlines: list of lines to join in the plot
#     - line_width: width of the lines (default = 1)
#     - plot_title: title of the plot (default = NULL)
#     - width: width of the plot (default = 500)
#     - height: height of the plot (default = 500)
# Output:
#     - plotly plot
plot_mshapes <- function(shapes, color, shape_labels, joinlines, line_width = 1, plot_title = NULL, width = 500, height = 500){
  # Initialize plot
  p <- plot_ly(width = width, height = height)
  # Determine x and y axis range
  all_coords <- do.call(rbind, shapes)
  x_y_range <- range(all_coords) + c(-10, 10)


  for (f in seq_along(shapes)){
    for (i in seq_along(joinlines)){
      showlegend = FALSE #if (i==1){showlegend=TRUE} else{showlegend=FALSE}
      mshape <- shapes[[f]]
      p <- add_trace(p,
        x = mshape[joinlines[[i]], 1],
        y = mshape[joinlines[[i]], 2],
        type = "scatter",
        mode = "lines + markers",
        marker = list(color = color[f], size = line_width * 2, symbol = "cross"),
        line = list(color = color[f], width = line_width),
        showlegend = showlegend,
        name = shape_labels[f])
        }
  }
  p <- layout(p,
            font = list(size = 19, family = "Helvetica"),
            title = list(text = plot_title, font = list(size = 40, family = "Helvetica"), x = 0.5, y=-0.5, xanchor = "center", yanchor = "center"),
            xaxis = list(range = x_y_range, title = "X", titlefont = list(size = 19, family = "Helvetica"), tickfont = list(size = 16, family = "Helvetica")),
            yaxis = list(range = x_y_range, title = "Y",  titlefont = list(size = 19, family = "Helvetica"), tickfont = list(size = 16, family = "Helvetica")))
  # Display the plot
  return(p)
}

# add_sd_polygons function: add standard deviation polygons to the plot
# Input:
#       - plot: plotly plot
#       - rotated coordinates: array of rotated coordinates 
#       - fill_color: color of the polygons (default = "rgba(128, 18, 12, 0.3)")
# Output:
#       - plotly plot with polygons
add_sd_polygons <- function(plot, rotated_array, fill_color = "rgba(128, 18, 12, 0.3)") {
  # Compute mean shape for the selected cluster
  mean_shape <- apply(rotated_array, c(1, 2), mean)
  # Compute standard deviation across the cluster shapes
  rotated_sd <- apply(rotated_array, c(1, 2), sd)
  
  # Construct polygons for each landmark
  polygon_x <- cbind(mean_shape[, 1] - rotated_sd[, 1],
                     mean_shape[, 1] + rotated_sd[, 1],
                     mean_shape[, 1] + rotated_sd[, 1],
                     mean_shape[, 1] - rotated_sd[, 1])
  
  polygon_y <- cbind(mean_shape[, 2] - rotated_sd[, 2],
                     mean_shape[, 2] - rotated_sd[, 2],
                     mean_shape[, 2] + rotated_sd[, 2],
                     mean_shape[, 2] + rotated_sd[, 2])
  
  # Add each polygon to the plot
  for (r in 1:nrow(mean_shape)) {
    plot <- plot %>% 
      plotly::add_polygons(
        x = polygon_x[r, ],
        y = polygon_y[r, ],
        type = "scatter",
        mode = "lines",
        fill = "toself",
        fillcolor = fill_color,
        line = list(color = "transparent"),
        showlegend = FALSE
      )
  }
  return(plot)
}



# compute_distance_centroid function: compute the distance of each sample to each corresponding cluster centroid
# Input:
#       - res_cluster: eclust object with values; $cluster; $data; $centers
# Output:
#       - vector with the distance to centroid per sample
compute_distance_centroid <- function(cluster_result) {
  distances <- c()
  cluster_assignments <- cluster_result$cluster

  for (i in seq_along(cluster_assignments)) {
    cluster <- cluster_assignments[i]
    point <- cluster_result$data[i, ]
    centroid <- cluster_result$centers[cluster, ]
    euclidean_distance <- sqrt(sum((point - centroid)^2))
    distances <- c(distances, euclidean_distance)
  }

  return(distances)
}

# create_cv_folds function: create cross-validation folds
# Input:
#       - cv_stratify: stratify the sampling based on a variable
#       - set_seed: set the random seed for reproducibility
#       - n_folds: number of folds
#       - df: dataframe to sample from
#       - verbose: print fold distribution for checking
# Output:
#       - vector with the fold number for each sample
create_cv_folds <- function(cv_stratify = NULL, set_seed = 123, n_folds = 5, df = NULL, verbose = FALSE) {
    # Set the random seed for reproducibility
    set.seed(set_seed)
    
    if (!is.null(cv_stratify)) {
        # Stratified sampling using createFolds
        folds <- createFolds(cv_stratify, k = n_folds, list = FALSE)
        if (verbose) {
            print("--Sampling--")
            # Print the first 5 fold counts for checking
            print(table(folds))
        }
    } else {
        # Random sampling if cv_stratify is NULL
        folds <- sample(1:n_folds, nrow(df), replace = TRUE)
        if (verbose) {
            print("--Sampling--")
            # Print fold distribution for checking
            print(table(folds))
        }
    }
    
    return(folds)
}

# cluster_matching function: match clusters from two different clusterings using the Hungarian algorithm
# Input:
#       - cluster_from: vector of cluster assignments
#       - cluster_to: vector of cluster assignments
#       - verbose: print the matching results if TRUE
# Output:
#       - list with the matching results
cluster_matching <- function(cluster_from, cluster_to, verbose = FALSE) {
  # Ensure input vectors have names for matching
  if (is.null(names(cluster_from)) || is.null(names(cluster_to))) {
    stop("Both cluster vectors must have names corresponding to sample identifiers.")
  }

  # Merge clusters based on names
  common_names <- intersect(names(cluster_from), names(cluster_to))
  if (length(common_names) == 0) {
    stop("No matching sample names found between the two cluster vectors.")
  }

  merged_clusters <- data.frame(
    cluster_1 = cluster_from[common_names],
    cluster_2 = cluster_to[common_names]
  )

  # Create contingency table (cross-tabulation)
  cont_table <- table(merged_clusters$cluster_1, merged_clusters$cluster_2)

  # Convert contingency table to a matrix
  cost_matrix <- as.matrix(cont_table)

  # Use Hungarian algorithm (solve_LSAP) to find the optimal matching
  assignment <- solve_LSAP(cost_matrix, maximum = TRUE)  # Maximizing the match

  # Calculate the matching score
  matching_score <- sum(cost_matrix[cbind(assignment, 1:length(assignment))]) / sum(cost_matrix)

  # Print the results if verbose is TRUE
  if (verbose) {
    cat("Matching from cluster 1 to cluster 2:\n")
    print(assignment)
    cat("\nMatching score: ", matching_score, "\n")
  }

  # Return results
  return(list(
    map = assignment,
    cont_table = cont_table,
    matching_score = matching_score))
}

# cluster_mapping function: map clusters from one clustering to another using the matching results
# Input:
#       - clusters: vector of cluster assignments
#       - map: matching results from cluster_matching function
# Output:
#       - vector of mapped cluster assignments
cluster_mapping <- function(clusters, map) {
  # Use map to map clusters in 'clusters' from cluster2 to cluster1
  transformed_clusters <- sapply(clusters, function(x) map[x])

  # Preserve the original names of the vector (indexing is kept intact)
  names(transformed_clusters) <- names(clusters)
  return(transformed_clusters)
}