library(plotly)
library(grid)
library(clue) # for cluster matching (solve lsap)

# euclidean_distance function: compute the euclidean distance between two points
# Input:
#       - x: first point
#       - y: second point
# Output:
#       - euclidean distance
euclidean_distance <- function(x, y) {
  sqrt(sum((x - y)^2, na.rm = TRUE))
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
