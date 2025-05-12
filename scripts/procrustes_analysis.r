library(shapes)

# proc_gpa function: performs Generalized Procrustes Analysis (GPA) on landmark coordinates dataframe
# Input:    - df_coordinates - dataframe with landmark coordinates
# Output:   - proc_gpa - list with the results of the GPA
generalized_procrustes <- function(df_coordinates) {
  patient_labels <- rownames(df_coordinates)
  coordinate_labels <- colnames(df_coordinates)

  ##### transform into an array of matrices -> coordinates_data
  list_dimnames <- list(unique(gsub("^X_|^Y_", "", colnames(df_coordinates))), c("x", "y"), patient_labels)
  matrix_coordinates <- df_to_array_matrix(df_coordinates, list_dimnames)
 
  ###### generalized Procrustes Analysis
  proc_gpa <- procGPA(matrix_coordinates)
  # defining labels for the output of the GPA
  # coordinates_labels are arranged in x1, x2, x3, ..., y1, y2, y3, ...
  # $pcar: eigenvectors of the PCA having the landmarks as rows
  # $tan: tangent space coordinates of the GPA, which are the residuals
  row.names(proc_gpa$pcar) <- c(coordinate_labels[seq(1, length(coordinate_labels), 2)], coordinate_labels[seq(2, length(coordinate_labels), 2)])# coordinates_labels
  row.names(proc_gpa$tan) <- c(coordinate_labels[seq(1, length(coordinate_labels), 2)], coordinate_labels[seq(2, length(coordinate_labels), 2)])# coordinates_labels
  colnames(proc_gpa$tan) <- patient_labels

  row.names(proc_gpa$stdscores) <- patient_labels
  row.names(proc_gpa$rawscores) <- patient_labels
  row.names(proc_gpa$scores) <- patient_labels

  proc_gpa$matrix_coordinates <- matrix_coordinates

  proc_gpa$tan <- t(proc_gpa$tan) # transpose the tangent space coordinates so that it is in x1, x2, y1, y2, ... format (columns) and patients as rows
  ######### save rotated coordinates into dataframe
  proc_gpa$df_rotated <- array_matrix_to_df(proc_gpa$rotated, coordinate_labels, patient_labels)

  return(proc_gpa)
}

# df_to_array_matrix function: transforms a dataframe with landmark coordinates into an array of matrices [,,]
# Input:    - df - dataframe with landmark coordinates
#           - dimnames - names of the dimensions of the array (type: list)
# Output:   - matrix - array of matrices [,,]
df_to_array_matrix <- function(df, dimnames = NULL) {
  ##### transform into an array of matrices -> coordinates_data
  colnames(df) <- NULL
  dimensions <- dim(df)
  k <- dimensions[2] / 2 #nb landmarks
  m <- 2 #2D coordinates (x,y)
  n <- dimensions[1] #nb of patients
  matrix <- array(0, c(k, m, n))
  matrix_data <- as.matrix(df)
  for (i in 1:n){
    matrix[, , i] <- matrix(matrix_data[i, ], ncol = 2, byrow = TRUE)
  }
  if (!is.null(dimnames)) {
    dimnames(matrix) <- dimnames
  }
  return(matrix)
}

# array_matrix_to_df function: transforms an array of matrices into a dataframe
# Input:    - array_matrix - array of matrices [,,]
#           - colnames - column names of the dataframe
#           - rownames - row names of the dataframe
# Output:   - df - dataframe
array_matrix_to_df <- function(array_matrix, colnames, rownames) {
  matrix_rotated_coordinates <- as.vector(t(array_matrix[, , 1]))
  # Loop over matrices and add rows to the matrix
  for (n in 2:dim(array_matrix)[3]){
    matrix_rotated_coordinates <- rbind(matrix_rotated_coordinates, as.vector(t(array_matrix[, , n]))) # x1, y1, x2, y2, ...
  }
  colnames(matrix_rotated_coordinates) <- colnames
  rownames(matrix_rotated_coordinates) <- rownames
  df <- data.frame(matrix_rotated_coordinates)
  return(df)
}


# shape_ordinary_procrustes function: performs Ordinary Procrustes Analysis (OPA) and computes the PCA scores for a new data point (after computing the tangent space coordinates)
# Input:    - reference_shape - the reference shape (mean shape of the cluster) in the form of a matrix
#           - eigenvectors - eigenvectors of the PCA
#           - shape - the shape of the new data point
# Output:   - opa - list with the results of the OPA, TAN, and the PCA scores
shape_ordinary_procrustes <- function(reference_shape, eigenvectors = NULL, shape) {
  # perform OPA
  opa <- procOPA(reference_shape, shape, scale = TRUE, reflect = FALSE)












  # Save the rotated shape (matrix)
  opa$rotated <- opa$Bhat

  # Compute tangent space coordinates (which are the residuals)
  opa$tan <- opa$Bhat - reference_shape

  # Linearize rotated shape
  opa$df_rotated <- as.vector(t(opa$Bhat)) # df_rotated should be in x1, y1, x2, y2, ... format

  opa$tan <- as.vector(opa$tan) # linearize the tangent space coordinates so that it is in x1, x2, y1, y2, ... format
  # Compute the pca scores for the new data point
  # the eigenvectors are orded by x1, x2, x3, ..., y1, y2, y3, ...
  if (!is.null(eigenvectors)) {
    opa$scores <- as.vector(as.vector(opa$tan) %*% eigenvectors)
  }

  return(opa)
}

# ordinary_procrustes function: performs the Ordinary Procrustes Analysis (OPA) for multiple patients (given an input a dataframe with the coordinates of the landmarks) and the reference shape
# Input:    - df_coordinates - dataframe with the coordinates of the landmarks (in the format x1, y1, x2, y2, ...)
#           - eigenvectors - eigenvectors of the PCA
#           - reference_shape - the reference shape (mean shape of the cluster) in the form of a matrix
#           - n_pca - number of PCA components to compute
# Output:   - opa - list with the results of the OPA, TAN, and the PCA scores
ordinary_procrustes <- function(df_coordinates, eigenvectors, reference_shape, n_pca = NULL) {
  # Transform the dataframe into an array of matrices
  array_matrix_coordinates_test <- df_to_array_matrix(df_coordinates)

  opa_list <- list()
  rotated = array(NA, dim = c(ncol(df_coordinates)/2, 2, nrow(df_coordinates))) # considering this is a 2D problem
  # Per sample perform the OPA
  for (idx in seq_len(nrow(df_coordinates))) {
      # Perform ordinary procrustes analysis, compute the residuals and transform the coordinates into the training PC space
      shape_opa <- shape_ordinary_procrustes(reference_shape = reference_shape,
                                       eigenvectors = eigenvectors,
                                       shape = array_matrix_coordinates_test[, , idx])
      # Store each result type in the opa_list
      rotated[, ,idx] <- shape_opa[["rotated"]]
      opa_list$tan[[idx]] <- shape_opa[["tan"]]
      opa_list$scores[[idx]] <- shape_opa[["scores"]]
      opa_list$df_rotated[[idx]] <- shape_opa[["df_rotated"]]
      if (!is.null(n_pca)) {
        opa_list$scores_slice[[idx]] <- shape_opa$scores[1:n_pca]
      }
  }
  opa_list$rotated <- rotated
  dimnames(opa_list$rotated) <- list(NULL, c("x", "y"), rownames(df_coordinates))

  # Convert the lists into data frames,include the colnames and rownames
  opa_list$tan <- data.frame(do.call(rbind, opa_list$tan), row.names = rownames(df_coordinates))
  opa_list$scores <- data.frame(do.call(rbind, opa_list$scores), row.names = rownames(df_coordinates))
  opa_list$df_rotated <- data.frame(do.call(rbind, opa_list$df_rotated), row.names = rownames(df_coordinates))
  if (!is.null(n_pca)) {
    opa_list$scores_slice <- data.frame(do.call(rbind, opa_list$scores_slice), row.names = rownames(df_coordinates))
    colnames(opa_list$scores_slice) <- paste0("PC", seq_len(ncol(opa_list$scores_slice)))
  }
  colnames(opa_list$tan) <- c(colnames(df_coordinates)[seq(1, ncol(df_coordinates), 2)], 
                              colnames(df_coordinates)[seq(2, ncol(df_coordinates), 2)])  # x1, x2, x3, ..., y1, y2, y3, ...
  colnames(opa_list$scores) <- paste0("PC", seq_len(ncol(opa_list$scores)))
  colnames(opa_list$df_rotated) <- colnames(df_coordinates)

  return(opa_list)
}