# Import Libraries
library(dplyr)
library(tidyr)
library(ggplot2)

library(caret) # for createFolds
library(gtools) # for permutations
library(class) # for the knn function
library(kknn) # for the kknn function
library(stringr) # for the str_replace
library(factoextra) # for the eclust function

# Load script
source("scripts/utils.R")
source("scripts/procrustes_analysis_utils.R")
set.seed(123)

# create_reference_clustering function: create a reference clustering using the coordinates data
# Input:  
#         - df_coordinates: data frame with coordinates
#         - fun_cluster: clustering function to use; "kmean" or "fuzzy"
#         - hc_metric: distance metric for hierarchical clustering
#         - hc_method: method for hierarchical clustering
#         - fuzzy_membership: membership parameter for fuzzy clustering
#         - k: number of clusters
#         - ground_truth_metric: for clustering (except for fuzzy) ground_truth_metric can be  "distance_to_centroid" or "silhouette_width"
#         - type: type of data to use for clustering; "tan" (residuals), "scores" (PCA scores on residuals data), "df_rotated" (coordinates after Procrustes Analysis), "scores_slice" (selected PCA scores on residuals data)
#         - n_pca: number of PCA components to use for clustering used in case of type = "scores_slice"
# Output:
#         - ground_truth: list with the ground truth clustering
create_reference_clustering <- function(df_coordinates,
                                        fun_cluster = "kmean",
                                        hc_metric = "euclidean",
                                        hc_method = NULL,
                                        fuzzy_membership = 1.2,
                                        k = 6,
                                        ground_truth_metric = "distance_centroid",
                                        type = "tan",  # "tan", "scores", "df_rotated", "scores_slice" -> type of data to use for clustering
                                        n_pca = 10) {
  # Check the type of data to use for clustering is valid
  stopifnot(type %in% c("tan", "scores", "df_rotated", "scores_slice"))

  # Reference clustering with all data
  gpa <- generalized_procrustes(df_coordinates)  # perform GPA on the data (assuming generalized_procrustes is a function from utils.R)
  
  if (type == "scores_slice") {
    gpa$scores_slice <- gpa$scores[, 1:n_pca]  # select the first n_pca scores
  }

  # Perform clustering (fuzzy or not) on the reference data
  if (fun_cluster == "fuzzy") {
    fuzzy_ground_truth <- FKM(X = gpa[[type]], m = fuzzy_membership, RS = 123, stand = 0, k = k)
    ground_truth <- list()
    ground_truth$cluster <- fuzzy_ground_truth$clus[, 1]  # extract the clusters from fuzzy output
    ground_truth$metric <- fuzzy_ground_truth$clus[, 2]  # extract fuzzy membership metric
  } else {
    ground_truth <- eclust(gpa[[type]], FUNcluster = fun_cluster, hc_metric = hc_metric, hc_method = hc_method, k = k, graph = FALSE)
    
    if (ground_truth_metric == "silhouette_width") {
      ground_truth$metric <- ground_truth$silinfo$widths$sil_width
    } else if (ground_truth_metric == "distance_centroid") {
      # Compute the distance to the centroid as ground truth metric
      ground_truth$metric <- compute_distance_centroid(ground_truth) # utils.R function
    }
  }
  ground_truth$mshape <- gpa$mshape  # Save the mean shape
  ground_truth$type <- gpa[[type]] # Save the data used for clustering

  # Return the ground truth clustering
  return(ground_truth)
}

# get_training_clusters function: get the training clusters
# This function returns the training clusters obtained by performing GPA on the training data and then clustering the data
# Input:
#         - df: data frame with coordinates
#         - type: type of data to use for clustering; "tan" (residuals), "scores" (PCA scores on residuals data), "df_rotated" (coordinates after Procrustes Analysis), "scores_slice" (selected PCA scores on residuals data)
#         - set_seed: seed for reproducibility
#         - k: number of clusters
#         - fun_cluster: clustering function to use
#         - hc_metric: distance metric for hierarchical clustering
#         - hc_method: method for hierarchical clustering
#         - n_pca: number of PCA components to use for clustering used in case of type = "scores_slice"
#         - fuzzy_membership: membership parameter for fuzzy clustering
get_training_clusters <- function(df, type, set_seed, k, fun_cluster, n_pca, fuzzy_membership,
                       hc_metric, hc_method) {

  training_results <- list()
  # Perform GPA on the training data
  training_results$gpa <- generalized_procrustes(df)
  if (type == "scores_slice") {
    gpa_training_data$gpa$scores_slice <- gpa_training_data$scores[, 1:n_pca]
  }
  # Get the training data
  training_results$train_data <- training_results$gpa[[type]]

  # Perform clustering
  if (fun_cluster == "fuzzy") {
    fuzzy_result <- FKM(X = training_results$train_data, m = fuzzy_membership, RS = set_seed, stand = 0, k = k)
    training_results$cluster <- fuzzy_result$clus[, 1] # Extract the hard cluster of fuzzy clustering
  } else {
    eclust_result <- eclust(training_results$train_data, FUNcluster = fun_cluster, hc_metric = hc_metric, 
                            hc_method = hc_method, k = k, graph = FALSE)
    training_results$cluster <- eclust_result$cluster
  }
  return(training_results)
}

# get_testing_clusters function: get the testing clusters: this function returns the testing clusters obtained by aligning the test data to the reference shape of the training data through 
# ordinary procrustes analysis and then assign a cluster through KNN
# Input:
#         - indices: indices of the testing data (optional) (this will filter the df)
#         - df: data frame with coordinates
#         - gpa: GPA object from the training data
#         - type: type of data to use for clustering; "tan" (residuals), "scores" (PCA scores on residuals data), "df_rotated" (coordinates after Procrustes Analysis), "scores_slice" (selected PCA scores on residuals data)
#         - set_seed: seed for reproducibility
#         - knn_params: list with the parameters for KNN (parameters for the kknn function) -> examples: k, distance, kerner, scale
#         - n_pca: number of PCA components to use for clustering used in case of type = "scores_slice"
#         - X_train: training data (used for KNN)
#         - y_train: training labels (used for KNN)
# Output:
#         - testing_results: list with the results of the testing data, including the result from ordinary procrustes analysis and the cluster assignment
get_testing_clusters <- function(indices = NULL, df, gpa, type = "tan", set_seed = 123, knn_params = list(
                                                                                                    k = 25,
                                                                                                    distance = 2,
                                                                                                    kernel = "rectangular",
                                                                                                    scale = FALSE
                                                                                                  ), n_pca = NULL, X_train, y_train) {
  # Default parameters for KNN
  default_knn_params <- list(
  k = 25,
  distance = 2,
  kernel = "rectangular",
  scale = FALSE
)

  testing_results <- list()

  # Filter the data if indices are provided
  if (!is.null(indices)) {
    df <- df[indices, ]
  }

  # Perform OPA on all testing data (from procrustes_analysis.r)
  testing_results$opa <- ordinary_procrustes(df_coordinates = df, eigenvectors = gpa$pcar, reference_shape = gpa$mshape, n_pca = n_pca)
  testing_results$test_data <- testing_results$opa[[type]] # Get the data used for clustering

  # Perform knn on the testing data - to assign the cluster of new data points
  # Knn is trained with the X_train and y_train
  #testing_results$cluster <- knn(train = X_train, test = testing_results$test_data, cl = y_train, k = n_nearst_neighbours)

  df_train <- as.data.frame(X_train)
  df_train$Cluster <- factor(y_train)
  testing_results$test_data <- as.data.frame(testing_results$test_data)

  knn_params <- modifyList(default_knn_params, knn_params)
  knn_input <- modifyList(knn_params, list(
    train = df_train,
    test = testing_results$test_data,
    formula = formula(Cluster ~ .)
  ))
  
  model_kknn <- do.call(kknn::kknn, knn_input)
  testing_results$cluster  <- fitted(model_kknn)
  testing_results$probability <- model_kknn$prob
  rownames(testing_results$probability) <- rownames(testing_results$test_data)
  return(testing_results)
}

# print_cross_validation_summary function: print the summary of the cross-validation results
# Input: cross_validation_results (list with the results of the cross-validation)
# Output: print the summary of the cross-validation results (void)
print_cross_validation_summary <- function(cross_validation_results) {
  # Extract values for training and testing accuracies
  mean_training <- round(cross_validation_results$mean_training_accuracy, 2)
  std_training <- round(cross_validation_results$std_training_accuracy, 2)
  mean_testing <- round(cross_validation_results$mean_testing_accuracy, 2)
  std_testing <- round(cross_validation_results$std_testing_accuracy, 2)

  # Format fold accuracies
  fold_training_accuracies <- paste(sapply(cross_validation_results$fold_accuracies_training, function(x) round(x, 2)), collapse = '; ')
  fold_testing_accuracies <- paste(sapply(cross_validation_results$fold_accuracies_testing, function(x) round(x, 2)), collapse = '; ')

  # Print results
  print("Cross-Validation Summary:")
  cat(sprintf("Training: (%.2f ± %.2f)%%\n", mean_training, std_training))
  cat("Fold Training Accuracies:", fold_training_accuracies, "(in %)\n")
  cat(sprintf("Testing: (%.2f ± %.2f)%%\n", mean_testing, std_testing))
  cat("Fold Testing Accuracies:", fold_testing_accuracies, "(in %)\n")
}


# clustering_cross_validation function: cross-validation for clustering using coordinates data
# Input:
#   df_coordinates: data frame with coordinates
#   cv_stratify: vector with labels to stratify the cross-validation
#   k: number of clusters
#   fun_cluster: clustering function to use
#   hc_metric: distance metric for hierarchical clustering
#   hc_method: method for hierarchical clustering
#   n_folds: number of folds
#   knn_params: list with the parameters for KNN (parameters for the kknn function) -> examples: k, distance, kerner, scale
#   type: type of data to use for clustering; "tan" (residuals), "scores" (PCA scores on residuals data), "df_rotated" (coordinates after Procrustes Analysis), "scores_slice" (selected PCA scores on residuals data)
#   n_pca: number of PCA components to use for clustering used in case of type = "scores_slice"
#   ground_truth_metric: for clustering (except for fuzzy) ground_truth_metric can be  "distance_to_centroid" or "silhouette_width"
#   verbose: print messages
#   summarize: summarize the results
#   set_seed: seed for reproducibility
# Output:
#   cross_validation_results: list with the results of the cross-validation
clustering_cross_validation <- function(df_coordinates,
                            cv_stratify = NULL,
                            k = 6,
                            fun_cluster = "kmean",
                            hc_metric = "euclidean",
                            hc_method = NULL,
                            fuzzy_membership = 1.2,
                            n_folds = 5,
                            knn_params = list(k = 25,
                                              distance = 2,
                                              kernel = "rectangular",
                                              scale = FALSE),
                            type = "tan", # "tan", "scores", "df_rotated", "scores_slice" -> type of data to use for clustering
                            n_pca = 10,
                            ground_truth_metric = "distance_centroid",
                            verbose = FALSE,
                            summarize = TRUE,
                            set_seed = 123) {
    # Check the type of data to use for clustering is valid                              
    stopifnot(type %in% c("tan", "scores", "df_rotated", "scores_slice"))

    if (type == "scores_slice") {
    stopifnot(!is.null(n_pca), n_pca <= ncol(df_coordinates), n_pca > 0)
    }

    patient_labels <- row.names(df_coordinates)

    # Create list - output object
    cross_validation_results <- list()
    cross_validation_results$fold_accuracies_testing <- vector("list", n_folds)
    cross_validation_results$fold_accuracies_training <- vector("list", n_folds)

    # Create cross_validation results matrix
    # Col fold_{n} with n the number of the fold -> "train" or "test"
    # Col partition_fold_{n} with n the number of the fold -> predicted cluster
    matrix_cross_validation <- matrix(nrow = length(patient_labels), ncol = (2 + n_folds + n_folds),
                                dimnames = list(patient_labels, c("ground_true", "ground_truth_metric", paste("fold", 1:n_folds, sep = "_"), paste("partition_fold",  1:n_folds, sep = "_")))
                                )

    ##################################
    # reference/ground_truth clustering with all data
    ground_truth <- create_reference_clustering(df_coordinates = df_coordinates,
                                              fun_cluster = fun_cluster,
                                              fuzzy_membership = fuzzy_membership,
                                              hc_metric = hc_metric,
                                              hc_method = hc_method,
                                              k = k,
                                              type = type,
                                              n_pca = n_pca,
                                              ground_truth_metric = ground_truth_metric)
   
    cross_validation_results$ground_truth <- ground_truth$cluster
    cross_validation_results$reference_data <- ground_truth[[type]]
    cross_validation_results$reference_mshape <- ground_truth$mshape
    matrix_cross_validation[patient_labels, "ground_true"] <- ground_truth$cluster
    matrix_cross_validation[patient_labels, "ground_truth_metric"] <- ground_truth$metric

    ##################################
    # Sampling for cross-validation
    folds <- create_cv_folds(cv_stratify = cv_stratify, set_seed = set_seed, n_folds = n_folds, df = df_coordinates, verbose = verbose) # function from utils.R

    ##################################
    # Cross-validation - loop over folds
    for (fold in 1:n_folds) {
      if (verbose) {print(paste("-----Fold ", fold, "-----"))}

      # Get training clusters
      if (verbose) {print("-> Training")}
      train_indices <- which(folds != fold)

      # Perform GPA and clustering on the training data
      train_results <- get_training_clusters(df = df_coordinates[train_indices, ], n_pca = n_pca, type = type,
                          fuzzy_membership = fuzzy_membership, set_seed = set_seed, k = k, fun_cluster = fun_cluster,
                          hc_metric = hc_metric, hc_method = hc_method)

      # Compute the matching between the training clusters and the ground truth clusters (with functions from utils.R)
      cluster_map <- cluster_matching(cluster_to = ground_truth$cluster[train_indices], cluster_from = train_results$cluster, verbose = FALSE)

      y_train <- train_results$cluster
      # Reassign the training clusters based on the matching
      train_results$cluster <- cluster_mapping(train_results$cluster, cluster_map$map)  # Map the training clusters to the ground truth clusters
 
      train_acc <- sum(ground_truth$cluster[train_indices] == train_results$cluster) * 100 / length(train_indices)

      if (verbose) {
        print(paste("Training Accuracy: ", train_acc, "%"))
      }

      # Get testing clusters
      if (verbose) {print("-> Testing")}
      test_indices <- which(folds == fold)

      test_results <- get_testing_clusters(indices = test_indices, df = df_coordinates, gpa = train_results$gpa, type = type,
                        set_seed = set_seed, knn_params = knn_params, n_pca = n_pca, X_train = train_results$train_data,
                        y_train = y_train)
      # Reassign the testing clusters based on the matching
      test_results$cluster <- cluster_mapping(test_results$cluster, cluster_map$map)  # Map the testing clusters to the ground truth clusters

      # Compute Accuracy
      test_acc <- sum(ground_truth$cluster[test_indices] == test_results$cluster) * 100 / length(test_indices)
      if (verbose) {
        print(paste("Testing Accuracy: ", test_acc, "%"))
      }


      # Save fold results
      matrix_cross_validation[patient_labels[train_indices], paste("fold", fold, sep = "_")] <- "train"
      matrix_cross_validation[patient_labels[test_indices], paste("fold", fold, sep = "_")] <- "test"
      
      matrix_cross_validation[patient_labels[train_indices], paste('partition_fold', fold, sep = '_')] <- train_results$cluster
      matrix_cross_validation[patient_labels[test_indices], paste('partition_fold', fold, sep = '_')] <- test_results$cluster

      cross_validation_results$fold_accuracies_training[[fold]] <- train_acc
      cross_validation_results$fold_accuracies_testing[[fold]] <- test_acc
    }
    # Save the cross-validation results
    cross_validation_results$fold_accuracies_training <- unlist(cross_validation_results$fold_accuracies_training)
    cross_validation_results$fold_accuracies_testing <- unlist(cross_validation_results$fold_accuracies_testing)
    # Compute overall accuracy mean and std
    cross_validation_results$mean_training_accuracy <- mean(cross_validation_results$fold_accuracies_training)
    cross_validation_results$mean_testing_accuracy <- mean(cross_validation_results$fold_accuracies_testing)
    cross_validation_results$std_training_accuracy <- sd(cross_validation_results$fold_accuracies_training)
    cross_validation_results$std_testing_accuracy <- sd(cross_validation_results$fold_accuracies_testing)
    cross_validation_results$df_cross_validation <- data.frame(matrix_cross_validation)

    if (summarize) {
      # Print the cross-validation summary
      print_cross_validation_summary(cross_validation_results)
    }

  return (cross_validation_results)
}

# violin_plot_training_metric function: plot violin plot of the ground truth metric values for the training data (correct and incorrect assignments)
# Input: df_cross_validation (data frame with the cross-validation results, ground_true, ground_truth_metric, partition_fold_n, fold_n)
#        which_folds (list of folds to plot)
#        name_metric (name of the metric to plot)
#        wilcox_alternative (alternative hypothesis for the Wilcoxon test) - default is "less"; other options are "greater" and "two.sided"
#        output_path (type: string) - path to save the plot (optional)
# Output: ggplot object
violin_plot_training_metric <- function(df_cross_validation, which_folds, name_metric = "ground_truth_metric", wilcox_alternative = "less", output_path = NULL) {
  # Data Preparation - df_cross_validation is a dataframe with ground_true,	ground_truth_metric,	fold_,	partition_fold_
  df_long <- df_cross_validation %>%
    pivot_longer(cols = starts_with("fold_"), names_to = "fold", values_to = "train_test") %>%
    pivot_longer(cols = starts_with("partition_fold_"), names_to = "partition_fold", values_to = "partition") %>%
    filter(str_replace(fold, "fold_", "") == str_replace(partition_fold, "partition_fold_", "")) %>%  # Compare the same fold
    filter(fold %in% paste0("fold_", which_folds)) %>%  # Filter by the selected folds (through the which_folds argument)
    mutate(fold = paste0("Fold ", str_replace(fold, "fold_", ""))) %>%  # Change to "Fold X"
    filter(train_test == "train") %>%  # Filter only the training data
    mutate(correct = as.factor(ifelse(partition == ground_true, "yes", "no")),  # Create a variable to distinguish correct and incorrect assignments in the training step of the CV
           # Convert the column referenced by name_metric to numeric
            !!name_metric := as.numeric(.data[[name_metric]]))

   # Wilcoxon test per fold and prepare annotation data
    wilcox_results <- df_long %>%
        group_by(fold) %>%
        summarise(
            # Wilcox test to compare the ground truth metric values for correct and incorrect assignments
            p_value = wilcox.test(
              .data[[name_metric]][correct == "yes"],
              .data[[name_metric]][correct == "no"],
                alternative = wilcox_alternative
            )$p.value
        ) %>%
        mutate(p_label = paste0("p = 10^", round(log10(p_value), 1)))

     df_long <- df_long %>% left_join(wilcox_results, by = "fold") # Join the wilcox_results to the main dataframe

    # Create violin plot with the ground truth metric values for the training data and the correct and incorrect assignments
     plot <- ggplot(df_long, aes(x = correct, y = .data[[name_metric]], fill = correct)) +
        geom_violin(color = "black", alpha = 0.7) +
        scale_fill_manual(values = c("red", "blue")) +
        labs(
            title = "",
            x = "Correct Assignment",
            y = name_metric
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(size = 0),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            strip.text = element_text(size = 20)  # Facet label text size
        ) +

        facet_wrap(~ fold, scales = "free", ncol = 2) +
        geom_text(
            data = wilcox_results,
            aes(x = 1.5, y = Inf, label = p_label),
            hjust = 0.5, vjust = 1.5,
            size = 5, color = "black",
            inherit.aes = FALSE
        )

      # Save plot to pdf if output_path is provided
      if (!is.null(output_path)) {
        ggsave(paste0(output_path, "violin_cv_distance_centroid.pdf"), plot = plot, device = "pdf", width = 20, height = 10)
      }
    return (plot)
}

# hist_fold_accuracy function: plot histogram of training and test accuracies per fold
# Input: training_accuracy (type : list), test_accuracy (type: list)
#        output_path (type: character) - path to save the plot (optional)
# Output: ggplot object
# Example 1 of usage: hist_fold_accuracy(training_accuracy = c(0.8, 0.9, 0.7), test_accuracy = c(0.7, 0.8, 0.6))
# Example 2 of usage: hist_fold_accuracy(res.cross_validation$fold_accuracies_training, res.cross_validation$fold_accuracies_testing)
hist_fold_accuracy <- function(training_accuracy, test_accuracy, output_path = NULL) {

  n_folds <- length(training_accuracy)
  list_folds <- 1:n_folds

  # create data frame with training data first, then test data
  data <- data.frame(
    Index = rep(list_folds, 2),
    Accuracy = c(training_accuracy, test_accuracy),
   Type = factor(rep(c("Training", "Testing"), each = n_folds), levels = c("Training", "Testing"))
  )
  
  # calculate mean and standard deviation
  mean_training <- mean(training_accuracy)
  sd_training <- sd(training_accuracy)
  mean_test <- mean(test_accuracy)
  sd_test <- sd(test_accuracy)
  
  # subtitle text 
  annotation_text <- paste(
    sprintf("Training: %.2f \u00b1 %.2f", mean_training, sd_training),
    sprintf("Testing: %.2f \u00b1 %.2f", mean_test, sd_test),
    sep = "\n"
  )
  
  # plotting
  plot <- ggplot(data, aes(x = as.factor(Index), y = Accuracy, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "Fold", y = "Accuracy", title = '', #paste0("Cross-Validation ", n_folds,"-folds - Accuracy"),
              subtitle = annotation_text) +
    scale_fill_manual(values = c("Training" = "#215F9B", "Testing" = "#DBB300")) +
    geom_text(
    aes(label = sprintf("%.2f", Accuracy)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4) +
    theme_classic() +  # Use classic theme to match the plain white background
    theme(
      plot.title = element_text(size = 18, face = "bold", family = "Helvetica"), 
      plot.subtitle = element_text(size = 18, family = "Helvetica"),
      axis.title.x = element_text(size = 16, family = "Helvetica"),
      axis.title.y = element_text(size = 16, family = "Helvetica"),
      axis.text.x = element_text(size = 14, family = "Helvetica"),
      axis.text.y = element_text(size = 14, family = "Helvetica"),
      legend.title = element_text(size = 14, family = "Helvetica", face = "bold"),
      legend.text = element_text(size = 14, family = "Helvetica")
    )

  # Save plot if output_path is provided
  if (!is.null(output_path)) {
    ggsave(paste0(output_path, "hist_fold_accuracy.pdf"), plot = plot, device = "pdf", width = 8, height = 8)
  }
  return (plot)
}

# perform_cross_validation_knn function: Perform cross-validation for different values of k nearest neighbors
# Input:
#       - df: data frame with coordinates
perform_cross_validation_knn <- function(df, cv_stratify_factors = NULL, fun_cluster = 'kmean', k = 6, n_folds, knn_param_to_combine, type = 'tan') {

  # Check if the grid_param is a list
  if (!is.list(knn_param_to_combine)) {
    stop("knn_param_to_combine should be a list.")
  }
  # Check if is not null
  if (is.null(knn_param_to_combine)) {
    stop("knn_param_to_combine should not be NULL.")
  }

  # Include all default parameters for knn
  default_knn_params <- list(
  k = 25,
  distance = 2,
  kernel = "rectangular",
  scale = FALSE
 )
 knn_param_to_combine <- modifyList(default_knn_params, knn_param_to_combine)
 

  # Create a dt with the parameters to combine
  grid_param <- expand.grid(knn_param_to_combine, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)   
  df_result <- data.frame()

  for (i in seq_len(nrow(grid_param))) {
    current_params <- as.list(grid_param[i, ])

    # Perform cross-validation (with a function from cross_validation.r)
    res <- clustering_cross_validation(
      df_coordinates = df,
      fun_cluster = fun_cluster,
      k = k,
      n_folds = n_folds,
      verbose = FALSE,
      summarize = FALSE,
      cv_stratify = cv_stratify_factors,
      knn_params = current_params,
      type = type
    )
    
    # Save row with results
    # Save results + hyperparameters
    fold_test_cols <- setNames(as.list(res$fold_accuracies_testing), 
                           paste0("fold_test_", seq_along(res$fold_accuracies_testing)))
    fold_train_cols <- setNames(as.list(res$fold_accuracies_training),
                            paste0("fold_train_", seq_along(res$fold_accuracies_training)))
    knn_param_cols <- as.vector(grid_param[i, ])
                      

    df_result <- rbind(df_result, data.frame(
      mean_training = res$mean_training_accuracy,
      sd_training = res$std_training_accuracy,
      mean_testing = res$mean_testing_accuracy,
      sd_testing = res$std_testing_accuracy,
      fold_test_cols,
      fold_train_cols,
      knn_param_cols,
      stringsAsFactors = FALSE
    ))}
  return(df_result)
}

# format_cross_validation_results function: format the cross-validation results
# Input:
#       - df_cv_knn: data frame with cross-validation results (obtained from perform_cross_validation_knn function or equivalent)
#       - output_path: path to save the formatted results
# Output:
#       - data frame with formatted cross-validation results
format_cross_validation_results <- function(df_cv_knn, output_path = NULL) {
  df_cv_knn_formatted <- data.frame()
  for (r in seq_len(nrow(df_cv_knn))) {
    Training <- paste0(round(df_cv_knn$mean_training[r], 2), " ± ", round(df_cv_knn$sd_training[r], 2))
    Testing <- paste0(round(df_cv_knn$mean_testing[r], 2), " ± ", round(df_cv_knn$sd_testing[r], 2))
    df_cv_knn_formatted <- rbind(df_cv_knn_formatted, data.frame(K_KNN = df_cv_knn$k[r], Training = Training, Testing = Testing))
  }
  
  # Save formatted results to a file
  if (!is.null(output_path)) {
    write.csv(df_cv_knn_formatted, paste0(output_path, "cross_validation_knn_results.csv"), row.names = FALSE)
  }
  return(df_cv_knn_formatted)
}

# plot_cross_validation_knn function: generate a plot of cross-validation results for KNN
# Input:
#       - df_cv_knn: data frame with cross-validation results for KNN (obtained from perform_cross_validation_knn function or equivalent)
#       - output_path: path to save the plot
# Output:
#       - ggplot object of the cross-validation plot
plot_cross_validation_knn <- function(df_cv_knn, output_path = NULL) {
  # Create a df with the cv_knn results in the format mean±sd (using function from cross_validation.r)
  df_annotation <- format_cross_validation_results(df_cv_knn)

  p <- df_cv_knn %>%
    ggplot(aes(x = factor(k))) +  # Convert x-axis to categorical
    geom_line(aes(y = mean_training, color = "Training", group = 1)) +
    geom_point(aes(y = mean_training, color = "Training")) +
    geom_errorbar(aes(ymin = mean_training - sd_training, ymax = mean_training + sd_training, color = "Training"), width = 0.1) +
    geom_line(aes(y = mean_testing, color = "Testing", group = 1)) +
    geom_point(aes(y = mean_testing, color = "Testing")) +
    geom_text(aes(y = mean_testing, label = df_annotation$Testing, color = "Testing"),
              hjust = 1, vjust = -1.3, size = 5, angle = 90) +
    geom_errorbar(aes(ymin = mean_testing - sd_testing, ymax = mean_testing + sd_testing, color = "Testing"), width = 0.1) +
    labs(title = "Cross-Validation for KNN", x = "Number of KNN Neighbors", y = "Accuracy") +
    scale_color_manual(name = "Metric", values = c("Training" = "blue", "Testing" = "red")) +
    theme_minimal() +
    theme(text = element_text(size = 18))
  
  # Save plot if output path is provided
  if (!is.null(output_path)) {
    ggsave(paste0(output_path, "cross_validation_knn_plot.pdf"), plot = p, width = 20, height = 10, dpi = 300)
  }
  
  return(p)
}