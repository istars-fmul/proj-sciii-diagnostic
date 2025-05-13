library(clValid)  # compute connectivity
library(fpc)
library(factoextra) # eclust function
source("scripts/utils.R") # load utils functions

# compute_clustering_selection function: this function computes clustering metrics for a given dataset and clustering method.
# Input: 
# - df: the dataset to be clustered
# - fun_cluster: the clustering function to be used (e.g., kmeans, pam, hclust, agnes, diana)
# - hc_metric: the distance metric to be used (e.g., euclidean, manhattan)
# - hc_method: the method to be used for hierarchical clustering (e.g., single, complete, average, centroid, median, mcquitty, ward.D, ward.D2)
# - k: the number of clusters to be formed
# - balanced_membership_threshold: the threshold for balanced membership (default is 80 -> no cluster should have more than 80% of the total membership)
# Output: a data frame containing various clustering metrics, including connectivity, Dunn index, silhouette width, average between and within cluster distances, entropy, Pearson gamma, and membership percentages.
# df with metrics for the considering clustering method 
compute_clustering <- function(df, fun_cluster, hc_metric, hc_method = NULL, k, balanced_membership_threshold = 80) {
  # Compute distance matrix
   distance_res <- dist(df, method = 'euclidean')
  # Perform clustering
  clust_res <- eclust(df, FUNcluster = fun_cluster, hc_metric = hc_metric, hc_method = hc_method, k = k, graph = FALSE)
  
  # Compute membership percentages
  perc_membership <- table(clust_res$cluster) * 100 / length(clust_res$cluster)
  min_max_membership_ratio <- min(perc_membership) / max(perc_membership)
  
  # Check balanced membership
  balanced_membership <- max(perc_membership) <= balanced_membership_threshold
  
  # Compute clustering statistics
  cluster_stat <- cluster.stats(distance_res, clust_res$cluster)
  
  # Construct clustering method name
  fun_cluster_name <- ifelse(is.null(hc_method), fun_cluster, paste0(fun_cluster, "-", hc_method))
  
  # Return results as a named list
  return(list(metrics = data.frame(
    hc_metric = hc_metric, fun_cluster = fun_cluster_name, k = k,
    connectivity = connectivity(distance = distance_res, clusters = clust_res$cluster, neighbSize = 25),
    dunn = cluster_stat$dunn, silhouette = cluster_stat$avg.silwidth,
    average_between = cluster_stat$average.between, average_within = cluster_stat$average.within,
    entropy = cluster_stat$entropy, persongamma = cluster_stat$pearsongamma, wb_ratio = cluster_stat$wb.ratio, gini_index = gini(perc_membership),
    max_perc_membership = max(perc_membership), min_perc_membership = min(perc_membership), 
    median_perc_membership = median(perc_membership), check_balance_treshold = balanced_membership, min_max_membership_ratio = min_max_membership_ratio
  ), cluster = clust_res$cluster))
}

# evaluate_clustering_types function: this function evaluates different clustering methods and metrics for a given dataset.
# Input:
# - df: the dataset to be clustered
# - n_clusters: a vector of the number of clusters to be formed (default is c(6))
# - balanced_membership_threshold: the threshold for balanced membership (default is 80 -> no cluster should have more than 80% of the total membership)
# Output: a data frame containing the results of clustering evaluations for different methods and metrics.
# df with metrics for a wide range of clustering method (distance metric, clustering method, number of clusters)
evaluate_clustering_types <- function(df, n_clusters = c(6), balanced_membership_threshold = 80) {
  # Loop through clustering methods and metrics
  res_metrics <- data.frame()
  res_clust <- list()
  for (hc_metric in c("euclidean", "manhattan")) {
    for (k in n_clusters) {
      for (fun_cluster in c("kmeans", "pam", "hclust", "agnes", "diana")) {
        if (fun_cluster == "hclust" || fun_cluster == "agnes" || fun_cluster == "diana") {
          for (hc_method in c("single", "complete", "average", "centroid", "median", "mcquitty", "ward.D", "ward.D2")) {
            # Compute clusters and associated metrics
            if (!(fun_cluster == "agnes" && (hc_method == "median" || hc_method == "mcquitty" || hc_method == "centroid"))) { # skip invalid methods
                res_clustering_selection <- compute_clustering(df, fun_cluster, hc_metric, hc_method = hc_method, k = k, balanced_membership_threshold = balanced_membership_threshold)
              res_metrics <- rbind(res_metrics, res_clustering_selection$metrics)
              res_clust[[paste0(fun_cluster, "_", hc_method,"_", hc_metric)]] <- res_clustering_selection$cluster}
          }
        } else {
          # Compute clusters and associated metrics
          res_clustering_selection <- compute_clustering(df, fun_cluster, hc_metric, hc_method = NULL, k = k, balanced_membership_threshold = balanced_membership_threshold)
          res_metrics <- rbind(res_metrics, res_clustering_selection$metrics)
          res_clust[[paste0(fun_cluster,"_", hc_metric)]] <- res_clustering_selection$cluster
        }
      }
    }
  }
  #res_clust <- data.frame(res_clust)
  #rownames(res_clust) <- rownames(df)

  return(list(df_metrics = res_metrics, df_clusters = data.frame(do.call(cbind, res_clust), stringsAsFactors = FALSE)))
}