#' Network Comparison Functions
#'
#' @description
#' Functions for comparing TNA networks and evaluating recovery performance.
#'
#' @name network_comparison
#' @keywords internal
NULL

#' Compare Two TNA Networks
#'
#' @description
#' Compare two TNA network models using multiple metrics including correlation,
#' RMSE, and edge-level differences.
#'
#' @param model1 A TNA model object (the reference/original model).
#' @param model2 A TNA model object (the comparison/simulated model).
#' @param metrics Character vector. Metrics to compute. Options:
#'   "correlation" (Pearson correlation), "rmse" (root mean square error),
#'   "mae" (mean absolute error), "edge_diff" (proportion with large differences),
#'   "cosine" (cosine similarity), or "all" to compute everything.
#'   Default: c("correlation", "rmse", "edge_diff").
#' @param scaling Character. How to scale edge weights before comparison.
#'   Options: "none" (no scaling), "minmax" (min-max normalization to 0-1),
#'   "zscore" (z-score standardization). Default: "none".
#' @param include_self Logical. Whether to include self-loops in comparison.
#'   Default: TRUE.
#' @param threshold Numeric. Threshold for considering an edge as "different"
#'   in edge_diff metric. Default: 0.05.
#'
#' @return A list with class "tna_comparison" containing:
#' \itemize{
#'   \item metrics: Named list of computed metrics.
#'   \item edge_comparison: Data frame comparing edges from both models.
#'   \item summary: Character summary of the comparison.
#' }
#'
#' @details
#' This function extracts edge weights from both models and computes various
#' comparison metrics. Edge weights are extracted from the transition matrices
#' stored in the TNA model objects.
#'
#' @examples
#' \dontrun{
#' # Generate original data and fit model
#' original_data <- simulate_sequences(trans_mat, init_probs, 20, 200)
#' model_original <- tna::tna(original_data)
#'
#' # Simulate from fitted model and fit new model
#' sim_data <- simulate_sequences(
#'   transition_matrix = model_original$weights,
#'   initial_probabilities = model_original$initial,
#'   max_seq_length = 20, num_rows = 200
#' )
#' model_sim <- tna::tna(sim_data)
#'
#' # Compare the models
#' comparison <- compare_networks(model_original, model_sim)
#' print(comparison$metrics)
#' }
#'
#' @seealso [compare_centralities()] for comparing centrality profiles,
#'   [compare_edge_recovery()] for edge recovery metrics.
#'
#' @importFrom stats cor sd
#' @export
compare_networks <- function(model1,
                             model2,
                             metrics = c("correlation", "rmse", "edge_diff"),
                             scaling = c("none", "minmax", "zscore"),
                             include_self = TRUE,
                             threshold = 0.05) {
  # Input validation
  if (!inherits(model1, "tna") && !is.list(model1)) {
    stop("model1 must be a TNA model object or a list with 'weights'.")
  }
  if (!inherits(model2, "tna") && !is.list(model2)) {
    stop("model2 must be a TNA model object or a list with 'weights'.")
  }

  scaling <- match.arg(scaling)
  if ("all" %in% metrics) {
    metrics <- c("correlation", "rmse", "mae", "edge_diff", "cosine")
  }

  # Extract weight matrices
  weights1 <- extract_transition_matrix(model1)
  weights2 <- extract_transition_matrix(model2)

  # Check dimensions match
  if (!identical(dim(weights1), dim(weights2))) {
    stop("Models have different numbers of states.")
  }

  # Check state names match
  states1 <- rownames(weights1)
  states2 <- rownames(weights2)
  if (!is.null(states1) && !is.null(states2) && !identical(sort(states1), sort(states2))) {
    warning("Models have different state names. Results may not be meaningful.")
    # Align by matching names
    common_states <- intersect(states1, states2)
    weights1 <- weights1[common_states, common_states]
    weights2 <- weights2[common_states, common_states]
  }

  # Remove self-loops if requested
  if (!include_self) {
    diag(weights1) <- NA
    diag(weights2) <- NA
  }

  # Apply scaling
  if (scaling == "minmax") {
    w1_range <- range(weights1, na.rm = TRUE)
    w2_range <- range(weights2, na.rm = TRUE)
    if (diff(w1_range) > 0) {
      weights1 <- (weights1 - w1_range[1]) / diff(w1_range)
    }
    if (diff(w2_range) > 0) {
      weights2 <- (weights2 - w2_range[1]) / diff(w2_range)
    }
  } else if (scaling == "zscore") {
    w1_mean <- mean(weights1, na.rm = TRUE)
    w1_sd <- sd(as.vector(weights1), na.rm = TRUE)
    w2_mean <- mean(weights2, na.rm = TRUE)
    w2_sd <- sd(as.vector(weights2), na.rm = TRUE)
    if (w1_sd > 0) weights1 <- (weights1 - w1_mean) / w1_sd
    if (w2_sd > 0) weights2 <- (weights2 - w2_mean) / w2_sd
  }

  # Flatten to vectors for comparison
  vec1 <- as.vector(weights1)
  vec2 <- as.vector(weights2)
  valid_idx <- !is.na(vec1) & !is.na(vec2)
  vec1 <- vec1[valid_idx]
  vec2 <- vec2[valid_idx]

  # Compute metrics
  results <- list()

  if ("correlation" %in% metrics) {
    results$correlation <- if (length(vec1) > 2 && sd(vec1) > 0 && sd(vec2) > 0) {
      cor(vec1, vec2)
    } else {
      NA_real_
    }
  }

  if ("rmse" %in% metrics) {
    results$rmse <- sqrt(mean((vec1 - vec2)^2))
  }

  if ("mae" %in% metrics) {
    results$mae <- mean(abs(vec1 - vec2))
  }

  if ("edge_diff" %in% metrics) {
    # Proportion of edges with absolute difference > threshold
    results$edge_diff <- mean(abs(vec1 - vec2) > threshold)
  }

  if ("cosine" %in% metrics) {
    norm1 <- sqrt(sum(vec1^2))
    norm2 <- sqrt(sum(vec2^2))
    results$cosine <- if (norm1 > 0 && norm2 > 0) {
      sum(vec1 * vec2) / (norm1 * norm2)
    } else {
      NA_real_
    }
  }

  # Create edge comparison data frame
  states <- rownames(weights1)
  if (is.null(states)) states <- paste0("S", seq_len(nrow(weights1)))

  edge_comparison <- expand.grid(from = states, to = states, stringsAsFactors = FALSE)
  edge_comparison$weight1 <- as.vector(weights1)
  edge_comparison$weight2 <- as.vector(weights2)
  edge_comparison$diff <- edge_comparison$weight2 - edge_comparison$weight1
  edge_comparison$abs_diff <- abs(edge_comparison$diff)

  if (!include_self) {
    edge_comparison <- edge_comparison[edge_comparison$from != edge_comparison$to, ]
  }

  # Create summary
  summary_text <- sprintf(
    "Network Comparison Summary:\n- States: %d\n- Edges compared: %d%s\n%s",
    nrow(weights1),
    sum(valid_idx),
    if (!include_self) " (excluding self-loops)" else "",
    paste(
      sprintf("- %s: %.4f", names(results), unlist(results)),
      collapse = "\n"
    )
  )

  # Create result object
  result <- list(
    metrics = results,
    edge_comparison = edge_comparison,
    summary = summary_text
  )
  class(result) <- c("tna_comparison", "list")

  return(result)
}


#' Compare Centrality Profiles
#'
#' @description
#' Compare centrality measures between two TNA networks using correlation
#' or rank-based methods.
#'
#' @param model1 A TNA model object (the reference/original model).
#' @param model2 A TNA model object (the comparison/simulated model).
#' @param measures Character vector. Centrality measures to compare. Options:
#'   "OutStrength", "InStrength", "Betweenness", "Closeness", or "all".
#'   Default: c("OutStrength", "InStrength", "Betweenness").
#' @param method Character. Comparison method: "correlation" (Pearson),
#'   "rank" (Spearman), or "both". Default: "both".
#'
#' @return A list containing:
#' \itemize{
#'   \item correlations: Named list of correlation values by measure.
#'   \item centrality_comparison: Data frame with centrality values for each state.
#'   \item summary: Character summary of the comparison.
#' }
#'
#' @details
#' Centrality measures are extracted from the TNA model objects using
#' `tna::centralities()`. The function compares how well the centrality
#' profiles are preserved between the original and simulated/comparison model.
#'
#' @examples
#' \dontrun{
#' # Compare centrality profiles between models
#' comparison <- compare_centralities(model_original, model_simulated)
#' print(comparison$correlations)
#' }
#'
#' @seealso [compare_networks()] for full network comparison.
#'
#' @importFrom stats cor
#' @export
compare_centralities <- function(model1,
                                 model2,
                                 measures = c("OutStrength", "InStrength", "Betweenness"),
                                 method = c("both", "correlation", "rank")) {
  # Input validation
  method <- match.arg(method)

  if ("all" %in% measures) {
    measures <- c("OutStrength", "InStrength", "Betweenness", "Closeness")
  }

  # Extract centralities
  cent1 <- tryCatch(
    tna::centralities(model1),
    error = function(e) {
      stop("Failed to extract centralities from model1: ", e$message)
    }
  )
  cent2 <- tryCatch(
    tna::centralities(model2),
    error = function(e) {
      stop("Failed to extract centralities from model2: ", e$message)
    }
  )

  # Align states
  states1 <- rownames(cent1)
  states2 <- rownames(cent2)
  common_states <- intersect(states1, states2)

  if (length(common_states) < 3) {
    warning("Fewer than 3 common states; correlation may not be meaningful.")
  }

  cent1 <- cent1[common_states, , drop = FALSE]
  cent2 <- cent2[common_states, , drop = FALSE]

  # Filter to requested measures
  available_measures <- intersect(measures, colnames(cent1))
  missing_measures <- setdiff(measures, colnames(cent1))
  if (length(missing_measures) > 0) {
    warning("Measures not available: ", paste(missing_measures, collapse = ", "))
  }

  # Compute correlations
  correlations <- list()
  for (measure in available_measures) {
    v1 <- cent1[, measure]
    v2 <- cent2[, measure]

    if (method %in% c("correlation", "both")) {
      pearson <- if (sd(v1) > 0 && sd(v2) > 0) cor(v1, v2) else NA_real_
      correlations[[paste0(measure, "_pearson")]] <- pearson
    }
    if (method %in% c("rank", "both")) {
      spearman <- if (sd(v1) > 0 && sd(v2) > 0) {
        cor(v1, v2, method = "spearman")
      } else {
        NA_real_
      }
      correlations[[paste0(measure, "_spearman")]] <- spearman
    }
  }

  # Create comparison data frame
  comparison_df <- data.frame(
    state = common_states,
    stringsAsFactors = FALSE
  )
  for (measure in available_measures) {
    comparison_df[[paste0(measure, "_1")]] <- cent1[, measure]
    comparison_df[[paste0(measure, "_2")]] <- cent2[, measure]
  }

  # Create summary
  summary_text <- sprintf(
    "Centrality Comparison Summary:\n- States: %d\n- Measures: %s\n%s",
    length(common_states),
    paste(available_measures, collapse = ", "),
    paste(
      sprintf("- %s: %.4f", names(correlations), unlist(correlations)),
      collapse = "\n"
    )
  )

  result <- list(
    correlations = correlations,
    centrality_comparison = comparison_df,
    summary = summary_text
  )
  class(result) <- c("centrality_comparison", "list")

  return(result)
}


#' Calculate Edge Recovery Metrics
#'
#' @description
#' Calculate metrics for how well edges from an original network are recovered
#' in a simulated/comparison network.
#'
#' @param original A TNA model object (the original/ground truth model).
#' @param simulated A TNA model object (the simulated/recovered model).
#' @param threshold Numeric. Minimum edge weight to consider an edge as
#'   "present". Default: 0.01.
#' @param return_edges Logical. Whether to return detailed edge-level results.
#'   Default: FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item true_positives: Number of edges correctly present in both.
#'   \item false_positives: Number of edges present in simulated but not original.
#'   \item false_negatives: Number of edges present in original but not simulated.
#'   \item true_negatives: Number of edges correctly absent in both.
#'   \item precision: TP / (TP + FP).
#'   \item recall: TP / (TP + FN), also known as sensitivity.
#'   \item f1_score: Harmonic mean of precision and recall.
#'   \item accuracy: (TP + TN) / total edges.
#'   \item jaccard: TP / (TP + FP + FN), Jaccard similarity.
#'   \item edges: (Optional) Data frame with edge-level results.
#' }
#'
#' @details
#' This function treats edge recovery as a binary classification problem:
#' \itemize{
#'   \item True Positive: Edge present in both original and simulated.
#'   \item False Positive: Edge present in simulated but not original.
#'   \item False Negative: Edge present in original but not simulated.
#'   \item True Negative: Edge absent in both.
#' }
#'
#' An edge is considered "present" if its weight exceeds the threshold.
#'
#' @examples
#' \dontrun{
#' # Calculate edge recovery
#' recovery <- compare_edge_recovery(model_original, model_simulated)
#' print(sprintf("Precision: %.2f, Recall: %.2f", recovery$precision, recovery$recall))
#' }
#'
#' @seealso [compare_networks()] for full network comparison,
#'   [run_bootstrap_iteration()] for bootstrap evaluation.
#'
#' @export
compare_edge_recovery <- function(original,
                                    simulated,
                                    threshold = 0.01,
                                    return_edges = FALSE) {
  # Extract weight matrices
  weights_orig <- extract_transition_matrix(original)
  weights_sim <- extract_transition_matrix(simulated)

  # Align states
  states_orig <- rownames(weights_orig)
  states_sim <- rownames(weights_sim)
  common_states <- intersect(states_orig, states_sim)

  if (length(common_states) == 0) {
    stop("No common states between original and simulated models.")
  }

  weights_orig <- weights_orig[common_states, common_states]
  weights_sim <- weights_sim[common_states, common_states]

  # Binarize based on threshold
  present_orig <- weights_orig >= threshold
  present_sim <- weights_sim >= threshold

  # Calculate confusion matrix
  tp <- sum(present_orig & present_sim)
  fp <- sum(!present_orig & present_sim)
  fn <- sum(present_orig & !present_sim)
  tn <- sum(!present_orig & !present_sim)

  # Calculate metrics
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }
  accuracy <- (tp + tn) / (tp + fp + fn + tn)
  jaccard <- if ((tp + fp + fn) > 0) tp / (tp + fp + fn) else NA_real_

  result <- list(
    true_positives = tp,
    false_positives = fp,
    false_negatives = fn,
    true_negatives = tn,
    precision = precision,
    recall = recall,
    f1_score = f1,
    accuracy = accuracy,
    jaccard = jaccard
  )

  # Add edge-level details if requested
  if (return_edges) {
    edges_df <- expand.grid(
      from = common_states,
      to = common_states,
      stringsAsFactors = FALSE
    )
    edges_df$weight_original <- as.vector(weights_orig)
    edges_df$weight_simulated <- as.vector(weights_sim)
    edges_df$present_original <- as.vector(present_orig)
    edges_df$present_simulated <- as.vector(present_sim)
    edges_df$status <- ifelse(
      edges_df$present_original & edges_df$present_simulated, "TP",
      ifelse(!edges_df$present_original & edges_df$present_simulated, "FP",
             ifelse(edges_df$present_original & !edges_df$present_simulated, "FN", "TN"))
    )
    result$edges <- edges_df
  }

  class(result) <- c("edge_recovery", "list")
  return(result)
}


# Backward compatibility alias (silent - no warning)

#' @rdname compare_edge_recovery
#' @param ... Arguments passed to \code{compare_edge_recovery}.
#' @export
calculate_edge_recovery <- function(...) compare_edge_recovery(...)
