# Minimal simulation helpers for Nestimate tests
# These replace Saqrlab simulation functions with standalone versions

simulate_mtna <- function(n_nodes = 5, n_types = 5,
                          type_names = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(type_names)) {
    all_types <- c("Metacognitive", "Cognitive", "Behavioral",
                   "Social", "Motivational", "Affective", "GroupRegulation")
    type_names <- all_types[seq_len(min(n_types, length(all_types)))]
    if (n_types > length(all_types)) {
      type_names <- c(type_names, paste0("Type", (length(all_types) + 1):n_types))
    }
  }
  total <- n_nodes * n_types
  node_names <- character(0)
  node_types_list <- list()
  for (i in seq_len(n_types)) {
    nms <- paste0(type_names[i], "_", seq_len(n_nodes))
    node_names <- c(node_names, nms)
    node_types_list[[type_names[i]]] <- nms
  }
  mat <- matrix(runif(total * total), nrow = total)
  diag(mat) <- 0
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  mat <- mat / row_sums
  rownames(mat) <- colnames(mat) <- node_names
  list(
    matrix = round(mat, 4),
    node_types = node_types_list,
    type_names = type_names,
    n_nodes_per_type = stats::setNames(rep(n_nodes, n_types), type_names)
  )
}

simulate_data <- function(type = "mlvar", seed = NULL, n_subjects = 20,
                          d = 6, n_obs = 40, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (type == "mlvar") {
    vars <- paste0("V", seq_len(d))
    rows <- n_subjects * n_obs

    # Generate known temporal and contemporaneous matrices
    true_B <- matrix(0, d, d)
    diag(true_B) <- runif(d, 0.2, 0.5)
    # Add a few off-diagonal temporal effects
    for (i in seq_len(min(d - 1, 3))) {
      true_B[i, i + 1] <- runif(1, 0.1, 0.3)
    }

    true_Theta <- diag(d)
    for (i in seq_len(min(d - 1, 2))) {
      true_Theta[i, i + 1] <- true_Theta[i + 1, i] <- runif(1, 0.1, 0.3)
    }

    # Generate VAR data
    between_means <- matrix(rnorm(n_subjects * d, 0, 1), n_subjects, d)
    df <- data.frame(
      id = rep(seq_len(n_subjects), each = n_obs),
      day = rep(rep(seq_len(n_obs %/% 10), each = 10), n_subjects),
      beep = rep(rep(seq_len(10), n_obs %/% 10), n_subjects)
    )
    y_mat <- matrix(0, rows, d)
    for (s in seq_len(n_subjects)) {
      idx <- ((s - 1) * n_obs + 1):(s * n_obs)
      y_mat[idx[1], ] <- between_means[s, ] + rnorm(d)
      for (t in 2:n_obs) {
        innovation <- rnorm(d)
        y_mat[idx[t], ] <- between_means[s, ] +
          true_B %*% (y_mat[idx[t - 1], ] - between_means[s, ]) + innovation
      }
    }
    for (j in seq_len(d)) df[[vars[j]]] <- y_mat[, j]

    rownames(true_B) <- colnames(true_B) <- vars
    rownames(true_Theta) <- colnames(true_Theta) <- vars

    attr(df, "vars") <- vars
    attr(df, "type") <- "mlvar"
    attr(df, "true_temporal") <- true_B
    attr(df, "true_contemporaneous") <- true_Theta
    attr(df, "info") <- list()
    return(df)
  }
  stop(sprintf("simulate_data type '%s' not supported in Nestimate test helper", type))
}
