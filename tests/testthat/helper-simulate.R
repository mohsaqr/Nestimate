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


# ---- Equivalence test infrastructure ----

#' Skip equivalence tests unless explicitly enabled
#' @noRd
skip_equiv_tests <- function() {
  run <- Sys.getenv("NESTIMATE_EQUIV_TESTS", unset = "false")
  if (!identical(run, "true")) {
    skip("Equivalence tests skipped (set NESTIMATE_EQUIV_TESTS=true)")
  }
}

#' Generate random sequence data for equivalence testing
#' @noRd
simulate_sequences <- function(n_actors = 10, n_states = 5,
                               seq_length = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- lapply(seq_len(n_actors), function(i) {
    sample(states, seq_length, replace = TRUE)
  })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

#' Generate random continuous data for association method equivalence testing
#' @noRd
simulate_continuous <- function(n = 100, p = 5, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- diag(p)
  for (i in seq_len(p - 1L)) {
    Sigma[i, i + 1L] <- rho
    Sigma[i + 1L, i] <- rho
  }
  L <- chol(Sigma)
  mat <- matrix(rnorm(n * p), n, p) %*% L
  df <- as.data.frame(mat)
  colnames(df) <- paste0("V", seq_len(p))
  df
}

#' Equivalence report logger
#' @noRd
equiv_report <- function() {
  env <- new.env(parent = emptyenv())
  env$rows <- list()

  env$log <- function(func, config, n_checked, n_failed,
                      max_abs_err, mean_abs_err, median_abs_err,
                      p95_abs_err, reference, notes = "") {
    env$rows[[length(env$rows) + 1L]] <- data.frame(
      func = func, config = config,
      n_checked = n_checked, n_failed = n_failed,
      max_abs_err = max_abs_err, mean_abs_err = mean_abs_err,
      median_abs_err = median_abs_err, p95_abs_err = p95_abs_err,
      reference = reference, notes = notes,
      stringsAsFactors = FALSE
    )
  }

  env$write_csv <- function(module) {
    if (length(env$rows) == 0L) return(invisible(NULL))
    df <- do.call(rbind, env$rows)
    dir.create("../../tmp", showWarnings = FALSE, recursive = TRUE)
    path <- sprintf("../../tmp/%s_equivalence_report.csv", module)
    write.csv(df, path, row.names = FALSE)
    message(sprintf("Equivalence report: %s (%d checks)", path, sum(df$n_checked)))
  }

  env$write_cvs <- function(module, test_file = NULL) {
    if (length(env$rows) == 0L) return(invisible(NULL))
    df <- do.call(rbind, env$rows)

    # Build vitest-compatible assertion results
    assertions <- lapply(seq_len(nrow(df)), function(i) {
      r <- df[i, ]
      passed <- r$n_failed == 0
      title <- sprintf("%s: %s delta=%.2e", r$func, r$config, r$max_abs_err)
      list(
        ancestorTitles = list(paste0(module, " equivalence")),
        title = title,
        fullName = sprintf("%s equivalence > %s", module, title),
        status = if (passed) "passed" else "failed",
        duration = 0,
        failureMessages = if (passed) list() else list(
          sprintf("max delta %.2e >= tolerance, %d/%d values failed",
                  r$max_abs_err, r$n_failed, r$n_checked)
        ),
        `_cvs` = list(
          delta = r$max_abs_err,
          tolerance = 1e-10,
          rFunction = r$func,
          rPackage = r$reference,
          module = module,
          target = "nestimate"
        )
      )
    })

    n_passed <- sum(df$n_failed == 0)
    n_failed <- sum(df$n_failed > 0)
    result <- list(
      numTotalTestSuites = 1L,
      numPassedTestSuites = if (n_failed == 0) 1L else 0L,
      numFailedTestSuites = if (n_failed > 0) 1L else 0L,
      numTotalTests = nrow(df),
      numPassedTests = n_passed,
      numFailedTests = n_failed,
      testResults = list(list(
        name = if (!is.null(test_file)) test_file
               else sprintf("tests/testthat/test-equiv-%s.R", module),
        assertionResults = assertions
      ))
    )

    inbox <- file.path("..", "..", "..", "validation", "data", "inbox")
    if (!dir.exists(inbox)) inbox <- "../../validation/data/inbox"
    if (!dir.exists(inbox)) {
      # Try absolute path
      inbox <- "/Users/mohammedsaqr/Documents/Github/validation/data/inbox"
    }
    if (dir.exists(inbox)) {
      ts <- format(Sys.time(), "%Y%m%dT%H%M%S")
      path <- file.path(inbox, sprintf("nestimate-%s-%s.json", module, ts))
      writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), path)
      message(sprintf("CVS report: %s", path))
    }
  }

  env$summary <- function() {
    if (length(env$rows) == 0L) return("No results logged.")
    df <- do.call(rbind, env$rows)
    sprintf("Total: %d values checked, %d failed, max delta %.2e",
            sum(df$n_checked), sum(df$n_failed), max(df$max_abs_err))
  }

  env
}
