# build_mogen(), summary.net_mogen(), path_counts(), mogen_transitions(): the multi-order generative model of
# paths (Scholtes 2017; Gote and Scholtes 2023). References: pathpy 2.2.0 (the authors' own package: higher-order
# count matrices and the multi-order model's likelihood), a Python reimplementation of the model selection written
# for these tests (data/pymogen.py: log-likelihoods, degrees of freedom, AIC, BIC per order), and k-gram counts
# from their definition. Re-expressed from local_testing_and_equivalence/test-equiv-mogen.R and
# test-mogen-pymogen-equiv.R (which hard-coded /opt/homebrew/bin/python3). Python runs through reticulate with the
# interpreter of RETICULATE_PYTHON. The count matrices and per-order statistics are read from the result's fields
# (API gap: summary() rounds them to two decimals and there is no accessor for the matrices).

ho_mogen_inputs <- function() {
  set.seed(202)
  second_order <- lapply(seq_len(80L), \(i) {
    start <- sample(c("A", "B", "C"), 2L, replace = TRUE)
    Reduce(\(path, t) {
      tail2 <- paste(utils::tail(path, 2L), collapse = "")
      c(path, if (tail2 == "AB") "C" else if (tail2 == "CB") "A" else sample(c("A", "B", "C"), 1L))
    }, seq_len(10L), init = start)
  })
  wide <- ho_sequences(30L, 20L, ho_stochastic(4L, 901L), 902L)
  list(list(name = "simple 3 states", max_order = 3L,
            trajectories = list(c("A", "B", "C", "A", "B", "C"), c("B", "C", "A", "B", "C", "A"),
                                c("C", "A", "B", "C", "A", "B"), c("A", "B", "C", "A", "B", "C"),
                                c("A", "C", "B", "A", "C", "B"))),
       list(name = "second-order dependency", max_order = 3L, trajectories = second_order),
       list(name = "first-order chain, 4 states", max_order = 2L,
            trajectories = lapply(seq_len(nrow(wide)), \(i) unlist(wide[i, ], use.names = FALSE))))
}

# A count matrix as one row per non-zero cell, node names joined by ",".
ho_cells <- function(counts) {
  cells <- which(counts != 0, arr.ind = TRUE)
  data.frame(from = gsub("\x01", ",", rownames(counts)[cells[, 1L]], fixed = TRUE),
             to = gsub("\x01", ",", colnames(counts)[cells[, 2L]], fixed = TRUE), count = counts[cells])
}

test_that("per-order count matrices and log-likelihoods equal pathpy 2.2.0", {
  skip_without_reference("reticulate")
  reticulate::py_run_file(reference_data("data/pathpy_mogen.py"))
  invisible(lapply(ho_mogen_inputs(), \(input) {
    found <- build_mogen(input$trajectories, max_order = input$max_order)
    paths <- reticulate::py$`_paths_from_sequences`(input$trajectories)
    invisible(lapply(seq_len(input$max_order), \(k) {
      reference <- reticulate::py$mogen_hon_transition_matrix(paths, as.integer(k))
      counts <- do.call(rbind, lapply(reference$counts, unlist))
      dimnames(counts) <- list(unlist(reference$nodes), unlist(reference$nodes))
      expect_reference(ho_cells(found$count_matrices[[k + 1L]]), ho_cells(counts), reference = "pathpy",
                       covers = "build_mogen()", quantity = sprintf("%s: order %d counts", input$name, k), tolerance = 0,
                       key = c("from", "to"))
      # pathpy's MultiOrderModel(max_order = k).likelihood(); 1e-6, the original test's tolerance.
      expect_reference(found$log_likelihood[[k + 1L]],
                       reticulate::py$mogen_cumulative_likelihood(paths, as.integer(k)), reference = "pathpy",
                       covers = "build_mogen()", quantity = sprintf("%s: log-likelihood up to order %d", input$name, k),
                       tolerance = 1e-6)
    }))
  }))
})

test_that("log-likelihoods, degrees of freedom, AIC and BIC per order equal the Python model selection", {
  skip_without_reference("reticulate")
  reticulate::py_run_file(reference_data("data/pymogen.py"))
  invisible(lapply(ho_mogen_inputs(), \(input) {
    found <- build_mogen(input$trajectories, max_order = input$max_order, criterion = "aic", dof = "observed")
    reference <- reticulate::py$build_mogen_ref(reticulate::r_to_py(input$trajectories),
                                                max_order = as.integer(input$max_order))
    orders <- paste0("order_", 0:input$max_order)
    expected <- data.frame(order = orders, loglik = unlist(reference$logliks), layer_dof = unlist(reference$layer_dofs),
                           cum_dof = unlist(reference$cum_dofs), aic = unlist(reference$aics), bic = unlist(reference$bics))
    actual <- data.frame(order = orders, loglik = unname(found$log_likelihood), layer_dof = unname(found$layer_dof),
                         cum_dof = unname(found$dof), aic = unname(found$aic), bic = unname(found$bic))
    expect_reference(actual, expected, reference = "independent: MOGen model selection reimplemented in Python (data/pymogen.py)",
                     covers = "build_mogen()", quantity = sprintf("%s: per-order statistics", input$name),
                     tolerance = 1e-10, measure = "relative", key = "order")
    printed <- utils::capture.output(table <- summary(found))
    # summary() rounds the log-likelihood, AIC and BIC to two decimals.
    expect_reference(data.frame(order = orders, subset(table, select = c(loglik, aic, bic))),
                     subset(expected, select = c(order, loglik, aic, bic)),
                     reference = "independent: MOGen model selection reimplemented in Python (data/pymogen.py)",
                     covers = c("build_mogen()", "summary.net_mogen()"),
                     quantity = sprintf("%s: summary, rounded", input$name), tolerance = 5e-3, key = "order")
  }))
})

test_that("path counts and first-order transitions equal counts of consecutive states", {
  invisible(lapply(ho_mogen_inputs(), \(input) {
    pairs <- unlist(lapply(input$trajectories, \(s) paste(utils::head(s, -1L), utils::tail(s, -1L), sep = " -> ")))
    counts <- table(pairs)
    expected <- data.frame(path = names(counts), count = as.vector(counts), proportion = as.vector(counts) / sum(counts))
    # path_counts() and mogen_transitions() round proportions and probabilities to four decimals.
    expect_reference(subset(path_counts(input$trajectories, k = 2L), select = c(path, count)),
                     subset(expected, select = c(path, count)), reference = "independent: counts of consecutive states, base R",
                     covers = "path_counts()", quantity = sprintf("%s: bigram counts", input$name), tolerance = 0, key = "path")
    expect_reference(subset(path_counts(input$trajectories, k = 2L), select = c(path, proportion)),
                     subset(expected, select = c(path, proportion)), reference = "independent: counts of consecutive states, base R",
                     covers = "path_counts()", quantity = sprintf("%s: bigram proportions", input$name), tolerance = 5e-5,
                     key = "path")
    from <- sub(" -> .*$", "", expected$path)
    transitions <- data.frame(path = expected$path, count = expected$count,
                              probability = expected$count / stats::ave(expected$count, from, FUN = sum))
    expect_reference(subset(mogen_transitions(build_mogen(input$trajectories, max_order = 1L), order = 1L),
                            select = c(path, count, probability)), transitions,
                     reference = "independent: counts of consecutive states, base R", covers = "mogen_transitions()",
                     quantity = sprintf("%s: first-order transitions", input$name), tolerance = 5e-5, key = "path")
  }))
})

# The trajectories as a wide table, NA after each one's end (the input tna takes).
ho_wide <- function(trajectories) {
  width <- max(lengths(trajectories))
  as.data.frame(do.call(rbind, lapply(trajectories, \(s) c(s, rep(NA_character_, width - length(s))))))
}

test_that("path counts and first-order transitions equal tna::ftna() and tna::tna()", {
  skip_without_reference("tna")
  invisible(lapply(ho_mogen_inputs(), \(input) {
    wide <- ho_wide(input$trajectories)
    counts <- out_edges(tna::ftna(wide)$weights, nonzero = TRUE)
    probabilities <- out_edges(tna::tna(wide)$weights, nonzero = TRUE)
    bigrams <- data.frame(path = paste(counts$from, counts$to, sep = " -> "), count = counts$weight,
                          proportion = counts$weight / sum(counts$weight))
    # path_counts() and mogen_transitions() round proportions and probabilities to four decimals.
    expect_reference(subset(path_counts(input$trajectories, k = 2L), select = c(path, count, proportion)), bigrams,
                     reference = "tna", covers = "path_counts()", quantity = sprintf("%s: bigram counts, tna", input$name),
                     tolerance = 5e-5, key = "path")
    transitions <- data.frame(path = paste(probabilities$from, probabilities$to, sep = " -> "),
                              count = counts$weight, probability = probabilities$weight)
    expect_reference(subset(mogen_transitions(build_mogen(input$trajectories, max_order = 1L), order = 1L),
                            select = c(path, count, probability)), transitions, reference = "tna",
                     covers = "mogen_transitions()", quantity = sprintf("%s: first-order transitions, tna", input$name),
                     tolerance = 5e-5, key = "path")
  }))
})

test_that("bigram counts equal the edge weights of pathpy's first-order network", {
  skip_without_reference("reticulate")
  reticulate::py_run_file(reference_data("data/pathpy_mogen.py"))
  invisible(lapply(ho_mogen_inputs(), \(input) {
    paths <- reticulate::py$`_paths_from_sequences`(input$trajectories)
    edges <- reticulate::py$hon_edges_df(paths, 1L)
    pathpy <- data.frame(path = vapply(edges, \(e) paste(e[[1L]], e[[2L]], sep = " -> "), character(1)),
                         count = vapply(edges, \(e) e[[3L]], numeric(1)))
    expect_reference(subset(path_counts(input$trajectories, k = 2L), select = c(path, count)), pathpy,
                     reference = "pathpy", covers = "path_counts()",
                     quantity = sprintf("%s: bigram counts, pathpy first-order network", input$name), tolerance = 0,
                     key = "path")
  }))
})

test_that("the summary's log-likelihoods equal pathpy's multi-order model of each maximal order", {
  skip_without_reference("reticulate")
  reticulate::py_run_file(reference_data("data/pathpy_mogen.py"))
  invisible(lapply(ho_mogen_inputs(), \(input) {
    found <- build_mogen(input$trajectories, max_order = input$max_order)
    paths <- reticulate::py$`_paths_from_sequences`(input$trajectories)
    orders <- seq_len(input$max_order)
    expected <- data.frame(order = orders,
                           loglik = vapply(orders, \(k) reticulate::py$mogen_cumulative_likelihood(paths, as.integer(k)),
                                           numeric(1)))
    printed <- utils::capture.output(table <- summary(found))
    # summary() rounds the log-likelihood to two decimals; its rows are the orders 0 to max_order.
    found_table <- subset(cbind(order = 0:input$max_order, subset(table, select = loglik)), order >= 1L)
    expect_reference(found_table, expected, reference = "pathpy",
                     covers = "summary.net_mogen()", quantity = sprintf("%s: summary log-likelihoods, pathpy", input$name),
                     tolerance = 5e-3, key = "order")
  }))
})
