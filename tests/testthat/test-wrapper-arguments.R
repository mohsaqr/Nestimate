# Shortcut verbs name the arguments a user needs instead of hiding them in `...`.
testthat::skip_on_cran()

test_that("transition wrappers name the column roles and match build_network()", {
  for (fn in list(build_tna, build_ftna, build_atna, build_cna)) {
    expect_true(all(c("actor", "action", "time", "session", "start", "end") %in%
                      names(formals(fn))))
  }
  g <- group_regulation_long
  expect_identical(
    build_tna(g, actor = "Actor", action = "Action", time = "Time"),
    build_network(g, method = "relative", actor = "Actor", action = "Action",
                  time = "Time"))
  expect_identical(build_ftna(trajectories),
                   build_network(trajectories, method = "frequency"))
  expect_identical(build_cna(trajectories, start = TRUE),
                   build_network(trajectories, method = "co_occurrence",
                                 start = TRUE))
})

test_that("transition wrappers preserve positional start and end arguments", {
  seqs <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "B"))
  wrappers <- list(build_tna, build_ftna, build_atna, build_cna)
  methods <- c("relative", "frequency", "attention", "co_occurrence")
  invisible(Map(function(fn, method) {
    net <- fn(seqs, TRUE, TRUE)
    expect_identical(net, build_network(seqs, method = method,
                                       start = TRUE, end = TRUE))
    expect_true(all(c("Start", "End") %in% rownames(net$weights)))
  }, wrappers, methods))
})

test_that("association wrappers name actor (and level where it applies)", {
  for (fn in list(build_cor, build_pcor, build_glasso)) {
    expect_true(all(c("actor", "level") %in% names(formals(fn))))
  }
  expect_true("actor" %in% names(formals(build_ising)))
  expect_identical(build_glasso(srl_strategies),
                   build_network(srl_strategies, method = "glasso"))
  expect_identical(build_cor(srl_strategies, threshold = 0.1),
                   build_network(srl_strategies, method = "cor",
                                 threshold = 0.1))
})

test_that("a wrapper refuses `method` with a classed error", {
  expect_error(build_tna(trajectories, method = "frequency"),
               class = "nestimate_bad_input")
  expect_error(build_glasso(srl_strategies, method = "cor"),
               class = "nestimate_bad_input")
  expect_error(build_tna(trajectories, method = "frequency"),
               "always uses method = \"relative\"")
})

test_that("sequence verbs take actor/action/time by name", {
  for (fn in list(sequence_patterns, sequence_metrics, pattern_outcome)) {
    expect_true(all(c("actor", "action", "time") %in% names(formals(fn))))
  }
  g <- group_regulation_long
  seqs <- extract_sequences(g, actor = "Actor", action = "Action",
                            time = "Time")
  expect_identical(
    sequence_metrics(g, actor = "Actor", action = "Action", time = "Time"),
    sequence_metrics(seqs))
  expect_identical(
    sequence_patterns(g, actor = "Actor", action = "Action", time = "Time"),
    sequence_patterns(seqs))
})

test_that("compare_mmm() names build_mmm()'s settings and passes them on", {
  expect_true(all(c("n_starts", "seed", "covariates", "estimator") %in%
                    names(formals(compare_mmm))))
  cm <- compare_mmm(trajectories, k = 2:3, seed = 1, n_starts = 3)
  bic <- vapply(2:3, \(m) build_mmm(trajectories, k = m, seed = 1,
                                    n_starts = 3)$BIC, numeric(1))
  expect_identical(summary(cm)$BIC, bic)
})

test_that("macro_network() passes threshold on", {
  seqs <- data.frame(t1 = c("A", "C", "A", "B"), t2 = c("B", "D", "C", "A"),
                     t3 = c("C", "A", "D", "C"))
  mc <- build_mcml(seqs, clusters = list(G1 = c("A", "B"), G2 = c("C", "D")),
                   cols = c("t1", "t2", "t3"))
  w <- macro_network(mc, threshold = 0.5)$weights
  expect_equal(w[, "G1"], c(G1 = 0, G2 = 0))
  expect_equal(unname(w[, "G2"]), c(0.6, 2 / 3))
})
