# mosaic_test() and mosaic_analysis(): association tests of a two-way table, plain and with observations nested in
# actors. References: stats::chisq.test() (Pearson), survey::svychisq() (Rao-Scott F), stats::mantelhaen.test()
# (Cochran-Mantel-Haenszel), geepack::geeglm() (GEE Wald, identity link, independence) and aod::donner()
# (Donner-Klar). Re-expressed from local_testing_and_equivalence/test-equiv-nested-association.R (random designs
# compared by maxima) and tests/testthat/test-nested-association.R ("Rao-Scott residuals equal chisq.test()'s with
# one row per cluster", moved here). mosaic_analysis() keeps its unrounded test as the field `test` (API gap:
# summary() rounds).

# K clusters of 3-15 observations; rows x (R levels), columns y (C levels) with a cluster-specific shift.
sq_nested <- function(K, R, C, seed, cluster_rows = FALSE) {
  set.seed(seed)
  rows <- lapply(seq_len(K), \(k) {
    m <- sample(3:15, 1L)
    x <- if (cluster_rows) rep(sample(letters[seq_len(R)], 1L), m) else sample(letters[seq_len(R)], m, replace = TRUE)
    shift <- stats::rnorm(C, sd = 0.8)
    y <- vapply(x, \(r) sample(LETTERS[seq_len(C)], 1L, prob = exp(shift + if (r == "a") c(0.6, rep(0, C - 1L)) else 0)), "")
    data.frame(cl = sprintf("k%03d", k), x = x, y = y)
  })
  do.call(rbind, rows)
}

# One mosaic_test() row as the reference shapes it.
sq_test_row <- function(method, statistic, df, p_value, df2 = NULL) {
  out <- data.frame(method = method, statistic = unname(statistic), df = unname(df), p_value = unname(p_value))
  if (!is.null(df2)) out$df2 <- unname(df2)
  out
}

test_that("Pearson chi-square of a table equals stats::chisq.test()", {
  data <- sq_nested(30L, 3L, 3L, 11L)
  reference <- stats::chisq.test(table(data$x, data$y))
  expected <- sq_test_row("chisq", reference$statistic, reference$parameter, reference$p.value)
  expect_reference(mosaic_test(table(data$x, data$y, data$cl), method = "chisq"), expected, reference = "stats",
                   covers = c("mosaic_test()", "mosaic_test.table()"), quantity = "Pearson test, table", tolerance = 1e-10,
                   measure = "relative", key = "method")
  expect_reference(mosaic_test(data, var1 = "x", var2 = "y", actor = "cl", method = "chisq"), expected, reference = "stats",
                   covers = c("mosaic_test()", "mosaic_test.data.frame()"), quantity = "Pearson test, data frame",
                   tolerance = 1e-10, measure = "relative", key = "method")
  analysis <- mosaic_analysis(data, "x", "y")
  expect_reference(data.frame(test = "Pearson", statistic = unname(analysis$test$statistic),
                              df = unname(analysis$test$parameter), p = analysis$test$p.value),
                   data.frame(test = "Pearson", statistic = unname(reference$statistic), df = unname(reference$parameter),
                              p = reference$p.value),
                   reference = "stats", covers = "mosaic_analysis()", quantity = "Pearson test, mosaic_analysis()",
                   tolerance = 1e-10, measure = "relative", key = "test")
  observed <- reference$observed
  expect_reference(summary(analysis),
                   data.frame(x = rep(rownames(observed), ncol(observed)), y = rep(colnames(observed), each = nrow(observed)),
                              observed = as.vector(observed)),
                   reference = "stats", covers = c("mosaic_analysis()", "summary.mosaic_analysis()"),
                   quantity = "observed counts, summary()", tolerance = 0, key = c("x", "y"))
})

test_that("a network's Pearson test equals stats::chisq.test() of its transition counts", {
  network <- build_network(group_regulation_long, method = "frequency", actor = "Actor", action = "Action", time = "Time")
  reference <- suppressWarnings(stats::chisq.test(extract_transition_matrix(network)))
  expect_reference(mosaic_test(network, method = "chisq"),
                   sq_test_row("chisq", reference$statistic, reference$parameter, reference$p.value), reference = "stats",
                   covers = c("mosaic_test()", "mosaic_test.netobject()"), quantity = "Pearson test of transition counts",
                   tolerance = 1e-10, measure = "relative", key = "method")
  skip_without_reference("tna")
  model <- tna::ftna(tna::group_regulation)
  counts_test <- suppressWarnings(stats::chisq.test(model$weights))
  expect_reference(mosaic_test(model, method = "chisq"),
                   sq_test_row("chisq", counts_test$statistic, counts_test$parameter, counts_test$p.value), reference = "tna",
                   covers = c("mosaic_test()", "mosaic_test.default()"), quantity = "Pearson test of a tna frequency model",
                   tolerance = 1e-10, measure = "relative", key = "method",
                   note = "stats::chisq.test() on tna::ftna()'s counts")
})

test_that("the Rao-Scott F test equals survey::svychisq(statistic = 'F')", {
  skip_without_reference("survey")
  invisible(lapply(list(c(40L, 3L, 3L, 101L), c(60L, 2L, 4L, 102L), c(25L, 4L, 2L, 103L)), \(input) {
    data <- sq_nested(input[[1L]], input[[2L]], input[[3L]], input[[4L]])
    design <- survey::svydesign(ids = ~cl, weights = ~w, data = transform(data, w = 1))
    reference <- survey::svychisq(~x + y, design, statistic = "F")
    expect_reference(mosaic_test(table(data$x, data$y, data$cl), method = "rao_scott"),
                     sq_test_row("rao_scott", reference$statistic, reference$parameter[[1L]], reference$p.value,
                                 df2 = reference$parameter[[2L]]),
                     reference = "survey", covers = c("mosaic_test()", "mosaic_test.table()"),
                     quantity = sprintf("Rao-Scott F, %d clusters, %d x %d", input[[1L]], input[[2L]], input[[3L]]),
                     tolerance = 1e-10, measure = "relative", key = "method")
  }))
})

test_that("Rao-Scott residuals with one observation per cluster equal stats::chisq.test()'s standardised residuals", {
  table_2x3 <- matrix(c(12, 5, 7, 9, 4, 11), 2, 3, dimnames = list(c("a", "b"), c("x", "y", "z")))
  cells <- which(table_2x3 > 0, arr.ind = TRUE)
  replicate_rows <- rep(seq_len(nrow(cells)), table_2x3[cells])
  singletons <- array(0, c(2L, 3L, length(replicate_rows)), dimnames = c(dimnames(table_2x3), list(NULL)))
  singletons[cbind(cells[replicate_rows, , drop = FALSE], seq_along(replicate_rows))] <- 1
  residuals <- stats::chisq.test(table_2x3)$stdres
  expect_reference(mosaic_test(singletons, method = "rao_scott", level = "cells"),
                   data.frame(row = rep(rownames(residuals), 3L), column = rep(colnames(residuals), each = 2L),
                              residual = as.vector(residuals)),
                   reference = "stats", covers = c("mosaic_test()", "mosaic_test.array()"),
                   quantity = "Rao-Scott residual, singleton clusters", tolerance = 1e-12, key = c("row", "column"))
})

test_that("the Cochran-Mantel-Haenszel test equals stats::mantelhaen.test(correct = FALSE)", {
  invisible(lapply(list(c(6L, 3L, 3L, 301L), c(8L, 2L, 4L, 302L)), \(input) {
    data <- sq_nested(input[[1L]], input[[2L]], input[[3L]], input[[4L]])
    counts <- table(data$x, data$y, data$cl)
    reference <- stats::mantelhaen.test(counts, correct = FALSE)
    expect_reference(mosaic_test(counts, method = "cmh"),
                     sq_test_row("cmh", reference$statistic, reference$parameter, reference$p.value),
                     reference = "stats", covers = c("mosaic_test()", "mosaic_test.table()"),
                     quantity = sprintf("CMH, %d clusters", input[[1L]]), tolerance = 1e-10, measure = "relative",
                     key = "method")
  }))
})

test_that("the GEE Wald test equals geepack::geeglm() (identity link, independence)", {
  skip_without_reference("geepack")
  data <- sq_nested(40L, 3L, 2L, 202L)
  data <- data[order(data$cl), ]
  fit <- geepack::geeglm(yb ~ x, data = transform(data, yb = as.numeric(y == "A"), cl = factor(cl)), id = cl,
                         family = stats::gaussian, corstr = "independence")
  wald <- stats::anova(fit)
  expect_reference(subset(mosaic_test(table(data$x, data$y, data$cl), method = "gee"), select = c(method, statistic, df)),
                   data.frame(method = "gee", statistic = wald$X2[[1L]], df = wald$Df[[1L]]), reference = "geepack",
                   covers = c("mosaic_test()", "mosaic_test.table()"), quantity = "GEE Wald, 40 clusters",
                   tolerance = 1e-8, measure = "relative", key = "method")
})

test_that("the Donner-Klar test and intraclass correlation equal aod::donner()", {
  skip_without_reference("aod")
  data <- sq_nested(40L, 3L, 2L, 404L, cluster_rows = TRUE)
  totals <- stats::aggregate(cbind(n = 1, yA = as.numeric(data$y == "A")) ~ cl + x, data = data, FUN = sum)
  reference <- aod::donner(cbind(yA, n - yA) ~ x, data = totals)
  expect_reference(subset(mosaic_test(table(data$x, data$y, data$cl), method = "donner_klar"), select = c(method, statistic, icc)),
                   data.frame(method = "donner_klar", statistic = reference@X2, icc = reference@rho), reference = "aod",
                   covers = c("mosaic_test()", "mosaic_test.table()"), quantity = "Donner-Klar, 40 clusters",
                   tolerance = 1e-8, measure = "relative", key = "method")
})

test_that("an htna network's Pearson and CMH tests equal stats::chisq.test() and stats::mantelhaen.test() of its transitions", {
  wide <- tn_simulate(n_actors = 25L, n_states = 4L, seq_length = 10L, seed = 421L)
  sequences <- tn_sequences(wide)
  states <- sort(unique(unlist(sequences)))
  network <- as_htna(wide, clusters = list(X = states[1:2], Y = states[3:4]))
  # Mosaic rows are the targets and columns the sources: the transpose of the pair table, one stratum per sequence.
  per_sequence <- vapply(sequences, \(s) t(tn_pair_counts(list(s), states)), numeric(length(states)^2))
  strata <- array(per_sequence, c(length(states), length(states), length(sequences)))
  pearson <- stats::chisq.test(t(tn_pair_counts(sequences, states)))
  cmh <- stats::mantelhaen.test(strata, correct = FALSE)
  expect_reference(mosaic_test(network, method = c("chisq", "cmh")),
                   rbind(sq_test_row("chisq", pearson$statistic, pearson$parameter, pearson$p.value),
                         sq_test_row("cmh", cmh$statistic, cmh$parameter, cmh$p.value)),
                   reference = "stats", covers = c("mosaic_test()", "mosaic_test.htna()"),
                   quantity = "htna transitions, Pearson and CMH", tolerance = 1e-10, measure = "relative", key = "method")
})

test_that("a grouped network's tests of group by state equal stats::chisq.test() and survey::svychisq()", {
  skip_without_reference("survey")
  wide <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 422L)
  group <- rep(c("g1", "g2"), length.out = nrow(wide))
  network <- build_network(cbind(wide, grp = group), method = "relative", group = "grp")
  # One row per event: its sequence (the cluster), its group and its state.
  events <- do.call(rbind, lapply(seq_len(nrow(wide)), \(i) {
    s <- as.character(stats::na.omit(unlist(wide[i, ], use.names = FALSE)))
    data.frame(cl = i, grp = group[[i]], state = s)
  }))
  pearson <- stats::chisq.test(table(events$grp, events$state))
  design <- survey::svydesign(ids = ~cl, weights = ~w, data = transform(events, w = 1))
  rao_scott <- survey::svychisq(~grp + state, design, statistic = "F")
  expect_reference(mosaic_test(network, method = "chisq"),
                   sq_test_row("chisq", pearson$statistic, pearson$parameter, pearson$p.value), reference = "stats",
                   covers = c("mosaic_test()", "mosaic_test.netobject_group()"), quantity = "group by state, Pearson",
                   tolerance = 1e-10, measure = "relative", key = "method")
  expect_reference(mosaic_test(network, method = "rao_scott"),
                   sq_test_row("rao_scott", rao_scott$statistic, rao_scott$parameter[[1L]], rao_scott$p.value,
                               df2 = rao_scott$parameter[[2L]]),
                   reference = "survey", covers = c("mosaic_test()", "mosaic_test.netobject_group()"),
                   quantity = "group by state, Rao-Scott F", tolerance = 1e-10, measure = "relative", key = "method")
})
