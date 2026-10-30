# Regression tests for the defects found by the review of everything changed
# since CRAN 0.9.24 (2026-10-30). One block per defect.
testthat::skip_on_cran()

test_that("numeric-coded wide states with an id column build the real network", {
  d <- data.frame(id = c("a", "b", "c"), T1 = c(1, 2, 3), T2 = c(2, 3, 1),
                  T3 = c(3, 1, 2))
  net <- build_network(d, method = "relative")
  expect_equal(sort(net$nodes$label), c("1", "2", "3"))
  expect_equal(unname(net$weights["1", "2"]), 1)
  expect_equal(names(net$metadata), "id")
})

test_that("sequence verbs refuse networks whose data are not sequences", {
  set.seed(1)
  df <- data.frame(a = stats::rnorm(30), b = stats::rnorm(30),
                   c = stats::rnorm(30))
  expect_error(sequence_patterns(build_network(df, method = "cor"), len = 2),
               class = "nestimate_bad_input")
  expect_error(sequence_metrics(build_network(df, method = "glasso")),
               class = "nestimate_bad_input")
})

test_that("a network grouped by two columns works with the sequence verbs", {
  ng <- build_network(group_regulation_long, method = "relative",
                      actor = "Actor", action = "Action", time = "Time",
                      group = c("Achiever", "Course"))
  expect_equal(nrow(extract_sequences(ng)), 2000L)
  sp <- sequence_patterns(ng, len = 2)
  expect_true("support_High-A" %in% names(sp))
  expect_equal(nrow(sequence_metrics(ng)), 2000L)
})

test_that("sequence_metrics: a row with no observed state is NA throughout", {
  er <- data.frame(T1 = c("A", NA), T2 = c("B", NA), T3 = c("A", NA))
  m <- sequence_metrics(er)
  expect_true(is.na(m$valid_proportion[2]))
  expect_true(is.na(m$cyclic_feedback_strength[2]))
})

test_that("pattern templates accept the spaced HON notation", {
  expect_equal(
    nrow(sequence_patterns(trajectories, pattern = "Active -> * -> Disengaged")),
    nrow(sequence_patterns(trajectories, pattern = "Active->*->Disengaged")))
})

test_that("nct() matches an unnamed matrix to a named one by position", {
  set.seed(1)
  m1 <- matrix(stats::rnorm(300), 100, 3, dimnames = list(NULL, c("a", "b", "c")))
  m2 <- matrix(stats::rnorm(300), 100, 3)
  fit <- nct(m1, m2, iter = 5)
  expect_equal(unique(c(fit$E$edge_names$from, fit$E$edge_names$to)),
               c("a", "b", "c"))
  m3 <- m2
  colnames(m3) <- c("a", "c", "b")
  expect_error(nct(m1, m3, iter = 5), class = "nestimate_bad_input")
})

test_that("a missing value inside a sequence splits it (HON, MOGen, path_counts)", {
  d <- data.frame(t1 = c("A", "B", "A"), t2 = c(NA, "A", "B"),
                  t3 = c("B", "B", NA), t4 = c("A", "A", NA))
  expect_s3_class(build_hon(d, max_order = 2), "net_hon")
  expect_s3_class(build_hon(d, max_order = 2, method = "hon"), "net_hon")
  m <- build_mogen(d, max_order = 2)
  expect_false("NA" %in% colnames(m$transition_matrices[[2]]))
  runs <- build_mogen(list(c("B", "A"), c("B", "A", "B", "A"), c("A", "B")),
                      max_order = 2)
  expect_equal(m$log_likelihood, runs$log_likelihood, tolerance = 1e-12)
  pc <- path_counts(data.frame(t1 = c("A", "B"), t2 = c("B", NA),
                               t3 = c("C", "A"), t4 = c("A", "C")), k = 2)
  expect_false("B -> A" %in% pc$path)
  expect_equal(sum(pc$count), 4L)
})

test_that("association_rules() reads factor columns by their labels", {
  df <- data.frame(t1 = factor(c("x", "y", "x", "y")), t2 = c("y", "z", "z", "z"))
  expect_setequal(association_rules(df, min_support = 0.1,
                                    min_confidence = 0.1)$items,
                  c("x", "y", "z"))
})

test_that("HYPA pathway simplices leave out paths that never occur", {
  tr <- c(rep(list(c("A", "B", "C")), 20), rep(list(c("X", "B", "Y")), 20))
  sc <- build_simplicial(build_hypa(tr), type = "pathway")
  tri <- lapply(Filter(\(s) length(s) == 3, sc$simplices), \(s) sort(sc$nodes[s]))
  expect_setequal(vapply(tri, paste, character(1), collapse = ""),
                  c("ABC", "BXY"))
})

test_that("persistent_homology() row names have no gaps", {
  pts <- as.matrix(stats::dist(data.frame(
    x = cos(seq(0, 2 * pi, length.out = 9)[-9]),
    y = sin(seq(0, 2 * pi, length.out = 9)[-9]))))
  ph <- persistent_homology(pts, max_dim = 1)$persistence
  expect_equal(rownames(ph), as.character(seq_len(nrow(ph))))
})

test_that("labels = on numeric-coded states keeps the bootstrap correct", {
  set.seed(7)
  d <- data.frame(T1 = sample(c(1, 2, 10), 30, TRUE),
                  T2 = sample(c(2, 10), 30, TRUE),
                  T3 = sample(c(1, 2, 10), 30, TRUE))
  net <- build_network(d, method = "frequency", labels = c("1" = "one"))
  b <- summary(bootstrap_network(net, iter = 100, seed = 1))
  expect_lt(max(abs(b$weight - b$mean)), 1.5)
})

test_that("landscape surfaces survive flat heights and refuse degenerate axes", {
  set.seed(3)
  gen <- function() {
    s <- character(12)
    s[1] <- sample(c("A", "B", "C"), 1)
    for (i in 2:12) s[i] <- sample(setdiff(c("A", "B", "C"), s[i - 1]), 1)
    s
  }
  m <- t(replicate(60, gen()))       # no self-transitions: every depth is 0
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  surf <- suppressWarnings(landscape_surface(m))
  expect_identical(plot(surf), surf)
  tr <- landscape_trend(m, iter = 5, seed = 1)
  expect_identical(plot(tr, type = "surface"), tr)
  const <- matrix(rep(c("A", "B", "C"), 200), 60, byrow = TRUE)
  expect_error(landscape_surface(const, measure = "transition_rate"),
               class = "nestimate_bad_input")
  expect_error(stability_landscape(matrix(c("A", NA, "B", NA, NA, "A"),
                                          ncol = 2, byrow = TRUE)),
               class = "nestimate_landscape_unestimable")
})

test_that("landscape_surface warns when slices rest on very few sequences", {
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Course")
  expect_warning(landscape_surface(net), class = "nestimate_thin_slice")
})

test_that("two-state occupancy landscapes report their one ridge once", {
  set.seed(2)
  m2 <- matrix(sample(c("A", "B"), 200, TRUE), 20)
  s <- summary(landscape_trend(m2, depth = "occupancy", iter = 0))
  expect_equal(sum(s$feature == "ridge" & s$value == min(s$value)), 1L)
})
