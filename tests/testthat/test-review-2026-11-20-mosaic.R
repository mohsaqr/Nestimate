# Review 2026-11-20, mosaic area: mosaic_test(), mosaic_plot(),
# mosaic_analysis() and the nested-association engine. Expected values come
# from hand computation or from an independent route through the package.
testthat::skip_on_cran()

# ---- H2: metadata named like an internal column ----

test_that("H2: an actor column called state/group/key/count gives the same clusters", {
  d <- data.frame(g = rep(c("x", "y"), each = 4),
                  team = rep(letters[1:4], each = 2), id = rep(1:2, 4),
                  T1 = rep(c("A", "B"), each = 2, times = 2),
                  T2 = c("A", "A", "A", "B", "A", "B", "B", "B"))
  run <- function(d, outer) {
    x <- build_network(d, method = "relative", group = "g",
                       state_cols = c("T1", "T2"))
    mosaic_test(x, method = "rao_scott", actor = c(outer, "id"))
  }
  # Independent route: the group x state x team table by hand.
  long <- data.frame(g = rep(d$g, 2), state = c(d$T1, d$T2),
                     team = rep(d$team, 2))
  ref <- mosaic_test(table(long$g, long$state, long$team),
                     method = "rao_scott")
  team <- run(d, "team")
  expect_identical(team$clusters, 4L)
  expect_equal(team$statistic, ref$statistic, tolerance = 1e-12)
  invisible(lapply(c("state", "group", "key", "count"), function(nm) {
    renamed <- d
    names(renamed)[names(renamed) == "team"] <- nm
    r <- run(renamed, nm)
    expect_identical(r$clusters, 4L, info = nm)
    expect_equal(r$statistic, ref$statistic, tolerance = 1e-12, info = nm)
    expect_equal(r$p_value, ref$p_value, tolerance = 1e-12, info = nm)
  }))
})

# ---- H3 / P3: undirected co-occurrence networks ----

test_that("H3/P3: the tests refuse a co-occurrence network", {
  d <- data.frame(T1 = c("A", "A", "A", "B"), T2 = c("A", "B", "B", "B"))
  net <- build_network(d, method = "co_occurrence")
  expect_error(mosaic_test(net), class = "nestimate_mosaic_undirected")
  expect_error(mosaic_test(net, method = "chisq"),
               class = "nestimate_mosaic_undirected")
  expect_error(mosaic_plot(net, report = TRUE),
               class = "nestimate_mosaic_undirected")
  expect_error(mosaic_plot(net, nested = TRUE),
               class = "nestimate_mosaic_undirected")
})

test_that("H3/P3: the co-occurrence mosaic counts each unordered pair once", {
  # Pairs AA, AB, AB, BB: unordered counts AA 1, AB 2, BB 1, four pairs.
  d <- data.frame(T1 = c("A", "A", "A", "B"), T2 = c("A", "B", "B", "B"))
  net <- build_network(d, method = "co_occurrence")
  p <- mosaic_plot(net, residuals = "asymptotic")
  freq <- stats::setNames(p$data$freq, paste0(p$data$row, p$data$col))
  expect_equal(freq[c("AA", "AB", "BA", "BB")],
               c(AA = 1, AB = 1, BA = 1, BB = 1))
  expect_equal(sum(p$data$freq), nrow(d))
  # Ends: A 4, B 4, so p = 1/2; expected unordered pairs 4 x (1/4, 1/2,
  # 1/4) = (1, 2, 1) = observed, so every residual is 0.
  expect_equal(p$data$stdres, rep(0, 4))
})

test_that("H3/P3: co-occurrence residuals follow random pairing by hand", {
  # Pairs AA, AA, BB, BB: U = (AA 2, AB 0, BB 2), M = 4, p = (1/2, 1/2).
  # AA: (2 - 1) / sqrt(4 * 1/4 * 3/4); AB: (0 - 2) / sqrt(4 * 1/2 * 1/2).
  d <- data.frame(T1 = c("A", "A", "B", "B"), T2 = c("A", "A", "B", "B"))
  net <- build_network(d, method = "co_occurrence")
  p <- mosaic_plot(net, residuals = "asymptotic")
  res <- stats::setNames(p$data$stdres, paste0(p$data$row, p$data$col))
  expect_equal(unname(res[c("AA", "BB")]), rep(1 / sqrt(0.75), 2),
               tolerance = 1e-12)
  expect_equal(unname(res[c("AB", "BA")]), c(-2, -2), tolerance = 1e-12)
  # Permutation residuals are symmetric across the mirrored cells.
  q <- mosaic_plot(net, seed = 1, n_perm = 200)
  rq <- stats::setNames(q$data$stdres, paste0(q$data$row, q$data$col))
  expect_identical(rq[["AB"]], rq[["BA"]])
})

test_that("H3: the pair table of 40 sequences of 4 holds 40 x choose(4, 2) pairs", {
  set.seed(3)
  m <- matrix(sample(c("A", "B", "C"), 160, TRUE), 40, 4)
  d <- stats::setNames(as.data.frame(m), paste0("T", 1:4))
  net <- build_network(d, method = "co_occurrence")
  p <- mosaic_plot(net, residuals = "asymptotic")
  expect_equal(sum(p$data$freq), 40 * choose(4, 2))
})

test_that("H3: undirected networks whose diagonal counts occurrences have no mosaic", {
  co <- cooccurrence(list(c("A", "B"), c("A", "C"), c("A", "B", "C")))
  expect_error(mosaic_plot(co), class = "nestimate_mosaic_undirected")
  d <- data.frame(T1 = c("A", "A"), T2 = c("A", "B"), T3 = c("B", "B"),
                  T4 = c("A", "B"))
  windowed <- build_network(d, method = "co_occurrence",
                            params = list(window_size = 2))
  expect_error(mosaic_plot(windowed), class = "nestimate_mosaic_undirected")
})

# ---- P2 / B4: probabilities are never read as counts ----

test_that("P2/B4: a relative network with 0/1 probabilities is tested on its counts", {
  d <- data.frame(T1 = rep(c("A", "B"), each = 4),
                  T2 = rep(c("B", "A"), each = 4))
  net <- build_network(d, method = "relative")
  # Counts A->B 4, B->A 4: E = 2 per cell, X2 = 4 x 2^2 / 2 = 8, n = 8.
  r <- mosaic_test(net, method = "chisq")
  expect_equal(r$statistic, 8)
  expect_identical(r$n, 8)
  # Adjusted residual (4 - 2) / sqrt(2 x 1/2 x 1/2) = 2 sqrt(2).
  p <- mosaic_plot(net, residuals = "asymptotic")
  expect_equal(abs(p$data$stdres), rep(2 * sqrt(2), 4), tolerance = 1e-12)
  expect_s3_class(mosaic_plot(net, nested = TRUE), "ggplot")
})

test_that("B4: mosaic_test() on 0/1 transitions with chisq alone runs on counts", {
  d <- data.frame(T1 = c("A", "A", "B", "B"), T2 = c("B", "B", "A", "A"))
  net <- build_network(d, method = "relative")
  r <- mosaic_test(net, method = "chisq")
  # Counts 0/2/2/0: E = 1 per cell, X2 = 4 x 1 / 1 = 4.
  expect_equal(r$statistic, 4)
  expect_identical(r$n, 4)
})

test_that("P2: a thresholded frequency network is drawn from its full counts", {
  d <- data.frame(T1 = c("A", "A", "A", "B", "B"),
                  T2 = c("A", "B", "B", "A", "B"))
  full <- build_network(d, method = "frequency")
  cut <- build_network(d, method = "frequency", threshold = 1.5)
  expect_identical(mosaic_test(cut, method = "chisq")$n, 5)
  expect_equal(mosaic_test(cut, method = "chisq")$statistic,
               mosaic_test(full, method = "chisq")$statistic)
})

test_that("P2: a probability network without counts or sequences is refused", {
  d <- data.frame(T1 = c("A", "A", "B", "B"), T2 = c("B", "B", "A", "A"))
  net <- build_network(d, method = "relative")
  net$frequency_matrix <- NULL
  net$data <- NULL
  expect_error(mosaic_plot(net), class = "nestimate_bad_input")
})

# ---- H9: Cramer's V ----

test_that("H9: Cramer's V is 1 for a perfect 2 x 2 association", {
  d <- data.frame(x = rep(c("A", "B"), each = 10),
                  y = rep(c("C", "D"), each = 10))
  expect_equal(mosaic_analysis(d, "x", "y")$cramers_v, 1)
})

test_that("H9: Cramer's V uses the uncorrected Pearson statistic", {
  # Table (8, 2; 3, 7), n = 20, all margins 10/10 except columns 11/9.
  d <- data.frame(x = rep(c("A", "B"), each = 10),
                  y = c(rep("C", 8), rep("D", 2), rep("C", 3), rep("D", 7)))
  tab <- matrix(c(8, 3, 2, 7), 2)
  e <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  v <- sqrt(sum((tab - e)^2 / e) / sum(tab))
  expect_equal(mosaic_analysis(d, "x", "y", min_count = 1)$cramers_v, v,
               tolerance = 1e-12)
})

# ---- B5: one outer cluster ----

test_that("B5: one outer cluster is named as the reason", {
  d <- data.frame(team = "T", student = rep(1:2, each = 4),
                  a = rep(c("A", "A", "B", "B"), 2),
                  b = c("X", "X", "Y", "Y", "X", "Y", "X", "Y"))
  expect_error(mosaic_test(d, var1 = "a", var2 = "b",
                           actor = c("team", "student"), method = "rao_scott"),
               class = "nestimate_bad_input")
  r <- mosaic_test(d, var1 = "a", var2 = "b", actor = c("team", "student"),
                   method = c("chisq", "rao_scott"))
  expect_match(r$note[r$method == "rao_scott"], "Fewer than two clusters")
})

# ---- P11: empty rows keep the other residuals ----

test_that("P11: an empty row is drawn at 0 and the other cells keep theirs", {
  # Rows (8, 2), (2, 8): E = 5, (8 - 5) / sqrt(5 x 1/2 x 1/2) = 2.683.
  p <- mosaic_plot(rbind(c(8, 2), c(2, 8), c(0, 0)),
                   residuals = "asymptotic")
  expect_equal(abs(p$data$stdres[p$data$freq > 0]), rep(3 / sqrt(1.25), 4),
               tolerance = 1e-12)
  expect_equal(p$data$stdres[p$data$row == "C"], c(0, 0))
})
