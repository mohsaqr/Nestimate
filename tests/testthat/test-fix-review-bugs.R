# Regression tests for the bugs left open after the 2026-11-08 review; each
# expectation failed on the code before the fix.
testthat::skip_on_cran()

grl <- group_regulation_long
roles <- function(...) build_network(grl, actor = "Actor", action = "Action",
                                     time = "Time", ...)

test_that("entropy_network() refuses scalings with a zero denominator", {
  det <- matrix(c(0, 1, 1, 0), 2, dimnames = list(c("a", "b"), c("a", "b")))
  expect_error(entropy_network(det, scaling = "share"),
               "positive entropy rate", class = "nestimate_bad_input")
  one <- matrix(1, 1, 1, dimnames = list("A", "A"))
  expect_error(entropy_network(one, weight = "surprisal", scaling = "chance"),
               "at least 2 states", class = "nestimate_bad_input")
})

test_that("entropy_network() keeps a bare cograph_network a netobject", {
  cg <- build_network(trajectories, method = "relative")
  class(cg) <- "cograph_network"
  expect_true(inherits(entropy_network(cg), "netobject"))
})

test_that("entropy_network() warns on a non-ergodic chain, as transition_entropy()", {
  P3 <- rbind(c(1, 0, 0), c(.5, .5, 0), c(0, .5, .5))
  dimnames(P3) <- list(letters[1:3], letters[1:3])
  expect_warning(entropy_network(P3), class = "nestimate_not_ergodic")
  expect_warning(transition_entropy(P3), class = "nestimate_not_ergodic")
})

test_that("entropy_bayes() networks declare the TNA geometry", {
  eb <- entropy_bayes(roles(method = "frequency"), draws = 200, seed = 1)
  expect_true(isTRUE(eb$network$meta$splot$defaults$tna_styling))
  expect_false(eb$network$meta$splot$defaults$psych_styling)
  expect_true(isTRUE(eb$model$meta$splot$defaults$tna_styling))
})

test_that("entropy_bayes() reads the output of wtna()", {
  oh <- data.frame(actor = rep(c("u1", "u2"), each = 6),
                   A = rep(c(1, 0, 1), 4), B = rep(c(0, 1, 1), 4),
                   C = rep(c(1, 1, 0), 4))
  w <- wtna(oh, actor = "actor")
  expect_true(inherits(entropy_bayes(w, draws = 200, seed = 1),
                       "net_entropy_bayes"))
})

test_that("grouped entropy results have summary() and plot()", {
  te <- transition_entropy(roles(method = "relative", group = "Achiever"))
  s <- summary(te)
  expect_true(is.data.frame(s))
  expect_identical(s$group, c("High", "Low"))
  expect_equal(s$entropy_rate, c(te$High$entropy_rate, te$Low$entropy_rate))
  expect_true(inherits(plot(te), "ggplot"))

  eb <- entropy_bayes(roles(method = "frequency", group = "Achiever"),
                      draws = 200, seed = 1)
  sb <- summary(eb)
  expect_true(is.data.frame(sb))
  expect_equal(nrow(sb), 6L)
  expect_true(inherits(plot(eb), "ggplot"))
})

test_that("entropy prints align their lines and never show -0.000", {
  eb <- entropy_bayes(roles(method = "frequency"), draws = 200, seed = 1)
  out <- utils::capture.output(print(eb))
  lines <- grep(":", out[3:5], value = TRUE)
  expect_true(all(startsWith(lines, "  ") & !startsWith(lines, "   ")))
  one <- utils::capture.output(print(
    transition_entropy(matrix(1, 1, 1, dimnames = list("A", "A")))))
  expect_false(any(grepl("-0.000", one, fixed = TRUE)))
})

test_that("entropy_trajectory() checks the time column and labels its axes", {
  d <- grl
  d$Tc <- format(d$Time)
  expect_error(entropy_trajectory(d, action = "Action", actor = "Actor",
                                  time = "Tc", window = 500),
               class = "nestimate_bad_input")
  tr <- entropy_trajectory(grl, action = "Action", actor = "Actor",
                           window = 500, base = exp(1))
  labs <- plot(tr)$labels
  expect_identical(labs$x, "Row of the event ending the transition")
  expect_identical(labs$y, "Entropy (nats per transition)")
  d$tn <- as.numeric(d$Time)
  trn <- entropy_trajectory(d, action = "Action", actor = "Actor",
                            time = "tn", window = 500)
  expect_identical(plot(trn)$labels$x, "tn")
})

test_that("asymptotic mosaic residuals survive an empty column", {
  tab <- as.table(matrix(c(10, 2, 0, 3, 12, 0, 5, 4, 0), 3,
                         dimnames = list(c("a", "b", "c"), c("x", "y", "z"))))
  res <- mosaic_plot(tab, residuals = "asymptotic")$data$stdres
  expect_gt(max(abs(res)), 1)
  # the six defined residuals are kept; the empty column's three cells are 0
  full <- as.vector(suppressWarnings(stats::chisq.test(tab))$stdres)
  expect_equal(sort(res), sort(c(full[is.finite(full)], 0, 0, 0)))
})

test_that("state counts and mosaics refuse association and one-hot networks", {
  dd <- as.data.frame(matrix(stats::rnorm(300), 60, 5))
  dd$g <- rep(c("a", "b"), 30)
  cg <- build_network(dd, method = "cor", group = "g")
  expect_error(state_distribution(cg), class = "nestimate_bad_input")
  expect_error(mosaic_plot(cg), class = "nestimate_bad_input")
  expect_error(mosaic_plot(build_network(dd[1:5], method = "cor")),
               "holds variables", class = "nestimate_bad_input")
  oh <- data.frame(A = rep(c(1, 0, 1), 4), B = rep(c(0, 1, 1), 4),
                   C = rep(c(1, 1, 0), 4))
  expect_error(state_distribution(build_network(oh, method = "wtna")),
               class = "nestimate_bad_input")
})

test_that("compare_model() reports undefined edge metrics as NA, silently", {
  r <- compare_model(roles(method = "relative", group = "Achiever"))
  num <- vapply(r$edge_metrics, is.numeric, logical(1))
  expect_false(any(vapply(r$edge_metrics[num],
                          function(v) any(is.nan(v) | is.infinite(v)),
                          logical(1))))
  expect_true(anyNA(r$edge_metrics$similarity_strength_index))
  set.seed(3)
  a <- matrix(stats::rexp(25), 5, 5, dimnames = list(letters[1:5], letters[1:5]))
  b <- matrix(stats::rexp(25), 5, 5, dimnames = dimnames(a))
  expect_no_warning(z <- compare_model(a, b, scaling = "zscore"))
  expect_true(anyNA(z$edge_metrics$logarithmic_ratio))
})

test_that("build_mcml(method = 'density') agrees between the two levels", {
  s <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"), T3 = c("C", "A"))
  m <- build_mcml(s, list(G1 = c("A", "B"), G2 = "C"), method = "density",
                  type = "raw", cols = c("T1", "T2", "T3"))
  within <- m$clusters$G1$weights
  # macro density = within-cluster sum / (n_i * n_j)
  expect_equal(m$macro$weights["G1", "G1"], sum(within) / 4)
})

test_that("build_mcml() refuses a weight matrix with exclude", {
  sq <- matrix(c(0, 1, 2, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(build_mcml(sq, clusters = list(G1 = "A", G2 = "B"),
                          exclude = "A"),
               "weight matrix", class = "nestimate_bad_input")
})

test_that("persistent homology prints VR essentials and ascending range", {
  dd <- as.matrix(stats::dist(matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2,
                                     byrow = TRUE)))
  out <- utils::capture.output(print(persistent_homology(dd, type = "vr")))
  expect_true(any(grepl("b0: 4 (1 persistent)", out, fixed = TRUE)))
  expect_true(any(grepl("[0.0000 → 1.4142]", out, fixed = TRUE)))
  d1 <- matrix(1, 3, 3)
  diag(d1) <- 0
  expect_error(persistent_homology(d1, type = "vr", max_scale = 0),
               "No edge enters", class = "nestimate_bad_input")
})

test_that("a (0, 0) pair is finite in bottleneck_distance(), as in landscapes", {
  bd <- bottleneck_distance(data.frame(dimension = 0, birth = 0, death = 0),
                            data.frame(dimension = 0, birth = 0.3,
                                       death = 0.3))
  expect_true(all(is.finite(unlist(bd[vapply(bd, is.numeric, TRUE)]))))
})

test_that("evaluate_links() validates k and returns a table, never NULL", {
  rn <- build_network(trajectories, method = "relative")
  pr <- predict_links(rn)
  expect_error(evaluate_links(pr, extract_edges(rn), k = 0),
               class = "nestimate_bad_input")
  ev <- evaluate_links(pr, extract_edges(rn))
  expect_true(is.data.frame(ev))
  expect_named(ev, c("method", "auc", "average_precision", "precision_at_5",
                     "precision_at_10", "precision_at_20"))
})

test_that("a constant normalised centrality warns instead of a silent NaN", {
  two <- build_network(matrix(c(0, .5, .5, 0), 2,
                              dimnames = list(c("A", "B"), c("A", "B"))),
                       method = "relative")
  expect_warning(net_centrality(two), class = "nestimate_constant_measure")
})

test_that("landscape_trend() and network_reliability() validate their inputs", {
  expect_error(landscape_trend(trajectories, iter = Inf),
               class = "nestimate_bad_input")
  tiny <- build_network(data.frame(V1 = c("A", "B", "C", "A"),
                                   V2 = c("B", "C", "A", "B")),
                        method = "relative")
  expect_error(network_reliability(tiny, iter = 10, split = 0.9, seed = 1),
               "second half is empty", class = "nestimate_bad_input")
})

test_that("boot_glasso() and path_dependence() refuse inputs they misread", {
  expect_error(boot_glasso(stats::cor(srl_strategies), iter = 5),
               "raw data", class = "nestimate_bad_input")
  expect_error(path_dependence(matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1), 3),
                               min_count = 1),
               class = "nestimate_bad_input")
})

test_that("chain_structure() explains a state with no outgoing transitions", {
  net <- build_network(trajectories, method = "relative", end = TRUE)
  expect_error(chain_structure(net), "mark_terminal_state",
               class = "nestimate_bad_input")
})

test_that("mogen_transitions(min_count = 0) lists only continuing contexts", {
  mg <- mogen_transitions(build_mogen(trajectories, max_order = 2),
                          order = 2, min_count = 0)
  expect_false(anyDuplicated(mg$path) > 0)
  # 9 two-state contexts x 3 next states
  expect_equal(nrow(mg), 27L)
})

test_that("logical flags refuse NA; small message fixes", {
  expect_error(sequence_plot(trajectories, combined = NA),
               class = "nestimate_bad_input")
  err <- tryCatch(build_hypa(trajectories, p_adjust = "zzz"),
                  error = conditionMessage)
  expect_identical(lengths(regmatches(err, gregexpr("none", err))), 1L)
  rn <- build_network(trajectories, method = "relative")
  expect_error(compare_networks(list(a = rn)), class = "nestimate_bad_input")
})
