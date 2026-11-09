# House-rule fixes in the resampling and centrality files (group F):
# bootstrap_network.R, boot_glasso.R, vertex_bootstrap.R,
# centrality_stability.R, centrality_measures.R, casedrop_reliability.R,
# reliability.R, certainty.R.

testthat::skip_on_cran()

.f_net <- function(...) {
  build_network(group_regulation_long, method = "relative",
                action = "Action", actor = "Actor", time = "Time", ...)
}

test_that("bootstrap and certainty prints label the interval with its level", {
  net <- .f_net()
  out_boot <- capture.output(
    print(bootstrap_network(net, iter = 20, seed = 1, ci_level = 0.1)))
  out_cert <- capture.output(print(certainty(net, ci_level = 0.1)))
  expect_true(any(grepl("90% CI", out_boot, fixed = TRUE)))
  expect_false(any(grepl("95% CI", out_boot, fixed = TRUE)))
  expect_true(any(grepl("90% CI", out_cert, fixed = TRUE)))
  expect_false(any(grepl("95% CI", out_cert, fixed = TRUE)))
})

test_that("edge rows of the bootstrap and certainty prints line up", {
  net <- .f_net()
  aligned <- function(out) {
    rows <- grep("→", out, value = TRUE)
    length(unique(regexpr("[", rows, fixed = TRUE))) == 1L
  }
  expect_true(aligned(capture.output(
    print(bootstrap_network(net, iter = 20, seed = 1)))))
  expect_true(aligned(capture.output(print(certainty(net)))))
})

test_that("a certainty group is not printed as a bootstrap with NA iterations", {
  grp <- .f_net(group = "Achiever")
  out <- capture.output(print(certainty(grp)))
  expect_true(any(grepl("Grouped Certainty", out, fixed = TRUE)))
  expect_false(any(grepl("NA iterations", out, fixed = TRUE)))
  out_boot <- capture.output(print(bootstrap_network(grp, iter = 10, seed = 1)))
  expect_true(any(grepl("Grouped Bootstrap  [2 groups | 10 iterations",
                        out_boot, fixed = TRUE)))
})

test_that("the bootstrap fast path refuses long data with a classed error", {
  expect_error(
    .precompute_per_sequence(data.frame(a = 1), "relative",
                             list(format = "long"), "A"),
    class = "nestimate_bad_input")
})

test_that("boot_glasso() refuses no more rows than variables by class", {
  wide <- srl_strategies[seq_len(4), seq_len(4)] + diag(4)
  expect_error(boot_glasso(wide, iter = 2, cs_iter = 1),
               class = "nestimate_bad_input")
})

test_that("boot_glasso() prints and labels what it computed", {
  bg <- boot_glasso(srl_strategies[, seq_len(4)], iter = 10, cs_iter = 10,
                    seed = 1)
  out <- capture.output(print(bg))
  expect_true(any(grepl("10 case-drop iterations", out, fixed = TRUE)))
  expect_false(any(grepl("per proportion", out, fixed = TRUE)))
  p <- plot(bg, type = "edge_diff")
  expect_match(p$labels$subtitle, "Red: column > row", fixed = TRUE)
  # the fill is the x-axis (column) edge minus the y-axis (row) edge
  cell <- subset(p$data, as.integer(edge1) != as.integer(edge2) &
                   !is.na(fill_val))
  means <- colMeans(bg$boot_edges)
  expect_equal(cell$fill_val,
               unname(means[as.character(cell$edge1)] -
                        means[as.character(cell$edge2)]))
})

test_that("the thresholded boot_glasso network zeroes both triangles", {
  pcor <- matrix(c(0, 0.3, 0.2, 0.3, 0, 0.1, 0.2, 0.1, 0), 3)
  ut <- .bg_upper_tri_indices(3)
  out <- .bg_threshold_network(pcor, c(0.1, -0.1, 0), c(0.5, 0.4, 0.2), ut)
  expect_equal(out, matrix(c(0, 0.3, 0, 0.3, 0, 0, 0, 0, 0), 3))
  expect_true(isSymmetric(out))
})

test_that("centrality_stability() warns by class when no measure varies", {
  net <- .f_net()
  expect_warning(
    res <- centrality_stability(
      net, measures = "foo", iter = 5,
      centrality_fn = function(m) list(foo = rep(1, nrow(m)))),
    class = "nestimate_constant_measure")
  # the all-constant branch still returns every requested name with cs = 0
  expect_identical(res$measures, "foo")
  expect_identical(unname(res$cs), 0)
})

test_that("centrality_stability() keeps its default measures (htna contract)", {
  expect_identical(eval(formals(centrality_stability)$measures),
                   c("InStrength", "OutStrength", "Betweenness"))
})

test_that("the stability prints show the threshold as given", {
  net <- .f_net()
  cs <- centrality_stability(net, iter = 5, seed = 1, threshold = 0.75,
                             drop_prop = c(0.1, 0.5))
  expect_true(any(grepl("threshold = 0.75)", capture.output(print(cs)),
                        fixed = TRUE)))
  cs_default <- centrality_stability(net, iter = 5, seed = 1,
                                     drop_prop = c(0.1, 0.5))
  expect_true(any(grepl("threshold = 0.7)", capture.output(print(cs_default)),
                        fixed = TRUE)))
})

test_that("a stability group with one measure prints", {
  grp <- .f_net(group = "Achiever")
  csg <- centrality_stability(grp, iter = 5, seed = 1, measures = "InStrength",
                              drop_prop = c(0.1, 0.5))
  out <- capture.output(print(csg))
  expect_true(any(grepl("InStrength", out, fixed = TRUE)))
  expect_true(any(grepl("^High", out)))
})

test_that("net_centrality() refuses a non-network by class", {
  expect_error(suppressMessages(net_centrality(diag(3))),
               class = "nestimate_bad_input")
})

test_that("plot(type = 'delta') says that `palette` is ignored for two networks", {
  cents <- suppressMessages(net_centrality(.f_net(group = "Achiever")))
  expect_warning(plot(cents, type = "delta", palette = "Dark2"),
                 class = "nestimate_ignored_argument")
  expect_no_warning(plot(cents, type = "delta"))
})

test_that("casedrop_reliability() says which proportions drop no case", {
  seqs <- data.frame(V1 = c("A", "B", "C", "A", "B"),
                     V2 = c("B", "C", "A", "C", "A"),
                     V3 = c("C", "A", "B", "B", "C"))
  net <- build_network(seqs, method = "frequency")
  expect_warning(
    res <- casedrop_reliability(net, iter = 10, drop_prop = c(0.1, 0.5),
                                seed = 1),
    class = "nestimate_drop_prop_skipped")
  expect_true(all(is.na(res$metrics$correlation[, "p0.1"])))
  expect_no_warning(
    casedrop_reliability(net, iter = 10, drop_prop = 0.5, seed = 1))
})

test_that("casedrop_reliability() print names its rows without `$`", {
  res <- casedrop_reliability(.f_net(), iter = 5, drop_prop = c(0.1, 0.5),
                              seed = 1)
  out <- capture.output(print(res))
  expect_false(any(grepl("$", out, fixed = TRUE)))
  expect_true(any(grepl("median|diff|", out, fixed = TRUE)))
  expect_false(any(grepl("MAD", out, fixed = TRUE)))
})

test_that("network_reliability() warns by class on mixed methods", {
  rel <- .f_net()
  freq <- build_network(group_regulation_long, method = "frequency",
                        action = "Action", actor = "Actor", time = "Time")
  expect_warning(res <- network_reliability(rel, freq, iter = 5, seed = 1),
                 class = "nestimate_mixed_methods")
  expect_identical(unique(res$summary$model), c("relative", "frequency"))
})

test_that("network_reliability() labels named, unnamed and grouped networks", {
  rel <- .f_net()
  grp <- .f_net(group = "Achiever")
  res <- network_reliability(mine = rel, grp, iter = 5, seed = 1)
  expect_identical(unique(res$summary$model), c("mine", "High", "Low"))
  out <- capture.output(print(res))
  expect_true(any(grepl("^  mine:$", out)))
  expect_error(network_reliability(iter = 5), class = "nestimate_bad_input")
})

test_that("vertex_compare() prints both replicate counts when they differ", {
  net <- .f_net()
  few <- vertex_bootstrap(net, iter = 10, seed = 1)
  many <- vertex_bootstrap(net, iter = 20, seed = 2)
  expect_true(any(grepl("10 and 20 replicates",
                        capture.output(print(vertex_compare(few, many))),
                        fixed = TRUE)))
  expect_true(any(grepl("10 replicates each",
                        capture.output(print(vertex_compare(few, few))),
                        fixed = TRUE)))
})

test_that("seeded resampling verbs leave the session's random stream alone", {
  net <- .f_net()
  set.seed(1)
  expected <- runif(1)
  untouched <- function(call) {
    set.seed(1)
    force(call)
    identical(runif(1), expected)
  }
  expect_true(untouched(bootstrap_network(net, iter = 10, seed = 5)))
  expect_true(untouched(centrality_stability(net, iter = 5, seed = 5,
                                             drop_prop = 0.5)))
  expect_true(untouched(casedrop_reliability(net, iter = 5, seed = 5,
                                             drop_prop = 0.5)))
  expect_true(untouched(network_reliability(net, iter = 5, seed = 5)))
  expect_true(untouched(vertex_bootstrap(net, iter = 5, seed = 5)))
})

test_that("bootstrap_network() counts association replicates it could not use", {
  # 14 of 15 values of `c` are equal, so many resamples leave it constant
  dat <- data.frame(a = srl_strategies[[1]][seq_len(15)],
                    b = srl_strategies[[2]][seq_len(15)],
                    c = c(rep(0, 14), 1))
  net <- suppressMessages(build_network(dat, method = "pcor"))
  expect_warning(
    boot <- suppressMessages(bootstrap_network(net, iter = 30, seed = 1)),
    class = "nestimate_resample_failed")
  expect_s3_class(boot, "net_bootstrap")
  whole <- suppressMessages(build_network(srl_strategies, method = "cor"))
  expect_no_warning(bootstrap_network(whole, iter = 10, seed = 1))
})
