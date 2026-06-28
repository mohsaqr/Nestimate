# Tests for build_mcml_pc() (experimental MCML for psychometric networks)

# Block-structured Gaussian data: within-block r = .5, between A-B = .30,
# A-C = .10, B-C = .02 -- a planted ordering the aggregations must recover.
make_block_data <- function(n = 600, seed = 42) {
  set.seed(seed)
  p <- 9
  sigma <- matrix(0, p, p)
  blocks <- list(1:3, 4:6, 7:9)
  fill <- function(a, b, r) {
    sigma[blocks[[a]], blocks[[b]]] <<- r
    sigma[blocks[[b]], blocks[[a]]] <<- r
  }
  fill(1, 1, 0.5); fill(2, 2, 0.5); fill(3, 3, 0.5)
  fill(1, 2, 0.30); fill(1, 3, 0.10); fill(2, 3, 0.02)
  diag(sigma) <- 1
  z <- matrix(rnorm(n * p), n, p) %*% chol(sigma)
  df <- as.data.frame(z)
  names(df) <- c(paste0("a", 1:3), paste0("b", 1:3), paste0("c", 1:3))
  df
}

block_clusters <- function() {
  list(A = paste0("a", 1:3), B = paste0("b", 1:3), C = paste0("c", 1:3))
}

test_that("all three aggregations recover the planted block ordering", {
  df <- make_block_data()
  cl <- block_clusters()
  for (agg in c("average", "scaled", "loadings")) {
    fit <- build_mcml_pc(df, cl, aggregation = agg, method = "cor")
    W <- fit$macro$weights
    expect_true(W["A", "B"] > W["A", "C"],
                info = paste(agg, ": A-B should exceed A-C"))
    expect_true(W["A", "C"] > W["B", "C"],
                info = paste(agg, ": A-C should exceed B-C"))
  }
})

test_that("aggregation = 'composite' remains a scaled-score alias", {
  df <- make_block_data()
  cl <- block_clusters()

  scaled <- build_mcml_pc(df, cl, aggregation = "scaled", method = "cor")
  composite <- build_mcml_pc(df, cl, aggregation = "composite", method = "cor")

  expect_equal(composite$macro$weights, scaled$macro$weights)
  expect_identical(composite$meta$aggregation, "composite")
  expect_identical(composite$meta$scale, TRUE)
})

test_that("result structure is complete and undirected", {
  df <- make_block_data()
  fit <- build_mcml_pc(df, block_clusters(), aggregation = "scaled",
                       method = "pcor")
  expect_s3_class(fit, "mcml_pc")
  expect_true(inherits(fit$macro, "netobject"))
  expect_false(fit$macro$directed)
  expect_identical(fit$meta$directed, FALSE)
  expect_true(fit$meta$experimental)
  expect_identical(names(fit$clusters), c("A", "B", "C"))
  expect_true(all(vapply(fit$clusters,
                         function(n) inherits(n, "netobject"), logical(1))))
  expect_true(isSymmetric(unname(fit$macro$weights)))
})

test_that("netobject input reuses its data and method", {
  df <- make_block_data()
  net <- build_network(df, method = "glasso")
  fit <- build_mcml_pc(net, block_clusters(), aggregation = "scaled")
  expect_identical(fit$meta$method, "glasso")
})

test_that("average mode works without raw data; re-estimation refuses", {
  df <- make_block_data()
  net <- build_network(df, method = "pcor")
  net$data <- NULL
  fit <- build_mcml_pc(net, block_clusters(), aggregation = "average")
  expect_s3_class(fit, "mcml_pc")
  expect_identical(fit$meta$within, "subnetwork")
  expect_error(
    build_mcml_pc(net, block_clusters(), aggregation = "scaled"),
    "requires raw data"
  )
})

test_that("directed netobjects are rejected with guidance", {
  tnet <- build_network(
    data.frame(T1 = c("x", "y", "z"), T2 = c("y", "z", "x")),
    method = "relative"
  )
  expect_error(
    build_mcml_pc(tnet, list(G = c("x", "y"), H = "z")),
    "directed network"
  )
})

test_that("loadings are tidy diagnostics with weights summing to 1", {
  df <- make_block_data()
  fit <- build_mcml_pc(df, block_clusters(), aggregation = "loadings",
                       method = "cor")
  ld <- fit$loadings
  expect_s3_class(ld, "data.frame")
  expect_identical(
    names(ld),
    c("node", "cluster", "loading", "weight", "sign",
      "max_cross", "cross_cluster", "misfit")
  )
  expect_identical(nrow(ld), 9L)
  sums <- tapply(ld$weight, ld$cluster, sum)
  expect_equal(as.numeric(sums), rep(1, 3))
  expect_true(all(ld$weight >= 0))
  # well-assigned exchangeable blocks: no misfit, no flips
  expect_false(any(ld$misfit))
  expect_true(all(ld$sign == 1))
})

test_that("within = 'subnetwork' slices the node-level matrix", {
  df <- make_block_data()
  net <- build_network(df, method = "cor")
  fit <- build_mcml_pc(net, block_clusters(), within = "subnetwork")
  members <- block_clusters()$A
  expect_equal(fit$clusters$A$weights, net$weights[members, members])
})

test_that("singleton clusters get NULL within networks", {
  df <- make_block_data()
  cl <- list(A = paste0("a", 1:3), B = paste0("b", 1:3),
             C1 = "c1", C23 = c("c2", "c3"))
  # splitting true block C trips the misfit diagnostic (correctly):
  # c2 is as connected to c1 as to its assigned cluster
  fit <- suppressWarnings(
    build_mcml_pc(df, cl, aggregation = "scaled", method = "cor")
  )
  expect_null(fit$clusters$C1)
  expect_true(inherits(fit$clusters$C23, "netobject"))
})

test_that("membership vector input matches list input", {
  df <- make_block_data()
  vec <- setNames(rep(c("A", "B", "C"), each = 3), names(df))
  fit_vec <- build_mcml_pc(df, vec, aggregation = "average",
                           method = "cor")
  fit_list <- build_mcml_pc(df, block_clusters(), aggregation = "average",
                            method = "cor")
  expect_equal(fit_vec$macro$weights, fit_list$macro$weights)
})

test_that("cluster validation errors are clear", {
  df <- make_block_data()
  expect_error(build_mcml_pc(df, list(A = c("a1", "bogus"))),
               "not present in the network")
  expect_error(
    build_mcml_pc(df, list(A = paste0("a", 1:3), B = paste0("b", 1:3))),
    "without a cluster assignment"
  )
  dup <- list(A = c("a1", "a2"), B = c("a2", "a3"),
              REST = c(paste0("b", 1:3), paste0("c", 1:3)))
  expect_error(build_mcml_pc(df, dup), "exactly one cluster")
  expect_error(build_mcml_pc(df, list(ALL = names(df))), "At least 2")
})

test_that("print, summary, and plot methods work", {
  df <- make_block_data()
  fit <- build_mcml_pc(df, block_clusters(), aggregation = "average",
                       method = "cor")
  expect_output(print(fit), "experimental")
  s <- summary(fit)
  expect_identical(names(s), c("from", "to", "weight"))
  expect_identical(nrow(s), 3L)
  expect_s3_class(plot(fit), "ggplot")
})

test_that("composite and loadings agree when loadings are uniform", {
  # Exchangeable within-block correlations -> near-equal loadings ->
  # loadings composite ~= plain composite
  df <- make_block_data()
  cl <- block_clusters()
  f1 <- build_mcml_pc(df, cl, aggregation = "scaled", method = "cor")
  f2 <- build_mcml_pc(df, cl, aggregation = "loadings", method = "cor")
  expect_true(max(abs(f1$macro$weights - f2$macro$weights)) < 0.01)
})

# ---- comprehensive feature battery ----

test_that("misassigned items are flagged as misfit with a warning", {
  df <- make_block_data()
  # put a1 (an A item) into cluster B on purpose
  bad <- list(A = c("a2", "a3"), B = c("a1", paste0("b", 1:3)),
              C = paste0("c", 1:3))
  expect_warning(
    fit <- build_mcml_pc(df, bad, aggregation = "loadings",
                         method = "cor"),
    "misassignment"
  )
  ld <- fit$loadings
  expect_true(ld$misfit[ld$node == "a1"])
  expect_identical(ld$cross_cluster[ld$node == "a1"], "A")
  expect_identical(fit$meta$n_misfit, 1L)
})

test_that("reverse-keyed items are flipped and the macro recovers", {
  df <- make_block_data()
  df_rev <- df
  df_rev$a1 <- -df_rev$a1   # reverse-key one item
  cl <- block_clusters()

  expect_warning(
    fit_rev <- build_mcml_pc(df_rev, cl, aggregation = "loadings",
                             method = "cor"),
    "Reverse-keyed"
  )
  expect_identical(fit_rev$loadings$sign[fit_rev$loadings$node == "a1"], -1)
  expect_identical(fit_rev$meta$n_flipped, 1L)

  # with the flip, the macro should match the unreversed analysis closely
  fit_orig <- build_mcml_pc(df, cl, aggregation = "loadings",
                            method = "cor")
  expect_true(max(abs(fit_rev$macro$weights - fit_orig$macro$weights))
              < 0.02)

  # without sign correction the A-B edge is attenuated
  suppressWarnings(
    fit_unsigned <- build_mcml_pc(df_rev, cl, aggregation = "loadings",
                                  method = "cor", signed = FALSE)
  )
  expect_true(fit_unsigned$macro$weights["A", "B"] <
                fit_rev$macro$weights["A", "B"])
})

test_that("rv and canonical aggregations recover the planted ordering", {
  df <- make_block_data()
  cl <- block_clusters()
  for (agg in c("escoufier", "cancor")) {
    fit <- build_mcml_pc(df, cl, aggregation = agg)
    W <- fit$macro$weights
    expect_true(W["A", "B"] > W["A", "C"], info = agg)
    expect_true(W["A", "C"] > W["B", "C"], info = agg)
    expect_true(all(W >= 0 & W <= 1 + 1e-10), info = agg)
    expect_true(all(diag(W) == 1), info = agg)
    expect_true(is.na(fit$meta$method))
  }
})

test_that("canonical correlation upper-bounds the composite correlation", {
  df <- make_block_data()
  cl <- block_clusters()
  comp <- build_mcml_pc(df, cl, aggregation = "scaled",
                        method = "cor")$macro$weights
  can <- build_mcml_pc(df, cl, aggregation = "cancor")$macro$weights
  off <- upper.tri(comp)
  expect_true(all(can[off] >= abs(comp[off]) - 1e-10))
})

test_that("rv and canonical refuse weight-only input", {
  df <- make_block_data()
  net <- build_network(df, method = "cor")
  net$data <- NULL
  expect_error(build_mcml_pc(net, block_clusters(), aggregation = "escoufier"),
               "requires raw data")
})

test_that("composites tolerate missing data", {
  df <- make_block_data()
  set.seed(9)
  holes <- cbind(sample(nrow(df), 40), sample(ncol(df), 40, replace = TRUE))
  df[holes] <- NA
  suppressMessages(
    fit <- build_mcml_pc(df, block_clusters(), aggregation = "scaled",
                         method = "cor")
  )
  expect_s3_class(fit, "mcml_pc")
  W <- fit$macro$weights
  expect_true(all(is.finite(W)))
  expect_true(W["A", "B"] > W["B", "C"])
})

test_that("polychoric estimation works on ordinal data", {
  skip_if_not_installed("lavaan")
  df <- make_block_data()
  # discretize to 5-point Likert
  likert <- as.data.frame(lapply(df, function(x) {
    as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, 0.2)),
                   include.lowest = TRUE))
  }))
  fit <- build_mcml_pc(likert, block_clusters(), aggregation = "loadings",
                       method = "pcor", cor_method = "polychoric")
  expect_s3_class(fit, "mcml_pc")
  expect_identical(fit$meta$cor_method, "polychoric")
  W <- fit$macro$weights
  expect_true(W["A", "B"] > W["B", "C"])
})

test_that("bootstrap machinery works on the composite macro network", {
  df <- make_block_data(n = 200)
  fit <- build_mcml_pc(df, block_clusters(), aggregation = "scaled",
                       method = "cor")
  boot <- bootstrap_network(fit$macro, iter = 30, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  vb <- vertex_bootstrap(fit$macro, iter = 30, seed = 1)
  expect_s3_class(vb, "net_vertex_bootstrap")
  expect_false(vb$directed)
})

test_that("loading_stability returns tidy CIs and is reproducible", {
  df <- make_block_data(n = 200)
  fit <- build_mcml_pc(df, block_clusters(), aggregation = "loadings",
                       method = "cor")
  ls1 <- loading_stability(fit, iter = 20, seed = 7)
  ls2 <- loading_stability(fit, iter = 20, seed = 7)
  expect_s3_class(ls1, "pc_loading_stability")
  expect_identical(ls1$boot_weights, ls2$boot_weights)
  expect_identical(
    names(ls1$summary),
    c("node", "cluster", "weight", "boot_mean", "boot_sd",
      "ci_lower", "ci_upper", "sign_flips")
  )
  expect_identical(nrow(ls1$summary), 9L)
  expect_true(all(ls1$summary$ci_lower <= ls1$summary$ci_upper))
  expect_output(print(ls1), "Stability")
  expect_s3_class(plot(ls1), "ggplot")

  # refuses objects without data
  fit$data <- NULL
  expect_error(loading_stability(fit), "no raw data")
})

test_that("singleton clusters are not flagged misfit and carry weight 1", {
  df <- make_block_data()
  cl <- list(A = paste0("a", 1:3), B = paste0("b", 1:3),
             C1 = "c1", C23 = c("c2", "c3"))
  fit <- suppressWarnings(
    build_mcml_pc(df, cl, aggregation = "scaled", method = "cor")
  )
  ld <- fit$loadings
  expect_false(ld$misfit[ld$node == "c1"])
  expect_identical(ld$weight[ld$node == "c1"], 1)
})

# ---- weighting schemes ----

test_that("all five weightings recover the planted block ordering", {
  df <- make_block_data()
  cl <- block_clusters()
  for (w in c("equal", "strength", "eigen", "pca", "factor")) {
    fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = w, method = "cor")
    W <- fit$macro$weights
    expect_true(W["A", "B"] > W["A", "C"],
                info = paste(w, ": A-B should exceed A-C"))
    expect_true(W["A", "C"] > W["B", "C"],
                info = paste(w, ": A-C should exceed B-C"))
  }
})

test_that("aggregation = 'loadings' is composite with strength weighting", {
  df <- make_block_data()
  cl <- block_clusters()
  f_load <- build_mcml_pc(df, cl, aggregation = "loadings",
                          method = "cor")
  f_str <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = "strength", method = "cor")
  expect_identical(f_load$macro$weights, f_str$macro$weights)
  expect_identical(f_load$meta$weighting, "strength")
  f_comp <- build_mcml_pc(df, cl, aggregation = "scaled",
                          method = "cor")
  expect_identical(f_comp$meta$weighting, "equal")
})

test_that("every weighting gives non-negative weights summing to 1", {
  df <- make_block_data()
  cl <- block_clusters()
  for (w in c("equal", "strength", "eigen", "pca", "factor")) {
    fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = w, method = "cor")
    ld <- fit$loadings
    expect_true(all(ld$weight >= 0), info = w)
    sums <- tapply(ld$weight, ld$cluster, sum)
    expect_true(all(abs(as.numeric(sums) - 1) < 1e-12),
                info = paste(w, ": weights should sum to 1 per cluster"))
  }
})

test_that("data-driven weightings diverge from equal on a weak item", {
  d <- make_block_data()
  # weaken a3's connections to its own cluster
  d$a3 <- 0.5 * d$a3 + rnorm(nrow(d), sd = 0.9)
  cl <- block_clusters()
  for (w in c("strength", "pca", "factor")) {
    fit <- suppressWarnings(
      build_mcml_pc(d, cl, aggregation = "scaled", weighting = w,
                    method = "cor")
    )
    ld <- fit$loadings
    wa <- setNames(ld$weight[ld$cluster == "A"], ld$node[ld$cluster == "A"])
    expect_true(max(abs(wa - 1 / 3)) > 0.02,
                info = paste(w, ": should diverge from equal weights"))
    expect_true(names(which.min(wa)) == "a3",
                info = paste(w, ": weak item a3 should get the smallest",
                             "weight"))
  }
})

test_that("pca and factor weightings recover a reverse-keyed item", {
  df <- make_block_data()
  cl <- block_clusters()
  df_rev <- df
  df_rev$a1 <- -df_rev$a1
  for (w in c("pca", "factor")) {
    fit_orig <- suppressWarnings(
      build_mcml_pc(df, cl, aggregation = "scaled", weighting = w,
                    method = "cor")
    )
    fit_rev <- suppressWarnings(
      build_mcml_pc(df_rev, cl, aggregation = "scaled", weighting = w,
                    method = "cor")
    )
    ld <- fit_rev$loadings
    expect_true(ld$sign[ld$node == "a1"] == -1,
                info = paste(w, ": reversed item should carry sign -1"))
    expect_true(
      abs(fit_rev$macro$weights["A", "B"] -
            fit_orig$macro$weights["A", "B"]) < 0.05,
      info = paste(w, ": sign-corrected macro A-B should match unreversed")
    )
  }
})

test_that("factor weighting on a 2-item cluster warns and falls back", {
  df <- make_block_data()
  # a3 lands in B so all 9 nodes stay assigned; that misassignment warns
  # separately, hence the suppressWarnings wrapper around expect_warning
  cl2 <- list(A = c("a1", "a2"),
              B = c("a3", paste0("b", 1:3)),
              C = paste0("c", 1:3))
  suppressWarnings(expect_warning(
    fit <- build_mcml_pc(df, cl2, aggregation = "scaled",
                         weighting = "factor", method = "cor"),
    "fewer than 3"
  ))
  expect_s3_class(fit, "mcml_pc")
  expect_identical(fit$meta$weighting, "factor")
  ld <- fit$loadings
  expect_equal(sum(ld$weight[ld$cluster == "A"]), 1)
})

test_that("factor weighting is deterministic across identical calls", {
  df <- make_block_data()
  cl <- block_clusters()
  f1 <- build_mcml_pc(df, cl, aggregation = "scaled",
                      weighting = "factor", method = "cor")
  f2 <- build_mcml_pc(df, cl, aggregation = "scaled",
                      weighting = "factor", method = "cor")
  expect_identical(f1$macro$weights, f2$macro$weights)
})

test_that("id_col drops identifier columns from data.frame input", {
  data(group_regulation_long)
  prof <- convert_sequence_format(group_regulation_long, action = "Action",
                                  id_col = "Actor", format = "frequency")
  cl <- list(
    Regulation = c("monitor", "plan", "coregulate", "adapt"),
    Social     = c("cohesion", "consensus", "discuss", "emotion", "synthesis")
  )
  fit <- suppressWarnings(
    build_mcml_pc(prof, cl, aggregation = "scaled", method = "pcor",
                  id_col = c("rid", "Actor"))
  )
  expect_s3_class(fit, "mcml_pc")
  expect_identical(fit$meta$n_nodes, 9L)
  expect_error(
    build_mcml_pc(prof, cl, id_col = "bogus"),
    "not found in data"
  )
})

# ---- expanded weighting schemes ----

test_that("the five new weightings recover the planted block ordering", {
  df <- make_block_data()
  cl <- block_clusters()
  for (w in c("closeness", "betweenness", "expected_influence",
              "specificity", "item_total")) {
    fit <- suppressWarnings(
      build_mcml_pc(df, cl, aggregation = "scaled",
                    weighting = w, method = "cor")
    )
    W <- fit$macro$weights
    expect_true(W["A", "B"] > W["B", "C"],
                info = paste(w, ": A-B should exceed B-C"))
    ld <- fit$loadings
    expect_true(all(ld$weight >= 0), info = w)
    sums <- tapply(ld$weight, ld$cluster, sum)
    expect_true(all(abs(as.numeric(sums) - 1) < 1e-12),
                info = paste(w, ": weights should sum to 1 per cluster"))
  }
})

test_that("expected_influence equals strength without negative edges", {
  df <- make_block_data()
  cl <- block_clusters()
  f_ei <- build_mcml_pc(df, cl, aggregation = "scaled",
                        weighting = "expected_influence", method = "cor")
  f_str <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = "strength", method = "cor")
  expect_equal(f_ei$macro$weights, f_str$macro$weights)
})

test_that("specificity downweights a misassigned item to zero", {
  df <- make_block_data()
  # a1 belongs to A but is placed in B: its misfit margin is <= 0
  bad <- list(A = c("a2", "a3"), B = c("a1", paste0("b", 1:3)),
              C = paste0("c", 1:3))
  fit <- suppressWarnings(
    build_mcml_pc(df, bad, aggregation = "scaled",
                  weighting = "specificity", method = "cor")
  )
  ld <- fit$loadings
  expect_equal(ld$weight[ld$node == "a1"], 0)
})

test_that("betweenness falls back to equal weights with a warning", {
  df <- make_block_data()
  # exchangeable blocks -> complete within-cluster cor networks ->
  # no item ever lies on a shortest path -> all-zero betweenness; the
  # fallback warns once per cluster, hence the suppressWarnings wrapper
  suppressWarnings(expect_warning(
    fit <- build_mcml_pc(df, block_clusters(), aggregation = "scaled",
                         weighting = "betweenness", method = "cor"),
    "betweenness"
  ))
  ld <- fit$loadings
  expect_true(all(abs(ld$weight - 1 / 3) < 1e-12))
})

test_that("custom named weight vectors are honored and validated", {
  df <- make_block_data()
  cl <- block_clusters()
  w_vec <- setNames(rep(1, 9), names(df))
  w_vec["a3"] <- 0
  fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                       weighting = w_vec, method = "cor")
  ld <- fit$loadings
  expect_equal(ld$weight[ld$node == "a3"], 0)
  expect_equal(ld$weight[ld$node == "a1"], 0.5)
  expect_identical(fit$meta$weighting, "custom (vector)")
  expect_error(
    build_mcml_pc(df, cl, aggregation = "scaled",
                  weighting = w_vec[names(w_vec) != "a3"],
                  method = "cor"),
    "lacks entries"
  )
})

test_that("custom weighting functions are honored and validated", {
  df <- make_block_data()
  cl <- block_clusters()
  fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                       weighting = function(Wb, db, nodes) seq_along(nodes),
                       method = "cor")
  expect_s3_class(fit, "mcml_pc")
  expect_identical(fit$meta$weighting, "custom (function)")
  ld <- fit$loadings
  expect_equal(as.numeric(tapply(ld$weight, ld$cluster, sum)), rep(1, 3))
  expect_error(
    build_mcml_pc(df, cl, aggregation = "scaled",
                  weighting = function(Wb, db, nodes) 1,
                  method = "cor"),
    "finite numeric"
  )
})

test_that("item_total recovers a reverse-keyed item", {
  # The rest-mean of a 3-item cluster is half made of the reversed item,
  # so corrected item-total correlations of the good items collapse to
  # ~0 there. A 5-item block keeps the rest-mean informative -- the
  # setting where item-total weighting is defensible at all.
  set.seed(7)
  n <- 600
  sigma <- matrix(0.3, 8, 8)
  sigma[1:5, 1:5] <- 0.5
  sigma[6:8, 6:8] <- 0.5
  diag(sigma) <- 1
  z <- matrix(rnorm(n * 8), n, 8) %*% chol(sigma)
  df <- as.data.frame(z)
  names(df) <- c(paste0("a", 1:5), paste0("b", 1:3))
  cl <- list(A = paste0("a", 1:5), B = paste0("b", 1:3))
  df_rev <- df
  df_rev$a1 <- -df_rev$a1
  fit_orig <- build_mcml_pc(df, cl, aggregation = "scaled",
                            weighting = "item_total", method = "cor")
  fit_rev <- suppressWarnings(
    build_mcml_pc(df_rev, cl, aggregation = "scaled",
                  weighting = "item_total", method = "cor")
  )
  ld <- fit_rev$loadings
  expect_true(ld$sign[ld$node == "a1"] == -1)
  expect_true(
    abs(fit_rev$macro$weights["A", "B"] -
          fit_orig$macro$weights["A", "B"]) < 0.05
  )
})

test_that("item_total recovers a reverse-keyed item even in a 3-item cluster", {
  # Regression: before the eigen-sign pre-orientation, the reversed member
  # contaminated the rest-mean and collapsed good items' item-totals to ~0
  # in 3-item clusters (recovery only worked for larger blocks).
  df <- make_block_data()
  df_rev <- df
  df_rev$a1 <- -df_rev$a1
  cl <- block_clusters()

  fit_rev <- suppressWarnings(
    build_mcml_pc(df_rev, cl, aggregation = "scaled",
                  weighting = "item_total", method = "cor")
  )
  fit_orig <- suppressWarnings(
    build_mcml_pc(df, cl, aggregation = "scaled",
                  weighting = "item_total", method = "cor")
  )
  expect_identical(fit_rev$loadings$sign[fit_rev$loadings$node == "a1"], -1)
  expect_true(max(abs(fit_rev$macro$weights - fit_orig$macro$weights))
              < 0.05)
})

# ---- factor extraction methods (fa_method) ----

test_that("all fa_method extractors run and agree on clean 1-factor blocks", {
  df <- make_block_data()
  cl <- block_clusters()
  fits <- lapply(c("ml", "paf", "minres"), function(fm) {
    build_mcml_pc(df, cl, aggregation = "scaled", weighting = "factor",
                  method = "cor", fa_method = fm)
  })
  names(fits) <- c("ml", "paf", "minres")
  for (fm in names(fits)) {
    W <- fits[[fm]]$macro$weights
    expect_true(W["A", "B"] > W["B", "C"], info = fm)
    sums <- tapply(fits[[fm]]$loadings$weight, fits[[fm]]$loadings$cluster,
                   sum)
    expect_equal(as.numeric(sums), rep(1, 3))
    expect_identical(fits[[fm]]$meta$fa_method, fm)
  }
  # extraction methods agree on exchangeable unidimensional blocks
  expect_true(max(abs(fits$ml$macro$weights - fits$paf$macro$weights))
              < 0.02)
  expect_true(max(abs(fits$paf$macro$weights - fits$minres$macro$weights))
              < 0.005)
})

test_that("cfa extraction works and matches ml for continuous items", {
  skip_if_not_installed("lavaan")
  df <- make_block_data()
  cl <- block_clusters()
  f_cfa <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = "factor", method = "cor",
                         fa_method = "cfa")
  f_ml <- build_mcml_pc(df, cl, aggregation = "scaled",
                        weighting = "factor", method = "cor",
                        fa_method = "ml")
  expect_identical(f_cfa$meta$fa_method, "cfa")
  expect_true(max(abs(f_cfa$macro$weights - f_ml$macro$weights)) < 0.01)
})

test_that("polychoric CFA (categorical factor model) runs on ordinal data", {
  skip_if_not_installed("lavaan")
  df <- make_block_data()
  likert <- as.data.frame(lapply(df, function(x) {
    as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, 0.2)),
                   include.lowest = TRUE))
  }))
  fit <- build_mcml_pc(likert, block_clusters(), aggregation = "scaled",
                       weighting = "factor", fa_method = "cfa",
                       cor_method = "polychoric", method = "pcor")
  expect_identical(fit$meta$fa_method, "cfa")
  expect_identical(fit$meta$cor_method, "polychoric")
  W <- fit$macro$weights
  expect_true(W["A", "B"] > W["B", "C"])
})

test_that("fa_method is validated and recorded only for factor weighting", {
  df <- make_block_data()
  cl <- block_clusters()
  expect_error(
    build_mcml_pc(df, cl, aggregation = "scaled", weighting = "factor",
                  method = "cor", fa_method = "bogus"),
    "should be one of"
  )
  fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                       weighting = "equal", method = "cor")
  expect_true(is.na(fit$meta$fa_method))
})

test_that("lavaan arguments pass through ... verbatim", {
  skip_if_not_installed("lavaan")
  df <- make_block_data()
  likert <- as.data.frame(lapply(df, function(x) {
    as.integer(cut(x, breaks = quantile(x, probs = seq(0, 1, 0.2)),
                   include.lowest = TRUE))
  }))
  cl <- block_clusters()

  # estimator choice reaches lavaan: WLSMV vs ULSMV give different
  # (but close) loadings on ordered items
  f_wlsmv <- build_mcml_pc(likert, cl, aggregation = "scaled",
                           weighting = "factor", fa_method = "cfa",
                           cor_method = "polychoric", method = "cor",
                           estimator = "WLSMV")
  f_ulsmv <- build_mcml_pc(likert, cl, aggregation = "scaled",
                           weighting = "factor", fa_method = "cfa",
                           cor_method = "polychoric", method = "cor",
                           estimator = "ULSMV")
  expect_identical(f_wlsmv$meta$fa_args, list(estimator = "WLSMV"))
  # WLSMV vs ULSMV ordered-CFA *loadings* coincide under Windows LAPACK (the
  # estimator there changes only SEs/fit, not point estimates), so the
  # "estimators differ numerically" checks below are platform-sensitive.
  # The argument pass-through above is verified on every platform.
  skip_on_os("windows")
  expect_false(identical(f_wlsmv$loadings$weight, f_ulsmv$loadings$weight))
  expect_true(max(abs(f_wlsmv$macro$weights - f_ulsmv$macro$weights)) < 0.05)
})

test_that("protected fa_args are ignored with a warning", {
  df <- make_block_data()
  cl <- block_clusters()
  expect_warning(
    fit <- build_mcml_pc(df, cl, aggregation = "scaled",
                         weighting = "factor", method = "cor",
                         factors = 3),
    "managed internally"
  )
  # result identical to the unmolested one-factor fit
  ref <- build_mcml_pc(df, cl, aggregation = "scaled",
                       weighting = "factor", method = "cor")
  expect_equal(fit$macro$weights, ref$macro$weights)
})

test_that("paf convergence controls pass through fa_args", {
  df <- make_block_data()
  cl <- block_clusters()
  f1 <- build_mcml_pc(df, cl, aggregation = "scaled",
                      weighting = "factor", method = "cor",
                      fa_method = "paf", tol = 1e-10)
  f2 <- build_mcml_pc(df, cl, aggregation = "scaled",
                      weighting = "factor", method = "cor",
                      fa_method = "paf")
  expect_true(max(abs(f1$macro$weights - f2$macro$weights)) < 1e-4)
  expect_true(is.na(build_mcml_pc(df, cl, aggregation = "scaled",
                                  weighting = "equal",
                                  method = "cor")$meta$fa_method))
})

test_that("stray ... arguments error unless weighting = 'factor'", {
  df <- make_block_data()
  cl <- block_clusters()
  expect_error(
    build_mcml_pc(df, cl, aggregation = "scaled", weighting = "equal",
                  method = "cor", missing = "fiml"),
    "Unused argument"
  )
  expect_error(
    build_mcml_pc(df, cl, aggregation = "scaled", weighting = "factor",
                  method = "cor", fa_method = "ml",
                  estimator = "WLSMV"),
    "applies only"
  )
})

# ---- as_networks(): promote an mcml_pc to a netobject_group ----

test_that("as_networks.mcml_pc returns a netobject_group of macro + clusters", {
  fit  <- build_mcml_pc(make_block_data(), block_clusters(),
                        aggregation = "scaled", method = "cor")
  nets <- as_networks(fit)

  expect_true(inherits(nets, "netobject_group"),
              info = "as_networks must return a netobject_group")
  expect_identical(names(nets), c("macro", "A", "B", "C"))
  expect_length(nets, 4L)  # 1 macro + 3 clusters
})

test_that("as_networks preserves psychometric semantics (undirected, estimator, data)", {
  fit  <- build_mcml_pc(make_block_data(), block_clusters(),
                        aggregation = "scaled", method = "cor")
  nets <- as_networks(fit)

  expect_true(inherits(nets$macro, "netobject"))
  expect_false(nets$macro$directed)          # psychometric nets are undirected
  expect_false(nets$A$directed)
  expect_identical(nets$macro$method, "cor")
  expect_identical(nets$A$method, "cor")
  expect_false(is.null(nets$macro$data))     # data carried through, not dropped
  # assembled, not re-wrapped: the macro netobject is the same object
  expect_identical(nets$macro$weights, fit$macro$weights)
})

test_that("as_networks output flows into downstream verbs (centrality)", {
  fit  <- build_mcml_pc(make_block_data(), block_clusters(),
                        aggregation = "scaled", method = "cor")
  nets <- as_networks(fit)
  ce   <- net_centrality(nets$macro)
  expect_true(is.data.frame(ce))
  expect_identical(sort(rownames(ce)), c("A", "B", "C"))
})

test_that("as_networks drops singleton clusters with a warning", {
  df <- make_block_data()
  # C becomes a singleton (one node) -> no within-network -> dropped
  cl <- list(A = paste0("a", 1:3), B = paste0("b", 1:3),
             Bx = paste0("c", 1:2), C = "c3")
  fit <- build_mcml_pc(df, cl, aggregation = "scaled", method = "cor")
  expect_warning(nets <- as_networks(fit),
                 "singleton clusters")
  expect_false("C" %in% names(nets))
  expect_true(all(c("macro", "A", "B", "Bx") %in% names(nets)))
})

test_that("as_networks.default passes through a group and errors otherwise", {
  fit  <- build_mcml_pc(make_block_data(), block_clusters(),
                        aggregation = "scaled", method = "cor")
  nets <- as_networks(fit)
  expect_identical(as_networks(nets), nets)             # passthrough
  expect_error(as_networks(list(a = 1)), "Cannot convert")
})

test_that("build_mcml_pc accepts clusters as a two-column data.frame", {
  set.seed(1)
  items <- as.data.frame(matrix(round(runif(200 * 6, 1, 5)), ncol = 6))
  names(items) <- paste0("v", 1:6)

  cl_list <- list(A = c("v1", "v2", "v3"), B = c("v4", "v5", "v6"))
  cl_df   <- data.frame(item  = paste0("v", 1:6),
                        group = rep(c("A", "B"), each = 3))

  f_list <- build_mcml_pc(items, cl_list, aggregation = "scaled", method = "cor")
  f_df   <- build_mcml_pc(items, cl_df,   aggregation = "scaled", method = "cor")

  # df read by position -> same partition -> identical macro network.
  expect_equal(f_df$macro$weights, f_list$macro$weights)
  expect_setequal(names(f_df$clusters), names(f_list$clusters))

  # column names are irrelevant (read by position).
  cl_df2 <- stats::setNames(cl_df, c("node", "cluster"))
  expect_no_error(build_mcml_pc(items, cl_df2, aggregation = "scaled", method = "cor"))

  # a one-column data.frame is rejected with a clear message.
  expect_error(build_mcml_pc(items, data.frame(x = paste0("v", 1:6))),
               "at least two columns")
})
