# Numerical equivalence: build_mcml() matrix path vs manual aggregation
#
# build_mcml(mat, clusters, method, directed, compute_within) splits a weight
# matrix into (a) a between-cluster n_clusters x n_clusters matrix of
# aggregated weights and (b) per-cluster within-cluster submatrices. The math
# on the matrix path is elementary — elementwise aggregation over row/col
# masks, no post-processing. The sequence and edgelist paths keep `type`, but
# the matrix path is aggregation-only.
#
# This file reimplements that math with base R (outer / colSums / mapply)
# and compares cell-by-cell to build_mcml() across a random grid of
# (n_nodes, n_clusters, method) configurations.

set.seed(20260422)
N_CONFIGS <- 100L
TOL       <- 1e-10

# ---- Inline reference implementation ------------------------------------

.aggregate_ref <- function(w, method, n_possible) {
  w <- w[!is.na(w) & w != 0]
  if (length(w) == 0L) return(0)
  switch(method,
    "sum"     = sum(w),
    "mean"    = mean(w),
    "median"  = stats::median(w),
    "max"     = max(w),
    "min"     = min(w),
    "density" = if (!is.null(n_possible) && n_possible > 0) sum(w) / n_possible
                else sum(w) / length(w),
    "geomean" = { pos <- w[w > 0]; if (!length(pos)) 0 else exp(mean(log(pos))) },
    stop("unknown method")
  )
}

.mcml_ref <- function(mat, clusters, method) {
  cluster_names <- sort(unique(clusters))
  k <- length(cluster_names)
  node_by_cluster <- lapply(cluster_names, function(c_) which(clusters == c_))
  names(node_by_cluster) <- cluster_names

  between_raw <- matrix(0, k, k, dimnames = list(cluster_names, cluster_names))
  pair_grid <- expand.grid(i = seq_len(k), j = seq_len(k), KEEP.OUT.ATTRS = FALSE)
  between_raw[cbind(pair_grid$i, pair_grid$j)] <- mapply(
    function(i, j) {
      idx_i <- node_by_cluster[[i]]; idx_j <- node_by_cluster[[j]]
      .aggregate_ref(as.vector(mat[idx_i, idx_j]), method,
                     length(idx_i) * length(idx_j))
    },
    pair_grid$i, pair_grid$j
  )

  col_sums_b <- colSums(between_raw)
  total_b    <- sum(col_sums_b)
  between_inits <- if (total_b > 0) col_sums_b / total_b else rep(1 / k, k)
  names(between_inits) <- cluster_names

  within <- lapply(cluster_names, function(c_) {
    idx <- node_by_cluster[[c_]]
    sub <- mat[idx, idx, drop = FALSE]
    n_i <- length(idx)
    col_sums_w <- colSums(sub, na.rm = TRUE)
    total_w    <- sum(col_sums_w, na.rm = TRUE)
    inits_w <- if (n_i <= 1L) {
      setNames(1, dimnames(sub)[[1]])
    } else if (!is.na(total_w) && total_w > 0) {
      col_sums_w / total_w
    } else {
      setNames(rep(1 / n_i, n_i), dimnames(sub)[[1]])
    }
    list(weights = sub, inits = inits_w)
  })
  names(within) <- cluster_names

  list(between_weights = between_raw,
       between_inits   = between_inits,
       within          = within)
}

# ---- Config grid --------------------------------------------------------

methods <- c("sum", "mean", "median", "max", "min", "density", "geomean")

configs <- replicate(N_CONFIGS, list(
  n_nodes    = sample(6:15, 1L),
  n_clusters = sample(2:4, 1L),
  method     = sample(methods, 1L),
  seed       = sample.int(100000L, 1L)
), simplify = FALSE)

.gen_mat <- function(n, seed) {
  set.seed(seed)
  m <- matrix(abs(stats::rnorm(n * n, 0.4, 0.3)), n, n)
  m[stats::runif(n * n) < 0.25] <- 0
  dimnames(m) <- list(paste0("n", seq_len(n)), paste0("n", seq_len(n)))
  m
}

.gen_clusters <- function(n, k, seed) {
  set.seed(seed + 7L)
  repeat {
    cl <- paste0("C", sample.int(k, n, replace = TRUE))
    if (length(unique(cl)) == k) return(cl)
  }
}

# ---- Tests --------------------------------------------------------------

test_that("build_mcml between-cluster weights match manual aggregation", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    mat <- .gen_mat(cfg$n_nodes, cfg$seed)
    clu <- .gen_clusters(cfg$n_nodes, cfg$n_clusters, cfg$seed)
    mcml <- build_mcml(mat, clusters = clu, method = cfg$method,
                       directed = TRUE, compute_within = TRUE)
    ref  <- .mcml_ref(mat, clu, cfg$method)

    cluster_names <- sort(unique(clu))
    got <- mcml$macro$weights[cluster_names, cluster_names]
    exp <- ref$between_weights[cluster_names, cluster_names]
    max(abs(got - exp))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})

test_that("build_mcml between-cluster inits match manual column-sum ratios", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    mat <- .gen_mat(cfg$n_nodes, cfg$seed)
    clu <- .gen_clusters(cfg$n_nodes, cfg$n_clusters, cfg$seed)
    mcml <- build_mcml(mat, clusters = clu, method = cfg$method,
                       directed = TRUE, compute_within = TRUE)
    ref  <- .mcml_ref(mat, clu, cfg$method)

    cluster_names <- sort(unique(clu))
    got <- mcml$macro$inits[cluster_names]
    exp <- ref$between_inits[cluster_names]
    max(abs(got - exp))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})

test_that("build_mcml within-cluster weights match manual submatrices", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    mat <- .gen_mat(cfg$n_nodes, cfg$seed)
    clu <- .gen_clusters(cfg$n_nodes, cfg$n_clusters, cfg$seed)
    mcml <- build_mcml(mat, clusters = clu, method = cfg$method,
                       directed = TRUE, compute_within = TRUE)
    ref  <- .mcml_ref(mat, clu, cfg$method)

    cluster_names <- sort(unique(clu))
    per_cluster <- vapply(cluster_names, function(c_) {
      got <- mcml$clusters[[c_]]$weights
      exp <- ref$within[[c_]]$weights
      got <- got[rownames(exp), colnames(exp), drop = FALSE]
      max(abs(got - exp))
    }, numeric(1))
    max(per_cluster)
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})

test_that("build_mcml within-cluster inits match manual column-sum ratios", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    mat <- .gen_mat(cfg$n_nodes, cfg$seed)
    clu <- .gen_clusters(cfg$n_nodes, cfg$n_clusters, cfg$seed)
    mcml <- build_mcml(mat, clusters = clu, method = cfg$method,
                       directed = TRUE, compute_within = TRUE)
    ref  <- .mcml_ref(mat, clu, cfg$method)

    cluster_names <- sort(unique(clu))
    per_cluster <- vapply(cluster_names, function(c_) {
      got <- mcml$clusters[[c_]]$inits
      exp <- ref$within[[c_]]$inits
      got <- got[names(exp)]
      max(abs(got - exp))
    }, numeric(1))
    max(per_cluster)
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})

# ---- Delta report + validation dashboard emit --------------------------

test_that("build_mcml equivalence report (CSV + CVS JSON)", {
  skip_on_cran()

  report <- equiv_report()

  invisible(lapply(configs, function(cfg) {
    mat <- .gen_mat(cfg$n_nodes, cfg$seed)
    clu <- .gen_clusters(cfg$n_nodes, cfg$n_clusters, cfg$seed)
    mcml <- build_mcml(mat, clusters = clu, method = cfg$method,
                       directed = TRUE, compute_within = TRUE)
    ref  <- .mcml_ref(mat, clu, cfg$method)

    cluster_names <- sort(unique(clu))
    between_err <- abs(mcml$macro$weights[cluster_names, cluster_names] -
                       ref$between_weights[cluster_names, cluster_names])
    within_err  <- unlist(lapply(cluster_names, function(c_) {
      got <- mcml$clusters[[c_]]$weights
      exp <- ref$within[[c_]]$weights
      got <- got[rownames(exp), colnames(exp), drop = FALSE]
      as.vector(abs(got - exp))
    }))
    all_err <- c(as.vector(between_err), within_err)

    report$log(
      func = "build_mcml",
      config = sprintf("n=%d k=%d method=%s",
                       cfg$n_nodes, cfg$n_clusters, cfg$method),
      n_checked = length(all_err),
      n_failed  = sum(all_err >= TOL),
      max_abs_err    = max(all_err),
      mean_abs_err   = mean(all_err),
      median_abs_err = stats::median(all_err),
      p95_abs_err    = stats::quantile(all_err, 0.95, names = FALSE),
      reference = "manual base-R aggregation (mapply over cluster pairs)",
      notes = "matrix path is aggregation-only; no type post-processing"
    )
  }))

  report$write_csv("mcml")
  report$write_cvs("mcml")
  expect_true(length(report$rows) > 0L)
})
