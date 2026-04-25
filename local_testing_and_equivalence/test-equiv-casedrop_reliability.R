# Numerical equivalence: casedrop_reliability() vs manual case-drop rewrite
#
# casedrop_reliability() uses a "transition fast path" for method = relative /
# frequency / co_occurrence — it precomputes per-sequence count matrices once
# and resamples sessions via column sums rather than rebuilding the TNA each
# iteration. For equivalence we validate that this fast path produces the same
# per-iteration edge metrics as rebuilding the network from scratch on each
# subset.
#
# Reference per iteration:
#   1. sample keep_n session indices without replacement
#   2. rebuild net_sub = build_network(data[idx, ], method = "relative")
#   3. compute edge vector off-diagonal
#   4. compare to the original net's edge vector: correlation, mean/median/max
#      absolute deviation
#
# Seed alignment: set.seed(seed) once at the top of both paths, then iterate
# drop_prop outer / iter inner — matches casedrop_reliability's loop order.

set.seed(20260422)
N_CONFIGS <- 10L
ITER      <- 20L
DROP_PROP <- c(0.3, 0.6)
TOL       <- 1e-10

skip_if_not_installed("tna")

.gen_seqs <- function(n_actors, n_states, seq_length, seed) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- replicate(n_actors, sample(states, seq_length, replace = TRUE),
                    simplify = FALSE)
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

.off_diag_vec <- function(mat) {
  mask <- matrix(TRUE, nrow(mat), ncol(mat))
  diag(mask) <- FALSE
  as.vector(mat[mask])
}

.ref_casedrop <- function(data, method, iter, drop_prop, states, seed) {
  # Original network's edge vector (reference)
  orig <- build_network(data, method = method)
  states_ord <- orig$nodes$label
  orig_mat <- orig$weights[states_ord, states_ord]
  orig_vec <- .off_diag_vec(orig_mat)
  n_cases <- nrow(data)

  # Replicate the outer/inner loop order and RNG stream
  set.seed(seed)
  per_prop <- lapply(drop_prop, function(p) {
    n_drop <- floor(n_cases * p)
    keep_n <- n_cases - n_drop

    vapply(seq_len(iter), function(it) {
      idx <- sample(seq_len(n_cases), keep_n, replace = FALSE)
      sub <- data[idx, , drop = FALSE]
      sub_net <- build_network(sub, method = method)
      sub_mat <- sub_net$weights[states_ord, states_ord]
      sub_vec <- .off_diag_vec(sub_mat)
      diffs   <- abs(orig_vec - sub_vec)
      r <- if (stats::sd(sub_vec) > 0) {
        stats::cor(orig_vec, sub_vec, method = "spearman")
      } else NA_real_
      c(mean_abs_dev   = mean(diffs),
        median_abs_dev = stats::median(diffs),
        correlation    = r,
        max_abs_dev    = max(diffs))
    }, numeric(4))
  })
  names(per_prop) <- paste0("p", drop_prop)
  per_prop
}

configs <- replicate(N_CONFIGS, list(
  n_actors   = sample(c(40L, 60L, 100L), 1L),
  n_states   = sample(3:5, 1L),
  seq_length = sample(c(12L, 20L), 1L),
  seed       = sample.int(1e5, 1L)
), simplify = FALSE)


# ---- Transition fast path vs full rebuild ------------------------------

test_that("casedrop_reliability metrics match per-iteration full rebuild", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    net  <- build_network(seqs, method = "relative")

    cdr <- casedrop_reliability(
      net, iter = ITER, drop_prop = DROP_PROP,
      method = "spearman", seed = cfg$seed
    )
    ref <- .ref_casedrop(seqs, method = "relative",
                         iter = ITER, drop_prop = DROP_PROP,
                         states = net$nodes$label, seed = cfg$seed)

    # Compare each metric matrix (iter x n_prop) cell-by-cell
    metric_deltas <- vapply(
      c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev"),
      function(m) {
        got <- cdr$metrics[[m]]
        exp <- do.call(cbind, lapply(ref, function(mat) mat[m, ]))
        max(abs(got - exp), na.rm = TRUE)
      }, numeric(1)
    )
    max(metric_deltas)
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- CS-coefficient is derived deterministically from metrics ------------

test_that("casedrop_reliability CS-coefficient matches reference post-hoc derivation", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    net  <- build_network(seqs, method = "relative")

    cdr <- casedrop_reliability(
      net, iter = ITER, drop_prop = DROP_PROP,
      threshold = 0.5, certainty = 0.8,
      method = "spearman", seed = cfg$seed
    )
    # Reference CS from the same correlations matrix
    prop_above <- colMeans(cdr$metrics$correlation >= 0.5, na.rm = TRUE)
    qualifying <- which(prop_above >= 0.8)
    cs_ref <- if (length(qualifying)) max(DROP_PROP[qualifying]) else 0
    abs(cdr$cs - cs_ref)
  }, numeric(1))

  expect_true(all(deltas < TOL))
})


# ---- Delta report + validation dashboard emit --------------------------

test_that("casedrop_reliability equivalence report (CSV + CVS JSON)", {
  skip_on_cran()

  report <- equiv_report()

  invisible(lapply(configs, function(cfg) {
    seqs <- .gen_seqs(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    net  <- build_network(seqs, method = "relative")
    cdr <- casedrop_reliability(
      net, iter = ITER, drop_prop = DROP_PROP,
      method = "spearman", seed = cfg$seed
    )
    ref <- .ref_casedrop(seqs, method = "relative",
                         iter = ITER, drop_prop = DROP_PROP,
                         states = net$nodes$label, seed = cfg$seed)
    errs <- unlist(lapply(
      c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev"),
      function(m) {
        got <- cdr$metrics[[m]]
        exp <- do.call(cbind, lapply(ref, function(mat) mat[m, ]))
        as.vector(abs(got - exp))
      }
    ))
    errs <- errs[!is.na(errs)]

    report$log(
      func = "casedrop_reliability",
      config = sprintf("n_actors=%d n_states=%d seq_len=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length),
      n_checked = length(errs),
      n_failed  = sum(errs >= TOL),
      max_abs_err    = max(errs),
      mean_abs_err   = mean(errs),
      median_abs_err = stats::median(errs),
      p95_abs_err    = stats::quantile(errs, 0.95, names = FALSE),
      reference = "full build_network() rebuild per subset (no fast-path precompute)",
      notes = "transition fast path vs full rebuild, per-iteration metrics"
    )
  }))

  report$write_csv("casedrop_reliability")
  report$write_cvs("casedrop_reliability")
  expect_true(length(report$rows) > 0L)
})
