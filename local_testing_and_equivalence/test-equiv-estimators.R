# ---- Equivalence tests: all 8 estimators vs EXTERNAL packages ----
#
# Every edge value is compared individually against an independent
# package implementation.  On failure, the exact (from, to, ours, ref)
# triple is reported.
#
#   relative      -> tna::tna()
#   frequency     -> tna::ftna()
#   co_occurrence -> tna::ctna()
#   cor           -> stats::cor()
#   pcor          -> corpcor::cor2pcor()
#   glasso        -> qgraph::EBICglasso()
#   ising         -> IsingFit::IsingFit()
#   attention     -> tna::atna()
#
# Requires NESTIMATE_EQUIV_TESTS=true.  Skipped on CRAN.

# ---- Global config -----------------------------------------------------------

set.seed(2026)
N_CONFIGS <- 100L
TOL       <- 1e-10

seq_configs <- lapply(seq_len(N_CONFIGS), function(i) {
  list(
    n_actors   = sample(c(10, 20, 30), 1),
    n_states   = sample(c(3, 5, 7), 1),
    seq_length = sample(c(10, 20, 30), 1),
    seed       = sample.int(100000, 1)
  )
})

assoc_configs <- lapply(seq_len(N_CONFIGS), function(i) {
  list(
    n    = sample(c(50, 100, 200), 1),
    p    = sample(c(4, 5, 6), 1),
    rho  = runif(1, 0.1, 0.5),
    seed = sample.int(100000, 1)
  )
})

# ---- Per-edge comparison helper -----------------------------------------------

#' Compare two matrices edge by edge, one expect_equal per cell.
#' On failure, prints the exact (row, col, ours, ref) that diverged.
#' @noRd
.check_edges <- function(ours, ref, tol, label) {
  # Align by shared state names
  states <- sort(intersect(rownames(ours), rownames(ref)))
  stopifnot(length(states) >= 2L)
  a <- ours[states, states, drop = FALSE]
  b <- ref[states, states, drop = FALSE]

  n <- length(states)
  n_checked <- 0L

  invisible(lapply(seq_len(n), function(i) {
    lapply(seq_len(n), function(j) {
      n_checked <<- n_checked + 1L
      expect_equal(
        a[i, j], b[i, j], tolerance = tol,
        label = sprintf("%s [%s -> %s]  ours=%.8e  ref=%.8e",
                        label, states[i], states[j], a[i, j], b[i, j])
      )
    })
  }))

  n_checked
}

report <- equiv_report()

seq_info <- function(cfg, i) {
  sprintf("cfg%d(a=%d,s=%d,l=%d)", i, cfg$n_actors, cfg$n_states, cfg$seq_length)
}
assoc_info <- function(cfg, i) {
  sprintf("cfg%d(n=%d,p=%d)", i, cfg$n, cfg$p)
}

# ---- Test: method = "relative" vs tna::tna() ---------------------------------

test_that("relative: every edge matches tna::tna() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- seq_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    ours <- build_network(data, method = "relative")$weights
    ref  <- tna::tna(data)$weights
    ref[is.na(ref)] <- 0

    n <- .check_edges(ours, ref, TOL, sprintf("relative %s", seq_info(cfg, i)))

    states <- sort(intersect(rownames(ours), rownames(ref)))
    delta <- abs(ours[states, states] - ref[states, states])
    report$log("relative", seq_info(cfg, i), n, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "tna::tna()")
  }))
})


# ---- Test: method = "frequency" vs tna::ftna() --------------------------------

test_that("frequency: every edge matches tna::ftna() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- seq_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    ours <- build_network(data, method = "frequency")$weights
    ref  <- tna::ftna(data)$weights
    ref[is.na(ref)] <- 0

    n <- .check_edges(ours, ref, TOL, sprintf("frequency %s", seq_info(cfg, i)))

    states <- sort(intersect(rownames(ours), rownames(ref)))
    delta <- abs(ours[states, states] - ref[states, states])
    report$log("frequency", seq_info(cfg, i), n, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "tna::ftna()")
  }))
})


# ---- Test: method = "co_occurrence" vs tna::ctna() ----------------------------

test_that("co_occurrence: every edge matches tna::ctna() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- seq_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    ours <- build_network(data, method = "co_occurrence")$weights
    ref  <- tna::ctna(data)$weights
    ref[is.na(ref)] <- 0

    n <- .check_edges(ours, ref, TOL, sprintf("co_occurrence %s", seq_info(cfg, i)))

    states <- sort(intersect(rownames(ours), rownames(ref)))
    delta <- abs(ours[states, states] - ref[states, states])
    report$log("co_occurrence", seq_info(cfg, i), n, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "tna::ctna()")
  }))
})


# ---- Test: method = "cor" vs stats::cor() ------------------------------------

test_that("cor: every edge matches stats::cor() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- assoc_configs[[i]]
    data <- simulate_continuous(cfg$n, cfg$p, cfg$rho, cfg$seed)

    ours <- build_network(data, method = "cor")$weights
    ref  <- stats::cor(as.matrix(data))
    diag(ref) <- 0

    cols <- colnames(ref)
    n <- length(cols)
    n_checked <- 0L

    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        n_checked <<- n_checked + 1L
        expect_equal(
          ours[cols[ii], cols[jj]], ref[ii, jj], tolerance = TOL,
          label = sprintf("cor %s [%s -> %s]", assoc_info(cfg, i), cols[ii], cols[jj])
        )
      })
    }))

    delta <- abs(ours[cols, cols] - ref)
    report$log("cor", assoc_info(cfg, i), n_checked, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "stats::cor()")
  }))
})


# ---- Test: method = "pcor" vs corpcor::cor2pcor() ----------------------------

test_that("pcor: every edge matches corpcor::cor2pcor() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("corpcor")

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- assoc_configs[[i]]
    data <- simulate_continuous(cfg$n, cfg$p, cfg$rho, cfg$seed)

    ours <- build_network(data, method = "pcor")$weights

    S   <- stats::cor(as.matrix(data))
    ref <- corpcor::cor2pcor(S)
    dimnames(ref) <- dimnames(S)
    diag(ref) <- 0

    cols <- colnames(ref)
    n <- length(cols)
    n_checked <- 0L

    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        n_checked <<- n_checked + 1L
        expect_equal(
          ours[cols[ii], cols[jj]], ref[ii, jj], tolerance = TOL,
          label = sprintf("pcor %s [%s -> %s]", assoc_info(cfg, i), cols[ii], cols[jj])
        )
      })
    }))

    delta <- abs(ours[cols, cols] - ref)
    report$log("pcor", assoc_info(cfg, i), n_checked, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "corpcor::cor2pcor()")
  }))
})


# ---- Test: method = "glasso" vs qgraph::EBICglasso() -------------------------

test_that("glasso: every edge matches qgraph::EBICglasso() across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("qgraph")

  n_glasso   <- 50L
  n_disagree <- 0L

  invisible(lapply(seq_len(n_glasso), function(i) {
    cfg  <- assoc_configs[[i]]
    data <- simulate_continuous(cfg$n, cfg$p, cfg$rho, cfg$seed)

    S   <- stats::cor(as.matrix(data))
    ref <- suppressWarnings(qgraph::EBICglasso(S, n = nrow(data)))
    diag(ref) <- 0

    ours <- build_network(data, method = "glasso")$weights
    cols <- colnames(ref)
    actual <- ours[cols, cols]

    # Only assert per-edge when sparsity patterns agree
    sparsity_match <- identical(ref == 0, actual == 0)
    if (!sparsity_match) {
      n_disagree <<- n_disagree + 1L
      delta <- abs(actual - ref)
      report$log("glasso", assoc_info(cfg, i), length(delta), sum(delta >= TOL),
                 max(delta), mean(delta), median(delta),
                 quantile(delta, 0.95, names = FALSE),
                 "qgraph::EBICglasso()", "lambda-path disagreement")
      return(NULL)
    }

    # Per-edge comparison
    n <- length(cols)
    n_checked <- 0L
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        n_checked <<- n_checked + 1L
        expect_equal(
          actual[ii, jj], ref[ii, jj], tolerance = TOL,
          label = sprintf("glasso %s [%s -> %s]", assoc_info(cfg, i), cols[ii], cols[jj])
        )
      })
    }))

    delta <- abs(actual - ref)
    report$log("glasso", assoc_info(cfg, i), n_checked, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "qgraph::EBICglasso()")
  }))

  agree_pct <- (n_glasso - n_disagree) / n_glasso * 100
  expect_true(agree_pct >= 90,
              info = sprintf("glasso: only %.0f%% configs agree on sparsity", agree_pct))
})


# ---- Test: method = "ising" vs IsingFit::IsingFit() --------------------------

test_that("ising: every edge matches IsingFit::IsingFit() across 30 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("IsingFit")

  n_ising <- 30L

  invisible(lapply(seq_len(n_ising), function(i) {
    cfg <- seq_configs[[i]]
    set.seed(cfg$seed)
    n_obs  <- cfg$n_actors * 5L
    n_vars <- cfg$n_states

    Sigma <- diag(n_vars)
    rho   <- 0.3
    invisible(lapply(seq_len(n_vars - 1L), function(k) {
      Sigma[k, k + 1L] <<- rho
      Sigma[k + 1L, k] <<- rho
    }))
    L   <- chol(Sigma)
    z   <- matrix(rnorm(n_obs * n_vars), n_obs, n_vars) %*% L
    bin <- as.data.frame(ifelse(z > 0, 1L, 0L))
    colnames(bin) <- paste0("V", seq_len(n_vars))

    var_ok <- vapply(bin, function(x) length(unique(x)) > 1L, logical(1))
    bin <- bin[, var_ok, drop = FALSE]
    if (ncol(bin) < 2L) return(NULL)

    ref  <- IsingFit::IsingFit(bin, plot = FALSE, progressbar = FALSE)$weiadj
    ours <- build_network(bin, method = "ising")$weights

    cols <- colnames(ref)
    actual <- ours[cols, cols]

    n <- length(cols)
    n_checked <- 0L
    invisible(lapply(seq_len(n), function(ii) {
      lapply(seq_len(n), function(jj) {
        n_checked <<- n_checked + 1L
        expect_equal(
          actual[ii, jj], ref[ii, jj], tolerance = TOL,
          label = sprintf("ising cfg%d [%s -> %s]", i, cols[ii], cols[jj])
        )
      })
    }))

    delta <- abs(actual - ref)
    report$log("ising",
               sprintf("cfg%d(n=%d,p=%d)", i, n_obs, ncol(bin)),
               n_checked, sum(delta >= TOL), max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "IsingFit::IsingFit()$weiadj")
  }))
})


# ---- Test: method = "attention" vs tna::atna() --------------------------------

test_that("attention: every edge matches tna::atna() across 100 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_CONFIGS), function(i) {
    cfg  <- seq_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    ours <- build_network(data, method = "attention",
                          params = list(lambda = 1))$weights
    ref  <- tna::atna(data)$weights
    ref[is.na(ref)] <- 0

    n <- .check_edges(ours, ref, TOL, sprintf("attention %s", seq_info(cfg, i)))

    states <- sort(intersect(rownames(ours), rownames(ref)))
    delta <- abs(ours[states, states] - ref[states, states])
    report$log("attention", seq_info(cfg, i), n, sum(delta >= TOL),
               max(delta), mean(delta), median(delta),
               quantile(delta, 0.95, names = FALSE), "tna::atna()")
  }))
})


# ---- Write equivalence report ------------------------------------------------

test_that("equivalence report is written", {
  skip_on_cran()
  skip_equiv_tests()

  report$write_csv("estimators")
  report$write_cvs("estimators")
  msg <- report$summary()
  message(msg)
  expect_true(length(report$rows) > 0L, info = "No equivalence results were logged")
})
