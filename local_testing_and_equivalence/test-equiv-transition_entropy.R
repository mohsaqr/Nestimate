# Numerical equivalence: transition_entropy() vs reference composition
#
# Validates Nestimate's transition_entropy() against an independent
# composition of two well-known packages:
#   - per-row Shannon entropy : entropy::entropy.plugin(p, unit = "log2")
#   - stationary distribution : markovchain::steadyStates(mc)
#   - entropy rate            : sum(pi * row_entropy)
#
# The reference therefore uses a different stationary-distribution solver
# (markovchain's iterative method via eigen on t(P) is similar in spirit but
# implemented separately) and a fundamentally independent per-row entropy
# implementation, so equivalence at machine precision exercises the whole
# pipeline.
#
# Sweeps 100 random irreducible aperiodic configurations spanning n in
# {3, ..., 7}. Bases tested: 2 (bits), e (nats), 10 (hartleys).
# TOL: 1e-10 for h(P), H(pi), and per-row H values.

set.seed(20260427)
N_CONFIGS <- 100L
TOL       <- 1e-10


# ---- Random irreducible aperiodic generator -----------------------------
# Lower bound 0.05 on each entry guarantees full support, hence
# irreducible + aperiodic; eigendecomposition for pi has a clean unique
# leading eigenvector.

.gen_regular <- function(n_states, seed) {
  set.seed(seed)
  M <- matrix(stats::runif(n_states^2, 0.05, 1), n_states, n_states)
  M / rowSums(M)
}


# ---- Reference: row entropies via entropy package -----------------------

.ref_row_entropy <- function(P, base) {
  unit <- if (isTRUE(all.equal(base, 2)))      "log2"
          else if (isTRUE(all.equal(base, exp(1)))) "log"
          else if (isTRUE(all.equal(base, 10))) "log10"
          else stop("Unsupported base for entropy::entropy.plugin")
  vapply(seq_len(nrow(P)),
         function(i) entropy::entropy.plugin(P[i, ], unit = unit),
         numeric(1))
}

.ref_stationary <- function(P) {
  state_names <- rownames(P)
  if (is.null(state_names)) state_names <- paste0("s", seq_len(nrow(P)))
  rownames(P) <- colnames(P) <- state_names
  mc <- methods::new("markovchain", states = state_names, byrow = TRUE,
                     transitionMatrix = P)
  ss <- markovchain::steadyStates(mc)
  drop(ss[1L, ])
}

.ref_transition_entropy <- function(P, base = 2) {
  H_row <- .ref_row_entropy(P, base = base)
  pi    <- .ref_stationary(P)
  pi    <- pi[colnames(P)]                  # align ordering
  list(
    row_entropy        = setNames(H_row, colnames(P)),
    stationary         = pi,
    entropy_rate       = sum(pi * H_row),
    stationary_entropy = {
      p <- pi[pi > 0]
      -sum(p * log(p, base = base))
    }
  )
}


# ---- Configurations -----------------------------------------------------

configs <- lapply(seq_len(N_CONFIGS), function(i) {
  list(n_states = sample(3:7, 1L),
       seed     = sample.int(1e6, 1L),
       base     = sample(c(2, exp(1), 10), 1L))
})

.gen_P <- function(cfg) {
  P <- .gen_regular(cfg$n_states, cfg$seed)
  rownames(P) <- colnames(P) <- paste0("s", seq_len(nrow(P)))
  P
}


# ---- Tests --------------------------------------------------------------

test_that("transition_entropy: entropy rate matches reference (1e-10)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  skip_if_not_installed("entropy")

  diffs <- vapply(configs, function(cfg) {
    P    <- .gen_P(cfg)
    ours <- Nestimate::transition_entropy(P, base = cfg$base)
    ref  <- .ref_transition_entropy(P, base = cfg$base)
    abs(ours$entropy_rate - ref$entropy_rate)
  }, numeric(1))

  expect_true(all(diffs < TOL),
              info = sprintf("max abs error = %.3e (tol %.0e)",
                             max(diffs), TOL))
  cat(sprintf("\n  entropy_rate       max abs err: %.3e (n=%d)\n",
              max(diffs), N_CONFIGS))
})

test_that("transition_entropy: stationary entropy matches reference (1e-10)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  skip_if_not_installed("entropy")

  diffs <- vapply(configs, function(cfg) {
    P    <- .gen_P(cfg)
    ours <- Nestimate::transition_entropy(P, base = cfg$base)
    ref  <- .ref_transition_entropy(P, base = cfg$base)
    abs(ours$stationary_entropy - ref$stationary_entropy)
  }, numeric(1))

  expect_true(all(diffs < TOL),
              info = sprintf("max abs error = %.3e (tol %.0e)",
                             max(diffs), TOL))
  cat(sprintf("  stationary_entropy max abs err: %.3e (n=%d)\n",
              max(diffs), N_CONFIGS))
})

test_that("transition_entropy: per-row entropies match reference (1e-10)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  skip_if_not_installed("entropy")

  diffs <- vapply(configs, function(cfg) {
    P    <- .gen_P(cfg)
    ours <- Nestimate::transition_entropy(P, base = cfg$base)
    ref  <- .ref_transition_entropy(P, base = cfg$base)
    max(abs(unname(ours$row_entropy) - unname(ref$row_entropy)))
  }, numeric(1))

  expect_true(all(diffs < TOL),
              info = sprintf("max abs error = %.3e (tol %.0e)",
                             max(diffs), TOL))
  cat(sprintf("  row_entropy        max abs err: %.3e (n=%d)\n",
              max(diffs), N_CONFIGS))
})

test_that("transition_entropy: stationary distribution matches reference (1e-8)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")

  diffs <- vapply(configs, function(cfg) {
    P    <- .gen_P(cfg)
    ours <- Nestimate::transition_entropy(P, base = cfg$base)
    ref  <- .ref_transition_entropy(P, base = cfg$base)
    max(abs(unname(ours$stationary) - unname(ref$stationary)))
  }, numeric(1))

  expect_true(all(diffs < 1e-8),
              info = sprintf("max abs error = %.3e", max(diffs)))
  cat(sprintf("  stationary         max abs err: %.3e (n=%d)\n",
              max(diffs), N_CONFIGS))
})


# ---- Validation report -------------------------------------------------

if (identical(Sys.getenv("NESTIMATE_EQUIV_TESTS"), "true")) {
  report <- data.frame(
    config = seq_len(N_CONFIGS),
    n      = vapply(configs, `[[`, integer(1), "n_states"),
    base   = vapply(configs, `[[`, numeric(1), "base"),
    seed   = vapply(configs, `[[`, integer(1), "seed"),
    stringsAsFactors = FALSE
  )

  per_config <- t(vapply(configs, function(cfg) {
    P    <- .gen_P(cfg)
    ours <- Nestimate::transition_entropy(P, base = cfg$base)
    ref  <- .ref_transition_entropy(P, base = cfg$base)
    c(
      h_err   = abs(ours$entropy_rate       - ref$entropy_rate),
      Hpi_err = abs(ours$stationary_entropy - ref$stationary_entropy),
      Hrow_err = max(abs(unname(ours$row_entropy) - unname(ref$row_entropy))),
      pi_err  = max(abs(unname(ours$stationary)  - unname(ref$stationary)))
    )
  }, numeric(4)))

  report <- cbind(report, as.data.frame(per_config))

  if (!dir.exists("tmp")) dir.create("tmp", recursive = TRUE)
  utils::write.csv(report,
                   "tmp/transition_entropy_equivalence_report.csv",
                   row.names = FALSE)
  cat(sprintf(
    "\n  Equivalence report written: tmp/transition_entropy_equivalence_report.csv (%d rows)\n",
    nrow(report)))
  cat(sprintf("  h_err    median %.2e | p95 %.2e | max %.2e\n",
              median(report$h_err),  stats::quantile(report$h_err,  0.95),
              max(report$h_err)))
  cat(sprintf("  Hpi_err  median %.2e | p95 %.2e | max %.2e\n",
              median(report$Hpi_err), stats::quantile(report$Hpi_err, 0.95),
              max(report$Hpi_err)))
  cat(sprintf("  Hrow_err median %.2e | p95 %.2e | max %.2e\n",
              median(report$Hrow_err), stats::quantile(report$Hrow_err, 0.95),
              max(report$Hrow_err)))
  cat(sprintf("  pi_err   median %.2e | p95 %.2e | max %.2e\n",
              median(report$pi_err),  stats::quantile(report$pi_err,  0.95),
              max(report$pi_err)))
}
