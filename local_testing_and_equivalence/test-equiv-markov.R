# ---- Equivalence: passage_time() vs markovchain::meanFirstPassageTime() ----
#
# markovchain sets diagonal to 0 (strict first-passage definition).
# passage_time() sets diagonal to 1/pi (mean recurrence time), matching
# the Kemeny-Snell convention and the paper's Figure.
# Equivalence is therefore checked on off-diagonal elements only.

skip_equiv_tests <- function() {
  run <- Sys.getenv("NESTIMATE_EQUIV_TESTS", unset = "false")
  if (!identical(run, "true"))
    skip("Equivalence tests skipped (set NESTIMATE_EQUIV_TESTS=true)")
}

TOL <- 1e-8

# ---- helpers ----

random_stochastic <- function(n, seed) {
  set.seed(seed)
  P <- matrix(runif(n * n), n, n)
  # Ensure irreducibility: add small floor so no row is all-zero to one state
  P <- P + 0.05
  P / rowSums(P)
}

off_diag <- function(M) {
  M[row(M) != col(M)]
}

# ---- 100-config equivalence sweep ----

test_that("passage_time off-diagonal matches markovchain::meanFirstPassageTime (100 configs)", {
  skip_equiv_tests()
  skip_if_not_installed("markovchain")

  rep <- equiv_report()

  sizes <- c(3L, 4L, 5L, 6L, 8L)
  seeds <- seq_len(20L)

  for (n in sizes) {
    for (seed in seeds) {
      P      <- random_stochastic(n, seed = seed * 1000L + n)
      states <- paste0("S", seq_len(n))
      colnames(P) <- rownames(P) <- states

      # Nestimate
      pt <- passage_time(P)

      # markovchain reference (off-diagonal only)
      mc  <- methods::new("markovchain", transitionMatrix = P, states = states)
      ref <- markovchain::meanFirstPassageTime(mc)

      our_off <- off_diag(pt$matrix)
      ref_off <- off_diag(ref)

      diffs   <- abs(our_off - ref_off)
      n_fail  <- sum(diffs > TOL)

      rep$log(
        func           = "passage_time",
        config         = sprintf("n=%d seed=%d", n, seed),
        n_checked      = length(diffs),
        n_failed       = n_fail,
        max_abs_err    = max(diffs),
        mean_abs_err   = mean(diffs),
        median_abs_err = stats::median(diffs),
        p95_abs_err    = stats::quantile(diffs, 0.95),
        reference      = "markovchain::meanFirstPassageTime",
        notes          = "off-diagonal only; diagonal convention differs"
      )

      expect_true(
        n_fail == 0L,
        info = sprintf("n=%d seed=%d: max delta=%.2e", n, seed, max(diffs))
      )
    }
  }

  rep$write_csv("markov")
  message(rep$summary())
})

# ---- diagonal convention is recurrence time ----

test_that("passage_time diagonal equals 1/pi (recurrence time), not 0", {
  skip_equiv_tests()

  P      <- random_stochastic(4L, seed = 42L)
  states <- paste0("S", seq_len(4L))
  colnames(P) <- rownames(P) <- states

  pt <- passage_time(P)
  expect_equal(diag(pt$matrix), pt$return_times,     tolerance = TOL)
  expect_equal(pt$return_times,  1 / pt$stationary,  tolerance = TOL)
  expect_true(all(diag(pt$matrix) > 0))
})

# ---- stationary distribution matches markovchain::steadyStates ----

test_that("stationary distribution matches markovchain::steadyStates", {
  skip_equiv_tests()
  skip_if_not_installed("markovchain")

  for (seed in c(1L, 7L, 13L, 42L, 99L)) {
    P      <- random_stochastic(5L, seed = seed)
    states <- paste0("S", seq_len(5L))
    colnames(P) <- rownames(P) <- states

    pt  <- passage_time(P)
    mc  <- methods::new("markovchain", transitionMatrix = P, states = states)
    ref <- as.vector(markovchain::steadyStates(mc))

    expect_equal(sort(unname(pt$stationary)), sort(unname(ref)), tolerance = TOL,
                 info = sprintf("seed=%d", seed))
  }
})
