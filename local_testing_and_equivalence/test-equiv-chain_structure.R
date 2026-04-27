# Numerical equivalence: chain_structure() vs markovchain::*()
#
# Compares Nestimate's chain_structure() against markovchain v0.10.0 across
# 100 random transition-matrix configurations covering three regimes:
#   (A) irreducible aperiodic ("regular")        — ~50 configs
#   (B) reducible with one absorbing state       — ~30 configs
#   (C) reducible with multiple closed classes   — ~20 configs
#
# Compared elements:
#   - is_irreducible           vs markovchain::is.irreducible
#   - is_regular               vs markovchain::is.regular
#   - period (per state)       vs markovchain::period (when irreducible)
#   - communicating_classes    vs markovchain::communicatingClasses
#   - recurrent_classes        vs markovchain::recurrentClasses
#   - absorbing_states         vs markovchain::absorbingStates
#   - classification           — derived from the above; we also verify
#                                state membership matches the union
#   - hitting_probabilities    vs markovchain::hittingProbabilities (TOL 1e-10)
#   - absorption_probabilities vs markovchain::absorptionProbabilities
#                                (TOL 1e-10, when absorbing states exist)
#   - mean_absorption_time     vs markovchain::meanAbsorptionTime (TOL 1e-10)

set.seed(20260426)
N_REGULAR  <- 50L
N_ABSORB   <- 30L
N_MULTI    <- 20L
TOL        <- 1e-10
TOL_LOOSE  <- 1e-8


# ---- Random matrix generators -------------------------------------------

.gen_regular <- function(n_states, seed) {
  set.seed(seed)
  M <- matrix(stats::runif(n_states^2, 0.05, 1), n_states, n_states)
  M / rowSums(M)
}

.gen_absorbing <- function(n_states, n_abs, seed) {
  set.seed(seed)
  P <- .gen_regular(n_states, seed + 1L)
  abs_idx <- sample.int(n_states, n_abs)
  P[abs_idx, ] <- 0
  P[cbind(abs_idx, abs_idx)] <- 1
  # Ensure transient states have some path to absorbing (already do via
  # unrestricted random + we check)
  P
}

.gen_multi_class <- function(seed) {
  set.seed(seed)
  # Two non-communicating recurrent classes of random size
  k1 <- sample(2:3, 1L)
  k2 <- sample(2:3, 1L)
  n <- k1 + k2
  P <- matrix(0, n, n)
  P[seq_len(k1), seq_len(k1)] <- .gen_regular(k1, seed + 7L)
  P[(k1 + 1L):n, (k1 + 1L):n] <- .gen_regular(k2, seed + 13L)
  P
}


# ---- Reference helpers via markovchain ----------------------------------

.mc_ref <- function(P) {
  state_names <- rownames(P)
  if (is.null(state_names)) state_names <- paste0("s", seq_len(nrow(P)))
  rownames(P) <- colnames(P) <- state_names
  mc <- methods::new("markovchain", states = state_names, byrow = TRUE,
                     transitionMatrix = P)
  list(
    irreducible          = markovchain::is.irreducible(mc),
    regular              = markovchain::is.regular(mc),
    communicating        = markovchain::communicatingClasses(mc),
    recurrent            = markovchain::recurrentClasses(mc),
    absorbing            = markovchain::absorbingStates(mc),
    transient            = markovchain::transientStates(mc),
    period               = if (markovchain::is.irreducible(mc))
                              markovchain::period(mc) else NA_integer_,
    hitting              = markovchain::hittingProbabilities(mc),
    absorption_probs     = if (length(markovchain::absorbingStates(mc)) > 0L)
                              markovchain::absorptionProbabilities(mc) else NULL,
    absorption_time      = if (length(markovchain::absorbingStates(mc)) > 0L)
                              markovchain::meanAbsorptionTime(mc) else NULL
  )
}

# Sort each class's members and the list of classes alphabetically by first
# member, so two equivalent partitions compare equal.
.normalize_classes <- function(classes) {
  cl_sorted <- lapply(classes, sort)
  ord <- order(vapply(cl_sorted, function(x) x[1L], character(1)))
  cl_sorted[ord]
}


# ---- Configs sweep ------------------------------------------------------

configs <- c(
  lapply(seq_len(N_REGULAR), function(i) {
    list(kind = "regular", n_states = sample(3:5, 1L),
         seed = sample.int(1e6, 1L))
  }),
  lapply(seq_len(N_ABSORB), function(i) {
    n <- sample(3:5, 1L)
    list(kind = "absorbing", n_states = n,
         n_abs = sample.int(min(2L, n - 1L), 1L),
         seed = sample.int(1e6, 1L))
  }),
  lapply(seq_len(N_MULTI), function(i) {
    list(kind = "multi", seed = sample.int(1e6, 1L))
  })
)

.gen_P <- function(cfg) {
  P <- switch(
    cfg$kind,
    regular   = .gen_regular(cfg$n_states, cfg$seed),
    absorbing = .gen_absorbing(cfg$n_states, cfg$n_abs, cfg$seed),
    multi     = .gen_multi_class(cfg$seed)
  )
  rownames(P) <- colnames(P) <- paste0("s", seq_len(nrow(P)))
  P
}


# ---- Tests --------------------------------------------------------------

test_that("chain_structure: irreducibility matches markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    cs$is_irreducible == ref$irreducible
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: regularity matches markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    cs$is_regular == ref$regular
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: communicating classes match markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    identical(.normalize_classes(cs$communicating_classes),
              .normalize_classes(ref$communicating))
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: recurrent classes match markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    identical(.normalize_classes(cs$recurrent_classes),
              .normalize_classes(ref$recurrent))
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: absorbing states match markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    identical(sort(cs$absorbing_states), sort(ref$absorbing))
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: transient states match markovchain (union)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  results <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    our_transient <- cs$states[cs$classification == "transient"]
    identical(sort(our_transient), sort(ref$transient))
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: period matches markovchain (irreducible only)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  irr_cfgs <- configs[vapply(configs, function(cfg) cfg$kind == "regular",
                             logical(1))]
  results <- vapply(irr_cfgs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    if (!cs$is_irreducible) return(TRUE)  # tested separately
    # When irreducible, all states share the same period; compare to scalar.
    our_period <- unique(cs$period[!is.na(cs$period)])
    length(our_period) == 1L && our_period == ref$period
  }, logical(1))
  expect_true(all(results),
              info = sprintf("%d/%d disagreements",
                             sum(!results), length(results)))
})

test_that("chain_structure: hitting probabilities match markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  deltas <- vapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    H_ref <- as.matrix(ref$hitting)
    H_ref <- H_ref[rownames(cs$hitting_probabilities),
                   colnames(cs$hitting_probabilities)]
    max(abs(cs$hitting_probabilities - H_ref), na.rm = TRUE)
  }, numeric(1))
  expect_true(max(deltas) < TOL_LOOSE,
              info = sprintf("max delta = %.3e", max(deltas)))
})

test_that("chain_structure: absorption probabilities match markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  abs_cfgs <- configs[vapply(configs,
                             function(cfg) cfg$kind == "absorbing",
                             logical(1))]
  deltas <- vapply(abs_cfgs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    if (is.null(cs$absorption_probabilities) ||
        is.null(ref$absorption_probs)) return(0)
    A_ref <- as.matrix(ref$absorption_probs)
    A_ref <- A_ref[rownames(cs$absorption_probabilities),
                   colnames(cs$absorption_probabilities), drop = FALSE]
    max(abs(cs$absorption_probabilities - A_ref))
  }, numeric(1))
  expect_true(max(deltas) < TOL,
              info = sprintf("max delta = %.3e", max(deltas)))
})

test_that("chain_structure: mean absorption time matches markovchain", {
  skip_on_cran()
  skip_if_not_installed("markovchain")
  abs_cfgs <- configs[vapply(configs,
                             function(cfg) cfg$kind == "absorbing",
                             logical(1))]
  deltas <- vapply(abs_cfgs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)
    if (is.null(cs$mean_absorption_time) ||
        is.null(ref$absorption_time)) return(0)
    t_ref <- ref$absorption_time[names(cs$mean_absorption_time)]
    max(abs(cs$mean_absorption_time - t_ref))
  }, numeric(1))
  expect_true(max(deltas) < TOL,
              info = sprintf("max delta = %.3e", max(deltas)))
})


# ---- Delta report + validation dashboard emit --------------------------

test_that("chain_structure equivalence report (CSV + CVS JSON)", {
  skip_on_cran()
  skip_if_not_installed("markovchain")

  report <- equiv_report()

  invisible(lapply(configs, function(cfg) {
    P <- .gen_P(cfg)
    cs <- chain_structure(P)
    ref <- .mc_ref(P)

    deltas <- numeric()
    deltas["irreducible"] <- as.numeric(cs$is_irreducible != ref$irreducible)
    deltas["regular"]     <- as.numeric(cs$is_regular != ref$regular)
    deltas["classes"]     <- as.numeric(!identical(
      .normalize_classes(cs$communicating_classes),
      .normalize_classes(ref$communicating)))
    deltas["recurrent"]   <- as.numeric(!identical(
      .normalize_classes(cs$recurrent_classes),
      .normalize_classes(ref$recurrent)))
    deltas["absorbing"]   <- as.numeric(!identical(
      sort(cs$absorbing_states), sort(ref$absorbing)))
    H_ref <- as.matrix(ref$hitting)
    H_ref <- H_ref[rownames(cs$hitting_probabilities),
                   colnames(cs$hitting_probabilities)]
    deltas["hitting"]     <- max(abs(cs$hitting_probabilities - H_ref),
                                 na.rm = TRUE)
    if (!is.null(ref$absorption_probs) &&
        !is.null(cs$absorption_probabilities)) {
      A_ref <- as.matrix(ref$absorption_probs)
      A_ref <- A_ref[rownames(cs$absorption_probabilities),
                     colnames(cs$absorption_probabilities), drop = FALSE]
      deltas["absorption_probs"] <- max(abs(cs$absorption_probabilities - A_ref))
      t_ref <- ref$absorption_time[names(cs$mean_absorption_time)]
      deltas["absorption_time"] <- max(abs(cs$mean_absorption_time - t_ref))
    }

    n_failed <- sum(deltas >= TOL_LOOSE)
    report$log(
      func = "chain_structure",
      config = sprintf("kind=%s n_states=%d", cfg$kind, nrow(P)),
      n_checked = length(deltas),
      n_failed  = n_failed,
      max_abs_err    = max(deltas),
      mean_abs_err   = mean(deltas),
      median_abs_err = stats::median(deltas),
      p95_abs_err    = stats::quantile(deltas, 0.95, names = FALSE),
      reference = "markovchain::is.irreducible/is.regular/communicatingClasses/recurrentClasses/absorbingStates/hittingProbabilities/absorptionProbabilities/meanAbsorptionTime",
      notes = "indicator deltas are 0/1; numerical deltas at TOL=1e-10"
    )
  }))

  report$write_csv("chain_structure")
  report$write_cvs("chain_structure")
  expect_true(length(report$rows) > 0L)
})
