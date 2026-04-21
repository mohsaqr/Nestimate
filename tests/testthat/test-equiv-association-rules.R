# ---- Equivalence: association_rules() vs arules::apriori() ----
#
# Validates that Nestimate::association_rules() produces a rule set numerically
# identical to arules::apriori(..., target = "rules") for matching thresholds,
# across a 108-config random sweep (n_trans in {50,100,200}, n_items in
# {5,8,12}, support in {0.05,0.1,0.2}, confidence in {0.3,0.5,0.7}, 3 seeds
# each). Comparison is restricted to Nestimate's single-consequent rules
# because arules::apriori() only emits single-consequent rules by default.
# Rule set identity is checked exact-set; for common rules, support /
# confidence / lift deltas must be <= TOL (1e-8).

skip_equiv_tests <- function() {
  run <- Sys.getenv("NESTIMATE_EQUIV_TESTS", unset = "false")
  if (!identical(run, "true"))
    skip("Equivalence tests skipped (set NESTIMATE_EQUIV_TESTS=true)")
}

TOL <- 1e-8

# ---- helpers ----

#' Generate random transactions as a list of character vectors.
#' density controls the Bernoulli probability of each item appearing in a
#' transaction (independent across items); transactions with 0 items get a
#' random singleton so arules doesn't error on empty rows.
simulate_transactions <- function(n_trans, n_items, density = 0.25,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  items <- paste0("I", sprintf("%02d", seq_len(n_items)))
  mat <- matrix(stats::runif(n_trans * n_items) < density,
                nrow = n_trans, ncol = n_items)
  trans <- lapply(seq_len(n_trans), function(i) sort(items[mat[i, ]]))
  # Fill empties so arules doesn't balk
  empty <- lengths(trans) == 0L
  if (any(empty)) {
    trans[empty] <- lapply(which(empty),
                           function(i) sample(items, 1L))
  }
  trans
}

#' Sort an item vector and build a tab-separated key for set-equality.
.key_sorted <- function(v) paste(sort(v), collapse = "\t")

#' Build canonical rule key "lhs_sorted => rhs_sorted" from a vector of
#' comma-separated item strings (as Nestimate returns them).
.ours_rule_key <- function(ante_str, cons_str) {
  a <- .key_sorted(strsplit(ante_str, ", ", fixed = TRUE)[[1]])
  b <- .key_sorted(strsplit(cons_str, ", ", fixed = TRUE)[[1]])
  paste(a, "=>", b)
}

#' Build canonical rule key from arules lhs/rhs lists of character vectors.
.arules_rule_key <- function(lhs_vec, rhs_vec) {
  paste(.key_sorted(lhs_vec), "=>", .key_sorted(rhs_vec))
}

# ---- 100+ config equivalence sweep ----

test_that("association_rules matches arules::apriori across 108 configs", {
  skip_equiv_tests()
  skip_if_not_installed("arules")

  rep <- equiv_report()

  n_trans_vals  <- c(50L, 100L, 200L)
  n_items_vals  <- c(5L, 8L, 12L)
  support_vals  <- c(0.05, 0.1, 0.2)
  conf_vals     <- c(0.3, 0.5, 0.7)
  seeds         <- c(11L, 47L, 113L)  # 3 seeds -> 3*3*3*3*3 = 243? no, we
  # iterate configs linearly; see grid below.

  grid <- expand.grid(
    n_trans = n_trans_vals,
    n_items = n_items_vals,
    support = support_vals,
    conf    = conf_vals,
    seed    = seeds,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  # Keep test runtime reasonable: take a random-but-deterministic sample of
  # the grid such that each (n_trans, n_items, support, conf) level appears
  # with at least one seed. We have 3*3*3*3*3 = 243 configs — cap to ~108 by
  # stratified thinning (drop one of the three seeds per cell randomly).
  set.seed(2026L)
  keep <- unlist(lapply(
    split(seq_len(nrow(grid)),
          interaction(grid$n_trans, grid$n_items,
                      grid$support, grid$conf, drop = TRUE)),
    function(idx) sample(idx, size = min(2L, length(idx)))
  ), use.names = FALSE)
  grid <- grid[sort(keep), , drop = FALSE]

  n_configs           <- nrow(grid)
  n_identity_match    <- 0L
  n_numeric_fail      <- 0L
  n_configs_skipped   <- 0L

  for (row_i in seq_len(n_configs)) {
    cfg <- grid[row_i, ]
    density <- 0.1 + stats::runif(1, 0, 0.3)  # 0.1-0.4 jitter per config
    trans <- simulate_transactions(cfg$n_trans, cfg$n_items,
                                   density = density,
                                   seed = cfg$seed * 7919L + row_i)

    # Nestimate
    ours <- association_rules(trans,
                              min_support    = cfg$support,
                              min_confidence = cfg$conf,
                              min_lift       = 0,
                              max_length     = 5L)
    r <- ours$rules
    # Restrict to single-consequent (arules default)
    if (nrow(r) > 0L) {
      n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
      r <- r[n_cons == 1L, , drop = FALSE]
    }

    # arules reference
    t_obj <- methods::as(trans, "transactions")
    theirs <- suppressWarnings(
      arules::apriori(
        t_obj,
        parameter = list(supp   = cfg$support,
                         conf   = cfg$conf,
                         minlen = 2L,
                         maxlen = 5L,
                         target = "rules"),
        control   = list(verbose = FALSE)
      )
    )

    # Handle both-empty case
    if (nrow(r) == 0L && length(theirs) == 0L) {
      n_identity_match <- n_identity_match + 1L
      rep$log(
        func           = "association_rules",
        config         = sprintf("n=%d k=%d sup=%.2f conf=%.2f seed=%d",
                                 cfg$n_trans, cfg$n_items, cfg$support,
                                 cfg$conf, cfg$seed),
        n_checked      = 0L,
        n_failed       = 0L,
        max_abs_err    = 0,
        mean_abs_err   = 0,
        median_abs_err = 0,
        p95_abs_err    = 0,
        reference      = "arules::apriori",
        notes          = "both-empty"
      )
      next
    }

    # Build rule keys
    our_keys <- if (nrow(r) > 0L) {
      vapply(seq_len(nrow(r)),
             function(i) .ours_rule_key(r$antecedent[i], r$consequent[i]),
             character(1))
    } else character(0)

    if (length(theirs) > 0L) {
      lhs_list <- methods::as(arules::lhs(theirs), "list")
      rhs_list <- methods::as(arules::rhs(theirs), "list")
      their_keys <- vapply(seq_along(lhs_list),
                           function(i) .arules_rule_key(lhs_list[[i]],
                                                        rhs_list[[i]]),
                           character(1))
      q <- arules::quality(theirs)
    } else {
      their_keys <- character(0)
      q <- NULL
    }

    # Set-identity check
    set_equal <- setequal(our_keys, their_keys)
    if (set_equal) n_identity_match <- n_identity_match + 1L

    # Numeric comparison on common keys
    common <- intersect(our_keys, their_keys)
    sup_d <- conf_d <- lift_d <- numeric(0)
    if (length(common) > 0L) {
      sup_d  <- numeric(length(common))
      conf_d <- numeric(length(common))
      lift_d <- numeric(length(common))
      for (i in seq_along(common)) {
        oi <- which(our_keys   == common[i])[1]
        ai <- which(their_keys == common[i])[1]
        sup_d[i]  <- abs(r$support[oi]    - q$support[ai])
        conf_d[i] <- abs(r$confidence[oi] - q$confidence[ai])
        lift_d[i] <- abs(r$lift[oi]       - q$lift[ai])
      }
    }

    all_d <- c(sup_d, conf_d, lift_d)
    missing_from_ours <- length(setdiff(their_keys, our_keys))
    extra_in_ours     <- length(setdiff(our_keys, their_keys))
    n_fail_numeric    <- sum(all_d > TOL)
    n_fail_total      <- n_fail_numeric + missing_from_ours + extra_in_ours
    if (n_fail_total > 0L) n_numeric_fail <- n_numeric_fail + 1L

    rep$log(
      func           = "association_rules",
      config         = sprintf("n=%d k=%d sup=%.2f conf=%.2f seed=%d",
                               cfg$n_trans, cfg$n_items, cfg$support,
                               cfg$conf, cfg$seed),
      n_checked      = length(all_d) + length(our_keys) + length(their_keys),
      n_failed       = n_fail_total,
      max_abs_err    = if (length(all_d)) max(all_d) else 0,
      mean_abs_err   = if (length(all_d)) mean(all_d) else 0,
      median_abs_err = if (length(all_d)) stats::median(all_d) else 0,
      p95_abs_err    = if (length(all_d)) as.numeric(stats::quantile(all_d, 0.95)) else 0,
      reference      = "arules::apriori",
      notes          = sprintf("set_equal=%s missing=%d extra=%d n_our=%d n_their=%d",
                               set_equal, missing_from_ours, extra_in_ours,
                               length(our_keys), length(their_keys))
    )

    expect_true(
      n_fail_total == 0L,
      info = sprintf(
        "n=%d k=%d sup=%.2f conf=%.2f seed=%d: max delta=%.2e missing=%d extra=%d",
        cfg$n_trans, cfg$n_items, cfg$support, cfg$conf, cfg$seed,
        if (length(all_d)) max(all_d) else 0,
        missing_from_ours, extra_in_ours
      )
    )
  }

  rep$write_csv("association_rules")
  message(sprintf("Rule-set identity hit rate: %d/%d",
                  n_identity_match, n_configs))
  message(rep$summary())

  expect_true(n_numeric_fail == 0L,
              info = sprintf("%d/%d configs failed equivalence",
                             n_numeric_fail, n_configs))
  expect_true(n_identity_match == n_configs,
              info = sprintf("%d/%d configs had identical rule sets",
                             n_identity_match, n_configs))
})
