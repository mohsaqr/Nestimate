# ---- Equivalence Tests: Higher-Order Networks ----
#
# Verifies structural and numerical correctness of four higher-order
# network functions against manual ground-truth implementations:
#   1. build_hon()   — edge counts match manual k-gram counting
#   2. build_hypa()  — p-values match manual hypergeometric computation
#   3. build_honem() — embedding properties (SVD invariants)
#   4. build_mogen() — log-likelihood monotonicity, DOF, IC validity
#
# Requires: NESTIMATE_EQUIV_TESTS=true

# ---- Ground-truth helpers ----

#' Manual k-gram counting for HON verification
#'
#' For order-1 HON (max_order = 1), each node is a single state.
#' For each (from, to) pair, count how many times state `to` immediately
#' follows state `from` across all sequences.
#'
#' @param seqs List of character vectors.
#' @return data.frame with columns: from, to, count
#' @noRd
.manual_order1_counts <- function(seqs) {
  pairs <- do.call(rbind, lapply(seqs, function(s) {
    n <- length(s)
    if (n < 2L) return(data.frame(from = character(0), to = character(0),
                                  stringsAsFactors = FALSE))
    data.frame(from = s[-n], to = s[-1L], stringsAsFactors = FALSE)
  }))
  if (nrow(pairs) == 0L) {
    return(data.frame(from = character(0), to = character(0),
                      count = integer(0), stringsAsFactors = FALSE))
  }
  agg <- stats::aggregate(
    list(count = rep(1L, nrow(pairs))),
    by = list(from = pairs$from, to = pairs$to),
    FUN = sum
  )
  agg
}

#' Manual hypergeometric p-value computation for HYPA verification
#'
#' Replicates the De Bruijn graph construction and hypergeometric null model.
#'
#' @param seqs List of character vectors.
#' @param k Integer order.
#' @return data.frame with from_node, to_node, observed, p_value columns.
#' @noRd
.manual_hypa_pvalues <- function(seqs, k) {
  sep <- "\x01"

  # Build k-gram nodes and transitions
  all_from <- character(0)
  all_to <- character(0)

  from_to <- lapply(seqs, function(s) {
    n <- length(s)
    if (n < k) return(list(from = character(0), to = character(0)))
    kgrams <- vapply(seq_len(n - k + 1L), function(i) {
      paste(s[i:(i + k - 1L)], collapse = sep)
    }, character(1L))
    if (length(kgrams) > 1L) {
      list(from = kgrams[-length(kgrams)], to = kgrams[-1L])
    } else {
      list(from = character(0), to = character(0))
    }
  })

  all_from <- unlist(lapply(from_to, `[[`, "from"))
  all_to <- unlist(lapply(from_to, `[[`, "to"))

  if (length(all_from) == 0L) {
    return(data.frame(from_node = character(0), to_node = character(0),
                      observed = integer(0), p_value = numeric(0),
                      stringsAsFactors = FALSE))
  }

  # Build adjacency matrix
  edge_keys <- paste(all_from, all_to, sep = "\x02")
  edge_tab <- table(edge_keys)
  edge_split <- strsplit(names(edge_tab), "\x02", fixed = TRUE)
  edge_from <- vapply(edge_split, `[`, character(1), 1L)
  edge_to <- vapply(edge_split, `[`, character(1), 2L)
  edge_weight <- as.integer(edge_tab)

  nodes <- sort(unique(c(edge_from, edge_to)))
  n <- length(nodes)
  adj <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  idx <- cbind(match(edge_from, nodes), match(edge_to, nodes))
  adj[idx] <- edge_weight

  # Xi = outer(out_strength, in_strength) * mask
  out_str <- rowSums(adj)
  in_str <- colSums(adj)
  mask <- adj > 0
  xi <- outer(out_str, in_str) * mask

  # Hypergeometric parameters
  N_total <- round(sum(xi))
  n_draws <- sum(adj)

  # Compute p-values for each edge
  edge_idx <- which(adj > 0, arr.ind = TRUE)
  results <- lapply(seq_len(nrow(edge_idx)), function(ei) {
    i <- edge_idx[ei, 1L]
    j <- edge_idx[ei, 2L]
    f_obs <- adj[i, j]
    K <- min(max(round(xi[i, j]), 0L), N_total)
    n_clamp <- min(n_draws, N_total)
    p_val <- stats::phyper(f_obs, K, N_total - K, n_clamp)

    data.frame(
      from_node = nodes[i], to_node = nodes[j],
      observed = as.integer(f_obs), p_value = p_val,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}


# ---- Config generation ----
set.seed(4242)
N_HON <- 50L
hon_configs <- lapply(seq_len(N_HON), function(i) {
  list(n_actors = sample(c(10L, 15L, 20L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(15L, 20L, 30L), 1),
       seed = sample.int(100000, 1))
})


# ---- 1. HON edge counts match manual k-gram counting (50 configs) ----

test_that("HON edges have valid counts and probabilities across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()

  TOL <- 1e-10
  report <- equiv_report()

  results <- lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    hon <- build_hon(seqs, max_order = 2)
    edges <- hon$ho_edges

    # (a) All counts > 0
    all_pos <- all(edges$count > 0L, na.rm = TRUE)

    # (b) All probabilities in [0, 1]
    prob_valid <- all(edges$probability >= 0 - TOL &
                      edges$probability <= 1 + TOL)

    # (c) Probabilities from same source sum to <= 1 (+ tolerance)
    prob_sums <- tapply(edges$probability, edges$from, sum)
    prob_sums_ok <- all(prob_sums <= 1 + TOL)

    # (d) Cross-check with manual order-1 counts: for edges where from_order==1,
    #     the count should match the manual pair counting
    manual <- .manual_order1_counts(seqs)
    order1_edges <- edges[edges$from_order == 1L, , drop = FALSE]

    errs <- vapply(seq_len(nrow(order1_edges)), function(j) {
      e <- order1_edges[j, ]
      match_row <- manual$from == e$from & manual$to == e$to
      if (any(match_row)) {
        abs(e$count - manual$count[match_row][1L])
      } else {
        # Edge in HON but not in manual: this should not happen for order-1
        as.numeric(e$count)
      }
    }, numeric(1))

    n_failed <- sum(errs > TOL)
    max_err <- if (length(errs) > 0) max(errs) else 0

    report$log(
      func = "build_hon",
      config = sprintf("actors=%d,states=%d,len=%d,seed=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed),
      n_checked = nrow(order1_edges),
      n_failed = n_failed,
      max_abs_err = max_err,
      mean_abs_err = if (length(errs) > 0) mean(errs) else 0,
      median_abs_err = if (length(errs) > 0) stats::median(errs) else 0,
      p95_abs_err = if (length(errs) > 0) stats::quantile(errs, 0.95) else 0,
      reference = "manual_kgram",
      notes = sprintf("n_edges=%d,counts_pos=%s,probs_valid=%s,sums_ok=%s",
                       nrow(edges), all_pos, prob_valid, prob_sums_ok)
    )

    list(all_pos = all_pos, prob_valid = prob_valid,
         prob_sums_ok = prob_sums_ok, n_failed = n_failed,
         max_err = max_err, config_idx = i)
  })

  report$write_csv("hon")

  # Assertions
  lapply(results, function(r) {
    expect_true(r$all_pos,
                info = sprintf("Config %d: not all counts > 0", r$config_idx))
    expect_true(r$prob_valid,
                info = sprintf("Config %d: probabilities out of [0,1]",
                               r$config_idx))
    expect_true(r$prob_sums_ok,
                info = sprintf("Config %d: probability sums exceed 1",
                               r$config_idx))
    expect_true(r$n_failed == 0L,
                info = sprintf("Config %d: %d order-1 edges mismatch manual (max_err=%.2e)",
                               r$config_idx, r$n_failed, r$max_err))
  })
})


# ---- 2. HYPA p-values match manual hypergeometric (50 configs) ----

test_that("HYPA p-values match manual hypergeometric across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()

  TOL <- 1e-10
  report <- equiv_report()

  results <- lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    hypa <- build_hypa(seqs, k = 2, min_count = 1L, p_adjust = "none")
    manual <- .manual_hypa_pvalues(seqs, k = 2)

    scores <- hypa$scores

    # (a) All p-values in [0, 1]
    pval_valid <- all(scores$p_value >= 0 - TOL &
                      scores$p_value <= 1 + TOL)

    # (b) All observed counts are non-negative integers
    obs_valid <- all(scores$observed >= 0L) &&
                 all(scores$observed == as.integer(scores$observed))

    # (c) Classification is one of "over", "under", "normal"
    class_valid <- all(scores$anomaly %in% c("over", "under", "normal"))

    # (d) Match p-values against manual computation
    # Build lookup key for manual results
    sep <- "\x01"
    manual_key <- paste(manual$from_node, manual$to_node, sep = "||")
    manual_pval <- stats::setNames(manual$p_value, manual_key)
    manual_obs <- stats::setNames(manual$observed, manual_key)

    # Reconstruct the De Bruijn node names for the package scores
    # scores$from is "A -> B" format, need to convert to sep-joined for lookup
    score_from_parts <- strsplit(scores$from, " -> ", fixed = TRUE)
    score_from_keys <- vapply(score_from_parts, function(p) {
      paste(p, collapse = sep)
    }, character(1))

    # Reconstruct to-node key (from_parts[-1] + to)
    score_to_keys <- vapply(seq_len(nrow(scores)), function(j) {
      from_parts <- score_from_parts[[j]]
      # target node in De Bruijn = from_parts[-1] + scores$to[j]
      paste(c(from_parts[-1L], scores$to[j]), collapse = sep)
    }, character(1))

    score_keys <- paste(score_from_keys, score_to_keys, sep = "||")

    # Compare p-values where keys match
    matched <- score_keys %in% names(manual_pval)
    if (any(matched)) {
      pkg_pval <- scores$p_value[matched]
      ref_pval <- manual_pval[score_keys[matched]]
      errs <- abs(pkg_pval - ref_pval)
      n_failed <- sum(errs > TOL)
      max_err <- max(errs)
      mean_err <- mean(errs)
    } else {
      errs <- numeric(0)
      n_failed <- 0L
      max_err <- 0
      mean_err <- 0
    }

    report$log(
      func = "build_hypa",
      config = sprintf("actors=%d,states=%d,len=%d,seed=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed),
      n_checked = sum(matched),
      n_failed = n_failed,
      max_abs_err = max_err,
      mean_abs_err = mean_err,
      median_abs_err = if (length(errs) > 0) stats::median(errs) else 0,
      p95_abs_err = if (length(errs) > 0) stats::quantile(errs, 0.95) else 0,
      reference = "manual_hypergeometric",
      notes = sprintf("n_edges=%d,pval_ok=%s,obs_ok=%s,class_ok=%s",
                       nrow(scores), pval_valid, obs_valid, class_valid)
    )

    list(pval_valid = pval_valid, obs_valid = obs_valid,
         class_valid = class_valid, n_failed = n_failed,
         max_err = max_err, n_matched = sum(matched), config_idx = i)
  })

  report$write_csv("hon")

  # Assertions
  lapply(results, function(r) {
    expect_true(r$pval_valid,
                info = sprintf("Config %d: p-values out of [0,1]",
                               r$config_idx))
    expect_true(r$obs_valid,
                info = sprintf("Config %d: observed counts invalid",
                               r$config_idx))
    expect_true(r$class_valid,
                info = sprintf("Config %d: invalid anomaly classification",
                               r$config_idx))
    expect_true(r$n_failed == 0L,
                info = sprintf("Config %d: %d/%d p-values mismatch manual (max_err=%.2e)",
                               r$config_idx, r$n_failed, r$n_matched,
                               r$max_err))
  })
})


# ---- 3. HONEM embedding properties (50 configs) ----

test_that("HONEM embeddings satisfy SVD invariants across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()

  TOL <- 1e-10
  report <- equiv_report()

  results <- lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    hon <- build_hon(seqs, max_order = 2)
    n_nodes <- hon$n_nodes
    emb_dim <- min(10L, n_nodes - 1L)

    if (n_nodes < 2L) {
      return(list(sv_sorted = TRUE, ev_valid = TRUE, dim_ok = TRUE,
                  no_nan = TRUE, config_idx = i, n_nodes = n_nodes))
    }

    honem <- build_honem(hon, dim = emb_dim)

    # (a) Singular values sorted descending
    sv <- honem$singular_values
    sv_sorted <- all(diff(sv) <= TOL)

    # (b) Explained variance in [0, 1]
    ev_valid <- honem$explained_variance >= 0 - TOL &&
                honem$explained_variance <= 1 + TOL

    # (c) Embeddings have correct dimensions (n_nodes x dim)
    actual_dim <- dim(honem$embeddings)
    dim_ok <- actual_dim[1L] == n_nodes && actual_dim[2L] == emb_dim

    # (d) No NaN or Inf in embeddings
    no_nan <- !any(is.nan(honem$embeddings) | is.infinite(honem$embeddings))

    # Compute error metrics from singular value sorting
    sv_diffs <- diff(sv)
    sv_violations <- sv_diffs[sv_diffs > TOL]

    report$log(
      func = "build_honem",
      config = sprintf("actors=%d,states=%d,len=%d,seed=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed),
      n_checked = length(sv) + 1L + 1L + length(honem$embeddings),
      n_failed = sum(!sv_sorted, !ev_valid, !dim_ok, !no_nan),
      max_abs_err = if (length(sv_violations) > 0) max(sv_violations) else 0,
      mean_abs_err = if (length(sv_violations) > 0) mean(sv_violations) else 0,
      median_abs_err = 0,
      p95_abs_err = 0,
      reference = "svd_invariants",
      notes = sprintf("n_nodes=%d,dim=%d,ev=%.4f,sv_sorted=%s,dim_ok=%s,no_nan=%s",
                       n_nodes, emb_dim, honem$explained_variance,
                       sv_sorted, dim_ok, no_nan)
    )

    list(sv_sorted = sv_sorted, ev_valid = ev_valid,
         dim_ok = dim_ok, no_nan = no_nan,
         config_idx = i, n_nodes = n_nodes)
  })

  report$write_csv("hon")

  # Assertions
  lapply(results, function(r) {
    expect_true(r$sv_sorted,
                info = sprintf("Config %d: singular values not sorted descending",
                               r$config_idx))
    expect_true(r$ev_valid,
                info = sprintf("Config %d: explained_variance out of [0,1]",
                               r$config_idx))
    expect_true(r$dim_ok,
                info = sprintf("Config %d: embedding dimensions incorrect (n_nodes=%d)",
                               r$config_idx, r$n_nodes))
    expect_true(r$no_nan,
                info = sprintf("Config %d: NaN/Inf in embeddings",
                               r$config_idx))
  })
})


# ---- 4. MOGen order selection and log-likelihood monotonicity (50 configs) ----

test_that("MOGen log-likelihoods are monotonically increasing across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()

  TOL <- 1e-10
  report <- equiv_report()

  results <- lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    mogen <- build_mogen(seqs, max_order = 3)

    ll <- mogen$log_likelihood
    dof <- mogen$dof
    aic <- mogen$aic
    bic <- mogen$bic
    opt <- mogen$optimal_order

    # (a) Log-likelihoods are monotonically non-decreasing
    ll_diffs <- diff(ll)
    ll_mono <- all(ll_diffs >= -TOL)

    # (b) DOF is monotonically increasing
    dof_diffs <- diff(dof)
    dof_mono <- all(dof_diffs >= 0L)

    # (c) AIC and BIC are finite
    aic_finite <- all(is.finite(aic))
    bic_finite <- all(is.finite(bic))

    # (d) optimal_order is in [0, max_order]
    max_order_actual <- length(ll) - 1L
    opt_valid <- opt >= 0L && opt <= max_order_actual

    # Track monotonicity violations
    ll_violations <- ll_diffs[ll_diffs < -TOL]
    max_violation <- if (length(ll_violations) > 0) max(abs(ll_violations)) else 0

    report$log(
      func = "build_mogen",
      config = sprintf("actors=%d,states=%d,len=%d,seed=%d",
                       cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed),
      n_checked = length(ll) + length(dof) + length(aic) + length(bic) + 1L,
      n_failed = sum(!ll_mono, !dof_mono, !aic_finite, !bic_finite, !opt_valid),
      max_abs_err = max_violation,
      mean_abs_err = if (length(ll_violations) > 0) mean(abs(ll_violations)) else 0,
      median_abs_err = 0,
      p95_abs_err = 0,
      reference = "statistical_invariants",
      notes = sprintf("ll_mono=%s,dof_mono=%s,aic_finite=%s,bic_finite=%s,opt=%d",
                       ll_mono, dof_mono, aic_finite, bic_finite, opt)
    )

    list(ll_mono = ll_mono, dof_mono = dof_mono,
         aic_finite = aic_finite, bic_finite = bic_finite,
         opt_valid = opt_valid, config_idx = i)
  })

  report$write_csv("hon")

  # Assertions
  lapply(results, function(r) {
    expect_true(r$ll_mono,
                info = sprintf("Config %d: log-likelihoods not monotonically increasing",
                               r$config_idx))
    expect_true(r$dof_mono,
                info = sprintf("Config %d: DOF not monotonically increasing",
                               r$config_idx))
    expect_true(r$aic_finite,
                info = sprintf("Config %d: non-finite AIC values",
                               r$config_idx))
    expect_true(r$bic_finite,
                info = sprintf("Config %d: non-finite BIC values",
                               r$config_idx))
    expect_true(r$opt_valid,
                info = sprintf("Config %d: optimal_order out of valid range",
                               r$config_idx))
  })
})
