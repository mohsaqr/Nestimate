# ---- Equivalence: predict_links() vs clean-room matrix-algebra refs ----
#
# Validates Nestimate::predict_links() scores for all six structural
# similarity methods against independent base-R reimplementations on
# random undirected binary graphs. Methods: (1) common_neighbors as
# (A %*% A)[i,j] (2) jaccard CN / (k_i + k_j - CN) (3) adamic_adar
# sum_z A[i,z]A[j,z]/log(k_z) with k_z>1 (4) resource_allocation
# sum_z A[i,z]A[j,z]/k_z with k_z>0 (5) preferential_attachment k_i*k_j
# (6) katz closed-form (I - beta*A)^{-1} - I with beta = 0.05. Cross-
# checked where applicable against igraph::similarity (jaccard, adamic
# invlogweighted) as a secondary reference. Clean-room is primary. We
# pass weighted=FALSE to predict_links so M = A matches the binary-A
# reference formulas byte-for-byte.

skip_equiv_tests <- function() {
  run <- Sys.getenv("NESTIMATE_EQUIV_TESTS", unset = "false")
  if (!identical(run, "true"))
    skip("Equivalence tests skipped (set NESTIMATE_EQUIV_TESTS=true)")
}

TOL <- 1e-8
KATZ_BETA <- 0.05

# ---- Random undirected binary graph generator ----
random_undirected_graph <- function(n, density, seed) {
  set.seed(seed)
  A <- matrix(0L, n, n)
  upper_idx <- which(upper.tri(A))
  n_edges <- max(1L, round(length(upper_idx) * density))
  chosen <- sample(upper_idx, n_edges)
  A[chosen] <- 1L
  A <- A + t(A)
  dimnames(A) <- list(paste0("n", seq_len(n)), paste0("n", seq_len(n)))
  A
}

# ---- Clean-room reference formulas (all ≤ 5 lines, vectorized) ----
ref_common_neighbors <- function(A) { s <- A %*% A; diag(s) <- 0; s }

ref_jaccard <- function(A) {
  CN <- A %*% A; k <- rowSums(A); U <- outer(k, k, "+") - CN
  s <- ifelse(U > 0, CN / U, 0); diag(s) <- 0; s
}

ref_adamic_adar <- function(A) {
  k <- rowSums(A); w <- ifelse(k > 1, 1 / log(k), 0)
  s <- A %*% diag(w) %*% A; diag(s) <- 0; s
}

ref_resource_allocation <- function(A) {
  k <- rowSums(A); w <- ifelse(k > 0, 1 / k, 0)
  s <- A %*% diag(w) %*% A; diag(s) <- 0; s
}

ref_preferential_attachment <- function(A) {
  k <- rowSums(A); s <- outer(k, k, "*"); diag(s) <- 0; s
}

ref_katz <- function(A, beta) {
  n <- nrow(A); I <- diag(n); s <- solve(I - beta * A) - I
  diag(s) <- 0; s
}

# ---- Wrap a binary adjacency as a netobject (undirected) ----
wrap_binary_netobject <- function(A) {
  n <- nrow(A)
  states <- rownames(A)
  nodes_df <- data.frame(
    id = seq_len(n), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )
  # edges data.frame (not strictly needed for predict_links, but mimic netobject)
  idx <- which(upper.tri(A) & A != 0, arr.ind = TRUE)
  edges_df <- data.frame(
    from = idx[, 1], to = idx[, 2],
    weight = A[idx], stringsAsFactors = FALSE
  )
  structure(
    list(
      data = NULL, weights = A * 1.0, nodes = nodes_df, edges = edges_df,
      directed = FALSE, meta = list(),
      method = "binary_ref", params = list(), scaling = "none",
      threshold = 0, level = NULL, n_nodes = n, n_edges = nrow(edges_df),
      node_groups = NULL, metadata = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}

# ---- Compare two matrices by zeroing diagonal and coercing to numeric ----
compare_mats <- function(ours, theirs) {
  ours <- as.matrix(ours); theirs <- as.matrix(theirs)
  diag(ours) <- 0; diag(theirs) <- 0
  abs(ours - theirs)
}

# ---- 100-config sweep (4 sizes * 3 densities * 9 seeds = 108 configs) ----

test_that("predict_links scores match clean-room refs across 6 methods (100+ configs)", {
  skip_equiv_tests()

  rep <- equiv_report()

  sizes    <- c(8L, 12L, 16L, 20L)
  dens_set <- c(0.2, 0.3, 0.4)
  seeds    <- seq_len(9L)

  methods_list <- c("common_neighbors", "resource_allocation", "adamic_adar",
                    "jaccard", "preferential_attachment", "katz")

  configs <- expand.grid(n = sizes, density = dens_set, seed = seeds,
                         KEEP.OUT.ATTRS = FALSE)

  n_fail_total <- 0L
  n_check_total <- 0L
  max_delta_by_method <- setNames(rep(0, length(methods_list)), methods_list)

  for (cfg_i in seq_len(nrow(configs))) {
    n       <- configs$n[cfg_i]
    density <- configs$density[cfg_i]
    seed    <- configs$seed[cfg_i]

    A   <- random_undirected_graph(n, density, seed = seed * 1000L + n * 10L + round(density * 100))
    net <- wrap_binary_netobject(A)

    # Single predict_links call with all 6 methods, binary (weighted = FALSE),
    # and an explicit Katz beta. Auto-damping is disabled by passing katz_damping.
    pred <- predict_links(net,
                          methods = methods_list,
                          weighted = FALSE,
                          exclude_existing = FALSE,
                          include_self = FALSE,
                          katz_damping = KATZ_BETA)

    # Reference scores
    refs <- list(
      common_neighbors        = ref_common_neighbors(A),
      resource_allocation     = ref_resource_allocation(A),
      adamic_adar             = ref_adamic_adar(A),
      jaccard                 = ref_jaccard(A),
      preferential_attachment = ref_preferential_attachment(A),
      katz                    = ref_katz(A, KATZ_BETA)
    )

    # Per-method delta
    deltas <- lapply(methods_list, function(m) {
      compare_mats(pred$scores[[m]], refs[[m]])
    })
    names(deltas) <- methods_list

    for (m in methods_list) {
      d     <- deltas[[m]]
      nfail <- sum(d > TOL)
      n_fail_total  <- n_fail_total + nfail
      n_check_total <- n_check_total + length(d)
      if (max(d) > max_delta_by_method[[m]]) max_delta_by_method[[m]] <- max(d)

      rep$log(
        func           = sprintf("predict_links[%s]", m),
        config         = sprintf("n=%d dens=%.2f seed=%d", n, density, seed),
        n_checked      = length(d),
        n_failed       = nfail,
        max_abs_err    = max(d),
        mean_abs_err   = mean(d),
        median_abs_err = stats::median(d),
        p95_abs_err    = stats::quantile(d, 0.95),
        reference      = "clean-room base-R matrix algebra",
        notes          = if (m == "katz")
                           sprintf("closed-form (I - beta*A)^-1 - I, beta=%.3f", KATZ_BETA)
                         else "binary undirected (weighted=FALSE)"
      )
    }

    expect_true(
      all(vapply(deltas, function(d) max(d) <= TOL, logical(1))),
      info = sprintf("n=%d dens=%.2f seed=%d: max deltas=%s",
                     n, density, seed,
                     paste(sprintf("%s=%.2e", methods_list,
                                   vapply(deltas, max, numeric(1))),
                           collapse = " "))
    )
  }

  rep$write_csv("link_prediction")
  message(sprintf(
    "link_prediction equivalence: %d configs * 6 methods = %d method-checks, %d cells, %d failed",
    nrow(configs), nrow(configs) * length(methods_list),
    n_check_total, n_fail_total
  ))
  message("Max delta by method:")
  for (m in methods_list) {
    message(sprintf("  %-25s %.3e", m, max_delta_by_method[[m]]))
  }
})

# ---- Secondary cross-check: igraph::similarity for jaccard and invlogweighted ----

test_that("predict_links jaccard and adamic_adar match igraph::similarity on 10 graphs", {
  skip_equiv_tests()
  skip_if_not_installed("igraph")

  rep <- equiv_report()

  configs <- expand.grid(n = c(10L, 15L, 20L),
                         density = c(0.25, 0.4),
                         seed = 1:4,
                         KEEP.OUT.ATTRS = FALSE)[1:12, ]

  for (cfg_i in seq_len(nrow(configs))) {
    n <- configs$n[cfg_i]; density <- configs$density[cfg_i]; seed <- configs$seed[cfg_i]
    A <- random_undirected_graph(n, density, seed = seed * 773L + n)
    net <- wrap_binary_netobject(A)

    pred <- predict_links(net,
                          methods = c("jaccard", "adamic_adar"),
                          weighted = FALSE, exclude_existing = FALSE,
                          include_self = FALSE)

    g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)

    # igraph jaccard: same definition as ours (|N(i) ∩ N(j)| / |N(i) ∪ N(j)|)
    ig_jac <- igraph::similarity(g, method = "jaccard")
    diag(ig_jac) <- 0

    # igraph invlogweighted ~ adamic_adar
    ig_aa <- igraph::similarity(g, method = "invlogweighted")
    diag(ig_aa) <- 0

    d_jac <- compare_mats(pred$scores$jaccard, ig_jac)
    d_aa  <- compare_mats(pred$scores$adamic_adar, ig_aa)

    rep$log(
      func = "predict_links[jaccard] vs igraph",
      config = sprintf("n=%d dens=%.2f seed=%d", n, density, seed),
      n_checked = length(d_jac), n_failed = sum(d_jac > TOL),
      max_abs_err = max(d_jac), mean_abs_err = mean(d_jac),
      median_abs_err = stats::median(d_jac),
      p95_abs_err = stats::quantile(d_jac, 0.95),
      reference = "igraph::similarity(method='jaccard')",
      notes = "secondary cross-check"
    )
    rep$log(
      func = "predict_links[adamic_adar] vs igraph",
      config = sprintf("n=%d dens=%.2f seed=%d", n, density, seed),
      n_checked = length(d_aa), n_failed = sum(d_aa > TOL),
      max_abs_err = max(d_aa), mean_abs_err = mean(d_aa),
      median_abs_err = stats::median(d_aa),
      p95_abs_err = stats::quantile(d_aa, 0.95),
      reference = "igraph::similarity(method='invlogweighted')",
      notes = "secondary cross-check"
    )

    expect_true(max(d_jac) <= TOL,
                info = sprintf("jaccard n=%d dens=%.2f seed=%d: max=%.2e",
                               n, density, seed, max(d_jac)))
    # adamic_adar: igraph may differ; expect with larger tolerance
    expect_true(max(d_aa) <= 1e-6 || max(d_aa) <= TOL,
                info = sprintf("adamic_adar n=%d dens=%.2f seed=%d: max=%.2e vs igraph",
                               n, density, seed, max(d_aa)))
  }

  rep$write_csv("link_prediction_igraph")
  message(rep$summary())
})
