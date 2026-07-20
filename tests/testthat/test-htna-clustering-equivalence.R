# Cell-level equivalence tests for clustering an HTNA-shaped netobject.
# Nestimate deliberately does not import htna; these fixtures encode the
# public HTNA contract directly and compare against ordinary netobjects and
# independently rebuilt per-assignment networks.

.htna_equiv_sequences <- function(n, width, n_human, n_ai, seed) {
  set.seed(seed)
  human <- paste0("H", seq_len(n_human))
  ai <- paste0("A", seq_len(n_ai))
  states <- c(human, ai)
  x <- matrix(
    sample(states, n * width, replace = TRUE),
    nrow = n,
    dimnames = list(paste0("S", seq_len(n)), paste0("T", seq_len(width)))
  )
  # Make state coverage deterministic rather than trusting a random draw.
  x[seq_along(states)] <- states
  as.data.frame(x, stringsAsFactors = FALSE)
}

.as_equiv_htna <- function(net) {
  labels <- as.character(net$nodes$label)
  groups <- ifelse(startsWith(labels, "H"), "Human", "AI")
  net$nodes$groups <- factor(groups, levels = c("Human", "AI"))
  net$node_groups <- data.frame(
    node = labels,
    group = net$nodes$groups,
    stringsAsFactors = FALSE
  )
  net$actor_levels <- c("Human", "AI")
  class(net) <- unique(c("htna", class(net)))
  net
}

.expect_equiv_matrix <- function(x, y, tolerance = 1e-13, info = NULL) {
  expect_setequal(rownames(x), rownames(y))
  expect_setequal(colnames(x), colnames(y))
  y <- y[rownames(x), colnames(x), drop = FALSE]
  expect_equal(unname(x), unname(y), tolerance = tolerance, info = info)
}

.expect_equiv_network <- function(x, y, info = NULL) {
  .expect_equiv_matrix(x$weights, y$weights, info = paste(info, "weights"))
  if (is.null(x$initial) || is.null(y$initial)) {
    expect_identical(is.null(x$initial), is.null(y$initial), info = info)
  } else {
    expect_setequal(names(x$initial), names(y$initial))
    expect_equal(unname(x$initial), unname(y$initial[names(x$initial)]),
                 tolerance = 1e-13, info = paste(info, "initial"))
  }
  expect_identical(x$method, y$method, info = info)
  expect_setequal(as.character(x$nodes$label), as.character(y$nodes$label))
  x_data <- x$data
  y_data <- y$data
  rownames(x_data) <- NULL
  rownames(y_data) <- NULL
  expect_equal(x_data, y_data, info = paste(info, "data"))
}

.expect_equiv_actor_partition <- function(net, expected_levels = c("Human", "AI")) {
  expect_s3_class(net, "htna")
  expect_identical(net$actor_levels, expected_levels)
  expect_identical(as.character(net$node_groups$node),
                   as.character(net$nodes$label))
  observed <- stats::setNames(
    as.character(net$node_groups$group),
    as.character(net$node_groups$node)
  )
  expected <- ifelse(startsWith(names(observed), "H"), "Human", "AI")
  expect_identical(unname(observed), expected)
  expect_identical(as.character(net$nodes$groups), unname(observed))
}

set.seed(20260720)
.htna_distance_configs <- lapply(seq_len(12L), function(i) {
  list(
    n = sample(c(18L, 24L, 32L), 1L),
    width = sample(c(6L, 8L, 11L), 1L),
    n_human = sample(2:4, 1L),
    n_ai = sample(2:4, 1L),
    k = sample(2:3, 1L),
    cluster_method = sample(c("pam", "ward.D2", "complete"), 1L),
    seed = sample.int(100000L, 1L)
  )
})

test_that("HTNA distance clustering is numerically identical to plain and manual paths", {
  methods <- c("relative", "frequency", "attention")
  checked_cells <- 0L

  for (cfg in .htna_distance_configs) {
    seqs <- .htna_equiv_sequences(
      cfg$n, cfg$width, cfg$n_human, cfg$n_ai, cfg$seed
    )
    plain <- build_network(seqs, method = "relative")
    htna <- .as_equiv_htna(plain)

    cl_plain <- build_clusters(
      plain, k = cfg$k, method = cfg$cluster_method, seed = cfg$seed
    )
    cl_htna <- build_clusters(
      htna, k = cfg$k, method = cfg$cluster_method, seed = cfg$seed
    )

    expect_identical(cl_htna$assignments, cl_plain$assignments)
    expect_identical(cl_htna$sizes, cl_plain$sizes)
    expect_equal(cl_htna$silhouette, cl_plain$silhouette, tolerance = 1e-14)
    expect_equal(as.matrix(cl_htna$distance), as.matrix(cl_plain$distance),
                 tolerance = 1e-14)
    expect_identical(cl_htna$htna_partition$actor_levels, htna$actor_levels)

    for (network_method in methods) {
      grouped_htna <- build_network(cl_htna, method = network_method)
      grouped_plain <- build_network(cl_plain, method = network_method)

      expect_s3_class(grouped_htna, "htna_group")
      expect_identical(attr(grouped_htna, "group_col"),
                       attr(grouped_plain, "group_col"))
      expect_identical(
        attr(grouped_htna, "clustering")$assignments,
        attr(grouped_plain, "clustering")$assignments
      )

      for (cluster_id in seq_len(cfg$k)) {
        manual_data <- cl_plain$data[
          cl_plain$assignments == cluster_id, , drop = FALSE
        ]
        manual <- build_network(manual_data, method = network_method)
        info <- paste0(
          "seed=", cfg$seed, " method=", network_method,
          " cluster=", cluster_id
        )

        .expect_equiv_network(grouped_htna[[cluster_id]],
                              grouped_plain[[cluster_id]], info)
        .expect_equiv_network(grouped_htna[[cluster_id]], manual, info)
        .expect_equiv_actor_partition(grouped_htna[[cluster_id]])
        checked_cells <- checked_cells + length(manual$weights)
      }
    }
  }

  expect_gte(checked_cells, 1500L)
})

test_that("cluster_network HTNA preservation changes metadata only", {
  checked_cells <- 0L
  configs <- .htna_distance_configs[seq_len(8L)]

  for (cfg in configs) {
    seqs <- .htna_equiv_sequences(
      cfg$n, cfg$width, cfg$n_human, cfg$n_ai, cfg$seed
    )
    plain <- build_network(seqs, method = "relative")
    htna <- .as_equiv_htna(plain)

    for (network_method in c("relative", "frequency", "attention")) {
      grouped_plain <- cluster_network(
        plain,
        k = cfg$k,
        cluster_by = cfg$cluster_method,
        seed = cfg$seed,
        method = network_method
      )
      grouped_htna <- cluster_network(
        htna,
        k = cfg$k,
        cluster_by = cfg$cluster_method,
        seed = cfg$seed,
        method = network_method
      )

      expect_s3_class(grouped_htna, "htna_group")
      expect_identical(attr(grouped_htna, "clustering")$assignments,
                       attr(grouped_plain, "clustering")$assignments)
      for (cluster_id in seq_len(cfg$k)) {
        .expect_equiv_network(
          grouped_htna[[cluster_id]], grouped_plain[[cluster_id]],
          paste(cfg$seed, network_method, cluster_id)
        )
        .expect_equiv_actor_partition(grouped_htna[[cluster_id]])
        checked_cells <- checked_cells + length(grouped_plain[[cluster_id]]$weights)
      }
    }
  }

  expect_gte(checked_cells, 900L)
})

test_that("MMM clustering of HTNA is numerically identical to plain input", {
  checked_cells <- 0L
  seeds <- c(104L, 209L, 313L, 419L, 521L, 631L)

  for (seed in seeds) {
    seqs <- .htna_equiv_sequences(30L, 9L, 3L, 3L, seed)
    plain <- build_network(seqs, method = "relative")
    htna <- .as_equiv_htna(plain)

    args <- list(k = 2L, n_starts = 2L, max_iter = 40L, seed = seed)
    grouped_plain <- do.call(cluster_mmm, c(list(data = plain), args))
    grouped_htna <- do.call(cluster_mmm, c(list(data = htna), args))
    fit_plain <- attr(grouped_plain, "clustering")
    fit_htna <- attr(grouped_htna, "clustering")

    expect_s3_class(grouped_htna, "htna_group")
    expect_identical(fit_htna$assignments, fit_plain$assignments)
    expect_equal(fit_htna$posterior, fit_plain$posterior, tolerance = 1e-12)
    expect_equal(fit_htna$mixing, fit_plain$mixing, tolerance = 1e-12)
    expect_equal(fit_htna$log_likelihood, fit_plain$log_likelihood,
                 tolerance = 1e-12)
    expect_equal(c(fit_htna$AIC, fit_htna$BIC, fit_htna$ICL),
                 c(fit_plain$AIC, fit_plain$BIC, fit_plain$ICL),
                 tolerance = 1e-12)
    expect_identical(fit_htna$converged, fit_plain$converged)

    for (cluster_id in seq_along(grouped_plain)) {
      .expect_equiv_network(
        grouped_htna[[cluster_id]], grouped_plain[[cluster_id]],
        paste("MMM", seed, cluster_id)
      )
      .expect_equiv_actor_partition(grouped_htna[[cluster_id]])
      checked_cells <- checked_cells + length(grouped_plain[[cluster_id]]$weights)
    }

    mmm_fit <- do.call(build_mmm, c(list(data = htna), args))
    rebuilt <- build_network(mmm_fit)
    expect_s3_class(rebuilt, "htna_group")
    expect_identical(attr(rebuilt, "clustering")$assignments,
                     mmm_fit$assignments)
    for (cluster_id in seq_along(rebuilt)) {
      .expect_equiv_network(rebuilt[[cluster_id]], mmm_fit$models[[cluster_id]],
                            paste("build_network(net_mmm)", seed, cluster_id))
      .expect_equiv_actor_partition(rebuilt[[cluster_id]])
    }
  }

  expect_gte(checked_cells, 300L)
})

test_that("HTNA actor levels survive clusters containing only one actor", {
  human <- data.frame(
    T1 = rep("H1", 8), T2 = rep("H2", 8), T3 = rep("H1", 8),
    stringsAsFactors = FALSE
  )
  ai <- data.frame(
    T1 = rep("A1", 8), T2 = rep("A2", 8), T3 = rep("A1", 8),
    stringsAsFactors = FALSE
  )
  seqs <- rbind(human, ai)
  htna <- .as_equiv_htna(build_network(seqs, method = "relative"))
  clustering <- build_clusters(htna, k = 2, seed = 1)
  clustering$assignments <- c(rep(1L, nrow(human)), rep(2L, nrow(ai)))
  clustering$sizes <- c(`1` = nrow(human), `2` = nrow(ai))

  grouped <- build_network(clustering)

  expect_s3_class(grouped, "htna_group")
  expect_identical(attr(grouped, "actor_levels"), c("Human", "AI"))
  expect_setequal(as.character(grouped[[1]]$node_groups$group), "Human")
  expect_setequal(as.character(grouped[[2]]$node_groups$group), "AI")
  expect_identical(grouped[[1]]$actor_levels, c("Human", "AI"))
  expect_identical(grouped[[2]]$actor_levels, c("Human", "AI"))
  .expect_equiv_network(
    grouped[[1]], build_network(human, method = "relative"), "Human-only"
  )
  .expect_equiv_network(
    grouped[[2]], build_network(ai, method = "relative"), "AI-only"
  )
})
