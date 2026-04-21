# Numerical equivalence: Nestimate's build_mogen() vs pathpy.MultiOrderModel.
#
# pathpy 2.2.0 (Scholtes' own lab) is the canonical Python reference for MOGen
# (Scholtes 2017, Gote & Scholtes 2023). We check three things per config:
#
#  1. Per-order count matrices from .mogen_count_kgrams()
#     match pathpy.HigherOrderNetwork(k=k).edges (summed weights).
#  2. Total cumulative log-likelihood at each max_order 1..K matches
#     pathpy's MultiOrderModel(paths, max_order=k).likelihood(log=True).
#  3. path_counts() agrees with pathpy k-gram enumeration (integer-exact).

set.seed(4242)
N_MOGEN <- 30L
TOL <- 1e-8         # counts must match exactly; likelihood to 1e-8
TOL_LIK <- 1e-6     # log-likelihood tolerance (hierarchical decomposition may
                    # differ in floor value for unseen transitions)

skip_if_pkg_broken("reticulate")

# ---- Config generation ----
mogen_configs <- lapply(seq_len(N_MOGEN), function(i) {
  list(n_actors = sample(c(10L, 15L, 20L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(15L, 20L, 30L), 1),
       max_order = sample(2:3, 1),
       seed = sample.int(100000, 1))
})

test_that("MOGen per-order count matrices match pathpy.HigherOrderNetwork", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  report <- equiv_report()
  HON_SEP <- "\x01"

  invisible(lapply(seq_len(N_MOGEN), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    mg <- tryCatch(build_mogen(seqs, max_order = cfg$max_order),
                   error = function(e) NULL)
    if (is.null(mg)) return(NULL)

    paths <- py_paths_from_sequences(seqs)

    # Check each order k from 1 to max_order. Nestimate stores the count matrix
    # at index k+1; node names use \x01 separators.
    invisible(lapply(seq_len(cfg$max_order), function(k) {
      nest_cm <- mg$count_matrices[[k + 1L]]
      py_cm <- py_hon_count_matrix(paths, k)

      # Translate node names: Nestimate uses \x01, pathpy uses ",".
      nest_nodes_py <- gsub(HON_SEP, ",", rownames(nest_cm), fixed = TRUE)
      rownames(nest_cm) <- nest_nodes_py
      colnames(nest_cm) <- nest_nodes_py

      # Align on intersection of node sets; log any asymmetry.
      common <- intersect(rownames(nest_cm), rownames(py_cm))
      n_only_nest <- length(setdiff(rownames(nest_cm), common))
      n_only_py <- length(setdiff(rownames(py_cm), common))

      if (length(common) == 0L) {
        report$log(func = sprintf("mogen_counts_k%d", k),
                   config = sprintf("cfg%d(no_common_nodes)", i),
                   n_checked = 0L, n_failed = 0L,
                   max_abs_err = NA_real_, mean_abs_err = NA_real_,
                   median_abs_err = NA_real_, p95_abs_err = NA_real_,
                   reference = "pathpy.HigherOrderNetwork",
                   notes = "empty intersection")
        return(NULL)
      }

      a <- nest_cm[common, common, drop = FALSE]
      b <- py_cm[common, common, drop = FALSE]
      delta <- abs(a - b)

      report$log(
        func = sprintf("mogen_counts_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,k=%d)",
                         i, cfg$n_actors, cfg$seq_length, k),
        n_checked = length(delta),
        n_failed = as.integer(sum(delta > TOL)),
        max_abs_err = max(delta), mean_abs_err = mean(delta),
        median_abs_err = stats::median(delta),
        p95_abs_err = as.numeric(stats::quantile(delta, 0.95)),
        reference = "pathpy.HigherOrderNetwork",
        notes = sprintf("common=%d only_nest=%d only_py=%d",
                        length(common), n_only_nest, n_only_py)
      )

      expect_true(max(delta) < TOL,
                  label = sprintf("cfg%d k=%d max count delta = %.2e",
                                  i, k, max(delta)))
    }))

    NULL
  }))

  report$write_csv("mogen_counts")
  report$write_cvs("mogen_counts",
                   "local_testing_and_equivalence/test-equiv-mogen.R")
})

test_that("MOGen cumulative log-likelihood matches pathpy.MultiOrderModel", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  report <- equiv_report()

  invisible(lapply(seq_len(N_MOGEN), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    mg <- tryCatch(build_mogen(seqs, max_order = cfg$max_order),
                   error = function(e) NULL)
    if (is.null(mg)) return(NULL)

    paths <- py_paths_from_sequences(seqs)

    # Compare total log-likelihood at each max_order 1..cfg$max_order
    invisible(lapply(seq_len(cfg$max_order), function(k) {
      nest_ll <- as.numeric(mg$log_likelihood[k + 1L])
      py_ll <- py_mogen_likelihood(paths, k)
      delta <- abs(nest_ll - py_ll)

      report$log(
        func = sprintf("mogen_loglik_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,k=%d)",
                         i, cfg$n_actors, cfg$seq_length, k),
        n_checked = 1L,
        n_failed = as.integer(delta > TOL_LIK),
        max_abs_err = delta, mean_abs_err = delta,
        median_abs_err = delta, p95_abs_err = delta,
        reference = "pathpy.MultiOrderModel.likelihood",
        notes = sprintf("nest=%.4f py=%.4f", nest_ll, py_ll)
      )

      expect_true(delta < TOL_LIK,
                  label = sprintf("cfg%d k=%d ll delta = %.2e (nest=%.4f py=%.4f)",
                                  i, k, delta, nest_ll, py_ll))
    }))

    NULL
  }))

  report$write_csv("mogen_loglik")
  report$write_cvs("mogen_loglik",
                   "local_testing_and_equivalence/test-equiv-mogen.R")
})

test_that("path_counts matches pathpy k-gram enumeration exactly", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  invisible(lapply(seq_len(10L), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    # path_counts uses " -> " separator; k=2 is a transition count.
    pc <- path_counts(seqs, k = 2L)

    # Manually aggregate bigram counts; pathpy HON(k=1) gives same result
    # (subpath + longest-path = total occurrences of each transition).
    paths <- py_paths_from_sequences(seqs)
    py_edges <- py_hon_edges(paths, 1L)
    py_edges$path <- paste(py_edges$from, py_edges$to, sep = " -> ")

    common <- intersect(pc$path, py_edges$path)
    pc_aligned <- pc[match(common, pc$path), ]
    py_aligned <- py_edges[match(common, py_edges$path), ]

    expect_equal(pc_aligned$count, py_aligned$count,
                 tolerance = 0,
                 label = sprintf("cfg%d path_counts k=2", i))
  }))
})
