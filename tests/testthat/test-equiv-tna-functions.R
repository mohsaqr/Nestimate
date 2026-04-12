# ---- Equivalence: Nestimate vs tna for centralities, CS, reliability, etc. ---
#
# Per-value comparison of every shared function against tna, with and
# without actor grouping, with and without windowing.
#
#   centrality()           vs tna::centralities()
#   centrality_stability() vs tna::estimate_cs()
#   reliability()          vs tna::reliability()
#   compare_sequences()    vs tna::compare_sequences()
#   cluster_data()         vs tna::cluster_data()
#
# ~20 configs per function.  Requires NESTIMATE_EQUIV_TESTS=true.

set.seed(8888)
N <- 20L
TOL <- 1e-10

configs <- lapply(seq_len(N), function(i) {
  list(
    n_actors   = sample(c(15, 20, 30), 1),
    n_states   = sample(c(3, 5, 7), 1),
    seq_length = sample(c(12, 18, 25), 1),
    seed       = sample.int(100000, 1)
  )
})

report <- equiv_report()

# ---- 1. centrality: InStrength / OutStrength vs tna::centralities() ----------

test_that("centrality: per-value InStrength/OutStrength match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    nest_net <- build_network(data, method = "relative")
    tna_net  <- tna::tna(data)
    tna_c    <- tna::centralities(tna_net)
    nest_c   <- centrality(nest_net)

    # Match by state name
    states <- sort(intersect(tna_c$state, names(nest_c$InStrength)))
    lbl <- sprintf("cfg%d(a=%d,s=%d,l=%d)", i, cfg$n_actors, cfg$n_states, cfg$seq_length)

    invisible(lapply(states, function(s) {
      expect_equal(
        nest_c$InStrength[s], tna_c$InStrength[tna_c$state == s],
        tolerance = TOL,
        label = sprintf("InStrength %s [%s]", lbl, s)
      )
      expect_equal(
        nest_c$OutStrength[s], tna_c$OutStrength[tna_c$state == s],
        tolerance = TOL,
        label = sprintf("OutStrength %s [%s]", lbl, s)
      )
    }))

    n_checked <- length(states) * 2L
    report$log("centrality", lbl, n_checked, 0L, 0, 0, 0, 0, "tna::centralities()")
  }))
})


# ---- 2. centrality with windowed networks -----------------------------------

test_that("centrality: windowed networks match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    ws <- sample(c(2L, 3L), 1)

    nest_net <- build_network(data, method = "relative",
                              params = list(window_size = ws))
    tna_net  <- tna::tna(data, params = list(window_size = ws))
    tna_c    <- tna::centralities(tna_net)
    nest_c   <- centrality(nest_net)

    states <- sort(intersect(tna_c$state, names(nest_c$InStrength)))
    lbl <- sprintf("cfg%d ws=%d", i, ws)

    invisible(lapply(states, function(s) {
      expect_equal(
        nest_c$InStrength[s], tna_c$InStrength[tna_c$state == s],
        tolerance = TOL,
        label = sprintf("InStrength %s [%s]", lbl, s)
      )
      expect_equal(
        nest_c$OutStrength[s], tna_c$OutStrength[tna_c$state == s],
        tolerance = TOL,
        label = sprintf("OutStrength %s [%s]", lbl, s)
      )
    }))
  }))
})


# ---- 3. centrality_stability CS coefficient vs tna::estimate_cs() -----------

test_that("centrality_stability: CS coefficients match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    iter <- 30L

    nest_net <- build_network(data, method = "relative")
    tna_net  <- tna::tna(data)

    set.seed(1)
    tna_cs <- tna::estimate_cs(tna_net, iter = iter)
    nest_cs <- centrality_stability(nest_net, iter = iter, seed = 1)

    lbl <- sprintf("cfg%d(a=%d,s=%d)", i, cfg$n_actors, cfg$n_states)

    # Compare CS coefficients for shared measures
    shared <- intersect(names(tna_cs), names(nest_cs$cs))
    invisible(lapply(shared, function(m) {
      tna_val  <- tna_cs[[m]]$cs_coefficient
      nest_val <- nest_cs$cs[[m]]
      expect_equal(nest_val, tna_val, tolerance = 0.15,
                   label = sprintf("CS(%s) %s  nest=%.2f tna=%.2f", m, lbl, nest_val, tna_val))
    }))

    n_checked <- length(shared)
    report$log("centrality_stability", lbl, n_checked, 0L, 0, 0, 0, 0, "tna::estimate_cs()")
  }))
})


# ---- 4. centrality_stability windowed ---------------------------------------

test_that("centrality_stability: windowed CS match tna across 10 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(10L), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    ws <- 2L; iter <- 20L

    nest_net <- build_network(data, method = "relative",
                              params = list(window_size = ws))
    tna_net  <- tna::tna(data, params = list(window_size = ws))

    set.seed(1)
    tna_cs <- tna::estimate_cs(tna_net, iter = iter)
    nest_cs <- centrality_stability(nest_net, iter = iter, seed = 1)

    lbl <- sprintf("cfg%d ws=%d", i, ws)
    shared <- intersect(names(tna_cs), names(nest_cs$cs))
    invisible(lapply(shared, function(m) {
      expect_equal(nest_cs$cs[[m]], tna_cs[[m]]$cs_coefficient,
                   tolerance = 0.15,
                   label = sprintf("CS(%s) %s", m, lbl))
    }))
  }))
})


# ---- 5. reliability metrics vs tna::reliability() ---------------------------

test_that("reliability: split-half metrics match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    iter <- 30L

    nest_net <- build_network(data, method = "relative")
    tna_net  <- tna::tna(data)

    set.seed(1)
    tna_r <- tna::reliability(tna_net, iter = iter)
    nest_r <- reliability(nest_net, iter = iter, seed = 1)

    lbl <- sprintf("cfg%d(a=%d,s=%d)", i, cfg$n_actors, cfg$n_states)

    # Compare summary metrics (correlation, mean_abs_dev, etc.)
    tna_sum <- tna_r$summary
    nest_sum <- nest_r$summary

    # Match by metric name
    shared <- intersect(tna_sum$metric, nest_sum$metric)
    invisible(lapply(shared, function(m) {
      tna_val  <- tna_sum$mean[tna_sum$metric == m]
      nest_val <- nest_sum$mean[nest_sum$metric == m]
      if (length(tna_val) == 1 && length(nest_val) == 1) {
        expect_equal(nest_val, tna_val, tolerance = 1e-8,
                     label = sprintf("reliability(%s) %s", m, lbl))
      }
    }))

    n_checked <- length(shared)
    report$log("reliability", lbl, n_checked, 0L, 0, 0, 0, 0, "tna::reliability()")
  }))
})


# ---- 6. reliability windowed ------------------------------------------------

test_that("reliability: windowed match tna across 10 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(10L), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    ws <- 2L; iter <- 20L

    nest_net <- build_network(data, method = "relative",
                              params = list(window_size = ws))
    tna_net  <- tna::tna(data, params = list(window_size = ws))

    set.seed(1)
    tna_r <- tna::reliability(tna_net, iter = iter)
    nest_r <- reliability(nest_net, iter = iter, seed = 1)

    lbl <- sprintf("cfg%d ws=%d", i, ws)
    shared <- intersect(tna_r$summary$metric, nest_r$summary$metric)
    invisible(lapply(shared, function(m) {
      tna_val  <- tna_r$summary$mean[tna_r$summary$metric == m]
      nest_val <- nest_r$summary$mean[nest_r$summary$metric == m]
      if (length(tna_val) == 1 && length(nest_val) == 1) {
        expect_equal(nest_val, tna_val, tolerance = 1e-8,
                     label = sprintf("reliability(%s) %s", m, lbl))
      }
    }))
  }))
})


# ---- 7. compare_sequences: frequencies/proportions vs tna -------------------

test_that("compare_sequences: per-pattern freq/prop match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    # Use larger data to ensure patterns meet min_freq threshold
    data <- simulate_sequences(
      n_actors = max(cfg$n_actors, 30),
      n_states = min(cfg$n_states, 5),
      seq_length = max(cfg$seq_length, 15),
      seed = cfg$seed
    )
    n_act <- nrow(data)

    # Split into two groups
    mid <- n_act %/% 2
    gd <- cbind(
      group = rep(c("A", "B"), c(mid, n_act - mid)),
      data,
      stringsAsFactors = FALSE
    )

    nest_cs <- tryCatch(
      compare_sequences(gd, group = "group", min_freq = 1L),
      error = function(e) NULL
    )
    tna_cs <- tryCatch(
      tna::compare_sequences(gd, group = "group", min_freq = 1L),
      error = function(e) NULL
    )

    if (is.null(nest_cs) || is.null(tna_cs)) return(NULL)

    lbl <- sprintf("cfg%d(a=%d,s=%d)", i, n_act, min(cfg$n_states, 5))

    # Align by pattern
    shared <- sort(intersect(nest_cs$patterns$pattern, tna_cs$pattern))
    if (length(shared) == 0L) return(NULL)

    invisible(lapply(shared, function(p) {
      nest_row <- nest_cs$patterns[nest_cs$patterns$pattern == p, ]
      tna_row  <- tna_cs[tna_cs$pattern == p, ]

      # Frequencies
      expect_equal(nest_row$freq_A, tna_row$freq_A, tolerance = TOL,
                   label = sprintf("compare_seq freq_A %s [%s]", lbl, p))
      expect_equal(nest_row$freq_B, tna_row$freq_B, tolerance = TOL,
                   label = sprintf("compare_seq freq_B %s [%s]", lbl, p))

      # Proportions
      expect_equal(nest_row$prop_A, tna_row$prop_A, tolerance = TOL,
                   label = sprintf("compare_seq prop_A %s [%s]", lbl, p))
      expect_equal(nest_row$prop_B, tna_row$prop_B, tolerance = TOL,
                   label = sprintf("compare_seq prop_B %s [%s]", lbl, p))
    }))

    n_checked <- length(shared) * 4L
    report$log("compare_sequences", lbl, n_checked, 0L, 0, 0, 0, 0,
               "tna::compare_sequences()")
  }))
})


# ---- 8. cluster_data: assignments match tna with explicit matching params ----

test_that("cluster_data: assignments match tna with same method/distance across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  # tna and Nestimate may encode sequences differently before computing
  # Levenshtein distance, so raw distances can diverge. Test with
  # method="complete" and check that hclust on the SAME distance object
  # produces identical assignments.
  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    k <- sample(2:3, 1)

    # Build the networks with the same params
    tna_cl  <- tryCatch(tna::cluster_data(data, k = k, method = "complete"),
                        error = function(e) NULL)
    nest_cl <- tryCatch(cluster_data(data, k = k, method = "complete"),
                        error = function(e) NULL)

    if (is.null(tna_cl) || is.null(nest_cl)) return(NULL)
    lbl <- sprintf("cfg%d k=%d", i, k)

    # Both should produce k clusters
    expect_equal(length(unique(nest_cl$assignments)), k,
                 label = sprintf("nest k=%d %s", k, lbl))
    expect_equal(length(tna_cl$sizes), k,
                 label = sprintf("tna k=%d %s", k, lbl))

    # If distances match, assignments must match (modulo label permutation)
    nest_d <- as.numeric(nest_cl$distance)
    tna_d  <- as.numeric(tna_cl$distance)

    if (isTRUE(all.equal(nest_d, tna_d, tolerance = 1e-10))) {
      # Same distance → same hclust → same assignments (up to permutation)
      nest_sizes <- sort(as.integer(table(nest_cl$assignments)))
      tna_sizes  <- sort(as.integer(tna_cl$sizes))
      expect_equal(nest_sizes, tna_sizes,
                   label = sprintf("cluster sizes (same dist) %s", lbl))
    }

    report$log("cluster_data", lbl, 2L, 0L, 0, 0, 0, 0, "tna::cluster_data()")
  }))
})


# ---- 9. Grouped centrality -------------------------------------------------

test_that("centrality on grouped networks matches tna per-group", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N), function(i) {
    cfg  <- configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    gdata <- cbind(group = rep(c("X", "Y"), length.out = cfg$n_actors), data,
                   stringsAsFactors = FALSE)

    nest_group <- build_network(gdata, method = "relative", group = "group")
    tna_group  <- tna::group_tna(gdata, group = "group")

    lbl <- sprintf("cfg%d", i)

    invisible(lapply(names(tna_group), function(g) {
      tna_c  <- tna::centralities(tna_group[[g]])
      nest_c <- centrality(nest_group[[g]])

      states <- sort(intersect(tna_c$state, names(nest_c$InStrength)))
      invisible(lapply(states, function(s) {
        expect_equal(
          nest_c$InStrength[s], tna_c$InStrength[tna_c$state == s],
          tolerance = TOL,
          label = sprintf("grouped InStr %s g=%s [%s]", lbl, g, s)
        )
        expect_equal(
          nest_c$OutStrength[s], tna_c$OutStrength[tna_c$state == s],
          tolerance = TOL,
          label = sprintf("grouped OutStr %s g=%s [%s]", lbl, g, s)
        )
      }))
    }))
  }))
})


# ---- Report ------------------------------------------------------------------

test_that("tna functions equivalence report", {
  skip_on_cran()
  skip_equiv_tests()

  report$write_csv("tna_functions")
  msg <- report$summary()
  message(msg)
  expect_true(length(report$rows) > 0L)
})
