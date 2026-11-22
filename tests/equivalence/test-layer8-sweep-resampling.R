# Layer 8, the argument sweep of the resampling verbs. bootstrap_network() and centrality_stability() draw from
# tna's random stream (test-bootstrap_network.R, test-centrality_stability.R), so under one seed every option is
# compared with tna::bootstrap() and tna::estimate_cs() at the same option, for relative, frequency and
# co-occurrence networks, scaled, weighted and with a start state: a setting dropped when the network is
# re-estimated shows as a difference. casedrop_reliability() draws its subsets with sample() after set.seed(seed)
# (?casedrop_reliability, Details); the reference replays those draws and recomputes every subset's network and
# metrics in base R, so the per-draw arithmetic (target, include_diag, method, scaling, threshold) is checked, not
# the random numbers. permutation() is checked on its deterministic parts: the observed differences and the
# multiplicity adjustment. The defaults are checked written out.

rs_wide <- tn_simulate(n_actors = 40L, n_states = 4L, seq_length = 10L, seed = 821L, ragged = TRUE)
rs_models <- list(relative = "tna", frequency = "ftna", co_occurrence = "ctna")

# tna's bootstrap summary with its bootstrap means and SDs, one row per edge Nestimate reports (an undirected
# network lists each pair once, from the earlier state in the matrix order).
rs_tna_boot <- function(boot, directed = TRUE) {
  summary <- sw_tna_bootstrap(boot)
  if (!directed) {
    states <- rownames(boot$weights_orig)
    summary <- summary[match(summary$from, states) <= match(summary$to, states), , drop = FALSE]
  }
  summary
}

# A tna model of the same kind and settings as a Nestimate network.
rs_tna_model <- function(wide, method, scaling = NULL, params = list(), start = NULL) {
  args <- c(list(wide), if (!is.null(scaling)) list(scaling = scaling), list(params = params),
            if (!is.null(start)) list(begin_state = start))
  do.call(getExportedValue("tna", rs_models[[method]]), args)
}

test_that("every inference, range, edge threshold and level of bootstrap_network() equals tna::bootstrap()", {
  skip_without_reference("tna")
  settings <- list(
    list(inference = "stability", consistency_range = c(0.75, 1.25), edge_threshold = NULL, ci_level = 0.05),
    list(inference = "stability", consistency_range = c(0.5, 1.5), edge_threshold = NULL, ci_level = 0.05),
    list(inference = "stability", consistency_range = c(0.9, 1.1), edge_threshold = NULL, ci_level = 0.2),
    list(inference = "threshold", consistency_range = c(0.75, 1.25), edge_threshold = NULL, ci_level = 0.05),
    list(inference = "threshold", consistency_range = c(0.75, 1.25), edge_threshold = 0.15, ci_level = 0.1)
  )
  cases <- expand.grid(method = names(rs_models), setting = seq_along(settings), stringsAsFactors = FALSE)
  invisible(Map(\(method, setting_index) {
    s <- settings[[setting_index]]
    model <- rs_tna_model(rs_wide, method)
    # A frequency threshold is on the scale of the counts.
    edge_threshold <- if (!is.null(s$edge_threshold) && !identical(method, "relative")) s$edge_threshold * 20 else s$edge_threshold
    tna_args <- c(list(model, iter = 60L, level = s$ci_level, method = s$inference,
                       consistency_range = s$consistency_range),
                  if (!is.null(edge_threshold)) list(threshold = edge_threshold))
    set.seed(7L)
    reference <- do.call(tna::bootstrap, tna_args)
    found <- bootstrap_network(build_network(rs_wide, method = method), iter = 60L, ci_level = s$ci_level,
                               inference = s$inference, consistency_range = s$consistency_range,
                               edge_threshold = edge_threshold, seed = 7L)
    expected <- rs_tna_boot(reference, directed = !identical(method, "co_occurrence"))
    if (identical(s$inference, "threshold")) expected <- expected[setdiff(names(expected), c("cr_lower", "cr_upper"))]
    expect_reference(summary(found), expected, reference = "tna",
                     covers = c("bootstrap_network()", "summary.net_bootstrap()"),
                     quantity = sprintf("%s, %s, range %s, edge threshold %s, level %s", method, s$inference,
                                        paste(s$consistency_range, collapse = "-"), edge_threshold %||% "default",
                                        s$ci_level),
                     tolerance = 1e-10, key = c("from", "to"))
  }, cases$method, cases$setting))
})

test_that("a scaled, weighted or started network is bootstrapped with its own settings, as tna::bootstrap()", {
  skip_without_reference("tna")
  variants <- list(
    list(method = "relative", scaling = "minmax", params = list(), start = NULL, label = "relative, minmax"),
    list(method = "frequency", scaling = "max", params = list(), start = NULL, label = "frequency, max"),
    list(method = "frequency", scaling = NULL, params = list(weighted = TRUE), start = NULL, label = "frequency, weighted"),
    list(method = "co_occurrence", scaling = NULL, params = list(weighted = TRUE), start = NULL,
         label = "co-occurrence, weighted"),
    list(method = "relative", scaling = NULL, params = list(), start = "Start", label = "relative, start state"),
    # rank scaling is left out: tna ranks the zero cells too (test-build_network.R, a known difference).
    list(method = "frequency", scaling = "minmax", params = list(), start = "Start", label = "frequency, minmax, start state")
  )
  invisible(lapply(variants, \(v) {
    set.seed(11L)
    reference <- tna::bootstrap(rs_tna_model(rs_wide, v$method, v$scaling, v$params, v$start), iter = 60L)
    net <- build_network(rs_wide, method = v$method, scaling = v$scaling, params = v$params,
                         start = if (is.null(v$start)) FALSE else TRUE)
    found <- bootstrap_network(net, iter = 60L, seed = 11L)
    expect_reference(summary(found), rs_tna_boot(reference, directed = !identical(v$method, "co_occurrence")),
                     reference = "tna", covers = c("bootstrap_network()", "summary.net_bootstrap()"),
                     quantity = sprintf("bootstrap of a network built with %s", v$label), tolerance = 1e-10,
                     key = c("from", "to"))
  }))
})

test_that("basic intervals, strict boundaries and the pruned network follow their definitions", {
  skip_without_reference("tna")
  invisible(lapply(names(rs_models), \(method) {
    set.seed(5L)
    reference <- tna::bootstrap(rs_tna_model(rs_wide, method), iter = 60L)
    expected <- rs_tna_boot(reference, directed = !identical(method, "co_occurrence"))
    basic <- bootstrap_network(build_network(rs_wide, method = method), iter = 60L, seed = 5L, ci_method = "basic")
    # Davison & Hinkley (1997) eq. 5.6: the percentile bounds reflected around the observed weight.
    expect_reference(summary(basic), data.frame(from = expected$from, to = expected$to,
                                                ci_lower = 2 * expected$weight - expected$ci_upper,
                                                ci_upper = 2 * expected$weight - expected$ci_lower),
                     reference = "independent: 2 w - q(1 - a/2), 2 w - q(a/2) from tna's percentile bounds",
                     covers = "bootstrap_network()", quantity = sprintf("basic interval, %s", method),
                     tolerance = 1e-10, key = c("from", "to"))
    inclusive <- bootstrap_network(build_network(rs_wide, method = method), iter = 60L, seed = 5L)
    strict <- bootstrap_network(build_network(rs_wide, method = method), iter = 60L, seed = 5L, boundary = "strict")
    # Strict counts a subset of the replicates inclusive counts: its p-value is never larger.
    expect_reference(as.numeric(sum(summary(strict)$p_value > summary(inclusive)$p_value + 1e-12)), 0,
                     reference = "independent: strict outside-the-range counts are a subset of inclusive ones",
                     covers = "bootstrap_network()", quantity = sprintf("strict p-values at most inclusive, %s", method),
                     tolerance = 0)
    sig <- summary(inclusive)
    expect_reference(data.frame(from = sig$from, to = sig$to, sig = sig$sig),
                     data.frame(from = sig$from, to = sig$to, sig = sig$p_value < 0.05),
                     reference = "independent: ?bootstrap_network, significant when the p-value is below ci_level",
                     covers = "bootstrap_network()", quantity = sprintf("significance flags, %s", method),
                     tolerance = 0, key = c("from", "to"))
    expected_pruned <- inclusive$original$weights * (inclusive$p_values < 0.05)
    expect_reference(extract_transition_matrix(inclusive$model), expected_pruned,
                     reference = "independent: ?bootstrap_network, original weights where p < ci_level, 0 elsewhere",
                     covers = "bootstrap_network()", quantity = sprintf("pruned network, %s", method), tolerance = 0)
  }))
})

test_that("every default of bootstrap_network() equals the value written out", {
  net <- build_network(rs_wide, method = "relative")
  omitted <- bootstrap_network(net, iter = 40L, seed = 3L)
  written <- bootstrap_network(net, iter = 40L, ci_level = 0.05, inference = "stability",
                               consistency_range = c(0.75, 1.25), edge_threshold = NULL, seed = 3L,
                               boundary = "inclusive", ci_method = "percentile", actor = NULL)
  expect_reference(summary(written), summary(omitted), reference = "independent: the documented defaults written out",
                   covers = "bootstrap_network()", quantity = "bootstrap defaults", tolerance = 0, key = c("from", "to"))
  threshold <- bootstrap_network(build_network(rs_wide, method = "frequency"), iter = 40L, seed = 3L,
                                 inference = "threshold")
  weights <- extract_transition_matrix(build_network(rs_wide, method = "frequency"))
  expect_reference(threshold$edge_threshold, unname(stats::quantile(abs(weights[weights != 0]), 0.1)),
                   reference = "independent: ?bootstrap_network, the 10th percentile of the non-zero absolute weights",
                   covers = "bootstrap_network()", quantity = "default edge threshold", tolerance = 1e-12)
})

# centrality_stability() against tna::estimate_cs() under one seed. Two divergences are kept out of the sweep
# and declared on their own below: rank correlations (spearman, kendall) of centralities whose ties tna computes
# exactly and Nestimate breaks by rounding error, and draws where a node's centrality is undefined (drop 0.9:
# unreachable nodes, clustering with no ties). Closeness with invert = FALSE is checked once, on its own.
cs_compare <- function(method, measures, s, o, label, defect = NULL) {
  set.seed(2L)
  reference <- tna::estimate_cs(rs_tna_model(rs_wide, method), loops = o$loops, normalize = o$normalize,
                                invert = o$invert, measures = measures, iter = 20L, method = s$method,
                                drop_prop = s$drop_prop, threshold = s$threshold, certainty = s$certainty)
  found <- withCallingHandlers(
    centrality_stability(build_network(rs_wide, method = method), measures = measures, iter = 20L,
                         drop_prop = s$drop_prop, threshold = s$threshold, certainty = s$certainty,
                         method = s$method, loops = o$loops, normalize = o$normalize, invert = o$invert, seed = 2L),
    nestimate_constant_measure = \(w) invokeRestart("muffleWarning"))
  list(found = found, reference = reference)
}

test_that("every measure, drop grid, threshold, certainty and correlation of centrality_stability() equals tna", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  linear_sets <- list(c("InStrength", "OutStrength", "Betweenness"), c("Closeness", "ClosenessIn", "ClosenessOut"),
                      c("Diffusion", "Clustering", "BetweennessRSP"))
  # Measures without rounding-error ties in these networks: integer betweenness, in-strength, diffusion, RSP.
  rank_set <- c("InStrength", "Betweenness", "Diffusion", "BetweennessRSP")
  settings <- list(
    list(drop_prop = seq(0.1, 0.8, by = 0.1), threshold = 0.7, certainty = 0.95, method = "pearson"),
    list(drop_prop = c(0.25, 0.5, 0.75), threshold = 0.5, certainty = 0.8, method = "pearson"),
    list(drop_prop = c(0.25, 0.5, 0.75), threshold = 0.5, certainty = 0.8, method = "spearman"),
    list(drop_prop = c(0.2, 0.4, 0.6), threshold = 0.9, certainty = 0.99, method = "kendall")
  )
  options <- list(list(loops = FALSE, normalize = FALSE, invert = TRUE), list(loops = TRUE, normalize = TRUE, invert = TRUE),
                  list(loops = FALSE, normalize = FALSE, invert = FALSE))
  cases <- expand.grid(set = seq_along(linear_sets), setting = seq_along(settings), option = seq_along(options),
                       method = c("relative", "frequency"), stringsAsFactors = FALSE)
  # The rank methods use one measure set: keep one copy of their cases.
  rank_case <- vapply(cases$setting, \(i) !identical(settings[[i]]$method, "pearson"), logical(1))
  cases <- cases[!rank_case | cases$set == 1L, , drop = FALSE]
  invisible(Map(\(set_index, setting_index, option_index, method) {
    s <- settings[[setting_index]]
    o <- options[[option_index]]
    measures <- if (identical(s$method, "pearson")) linear_sets[[set_index]] else rank_set
    if (!o$invert) measures <- setdiff(measures, "Closeness")
    pair <- cs_compare(method, measures, s, o)
    label <- sprintf("%s, %s, drop %s, threshold %s, certainty %s, %s, loops %s, normalize %s, invert %s", method,
                     paste(measures, collapse = "+"), paste(s$drop_prop, collapse = "/"), s$threshold, s$certainty,
                     s$method, o$loops, o$normalize, o$invert)
    # The CS coefficients are a field of the result (API gap, test-centrality_stability.R).
    expect_reference(pair$found$cs, vapply(pair$reference, \(m) m$cs_coefficient, numeric(1)), reference = "tna",
                     covers = "centrality_stability()", quantity = sprintf("CS coefficients, %s", label),
                     tolerance = 1e-10)
    expect_reference(summary(pair$found), sw_stability_summary(pair$reference, s$threshold), reference = "tna",
                     covers = c("centrality_stability()", "summary.net_stability()"),
                     quantity = sprintf("per drop proportion, %s", label), tolerance = 1e-10,
                     key = c("measure", "drop_prop"))
  }, cases$set, cases$setting, cases$option, cases$method))
})

test_that("centrality_stability() differs from tna where ties or undefined nodes enter, and agrees on invert = FALSE closeness", {
  skip_without_reference("tna")
  sw_skip_old_tna()
  base <- list(loops = FALSE, normalize = FALSE, invert = TRUE)
  draws <- \(pair, m) pair$found$correlations[[m]]
  tied <- cs_compare("relative", c("InStrength", "OutStrength"),
                     list(drop_prop = c(0.25, 0.5, 0.75), threshold = 0.5, certainty = 0.8, method = "spearman"), base)
  expect_known_difference(draws(tied, "OutStrength"), tied$reference$OutStrength$correlations, reference = "tna",
                          covers = "centrality_stability()", quantity = "spearman correlations per draw, OutStrength",
                          tolerance = 1e-10,
                          reason = "Out-strengths that are equal in exact arithmetic (1 - P_ii) differ in the last bits in Nestimate, so rank() breaks a tie tna keeps (one draw: -0.80 against tna's -0.74).")
  undefined <- cs_compare("relative", c("Closeness", "ClosenessIn", "ClosenessOut"),
                          list(drop_prop = c(0.5, 0.9), threshold = 0.7, certainty = 0.95, method = "pearson"), base)
  expect_known_difference(summary(undefined$found), sw_stability_summary(undefined$reference, 0.7), reference = "tna",
                          covers = "centrality_stability()", quantity = "closeness at drop 0.9, per drop proportion",
                          tolerance = 1e-10, key = c("measure", "drop_prop"),
                          reason = "At drop 0.9 some nodes are unreachable in a draw; Nestimate correlates the nodes whose centrality is finite in both networks, tna's cor(use = \"complete.obs\") after an sd() check that is NA when any node is undefined.")
  inverted <- cs_compare("relative", c("Closeness", "InStrength"),
                         list(drop_prop = c(0.3, 0.6), threshold = 0.7, certainty = 0.95, method = "pearson"),
                         list(loops = FALSE, normalize = FALSE, invert = FALSE))
  expect_reference(draws(inverted, "Closeness"), inverted$reference$Closeness$correlations, reference = "tna",
                   covers = "centrality_stability()", quantity = "Closeness correlations per draw, invert = FALSE",
                   tolerance = 1e-10)
})
# casedrop_reliability()'s documented procedure, replayed with its draws: set.seed(seed), then for each drop
# proportion in order and each iteration, sample(n, n - floor(n p)) sequences; every subset's network is
# recounted from the sequences here and compared with the original over the cells `target` and `include_diag`
# select.
cd_replay <- function(wide, method, seed, iter, drop_prop, target, include_diag, cor_method, scaling = NULL,
                      threshold = 0) {
  mat <- sw_positions(wide)
  states <- sw_states(mat)
  kind <- if (identical(method, "co_occurrence")) "co_occurrence" else "frequency"
  per_sequence <- lapply(seq_len(nrow(mat)), \(i) sw_reference(mat[i, , drop = FALSE], states, kind))
  finish <- function(counts) {
    m <- if (identical(method, "relative")) sw_row_normalise(counts) else counts
    if (!is.null(scaling)) m <- sw_scale(m, scaling)
    sw_threshold(m, threshold)
  }
  original <- finish(Reduce(`+`, per_sequence))
  target <- if (identical(target, "auto")) if (is.null(scaling) && threshold == 0) "existing" else "all" else target
  cells <- (include_diag | row(original) != col(original)) & (identical(target, "all") | original != 0)
  n <- length(per_sequence)
  set.seed(seed)
  draws <- lapply(drop_prop, \(p) {
    keep <- n - floor(n * p)
    vapply(seq_len(iter), \(i) {
      sub <- finish(Reduce(`+`, per_sequence[sample(seq_len(n), keep)]))
      change <- sub[cells] - original[cells]
      c(mean_abs_dev = mean(abs(change)), median_abs_dev = stats::median(abs(change)),
        # A subset whose compared cells are all equal has no correlation (?casedrop_reliability: undefined).
        correlation = if (stats::sd(sub[cells]) > 0) stats::cor(original[cells], sub[cells], method = cor_method) else NA_real_,
        max_abs_dev = max(abs(change)))
    }, numeric(4))
  })
  metrics <- lapply(c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev"), \(m) {
    out <- vapply(draws, \(d) d[m, ], numeric(iter))
    matrix(out, iter, length(drop_prop), dimnames = list(NULL, paste0("p", drop_prop)))
  })
  names(metrics) <- c("mean_abs_dev", "median_abs_dev", "correlation", "max_abs_dev")
  list(metrics = metrics, n_edges = sum(cells))
}

test_that("every target, include_diag, correlation, threshold and certainty of casedrop_reliability() equals the replayed procedure", {
  configs <- list(list(method = "relative", scaling = NULL, threshold = 0),
                  list(method = "frequency", scaling = NULL, threshold = 0),
                  list(method = "co_occurrence", scaling = NULL, threshold = 0),
                  list(method = "relative", scaling = "max", threshold = 0),
                  list(method = "frequency", scaling = NULL, threshold = 3))
  cases <- expand.grid(config = seq_along(configs), target = c("auto", "existing", "all"),
                       include_diag = c(FALSE, TRUE), cor_method = c("spearman", "pearson", "kendall"),
                       stringsAsFactors = FALSE)
  drop_prop <- c(0.2, 0.5, 0.8)
  invisible(Map(\(config_index, target, include_diag, cor_method) {
    config <- configs[[config_index]]
    net <- build_network(rs_wide, method = config$method, scaling = config$scaling, threshold = config$threshold)
    found <- casedrop_reliability(net, iter = 15L, drop_prop = drop_prop, method = cor_method,
                                  include_diag = include_diag, target = target, seed = 13L)
    reference <- cd_replay(rs_wide, config$method, seed = 13L, iter = 15L, drop_prop = drop_prop, target = target,
                           include_diag = include_diag, cor_method = cor_method, scaling = config$scaling,
                           threshold = config$threshold)
    label <- sprintf("%s, scaling %s, threshold %s, target %s, include_diag %s, %s", config$method,
                     config$scaling %||% "none", config$threshold, target, include_diag, cor_method)
    lapply(names(reference$metrics), \(m) {
      # The per-iteration metrics are a field of the result (API gap: summary() reports their moments only).
      expect_reference(found$metrics[[m]], reference$metrics[[m]],
                       reference = "independent: case-dropping replayed with the same draws, base R",
                       covers = "casedrop_reliability()", quantity = sprintf("%s per draw, %s", m, label),
                       tolerance = 1e-10)
    })
    expect_reference(found$n_edges, reference$n_edges,
                     reference = "independent: cells selected by target and include_diag, base R",
                     covers = "casedrop_reliability()", quantity = sprintf("cells compared, %s", label), tolerance = 0)
  }, cases$config, cases$target, cases$include_diag, cases$cor_method))
  # threshold and certainty only decide the CS coefficient from the same correlations.
  grid <- expand.grid(threshold = c(0.3, 0.7, 0.95), certainty = c(0.5, 0.95), stringsAsFactors = FALSE)
  reference <- cd_replay(rs_wide, "relative", seed = 17L, iter = 25L, drop_prop = c(0.1, 0.3, 0.5, 0.7),
                         target = "auto", include_diag = FALSE, cor_method = "spearman")
  invisible(Map(\(threshold, certainty) {
    found <- casedrop_reliability(build_network(rs_wide, method = "relative"), iter = 25L,
                                  drop_prop = c(0.1, 0.3, 0.5, 0.7), threshold = threshold, certainty = certainty,
                                  seed = 17L)
    share <- colMeans(!is.na(reference$metrics$correlation) & reference$metrics$correlation >= threshold)
    passing <- c(0.1, 0.3, 0.5, 0.7)[share >= certainty]
    expect_reference(found$cs, if (length(passing)) max(passing) else 0,
                     reference = "independent: largest drop proportion whose share at or above threshold reaches certainty",
                     covers = "casedrop_reliability()",
                     quantity = sprintf("CS coefficient, threshold %s, certainty %s", threshold, certainty), tolerance = 0)
  }, grid$threshold, grid$certainty))
})

test_that("every default of casedrop_reliability() and centrality_stability() equals the value written out", {
  net <- build_network(rs_wide, method = "relative")
  omitted <- casedrop_reliability(net, iter = 10L, seed = 1L)
  written <- casedrop_reliability(net, iter = 10L, drop_prop = seq(0.1, 0.9, by = 0.1), threshold = 0.7,
                                  certainty = 0.95, method = "spearman", include_diag = FALSE, seed = 1L,
                                  target = "auto", actor = NULL)
  expect_reference(summary(written), summary(omitted), reference = "independent: the documented defaults written out",
                   covers = "casedrop_reliability()", quantity = "case-dropping defaults", tolerance = 0,
                   key = c("metric", "drop_prop"))
  omitted <- centrality_stability(net, iter = 10L, seed = 1L)
  written <- centrality_stability(net, measures = c("InStrength", "OutStrength", "Betweenness"), iter = 10L,
                                  drop_prop = seq(0.1, 0.9, by = 0.1), threshold = 0.7, certainty = 0.95,
                                  method = "pearson", centrality_fn = NULL, loops = FALSE, normalize = FALSE,
                                  invert = TRUE, normalize_diffusion = TRUE, seed = 1L, actor = NULL)
  expect_reference(summary(written), summary(omitted), reference = "independent: the documented defaults written out",
                   covers = "centrality_stability()", quantity = "centrality-stability defaults", tolerance = 0,
                   key = c("measure", "drop_prop"))
})

test_that("observed differences of permutation() equal the difference of the two networks, every method and setting", {
  x <- tn_simulate(n_actors = 25L, n_states = 4L, seq_length = 9L, seed = 822L, ragged = TRUE)
  y <- tn_simulate(n_actors = 25L, n_states = 4L, seq_length = 9L, seed = 823L, ragged = TRUE)
  configs <- list(
    list(method = "relative", params = list(), scaling = NULL, start = FALSE),
    list(method = "frequency", params = list(), scaling = NULL, start = FALSE),
    list(method = "co_occurrence", params = list(), scaling = NULL, start = FALSE),
    list(method = "relative", params = list(), scaling = "max", start = TRUE),
    list(method = "frequency", params = list(weighted = TRUE), scaling = NULL, start = FALSE),
    list(method = "attention", params = list(lambda = 2), scaling = NULL, start = FALSE),
    list(method = "gap", params = list(max_gap = 2L), scaling = NULL, start = FALSE)
  )
  invisible(lapply(configs, \(config) {
    fit <- \(data) build_network(data, method = config$method, params = config$params, scaling = config$scaling,
                                 start = config$start)
    reference <- \(data) {
      start <- if (isTRUE(config$start)) "Start" else NULL
      mat <- sw_positions(data, start)
      weights <- sw_reference(mat, sw_states(sw_positions(rbind(x, y), start), start), config$method, config$params)
      if (is.null(config$scaling)) weights else sw_scale(weights, config$scaling)
    }
    found <- permutation(fit(x), fit(y), iter = 20L, seed = 1L)
    diff <- reference(x) - reference(y)
    table <- rs_edge_table(weight_x = reference(x), weight_y = reference(y), diff = diff)
    table <- table[table$weight_x != 0 | table$weight_y != 0, , drop = FALSE]
    if (identical(config$method, "co_occurrence")) {
      states <- rownames(diff)
      table <- table[match(table$from, states) <= match(table$to, states), , drop = FALSE]
    }
    expect_reference(summary(found), table, reference = "independent: weights of each sequence set, base R",
                     covers = c("permutation()", "summary.net_permutation()"),
                     quantity = sprintf("observed differences, %s, %s, scaling %s, start %s", config$method,
                                        paste(names(config$params), unlist(config$params), collapse = " "),
                                        config$scaling %||% "none", config$start),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("every adjust method of permutation() equals stats::p.adjust() of the raw p-values, paired and not", {
  x <- tn_simulate(n_actors = 20L, n_states = 3L, seq_length = 8L, seed = 824L)
  y <- tn_simulate(n_actors = 20L, n_states = 3L, seq_length = 8L, seed = 825L)
  cases <- expand.grid(method = c("relative", "frequency", "co_occurrence"), paired = c(FALSE, TRUE),
                       adjust = stats::p.adjust.methods, stringsAsFactors = FALSE)
  invisible(Map(\(method, paired, adjust) {
    net_x <- build_network(x, method = method)
    net_y <- build_network(y, method = method)
    raw <- permutation(net_x, net_y, iter = 30L, paired = paired, adjust = "none", seed = 9L)
    adjusted <- permutation(net_x, net_y, iter = 30L, paired = paired, adjust = adjust, seed = 9L)
    # The p-value matrices are fields of the result (API gap, test-permutation.R); the adjustment runs over every
    # cell the test reports.
    cells <- if (identical(method, "co_occurrence")) upper.tri(raw$p_values, diag = TRUE) else
      matrix(TRUE, nrow(raw$p_values), ncol(raw$p_values))
    expected <- raw$p_values
    expected[cells] <- stats::p.adjust(raw$p_values[cells], method = adjust)
    expect_reference(adjusted$p_values[cells], expected[cells], reference = "stats", covers = "permutation()",
                     quantity = sprintf("%s, paired %s, adjust %s", method, paired, adjust), tolerance = 1e-12)
  }, cases$method, cases$paired, cases$adjust))
})
