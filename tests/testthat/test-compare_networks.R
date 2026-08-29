# ---- compare_networks(): N-way descriptive comparison + optional inference ----

.cn_fixture <- function() {
  s1 <- simulate_sequences(12, 4, 15, seed = 1)
  s2 <- simulate_sequences(12, 4, 15, seed = 2)
  s3 <- simulate_sequences(12, 4, 15, seed = 3)
  list(
    seqs = list(s1, s2, s3),
    early = build_network(s1, method = "relative"),
    mid = build_network(s2, method = "relative"),
    late = build_network(s3, method = "relative")
  )
}

.cn_edge_cols <- c("pair", "network_a", "network_b", "from", "to", "weight_a",
                   "weight_b", "diff", "abs_diff", "rel_diff", "ratio",
                   "log_ratio", "rank_a", "rank_b", "rank_diff",
                   "percentile_diff", "status", "higher")
.cn_node_cols <- c("pair", "network_a", "network_b", "node", "measure",
                   "value_a", "value_b", "diff", "abs_diff", "rank_a",
                   "rank_b", "higher")
.cn_global_cols <- c("pair", "network_a", "network_b", "category", "metric",
                     "key", "value")

.no_nonfinite <- function(df) {
  num <- vapply(df, is.numeric, logical(1))
  !any(vapply(df[num], function(v) any(is.nan(v) | is.infinite(v)), logical(1)))
}

# 1. ---------------------------------------------------------------------------
test_that("class stack, slots and exact table columns", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  expect_s3_class(cmp, "net_network_comparison")
  expect_named(cmp, c("networks", "matrices", "pairs", "edges", "nodes",
                      "global", "network_metrics", "differences", "scaling",
                      "reference", "measures", "test", "iter", "alpha",
                      "adjust", "paired", "rope", "directed", "n_networks",
                      "n_pairs"))
  expect_named(cmp$pairs, c("pair", "network_a", "network_b"))
  expect_named(cmp$edges, .cn_edge_cols)
  expect_named(cmp$nodes, .cn_node_cols)
  expect_named(cmp$global, .cn_global_cols)
  expect_named(cmp$network_metrics, c("network", "metric", "value"))
  expect_equal(cmp$n_pairs, choose(3, 2))
  expect_equal(nrow(cmp$edges), cmp$n_pairs * 4^2)
  expect_equal(nrow(cmp$nodes), cmp$n_pairs * 4 * 3)
  expect_equal(nrow(cmp$global), cmp$n_pairs * 22)
  expect_equal(nrow(cmp$network_metrics), 3 * 13)
  expect_true(all(vapply(cmp$differences, inherits, logical(1), "netdifference")))
  expect_equal(cmp$test, "none")
  # `higher` names the side with the larger weight (colour channel contract).
  e <- cmp$edges
  tol <- sqrt(.Machine$double.eps)
  expect_equal(as.character(e$higher),
               ifelse(e$diff > tol, e$network_a,
                      ifelse(e$diff < -tol, e$network_b, "equal")))
  # rbind() across pairs unions the factor levels; every pair's levels are
  # its own two names plus "equal".
  expect_true(all(levels(e$higher) %in% c("early", "mid", "late", "equal")))
  nd <- cmp$nodes
  expect_equal(as.character(nd$higher),
               ifelse(nd$diff > tol, nd$network_a,
                      ifelse(nd$diff < -tol, nd$network_b, "equal")))
})

# 2. ---------------------------------------------------------------------------
test_that("all pairs: choose(n, 2) pairs in combn order with 'a vs b' names", {
  f <- .cn_fixture()
  cmp <- compare_networks(a = f$early, b = f$mid, c = f$late, d = f$early)
  expect_equal(nrow(cmp$pairs), 6L)
  expect_equal(cmp$pairs$pair,
               c("a vs b", "a vs c", "a vs d", "b vs c", "b vs d", "c vs d"))
})

# 3. ---------------------------------------------------------------------------
test_that("reference mode by name and by index, reference is always network_a", {
  f <- .cn_fixture()
  by_name <- compare_networks(early = f$early, mid = f$mid, late = f$late,
                              reference = "mid")
  by_index <- compare_networks(early = f$early, mid = f$mid, late = f$late,
                               reference = 2L)
  expect_equal(nrow(by_name$pairs), 2L)
  expect_true(all(by_name$pairs$network_a == "mid"))
  expect_equal(by_name$reference, "mid")
  expect_equal(by_name$edges, by_index$edges)
  expect_error(compare_networks(f$early, f$mid, reference = "zzz"),
               class = "nestimate_compare_reference_unknown")
})

# 4. ---------------------------------------------------------------------------
test_that("swapping a pair flips signs and leaves symmetric metrics unchanged", {
  f <- .cn_fixture()
  ab <- compare_networks(early = f$early, late = f$late)
  ba <- compare_networks(late = f$late, early = f$early)
  e1 <- ab$edges
  e2 <- ba$edges
  expect_equal(e1$diff, -e2$diff)
  expect_equal(e1$rank_diff, -e2$rank_diff)
  expect_equal(e1$percentile_diff, -e2$percentile_diff)
  expect_equal(e1$log_ratio, -e2$log_ratio)
  expect_equal(e1$abs_diff, e2$abs_diff)
  expect_equal(e1$rel_diff, e2$rel_diff)
  finite <- is.finite(e1$ratio) & is.finite(e2$ratio)
  expect_equal(e1$ratio[finite], 1 / e2$ratio[finite])
  expect_equal(as.character(e1$higher), as.character(e2$higher))
  # cv_ratio and rel_mean_abs are asymmetric by definition and excluded.
  sym <- setdiff(unique(ab$global$key), c("cv_ratio", "rel_mean_abs"))
  g1 <- ab$global$value[match(sym, ab$global$key)]
  g2 <- ba$global$value[match(sym, ba$global$key)]
  expect_equal(g1, g2)
})

# 5. ---------------------------------------------------------------------------
test_that("two-network case reproduces compare_model() numerically", {
  f <- .cn_fixture()
  for (sc in c("none", "minmax")) {
    cmp <- compare_networks(x = f$early, y = f$late, scaling = sc)
    # compare_model() is the oracle only; its own loop message and the
    # zero-sd warning from its standardized-score ratio are not under test.
    old <- suppressMessages(suppressWarnings(
      compare_model(f$early, f$late, scaling = sc,
                    measures = c("InStrength", "OutStrength", "Betweenness"))
    ))
    e <- cmp$edges
    o <- old$edge_metrics
    # compare_model() is column-major (source varies fastest); align by key.
    key_new <- paste(e$from, e$to)
    key_old <- paste(o$source, o$target)
    idx <- match(key_new, key_old)
    expect_equal(e$weight_a, o$weight_x[idx], tolerance = 1e-12)
    expect_equal(e$weight_b, o$weight_y[idx], tolerance = 1e-12)
    expect_equal(e$diff, o$raw_difference[idx], tolerance = 1e-12)
    expect_equal(e$abs_diff, o$absolute_difference[idx], tolerance = 1e-12)
    expect_equal(abs(e$rank_diff), o$rank_difference[idx], tolerance = 1e-12)
    expect_equal(abs(e$percentile_diff), o$percentile_difference[idx],
                 tolerance = 1e-12)
    nonneg <- e$weight_a >= 0 & e$weight_b >= 0
    expect_equal(e$log_ratio[nonneg], o$logarithmic_ratio[idx][nonneg],
                 tolerance = 1e-12)
    expect_equal(cmp$global$value, old$summary_metrics$value, tolerance = 1e-12)
    nm <- cmp$network_metrics
    expect_equal(nm$value[nm$network == "x"], old$network_metrics$x,
                 tolerance = 1e-12)
    expect_equal(nm$value[nm$network == "y"], old$network_metrics$y,
                 tolerance = 1e-12)
  }
})

# 6. ---------------------------------------------------------------------------
test_that("ratios are guarded: no Inf/NaN, NA exactly where undefined", {
  A <- matrix(c(0, 2, 0, 0, 1, 3, 0, 0, 0), 3, 3,
              dimnames = list(letters[1:3], letters[1:3]))
  B <- matrix(c(1, 0, 0, 2, 1, 0, 0, 0, 0), 3, 3,
              dimnames = list(letters[1:3], letters[1:3]))
  cmp <- compare_networks(A = A, B = B, measures = NULL)
  e <- cmp$edges
  expect_true(.no_nonfinite(e))
  expect_true(.no_nonfinite(cmp$global))
  expect_equal(is.na(e$ratio), e$weight_b == 0)
  expect_equal(is.na(e$rel_diff), e$weight_a == 0 & e$weight_b == 0)
  expect_equal(as.character(e$status[e$from == "a" & e$to == "b"]), "only_b")
  expect_equal(as.character(e$status[e$from == "b" & e$to == "a"]), "only_a")
  expect_equal(as.character(e$status[e$from == "c" & e$to == "c"]), "neither")
})

# 7. ---------------------------------------------------------------------------
test_that("undirected pair keeps from <= to cells including the diagonal", {
  set.seed(3)
  M1 <- matrix(runif(16), 4, 4); M1 <- (M1 + t(M1)) / 2
  M2 <- matrix(runif(16), 4, 4); M2 <- (M2 + t(M2)) / 2
  dimnames(M1) <- dimnames(M2) <- list(LETTERS[1:4], LETTERS[1:4])
  cmp <- compare_networks(M1, M2, measures = NULL)
  expect_false(any(cmp$directed))
  expect_equal(nrow(cmp$edges), 4 * 5 / 2)
  expect_equal(sum(cmp$edges$from == cmp$edges$to), 4L)
})

# 8. ---------------------------------------------------------------------------
test_that("inputs: matrices, groups, lists, labels", {
  f <- .cn_fixture()
  M <- unname(f$early$weights)
  cmp <- compare_networks(M, unname(f$late$weights), measures = NULL)
  expect_equal(rownames(cmp$matrices[[1]]), paste0("V", 1:4))
  expect_equal(names(cmp$networks), c("network_1", "network_2"))

  long <- rbind(cbind(f$seqs[[1]], grp = "x"), cbind(f$seqs[[2]], grp = "y"))
  grp <- build_network(long, method = "relative", group = "grp")
  cg <- compare_networks(grp)
  expect_equal(cg$pairs$pair, "x vs y")
  expect_equal(cg$matrices$x, grp$x$weights)

  cl <- compare_networks(list(f$early, f$mid, f$late))
  expect_equal(cl$n_networks, 3L)

  lab <- compare_networks(f$early, f$late, labels = c("pre", "post"))
  expect_equal(names(lab$networks), c("pre", "post"))
  expect_equal(lab$pairs$pair, "pre vs post")

  dup <- compare_networks(f$early, f$late)
  expect_equal(names(dup$networks), c("relative", "relative_1"))
})

test_that("tna and group_tna inputs match tna::tna() weights and support inference", {
  skip_if_not_installed("tna")
  f <- .cn_fixture()
  t1 <- tna::tna(f$seqs[[1]])
  t2 <- tna::tna(f$seqs[[2]])
  cmp <- compare_networks(one = t1, two = t2, measures = NULL)
  expect_equal(unname(cmp$matrices$one), unname(t1$weights), tolerance = 1e-12)
  expect_equal(unname(cmp$matrices$two), unname(t2$weights), tolerance = 1e-12)
  expect_true(is.data.frame(cmp$networks$one$data))
  ci <- compare_networks(one = t1, two = t2, measures = NULL,
                         test = "permutation", iter = 20, seed = 1)
  expect_true("perm_p" %in% names(ci$edges))
})

# 9. ---------------------------------------------------------------------------
test_that("classed error conditions", {
  f <- .cn_fixture()
  W <- f$early$weights
  expect_error(compare_networks(f$early), class = "nestimate_compare_too_few")
  expect_error(compare_networks(f$early, "not a network"),
               class = "nestimate_compare_bad_input")
  expect_error(compare_networks(W, W[1:3, 1:3]),
               class = "nestimate_compare_dim_mismatch")
  W2 <- W[c(2, 1, 3, 4), c(2, 1, 3, 4)]
  expect_error(compare_networks(W, W2), class = "nestimate_compare_node_mismatch")
  W3 <- W; W3[1, 1] <- NA
  expect_error(compare_networks(W, W3), class = "nestimate_compare_na_weights")
  expect_error(compare_networks(f$early, f$late, labels = "only_one"),
               class = "nestimate_compare_labels_length")
  Z <- W; Z[1, 2] <- 0
  expect_error(compare_networks(Z, W, scaling = "log"),
               class = "nestimate_compare_scaling_domain")
  expect_error(compare_networks(f$early, f$late, scaling = "minmax",
                                test = "permutation"),
               class = "nestimate_compare_scaling_inference")
  expect_error(compare_networks(W, W, test = "permutation"),
               class = "nestimate_compare_test_unsupported")
  expect_warning(compare_networks(f$early, f$late,
                                  measures = c("InStrength", "Bogus")),
                 class = "nestimate_compare_unknown_measure")
  no_nodes <- compare_networks(f$early, f$late, measures = NULL)
  expect_null(no_nodes$nodes)
  expect_error(node_differences(no_nodes),
               class = "nestimate_compare_no_nodes")
  expect_error(edge_differences(no_nodes, pair = "nope"),
               class = "nestimate_compare_unknown_pair")
  expect_error(plot(no_nodes, what = "nodes"), class = "nestimate_compare_no_nodes")
})

# 10. --------------------------------------------------------------------------
test_that("named table verbs, summary(), digits, print()", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  expect_named(edge_differences(cmp), .cn_edge_cols)
  expect_named(node_differences(cmp), .cn_node_cols)
  expect_named(global_differences(cmp), .cn_global_cols)
  expect_named(network_metrics(cmp), c("network", "metric", "value"))
  one <- edge_differences(cmp, pair = "early vs mid")
  expect_equal(unique(one$pair), "early vs mid")
  expect_equal(nrow(one), 16L)
  expect_identical(attr(one, "row.names"), seq_len(16L))
  nd <- node_differences(cmp, measure = "Betweenness")
  expect_equal(unique(nd$measure), "Betweenness")
  expect_error(node_differences(cmp, measure = "Nope"),
               class = "nestimate_compare_unknown_measure")
  # summary() is a tidy one-row-per-pair data.frame, like every other
  # Nestimate summary() method; the full tables are the four verbs.
  s <- summary(cmp)
  expect_s3_class(s, "net_table")
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 3L)
  expect_named(s, c("pair", "network_a", "network_b", "n_cells", "n_differing",
                    "share_higher_a", "share_higher_b", "mean_abs_diff",
                    "max_abs_diff", "pearson", "spearman", "cosine", "jaccard",
                    "top_edge", "top_edge_higher"))
  expect_equal(summary(cmp, pair = "mid vs late")$pair, "mid vs late")
  expect_true(all(s$share_higher_a + s$share_higher_b <= 1 + 1e-12))
  expect_error(edge_differences("not a comparison"), "net_network_comparison")
  # digits: 2 by default, p-values never rounded, zeros stay 0
  e2 <- edge_differences(cmp)
  expect_true(all(e2$diff == round(e2$diff, 2)))
  e4 <- edge_differences(cmp, digits = 4)
  expect_true(any(e4$diff != round(e4$diff, 2)))
  expect_error(edge_differences(cmp, digits = -1), "digits")
  expect_equal(.cn_fmt(c(0, 0.123456, -0.5), 2), c("0", "0.12", "-0.50"))
  expect_s3_class(e2, "data.frame")
  expect_s3_class(e2, "net_table")
  out <- capture.output(print(e2))
  expect_false(any(grepl("0\\.00", out)))
  expect_false(any(grepl("e[-+][0-9]", out)))
  expect_equal(.cn_fmt(c(0, 0.123456), 2, signed = TRUE), c("0", "+0.12"))
  out <- capture.output(print(cmp))
  expect_true(any(grepl("3 networks, 3 pairs", out)))
  expect_true(any(grepl("edge_differences", out)))
  expect_false(any(grepl("data.frame|\\$", out)))
})

# 11. --------------------------------------------------------------------------
test_that("every plot view builds and carries group names, never x/y", {
  f <- .cn_fixture()
  cmp2 <- compare_networks(early = f$early, late = f$late)
  cmp3 <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  views <- c("edges", "nodes", "global", "heatmap", "scatter")
  for (obj in list(cmp2, cmp3)) {
    for (w in views) {
      p <- plot(obj, what = w, labels = TRUE)
      expect_s3_class(p, "ggplot")
      expect_silent(ggplot2::ggplot_build(p))
    }
  }
  sc <- plot(cmp2, what = "scatter")
  labs <- ggplot2::get_labs(sc)
  expect_equal(labs$x, "early")
  expect_equal(labs$y, "late")
  # top_n caps rows in the edge panel
  b <- ggplot2::ggplot_build(plot(cmp2, what = "edges", top_n = 5))
  expect_equal(nrow(b$data[[1]]), 5L)
  # heatmap limits follow the data, not c(-1, 1)
  counts <- compare_networks(A = f$early$frequency_matrix * 1,
                             B = f$late$frequency_matrix * 1, measures = NULL)
  hb <- ggplot2::ggplot_build(plot(counts, what = "heatmap"))
  lim <- max(abs(counts$edges$diff))
  expect_equal(hb$plot$scales$get_scales("fill")$limits, c(-lim, lim))
  # colour contract: network_a maps to blue
  eb <- ggplot2::ggplot_build(plot(cmp2, what = "edges"))
  seg <- eb$data[[1]]
  higher_a <- cmp2$edges$higher == "early"
  expect_true(all(seg$colour[seg$colour != "grey65"] %in% c("#4A6FE3", "#D33F6A")))
  # measure filter on node panel
  expect_s3_class(plot(cmp2, what = "nodes", measure = "InStrength"), "ggplot")
  expect_error(plot(cmp2, what = "nodes", measure = "Nope"),
               class = "nestimate_compare_unknown_measure")
})

test_that("plot() default draws the network panels in base graphics", {
  f <- .cn_fixture()
  cmp3 <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot(cmp3))
  expect_null(plot(cmp3, what = "networks", pair = "early vs late"))
  # built-in base drawer (used when cograph is absent) draws without error
  expect_null(.cn_draw_network_base(cmp3$matrices$early, "early", FALSE, TRUE))
  expect_null(.cn_draw_network_base(cmp3$matrices$early - cmp3$matrices$late,
                                    "early - late", TRUE, TRUE))
})

test_that("metric tables print in fixed notation", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, late = f$late)
  out <- capture.output(print(network_metrics(cmp)))
  expect_false(any(grepl("e[-+][0-9]", out)))
  out <- capture.output(print(global_differences(cmp)))
  expect_false(any(grepl("e[-+][0-9]", out)))
})

# 12. inference ---------------------------------------------------------------
test_that("permutation backend: columns, equivalence with permutation(), M/S", {
  f <- .cn_fixture()
  ci <- compare_networks(early = f$early, late = f$late, test = "permutation",
                         iter = 60, seed = 10)
  expect_true(all(c("perm_effect", "perm_p", "perm_sig", "sig", "evidence")
                  %in% names(ci$edges)))
  expect_true(all(c("perm_effect", "perm_p", "perm_sig", "sig") %in% names(ci$nodes)))
  expect_true(all(ci$edges$perm_p > 0 & ci$edges$perm_p <= 1))
  expect_equal(ci$edges$sig, ci$edges$perm_p < 0.05)
  expect_equal(unique(ci$edges$evidence), "permutation")
  # pair 1 uses seed + 1
  direct <- permutation(f$early, f$late, iter = 60, seed = 11,
                        measures = c("InStrength", "OutStrength", "Betweenness"))
  expect_equal(ci$edges$perm_p, as.vector(t(direct$p_values)), tolerance = 1e-12)
  expect_equal(ci$edges$perm_effect, as.vector(t(direct$effect_size)),
               tolerance = 1e-12)
  st <- direct$centralities$stats
  key <- paste(st$state, st$centrality)
  expect_equal(ci$nodes$perm_p,
               st$p_value[match(paste(ci$nodes$node, ci$nodes$measure), key)],
               tolerance = 1e-12)
  g <- ci$global
  expect_equal(g$value[g$key == "perm_M"], sum(abs(ci$edges$diff)))
  expect_equal(g$value[g$key == "perm_S"], max(abs(ci$edges$diff)))
  expect_equal(g$perm_p[g$key == "perm_M"], direct$global$p_value[1])
  s <- summary(ci)
  expect_true(all(c("n_sig_edges", "n_sig_nodes", "m_stat", "m_p", "s_stat",
                    "s_p") %in% names(s)))
  expect_equal(s$m_stat, round(sum(abs(ci$edges$diff)), 2))
  expect_equal(summary(ci, digits = 6)$m_p, s$m_p)  # p-values never rounded
})

test_that("Bayesian backend: CI brackets the posterior difference, ROPE decisions", {
  f <- .cn_fixture()
  ci <- compare_networks(early = f$early, late = f$late, test = "bayes",
                         iter = 400, seed = 5, rope = 0.05, measures = NULL)
  cols <- c("bayes_diff", "bayes_ci_lower", "bayes_ci_upper", "bayes_pd",
            "bayes_p", "bayes_sig", "bayes_p_rope", "bayes_decision", "sig",
            "evidence")
  expect_true(all(cols %in% names(ci$edges)))
  e <- ci$edges
  expect_true(all(e$bayes_ci_lower <= e$bayes_diff + 1e-8))
  expect_true(all(e$bayes_ci_upper >= e$bayes_diff - 1e-8))
  expect_true(all(e$bayes_pd >= 0.5 & e$bayes_pd <= 1))
  expect_true(all(e$bayes_p_rope >= 0 & e$bayes_p_rope <= 1))
  expect_true(all(levels(e$bayes_decision) == c("different", "equivalent", "undecided")))
  expect_equal(e$sig, e$bayes_sig)
  expect_equal(unique(e$evidence), "bayes")
  direct <- bayes_compare(f$early, f$late, draws = 400, ci = 0.95, seed = 6)
  expect_equal(e$bayes_ci_lower, as.vector(t(direct$ci_lower)), tolerance = 1e-12)
  expect_equal(e$bayes_pd, as.vector(t(direct$p_difference)), tolerance = 1e-12)
})

test_that("bootstrap backend appends structural rows with CIs", {
  f <- .cn_fixture()
  ci <- compare_networks(early = f$early, late = f$late, test = "bootstrap",
                         iter = 50, seed = 2, measures = NULL)
  g <- ci$global
  boot <- g[g$category == "Structure (bootstrap)", ]
  expect_equal(boot$metric, c("density", "mean_weight", "centralization",
                              "reciprocity"))
  expect_true(all(c("boot_observed_a", "boot_observed_b", "boot_se",
                    "boot_ci_lower", "boot_ci_upper", "boot_z", "boot_p",
                    "boot_sig") %in% names(g)))
  expect_true(all(boot$boot_ci_lower <= boot$value + 1e-8))
  expect_true(all(boot$boot_ci_upper >= boot$value - 1e-8))
  expect_true(all(is.na(g$boot_p[g$category != "Structure (bootstrap)"])))
  direct <- vertex_compare(f$early, f$late, iter = 50, ci_level = 0.05,
                           seed = 3, labels = c("early", "late"))
  expect_equal(boot$value, direct$summary$diff, tolerance = 1e-12)
})

test_that("combined backends fill without collision; plots encode evidence", {
  f <- .cn_fixture()
  ci <- compare_networks(early = f$early, late = f$late,
                         test = c("permutation", "bayes", "bootstrap"),
                         iter = 40, seed = 1, rope = 0.05)
  expect_equal(ci$test, c("permutation", "bayes", "bootstrap"))
  expect_true(all(c("perm_p", "bayes_p", "sig") %in% names(ci$edges)))
  expect_equal(unique(ci$edges$evidence), "permutation")
  g <- ci$global
  expect_equal(sum(g$category == "Global Test (permutation)"), 2L)
  expect_equal(sum(g$category == "Structure (bootstrap)"), 4L)
  for (w in c("edges", "nodes", "global", "heatmap")) {
    p <- plot(ci, what = w)
    expect_silent(b <- ggplot2::ggplot_build(p))
  }
  eb <- ggplot2::ggplot_build(plot(ci, what = "edges"))
  alphas <- unique(eb$data[[1]]$alpha)
  expect_true(all(alphas %in% c(0.35, 1)))
  expect_true(grepl("Evidence: permutation", ggplot2::get_labs(plot(ci, what = "edges"))$caption))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot(ci))
})

test_that("paired permutation errors on unequal sizes", {
  f <- .cn_fixture()
  small <- build_network(f$seqs[[2]][1:6, ], method = "relative")
  expect_error(compare_networks(f$early, small, test = "permutation",
                                iter = 10, paired = TRUE),
               "equal number")
})


# ---- Revision: type =, combined =, per-facet ordering, signed zeros ---------

test_that("plot() takes type = and keeps what = as a working alias", {
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  by_type <- plot(cmp, type = "edges")
  by_what <- plot(cmp, what = "edges")
  expect_s3_class(by_type, "ggplot")
  expect_equal(by_type$data, by_what$data)
  expect_error(plot(cmp, type = "nope"), "should be one of")
})

test_that("each facet of the edge view keeps its own ordering", {
  # Regression: the previous implementation de-duplicated factor levels via
  # `levels<-`, which MERGES levels and forced every facet onto the first
  # facet's row order.
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  d <- plot(cmp, type = "edges")$data
  ordered_within <- vapply(split(d, d$pair), function(dd) {
    !is.unsorted(dd$abs_diff[order(as.integer(dd$label))])
  }, logical(1))
  expect_true(all(ordered_within))
  # The three facets do not share a row order.
  first_rows <- vapply(split(d, d$pair), function(dd) {
    as.character(dd$label[which.max(dd$abs_diff)])
  }, character(1))
  expect_gt(length(unique(first_rows)), 1L)
  # Axis labels are rendered without the pair prefix.
  expect_false(any(grepl("\r", .cn_strip_prefix(levels(d$label)))))
})

test_that("no signed zero is ever printed", {
  # Regression: .cn_fmt() tested for zero before rounding, so -0.0004 printed
  # as "-0.00".
  expect_equal(.cn_fmt(-4e-04, 2L, signed = TRUE), "0")
  expect_equal(.cn_fmt(4e-04, 2L, signed = TRUE), "0")
  expect_equal(.cn_fmt(-0.004, 2L), "0")
  expect_equal(.cn_fmt(-0.015, 2L, signed = TRUE), "-0.01")
  expect_equal(.cn_fmt(-4e-04, 4L, signed = TRUE), "-0.0004")
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  for (tp in c("edges", "heatmap")) {
    b <- ggplot2::ggplot_build(plot(cmp, type = tp))
    labs <- b$data[[length(b$data)]]$label
    expect_false(any(grepl("^[-+]0[.]0+$", labs)))
  }
})

test_that("combined = FALSE returns one plot per pair", {
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  for (tp in c("edges", "nodes", "global", "heatmap", "scatter")) {
    lst <- plot(cmp, type = tp, combined = FALSE)
    expect_type(lst, "list")
    expect_named(lst, cmp$pairs$pair)
    expect_true(all(vapply(lst, inherits, logical(1), "ggplot")))
    expect_true(all(vapply(lst, function(p) length(unique(p$data$pair)) == 1L,
                           logical(1))))
  }
  lst1 <- plot(cmp, type = "edges", pair = "mid vs late", combined = FALSE)
  expect_named(lst1, "mid vs late")
  expect_error(plot(cmp, type = "edges", combined = NA), "combined")
})

test_that("networks view draws each network once; difference is its own view", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  # Three networks, three pairs: the old triptych-per-pair layout drew nine
  # panels for three distinct networks.
  expect_equal(.cn_networks_in(cmp, cmp$pairs$pair), c("early", "mid", "late"))
  expect_equal(.cn_networks_in(cmp, "early vs mid"), c("early", "mid"))
  expect_equal(.cn_grid(3L), c(1L, 3L))
  expect_equal(.cn_grid(4L), c(2L, 2L))
  expect_equal(.cn_grid(6L), c(2L, 3L))
  expect_equal(.cn_grid(1L), c(1L, 1L))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(plot(cmp, type = "networks"))
  expect_null(plot(cmp, type = "difference"))
  expect_null(plot(cmp, type = "networks", combined = FALSE))
  expect_null(plot(cmp, type = "difference", pair = "early vs mid"))
})

test_that("summary() is the tidy pair table and stays under 20 printed lines", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  s <- summary(cmp)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), nrow(cmp$pairs))
  expect_lt(length(capture.output(print(s))), 20L)
})

test_that("type = 'inference' is a forest on the difference scale", {
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  expect_error(plot(cmp, type = "inference"),
               class = "nestimate_compare_no_test")

  perm <- compare_networks(early = f$early, mid = f$mid, late = f$late,
                           test = "permutation", iter = 60, seed = 1)
  p <- plot(perm, type = "inference", top_n = 6)
  expect_s3_class(p, "ggplot")
  # x is the difference, not the weight.
  expect_equal(ggplot2::get_labs(p)$x,
               "Difference (first - second network of pair)")
  expect_true(all(p$data$diff == p$data$weight_a - p$data$weight_b))
  expect_equal(nrow(p$data), 6L * nrow(perm$pairs))
  # Filled marker = significant, open marker = not.
  expect_setequal(unique(p$data$point_fill[!p$data$sig]), "white")
  expect_false(any(p$data$point_fill[p$data$sig] == "white"))
  # p-values are printed, and never rounded to the plot's `digits`.
  b <- ggplot2::ggplot_build(p)
  labs <- b$data[[length(b$data)]]$label
  expect_true(all(grepl("^p = [01][.][0-9]{3}$", labs)))

  bayes <- compare_networks(early = f$early, late = f$late,
                            test = c("permutation", "bayes"), iter = 60,
                            seed = 1)
  pb <- plot(bayes, type = "inference", top_n = 5)
  expect_true(any(vapply(pb$layers, function(l)
    inherits(l$geom, "GeomErrorbar"), logical(1))))
  expect_match(ggplot2::get_labs(pb)$subtitle, "credible interval")
  # A Bayesian figure names its own quantity: bayes_p is 2 * (1 - pd), not a
  # frequentist p-value.
  bayes_only <- compare_networks(early = f$early, late = f$late,
                                 test = "bayes", iter = 60, seed = 1)
  bo <- ggplot2::ggplot_build(plot(bayes_only, type = "inference", top_n = 4))
  expect_true(all(grepl("^p_bayes = ", bo$data[[length(bo$data)]]$label)))
  # No interval layer when only the permutation test ran.
  expect_false(any(vapply(plot(perm, type = "inference")$layers, function(l)
    inherits(l$geom, "GeomErrorbar"), logical(1))))

  lst <- plot(perm, type = "inference", combined = FALSE)
  expect_named(lst, perm$pairs$pair)
  expect_true(all(vapply(lst, inherits, logical(1), "ggplot")))
  expect_silent(ggplot2::ggplot_build(plot(perm, type = "inference")))
})

test_that("pair = accepts network names and indices, not only the composed name", {
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  one <- function(p) unique(edge_differences(cmp, pair = p)$pair)
  expect_equal(one("early vs late"), "early vs late")
  expect_equal(one(c("early", "late")), "early vs late")
  expect_equal(one(c("late", "early")), "early vs late")   # order-free
  expect_equal(one(2L), "early vs late")
  expect_setequal(one("early"), c("early vs mid", "early vs late"))
  expect_setequal(one(c(1L, 3L)), c("early vs mid", "mid vs late"))
  expect_setequal(one(c("early vs mid", "mid vs late")),
                  c("early vs mid", "mid vs late"))
  # Every pair-taking surface uses the same resolver.
  expect_equal(summary(cmp, pair = c("mid", "late"))$pair, "mid vs late")
  expect_equal(unique(node_differences(cmp, pair = c("mid", "late"))$pair),
               "mid vs late")
  expect_equal(unique(global_differences(cmp, pair = 3L)$pair), "mid vs late")
  skip_if_not_installed("ggplot2")
  expect_equal(unique(plot(cmp, type = "edges", pair = c("early", "late"))$data$pair),
               "early vs late")
  expect_error(edge_differences(cmp, pair = "nope"),
               class = "nestimate_compare_unknown_pair")
  expect_error(edge_differences(cmp, pair = 99L),
               class = "nestimate_compare_unknown_pair")
  expect_error(edge_differences(cmp, pair = TRUE),
               class = "nestimate_compare_unknown_pair")
})

test_that("inference and global views encode magnitude as a bar, not a bare dot", {
  skip_if_not_installed("ggplot2")
  f <- .cn_fixture()
  cmp <- compare_networks(early = f$early, mid = f$mid, late = f$late)
  perm <- compare_networks(early = f$early, mid = f$mid, late = f$late,
                           test = "permutation", iter = 60, seed = 1)
  has_seg <- function(p, cls) any(vapply(p$layers, function(l)
    inherits(l$geom, cls), logical(1)))
  # Inference: a bar anchored at zero.
  pi <- plot(perm, type = "inference", top_n = 5)
  expect_true(has_seg(pi, "GeomSegment"))
  expect_match(ggplot2::get_labs(pi)$subtitle, "Bar = difference from zero")
  # The global view stays a point plot against its reference lines: its
  # metrics live on five unrelated scales, so a bar length is not meaningful
  # across panels.
  pg <- plot(cmp, type = "global", pair = c("early", "mid"))
  expect_false(has_seg(pg, "GeomSegment"))
  expect_true(has_seg(pg, "GeomPoint"))
  expect_equal(.cn_global_ref(c("Correlations", "Similarities",
                                "Pattern Similarities", "Weight Deviations",
                                "Dissimilarities")),
               c(1, 1, 1, 0, 0))
  expect_silent(ggplot2::ggplot_build(plot(cmp, type = "global")))
})
