# ==============================================================================
# compare_networks(): descriptive N-way network comparison with optional
# inference. One tidy object, one colour contract (see compare_networks_plot.R).
# ==============================================================================

#' Compare two or more networks
#'
#' Compares any number of networks pairwise (all pairs, or every network
#' against one `reference`) and returns one tidy object: an edge table, a
#' node (centrality) table, a global-metric table and per-network structural
#' metrics, each returned by a named verb and carrying a `pair` column so
#' nothing ever needs list indexing.
#' Descriptive by default; `test =` adds permutation, Bayesian and bootstrap
#' evidence to the same tables. `plot()` draws one view per call.
#'
#' @param ... Two or more networks, in any mix of: `netobject`,
#'   `netobject_group` (members are flattened and keep their names),
#'   `cograph_network` / `psychnet`, `mcml`, `tna`, `group_tna`, square numeric
#'   matrices, or one unnamed `list` of these. Name the arguments to name the
#'   networks (`compare_networks(early = a, late = b)`).
#' @param reference `NULL` (default) compares all pairs. A single network name
#'   or index compares every other network against that one; the reference is
#'   always `network_a`, so `diff = reference - other`.
#' @param scaling Scaling applied to every network before comparison; one of
#'   `"none"` (default), `"minmax"`, `"max"`, `"rank"`, `"zscore"`, `"robust"`,
#'   `"log"`, `"log1p"`, `"softmax"`, `"quantile"`, `"frobenius"`, `"row"`.
#'   Inference (`test != "none"`) requires `"none"`: the tests are defined on
#'   the networks as estimated.
#' @param measures Centrality measures for the node table. Any of
#'   `OutStrength`, `InStrength`, `ClosenessIn`, `ClosenessOut`, `Closeness`,
#'   `Betweenness`, `BetweennessRSP`, `Diffusion`, `Clustering`. `NULL` or
#'   `character(0)` skips the node table.
#' @param labels Optional character vector naming the networks (one per
#'   network after flattening groups); overrides argument names.
#' @param test Character vector of inference backends, any of `"none"`
#'   (default), `"permutation"`, `"bayes"`, `"bootstrap"`. Several may be
#'   combined; each fills the columns it supports (see Details).
#' @param iter Number of permutations / posterior draws / bootstrap
#'   replicates per pair. Default `1000`.
#' @param alpha Significance level (permutation, bootstrap) and `1 - ci`
#'   (Bayesian credible interval). Default `0.05`.
#' @param adjust Multiplicity adjustment for permutation p-values, passed to
#'   [stats::p.adjust()] within each pair and table. Default `"none"`.
#' @param paired Logical; paired permutation (equal observation counts).
#' @param rope Optional half-width of a region of practical equivalence on the
#'   difference scale (Bayesian backend only). Adds `bayes_p_rope` and a
#'   `bayes_decision` of `"different"`, `"equivalent"` or `"undecided"`.
#' @param seed Optional integer seed. Each pair uses `seed + pair index`, so
#'   results are reproducible and independent of pair order.
#'
#' @details
#' **Guarding.** Ratios are never `Inf`/`NaN`: `ratio` is `NA` when
#' `weight_b == 0`, `rel_diff` is `NA` when both weights are 0, and
#' `log_ratio = log1p(a) - log1p(b)` is `NA` when either weight is negative.
#' No pseudo-counts are added.
#'
#' **Cells.** Every cell of the weight matrix is a row, including the
#' diagonal and edges absent from one network (weight 0); when both networks
#' are undirected only `from <= to` cells are kept.
#'
#' **Inference.** `"permutation"` (via [permutation()]) adds
#' `perm_effect`, `perm_p`, `perm_sig` to `edges` and `nodes`, and two rows
#' `M` (sum of absolute edge differences) and `S` (largest absolute edge
#' difference) to `global` with permutation p-values. `"bayes"` (via
#' [bayes_compare()]) adds `bayes_diff` (posterior mean difference),
#' `bayes_ci_lower`, `bayes_ci_upper`, `bayes_pd` (probability of direction),
#' `bayes_p`, `bayes_sig` to `edges`; with `rope`, `bayes_p_rope` (normal
#' approximation from the posterior mean and SD) and `bayes_decision`.
#' `"bootstrap"` (via [vertex_compare()]) appends structural rows
#' (density, mean weight, centralization, reciprocity) to `global` with
#' `boot_se`, `boot_ci_lower`, `boot_ci_upper`, `boot_z`, `boot_p`,
#' `boot_sig`. Unified `sig` and `evidence` columns on `edges`/`nodes` take
#' the permutation result when run, else the Bayesian one.
#' Permutation and Bayesian tests need networks that carry their data
#' (`build_network()` output, or `tna` objects, which are rebuilt); plain
#' matrices support `"bootstrap"` only.
#'
#' @return An object of class `net_network_comparison`: a list with
#' * `networks`: named list of the input networks as `netobject`s (scaled
#'   weights);
#' * `matrices`: named list of scaled weight matrices;
#' * `pairs`: data.frame, one row per comparison: `pair`, `network_a`,
#'   `network_b`;
#' * `edges`: data.frame, one row per pair x cell: `pair`, `network_a`,
#'   `network_b`, `from`, `to`, `weight_a`, `weight_b`, `diff`, `abs_diff`,
#'   `rel_diff`, `ratio`, `log_ratio`, `rank_a`, `rank_b`, `rank_diff`,
#'   `percentile_diff`, `status` (`both`/`only_a`/`only_b`/`neither`),
#'   `higher` (`network_a`, `network_b` or `equal`), plus inference columns;
#' * `nodes`: data.frame (or `NULL`), one row per pair x node x measure:
#'   `pair`, `network_a`, `network_b`, `node`, `measure`, `value_a`,
#'   `value_b`, `diff`, `abs_diff`, `rank_a`, `rank_b`, `higher`, plus
#'   inference columns;
#' * `global`: data.frame, one row per pair x metric (22 descriptive metrics
#'   in five categories, plus inference rows): `pair`, `network_a`,
#'   `network_b`, `category`, `metric`, `key`, `value`, plus inference
#'   columns;
#' * `network_metrics`: data.frame, one row per network x structural metric:
#'   `network`, `metric`, `value`;
#' * `differences`: named list of `netdifference` objects (one per pair);
#' * `scaling`, `reference`, `measures`, `test`, `iter`, `alpha`, `adjust`,
#'   `paired`, `rope`, `directed` (named logical), `n_networks`, `n_pairs`.
#'
#' `summary()` returns the one-row-per-pair overview table; the full tables
#' come from the named verbs [edge_differences()], [node_differences()],
#' [global_differences()] and [network_metrics()]. `plot()` draws one view
#' per call.
#'
#' @section Errors:
#' Classed conditions (`nestimate_compare_*`): `too_few`, `bad_input`,
#' `dim_mismatch`, `node_mismatch`, `na_weights`, `reference_unknown`,
#' `labels_length`, `scaling_domain`, `scaling_inference`,
#' `test_unsupported`, `unknown_pair`, `no_nodes`; warning
#' `unknown_measure`.
#'
#' @seealso [compare_model()] (two-network predecessor), [permutation()],
#'   [bayes_compare()], [vertex_compare()], [subtract_networks()].
#' @examples
#' # Regulation networks for the three courses, compared pairwise.
#' courses <- build_network(group_regulation_long, method = "relative",
#'                          actor = "Actor", action = "Action", time = "Time",
#'                          group = "Course")
#' cmp <- compare_networks(courses)
#' cmp
#' summary(cmp)
#' # `pair` takes the two network names, in either order.
#' edge_differences(cmp, pair = c("A", "B"))
#' global_differences(cmp, pair = c("A", "B"))
#' plot(cmp)                                    # each network once
#' plot(cmp, type = "difference")               # signed difference per pair
#' plot(cmp, type = "edges", pair = c("A", "B"))
#' \donttest{
#' # High against low achievers, with a permutation test on every edge.
#' achievers <- build_network(group_regulation_long, method = "relative",
#'                            actor = "Actor", action = "Action",
#'                            time = "Time", group = "Achiever")
#' cmp_perm <- compare_networks(achievers, test = "permutation",
#'                              iter = 100, seed = 1)
#' summary(cmp_perm)
#' plot(cmp_perm, type = "inference", top_n = 10)
#' }
#' @export
compare_networks <- function(...,
                             reference = NULL,
                             scaling = c("none", "minmax", "max", "rank",
                                         "zscore", "robust", "log", "log1p",
                                         "softmax", "quantile", "frobenius",
                                         "row"),
                             measures = c("InStrength", "OutStrength",
                                          "Betweenness"),
                             labels = NULL,
                             test = "none",
                             iter = 1000L,
                             alpha = 0.05,
                             adjust = "none",
                             paired = FALSE,
                             rope = NULL,
                             seed = NULL) {
  scaling <- match.arg(scaling)
  stopifnot(
    "`test` must be a character vector" = is.character(test) && length(test) >= 1L,
    "`iter` must be a single integer >= 2" =
      is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 2,
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1,
    "`adjust` must be a single string" =
      is.character(adjust) && length(adjust) == 1L,
    "`paired` must be TRUE or FALSE" = isTRUE(paired) || isFALSE(paired),
    "`rope` must be NULL or a single positive number" =
      is.null(rope) || (is.numeric(rope) && length(rope) == 1L && rope > 0),
    "`seed` must be NULL or a single number" =
      is.null(seed) || (is.numeric(seed) && length(seed) == 1L)
  )
  test <- match.arg(test, c("none", "permutation", "bayes", "bootstrap"),
                    several.ok = TRUE)
  test <- setdiff(test, "none")
  iter <- as.integer(iter)
  measures <- .cn_resolve_measures(measures)

  collected <- .cn_collect(list(...), labels)
  nets <- collected$networks
  net_names <- names(nets)
  n_networks <- length(nets)
  if (n_networks < 2L) {
    .cn_stop("`compare_networks()` needs at least two networks.", "too_few")
  }

  raw_mats <- lapply(nets, `[[`, "weights")
  .cn_validate_matrices(raw_mats, net_names)
  directed <- vapply(nets, function(z) isTRUE(z$directed), logical(1))
  names(directed) <- net_names

  if (length(test) > 0L && scaling != "none") {
    .cn_stop(paste0("Inference is defined on the networks as estimated; use ",
                    "scaling = \"none\" together with `test`."),
             "scaling_inference")
  }
  .cn_check_scaling_domain(raw_mats, scaling)
  mats <- lapply(raw_mats, .apply_compare_scaling, scaling = scaling)
  nets <- Map(function(net, W) { net$weights <- W; net }, nets, mats)

  pairs <- .cn_pairs(net_names, reference)
  pair_rows <- split(pairs, seq_len(nrow(pairs)))
  seeds <- if (is.null(seed)) rep(list(NULL), nrow(pairs)) else
    as.list(seed + seq_len(nrow(pairs)))

  if (length(test) > 0L) {
    .cn_check_test_support(nets, pairs, test)
  }

  edges <- do.call(rbind, lapply(pair_rows, function(pr) {
    .cn_edge_table(mats[[pr$network_a]], mats[[pr$network_b]], pr,
                   both_undirected = !directed[[pr$network_a]] &&
                     !directed[[pr$network_b]])
  }))
  rownames(edges) <- NULL

  global <- do.call(rbind, lapply(pair_rows, function(pr) {
    .cn_global_table(mats[[pr$network_a]], mats[[pr$network_b]], pr)
  }))
  rownames(global) <- NULL

  network_metrics <- do.call(rbind, Map(function(nm, W, d) {
    sm <- .summary_metrics_from_weights(W, directed = d)
    # Zero out floating-point noise (e.g. centralization ~1e-16) so the
    # column prints in fixed notation next to the integer counts.
    data.frame(network = nm, metric = sm$metric, value = .cn_zap(sm$value),
               stringsAsFactors = FALSE)
  }, net_names, mats, directed))
  rownames(network_metrics) <- NULL

  nodes <- NULL
  if (length(measures) > 0L) {
    cents <- Map(function(W, d) {
      cm <- .perm_centralities_mat(W, rownames(W), d, measures)
      dim(cm) <- c(nrow(W), length(measures))
      dimnames(cm) <- list(rownames(W), measures)
      cm
    }, mats, directed)
    nodes <- do.call(rbind, lapply(pair_rows, function(pr) {
      .cn_node_table(cents[[pr$network_a]], cents[[pr$network_b]], pr,
                     measures)
    }))
    rownames(nodes) <- NULL
  }

  differences <- lapply(pair_rows, function(pr) {
    subtract_networks(nets[[pr$network_a]], nets[[pr$network_b]])
  })
  names(differences) <- pairs$pair

  result <- list(
    networks = nets,
    matrices = mats,
    pairs = pairs,
    edges = edges,
    nodes = nodes,
    global = global,
    network_metrics = network_metrics,
    differences = differences,
    scaling = scaling,
    reference = if (is.null(reference)) NULL else pairs$network_a[1L],
    measures = measures,
    test = if (length(test) == 0L) "none" else test,
    iter = iter,
    alpha = alpha,
    adjust = adjust,
    paired = paired,
    rope = rope,
    directed = directed,
    n_networks = n_networks,
    n_pairs = nrow(pairs)
  )

  if (length(test) > 0L) {
    result <- .cn_add_inference(result, nets, pair_rows, seeds, test,
                                iter, alpha, adjust, paired, rope, measures)
  }

  structure(result, class = "net_network_comparison")
}


# ---- Conditions ----------------------------------------------------------

.cn_stop <- function(msg, class) {
  stop(errorCondition(msg, class = c(paste0("nestimate_compare_", class),
                                     "nestimate_compare_error"),
                      call = NULL))
}

.cn_warn <- function(msg, class) {
  warning(warningCondition(msg, class = c(paste0("nestimate_compare_", class),
                                          "nestimate_compare_warning"),
                           call = NULL))
}


# ---- Input collection ------------------------------------------------------

.cn_resolve_measures <- function(measures) {
  if (is.null(measures) || length(measures) == 0L) return(character(0))
  stopifnot("`measures` must be a character vector" = is.character(measures))
  if (length(measures) == 1L && identical(tolower(measures), "all")) {
    return(.centrality_all_measures())
  }
  valid <- .centrality_builtin_measures()
  bad <- setdiff(measures, valid)
  if (length(bad) > 0L) {
    .cn_warn(paste0("Unknown centrality measure(s) ignored: ",
                    paste(bad, collapse = ", "), ". Valid: ",
                    paste(valid, collapse = ", "), "."),
             "unknown_measure")
  }
  intersect(measures, valid)
}

# Flatten `...` into a named list of netobjects. Groups contribute their
# members; an unnamed plain list is recursed once.
.cn_collect <- function(dots, labels) {
  arg_names <- names(dots) %||% rep("", length(dots))
  arg_names[is.na(arg_names)] <- ""
  units <- Map(.cn_expand, dots, arg_names, seq_along(dots))
  units <- unlist(units, recursive = FALSE)
  nets <- lapply(units, `[[`, "net")
  auto_names <- vapply(units, `[[`, character(1), "name")
  if (!is.null(labels)) {
    if (!is.character(labels) || length(labels) != length(nets)) {
      .cn_stop(sprintf("`labels` must be a character vector of length %d (one per network); got length %d.",
                       length(nets), length(labels)), "labels_length")
    }
    auto_names <- labels
  }
  names(nets) <- make.unique(auto_names, sep = "_")
  list(networks = nets)
}

# One argument -> list of list(net = netobject, name = character)
.cn_expand <- function(obj, arg_name, position, depth = 0L) {
  is_group <- inherits(obj, "netobject_group") || inherits(obj, "group_tna")
  if (inherits(obj, "mcml")) {
    obj <- as_tna(obj)
    is_group <- inherits(obj, "netobject_group")
  }
  if (is_group) {
    member_names <- names(obj) %||% paste0(arg_name %||% "group", "_",
                                            seq_along(obj))
    member_names[!nzchar(member_names)] <- paste0("group_", which(!nzchar(member_names)))
    return(unlist(Map(function(member, nm) {
      list(list(net = .cn_as_netobject(member, nm), name = nm))
    }, unclass(obj), member_names), recursive = FALSE))
  }
  plain_list <- is.list(obj) && !is.data.frame(obj) &&
    is.null(attr(obj, "class"))
  if (plain_list) {
    if (depth > 0L) {
      .cn_stop("Nested lists of networks are not supported; pass networks directly.",
               "bad_input")
    }
    inner_names <- names(obj) %||% rep("", length(obj))
    inner_names[!nzchar(inner_names)] <- paste0("network_", which(!nzchar(inner_names)))
    return(unlist(Map(.cn_expand, obj, inner_names, seq_along(obj),
                      MoreArgs = list(depth = depth + 1L)),
                  recursive = FALSE))
  }
  nm <- if (nzchar(arg_name)) arg_name else NULL
  net <- .cn_as_netobject(obj, nm)
  nm <- nm %||% .cn_default_name(net, position)
  list(list(net = net, name = nm))
}

.cn_default_name <- function(net, position) {
  m <- net$method
  if (is.character(m) && length(m) == 1L && nzchar(m) &&
      !m %in% c("matrix", "difference")) {
    return(m)
  }
  paste0("network_", position)
}

# Coerce one input to a netobject.
.cn_as_netobject <- function(obj, name = NULL) {
  if (inherits(obj, "netobject")) return(obj)
  if (inherits(obj, "psychnet") || inherits(obj, "cograph_network")) {
    return(as_netobject(obj))
  }
  if (inherits(obj, "tna")) return(.cn_from_tna(obj))
  if (is.matrix(obj) && is.numeric(obj)) {
    if (nrow(obj) != ncol(obj)) {
      .cn_stop(sprintf("Network %s: matrix must be square (got %d x %d).",
                       name %||% "", nrow(obj), ncol(obj)), "bad_input")
    }
    if (anyNA(obj)) {
      .cn_stop(sprintf("Weights contain NA in network %s.", name %||% "(matrix)"),
               "na_weights")
    }
    W <- .weights_of(obj)
    return(.wrap_netobject(W, method = "matrix",
                           directed = !isSymmetric(unname(W))))
  }
  .cn_stop(paste0("Unsupported input", if (!is.null(name)) paste0(" `", name, "`"),
                  " of class <", paste(class(obj), collapse = "/"),
                  ">. Accepted: netobject, netobject_group, cograph_network, ",
                  "psychnet, mcml, tna, group_tna, square numeric matrix."),
           "bad_input")
}

# A tna object: rebuild through build_network() when its sequence data is
# available so permutation/Bayesian tests are possible; otherwise wrap the
# weights.
.cn_from_tna <- function(x) {
  W <- x$weights
  labels <- x$labels %||% rownames(W)
  if (is.null(rownames(W))) dimnames(W) <- list(labels, labels)
  type <- attr(x, "type") %||% ""
  method <- switch(type, relative = "relative", frequency = "frequency",
                   `co-occurrence` = "co_occurrence", NULL)
  unscaled <- length(attr(x, "scaling") %||% character(0)) == 0L
  seq_data <- x$data
  if (!is.null(method) && unscaled && is.matrix(seq_data)) {
    alphabet <- attr(seq_data, "labels") %||% attr(seq_data, "alphabet") %||% labels
    df <- as.data.frame(matrix(alphabet[as.integer(seq_data)],
                               nrow = nrow(seq_data)),
                        stringsAsFactors = FALSE)
    net <- build_network(df, method = method,
                         params = list(alphabet = labels))
    if (identical(dim(net$weights), dim(W)) &&
        isTRUE(all.equal(unname(net$weights[labels, labels]), unname(W),
                         tolerance = 1e-8))) {
      return(net)
    }
  }
  .wrap_netobject(W, method = if (is.null(method)) "tna" else method,
                  directed = TRUE, inits = x$inits)
}

.cn_validate_matrices <- function(mats, names) {
  dims <- vapply(mats, nrow, integer(1))
  if (length(unique(dims)) != 1L) {
    .cn_stop(paste0("All networks must have the same number of nodes; got ",
                    paste(sprintf("%s = %d", names, dims), collapse = ", "), "."),
             "dim_mismatch")
  }
  ref_nodes <- rownames(mats[[1L]])
  same <- vapply(mats, function(W) identical(rownames(W), ref_nodes) &&
                   identical(colnames(W), ref_nodes), logical(1))
  if (!all(same)) {
    .cn_stop(paste0("Node labels must be identical (same set and order) in ",
                    "every network; differs for: ",
                    paste(names[!same], collapse = ", "), "."),
             "node_mismatch")
  }
  has_na <- vapply(mats, anyNA, logical(1))
  if (any(has_na)) {
    .cn_stop(paste0("Weights contain NA in: ",
                    paste(names[has_na], collapse = ", "), "."),
             "na_weights")
  }
  invisible(TRUE)
}

.cn_check_scaling_domain <- function(mats, scaling) {
  if (scaling == "log" && any(vapply(mats, function(W) any(W <= 0), logical(1)))) {
    .cn_stop("scaling = \"log\" needs strictly positive weights.", "scaling_domain")
  }
  if (scaling == "log1p" && any(vapply(mats, function(W) any(W < -1), logical(1)))) {
    .cn_stop("scaling = \"log1p\" needs weights >= -1.", "scaling_domain")
  }
  invisible(TRUE)
}

.cn_pairs <- function(names, reference) {
  if (is.null(reference)) {
    idx <- utils::combn(length(names), 2L)
    a <- names[idx[1L, ]]
    b <- names[idx[2L, ]]
  } else {
    stopifnot("`reference` must be a single name or index" = length(reference) == 1L)
    ref <- if (is.character(reference)) match(reference, names) else
      as.integer(reference)
    if (is.na(ref) || ref < 1L || ref > length(names)) {
      .cn_stop(paste0("`reference` \"", reference, "\" is not one of: ",
                      paste(names, collapse = ", "), "."),
               "reference_unknown")
    }
    a <- rep(names[ref], length(names) - 1L)
    b <- names[-ref]
  }
  data.frame(pair = paste(a, "vs", b), network_a = a, network_b = b,
             stringsAsFactors = FALSE)
}


# ---- Tables ------------------------------------------------------------------

.cn_tol <- sqrt(.Machine$double.eps)

# Replace |v| < 1e-10 by exactly 0: such values are floating-point residue
# (e.g. a centralization of 1e-16) and would force a whole printed column
# into scientific notation. Real values are untouched.
.cn_zap <- function(v) {
  v[!is.na(v) & abs(v) < 1e-10] <- 0
  v
}

# Colour channel: which side is higher.
.cn_higher <- function(diff, name_a, name_b) {
  lev <- c(name_a, name_b, "equal")
  out <- ifelse(diff > .cn_tol, name_a, ifelse(diff < -.cn_tol, name_b, "equal"))
  factor(out, levels = lev)
}

.cn_edge_table <- function(A, B, pr, both_undirected = FALSE) {
  nodes <- rownames(A)
  n <- length(nodes)
  from_i <- rep(seq_len(n), each = n)
  to_i <- rep(seq_len(n), times = n)
  keep <- if (both_undirected) from_i <= to_i else rep(TRUE, n * n)
  a <- as.vector(t(A))[keep]
  b <- as.vector(t(B))[keep]
  d <- a - b
  both_zero <- abs(a) <= .cn_tol & abs(b) <= .cn_tol
  a_zero <- abs(a) <= .cn_tol
  b_zero <- abs(b) <= .cn_tol
  status <- ifelse(!a_zero & !b_zero, "both",
            ifelse(!a_zero & b_zero, "only_a",
            ifelse(a_zero & !b_zero, "only_b", "neither")))
  rank_a <- rank(a, ties.method = "average")
  rank_b <- rank(b, ties.method = "average")
  data.frame(
    pair = pr$pair, network_a = pr$network_a, network_b = pr$network_b,
    from = nodes[from_i[keep]], to = nodes[to_i[keep]],
    weight_a = a, weight_b = b,
    diff = d, abs_diff = abs(d),
    rel_diff = ifelse(both_zero, NA_real_, abs(d) / (abs(a) + abs(b))),
    ratio = ifelse(b_zero, NA_real_, a / b),
    log_ratio = ifelse(a < 0 | b < 0, NA_real_, log1p(a) - log1p(b)),
    rank_a = rank_a, rank_b = rank_b, rank_diff = rank_a - rank_b,
    percentile_diff = stats::ecdf(a)(a) - stats::ecdf(b)(b),
    status = factor(status, levels = c("both", "only_a", "only_b", "neither")),
    higher = .cn_higher(d, pr$network_a, pr$network_b),
    stringsAsFactors = FALSE
  )
}

.cn_global_labels <- function() {
  data.frame(
    category = c(rep("Weight Deviations", 6L), rep("Correlations", 4L),
                 rep("Dissimilarities", 5L), rep("Similarities", 5L),
                 rep("Pattern Similarities", 2L)),
    metric = c("Mean Abs. Diff.", "Median Abs. Diff.", "RMS Diff.",
               "Max Abs. Diff.", "Rel. Mean Abs. Diff.", "CV Ratio",
               "Pearson", "Spearman", "Kendall", "Distance",
               "Euclidean", "Manhattan", "Canberra", "Bray-Curtis", "Frobenius",
               "Cosine", "Jaccard", "Dice", "Overlap", "RV",
               "Rank Agreement", "Sign Agreement"),
    key = c("mean_abs_diff", "median_abs_diff", "rms_diff", "max_abs_diff",
            "rel_mean_abs", "cv_ratio",
            "pearson", "spearman", "kendall", "distance_cor",
            "euclidean", "manhattan", "canberra", "bray_curtis", "frobenius",
            "cosine", "jaccard", "dice", "overlap", "rv",
            "rank_agreement", "sign_agreement"),
    stringsAsFactors = FALSE
  )
}

.cn_global_table <- function(A, B, pr) {
  lab <- .cn_global_labels()
  sim <- .network_similarity(A, B)
  data.frame(
    pair = pr$pair, network_a = pr$network_a, network_b = pr$network_b,
    category = lab$category, metric = lab$metric, key = lab$key,
    value = .cn_zap(unname(sim[lab$key])),
    stringsAsFactors = FALSE
  )
}

.cn_node_table <- function(CA, CB, pr, measures) {
  nodes <- rownames(CA)
  do.call(rbind, lapply(measures, function(m) {
    a <- unname(CA[, m])
    b <- unname(CB[, m])
    d <- a - b
    data.frame(
      pair = pr$pair, network_a = pr$network_a, network_b = pr$network_b,
      node = nodes, measure = m,
      value_a = a, value_b = b, diff = d, abs_diff = abs(d),
      rank_a = rank(a, ties.method = "average"),
      rank_b = rank(b, ties.method = "average"),
      higher = .cn_higher(d, pr$network_a, pr$network_b),
      stringsAsFactors = FALSE
    )
  }))
}


# ---- Inference -------------------------------------------------------------

.cn_supports <- function(net, backend) {
  method <- .resolve_method_alias(net$method %||% "")
  has_data <- is.data.frame(net$data) && ncol(net$data) > 0L
  switch(backend,
    permutation = has_data,
    bayes = has_data && method %in% c("relative", "frequency") &&
      !is.null(net$frequency_matrix),
    bootstrap = TRUE
  )
}

.cn_check_test_support <- function(nets, pairs, test) {
  problems <- unlist(lapply(test, function(bk) {
    ok <- vapply(nets, .cn_supports, logical(1), backend = bk)
    if (all(ok)) return(character(0))
    sprintf("%s: %s", bk, paste(names(nets)[!ok], collapse = ", "))
  }))
  methods <- vapply(nets, function(z) .resolve_method_alias(z$method %||% ""),
                    character(1))
  if (any(test %in% c("permutation", "bayes")) && length(unique(methods)) > 1L) {
    problems <- c(problems, paste0("permutation/bayes need one estimation method; got ",
                                   paste(unique(methods), collapse = ", ")))
  }
  if (length(problems) > 0L) {
    .cn_stop(paste0("Requested inference is not supported by every network ",
                    "(need build_network() output with data; Bayesian needs ",
                    "relative/frequency counts):\n  ",
                    paste(problems, collapse = "\n  ")),
             "test_unsupported")
  }
  invisible(TRUE)
}

# rbind two data.frames, filling columns missing on either side with NA.
.cn_rbind_fill <- function(a, b) {
  cols <- union(names(a), names(b))
  fill <- function(d) {
    missing <- setdiff(cols, names(d))
    d[missing] <- NA
    d[cols]
  }
  out <- rbind(fill(a), fill(b))
  rownames(out) <- NULL
  out
}

# Row-major cell vector of a nodes x nodes matrix, aligned with .cn_edge_table
.cn_cells <- function(M, nodes, keep) as.vector(t(M[nodes, nodes]))[keep]

.cn_add_inference <- function(result, nets, pair_rows, seeds, test, iter,
                              alpha, adjust, paired, rope, measures) {
  edges <- result$edges
  nodes_tab <- result$nodes
  global <- result$global
  node_labels <- rownames(result$matrices[[1L]])
  edge_keep <- function(pr) {
    both_undirected <- !result$directed[[pr$network_a]] &&
      !result$directed[[pr$network_b]]
    n <- length(node_labels)
    if (both_undirected) rep(seq_len(n), each = n) <= rep(seq_len(n), times = n)
    else rep(TRUE, n * n)
  }

  if ("permutation" %in% test) {
    runs <- Map(function(pr, sd) {
      permutation(nets[[pr$network_a]], nets[[pr$network_b]], iter = iter,
                  alpha = alpha, paired = paired, adjust = adjust,
                  measures = if (length(measures) > 0L) measures else NULL,
                  seed = sd)
    }, pair_rows, seeds)
    edge_add <- do.call(rbind, Map(function(pr, r) {
      keep <- edge_keep(pr)
      data.frame(perm_effect = .cn_cells(r$effect_size, node_labels, keep),
                 perm_p = .cn_cells(r$p_values, node_labels, keep))
    }, pair_rows, runs))
    edges$perm_effect <- edge_add$perm_effect
    edges$perm_p <- edge_add$perm_p
    edges$perm_sig <- edges$perm_p < alpha
    if (!is.null(nodes_tab)) {
      node_add <- do.call(rbind, Map(function(pr, r) {
        st <- r$centralities$stats
        key <- paste(st$state, st$centrality)
        want <- paste(nodes_tab$node[nodes_tab$pair == pr$pair],
                      nodes_tab$measure[nodes_tab$pair == pr$pair])
        idx <- match(want, key)
        data.frame(perm_effect = st$effect_size[idx], perm_p = st$p_value[idx])
      }, pair_rows, runs))
      nodes_tab$perm_effect <- node_add$perm_effect
      nodes_tab$perm_p <- node_add$perm_p
      nodes_tab$perm_sig <- nodes_tab$perm_p < alpha
    }
    global_add <- do.call(rbind, Map(function(pr, r) {
      g <- r$global
      data.frame(pair = pr$pair, network_a = pr$network_a,
                 network_b = pr$network_b,
                 category = "Global Test (permutation)",
                 metric = c("M (global strength)", "S (max edge)"),
                 key = c("perm_M", "perm_S"),
                 value = g$observed, perm_p = g$p_value,
                 stringsAsFactors = FALSE)
    }, pair_rows, runs))
    global$perm_p <- NA_real_
    global <- .cn_rbind_fill(global, global_add)
    global$perm_sig <- ifelse(is.na(global$perm_p), NA, global$perm_p < alpha)
  }

  if ("bayes" %in% test) {
    runs <- Map(function(pr, sd) {
      bayes_compare(nets[[pr$network_a]], nets[[pr$network_b]], prior = 0.5,
                    draws = iter, ci = 1 - alpha, seed = sd)
    }, pair_rows, seeds)
    edge_add <- do.call(rbind, Map(function(pr, r) {
      keep <- edge_keep(pr)
      bdiff <- .cn_cells(r$diff, node_labels, keep)
      es <- .cn_cells(r$effect_size, node_labels, keep)
      post_sd <- ifelse(abs(es) > .cn_tol, bdiff / es, NA_real_)
      out <- data.frame(
        bayes_diff = bdiff,
        bayes_ci_lower = .cn_cells(r$ci_lower, node_labels, keep),
        bayes_ci_upper = .cn_cells(r$ci_upper, node_labels, keep),
        bayes_pd = .cn_cells(r$p_difference, node_labels, keep),
        bayes_p = .cn_cells(r$p_bayes, node_labels, keep),
        bayes_sig = .cn_cells(r$sig, node_labels, keep)
      )
      if (!is.null(rope)) {
        out$bayes_p_rope <- stats::pnorm(rope, bdiff, post_sd) -
          stats::pnorm(-rope, bdiff, post_sd)
        inside <- out$bayes_ci_lower >= -rope & out$bayes_ci_upper <= rope
        outside <- out$bayes_ci_lower > rope | out$bayes_ci_upper < -rope
        out$bayes_decision <- factor(
          ifelse(outside, "different", ifelse(inside, "equivalent", "undecided")),
          levels = c("different", "equivalent", "undecided"))
      }
      out
    }, pair_rows, runs))
    edges <- cbind(edges, edge_add)
  }

  if ("bootstrap" %in% test) {
    runs <- Map(function(pr, sd) {
      vertex_compare(nets[[pr$network_a]], nets[[pr$network_b]], iter = iter,
                     ci_level = alpha, seed = sd,
                     labels = c(pr$network_a, pr$network_b))
    }, pair_rows, seeds)
    global_add <- do.call(rbind, Map(function(pr, r) {
      s <- r$summary
      data.frame(pair = pr$pair, network_a = pr$network_a,
                 network_b = pr$network_b,
                 category = "Structure (bootstrap)",
                 metric = s$statistic, key = paste0("boot_", s$statistic),
                 value = s$diff,
                 boot_observed_a = s[[2L]], boot_observed_b = s[[3L]],
                 boot_se = s$se_diff, boot_ci_lower = s$ci_lower,
                 boot_ci_upper = s$ci_upper, boot_z = s$z, boot_p = s$p_value,
                 stringsAsFactors = FALSE)
    }, pair_rows, runs))
    global <- .cn_rbind_fill(global, global_add)
    global$boot_sig <- ifelse(is.na(global$boot_p), NA, global$boot_p < alpha)
  }

  # Unified evidence channel for plots
  if ("permutation" %in% test) {
    edges$sig <- edges$perm_sig
    edges$evidence <- "permutation"
    if (!is.null(nodes_tab)) {
      nodes_tab$sig <- nodes_tab$perm_sig
      nodes_tab$evidence <- "permutation"
    }
  } else if ("bayes" %in% test) {
    edges$sig <- edges$bayes_sig
    edges$evidence <- "bayes"
  }

  result$edges <- edges
  result$nodes <- nodes_tab
  result$global <- global
  result
}


# ---- print, summary, table verbs -------------------------------------------

.cn_header <- function(x) {
  test_txt <- if (identical(x$test, "none")) "descriptive" else
    sprintf("%s, iter = %d, alpha = %s, adjust = %s",
            paste(x$test, collapse = " + "), x$iter, format(x$alpha), x$adjust)
  net_lines <- vapply(names(x$networks), function(nm) {
    sprintf("%s (%d nodes, %s)", nm, nrow(x$matrices[[nm]]),
            if (x$directed[[nm]]) "directed" else "undirected")
  }, character(1))
  paste0(
    sprintf("Network comparison (%s): %d networks, %d pairs, scaling = %s\n",
            test_txt, x$n_networks, x$n_pairs, x$scaling),
    if (!is.null(x$reference)) sprintf("  reference: %s\n", x$reference) else "",
    "  networks: ", paste(net_lines, collapse = ", ")
  )
}

#' @rdname compare_networks
#' @param x A `net_network_comparison` object.
#' @param digits Decimals shown. Default `2`.
#' @return `print()` returns `x` invisibly.
#' @export
print.net_network_comparison <- function(x, digits = 2L, ...) {
  digits <- .cn_check_digits(digits)
  cat(.cn_header(x), "\n\n", sep = "")
  s <- .cn_summary(x)
  show <- data.frame(
    pair = s$pair,
    pearson = .cn_fmt(s$pearson, digits),
    `mean|diff|` = .cn_fmt(s$mean_abs_diff, digits),
    `max|diff|` = .cn_fmt(s$max_abs_diff, digits),
    `largest change` = sprintf("%s (%s higher)", s$top_edge, s$top_edge_higher),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if ("n_sig_edges" %in% names(s)) show$`sig edges` <- s$n_sig_edges
  print(show, row.names = FALSE, right = FALSE)
  cat("\nTables: summary(x), edge_differences(x), node_differences(x), ",
      "global_differences(x), network_metrics(x). Plot: plot(x, type = ...)\n",
      sep = "")
  invisible(x)
}

# One row per pair (backs summary() and print()).
.cn_summary <- function(object) {
  e <- object$edges
  g <- object$global
  per_pair <- lapply(split(e, factor(e$pair, levels = object$pairs$pair)),
                     function(d) {
    top <- order(-d$abs_diff, d$from, d$to)[1L]
    data.frame(
      n_cells = nrow(d),
      n_differing = sum(d$abs_diff > .cn_tol),
      share_higher_a = mean(d$higher == d$network_a[1L]),
      share_higher_b = mean(d$higher == d$network_b[1L]),
      mean_abs_diff = mean(d$abs_diff),
      max_abs_diff = max(d$abs_diff),
      top_edge = paste(d$from[top], "->", d$to[top]),
      top_edge_higher = as.character(d$higher[top]),
      stringsAsFactors = FALSE
    )
  })
  out <- cbind(object$pairs, do.call(rbind, per_pair))
  pick <- function(key) {
    v <- g$value[g$key == key]
    v[match(object$pairs$pair, g$pair[g$key == key])]
  }
  out$pearson <- pick("pearson")
  out$spearman <- pick("spearman")
  out$cosine <- pick("cosine")
  out$jaccard <- pick("jaccard")
  out <- out[c("pair", "network_a", "network_b", "n_cells", "n_differing",
               "share_higher_a", "share_higher_b", "mean_abs_diff",
               "max_abs_diff", "pearson", "spearman", "cosine", "jaccard",
               "top_edge", "top_edge_higher")]
  if ("sig" %in% names(e)) {
    out$n_sig_edges <- vapply(split(e$sig, factor(e$pair, levels = object$pairs$pair)),
                              function(v) sum(v, na.rm = TRUE), integer(1))
    if (!is.null(object$nodes) && "sig" %in% names(object$nodes)) {
      nd <- object$nodes
      out$n_sig_nodes <- vapply(split(nd$sig, factor(nd$pair, levels = object$pairs$pair)),
                                function(v) sum(v, na.rm = TRUE), integer(1))
    }
  }
  if ("perm_p" %in% names(g)) {
    out$m_stat <- pick("perm_M")
    out$m_p <- g$perm_p[g$key == "perm_M"][match(object$pairs$pair, g$pair[g$key == "perm_M"])]
    out$s_stat <- pick("perm_S")
    out$s_p <- g$perm_p[g$key == "perm_S"][match(object$pairs$pair, g$pair[g$key == "perm_S"])]
  }
  rownames(out) <- NULL
  out
}

#' Tables of a network comparison
#'
#' Named accessors for the tables inside a `net_network_comparison` object
#' (from [compare_networks()]). Each returns a plain data frame with one row
#' per unit and no row names. Numeric columns are rounded to `digits`
#' decimals; p-value columns are never rounded.
#'
#' @param x,object A `net_network_comparison` object (for `print.net_table()`,
#'   a table returned by one of these verbs).
#' @param pair Optional selection of comparisons. Any of: the pair name(s) as
#'   printed (`"A vs B"`); the two network names (`c("A", "B")`, either
#'   order); one network name (`"A"`, every pair it takes part in); or
#'   index/indices into the pair table (`1`, `c(1, 3)`). `NULL` keeps all.
#' @param measure Optional centrality measure name(s) to keep.
#' @param digits Decimals kept in numeric columns. Default `2`.
#' @param ... Ignored.
#' @return
#' * `summary()`: one row per pair -- `pair`, `network_a`, `network_b`,
#'   `n_cells`, `n_differing`, `share_higher_a`, `share_higher_b`,
#'   `mean_abs_diff`, `max_abs_diff`, `pearson`, `spearman`, `cosine`,
#'   `jaccard`, `top_edge`, `top_edge_higher`, and with inference
#'   `n_sig_edges`, `n_sig_nodes`, `m_stat`, `m_p`, `s_stat`, `s_p`. The
#'   per-edge, per-node, per-metric and per-network tables are the four
#'   verbs below.
#' * `edge_differences()`: one row per pair x transition: `pair`,
#'   `network_a`, `network_b`, `from`, `to`, `weight_a`, `weight_b`, `diff`,
#'   `abs_diff`, `rel_diff`, `ratio`, `log_ratio`, `rank_a`, `rank_b`,
#'   `rank_diff`, `percentile_diff`, `status`, `higher`, and the inference
#'   columns when `test` was used.
#' * `node_differences()`: one row per pair x state x measure: `pair`,
#'   `network_a`, `network_b`, `node`, `measure`, `value_a`, `value_b`,
#'   `diff`, `abs_diff`, `rank_a`, `rank_b`, `higher`, plus inference columns.
#'   Errors (class `nestimate_compare_no_nodes`) when `compare_networks()`
#'   was called with no `measures`.
#' * `global_differences()`: one row per pair x metric: `pair`, `network_a`,
#'   `network_b`, `category`, `metric`, `key`, `value`, plus inference rows
#'   and columns.
#' * `network_metrics()`: one row per network x structural metric:
#'   `network`, `metric`, `value`.
#'
#' Every table is a data frame of class `net_table` whose `print()` shows
#' whole numbers without decimals, exact zeros as `0`, other values with
#' `digits` decimals, and p-values with three decimals.
#' @examples
#' achievers <- build_network(group_regulation_long, method = "relative",
#'                            actor = "Actor", action = "Action",
#'                            time = "Time", group = "Achiever")
#' cmp <- compare_networks(achievers)
#' summary(cmp)
#' edge_differences(cmp)
#' edge_differences(cmp, digits = 4)
#' node_differences(cmp, measure = "InStrength")
#' global_differences(cmp)
#' network_metrics(cmp)
#' @name comparison_tables
NULL

#' @rdname comparison_tables
#' @export
summary.net_network_comparison <- function(object, pair = NULL, digits = 2L,
                                           ...) {
  .cn_check_object(object)
  .cn_tidy(.cn_pair_filter(.cn_summary(object), object, pair), digits)
}

#' @rdname comparison_tables
#' @export
edge_differences <- function(x, pair = NULL, digits = 2L) {
  .cn_check_object(x)
  .cn_tidy(.cn_pair_filter(x$edges, x, pair), digits)
}

#' @rdname comparison_tables
#' @export
node_differences <- function(x, pair = NULL, measure = NULL, digits = 2L) {
  .cn_check_object(x)
  if (is.null(x$nodes)) {
    .cn_stop("No node table: `compare_networks()` was called with no `measures`.",
             "no_nodes")
  }
  tab <- .cn_pair_filter(x$nodes, x, pair)
  if (!is.null(measure)) {
    unknown <- setdiff(measure, x$measures)
    if (length(unknown) > 0L) {
      .cn_stop(paste0("Unknown measure(s): ", paste(unknown, collapse = ", "),
                      ". Available: ", paste(x$measures, collapse = ", "), "."),
               "unknown_measure")
    }
    tab <- tab[tab$measure %in% measure, , drop = FALSE]
  }
  .cn_tidy(tab, digits)
}

#' @rdname comparison_tables
#' @export
global_differences <- function(x, pair = NULL, digits = 2L) {
  .cn_check_object(x)
  .cn_tidy(.cn_pair_filter(x$global, x, pair), digits)
}

#' @rdname comparison_tables
#' @export
network_metrics <- function(x, digits = 2L) {
  .cn_check_object(x)
  .cn_tidy(x$network_metrics, digits)
}

.cn_check_object <- function(x) {
  stopifnot("`x` must be a net_network_comparison (from compare_networks())" =
              inherits(x, "net_network_comparison"))
  invisible(TRUE)
}

.cn_check_digits <- function(digits) {
  stopifnot("`digits` must be a single non-negative integer" =
              is.numeric(digits) && length(digits) == 1L && !is.na(digits) &&
              digits >= 0)
  as.integer(digits)
}

# Round numeric columns to `digits` (p-value columns untouched), drop row
# names, and tag the frame so print() shows exact zeros as "0" rather than
# "0.00". The object is still a plain data.frame for every other purpose.
.cn_tidy <- function(tab, digits = 2L) {
  digits <- .cn_check_digits(digits)
  is_p <- grepl("(^|_)p($|_)|_pd$", names(tab))
  num <- vapply(tab, is.numeric, logical(1)) & !is_p
  tab[num] <- lapply(tab[num], round, digits = digits)
  rownames(tab) <- NULL
  attr(tab, "digits") <- digits
  class(tab) <- c("net_table", "data.frame")
  tab
}

#' @rdname comparison_tables
#' @export
print.net_table <- function(x, digits = attr(x, "digits") %||% 2L, ...) {
  digits <- .cn_check_digits(digits)
  shown <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  is_p <- grepl("(^|_)p($|_)|_pd$", names(shown))
  num <- vapply(shown, is.numeric, logical(1))
  fmt_value <- function(v) {
    whole <- !is.na(v) & abs(v - round(v)) < .cn_tol
    out <- .cn_fmt(v, digits)
    out[whole] <- format(round(v[whole]), scientific = FALSE, trim = TRUE)
    out
  }
  shown[num & !is_p] <- lapply(shown[num & !is_p], fmt_value)
  shown[num & is_p] <- lapply(shown[num & is_p], function(v)
    ifelse(is.na(v), NA_character_, formatC(v, digits = 3L, format = "f")))
  shown[] <- lapply(shown, function(v) { v <- as.character(v); v[is.na(v)] <- "NA"; v })
  print.data.frame(shown, row.names = FALSE, right = TRUE)
  invisible(x)
}

# Format numbers for print/plot labels: anything that rounds to zero at
# `digits` prints as "0" (never "-0.00"), others with `digits` decimals;
# `signed = TRUE` adds a leading sign. Rounding happens before the zero test
# so a value of -0.0004 is reported as 0 rather than as a signed zero.
.cn_fmt <- function(v, digits = 2L, signed = FALSE) {
  fmt <- paste0("%", if (signed) "+" else "", ".", digits, "f")
  r <- round(v, digits)
  ifelse(is.na(v), NA_character_,
         ifelse(abs(r) < .cn_tol, "0", sprintf(fmt, r)))
}

.cn_pair_filter <- function(tab, x, pair) {
  if (is.null(pair)) return(tab)
  .cn_pair_subset(tab, x, pair)
}

# Resolve `pair` to pair names. Accepts, in this order:
#   * NULL                      -> every pair
#   * the composed names        -> "A vs B", or several of them
#   * two network names         -> c("A", "C"), in either order
#   * one network name          -> every pair that network takes part in
#   * indices into x$pairs      -> 1, or c(1, 3)
# Composed names win, so a network named like a pair can never be ambiguous.
.cn_match_pairs <- function(x, pair) {
  if (is.null(pair)) return(x$pairs$pair)
  known <- x$pairs$pair
  nets <- names(x$networks)
  if (is.numeric(pair)) {
    idx <- as.integer(pair)
    if (anyNA(idx) || any(idx < 1L) || any(idx > length(known))) {
      .cn_stop(paste0("`pair` index out of range; there are ", length(known),
                      " pairs: ", paste(known, collapse = ", "), "."),
               "unknown_pair")
    }
    return(known[idx])
  }
  if (!is.character(pair)) {
    .cn_stop("`pair` must be pair name(s), network name(s), or index/indices.",
             "unknown_pair")
  }
  if (all(pair %in% known)) return(unique(pair))
  if (all(pair %in% nets)) {
    hit <- (x$pairs$network_a %in% pair & x$pairs$network_b %in% pair) |
      (length(pair) == 1L &
         (x$pairs$network_a %in% pair | x$pairs$network_b %in% pair))
    if (any(hit)) return(known[hit])
    .cn_stop(paste0("No comparison of ", paste(pair, collapse = " and "),
                    ". Available: ", paste(known, collapse = ", "), "."),
             "unknown_pair")
  }
  unknown <- setdiff(pair, c(known, nets))
  .cn_stop(paste0("Unknown pair(s): ", paste(unknown, collapse = ", "),
                  ". Available pairs: ", paste(known, collapse = ", "),
                  "; networks: ", paste(nets, collapse = ", "), "."),
           "unknown_pair")
}

# Filter a table to one or more pairs, however `pair` was written.
.cn_pair_subset <- function(tab, x, pair) {
  tab[tab$pair %in% .cn_match_pairs(x, pair), , drop = FALSE]
}
