# ---- Centrality Stability ----

#' Centrality Stability Coefficient (CS-coefficient)
#'
#' @description
#' `centrality_stability()` estimates how stable a network's centrality
#' indices are under case-dropping (Epskamp, Borsboom & Fried, 2018). For
#' each proportion in `drop_prop`, that share of the sequences (or, for
#' association networks, of the rows) is removed at random `iter` times,
#' the network is re-estimated on what is left, and its centralities are
#' correlated with those of the full network. The CS-coefficient of a
#' measure is the largest proportion that can be dropped while the
#' correlation stays at or above `threshold` in at least a `certainty`
#' share of the draws. It takes a `netobject` from [build_network()] of any
#' method, a bare `cograph_network`, or a `netobject_group` or `mcml`, whose
#' networks are assessed one by one.
#'
#' The CS-coefficient says how well the nodes' centralities (their values,
#' or with `method = "spearman"` or `"kendall"` their ordering) survive the
#' loss of cases; Epskamp et al. (2018) suggest that it should
#' not be below 0.25, and preferably above 0.5, before centrality
#' differences are interpreted. It does not say whether two nodes differ in
#' centrality, nor how precise single edges are: check the edges with
#' [bootstrap_network()], and for `glasso` networks the centrality
#' differences with [boot_glasso()].
#'
#' @details
#' For transition methods, per-sequence count matrices are computed once
#' and each draw only re-sums them; the network's `scaling` and `threshold`
#' are applied to every re-estimated network as to the original. For the
#' other methods the estimator is re-run on each subset with the network's
#' `params` (after the multilevel decomposition, when the network was built
#' with a `level` and an id column and its estimator is undirected).
#' The built-in measures are computed from the weight matrix
#' directly, without igraph, as in [net_centrality()]: strength from row
#' and column sums, betweenness and closeness from Floyd-Warshall shortest
#' paths. A network that rests on one sequence, or one-hot data with no
#' actor or session, has nothing to drop and raises an error of class
#' `nestimate_single_sequence`.
#'
#' Every draw counts towards `certainty`: a draw whose correlation is
#' undefined (the measure is constant in it), and a draw whose network
#' could not be estimated, count as not reaching `threshold`. Failed
#' estimates are reported in a warning of class
#' `nestimate_resample_failed`.
#'
#' Measures without variance are handled by two different rules, both
#' long-standing behaviour. When *some* requested measures do not vary
#' across the nodes of the original network (or are undefined on it), they
#' are dropped, and the result covers only the retained ones. When *every*
#' requested measure is so, a warning is issued and *all* requested names
#' are returned with a CS-coefficient of 0 and correlations that are all
#' `NA`. A measure counts as constant when its values span no more than
#' the square root of the machine epsilon relative to their size, so
#' `"OutStrength"` of a `"relative"` network with `loops = TRUE`, which is
#' 1 for every node up to rounding error, is dropped rather than correlated
#' on that error. The same rule makes a measure constant in a draw an
#' undefined correlation.
#'
#' @param x The network: a `netobject` from [build_network()], a bare
#'   `cograph_network` (converted to one), or a `netobject_group` or `mcml`,
#'   each of whose networks is assessed. Any other object, or a network
#'   built without its data, raises an error. For the `print()` methods: a
#'   `net_stability` or `net_stability_group`; for the `plot()` method: a
#'   `net_stability`.
#' @param measures The centrality measures to assess. The default is
#'   `"InStrength"`, `"OutStrength"` and `"Betweenness"`; `NULL` takes
#'   [net_centrality()]'s default (`"InStrength"`, `"Betweenness"`,
#'   `"Diffusion"`), and `"all"` every built-in measure: `"OutStrength"`,
#'   `"InStrength"`, `"ClosenessIn"`, `"ClosenessOut"`, `"Closeness"`,
#'   `"Betweenness"`, `"BetweennessRSP"`, `"Diffusion"` and `"Clustering"`
#'   ([net_centrality()] defines them). `"InCloseness"` and
#'   `"OutCloseness"` are accepted for `"ClosenessIn"` and `"ClosenessOut"`.
#'   Other names need `centrality_fn`; without it they raise an error.
#' @param iter The number of random subsets drawn at each drop proportion,
#'   at least 2. Default `1000`.
#' @param drop_prop The proportions of cases to drop, each strictly between
#'   0 and 1. Default `seq(0.1, 0.9, by = 0.1)`. A proportion that would
#'   leave fewer than 2 cases is skipped and counts as not reaching
#'   `threshold`. The CS-coefficient is the largest proportion that meets
#'   `certainty`, in whatever order they are given.
#' @param threshold The correlation with the full network's centralities
#'   that a draw must reach, between 0 and 1. Default `0.7`.
#' @param certainty The share of draws at a drop proportion that must reach
#'   `threshold`, between 0 and 1. Default `0.95`. It is a share of all
#'   draws: one with an undefined correlation counts as not reaching
#'   `threshold`.
#' @param method The correlation between the full and the subset
#'   centralities: `"pearson"` (default), `"spearman"` or `"kendall"`.
#' @param centrality_fn A function of your own that takes the weight matrix
#'   and returns a named list of centrality vectors, such as
#'   `list(my_metric = ...)`. It is called as `centrality_fn(mat)` and used
#'   only for the requested measures that are not built-in; list their
#'   names in `measures`. Default `NULL`.
#' @param loops Whether self-loops (staying in a state) count. With `FALSE`
#'   (default) they are left out of every measure; the stored matrix is not
#'   changed.
#' @param normalize Whether to range-normalise every measure to 0--1, as
#'   `tna::centralities(normalize = TRUE)` does. Default `FALSE`.
#' @param invert Whether a heavier edge is a shorter step in the
#'   shortest-path measures (distance = 1 / weight), as in `tna`. Default
#'   `TRUE`; with `FALSE` the weight is the distance. It affects betweenness
#'   and closeness only, not `"BetweennessRSP"`.
#' @param normalize_diffusion Whether to range-normalise `Diffusion` even
#'   when `normalize = FALSE`. Default `TRUE`.
#' @param seed An integer that makes the random subsets reproducible.
#'   Default `NULL`. The session's own random stream is left as it was; a
#'   group's networks
#'   are each started from the same seed.
#'
#' @return A `net_stability` object. `print()` shows the CS-coefficient of
#'   each measure; `summary()` returns one row per measure and drop
#'   proportion, with the mean and SD of the correlations (over the draws
#'   where it is defined) and `prop_above`, the share of all draws that
#'   reach `threshold`; `plot()` draws the mean correlation against the
#'   proportion dropped, one line per measure with a band of one SD and a
#'   dashed line at `threshold`, and returns the ggplot invisibly. For a
#'   `netobject_group` or `mcml`, the result is a `net_stability_group`, a
#'   named list of `net_stability` objects, one per network, whose `print()`
#'   shows the CS-coefficients as a networks-by-measures table and whose
#'   `summary()` stacks the per-network summaries with a leading `group`
#'   column. The `print()` methods return their input invisibly. The
#'   `net_stability`'s fields:
#' \describe{
#'   \item{cs}{The CS-coefficients, a named numeric vector with one per
#'     retained measure.}
#'   \item{correlations}{A named list, one per retained measure, of
#'     `iter` x `length(drop_prop)` matrices of correlations.}
#'   \item{measures}{The measures assessed (see the rule for measures
#'     without variance under *Details*).}
#'   \item{drop_prop, threshold, certainty, iter, method}{The settings
#'     used.}
#' }
#'
#' @examples
#' seqs <- data.frame(
#'   T1 = c("plan", "code", "debug", "plan", "test", "code"),
#'   T2 = c("code", "debug", "code", "plan", "code", "test"),
#'   T3 = c("debug", "code", "plan", "code", "debug", "plan"),
#'   T4 = c("test", "plan", "test", "debug", "plan", "code")
#' )
#' net <- build_network(seqs, method = "relative")
#' cs <- centrality_stability(net, iter = 10, drop_prop = 0.3, seed = 1)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' cs <- centrality_stability(net, iter = 100, seed = 42,
#'   measures = c("InStrength", "OutStrength"))
#' print(cs)
#' }
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods} 50(1), 195-212.
#' \doi{10.3758/s13428-017-0862-1}
#'
#' @seealso [net_centrality()], [bootstrap_network()],
#'   [network_reliability()], [build_network()]
#'
#' @importFrom stats cor sd
#' @export
centrality_stability <- function(x,
                                 measures = c("InStrength", "OutStrength",
                                              "Betweenness"),
                                 iter = 1000L,
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 threshold = 0.7,
                                 certainty = 0.95,
                                 method = "pearson",
                                 centrality_fn = NULL,
                                 loops = FALSE,
                                 normalize = FALSE,
                                 invert = TRUE,
                                 normalize_diffusion = TRUE,
                                 seed = NULL) {

  # ---- Input validation ----
  measures <- .resolve_measures(measures)   # "all" -> every built-in measure
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(x, "netobject_group")) {
    .stop_single_unit_group(x)
    out <- lapply(x, function(net) {
      centrality_stability(net, measures = measures, iter = iter,
                           drop_prop = drop_prop, threshold = threshold,
                           certainty = certainty, method = method,
                           centrality_fn = centrality_fn, loops = loops,
                           normalize = normalize, invert = invert,
                           normalize_diffusion = normalize_diffusion,
                           seed = seed)
    })
    class(out) <- c("net_stability_group", "list")
    return(out)
  }
  if (!inherits(x, "netobject")) {
    stop(errorCondition("'x' must be a netobject from build_network().",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (is.null(x$data)) {
    stop(errorCondition(
      "netobject does not contain $data. Rebuild with build_network().",
      class = "nestimate_bad_input", call = NULL))
  }
  .stop_single_unit(x)
  .assert(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(drop_prop), all(drop_prop > 0), all(drop_prop < 1),
    is.numeric(threshold), length(threshold) == 1,
    threshold >= 0, threshold <= 1,
    is.numeric(certainty), length(certainty) == 1,
    certainty >= 0, certainty <= 1
  )
  iter <- as.integer(iter)
  method <- match.arg(method, c("pearson", "spearman", "kendall"))

  if (!is.null(centrality_fn)) {
    .assert("centrality_fn must be a function" = is.function(centrality_fn))
  }

  valid_measures <- .centrality_builtin_measures()
  bad <- setdiff(measures, valid_measures)
  # A custom centrality_fn resolves non-builtin measure names (mirrors
  # net_centrality, which has no whitelist). Only reject unknown names
  # when no centrality_fn is supplied; otherwise .compute_centralities
  # routes them through centrality_fn.
  if (length(bad) > 0L && is.null(centrality_fn)) {
    stop(errorCondition(
      paste0(
        "Unknown measures: ",
        paste(bad, collapse = ", "),
        ". Options: ",
        paste(valid_measures, collapse = ", "),
        ". Custom measures require centrality_fn."),
      class = "nestimate_bad_input", call = NULL))
  }

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
  }

  # ---- Setup ----
  net_method <- .resampling_method(.resolve_method_alias(x$method), x$params)
  states <- x$nodes$label
  n_states <- length(states)
  is_transition <- net_method %in% c("relative", "frequency", "co_occurrence")
  is_relative <- .counts_are_relative(net_method, x$params)
  scaling <- x$scaling
  thresh <- x$threshold

  # ---- Original centralities ----
  orig_cents <- .compute_centralities(
    x$weights, states, x$directed, measures,
    loops = loops, normalize = normalize, invert = invert,
    normalize_diffusion = normalize_diffusion, centrality_fn = centrality_fn
  )

  # Drop measures with zero variance (e.g. OutStrength for relative networks)
  # or that are undefined on this network (e.g. Diffusion is NaN on a small
  # cyclic net): an all-NA sd() would otherwise poison `any()` with NA.
  keep <- vapply(measures,
                 function(m) !.is_constant_vector(orig_cents[[m]]),
                 logical(1))
  if (!any(keep)) {
    warning("All centrality measures have zero variance. ",
            "No stability can be assessed.", call. = FALSE)
    result <- list(
      cs = stats::setNames(rep(0, length(measures)), measures),
      correlations = stats::setNames(
        lapply(measures, function(m) {
          matrix(NA_real_, nrow = iter, ncol = length(drop_prop))
        }), measures),
      measures = measures,
      drop_prop = drop_prop,
      threshold = threshold,
      certainty = certainty,
      iter = iter,
      method = method
    )
    class(result) <- "net_stability"
    return(result)
  }
  measures <- measures[keep]
  orig_cents <- orig_cents[measures]

  # ---- Pre-compute for transitions ----
  if (is_transition) {
    transition_data <- .resampling_transition_data(x$data, x$metadata, x$params)
    trans_2d <- .precompute_per_sequence(transition_data, net_method, x$params,
                                         states)
    n_seq <- nrow(trans_2d)
  } else {
    data <- x$data
    n_seq <- nrow(data)
    estimator <- get_estimator(net_method)
    params <- x$params
    level <- x$level
    id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")
  }

  # ---- Build matrix from subset (transition fast path) ----
  build_matrix_transition <- function(idx) {
    counts <- colSums(trans_2d[idx, , drop = FALSE])
    mat <- matrix(counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) { # nocov start
      mat <- .apply_scaling(mat, scaling,
                            include_zeros = .scaling_includes_zeros(net_method))
    } # nocov end
    if (thresh > 0) mat[abs(mat) < thresh] <- 0 # nocov
    dimnames(mat) <- list(states, states)
    mat
  }

  # ---- Build matrix from subset (association path) ----
  build_matrix_association <- function(idx) {
    sub_data <- data[idx, , drop = FALSE]
    if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
      sub_data <- tryCatch( # nocov start
        .decompose_multilevel(sub_data, id_col = id_col, level = level),
        error = function(e) NULL
      )
      if (is.null(sub_data)) return(NULL) # nocov end
    }
    est <- tryCatch(
      do.call(estimator$fn, c(list(data = sub_data), params)),
      error = function(e) NULL
    )
    if (is.null(est)) return(NULL)
    # a state lost to the subset (a column constant in it) leaves no network to
    # correlate: an unbuilt draw, which counts as below the threshold
    if (!all(states %in% rownames(est$matrix))) return(NULL)
    mat <- est$matrix[states, states] # nocov start
    if (!is.null(scaling)) {
      mat <- .apply_scaling(mat, scaling,
                            include_zeros = .scaling_includes_zeros(net_method))
    }
    if (thresh > 0) mat[abs(mat) < thresh] <- 0
    mat # nocov end
  }

  build_matrix <- if (is_transition) build_matrix_transition
                  else build_matrix_association

  # ---- Run stability assessment ----
  n_prop <- length(drop_prop)
  n_measures <- length(measures)

  # Pre-allocate correlation storage
  corr_storage <- lapply(measures, function(m) {
    matrix(NA_real_, nrow = iter, ncol = n_prop)
  })
  names(corr_storage) <- measures

  failed <- new.env(parent = emptyenv())
  failed$n <- 0L
  failed$total <- 0L
  for (p_idx in seq_len(n_prop)) {
    prop <- drop_prop[p_idx]
    n_keep <- n_seq - floor(n_seq * prop)
    if (n_keep < 2L) next # nocov

    # vapply: each iteration returns correlation vector (one per measure)
    corr_mat <- vapply(seq_len(iter), function(i) {
      idx <- sample.int(n_seq, n_keep, replace = FALSE)
      mat <- build_matrix(idx)
      failed$total <- failed$total + 1L
      if (is.null(mat)) {
        failed$n <- failed$n + 1L
        return(rep(NA_real_, n_measures))
      }

      sub_cents <- .compute_centralities(
        mat, states, x$directed, measures,
        loops = loops, normalize = normalize, invert = invert,
        normalize_diffusion = normalize_diffusion, centrality_fn = centrality_fn
      )

      vapply(measures, function(m) {
        sc <- sub_cents[[m]]
        oc <- orig_cents[[m]]
        ok <- is.finite(sc) & is.finite(oc)
        # a resampled network can leave a measure constant or NA (e.g.
        # Diffusion / RSP betweenness on a degenerate draw): no correlation
        if (sum(ok) < 2L || .is_constant_vector(sc[ok])) return(NA_real_)
        cor(sc[ok], oc[ok], method = method)
      }, numeric(1))
    }, numeric(n_measures))

    # corr_mat is n_measures x iter, store each measure's column
    if (n_measures == 1L) {
      corr_storage[[1]][, p_idx] <- corr_mat
    } else {
      for (k in seq_len(n_measures)) {
        corr_storage[[measures[k]]][, p_idx] <- corr_mat[k, ]
      }
    }
  }

  .warn_failed_resamples(failed$n, failed$total, "centrality_stability()",
                         below_threshold = TRUE)

  # ---- Compute CS coefficients ----
  cs <- vapply(measures, function(m) {
    .calculate_cs(corr_storage[[m]], threshold, certainty, drop_prop)
  }, numeric(1))

  result <- list(
    cs = cs,
    correlations = corr_storage,
    measures = measures,
    drop_prop = drop_prop,
    threshold = threshold,
    certainty = certainty,
    iter = iter,
    method = method
  )
  class(result) <- "net_stability"
  result
}


# ---- Helpers ----

#' Compute centralities from a weight matrix
#'
#' Built-in measures are dependency-free and follow tna centrality semantics.
#'
#' @param mat Weight matrix.
#' @param states Character vector of state names.
#' @param directed Logical.
#' @param measures Character vector of requested measures.
#' @param loops Logical. If FALSE, zero out diagonal before computing.
#' @param normalize Logical. Range-normalize all requested built-in measures.
#' @param invert Logical. Invert weights for path-based measures.
#' @param normalize_diffusion Logical. Range-normalize Diffusion even when
#'   \code{normalize = FALSE}.
#' @param centrality_fn Optional function taking a weight matrix and
#'   returning a named list of centrality vectors.
#' @noRd
.compute_centralities <- function(mat, states, directed, measures,
                                  loops = FALSE, normalize = FALSE,
                                  invert = TRUE,
                                  normalize_diffusion = TRUE,
                                  centrality_fn = NULL) {
  n <- length(states)
  dimnames(mat) <- list(states, states)
  if (!loops) diag(mat) <- 0
  result <- list()

  builtin <- .centrality_builtin_measures()
  canonical <- vapply(measures, .centrality_canonical_measure, character(1))

  if (any(c("InStrength", "OutStrength") %in% canonical)) {
    if ("InStrength" %in% canonical) result[["InStrength"]] <- colSums(mat)
    if ("OutStrength" %in% canonical) result[["OutStrength"]] <- rowSums(mat)
  }

  abs_mat <- abs(mat)
  if (any(c("Betweenness", "ClosenessIn", "ClosenessOut", "Closeness") %in%
          canonical)) {
    if ("Betweenness" %in% canonical) {
      result[["Betweenness"]] <- .betweenness(abs_mat, directed = directed,
                                               invert = invert)
    }
    if (any(c("ClosenessIn", "ClosenessOut", "Closeness") %in% canonical)) {
      cl <- .closeness(abs_mat, directed = directed, invert = invert)
      for (m in intersect(c("ClosenessIn", "ClosenessOut", "Closeness"),
                          canonical)) {
        result[[m]] <- cl[[m]]
      }
    }
  }

  if ("BetweennessRSP" %in% canonical) {
    result[["BetweennessRSP"]] <- .rsp_betweenness(abs_mat)
  }
  if ("Diffusion" %in% canonical) {
    result[["Diffusion"]] <- .diffusion_centrality(mat)
  }
  if ("Clustering" %in% canonical) {
    result[["Clustering"]] <- .weighted_clustering(mat + t(mat))
  }

  # External centralities via centrality_fn (user-supplied measures only)
  external <- measures[!(measures %in% builtin)]
  if (length(external) > 0L) {
    if (is.null(centrality_fn)) {
      stop(errorCondition(paste0("centrality_fn is required for measures: ",
           paste(external, collapse = ", "),
           ". Provide a function that takes a weight matrix and returns ",
           "a named list of centrality vectors."),
                          class = "nestimate_bad_input", call = NULL))
    }
    custom <- centrality_fn(mat)
    if (!is.list(custom)) { # nocov start
      stop(errorCondition("centrality_fn must return a named list.",
                          class = "nestimate_bad_input", call = NULL))
    } # nocov end
    for (m in external) {
      if (is.null(custom[[m]])) { # nocov start
        stop(errorCondition(
          paste0("centrality_fn did not return measure '", m, "'."),
          class = "nestimate_bad_input", call = NULL))
      } # nocov end
      result[[m]] <- custom[[m]]
    }
  }

  for (m in names(result)) {
    if (isTRUE(normalize) ||
        (identical(m, "Diffusion") && isTRUE(normalize_diffusion))) {
      result[[m]] <- .range01_tna(result[[m]])
    }
  }

  out <- vector("list", length(measures))
  names(out) <- measures
  for (i in seq_along(measures)) {
    m <- measures[[i]]
    cm <- canonical[[i]]
    if (!is.null(result[[m]])) {
      out[[m]] <- result[[m]]
    } else {
      out[[m]] <- result[[cm]]
    }
    names(out[[m]]) <- states
  }
  out
}


#' Whether each resample's correlation reaches the threshold
#'
#' An undefined correlation (the values are constant in that resample) or a
#' resample that could not be built does not reach it, so it stays in the
#' denominator of the certainty rather than being dropped from it.
#' @noRd
.reaches_threshold <- function(correlations, threshold) {
  !is.na(correlations) & correlations >= threshold
}

#' Calculate CS-coefficient from correlation matrix
#' @noRd
.calculate_cs <- function(corr_mat, threshold, certainty, drop_prop) {
  prop_above <- colMeans(.reaches_threshold(corr_mat, threshold))
  # CS = the largest drop proportion that passes (Epskamp, Borsboom & Fried
  # 2018): the largest value, not the last position, so an unsorted
  # `drop_prop` gives the same CS as a sorted one.
  max(c(0, drop_prop[prop_above >= certainty]))
}

#' Is a vector constant up to rounding?
#'
#' A measure whose finite values span no more than sqrt(machine epsilon)
#' relative to their size has no variance to correlate: OutStrength of a
#' relative network with self-loops is 1 for every state up to 1e-16, and
#' an exact `sd() > 0` test kept it and correlated rounding noise.
#' @noRd
.is_constant_vector <- function(v) {
  v <- v[is.finite(v)]
  if (length(v) < 2L) return(TRUE)
  diff(range(v)) <= sqrt(.Machine$double.eps) * max(1, abs(v))
}


# ---- S3 Methods ----

#' @rdname centrality_stability
#' @param ... For the `print()`, `summary()` and `plot()` methods: not
#'   used.
#' @param object For the `summary()` methods: a `net_stability` or
#'   `net_stability_group`.
#' @export
print.net_stability <- function(x, ...) {
  cat(sprintf("Centrality Stability (%d iterations, threshold = %.1f)\n",
              x$iter, x$threshold))
  cat(sprintf("  Drop proportions: %s\n",
              paste(x$drop_prop, collapse = ", ")))
  cat("\n  CS-coefficients:\n")
  for (m in names(x$cs)) {
    cat(sprintf("    %-15s  %.2f\n", m, x$cs[m]))
  }
  invisible(x)
}

#' @rdname centrality_stability
#' @export
print.net_stability_group <- function(x, ...) {
  measures <- unique(unlist(lapply(x, function(e) names(e$cs))))
  cs_mat <- sapply(x, function(e) e$cs[measures])
  rownames(cs_mat) <- measures
  cat(sprintf("Centrality Stability (%d networks, threshold = %.2f)\n",
              length(x), x[[1]]$threshold))
  print(round(t(cs_mat), 2))
  invisible(x)
}


#' @rdname centrality_stability
#' @export
summary.net_stability_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df <- summary(object[[nm]])
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' @rdname centrality_stability
#' @export
summary.net_stability <- function(object, ...) {
  rows <- do.call(rbind, lapply(object$measures, function(m) {
    corr_mat <- object$correlations[[m]]
    do.call(rbind, lapply(seq_along(object$drop_prop), function(j) {
      vals <- corr_mat[, j]
      data.frame(
        measure = m,
        drop_prop = object$drop_prop[j],
        mean_cor = mean(vals, na.rm = TRUE),
        sd_cor = sd(vals, na.rm = TRUE),
        prop_above = mean(.reaches_threshold(vals, object$threshold)),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))
  rownames(rows) <- NULL
  rows
}


#' @rdname centrality_stability
#' @export
plot.net_stability <- function(x, ...) {
  summ <- summary(x)

  p <- ggplot2::ggplot(summ, ggplot2::aes(
    x = .data$drop_prop, y = .data$mean_cor,
    color = .data$measure)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$mean_cor - .data$sd_cor,
        ymax = pmin(.data$mean_cor + .data$sd_cor, 1),
        fill = .data$measure),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_hline(yintercept = x$threshold,
                        linetype = "dashed", color = "grey40") +
    ggplot2::annotate("text", x = max(x$drop_prop), y = x$threshold,
                      label = sprintf("threshold = %.1f", x$threshold),
                      hjust = 1, vjust = -0.5, size = 3, color = "grey40") +
    ggplot2::scale_x_continuous(
      breaks = x$drop_prop,
      labels = x$drop_prop
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Proportion dropped",
      y = sprintf("Mean correlation (%s)", x$method),
      title = "Centrality Stability",
      color = "Measure", fill = "Measure"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}
