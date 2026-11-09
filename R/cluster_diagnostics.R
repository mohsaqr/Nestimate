# ==============================================================================
# cluster_diagnostics() -- unified clustering quality surface.
#
# Single entry point that returns a structured `net_cluster_diagnostics`
# object regardless of whether the clustering came from build_clusters()
# (distance), build_mmm() (mmm), or one of the cluster_*() wrappers that
# return a netobject_group with an attr(, "clustering") of class
# `net_clustering` or `net_mmm_clustering`.
#
# This module only PACKAGES information that's already computed:
#  - .mmm_quality()    in R/mmm.R         -- avepp, entropy, classification_error
#  - cluster::silhouette() (already used inside plot.net_clustering)
#  - cl$silhouette / mmm$BIC / etc.       -- read straight off the source
#
# Helpers it reuses (do not reimplement):
#  - .cluster_table_lines() / .fmt_size_pct() in R/cluster_data.R
#  - .plot_mmm_posterior() / .plot_covariate_forest() if a method needs them
# ==============================================================================

#' Cluster Diagnostics
#'
#' @description
#' `cluster_diagnostics()` collects the quality measures of a clustering in
#' one `net_cluster_diagnostics` object, laid out the same way whichever
#' route produced the clusters: a distance clustering from
#' [build_clusters()], a mixed Markov model from [build_mmm()] or
#' [cluster_mmm()], the `netobject_group` that [cluster_network()] or
#' [build_network()] builds from one of them, or the `net_mmm_clustering`
#' attribute of such a group. It packages what the fit already computed and
#' adds only the per-cluster summaries the fit does not store (silhouette
#' means and mean within-cluster distances for a distance clustering, the
#' classification error within each cluster for a mixture); it does not
#' refit.
#'
#' For a distance clustering it reports the average silhouette width,
#' overall and per cluster, and each cluster's mean within-cluster
#' distance. For a mixture it reports each cluster's average posterior
#' probability (AvePP) and mixing proportion, the overall AvePP, the
#' normalised entropy of the posterior probabilities (0 when every sequence
#' is assigned with certainty, 1 when all assignments are equally likely;
#' the reverse of the relative entropy often reported for latent class
#' models), the classification error (the expected share of sequences
#' misclassified when each is assigned to its most likely cluster: one minus
#' the largest posterior, averaged over all sequences and over each
#' cluster's members, so it equals one minus the AvePP; Vermunt & Magidson,
#' 2005), and the log-likelihood with BIC, AIC and ICL.
#'
#' These measures describe how cleanly one solution separates the
#' sequences; none of them picks the number of clusters. Compare solutions
#' with [cluster_choice()] for distance clustering or [compare_mmm()] for
#' mixtures, and build the clusters' networks with [build_network()] on the
#' clustering to see how they differ, with [bootstrap_network()] for how
#' stable their edges are.
#'
#' @param x The clustering: a `net_clustering`, a `net_mmm`, a
#'   `netobject_group` carrying a clustering (built by [cluster_network()]
#'   or by [build_network()] from a `net_clustering` or `net_mmm`), or a
#'   `net_mmm_clustering`. A `netobject_group` without a clustering, and any
#'   other object, raise an error. For the `print()` and `plot()`
#'   methods: a `net_cluster_diagnostics`.
#' @param ... Not used: an extra argument raises an error, here and in the
#'   `print()` and `summary()` methods. For `plot()`: further
#'   arguments for the clustering's own `plot()` method, such as `combined`.
#' @return A `net_cluster_diagnostics` object. `print()` shows the numbers
#'   of sequences and clusters (with the method and dissimilarity of a
#'   distance clustering, or the number of states of a mixture), the overall
#'   quality measures, the information criteria of a mixture, one row per
#'   cluster, and the covariates of the fit, if any; `plot()` draws the clustering's own views. The object's
#'   fields:
#' \describe{
#'   \item{family}{\code{"distance"} or \code{"mmm"}.}
#'   \item{k, n, sizes}{Number of clusters, number of sequences, and the
#'     cluster sizes.}
#'   \item{per_cluster}{Data frame with one row per cluster; its columns
#'     depend on the family. Distance: \code{cluster}, \code{size},
#'     \code{pct}, \code{mean_within_dist}, \code{sil_mean}. Mixture:
#'     \code{cluster}, \code{size}, \code{pct}, \code{mix_pct},
#'     \code{avepp}, \code{class_err_pct}. \code{pct}, \code{mix_pct}
#'     and \code{class_err_pct} are in percent.}
#'   \item{overall}{List of the family's summary measures:
#'     \code{silhouette} for distance; \code{avepp_overall},
#'     \code{entropy} and \code{classification_error} (a proportion) for a
#'     mixture.}
#'   \item{ics}{For a mixture, a list with \code{BIC}, \code{AIC},
#'     \code{ICL} and \code{log_likelihood}; NULL for distance.}
#'   \item{metadata}{The settings of the fit. Distance: \code{method},
#'     \code{dissimilarity}, \code{weighted}, \code{lambda},
#'     \code{medoids}, \code{covariates}. Mixture: \code{states},
#'     \code{n_states}, \code{iterations}, \code{converged},
#'     \code{covariates}.}
#'   \item{source}{The clustering object itself, which \code{plot()} draws
#'     without recomputing anything.}
#' }
#' @references
#' Vermunt, J. K., & Magidson, J. (2005). \emph{Latent GOLD 4.0 User's
#' Guide}. Statistical Innovations.
#'
#' @seealso \code{\link{build_clusters}}, \code{\link{build_mmm}},
#'   \code{\link{cluster_choice}},
#'   \code{\link{compare_mmm}} for k-sweep model selection (MMM only).
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' cl <- build_clusters(seqs, k = 2, method = "ward.D2")
#' cluster_diagnostics(cl)
#' \donttest{
#' fit <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
#' cluster_diagnostics(fit)
#' summary(cluster_diagnostics(fit))
#' }
#' @export
cluster_diagnostics <- function(x, ...) {
  UseMethod("cluster_diagnostics")
}

#' @export
cluster_diagnostics.default <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  stop(errorCondition(paste0("cluster_diagnostics() has no method for class '",
       paste(class(x), collapse = "/"), "'. Supported inputs: ",
       "net_clustering (from build_clusters), net_mmm (from build_mmm), ",
       "netobject_group (from cluster_network / build_network), ",
       "or net_mmm_clustering (the clustering attribute itself)."),
       class = "nestimate_bad_input", call = NULL))
}

# ---------------------------------------------------------------------------
# Distance-based: net_clustering
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_clustering <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  k <- as.integer(x$k)
  assignments <- as.integer(x$assignments)
  sizes <- as.integer(x$sizes %||% tabulate(assignments, nbins = k))
  n_total <- sum(sizes)

  mean_within <- if (!is.null(x$distance)) {
    .per_cluster_within_dist(x$distance, assignments, k)
  } else {
    rep(NA_real_, k)
  }

  # Per-cluster silhouette mean -- compute once, aggregate. Skip when the
  # cluster package isn't available (it's a Suggests-style dep here, used
  # only inside plot.net_clustering already).
  sil_mean <- rep(NA_real_, k)
  overall_sil <- as.numeric(x$silhouette %||% NA_real_)
  if (!is.null(x$distance) && requireNamespace("cluster", quietly = TRUE)) {
    # A silhouette that cannot be computed leaves the columns NA, as before,
    # and now says why.
    sil <- tryCatch(
      cluster::silhouette(assignments, dist = x$distance),
      error = function(e) {
        warning(warningCondition(
          paste0("cluster_diagnostics(): the per-cluster silhouette could ",
                 "not be computed (", conditionMessage(e), "); it is NA."),
          class = "nestimate_no_silhouette", call = NULL))
        NULL
      }
    )
    # silhouette() returns NA, not a matrix, for a single cluster.
    if (is.matrix(sil) && nrow(sil) > 0L) {
      widths <- sil[, 3L]
      means  <- tapply(widths, sil[, 1L], mean)
      found  <- intersect(as.character(seq_len(k)), names(means))
      sil_mean[as.integer(found)] <- as.numeric(means[found])
      if (is.na(overall_sil)) overall_sil <- mean(widths)
    }
  }

  per_cluster <- data.frame(
    cluster          = seq_len(k),
    size             = sizes,
    pct              = if (n_total > 0L) sizes / n_total * 100 else
                        rep(0, k),
    mean_within_dist = mean_within,
    sil_mean         = sil_mean,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      family      = "distance",
      k           = k,
      n           = as.integer(n_total),
      sizes       = sizes,
      per_cluster = per_cluster,
      overall     = list(silhouette = overall_sil),
      ics         = NULL,
      metadata    = list(
        method        = x$method,
        dissimilarity = x$dissimilarity,
        weighted      = isTRUE(x$weighted),
        lambda        = x$lambda,
        medoids       = x$medoids,
        covariates    = if (is.null(x$covariates)) NULL else
                          setdiff(unique(x$covariates$coefficients$variable),
                                  "(Intercept)")
      ),
      source      = x
    ),
    class = "net_cluster_diagnostics"
  )
}

# ---------------------------------------------------------------------------
# Model-based: net_mmm
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_mmm <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  k <- as.integer(x$k)
  assignments <- as.integer(x$assignments)
  sizes <- as.integer(tabulate(assignments, nbins = k))
  n_total <- as.integer(x$n_sequences %||% sum(sizes))

  # Per-cluster classification error: the expected misclassification
  # 1 - max posterior averaged over the sequences assigned to cluster m
  # (Vermunt & Magidson, 2005), the breakdown of quality$classification_error.
  class_err_pct <- rep(NA_real_, k)
  if (!is.null(x$posterior)) {
    max_post <- apply(x$posterior, 1L, max)
    class_err_pct <- vapply(seq_len(k), \(cl) {
      members <- which(assignments == cl)
      if (length(members) == 0L) NA_real_ else mean(1 - max_post[members]) * 100
    }, numeric(1))
  }

  avepp <- as.numeric(x$quality$avepp %||% rep(NA_real_, k))
  mix   <- as.numeric(x$mixing %||% rep(NA_real_, k))

  per_cluster <- data.frame(
    cluster       = seq_len(k),
    size          = sizes,
    pct           = if (n_total > 0L) sizes / n_total * 100 else rep(0, k),
    mix_pct       = mix * 100,
    avepp         = avepp,
    class_err_pct = class_err_pct,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      family      = "mmm",
      k           = k,
      n           = n_total,
      sizes       = sizes,
      per_cluster = per_cluster,
      overall     = list(
        avepp_overall        = as.numeric(x$quality$avepp_overall),
        entropy              = as.numeric(x$quality$entropy),
        classification_error = as.numeric(x$quality$classification_error)
      ),
      ics         = list(BIC = as.numeric(x$BIC),
                         AIC = as.numeric(x$AIC),
                         ICL = as.numeric(x$ICL),
                         log_likelihood = as.numeric(x$log_likelihood)),
      metadata    = list(
        states     = x$states,
        n_states   = length(x$states),
        iterations = x$iterations,
        converged  = x$converged,
        covariates = if (is.null(x$covariates)) NULL else
                       setdiff(unique(x$covariates$coefficients$variable),
                               "(Intercept)")
      ),
      source      = x
    ),
    class = "net_cluster_diagnostics"
  )
}

# ---------------------------------------------------------------------------
# net_mmm_clustering -- the slim attribute object on a netobject_group.
# Same field shape as net_mmm minus $models, plus $data. We reuse the
# net_mmm method by re-classing in place (no copy of the heavy slots).
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_mmm_clustering <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  shadow <- structure(unclass(x), class = "net_mmm")
  out <- cluster_diagnostics.net_mmm(shadow)
  # Keep the source pointing at the actual input class so plot delegation
  # routes through plot.net_mmm_clustering rather than plot.net_mmm.
  out$source <- x
  out
}

# ---------------------------------------------------------------------------
# netobject_group -- delegate to attr(, "clustering")
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.netobject_group <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  cl <- attr(x, "clustering")
  if (is.null(cl)) {
    stop(errorCondition(
      paste0(
        "cluster_diagnostics() requires a clustering attribute on the ",
        "netobject_group. Build with cluster_network(), or with ",
        "build_network() on a build_clusters() or cluster_mmm() fit ",
        "(or attach an existing net_clustering / net_mmm_clustering as ",
        "attr(grp, \"clustering\"))."),
      class = "nestimate_bad_input", call = NULL))
  }
  cluster_diagnostics(cl)
}

.cluster_diag_check_unused_dots <- function(...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(errorCondition(paste0(
    "cluster_diagnostics() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", ")),
    class = "nestimate_bad_input", call = NULL)
  )
}

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

#' @rdname cluster_diagnostics
#' @param digits For `print()`: decimal places of the silhouette widths,
#'   distances, AvePP, entropy and information criteria, a non-negative
#'   whole number. Default `3`. Percentages are always shown to one
#'   decimal.
#' @param type For `plot()`: the view, passed to the clustering's own
#'   `plot()` method. Default `NULL`, which draws `"silhouette"` for a
#'   distance clustering and `"posterior"` for a mixture. A distance
#'   clustering also takes `"mds"`, `"heatmap"` and `"predictors"` (see
#'   [build_clusters()]); a mixture takes `"covariates"`, with
#'   `"predictors"` accepted for it when the source is a
#'   `net_mmm_clustering` (see [build_mmm()]).
#' @param object For the `summary()` method: a `net_cluster_diagnostics`.
#' @return For `print()`: the input, invisibly. For `plot()`: what the
#'   clustering's own `plot()` method returns, a `ggplot` object (invisibly)
#'   or, for the covariate views with `combined = FALSE`, a list of
#'   `ggplot` objects named by cluster (invisibly). For `summary()`:
#'   the per-cluster table described under `per_cluster`, a data frame with
#'   one row per cluster.
#' @export
print.net_cluster_diagnostics <- function(x, digits = 3L, ...) {
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop(errorCondition(
      paste0(
        "print.net_cluster_diagnostics() got unsupported argument",
        if (length(dots) == 1L) ": " else "s: ",
        paste(dot_names, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  .net_clustering_check_digits(digits)
  digits <- as.integer(digits)
  k <- x$k
  n <- x$n

  if (x$family == "distance") {
    cat(sprintf("Cluster Diagnostics (distance) [%s / %s]\n",
                x$metadata$method %||% "?", x$metadata$dissimilarity %||% "?"))
    cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n, k))
    if (!is.null(x$overall$silhouette) && !is.na(x$overall$silhouette)) {
      cat(sprintf("  Quality: silhouette = %.*f\n",
                  digits, x$overall$silhouette))
    }
  } else {
    cat(sprintf("Cluster Diagnostics (mmm) [k = %d]\n", k))
    cat(sprintf("  Sequences: %d  |  Clusters: %d  |  States: %d\n",
                n, k, x$metadata$n_states %||% NA_integer_))
    if (!is.null(x$overall$avepp_overall)) {
      cat(sprintf(
        "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
        digits, x$overall$avepp_overall,
        digits, x$overall$entropy,
        x$overall$classification_error * 100))
    }
    if (!is.null(x$ics)) {
      cat(sprintf(
        "  ICs: LL = %.*f  |  BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
        digits, x$ics$log_likelihood, digits, x$ics$BIC,
        digits, x$ics$AIC,            digits, x$ics$ICL))
    }
  }

  cat("\n")
  pc <- x$per_cluster
  if (x$family == "distance") {
    cols <- list(
      Cluster            = sprintf("%d", pc$cluster),
      N                  = .fmt_size_pct(pc$size, n),
      `Mean within-dist` = ifelse(is.na(pc$mean_within_dist), "--",
                                   sprintf(paste0("%.", digits, "f"),
                                           pc$mean_within_dist)),
      Silhouette         = ifelse(is.na(pc$sil_mean), "--",
                                   sprintf(paste0("%.", digits, "f"),
                                           pc$sil_mean))
    )
  } else {
    cols <- list(
      Cluster      = sprintf("%d", pc$cluster),
      N            = .fmt_size_pct(pc$size, n),
      `Mix%`       = ifelse(is.na(pc$mix_pct), "--",
                             sprintf("%4.1f%%", pc$mix_pct)),
      AvePP        = ifelse(is.na(pc$avepp), "--",
                             sprintf(paste0("%.", digits, "f"), pc$avepp)),
      `Class.Err%` = ifelse(is.na(pc$class_err_pct), "--",
                             sprintf("%4.1f%%", pc$class_err_pct))
    )
  }
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  # Read from the fitted object, not metadata$covariates: that field lists
  # the logit's design columns and is empty after a chi-square analysis.
  cov <- x$source$covariates
  if (!is.null(cov)) {
    label <- if (x$family == "distance") "post-hoc" else
      .mmm_covariate_label(x$source$covariate_effect)
    cat(.covariate_line(cov, label))
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# plot method -- delegates to the source's plot method
# ---------------------------------------------------------------------------

#' @rdname cluster_diagnostics
#' @export
plot.net_cluster_diagnostics <- function(x, type = NULL, ...) {
  if (is.null(type)) {
    type <- if (x$family == "distance") "silhouette" else "posterior"
  }
  plot(x$source, type = type, ...)
}

# ---------------------------------------------------------------------------
# summary method: the per-cluster table
# ---------------------------------------------------------------------------

#' @rdname cluster_diagnostics
#' @export
summary.net_cluster_diagnostics <- function(object, ...) {
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots) %||% rep("", length(dots))
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop(errorCondition(
      paste0(
        "summary.net_cluster_diagnostics() got unsupported argument",
        if (length(dots) == 1L) ": " else "s: ",
        paste(dot_names, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
  }
  object$per_cluster
}

# Kept so that code written against earlier versions still runs; the
# documented verb is summary().
#' @method as.data.frame net_cluster_diagnostics
#' @export
#' @noRd
as.data.frame.net_cluster_diagnostics <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
  summary(x, ...)
}
