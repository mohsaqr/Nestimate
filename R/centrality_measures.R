# ---- Built-in centrality measures (no external dependencies) ----
#
# All path-based measures (Betweenness, Closeness) are derived from
# all-pairs shortest paths via Floyd-Warshall. For transition/frequency
# networks, weights are inverted so higher weight = shorter distance.

#' All-pairs shortest paths via Floyd-Warshall (vectorized)
#'
#' @param W Square numeric weight matrix (zeros = no edge).
#' @param invert Logical. Convert weights to distances by 1/w? Default TRUE.
#' @return List with \code{D} (distance matrix) and \code{sigma} (shortest
#'   path count matrix).
#' @noRd
.floyd_warshall_sp <- function(W, invert = TRUE) {
  n   <- nrow(W)
  pos <- W > 0

  D              <- matrix(Inf, n, n)
  D[pos]         <- if (invert) 1 / W[pos] else W[pos]
  diag(D)        <- 0

  sigma          <- matrix(0L, n, n)
  sigma[pos]     <- 1L
  diag(sigma)    <- 1L

  Reduce(function(s, k) {
    D     <- s$D
    sigma <- s$sigma
    new_d <- outer(D[, k], D[k, ], "+")
    new_s <- outer(sigma[, k], sigma[k, ], "*")
    # Standard Floyd-Warshall only relaxes (s, t) when both s and t differ
    # from k. With diag(D) = 0, new_d[k, j] would otherwise match D[k, j]
    # under the equal-path test and double-count sigma.
    shorter        <- new_d < D & is.finite(new_d)
    equal          <- (new_d == D) & is.finite(new_d) & new_d > 0
    shorter[k, ]   <- FALSE; shorter[, k] <- FALSE
    equal[k, ]     <- FALSE; equal[, k]   <- FALSE
    sigma[shorter] <- new_s[shorter]
    sigma[equal]   <- sigma[equal] + new_s[equal]
    new_D          <- D
    new_D[shorter] <- new_d[shorter]
    list(D = new_D, sigma = sigma)
  }, seq_len(n), list(D = D, sigma = sigma))
}


#' Betweenness centrality (directed or undirected, weighted)
#'
#' Number of all shortest paths passing through each node. This matches
#' \code{igraph::betweenness(normalized = FALSE)}, which is the scale used by
#' \code{tna::centralities()}.
#'
#' @param W Square numeric weight matrix.
#' @param directed Logical.
#' @param invert Logical. Invert weights to distances? Default TRUE.
#' @return Named numeric vector of betweenness values.
#' @noRd
.betweenness <- function(W, directed = TRUE, invert = TRUE) {
  n  <- nrow(W)
  if (n < 3L) return(setNames(numeric(n), rownames(W)))
  sp <- .floyd_warshall_sp(W, invert)
  D  <- sp$D
  sg <- sp$sigma

  btw <- vapply(seq_len(n), function(v) {
    idx    <- seq_len(n)[-v]
    d_sv   <- D[idx, v]
    d_vt   <- D[v, idx]
    d_st   <- D[idx, idx]
    sg_sv  <- sg[idx, v]
    sg_vt  <- sg[v, idx]
    sg_st  <- sg[idx, idx]
    d_svt  <- outer(d_sv, d_vt, "+")
    on_path <- is.finite(d_st) & sg_st > 0L &
               abs(d_svt - d_st) < 1e-10
    diag(on_path) <- FALSE
    num  <- outer(sg_sv, sg_vt, "*")
    sum((num / sg_st)[on_path], na.rm = TRUE)
  }, numeric(1))

  if (!directed) btw <- btw / 2
  setNames(btw, rownames(W))
}


#' Edge betweenness (directed or undirected, weighted)
#'
#' Replaces each edge's weight with the number of shortest paths that
#' traverse it (fractional when shortest paths tie), matching
#' \code{tna::betweenness_network()} / \code{igraph::edge_betweenness()}.
#' For a probability/transition network, weights are inverted to distances
#' (\code{invert = TRUE}) so the geodesic between two states is the most
#' probable route, not the fewest hops.
#'
#' @param W Square numeric weight matrix (zeros = no edge).
#' @param invert Logical. Invert weights to distances by \code{1/w}?
#'   Default \code{TRUE}.
#' @return A numeric matrix the same shape as \code{W}: edge-betweenness on
#'   edge positions, zero elsewhere. Symmetric when \code{W} is symmetric.
#' @noRd
.edge_betweenness <- function(W, invert = TRUE) {
  n   <- nrow(W)
  EB  <- matrix(0, n, n, dimnames = dimnames(W))
  if (n < 2L) return(EB)

  sp  <- .floyd_warshall_sp(W, invert)
  D   <- sp$D
  sg  <- sp$sigma
  pos <- W > 0
  len <- matrix(Inf, n, n)
  len[pos] <- if (invert) 1 / W[pos] else W[pos]

  edges <- which(pos, arr.ind = TRUE)
  vals  <- vapply(seq_len(nrow(edges)), function(e) {
    a <- edges[e, 1L]
    b <- edges[e, 2L]
    # A shortest s->t path uses edge a->b iff
    #   d(s, a) + len(a, b) + d(b, t) == d(s, t).
    # Such paths number sigma(s, a) * sigma(b, t); divide by sigma(s, t).
    through <- outer(D[, a], D[b, ], "+") + len[a, b]
    on_path <- is.finite(D) & sg > 0L & abs(through - D) < 1e-9
    diag(on_path) <- FALSE                       # exclude s == t
    sum((outer(sg[, a], sg[b, ]) / sg)[on_path])
  }, numeric(1))

  EB[edges] <- vals
  EB
}


#' Closeness centrality (directed or undirected, weighted)
#'
#' Returns \code{ClosenessIn}, \code{ClosenessOut}, and \code{Closeness} on the
#' unnormalized \code{igraph::closeness()} scale used by \code{tna}.
#'
#' @param W Square numeric weight matrix.
#' @param directed Logical.
#' @param invert Logical. Invert weights to distances? Default TRUE.
#' @return For directed: named list with \code{InCloseness} and
#'   \code{OutCloseness} vectors. For undirected: named list with
#'   \code{Closeness} vector.
#' @noRd
.closeness <- function(W, directed = TRUE, invert = TRUE) {
  n   <- nrow(W)
  nms <- rownames(W)
  D   <- .floyd_warshall_sp(W, invert)$D
  W_all <- pmax(W, t(W))
  D_all <- .floyd_warshall_sp(W_all, invert)$D

  .cl <- function(d_vec) {
    finite_d <- d_vec[is.finite(d_vec) & d_vec > 0]
    if (length(finite_d) == 0L) NaN else 1 / sum(finite_d)
  }

  list(
    ClosenessIn = setNames(vapply(seq_len(n), function(v) .cl(D[, v]),
                                  numeric(1)), nms),
    ClosenessOut = setNames(vapply(seq_len(n), function(v) .cl(D[v, ]),
                                   numeric(1)), nms),
    Closeness = setNames(vapply(seq_len(n), function(v) .cl(D_all[v, ]),
                                numeric(1)), nms)
  )
}

#' Randomized shortest-path betweenness
#' @noRd
.rsp_betweenness <- function(mat, beta = 0.01) {
  n <- ncol(mat)
  d <- rowSums(mat)
  if (any(d == 0)) return(setNames(rep(NA_real_, n), rownames(mat)))
  p_ref <- diag(d^-1, n, n) %*% mat
  cost <- mat^-1
  cost[is.infinite(cost)] <- 0
  W <- p_ref * exp(-beta * cost)
  Z <- tryCatch(solve(diag(1, n, n) - W), error = function(e) NULL)
  if (is.null(Z)) return(setNames(rep(NA_real_, n), rownames(mat)))
  Z_recip <- Z^-1
  Z_recip[is.infinite(Z_recip)] <- 0
  Z_recip_diag <- diag(Z_recip) * diag(1, n, n)
  out <- diag(tcrossprod(Z, Z_recip - n * Z_recip_diag) %*% Z)
  out <- round(out)
  out <- out - min(out, na.rm = TRUE) + 1
  setNames(out, rownames(mat))
}

#' Diffusion centrality
#' @noRd
.diffusion_centrality <- function(mat) {
  n <- ncol(mat)
  p <- diag(1, n, n)
  s <- 0
  for (i in seq_len(n)) {
    p <- p %*% mat
    s <- s + p
  }
  setNames(rowSums(s), rownames(mat))
}

#' Weighted clustering coefficient used by tna
#' @noRd
.weighted_clustering <- function(mat) {
  diag(mat) <- 0
  n <- ncol(mat)
  num <- diag(mat %*% mat %*% mat)
  den <- colSums(mat)^2 - colSums(mat^2)
  setNames(num / den, rownames(mat))
}

#' Range-normalize on the same scale as tna::centralities(normalize = TRUE)
#' @noRd
.range01_tna <- function(x, na.rm = FALSE) {
  mi <- min(x, na.rm = na.rm)
  ma <- max(x, na.rm = na.rm)
  (x - mi) / (ma - mi)
}

#' @noRd
.centrality_builtin_measures <- function() {
  c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut",
    "Closeness", "Betweenness", "BetweennessRSP", "Diffusion",
    "Clustering", "InCloseness", "OutCloseness")
}

#' @noRd
.centrality_default_measures <- function() {
  c("InStrength", "Betweenness", "Diffusion")
}

#' @noRd
.centrality_all_measures <- function() {
  c("OutStrength", "InStrength", "ClosenessIn", "ClosenessOut",
    "Closeness", "Betweenness", "BetweennessRSP", "Diffusion",
    "Clustering")
}

#' @noRd
.resolve_measures <- function(measures) {
  if (is.null(measures)) return(.centrality_default_measures())
  if (length(measures) == 1L && identical(tolower(measures), "all")) {
    return(.centrality_all_measures())
  }
  measures
}

#' Format centrality values: blank zeros, integers without trailing zeros
#' @noRd
.fmt_centrality <- function(x, signed = FALSE) {
  out <- character(length(x))
  fin <- is.finite(x)
  zero <- fin & abs(x) < 1e-8            # drop 0 labels entirely
  whole <- fin & !zero & abs(x - round(x)) < 1e-8
  out[whole] <- sprintf(if (signed) "%+d" else "%d",
                        as.integer(round(x[whole])))
  frac <- fin & !zero & !whole
  out[frac] <- sprintf(if (signed) "%+.2f" else "%.2f", x[frac])
  out
}

#' @noRd
.centrality_canonical_measure <- function(measure) {
  switch(measure,
    InCloseness = "ClosenessIn",
    OutCloseness = "ClosenessOut",
    measure
  )
}

#' @noRd
.as_net_centrality_frame <- function(result, states, measures,
                                     normalized, invert,
                                     normalize_diffusion) {
  out <- data.frame(state = factor(states, levels = states),
                    check.names = FALSE)
  for (m in measures) out[[m]] <- unname(result[[m]])
  rownames(out) <- states
  attr(out, "measures") <- measures
  attr(out, "normalized") <- normalized
  attr(out, "invert") <- invert
  attr(out, "normalize_diffusion") <- normalize_diffusion
  class(out) <- c("net_centrality", class(out))
  out
}


# ---- Internal centrality() generic (S3 dispatch) ----

#' @noRd
centrality <- function(x, ...) {
  UseMethod("centrality")
}


# ---- Exported net_centrality() ----

#' Compute Centrality Measures for a Network
#'
#' Computes centrality measures from a \code{netobject},
#' \code{netobject_group}, \code{mcml}, or \code{cograph_network}. The built-in
#' measures match \code{tna::centralities()} without importing \code{tna} or
#' \code{igraph}. The only intentional default difference is that
#' \code{Diffusion} is range-normalized by default.
#'
#' @param x A \code{netobject}, \code{netobject_group}, or
#'   \code{cograph_network}.
#' @param measures Character vector. Centrality measures to compute.
#'   Defaults to \code{c("InStrength", "Betweenness", "Diffusion")}. Pass
#'   \code{"all"} for every built-in measure: \code{"OutStrength"},
#'   \code{"InStrength"}, \code{"ClosenessIn"}, \code{"ClosenessOut"},
#'   \code{"Closeness"}, \code{"Betweenness"}, \code{"BetweennessRSP"},
#'   \code{"Diffusion"}, and \code{"Clustering"}. The legacy aliases
#'   \code{"InCloseness"} and \code{"OutCloseness"} are also accepted.
#' @param loops Logical. Include self-loops (diagonal) in computation?
#'   Default: \code{FALSE}.
#' @param normalize Logical. Range-normalize all requested measures using the
#'   same transformation as \code{tna::centralities(normalize = TRUE)}.
#'   Default: \code{FALSE}.
#' @param invert Logical. Invert weights for shortest-path measures?
#'   Default: \code{TRUE}, matching \code{tna}.
#' @param normalize_diffusion Logical. Range-normalize \code{Diffusion} even
#'   when \code{normalize = FALSE}. Default: \code{TRUE}.
#' @param centrality_fn Optional function. Custom centrality function that
#'   takes a weight matrix and returns a named list of centrality vectors.
#' @param ... Additional arguments (ignored).
#'
#' @return For a \code{netobject}: a \code{net_centrality} data frame with
#'   node names as rows, a \code{state} column, and one column per centrality
#'   measure. For a \code{netobject_group}: a \code{net_centrality_group}
#'   list of such data frames.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' net_centrality(net)
#'
#' @export
net_centrality <- function(x, measures = NULL, loops = FALSE,
                            normalize = FALSE, invert = TRUE,
                            normalize_diffusion = TRUE,
                            centrality_fn = NULL, ...) {
  if (isFALSE(loops)) {
    message("centralities computed excluding loops (diagonal). ",
            "Pass `loops = TRUE` to include self-transitions.")
  }
  centrality(x, measures = measures, loops = loops,
             normalize = normalize, invert = invert,
             normalize_diffusion = normalize_diffusion,
             centrality_fn = centrality_fn, ...)
}


#' @noRd
centrality.netobject <- function(x, measures = NULL, loops = FALSE,
                                  normalize = FALSE, invert = TRUE,
                                  normalize_diffusion = TRUE,
                                  centrality_fn = NULL, ...) {
  mat      <- x$weights
  states   <- x$nodes$label
  directed <- x$directed

  measures <- .resolve_measures(measures)

  res <- .compute_centralities(mat, states, directed, measures,
                                loops = loops, normalize = normalize,
                                invert = invert,
                                normalize_diffusion = normalize_diffusion,
                                centrality_fn = centrality_fn)
  .as_net_centrality_frame(res, states, measures,
                           normalized = normalize, invert = invert,
                           normalize_diffusion = normalize_diffusion)
}


#' @noRd
centrality.netobject_group <- function(x, measures = NULL, loops = FALSE,
                                        normalize = FALSE, invert = TRUE,
                                        normalize_diffusion = TRUE,
                                        centrality_fn = NULL, ...) {
  out <- lapply(x, function(net) {
    centrality.netobject(net, measures = measures, loops = loops,
                         normalize = normalize, invert = invert,
                         normalize_diffusion = normalize_diffusion,
                         centrality_fn = centrality_fn)
  })
  names(out) <- names(x)
  attr(out, "measures") <- .resolve_measures(measures)
  class(out) <- c("net_centrality_group", "list")
  out
}


#' @noRd
centrality.cograph_network <- function(x, measures = NULL, loops = FALSE,
                                        normalize = FALSE, invert = TRUE,
                                        normalize_diffusion = TRUE,
                                        centrality_fn = NULL, ...) {
  centrality.netobject(.as_netobject(x), measures = measures, loops = loops,
                        normalize = normalize, invert = invert,
                        normalize_diffusion = normalize_diffusion,
                        centrality_fn = centrality_fn)
}


#' @noRd
centrality.mcml <- function(x, measures = NULL, loops = FALSE,
                             normalize = FALSE, invert = TRUE,
                             normalize_diffusion = TRUE,
                             centrality_fn = NULL, ...) {
  centrality.netobject_group(as_tna(x), measures = measures, loops = loops,
                              normalize = normalize, invert = invert,
                              normalize_diffusion = normalize_diffusion,
                              centrality_fn = centrality_fn)
}

#' Plot centrality measures
#'
#' @param x A \code{net_centrality} object returned by
#'   \code{\link{net_centrality}}.
#' @param reorder Logical. Reorder states within each centrality panel by
#'   centrality value. Default: \code{TRUE}.
#' @param type Plot type. \code{"bar"} shows one faceted horizontal bar chart
#'   per measure; \code{"line"} shows state profiles as lines across measures
#'   (the value \code{"profile"} is still accepted as an alias);
#'   \code{"heatmap"} shows a states-by-measures tile grid, each measure
#'   scaled to 0--1 for cross-measure comparability with the raw value
#'   printed in the tile.
#' @param ncol Integer. Number of facet columns.
#' @param scales Facet scale mode. \code{"free_x"} uses free centrality axes;
#'   \code{"fixed"} keeps a common centrality axis.
#' @param profile_scale Scaling used by \code{type = "line"}.
#'   \code{"measure"} rescales each centrality measure to 0--1 before drawing
#'   cross-measure profiles; \code{"none"} uses raw values.
#' @param labels Logical. Add compact value labels.
#' @param drop_zero Logical. Drop measures whose values are all (near) zero
#'   so empty panels do not waste space. Default: \code{FALSE} (every
#'   requested measure is shown).
#' @param ... Additional arguments ignored.
#' @return A \code{ggplot} object.
#' @export
plot.net_centrality <- function(x, reorder = TRUE, ncol = 3L,
                                type = c("bar", "line", "heatmap"),
                                scales = c("free_x", "fixed"),
                                profile_scale = c("measure", "none"),
                                labels = TRUE, drop_zero = FALSE, ...) {
  if (length(type) == 1L && identical(type, "profile")) type <- "line"
  type <- match.arg(type)
  scales <- match.arg(scales)
  profile_scale <- match.arg(profile_scale)
  if (type == "line") {
    profile_ncol <- if (missing(ncol)) NULL else ncol
    profile_labels <- if (missing(labels)) FALSE else labels
    return(.plot_net_centrality_profile(
      x, ncol = profile_ncol, profile_scale = profile_scale,
      labels = profile_labels
    ))
  }
  if (type == "heatmap") {
    return(.plot_net_centrality_heatmap(
      .centrality_plot_data(x), reorder = reorder, labels = labels
    ))
  }
  scales <- if (scales == "free_x") "free_x" else "fixed"
  d <- .centrality_plot_data(x)
  if (isTRUE(drop_zero)) d <- .centrality_drop_zero(d)
  d$state <- .centrality_state_factor(d, reorder = reorder)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value,
                                        y = .data$state,
                                        fill = .data$measure)) +
    ggplot2::geom_col(width = 0.72, alpha = 0.88, show.legend = FALSE,
                      na.rm = TRUE) +
    ggplot2::facet_wrap(~ measure, ncol = ncol, scales = scales) +
    ggplot2::scale_fill_manual(values = .centrality_measure_palette()) +
    ggplot2::labs(x = "Centrality", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      hjust = -0.08, size = 2.8, na.rm = TRUE
    ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.16))
      )
  }
  p
}

#' Plot grouped centrality measures
#'
#' @param x A \code{net_centrality_group} object returned by
#'   \code{\link{net_centrality}}.
#' @param reorder Logical. Reorder states by their mean value within each
#'   centrality panel. Default: \code{TRUE}.
#' @param type Plot type. \code{"bar"} shows grouped bars within each measure;
#'   \code{"line"} facets by state and draws one line per group across
#'   centrality measures (\code{"profile"} is accepted as an alias);
#'   \code{"delta"} draws a diverging bar of group
#'   differences. With two groups it is the per-state difference (second group
#'   minus first); with three or more groups it is each group's deviation from
#'   the per-state group mean, so the largest gaps stand out either way.
#' @param ncol Integer. Number of facet columns.
#' @param scales Facet scale mode. \code{"free_x"} uses free centrality axes;
#'   \code{"fixed"} keeps a common centrality axis.
#' @param palette Brewer palette for groups.
#' @param profile_scale Scaling used by \code{type = "line"}.
#'   \code{"measure"} rescales each centrality measure to 0--1 before drawing
#'   cross-measure profiles; \code{"none"} uses raw values.
#' @param labels Logical. Add compact value labels.
#' @param drop_zero Logical. Drop measures whose values are all (near) zero
#'   so empty panels do not waste space. Default: \code{FALSE}.
#' @param ... Additional arguments ignored.
#' @return A \code{ggplot} object.
#' @export
plot.net_centrality_group <- function(x, reorder = TRUE, ncol = 3L,
                                      type = c("bar", "line", "delta"),
                                      scales = c("free_x", "fixed"),
                                      palette = "Set2",
                                      profile_scale = c("measure", "none"),
                                      labels = FALSE, drop_zero = FALSE,
                                      ...) {
  if (length(type) == 1L && identical(type, "profile")) type <- "line"
  type <- match.arg(type)
  scales <- match.arg(scales)
  profile_scale <- match.arg(profile_scale)
  scales <- if (scales == "free_x") "free_x" else "fixed"
  d <- do.call(rbind, lapply(names(x), function(g) {
    z <- .centrality_plot_data(x[[g]])
    z$group <- g
    z
  }))
  rownames(d) <- NULL
  if (isTRUE(drop_zero)) d <- .centrality_drop_zero(d)
  if (type == "delta") {
    return(.plot_net_centrality_delta(d, ncol = ncol, reorder = reorder,
                                      labels = labels, palette = palette))
  }
  if (type == "line") {
    profile_ncol <- if (missing(ncol)) NULL else ncol
    return(.plot_net_centrality_group_profile(
      d, ncol = profile_ncol, palette = palette,
      profile_scale = profile_scale, labels = labels
    ))
  }
  d$state <- .centrality_state_factor(d, reorder = reorder)
  dodge <- ggplot2::position_dodge2(width = 0.76, preserve = "single",
                                    padding = 0.08)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value,
                                        y = .data$state,
                                        fill = .data$group,
                                        group = .data$group)) +
    ggplot2::geom_col(position = dodge, width = 0.68, alpha = 0.88,
                      colour = "white", linewidth = 0.25, na.rm = TRUE) +
    ggplot2::facet_wrap(~ measure, ncol = ncol, scales = scales) +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.16))
    ) +
    ggplot2::labs(x = "Centrality", y = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      position = dodge, hjust = -0.08, size = 2.5,
      show.legend = FALSE, na.rm = TRUE
    )
  }
  p
}

#' @noRd
.centrality_plot_data <- function(x) {
  d <- as.data.frame(x)
  measures <- setdiff(names(d), "state")
  out <- stats::reshape(
    d,
    varying = measures,
    v.names = "value",
    timevar = "measure",
    times = measures,
    idvar = "state",
    direction = "long"
  )
  out$state <- as.character(out$state)
  out$measure <- factor(out$measure, levels = measures)
  out$state_panel <- paste(out$state, out$measure, sep = "__")
  out$label <- .fmt_centrality(out$value)
  rownames(out) <- NULL
  out
}

#' @noRd
.plot_net_centrality_profile <- function(x, ncol = NULL,
                                         profile_scale = "measure",
                                         labels = FALSE) {
  .centrality_line_plot(.centrality_plot_data(x), group = FALSE, ncol = ncol,
                        profile_scale = profile_scale, labels = labels)
}

#' @noRd
.plot_net_centrality_group_profile <- function(d, ncol = NULL,
                                               palette = "Set2",
                                               profile_scale = "measure",
                                               labels = FALSE) {
  .centrality_line_plot(d, group = TRUE, ncol = ncol, palette = palette,
                        profile_scale = profile_scale, labels = labels)
}

#' Centrality line chart: states on x, one panel per measure, line per group
#' @noRd
.centrality_line_plot <- function(d, group = FALSE, ncol = NULL,
                                  palette = "Set2",
                                  profile_scale = "measure",
                                  labels = FALSE) {
  scaled <- identical(profile_scale, "measure")
  d <- .centrality_profile_data(d, profile_scale = "measure")  # value_plot, label
  # consistent state order across panels: highest mean centrality on the left
  ord <- stats::aggregate(value_plot ~ state, d, mean, na.rm = TRUE)
  ord <- ord[order(ord$value_plot, decreasing = TRUE), , drop = FALSE]
  d$state <- factor(as.character(d$state), levels = ord$state)
  if (!scaled) d$value_plot <- d$value      # display raw values, free y per panel
  if (is.null(ncol)) ncol <- min(3L, length(levels(d$measure)))
  y_lab <- if (scaled) "Centrality (within-measure scaled)" else "Centrality"
  facet_scales <- if (scaled) "fixed" else "free_y"

  if (group) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$state, y = .data$value_plot,
                                         colour = .data$group,
                                         group = .data$group)) +
      ggplot2::geom_line(linewidth = 0.85, alpha = 0.9, na.rm = TRUE) +
      ggplot2::geom_point(size = 2.1, na.rm = TRUE) +
      ggplot2::scale_colour_brewer(palette = palette)
  } else {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$state, y = .data$value_plot,
                                         group = 1L)) +
      ggplot2::geom_line(linewidth = 0.85, colour = "#2A6FBB", alpha = 0.9,
                         na.rm = TRUE) +
      ggplot2::geom_point(size = 2.1, colour = "#2A6FBB", na.rm = TRUE)
  }
  p <- p +
    ggplot2::facet_wrap(~ measure, ncol = ncol, scales = facet_scales) +
    ggplot2::labs(x = NULL, y = y_lab, colour = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = if (group) "bottom" else "none",
      panel.spacing = ggplot2::unit(1, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      vjust = -0.65, size = 2.3, show.legend = FALSE, na.rm = TRUE
    )
  }
  p
}

#' @noRd
.centrality_profile_data <- function(d, profile_scale = "measure") {
  d$value_plot <- d$value
  if (profile_scale == "measure") {
    split_idx <- split(seq_len(nrow(d)), d$measure)
    for (idx in split_idx) {
      vals <- d$value[idx]
      finite <- is.finite(vals)
      if (!any(finite)) {
        d$value_plot[idx] <- NA_real_
        next
      }
      mi <- min(vals[finite])
      ma <- max(vals[finite])
      if (!is.finite(mi) || !is.finite(ma) || isTRUE(all.equal(mi, ma))) {
        d$value_plot[idx] <- ifelse(finite, 0.5, NA_real_)
      } else {
        d$value_plot[idx] <- (vals - mi) / (ma - mi)
      }
    }
  }
  d$label <- .fmt_centrality(d$value)
  d
}

#' @noRd
.centrality_state_factor <- function(d, reorder = TRUE) {
  states <- unique(as.character(d$state))
  if (isTRUE(reorder)) {
    ord <- stats::aggregate(value ~ state, d, mean, na.rm = TRUE)
    ord <- ord[order(ord$value, na.last = TRUE), , drop = FALSE]
    states <- ord$state
  }
  factor(as.character(d$state), levels = states)
}

#' @noRd
.centrality_drop_zero <- function(d) {
  keep <- vapply(split(d$value, d$measure), function(v) {
    any(is.finite(v) & abs(v) > .Machine$double.eps^0.5)
  }, logical(1))
  kept <- names(keep)[keep]
  if (!length(kept)) return(d)
  d <- d[d$measure %in% kept, , drop = FALSE]
  d$measure <- factor(as.character(d$measure), levels = kept)
  rownames(d) <- NULL
  d
}

#' @noRd
.plot_net_centrality_heatmap <- function(d, reorder = TRUE, labels = TRUE) {
  d <- .centrality_profile_data(d, profile_scale = "measure")
  d$state <- .centrality_state_factor(d, reorder = reorder)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$measure, y = .data$state,
                                        fill = .data$value_plot)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6, na.rm = TRUE) +
    ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#2A6FBB",
                                 na.value = "grey92",
                                 limits = c(0, 1),
                                 name = "Scaled\n(per measure)") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label,
                   colour = .data$value_plot > 0.6),
      size = 2.8, show.legend = FALSE, na.rm = TRUE
    ) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white",
                                              `FALSE` = "grey20"))
  }
  p
}

#' @noRd
.plot_net_centrality_delta <- function(d, ncol = 3L, reorder = TRUE,
                                       labels = FALSE, palette = "Set2") {
  groups <- unique(d$group)
  if (length(groups) < 2L) {
    stop("type = \"delta\" needs at least two groups; got ", length(groups),
         ".", call. = FALSE)
  }
  if (length(groups) > 2L) {
    return(.plot_net_centrality_delta_multi(d, ncol = ncol, reorder = reorder,
                                            labels = labels, palette = palette))
  }
  wide <- stats::reshape(
    d[, c("state", "measure", "group", "value")],
    idvar = c("state", "measure"), timevar = "group",
    direction = "wide"
  )
  g1 <- paste0("value.", groups[1L])
  g2 <- paste0("value.", groups[2L])
  wide$delta <- wide[[g2]] - wide[[g1]]
  wide$measure <- factor(as.character(wide$measure), levels = levels(d$measure))
  wide$sign <- ifelse(wide$delta >= 0, groups[2L], groups[1L])
  wide$sign <- factor(wide$sign, levels = groups)
  ord <- stats::aggregate(delta ~ state, wide, mean, na.rm = TRUE)
  ord <- ord[order(ord$delta, na.last = TRUE), , drop = FALSE]
  wide$state <- factor(as.character(wide$state), levels = ord$state)
  wide$label <- .fmt_centrality(wide$delta, signed = TRUE)
  p <- ggplot2::ggplot(wide, ggplot2::aes(x = .data$delta, y = .data$state,
                                          fill = .data$sign)) +
    ggplot2::geom_col(width = 0.72, alpha = 0.9, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::facet_wrap(~ measure, ncol = ncol, scales = "free_x") +
    ggplot2::scale_fill_manual(
      values = stats::setNames(c("#D33F6A", "#4A6FE3"), groups),
      name = "Higher in"
    ) +
    ggplot2::labs(
      x = sprintf("Difference (%s − %s)", groups[2L], groups[1L]),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      hjust = ifelse(wide$delta >= 0, -0.1, 1.1),
      size = 2.5, na.rm = TRUE
    ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.14, 0.14))
      )
  }
  p
}

#' Diverging delta view for 3+ groups: deviation from the per-cell group mean
#' @noRd
.plot_net_centrality_delta_multi <- function(d, ncol = 3L, reorder = TRUE,
                                             labels = FALSE, palette = "Set2") {
  groups <- unique(d$group)
  # per (state, measure) mean across groups, then each group's deviation
  cell <- paste(d$state, d$measure, sep = "\r")
  mu <- stats::ave(d$value, cell, FUN = function(v) mean(v, na.rm = TRUE))
  d$dev <- d$value - mu
  d$state <- .centrality_state_factor(d, reorder = reorder)
  d$group <- factor(d$group, levels = groups)
  d$label <- .fmt_centrality(d$dev, signed = TRUE)
  dodge <- ggplot2::position_dodge2(width = 0.76, preserve = "single",
                                    padding = 0.08)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$dev, y = .data$state,
                                       fill = .data$group, group = .data$group)) +
    ggplot2::geom_col(position = dodge, width = 0.68, alpha = 0.9,
                      colour = "white", linewidth = 0.25, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    ggplot2::facet_wrap(~ measure, ncol = ncol, scales = "free_x") +
    ggplot2::scale_fill_brewer(palette = palette) +
    ggplot2::labs(x = "Deviation from group mean", y = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      panel.spacing = ggplot2::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      position = dodge, hjust = -0.08, size = 2.3,
      show.legend = FALSE, na.rm = TRUE
    ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.14, 0.14))
      )
  }
  p
}

#' @noRd
.centrality_measure_palette <- function() {
  c(
    OutStrength = "#2A6FBB",
    InStrength = "#D94F45",
    ClosenessIn = "#2F9E72",
    ClosenessOut = "#8E63B0",
    Closeness = "#D98C20",
    Betweenness = "#4F6D7A",
    BetweennessRSP = "#C44E7F",
    Diffusion = "#6F9E2F",
    Clustering = "#7B6D5F",
    InCloseness = "#2F9E72",
    OutCloseness = "#8E63B0"
  )
}


# ---- Exported net_edge_betweenness() ----

#' Edge Betweenness Network
#'
#' Builds a network in which each edge's weight is replaced by its
#' betweenness: the number of shortest paths between all node pairs that
#' traverse that edge (fractional when shortest paths tie). This is the
#' Nestimate counterpart of \code{tna::betweenness_network()} and produces
#' identical values for transition networks; the name differs to avoid a
#' clash with \code{tna::betweenness_network()} and
#' \code{igraph::edge_betweenness()}.
#'
#' For a probability/transition network the edge weights are transition
#' probabilities, so they are inverted to distances (\code{invert = TRUE})
#' before path-finding: the geodesic between two states is then the most
#' probable route rather than the one with the fewest hops. Pass
#' \code{invert = FALSE} when the weights already represent distances.
#'
#' Directedness is taken from the network itself. A directed network yields
#' an asymmetric betweenness matrix; an undirected (symmetric) network
#' yields a symmetric one.
#'
#' @param x A \code{netobject} or \code{netobject_group}.
#' @param invert Logical. Invert weights to distances by \code{1/w} before
#'   computing shortest paths? Default \code{TRUE} (correct for probability
#'   and frequency networks).
#' @param ... Additional arguments (ignored).
#'
#' @return For a \code{netobject}: a new \code{netobject} (class
#'   \code{c("netobject", "cograph_network")}) whose \code{$weights} are the
#'   edge-betweenness scores, with \code{method = "edge_betweenness"}. Call
#'   \code{extract_edges()} on it for a tidy per-edge table, or \code{plot()}
#'   to render it. For a \code{netobject_group}: a \code{netobject_group} of
#'   such networks, one per group.
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' net <- build_network(seqs, method = "relative")
#' eb  <- net_edge_betweenness(net)
#' extract_edges(eb)
#'
#' @export
net_edge_betweenness <- function(x, invert = TRUE, ...) {
  UseMethod("net_edge_betweenness")
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.netobject <- function(x, invert = TRUE, ...) {
  mat <- x$weights
  if (is.null(dimnames(mat)) && !is.null(x$nodes$name)) {
    dimnames(mat) <- list(x$nodes$name, x$nodes$name)
  }
  directed <- if (!is.null(x$directed)) isTRUE(x$directed) else TRUE
  eb <- .edge_betweenness(mat, invert = invert)
  out <- .wrap_netobject(eb, data = x$data, method = "edge_betweenness",
                         directed = directed, inits = x$inits)
  class(out) <- c("net_edge_betweenness", class(out))
  out
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.netobject_group <- function(x, invert = TRUE, ...) {
  result <- lapply(x, net_edge_betweenness.netobject, invert = invert)
  names(result) <- names(x)
  class(result) <- "netobject_group"
  result
}

#' @rdname net_edge_betweenness
#' @export
net_edge_betweenness.default <- function(x, invert = TRUE, ...) {
  stop("net_edge_betweenness() requires a 'netobject' or 'netobject_group'; ",
       "got '", class(x)[1L], "'.", call. = FALSE)
}

#' Plot edge-betweenness scores
#'
#' Draws the edges of a \code{\link{net_edge_betweenness}} network ranked by
#' their betweenness, as a horizontal bar chart. This is the tidy,
#' cograph-free companion to the node-link diagram: render the diagram with
#' \code{cograph::splot(eb)} and the ranking with \code{plot(eb)}.
#'
#' @param x A \code{net_edge_betweenness} network from
#'   \code{\link{net_edge_betweenness}}.
#' @param style Plot style. \code{"bar"} (default) draws one horizontal bar
#'   per edge; \code{"forest"} draws a forest/lollipop chart (a stem from
#'   zero to a point) with a dashed reference line at the mean betweenness;
#'   \code{"delta"} draws each edge's deviation from the mean edge betweenness
#'   as a diverging bar (above the mean in blue, below in red).
#' @param top_n Integer or \code{NULL}. Keep only the \code{top_n} highest
#'   edges. Default \code{NULL} (all edges with non-zero betweenness).
#' @param labels Logical. Print the betweenness value beside each edge.
#'   Default \code{TRUE}.
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object.
#' @examples
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"))
#' eb <- net_edge_betweenness(build_network(seqs, method = "relative"))
#' plot(eb)
#' plot(eb, style = "forest")
#' plot(eb, style = "delta")
#' @export
plot.net_edge_betweenness <- function(x, style = c("bar", "forest", "delta"),
                                      top_n = NULL, labels = TRUE, ...) {
  style <- match.arg(style)
  e <- extract_edges(x)
  e <- e[is.finite(e$weight) & e$weight > 0, , drop = FALSE]
  if (!nrow(e)) {
    stop("No edges with non-zero betweenness to plot.", call. = FALSE)
  }
  e <- e[order(e$weight, decreasing = TRUE), , drop = FALSE]
  if (!is.null(top_n)) e <- utils::head(e, top_n)
  e$edge <- factor(sprintf("%s -> %s", e$from, e$to),
                   levels = rev(sprintf("%s -> %s", e$from, e$to)))
  e$label <- .fmt_centrality(e$weight)

  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 14, 8, 8)
    )

  if (style == "delta") {
    ref <- mean(e$weight, na.rm = TRUE)
    e$delta <- e$weight - ref
    # re-order edges by the deviation (most-above at top)
    e$edge <- factor(sprintf("%s -> %s", e$from, e$to),
                     levels = sprintf("%s -> %s", e$from, e$to)[order(e$delta)])
    e$sign <- factor(ifelse(e$delta >= 0, "above", "below"),
                     levels = c("below", "above"))
    e$label <- .fmt_centrality(e$delta, signed = TRUE)
    p <- ggplot2::ggplot(e, ggplot2::aes(x = .data$delta, y = .data$edge,
                                         fill = .data$sign)) +
      ggplot2::geom_col(width = 0.72, alpha = 0.9, na.rm = TRUE) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
      ggplot2::scale_fill_manual(
        values = c(below = "#D33F6A", above = "#4A6FE3"),
        labels = c(below = "below mean", above = "above mean"), name = NULL
      ) +
      ggplot2::labs(x = sprintf("Edge betweenness − mean (%.2f)", ref),
                    y = NULL) +
      base_theme + ggplot2::theme(legend.position = "bottom")
    if (isTRUE(labels)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$label),
        hjust = ifelse(e$delta >= 0, -0.1, 1.1), size = 2.8, na.rm = TRUE
      ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0.14, 0.14))
        )
    }
    return(p)
  }

  if (style == "forest") {
    ref <- mean(e$weight, na.rm = TRUE)
    p <- ggplot2::ggplot(e, ggplot2::aes(x = .data$weight, y = .data$edge)) +
      ggplot2::geom_vline(xintercept = ref, linetype = "dashed",
                          colour = "grey55", linewidth = 0.4) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = .data$weight,
                     y = .data$edge, yend = .data$edge),
        colour = "#9DB4C0", linewidth = 1.1, na.rm = TRUE
      ) +
      ggplot2::geom_point(size = 3.1, colour = "#4F6D7A", na.rm = TRUE) +
      ggplot2::annotate("text", x = ref, y = Inf,
                        label = sprintf("mean = %.2f", ref),
                        vjust = 1.4, hjust = -0.05, size = 2.8,
                        colour = "grey45") +
      ggplot2::labs(x = "Edge betweenness", y = NULL) +
      base_theme
    if (isTRUE(labels)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = .data$label), hjust = -0.5, size = 2.8,
        na.rm = TRUE
      ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0.02, 0.18))
        )
    }
    return(p)
  }

  p <- ggplot2::ggplot(e, ggplot2::aes(x = .data$weight, y = .data$edge)) +
    ggplot2::geom_col(width = 0.72, alpha = 0.9, fill = "#4F6D7A",
                      na.rm = TRUE) +
    ggplot2::labs(x = "Edge betweenness", y = NULL) +
    base_theme
  if (isTRUE(labels)) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$label), hjust = -0.12, size = 2.8,
      na.rm = TRUE
    ) +
      ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(mult = c(0.02, 0.16))
      )
  }
  p
}
