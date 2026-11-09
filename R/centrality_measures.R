# ---- Built-in centrality measures (no external dependencies) ----
#
# All path-based measures (Betweenness, Closeness) are derived from
# all-pairs shortest paths via Floyd-Warshall. For transition/frequency
# networks, weights are inverted so higher weight = shorter distance.

# Tolerance for "same path length": relative, so a tie that floating-point
# addition breaks (1/10 + 1/15 vs 1/6) is still a tie. Shared by the path
# counting and the on-path tests so both agree on what a tie is.
.sp_tol <- function(d) 1e-10 * pmax(1, abs(d))


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
    tol            <- .sp_tol(ifelse(is.finite(D), D, 0))
    shorter        <- is.finite(new_d) & (!is.finite(D) | new_d < D - tol)
    equal          <- is.finite(new_d) & is.finite(D) &
                      abs(new_d - D) <= tol & new_d > 0
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
               abs(d_svt - d_st) <= .sp_tol(d_st)
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
    on_path <- is.finite(D) & sg > 0L & abs(through - D) <= .sp_tol(D)
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
  # reach[s, t]: t can be reached from s (or s == t). Only these pairs carry
  # randomized-shortest-path walks; the unreachable entries of Z are zero in
  # exact arithmetic and their solver noise must not enter.
  reach <- is.finite(.floyd_warshall_sp(1 * (mat > 0), invert = FALSE)$D)
  Z_recip <- Z^-1
  Z_recip[!reach | !is.finite(Z_recip)] <- 0
  # The correction for walks absorbed at target t is counted once per source
  # that reaches t: n of them in a strongly connected network, fewer when a
  # state is never entered or never left. Until 1.0.5 it was always n.
  n_sources <- .colSums(reach, n, n)
  Z_recip_diag <- diag(n_sources * diag(Z_recip), n, n)
  out <- diag(tcrossprod(Z, Z_recip - Z_recip_diag) %*% Z)
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
#'
#' A numerically constant vector has no min-max normalization: it returns
#' NaN, as tna's ranger() does, instead of stretching rounding error to 0-1.
#' @noRd
.range01_tna <- function(x, na.rm = FALSE) {
  mi <- min(x, na.rm = na.rm)
  ma <- max(x, na.rm = na.rm)
  span <- ma - mi
  if (!isTRUE(span > 1e-10 * max(abs(c(mi, ma))))) {
    return(rep(NaN, length(x)))
  }
  (x - mi) / span
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
#' @description
#' `net_centrality()` scores how central each node is in a network, for one
#' or more measures, and returns one row per node. Strength sums a node's
#' edge weights (incoming or outgoing); closeness and betweenness read the
#' shortest paths between nodes, where a heavier edge is a shorter step;
#' diffusion counts the weighted walks that start at a node; and clustering
#' measures how strongly a node's neighbours are tied to each other. It
#' takes a single network or a group of them (a `netobject_group`, the
#' clusters of an `mcml`), scoring each separately.
#'
#' Centralities describe the estimated network, not the population it was
#' estimated from. Rankings built on a few weak edges change from sample to
#' sample: check them with [centrality_stability()] before reporting which
#' node is most central.
#'
#' @details
#' Strength is read from the weight matrix directly, and closeness and
#' betweenness from all-pairs shortest paths (Floyd-Warshall), without
#' igraph; the values agree with `tna::centralities()`. Diffusion is
#' range-normalised by default, which `tna` does not do.
#'
#' Range normalisation of a measure that is the same for every node is
#' undefined and gives `NaN`, with a warning of class
#' `nestimate_constant_measure`. Diffusion on a `"relative"` network with
#' `loops = TRUE` is such a case: each row sums to 1, so every node has the
#' same diffusion; use `normalize_diffusion = FALSE` there.
#'
#' @param x The network: a `netobject`, `netobject_group`, `mcml` or
#'   `cograph_network`. For the `plot()` methods: the result of
#'   `net_centrality()`.
#' @param measures The measures to compute. The default is `"InStrength"`,
#'   `"Betweenness"` and `"Diffusion"`; `"all"` computes every built-in one:
#'   `"OutStrength"` and `"InStrength"` (sum of outgoing or incoming weights),
#'   `"ClosenessOut"`, `"ClosenessIn"` and `"Closeness"` (inverse of the
#'   summed shortest-path distances from, to, or ignoring the direction of
#'   edges), `"Betweenness"` (how many shortest paths pass through the node),
#'   `"BetweennessRSP"` (betweenness over randomised shortest paths, which
#'   also counts near-shortest routes), `"Diffusion"` (the summed weights of
#'   walks of length 1 to n from the node) and `"Clustering"` (weighted
#'   clustering coefficient). `"InCloseness"` and `"OutCloseness"` are
#'   accepted for `"ClosenessIn"` and `"ClosenessOut"`. With `centrality_fn`,
#'   name its measures here to have them returned.
#' @param loops Whether self-loops (staying in a state) count. With `FALSE`
#'   (default) they are left out of every measure, and a message says so.
#' @param normalize Whether to range-normalise every measure to 0--1, as
#'   `tna::centralities(normalize = TRUE)` does. Default `FALSE`.
#' @param invert Whether a heavier edge is a shorter step in the
#'   shortest-path measures (distance = 1 / weight), as in `tna`. With
#'   `FALSE`, the weight itself is the distance. Default `TRUE`.
#' @param normalize_diffusion Whether to range-normalise `Diffusion` even
#'   when `normalize = FALSE`. Default `TRUE`.
#' @param centrality_fn A function of your own that takes the weight matrix
#'   and returns a named list of centrality vectors, one per measure; list
#'   their names in `measures`. Default `NULL`.
#' @param ... Not used.
#' @return For a `netobject` or `cograph_network`: a `net_centrality` data
#'   frame with one row per node, a `state` column and one column per
#'   requested measure (the node names are also the row names). For a
#'   `netobject_group` or an `mcml`: a `net_centrality_group`, a list of such
#'   data frames, one per network. `plot()` draws the measures per node.
#'
#' @references
#' Freeman, L. C. (1978). Centrality in social networks: conceptual
#' clarification. \emph{Social Networks}, 1(3), 215--239. (betweenness,
#' closeness)
#'
#' Opsahl, T., Agneessens, F. & Skvoretz, J. (2010). Node centrality in
#' weighted networks: generalizing degree and shortest paths.
#' \emph{Social Networks}, 32(3), 245--251. (weighted strength and geodesics)
#'
#' Kivimaki, I., Lebichot, B., Saramaki, J. & Saerens, M. (2016). Two
#' betweenness centrality measures based on randomized shortest paths.
#' \emph{Scientific Reports}, 6, 19668. (\code{BetweennessRSP})
#'
#' Banerjee, A., Chandrasekhar, A. G., Duflo, E. & Jackson, M. O. (2013). The
#' diffusion of microfinance. \emph{Science}, 341(6144), 1236498.
#' (\code{Diffusion})
#'
#' Onnela, J.-P., Saramaki, J., Kertesz, J. & Kaski, K. (2005). Intensity and
#' coherence of motifs in weighted complex networks. \emph{Physical Review E},
#' 71, 065103. (\code{Clustering})
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
  # a measure equal for every node has no range to normalise: tna's NaN is
  # kept, and said
  normalised <- if (isTRUE(normalize)) names(res) else
    if (isTRUE(normalize_diffusion)) intersect(names(res), "Diffusion")
  undefined <- Filter(function(m) all(is.nan(res[[m]])), normalised)
  if (length(undefined)) {
    warning(warningCondition(sprintf(paste0(
      "%s is the same for every node, so it cannot be range-normalised and ",
      "is NaN. %s"), paste(undefined, collapse = ", "),
      if (identical(undefined, "Diffusion") && !isTRUE(normalize)) {
        "Use `normalize_diffusion = FALSE` for the raw values."
      } else "Use `normalize = FALSE` for the raw values."),
      class = "nestimate_constant_measure", call = NULL))
  }
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

#' @rdname net_centrality
#' @param reorder For the `plot()` methods: whether to sort the nodes in each
#'   measure's panel by their value (for a group, by their mean value across
#'   networks). Default `TRUE`.
#' @param type For the `plot()` methods: the kind of plot. `"bar"` (default)
#'   draws one panel of bars per measure, grouped bars for several networks.
#'   `"line"` draws each node's profile across the measures (for several
#'   networks, one panel per node and one line per network); `"profile"` is
#'   accepted for it. A single network also takes `"heatmap"`, a
#'   nodes-by-measures grid with each measure scaled to 0--1 for the colour
#'   and its raw value printed in the tile. Several networks also take
#'   `"delta"`, bars of the differences between networks: the second minus
#'   the first with two, each network's deviation from the per-node mean with
#'   more.
#' @param ncol For the `plot()` methods: the number of panel columns.
#'   Default `3`.
#' @param scales For the `plot()` methods: `"free_x"` (default) gives each
#'   measure its own value axis, `"fixed"` one shared axis.
#' @param profile_scale For `type = "line"`: `"measure"` (default) rescales
#'   each measure to 0--1 so the profiles compare across measures; `"none"`
#'   draws raw values.
#' @param labels For the `plot()` methods: whether to print the values on the
#'   plot. Default `TRUE` for one network, `FALSE` for several.
#' @param drop_zero For the `plot()` methods: whether to leave out measures
#'   that are zero for every node. Default `FALSE`.
#' @param palette For `plot()` of several networks: the ColorBrewer palette
#'   for the networks. Default `"Set2"`. Ignored by `type = "delta"` with two
#'   networks, whose bars are coloured by which network is higher.
#' @return For the `plot()` methods: a `ggplot` object.
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

#' @rdname net_centrality
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
    stop(errorCondition(
      paste0(
        "type = \"delta\" needs at least two groups; got ",
        length(groups),
        "."),
      class = "nestimate_bad_input", call = NULL))
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
      x = sprintf("Difference (%s - %s)", groups[2L], groups[1L]),
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
#' @description
#' `net_edge_betweenness()` returns a copy of a network in which each edge's
#' weight is replaced by its betweenness: the number of shortest paths
#' between ordered pairs of distinct nodes that pass through that edge, with
#' a fractional share when several shortest paths tie (Girvan & Newman
#' 2002). It
#' takes a single network or a `netobject_group`, which is scored one
#' network at a time.
#'
#' With `invert = TRUE` (the default) a heavier edge is a shorter step, so
#' the shortest path between two states of a transition network is the most
#' probable route, not the one with the fewest hops.
#' Only positive weights count as edges, and a self-loop never lies on a
#' shortest path, so its betweenness is 0. A directed network gives an
#' asymmetric betweenness matrix; an undirected one a symmetric matrix,
#' each cell holding the betweenness of that edge. The result keeps the
#' source network's directedness (directed when it records none).
#'
#' An edge with high betweenness carries many of the shortest routes
#' between states; Girvan & Newman (2002) read such edges as the bridges
#' between groups of states. To compare the edge betweenness of
#' two groups, pass two results to [permutation()] or [bayes_compare()],
#' which resample the source networks behind them; both need two results
#' with the same source method and `invert`, and [bayes_compare()] needs
#' `"relative"` source networks.
#'
#' @details
#' Shortest paths come from all-pairs Floyd-Warshall on the distance
#' matrix (1 / weight, or the weight with `invert = FALSE`), counting tied
#' paths with a relative tolerance so that floating-point rounding does not
#' break a tie. The values agree with `tna::betweenness_network()` for
#' transition networks, and on an undirected network each cell equals
#' `igraph::edge_betweenness()` for that edge; the function is named
#' differently to avoid masking either.
#'
#' @param x A `netobject` or `netobject_group`; any other input raises an
#'   error. For the `plot()` method: a `net_edge_betweenness` network.
#' @param invert Whether to turn weights into distances by `1 / weight`
#'   before finding shortest paths, which suits probability and frequency
#'   networks. With `FALSE`, the weight itself is the distance; use it when
#'   the weights already are distances. Default `TRUE`.
#' @param ... Not used.
#' @return For a `netobject`: a network of class
#'   `c("net_edge_betweenness", "netobject", "cograph_network")` with method
#'   `"edge_betweenness"`, whose weights are the edge betweenness scores.
#'   It prints like any network; [extract_edges()] gives one row per edge,
#'   `plot()` ranks the edges and `cograph::splot()` draws the network. The
#'   source network's data, metadata, parameters, weights and initial
#'   probabilities are kept on it, which is what lets [permutation()] and
#'   [bayes_compare()] test differences in edge betweenness between two
#'   such networks. For a `netobject_group`: a `netobject_group` of these
#'   networks, one per group.
#'
#' @references
#' Girvan, M. & Newman, M. E. J. (2002). Community structure in social and
#' biological networks. \emph{Proceedings of the National Academy of
#' Sciences}, 99(12), 7821--7826. \doi{10.1073/pnas.122653799}
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
  out$metadata <- x$metadata
  out$params <- x$params %||% list()
  out$scaling <- x$scaling
  out$threshold <- x$threshold %||% 0
  out$level <- x$level
  out$build_args <- x$build_args
  out$source_weights <- mat
  out$frequency_matrix <- x$frequency_matrix
  out$edge_betweenness <- list(
    source_method = x$method,
    invert = isTRUE(invert)
  )
  out$meta$edge_betweenness <- out$edge_betweenness
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
  stop(errorCondition(
    paste0(
      "net_edge_betweenness() requires a 'netobject' or 'netobject_group'; ",
      "got '",
      class(x)[1L],
      "'."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname net_edge_betweenness
#' @section Methods:
#' `plot()` ranks the edges with non-zero betweenness, highest first, one
#' row per edge (`from -> to`). It needs no cograph; `cograph::splot()`
#' draws the same object as a node-link diagram. A network whose edges all
#' have zero betweenness raises an error.
#' @param style For the `plot()` method: the kind of chart. `"bar"`
#'   (default) draws one horizontal bar per edge; `"forest"` a stem from
#'   zero to a point, with a dashed line at the mean betweenness of the
#'   edges shown; `"delta"` each edge's deviation from that mean as a bar
#'   pointing right (above the mean, blue) or left (below, red).
#' @param top_n For the `plot()` method: the number of highest edges to
#'   keep. Default `NULL` (every edge with non-zero betweenness).
#' @param labels For the `plot()` method: whether to print each edge's
#'   value (its deviation, for `"delta"`) beside it. Default `TRUE`.
#' @return For the `plot()` method: a `ggplot` object.
#' @export
plot.net_edge_betweenness <- function(x, style = c("bar", "forest", "delta"),
                                      top_n = NULL, labels = TRUE, ...) {
  style <- match.arg(style)
  e <- extract_edges(x)
  e <- e[is.finite(e$weight) & e$weight > 0, , drop = FALSE]
  if (!nrow(e)) {
    stop(errorCondition("No edges with non-zero betweenness to plot.",
                        class = "nestimate_bad_input", call = NULL))
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
      ggplot2::labs(x = sprintf("Edge betweenness - mean (%.2f)", ref),
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
