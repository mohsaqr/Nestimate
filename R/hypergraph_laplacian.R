# ---- Hypergraph Laplacians, Spectral Clustering, Transduction ----------
#
# Zhou, Huang & Scholkopf (2006): normalized hypergraph Laplacian
#   L = I - D_v^{-1/2} H W D_e^{-1} H^T D_v^{-1/2}
# on the BINARY incidence pattern H with hyperedge weights w(e).
#
# Hayashi, Aksoy, Park & Park (2020): random-walk Laplacian with
# edge-dependent vertex weights (EDVW). The WEIGHTED incidence cells
# gamma_e(v) drive the walk (pick e ~ w(e), then v ~ gamma_e(v)); the
# non-reversible walk is symmetrized through its stationary distribution
# with Chung's (2005) directed Laplacian. With a binary incidence and unit
# hyperedge weights the two Laplacians coincide (the walk is reversible
# with pi ~ d(v)) - pinned by an invariant test.
#
# Reference implementation used as the equivalence oracle: HyperNetX
# `hypernetx/algorithms/clustering/laplacians_clustering.py` (PNNL;
# author-adjacent). Its hyperedge-weight default - the population standard
# deviation of the edge's vertex weights plus one - is reproduced here as
# `edge_weights = NULL` for `type = "random_walk"`.

#' Normalised hypergraph Laplacian
#'
#' @description
#' `hypergraph_laplacian()` computes the normalised Laplacian of a hypergraph
#' and returns it as a symmetric node-by-node matrix, the operator that
#' [hypergraph_cluster()] and [hypergraph_transduction()] are built on. It
#' comes in two forms. `type = "zhou"` is the Laplacian of Zhou, Huang and
#' Sch\enc{ö}{o}lkopf (2006), computed on the binary incidence pattern, so
#' the values in the incidence cells play no part. `type = "random_walk"` is
#' the Laplacian of Hayashi et al. (2020) for a random walk with
#' edge-dependent vertex weights: the weighted incidence cells (for example
#' the summed weights [bipartite_groups()] produces with `weight`) set where
#' a walker that enters a hyperedge lands, and the walk, which is not
#' reversible, is symmetrised through its stationary distribution (Chung
#' 2005).
#'
#' Both Laplacians are symmetric positive semi-definite with eigenvalues in
#' \[0, 2\]. On a binary incidence with unit hyperedge weights the two forms
#' coincide. The hypergraph must be connected, so the eigenvalue 0 occurs
#' once; further eigenvalues close to 0 indicate groups of nodes that the
#' walk seldom leaves, which is what [hypergraph_cluster()] partitions by.
#'
#' @details
#' For `"zhou"`, with \eqn{H} the binary incidence, \eqn{W} the diagonal of
#' hyperedge weights, \eqn{D_e} the hyperedge sizes and \eqn{D_v} the node
#' degrees (the summed weights of the hyperedges a node belongs to),
#' \deqn{L = I - D_v^{-1/2} H W D_e^{-1} H^T D_v^{-1/2},}
#' and the stationary distribution is proportional to \eqn{D_v}. For
#' `"random_walk"`, a walker at node \eqn{v} picks a hyperedge containing
#' it with probability proportional to its weight, then a node \eqn{u} of
#' that hyperedge with probability proportional to the incidence cell
#' \eqn{\gamma_e(u)}. The stationary distribution \eqn{\pi} is the leading
#' left eigenvector of this transition matrix \eqn{P}, and
#' \eqn{L = I - (\Pi^{1/2} P \Pi^{-1/2} + \Pi^{-1/2} P^T \Pi^{1/2}) / 2}.
#' The random-walk Laplacian, its stationary distribution and its default
#' hyperedge weights follow the reference implementation in HyperNetX
#' (`laplacians_clustering.py`).
#'
#' @param hg The hypergraph: a `net_hypergraph` from [build_hypergraph()] or
#'   [bipartite_groups()], with at least two nodes and one hyperedge. A
#'   hypergraph that is not connected has no unique stationary distribution
#'   and raises an error of class `nestimate_hypergraph_disconnected`;
#'   analyse its components separately.
#' @param type The Laplacian. `"zhou"` (default) is the Zhou et al. (2006)
#'   Laplacian on the binary incidence pattern; `"random_walk"` is the
#'   Hayashi et al. (2020) random-walk Laplacian on the weighted incidence.
#' @param edge_weights Positive hyperedge weights, one per hyperedge in the
#'   order of the incidence columns. Default `NULL`, which gives unit weights
#'   for `"zhou"` and, for `"random_walk"`, the heuristic of Hayashi et al.:
#'   the population standard deviation of each hyperedge's non-zero
#'   incidence cells plus one, which is 1 on a binary incidence. A vector of
#'   the wrong length, or with a value that is not finite and positive,
#'   raises an error.
#'
#' @return A symmetric numeric matrix with one row and one column per node
#'   (the node names as dimnames). Its attributes record the `type`, `pi`
#'   (the stationary distribution of the walk, named by node) and
#'   `edge_weights` (the hyperedge weights used).
#'
#' @references
#' Zhou, D., Huang, J., & Sch\enc{ö}{o}lkopf, B. (2006). Learning with
#' hypergraphs: Clustering, classification, and embedding.
#' \emph{NeurIPS 19}.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{CIKM 2020}, 495-504.
#' \doi{10.1145/3340531.3412034}
#'
#' Chung, F. (2005). Laplacians and the Cheeger inequality for directed
#' graphs. \emph{Annals of Combinatorics}, 9(1), 1-19.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4"),
#'   hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1)
#' )
#' hg <- bipartite_groups(events, player = "person", group = "meeting",
#'                        weight = "hours")
#' L <- hypergraph_laplacian(hg, type = "random_walk")
#' range(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
#'
#' @export
hypergraph_laplacian <- function(hg,
                                 type = c("zhou", "random_walk"),
                                 edge_weights = NULL) {
  type <- match.arg(type)
  .hl_validate_hg(hg)
  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  L <- parts$L
  attr(L, "type") <- type
  attr(L, "pi") <- parts$pi
  attr(L, "edge_weights") <- parts$w
  L
}

#' Spectral clustering of hypergraph vertices
#'
#' @description
#' `hypergraph_cluster()` partitions the nodes of a hypergraph into `k`
#' clusters with the spectral algorithm of Hayashi et al. (2020,
#' "RDC-Spec"): it takes the eigenvectors of the `k` smallest eigenvalues of
#' the normalised Laplacian from [hypergraph_laplacian()], scales each
#' node's row of them to unit length, and clusters the rows with k-means.
#' It works with either Laplacian. With `type = "random_walk"` and a
#' weighted incidence (for example from [bipartite_groups()] with
#' `weight`), the edge-dependent vertex weights can change the partition;
#' with edge-independent weights the walk reduces to a random walk on a
#' graph (Chitra & Raphael 2019).
#'
#' The eigengap, the distance between the `k`-th and the next eigenvalue,
#' shows how well the spectrum separates `k` groups: a large gap supports
#' `k`, a small one means the spectrum does not single it out. `plot()`
#' draws the spectrum, so neighbouring values of `k` can be compared.
#' k-means is stochastic, so a partition that matters is worth refitting
#' with other seeds to see whether it holds.
#'
#' @param hg The hypergraph: a connected `net_hypergraph`, as for
#'   [hypergraph_laplacian()].
#' @param k The number of clusters, a whole number from 2 to one less than
#'   the number of nodes. Required; a value that is not a single whole
#'   number in that range raises an error.
#' @param type,edge_weights The Laplacian and its hyperedge weights, passed
#'   to [hypergraph_laplacian()]. Defaults `"zhou"` and `NULL`.
#' @param nstart The number of random starts k-means tries, keeping the
#'   best. Default `25`.
#' @param seed A whole number that makes the k-means starts reproducible.
#'   The session's own random stream is left as it was. Default `NULL`.
#'
#' @return A `net_hypergraph_cluster` object. `print()` reports the
#'   Laplacian, the cluster sizes and the eigengap; [summary()] returns one
#'   row per cluster; `plot()` draws the spectrum and the embedding. The
#'   object's fields:
#'   \describe{
#'     \item{clusters}{A data frame with one row per node: `node` and
#'       `cluster`. Clusters are labelled `"Cluster 1"`, `"Cluster 2"`, ...
#'       in the order they first appear among the nodes.}
#'     \item{embedding}{The node-by-`k` spectral embedding k-means
#'       clustered, each row of unit length, with columns `dim1` to
#'       `dimk`.}
#'     \item{k, type}{The number of clusters and the Laplacian type.}
#'     \item{eigenvalues}{The full spectrum of the Laplacian, increasing.}
#'     \item{eigengap}{The `(k + 1)`-th minus the `k`-th eigenvalue.}
#'     \item{sizes}{A data frame of `cluster` and `size`.}
#'     \item{pi}{The stationary distribution of the walk, named by node.}
#'     \item{n_nodes, n_hyperedges}{The size of the hypergraph.}
#'     \item{params}{The `edge_weights` used, `nstart`, `seed` and the
#'       k-means total within-cluster sum of squares (`tot_withinss`).}
#'   }
#'
#' @references
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{CIKM 2020}, 495-504.
#' \doi{10.1145/3340531.3412034}
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{ICML 2019}.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
#'              "d", "e", "f", "c", "d"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4", "m4", "m5", "m5")
#' )
#' hg <- bipartite_groups(events, player = "person", group = "meeting")
#' cl <- hypergraph_cluster(hg, k = 2, seed = 1)
#' cl
#' as.data.frame(cl)
#'
#' @export
hypergraph_cluster <- function(hg, k,
                               type = c("zhou", "random_walk"),
                               edge_weights = NULL,
                               nstart = 25L,
                               seed = NULL) {
  type <- match.arg(type)
  .hl_validate_hg(hg)
  .assert(
    "`k` must be a single whole number" =
      is.numeric(k) && length(k) == 1L && is.finite(k) && k == round(k),
    "`nstart` must be a single positive whole number" =
      is.numeric(nstart) && length(nstart) == 1L && nstart >= 1
  )
  k <- as.integer(k)
  n <- hg$n_nodes
  if (k < 2L || k > n - 1L) {
    stop(errorCondition(
      sprintf("`k` must be between 2 and n_nodes - 1 (= %d), got %d.",
        n - 1L, k),
      class = "nestimate_bad_input", call = NULL))
  }
  if (!is.null(seed)) .local_seed(as.integer(seed))

  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  eig <- eigen((parts$L + t(parts$L)) / 2, symmetric = TRUE)
  # eigen() returns decreasing order; take the k SMALLEST, increasing
  ord <- rev(seq_len(n))
  values <- eig$values[ord]
  U <- eig$vectors[, ord[seq_len(k)], drop = FALSE]

  # Row-normalize the spectral embedding to unit length (zero rows kept)
  row_norm <- sqrt(rowSums(U^2))
  nz <- row_norm > 0
  U[nz, ] <- U[nz, , drop = FALSE] / row_norm[nz]

  km <- stats::kmeans(U, centers = k, nstart = as.integer(nstart),
                      iter.max = 100L)

  # Deterministic labels: "Cluster 1" = first node's cluster, etc.
  relabel <- match(km$cluster, unique(km$cluster))
  cluster_lab <- paste("Cluster", relabel)
  clusters <- data.frame(node = hg$nodes, cluster = cluster_lab,
                         stringsAsFactors = FALSE)
  dimnames(U) <- list(hg$nodes, paste0("dim", seq_len(k)))
  sizes <- as.data.frame(table(cluster = cluster_lab),
                         stringsAsFactors = FALSE)
  names(sizes) <- c("cluster", "size")
  sizes <- sizes[order(sizes$cluster), , drop = FALSE]
  rownames(sizes) <- NULL

  structure(
    list(
      clusters    = clusters,
      embedding   = U,
      k           = k,
      type        = type,
      eigenvalues = values,
      eigengap    = if (k < n) values[k + 1L] - values[k] else NA_real_,
      sizes       = sizes,
      pi          = parts$pi,
      n_nodes     = n,
      n_hyperedges = hg$n_hyperedges,
      params = list(edge_weights = parts$w, nstart = as.integer(nstart),
                    seed = seed, tot_withinss = km$tot.withinss)
    ),
    class = "net_hypergraph_cluster"
  )
}

#' Transductive label spreading on a hypergraph
#'
#' @description
#' `hypergraph_transduction()` classifies the unlabelled nodes of a
#' hypergraph from a few labelled ones, by the label-spreading framework of
#' Zhou et al. (2006). The labels spread over the hypergraph as the scores
#' \deqn{F = (1 - \xi) (I - \xi S)^{-1} Y,}
#' where \eqn{S = I - L} is the similarity operator of the chosen Laplacian
#' (see [hypergraph_laplacian()]) and \eqn{Y} marks each labelled node's
#' class, and each node is assigned the class with the highest score. This
#' is the non-neural ancestor of hypergraph-attention text classifiers: with
#' documents as hyperedges over words (or the reverse), it classifies
#' unlabelled nodes from a handful of labelled ones.
#'
#' The `margin` of a node, its winning minus its runner-up score, says how
#' clearly it was assigned; assignments with a small margin are the ones to
#' distrust. The given labels are not held fixed, so a labelled node can be
#' predicted into another class, and a class's scores grow with its number
#' of labelled nodes. As `xi` approaches 1, every node's scores become
#' proportional to the same vector, so in the limit every node is assigned
#' the same class. On a small hypergraph the default `0.99` can come close
#' to that; compare the predictions under a smaller `xi`.
#'
#' @param hg The hypergraph: a connected `net_hypergraph`, as for
#'   [hypergraph_laplacian()].
#' @param labels The known class labels, in one of two forms: a vector
#'   named by node, covering some of the nodes, or an unnamed vector with one
#'   entry per node, in the hypergraph's node order, and `NA` for the
#'   unlabelled ones. Required. It must hold at least two distinct classes;
#'   a name that is not a node, or an unnamed vector of the wrong length,
#'   raises an error.
#' @param xi The spreading coefficient, a number strictly between 0 and 1.
#'   Larger values give the hypergraph structure more weight relative to the
#'   given labels. Default `0.99`.
#' @param type,edge_weights The Laplacian and its hyperedge weights, passed
#'   to [hypergraph_laplacian()]. Defaults `"zhou"` and `NULL`.
#'
#' @return A `net_hypergraph_transduction` object. `print()` reports the
#'   Laplacian, `xi`, the classes and how many nodes each was predicted
#'   for; [summary()] returns one row per class; `plot()` draws the
#'   node-by-class score heatmap. The object's fields:
#'   \describe{
#'     \item{predictions}{A data frame with one row per node: `node`,
#'       `label` (the given label, `NA` for an unlabelled node),
#'       `predicted`, `score` (the winning class's score) and `margin` (the
#'       winning minus the runner-up score).}
#'     \item{classes}{The classes, sorted.}
#'     \item{scores}{The node-by-class score matrix.}
#'     \item{xi, type}{The spreading coefficient and the Laplacian type.}
#'     \item{n_labeled, n_nodes}{The numbers of labelled and of all nodes.}
#'     \item{params}{The `edge_weights` used.}
#'   }
#'
#' @references
#' Zhou, D., Huang, J., & Sch\enc{ö}{o}lkopf, B. (2006). Learning with
#' hypergraphs: Clustering, classification, and embedding.
#' \emph{NeurIPS 19}.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
#'              "d", "e", "f", "c", "d"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4", "m4", "m5", "m5")
#' )
#' hg <- bipartite_groups(events, player = "person", group = "meeting")
#' tr <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"))
#' tr
#' as.data.frame(tr)
#'
#' @export
hypergraph_transduction <- function(hg, labels, xi = 0.99,
                                    type = c("zhou", "random_walk"),
                                    edge_weights = NULL) {
  type <- match.arg(type)
  .hl_validate_hg(hg)
  .assert(
    "`xi` must be a single number in (0, 1)" =
      is.numeric(xi) && length(xi) == 1L && xi > 0 && xi < 1
  )
  n <- hg$n_nodes
  nodes <- hg$nodes

  # Normalize `labels` to a full-length character vector with NAs
  lab <- if (!is.null(names(labels))) {
    unknown <- setdiff(names(labels), nodes)
    if (length(unknown) > 0L) {
      stop(errorCondition(paste0("Unknown node names in `labels`: ",
           paste(unknown, collapse = ", ")),
                          class = "nestimate_bad_input", call = NULL))
    }
    full <- rep(NA_character_, n)
    full[match(names(labels), nodes)] <- as.character(labels)
    full
  } else {
    if (length(labels) != n) {
      stop(errorCondition(
        sprintf(paste0("Unnamed `labels` must have length n_nodes ",
          "(= %d), got %d."), n, length(labels)),
        class = "nestimate_bad_input", call = NULL))
    }
    as.character(labels)
  }
  classes <- sort(unique(lab[!is.na(lab)]))
  if (length(classes) < 2L) {
    stop(errorCondition("`labels` must contain at least two distinct classes.",
         class = "nestimate_bad_input", call = NULL))
  }

  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  S <- diag(n) - (parts$L + t(parts$L)) / 2
  Y <- vapply(classes, function(cl) as.numeric(!is.na(lab) & lab == cl),
              numeric(n))
  F_scores <- (1 - xi) * solve(diag(n) - xi * S, Y)
  dimnames(F_scores) <- list(nodes, classes)

  win <- max.col(F_scores, ties.method = "first")
  score <- F_scores[cbind(seq_len(n), win)]
  runner <- vapply(seq_len(n), function(i) {
    max(F_scores[i, -win[i]])
  }, numeric(1L))

  predictions <- data.frame(
    node      = nodes,
    label     = lab,
    predicted = classes[win],
    score     = score,
    margin    = score - runner,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      predictions = predictions,
      classes     = classes,
      scores      = F_scores,
      xi          = xi,
      type        = type,
      n_labeled   = sum(!is.na(lab)),
      n_nodes     = n,
      params      = list(edge_weights = parts$w)
    ),
    class = "net_hypergraph_transduction"
  )
}

# ---- Internal machinery -------------------------------------------------

#' Validate a net_hypergraph for spectral work (connected, non-degenerate)
#' @noRd
.hl_validate_hg <- function(hg) {
  .assert(
    "`hg` must be a net_hypergraph (build_hypergraph/bipartite_groups)" =
      inherits(hg, "net_hypergraph"),
    "`hg` must have at least 2 nodes" = hg$n_nodes >= 2L,
    "`hg` must have at least 1 hyperedge" = hg$n_hyperedges >= 1L
  )
  if (!.hl_is_connected(hg$incidence > 0)) {
    stop(errorCondition(
      paste0("The hypergraph is not connected; the random walk has no ",
             "unique stationary distribution. Analyze components ",
             "separately."),
      class = "nestimate_hypergraph_disconnected", call = NULL
    ))
  }
  invisible(TRUE)
}

#' Connectivity of the hypergraph via BFS on the node co-membership graph
#' @noRd
.hl_is_connected <- function(pattern) {
  n <- nrow(pattern)
  if (n == 1L) return(TRUE)
  adj <- tcrossprod(pattern * 1.0) > 0
  reached <- c(TRUE, rep(FALSE, n - 1L))
  # BFS frontier expansion; bounded by n iterations (diameter <= n - 1)
  repeat {
    nxt <- reached | (as.vector(adj %*% reached) > 0)
    if (identical(nxt, reached)) break
    reached <- nxt
  }
  all(reached)
}

#' Population standard deviation (ddof = 0, numpy convention)
#' @noRd
.hl_pop_sd <- function(x) sqrt(mean((x - mean(x))^2))

#' Build Laplacian + stationary distribution + edge weights for a type
#' @return list(L, pi, w)
#' @noRd
.hl_build <- function(hg, type, edge_weights) {
  gamma <- hg$incidence * 1.0
  pattern <- (gamma > 0) * 1.0
  n <- nrow(gamma)
  m <- ncol(gamma)
  nodes <- rownames(gamma) %||% hg$nodes

  if (!is.null(edge_weights)) {
    .assert(
      "`edge_weights` must be a positive numeric vector, one per hyperedge" =
        is.numeric(edge_weights) && length(edge_weights) == m &&
        all(is.finite(edge_weights)) && all(edge_weights > 0)
    )
  }

  if (type == "zhou") {
    # Binary incidence pattern; default unit hyperedge weights
    w <- edge_weights %||% rep(1, m)
    delta_e <- colSums(pattern)
    d_v <- as.vector(pattern %*% w)
    Hs <- pattern / sqrt(d_v)
    Theta <- Hs %*% ((w / delta_e) * t(Hs))
    L <- diag(n) - (Theta + t(Theta)) / 2
    pi_v <- d_v / sum(d_v)
  } else {
    # Hayashi EDVW random walk on the weighted incidence
    w <- edge_weights %||% vapply(seq_len(m), function(j) {
      .hl_pop_sd(gamma[gamma[, j] > 0, j]) + 1
    }, numeric(1L))
    delta_e <- colSums(gamma)
    d_v <- as.vector(pattern %*% w)
    A <- sweep(pattern, 2L, w, "*") / d_v      # A[v,e] = w(e) 1[v in e]/d(v)
    B <- t(gamma) / delta_e                     # B[e,u] = gamma_e(u)/delta(e)
    P <- A %*% B

    # Stationary distribution: dominant left eigenvector of P
    eg <- eigen(t(P))
    lead <- which.max(Mod(eg$values))
    pi_v <- Re(eg$vectors[, lead])
    pi_v <- pi_v / sum(pi_v)
    if (min(pi_v) < -1e-8) {
      stop("Stationary distribution has negative mass; the walk appears ",
           "reducible despite the connectivity check.", call. = FALSE)
    }
    pi_v <- pmax(pi_v, 0)
    pi_v <- pi_v / sum(pi_v)

    # Chung (2005) directed Laplacian, symmetrized via pi
    G <- (P * sqrt(pi_v))                       # scale rows by sqrt(pi)
    G <- sweep(G, 2L, sqrt(pi_v), "/")          # scale cols by 1/sqrt(pi)
    L <- diag(n) - (G + t(G)) / 2
  }

  dimnames(L) <- list(nodes, nodes)
  names(pi_v) <- nodes
  list(L = L, pi = pi_v, w = w)
}

# ---- S3: net_hypergraph_cluster ----------------------------------------

#' @rdname hypergraph_cluster
#' @param x For the `print()`, `as.data.frame()` and `plot()` methods: a
#'   `net_hypergraph_cluster` object.
#' @param ... For the methods: ignored.
#' @param object For the `summary()` method: a `net_hypergraph_cluster`
#'   object.
#' @param what For the `plot()` method: the panel to draw. `"spectrum"` is
#'   a scree plot of the smallest eigenvalues with the `k` used marked, so
#'   the eigengap after `k` can be read; `"embedding"` places the nodes on
#'   the first two embedding dimensions, labelled, coloured and shaped by
#'   cluster and sized by stationary probability, which is the geometry
#'   k-means clustered; `"both"` (default) draws the two side by side.
#' @param n_values For the `plot()` method: how many of the smallest
#'   eigenvalues the spectrum panel shows. Default `NULL`, which shows
#'   `min(3 * k, n_nodes)`.
#' @return The `print()` method returns its input invisibly. The `summary()`
#'   method returns a data frame with one row per cluster: `cluster`, `size`
#'   and `share` (the proportion of nodes). The `as.data.frame()` method
#'   returns a data frame with one row per node: `node`, `cluster`, `pi`
#'   (the node's stationary probability under the Laplacian's walk) and the
#'   embedding coordinates `dim1` to `dimk`. The `plot()` method returns a
#'   `ggplot` for `"spectrum"` or `"embedding"`. For `"both"` it draws the
#'   two panels and returns the arranged `gtable` invisibly when
#'   \pkg{gridExtra} is installed, and otherwise draws them in grid
#'   viewports and returns the two `ggplot`s as a list, invisibly.
#' @export
print.net_hypergraph_cluster <- function(x, ...) {
  cat("Hypergraph spectral clustering (", x$type, " Laplacian)\n", sep = "")
  cat(sprintf("  Nodes: %d | Hyperedges: %d | k: %d\n",
              x$n_nodes, x$n_hyperedges, x$k))
  cat(sprintf("  Cluster sizes: %s\n",
              paste(sprintf("%s = %d", x$sizes$cluster, x$sizes$size),
                    collapse = ", ")))
  cat(sprintf("  Eigengap after k: %.4f\n", x$eigengap))
  invisible(x)
}

#' @rdname hypergraph_cluster
#' @export
summary.net_hypergraph_cluster <- function(object, ...) {
  out <- object$sizes
  out$share <- out$size / sum(out$size)
  out
}

#' @rdname hypergraph_cluster
#' @export
as.data.frame.net_hypergraph_cluster <- function(x, ...) {
  out <- x$clusters
  out$pi <- as.numeric(x$pi)
  cbind(out, as.data.frame(x$embedding), row.names = NULL)
}

#' @rdname hypergraph_cluster
#' @export
plot.net_hypergraph_cluster <- function(x,
                                        what = c("both", "spectrum",
                                                 "embedding"),
                                        n_values = NULL, ...) {
  what <- match.arg(what)
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#999999", "#000000")

  n_show <- as.integer(n_values %||% min(3L * x$k, x$n_nodes))
  df_s <- data.frame(
    index = seq_len(n_show),
    value = x$eigenvalues[seq_len(n_show)],
    used  = ifelse(seq_len(n_show) <= x$k, "used", "unused")
  )
  p_spec <- ggplot2::ggplot(df_s, ggplot2::aes(x = .data$index,
                                               y = .data$value,
                                               color = .data$used,
                                               shape = .data$used)) +
    ggplot2::geom_line(color = "#999999", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = x$k + 0.5, linetype = "dashed",
                        color = "#999999") +
    ggplot2::scale_color_manual(
      values = c(used = "#0072B2", unused = "#E69F00")
    ) +
    ggplot2::scale_shape_manual(values = c(used = 16, unused = 17)) +
    ggplot2::labs(
      title = sprintf("Laplacian spectrum (k = %d)", x$k),
      subtitle = sprintf("Eigengap after k: %.4f", x$eigengap),
      x = "Eigenvalue index (increasing)", y = "Eigenvalue",
      color = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (what == "spectrum") return(p_spec)

  df_e <- as.data.frame(x)
  p_emb <- ggplot2::ggplot(df_e, ggplot2::aes(x = .data$dim1,
                                              y = .data$dim2,
                                              color = .data$cluster,
                                              shape = .data$cluster)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$pi), alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = .data$node),
                       vjust = -1.1, size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = okabe) +
    ggplot2::scale_shape_manual(
      values = rep(c(16, 17, 15, 18, 8, 7, 3, 4, 6), length.out = x$k)
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6),
                                   guide = "none") +
    ggplot2::expand_limits(
      x = range(df_e$dim1) + c(-0.15, 0.15) * diff(range(df_e$dim1)),
      y = range(df_e$dim2) + c(-0.1, 0.25) * diff(range(df_e$dim2))
    ) +
    ggplot2::labs(
      title = "Spectral embedding (first two dimensions)",
      subtitle = "Point size = stationary probability of the node",
      x = "dim1", y = "dim2", color = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (what == "embedding") return(p_emb)

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(p_spec, p_emb, ncol = 2))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, 2L)))
  print(p_spec, vp = grid::viewport(layout.pos.row = 1L,
                                    layout.pos.col = 1L))
  print(p_emb, vp = grid::viewport(layout.pos.row = 1L,
                                   layout.pos.col = 2L))
  grid::popViewport()
  invisible(list(spectrum = p_spec, embedding = p_emb))
}

# ---- S3: net_hypergraph_transduction -----------------------------------

#' @rdname hypergraph_transduction
#' @param x For the `print()`, `as.data.frame()` and `plot()` methods: a
#'   `net_hypergraph_transduction` object.
#' @param ... For the methods: ignored.
#' @param object For the `summary()` method: a `net_hypergraph_transduction`
#'   object.
#' @param row.names For the `as.data.frame()` method: row names for the
#'   returned data frame. Default `NULL` keeps the default row numbers.
#' @param optional For the `as.data.frame()` method: ignored; it is there to
#'   match the [base::as.data.frame()] generic.
#' @param what For the `as.data.frame()` method: the table to return.
#'   `"predictions"` (default) gives one row per node; `"scores"` gives one
#'   row per node and class.
#' @return The `print()` method returns its input invisibly. The `summary()`
#'   method returns a data frame with one row per class: `class`,
#'   `n_labeled`, `n_predicted` and `mean_margin` (the mean margin of the
#'   nodes predicted into the class, `NA` when there are none). The
#'   `as.data.frame()` method returns, for `what = "predictions"`, one row
#'   per node with `node`, `label` (`NA` for an unlabelled node),
#'   `predicted`, `score` and `margin`, and for `what = "scores"` one row per
#'   node and class with `node`, `class` and `score`. The `plot()` method
#'   returns the score heatmap as a `ggplot`, visibly, so `plot(x)` draws
#'   it: one row per node, grouped by
#'   predicted class and ordered by margin within it, one column per class,
#'   the tiles shaded and labelled by score. Labelled nodes have a black
#'   border on their given class, and a dot marks each node's winning class,
#'   so seeds, scores and assignments can be compared in one panel.
#' @export
print.net_hypergraph_transduction <- function(x, ...) {
  cat("Hypergraph transductive label spreading (", x$type,
      " Laplacian, xi = ", format(x$xi), ")\n", sep = "")
  tab <- table(x$predictions$predicted)
  cat(sprintf("  Nodes: %d (%d labeled) | Classes: %s\n",
              x$n_nodes, x$n_labeled, paste(x$classes, collapse = ", ")))
  cat(sprintf("  Predicted: %s\n",
              paste(sprintf("%s = %d", names(tab), as.integer(tab)),
                    collapse = ", ")))
  invisible(x)
}

#' @rdname hypergraph_transduction
#' @export
summary.net_hypergraph_transduction <- function(object, ...) {
  p <- object$predictions
  out <- do.call(rbind, lapply(object$classes, function(cl) {
    sel <- p$predicted == cl
    data.frame(
      class       = cl,
      n_labeled   = sum(!is.na(p$label) & p$label == cl),
      n_predicted = sum(sel),
      mean_margin = if (any(sel)) mean(p$margin[sel]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

#' @rdname hypergraph_transduction
#' @export
as.data.frame.net_hypergraph_transduction <- function(
    x, row.names = NULL, optional = FALSE,
    what = c("predictions", "scores"), ...) {
  what <- match.arg(what)
  out <- if (what == "predictions") {
    x$predictions
  } else {
    data.frame(
      node  = rep(rownames(x$scores), times = ncol(x$scores)),
      class = rep(colnames(x$scores), each = nrow(x$scores)),
      score = as.vector(x$scores),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(row.names)) rownames(out) <- row.names
  out
}

#' @rdname hypergraph_transduction
#' @export
plot.net_hypergraph_transduction <- function(x, ...) {
  sc <- as.data.frame(x, what = "scores")
  p <- x$predictions
  sc$seeded <- !is.na(p$label[match(sc$node, p$node)]) &
    p$label[match(sc$node, p$node)] == sc$class
  sc$winner <- p$predicted[match(sc$node, p$node)] == sc$class
  # order rows by predicted class, then margin (most confident on top)
  ord <- p$node[order(p$predicted, -p$margin)]
  sc$node <- factor(sc$node, levels = rev(ord))

  ggplot2::ggplot(sc, ggplot2::aes(x = .data$class, y = .data$node)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$score),
                       color = "white", linewidth = 0.6) +
    ggplot2::geom_tile(
      data = sc[sc$seeded, , drop = FALSE],
      fill = NA, color = "#000000", linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = sc[sc$winner, , drop = FALSE],
      ggplot2::aes(x = .data$class, y = .data$node),
      shape = 21, size = 2.6, stroke = 0.9,
      fill = "white", color = "#D55E00",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$score)),
                       size = 3, vjust = -0.9) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#4A6FE3") +
    ggplot2::labs(
      title = "Transductive label spreading: score matrix",
      subtitle = sprintf(paste0("%d labeled of %d nodes, xi = %s. Black ",
                                "border = seed; dot = winning class."),
                         x$n_labeled, x$n_nodes, format(x$xi)),
      x = "Class", y = NULL, fill = "Score"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid = ggplot2::element_blank())
}
