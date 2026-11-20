# Inputs and reference computations shared by the hypergraph equivalence tests, written in base R. Nothing here
# calls Nestimate.

#' A random long table of memberships: `n_groups` groups drawing members from `n_people` people
#'
#' @param n_people,n_groups People and groups.
#' @param max_size Largest group.
#' @param seed Seed of the draw.
#' @param weighted Draw a membership weight in (0.5, 3) instead of 1.
#' @return A data frame `person`, `group`, `w`, one row per distinct membership.
hg_memberships <- function(n_people, n_groups, max_size, seed, weighted = FALSE) {
  set.seed(seed)
  rows <- lapply(seq_len(n_groups), \(group) {
    members <- sample(n_people, sample(2:max_size, 1L))
    data.frame(person = sprintf("p%02d", members), group = sprintf("g%02d", group),
               w = if (weighted) round(stats::runif(length(members), 0.5, 3), 3) else 1)
  })
  do.call(rbind, rows)
}

#' The 0/1 incidence matrix of a membership table: people in rows, groups in columns, both sorted
#'
#' @param memberships A data frame `person`, `group`.
#' @return An integer matrix.
hg_incidence <- function(memberships) {
  counts <- table(memberships$person, memberships$group)
  incidence <- matrix(as.integer(counts > 0), nrow(counts), ncol(counts), dimnames = dimnames(counts))
  incidence[sort(rownames(incidence)), sort(colnames(incidence)), drop = FALSE]
}

#' The weighted incidence matrix of a membership table (the membership weight in each cell)
#'
#' @param memberships A data frame `person`, `group`, `w`.
#' @return A numeric matrix, people by groups.
hg_weighted_incidence <- function(memberships) {
  weights <- tapply(memberships$w, list(memberships$person, memberships$group), sum)
  weights[is.na(weights)] <- 0
  weights[sort(rownames(weights)), sort(colnames(weights)), drop = FALSE]
}

#' A square matrix with its diagonal set to zero
#'
#' @param m A square matrix.
#' @return The matrix with a zero diagonal.
hg_offdiag <- function(m) {
  diag(m) <- 0
  m
}

#' Zhou, Huang and Schoelkopf (2006): Theta = Dv^-1/2 H W De^-1 H' Dv^-1/2 of a 0/1 incidence matrix, unit edge weights
#'
#' @param incidence A 0/1 matrix, nodes by hyperedges.
#' @return A list: `theta` (n by n), `laplacian` (I - theta), `pi` (d_v / sum d_v).
hg_zhou <- function(incidence) {
  degree_v <- rowSums(incidence)
  degree_e <- colSums(incidence)
  half <- diag(1 / sqrt(degree_v), length(degree_v))
  theta <- half %*% incidence %*% diag(1 / degree_e, length(degree_e)) %*% t(incidence) %*% half
  dimnames(theta) <- list(rownames(incidence), rownames(incidence))
  list(theta = theta, laplacian = diag(nrow(theta)) - theta, pi = degree_v / sum(degree_v))
}

#' A partition as one row per unordered pair of nodes: 1 when the two share a cluster
#'
#' @param nodes Node names.
#' @param cluster Cluster labels, one per node (any labels).
#' @return A data frame `from`, `to`, `same`.
hg_same_cluster <- function(nodes, cluster) {
  pairs <- utils::combn(order(nodes), 2L)
  data.frame(from = nodes[pairs[1L, ]], to = nodes[pairs[2L, ]],
             same = as.integer(cluster[pairs[1L, ]] == cluster[pairs[2L, ]]))
}

#' Every clique of a graph (the simplices of its clique complex), found by testing each subset with combn()
#'
#' @param adjacency A symmetric 0/1 matrix with names.
#' @return A list of character vectors (sorted node names), one per clique, single nodes included.
hg_cliques <- function(adjacency) {
  nodes <- rownames(adjacency)
  grow <- \(cliques, size) {
    candidates <- utils::combn(nodes, size, simplify = FALSE)
    found <- Filter(\(s) all(adjacency[s, s][upper.tri(diag(size))] > 0), candidates)
    if (length(found) == 0L) return(cliques)
    grow(c(cliques, found), size + 1L)
  }
  grow(as.list(nodes), 2L)
}

#' Betti numbers over the rationals of a simplicial complex: b_k = dim C_k - rank d_k - rank d_(k+1)
#'
#' @param simplices A list of character vectors, closed under faces.
#' @return An integer vector b0, b1, ... up to the top dimension.
hg_betti <- function(simplices) {
  keys <- vapply(simplices, paste, character(1), collapse = ",")
  dims <- lengths(simplices) - 1L
  top <- max(dims)
  boundary_rank <- \(k) {
    if (k < 1L || k > top) return(0L)
    rows <- which(dims == k - 1L)
    cols <- which(dims == k)
    entries <- vapply(cols, \(col) {
      faces <- vapply(seq_along(simplices[[col]]), \(drop) paste(simplices[[col]][-drop], collapse = ","), character(1))
      signs <- (-1)^(seq_along(faces) - 1L)
      column <- numeric(length(rows))
      column[match(faces, keys[rows])] <- signs
      column
    }, numeric(length(rows)))
    qr(matrix(entries, nrow = length(rows)))$rank
  }
  vapply(0:top, \(k) as.integer(sum(dims == k) - boundary_rank(k) - boundary_rank(k + 1L)), integer(1))
}
