# ---- Bipartite group hypergraph (EG-6) -----------------------------------
# Direct constructor: long-format event data with player + group columns
# becomes a net_hypergraph where each group is a hyperedge spanning all
# players that appeared in it.

#' Hypergraph from bipartite group / event data
#'
#' Constructs a [net_hypergraph][build_hypergraph] from long-format event
#' data in which each row records a `player` participating in a `group`
#' (a session, team, project, transaction, or any group context). Each
#' unique group becomes one hyperedge spanning the players that appeared in
#' it. Optional `weight` column produces a weighted incidence matrix.
#'
#' @param data Data frame in long format. Must contain `player` and `group`
#'   columns; optionally a `weight` column.
#' @param player Character. Name of the column whose values become the
#'   hypergraph's nodes (players, participants, actors).
#' @param group Character. Name of the column whose values become the
#'   hypergraph's hyperedges (groups, sessions, teams).
#' @param weight Character or `NULL`. If supplied, the column is summed per
#'   `(player, group)` pair to produce a weighted incidence matrix. Default
#'   `NULL` produces a 0/1 binary incidence matrix.
#'
#' @return A `net_hypergraph` object with the same structure produced by
#'   [build_hypergraph()] (`hyperedges`, `incidence`, `nodes`, `n_nodes`,
#'   `n_hyperedges`, `size_distribution`, `params`). The `params` list
#'   records `source = "bipartite_groups"` and the original column names.
#'
#' @details
#' The bipartite representation preserves the full group structure without
#' projecting to a pairwise network. A group of three players A, B, C
#' produces a single 3-hyperedge containing all three, not three pairwise
#' edges AB, AC, BC. This avoids information loss when group interactions are
#' the primary unit of analysis (Perc et al. 2013).
#'
#' Unlike [build_hypergraph()] (which derives hyperedges from a network's
#' clique structure), `bipartite_groups()` takes group memberships
#' directly. The two functions are complementary:
#' \itemize{
#'   \item `bipartite_groups()` — when group membership is observed
#'     (sessions, transactions, co-authorships).
#'   \item `build_hypergraph()` — when only pairwise interactions are
#'     observed and triadic structure must be inferred from triangles.
#' }
#'
#' Rows with `NA` in either the `player` or `group` column (or, when
#' supplied, the `weight` column) are dropped silently.
#'
#' @seealso [build_hypergraph()] for the clique-based constructor.
#'
#' @examples
#' df <- data.frame(
#'   player = c("Alice", "Bob", "Carol", "Alice", "Bob",
#'              "Dave", "Carol", "Dave", "Eve"),
#'   session = c("S1", "S1", "S1", "S2", "S2",
#'               "S3", "S3", "S3", "S3")
#' )
#' hg <- bipartite_groups(df, player = "player", group = "session")
#' print(hg)
#' summary(hg)
#'
#' @references
#' Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
#' (2013). Evolutionary dynamics of group interactions on structured
#' populations: a review. \emph{Journal of the Royal Society Interface}
#' 10(80), 20120997. \doi{10.1098/rsif.2012.0997}
#'
#' @export
bipartite_groups <- function(data, player, group, weight = NULL) {
  stopifnot(
    is.data.frame(data),
    is.character(player), length(player) == 1L,
    is.character(group),  length(group)  == 1L,
    player %in% names(data),
    group  %in% names(data),
    is.null(weight) ||
      (is.character(weight) && length(weight) == 1L && weight %in% names(data))
  )

  cols <- c(player, group, weight)
  d <- data[, cols, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No complete observations after dropping NAs.", call. = FALSE)
  }

  d[[player]] <- as.character(d[[player]])
  d[[group]]  <- as.character(d[[group]])

  player_levels <- sort(unique(d[[player]]))
  group_levels  <- sort(unique(d[[group]]))
  n_players <- length(player_levels)
  n_groups  <- length(group_levels)

  # Map values to row/col indices, then accumulate via tabulate
  pi <- match(d[[player]], player_levels)
  gj <- match(d[[group]],  group_levels)
  if (is.null(weight)) {
    cell <- (gj - 1L) * n_players + pi
    counts <- tabulate(cell, nbins = n_players * n_groups)
    incidence <- matrix(as.integer(counts > 0L), n_players, n_groups,
                        dimnames = list(player_levels, group_levels))
  } else {
    incidence <- matrix(0, n_players, n_groups,
                        dimnames = list(player_levels, group_levels))
    w <- as.numeric(d[[weight]])
    for (k in seq_along(pi)) {
      incidence[pi[k], gj[k]] <- incidence[pi[k], gj[k]] + w[k]
    }
  }

  # Drop hyperedges that ended up empty (e.g. all-zero weight)
  he_sizes_pre <- colSums(incidence > 0)
  keep <- he_sizes_pre > 0
  incidence <- incidence[, keep, drop = FALSE]
  group_levels <- group_levels[keep]
  n_groups <- length(group_levels)

  hyperedges <- lapply(seq_len(n_groups), function(j) {
    sort(which(incidence[, j] > 0))
  })

  he_sizes <- vapply(hyperedges, length, integer(1L))
  size_dist <- if (length(he_sizes)) {
    tab <- table(he_sizes)
    out <- as.integer(tab)
    names(out) <- paste0("size_", names(tab))
    out
  } else {
    integer(0L)
  }

  structure(
    list(
      hyperedges        = hyperedges,
      incidence         = incidence,
      nodes             = player_levels,
      n_nodes           = n_players,
      n_hyperedges      = n_groups,
      size_distribution = size_dist,
      params = list(
        source         = "bipartite_groups",
        player         = player,
        group          = group,
        weight         = weight,
        n_observations = nrow(d)
      )
    ),
    class = "net_hypergraph"
  )
}
