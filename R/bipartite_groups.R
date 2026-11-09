# ---- Bipartite group hypergraph (EG-6) -----------------------------------
# Direct constructor: long-format event data with player + group columns
# becomes a net_hypergraph where each group is a hyperedge spanning all
# players that appeared in it.

#' Build a Hypergraph from Group Memberships
#'
#' @description
#' `bipartite_groups()` builds a hypergraph from long-format membership
#' data, in which each row records that a `player` took part in a `group`
#' (a session, team, project, transaction or any other group context), and
#' returns a `net_hypergraph`. Each distinct group becomes one hyperedge
#' spanning the players that appeared in it. The incidence matrix of
#' players by groups is 0/1, or holds summed weights when `weight` is
#' given.
#'
#' @param data A data frame with one row per membership (a player in a
#'   group), holding the `player`, `group` and, when used, `weight`
#'   columns. Rows with a missing value in any of these columns are dropped
#'   without a message; when no row is left, an error is raised.
#' @param player The name of the column whose values become the nodes
#'   (players, participants, actors). Values are read as character and the
#'   nodes are sorted. A name that is not a single column of `data` raises
#'   an error, as for `group` and `weight`.
#' @param group The name of the column whose values become the hyperedges
#'   (groups, sessions, teams). Values are read as character, and the
#'   hyperedges follow their sorted order (as text, so `"10"` comes before
#'   `"2"`).
#' @param weight The name of a numeric column summed per player and group
#'   to give a weighted incidence matrix. With the default `NULL`, the
#'   matrix is 0/1, and a player listed more than once in a group counts
#'   once. With weights, a player belongs to a group when the summed weight
#'   is positive; a group with no such member is dropped, and a player left
#'   without a group stays a node.
#'
#' @return A `net_hypergraph`, the class [build_hypergraph()] returns. Print
#'   it for the numbers of nodes and hyperedges and the hyperedge sizes;
#'   `summary()` prints the counts with the mean and largest hyperedge size,
#'   and returns a data frame of node degrees (`node`, and `degree`, the
#'   number of hyperedges a player belongs to). [hypergraph_measures()],
#'   [hypergraph_centrality()], [hypergraph_laplacian()] and
#'   [hypergraph_cluster()] read it. The fields:
#' \describe{
#'   \item{hyperedges}{An unnamed list with one vector per hyperedge, in
#'     the sorted order of the groups: the indices of its players in
#'     `nodes`, sorted and named by player.}
#'   \item{incidence}{The players-by-groups matrix, one column per
#'     hyperedge, with player and group names as dimnames: integer 0/1, or
#'     the summed weights (a player whose sum is zero or negative keeps
#'     that sum here but is not a member of the hyperedge).}
#'   \item{nodes}{The player names, sorted.}
#'   \item{n_nodes, n_hyperedges}{The numbers of players and of hyperedges
#'     (after dropping empty groups).}
#'   \item{size_distribution}{A named integer vector with the number of
#'     hyperedges of each size (`size_1`, `size_2`, ...).}
#'   \item{params}{A list with `source = "bipartite_groups"`, the column
#'     names given as `player`, `group` and `weight`, and
#'     `n_observations`, the number of rows used.}
#' }
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
#'   \item `bipartite_groups()` - when group membership is observed
#'     (sessions, transactions, co-authorships).
#'   \item `build_hypergraph()` - when only pairwise interactions are
#'     observed and group structure must be inferred from cliques
#'     (triangles by default).
#' }
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
#' @note (experimental) Validated against a hand-computed `table()` incidence
#'   reference only; no independent R package exposes the
#'   long-format-to-binary-incidence primitive, because the operation is
#'   definitionally `table()`. The code path is a direct one-to-one
#'   restatement of its definition.
#'
#' @export
bipartite_groups <- function(data, player, group, weight = NULL) {
  .assert(
    "`data` must be a data frame" = is.data.frame(data),
    "`player` must name one column of `data`" =
      is.character(player) && length(player) == 1L && player %in% names(data),
    "`group` must name one column of `data`" =
      is.character(group) && length(group) == 1L && group %in% names(data),
    "`weight` must be NULL or name one column of `data`" =
      is.null(weight) ||
        (is.character(weight) && length(weight) == 1L &&
           weight %in% names(data))
  )

  cols <- c(player, group, weight)
  d <- data[, cols, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0L) {
    stop(errorCondition("No complete observations after dropping NAs.",
                        class = "nestimate_bad_input", call = NULL))
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
    # the weights of each (player, group) cell, added in row order
    cell <- (gj - 1L) * n_players + pi
    totals <- vapply(split(w, cell), function(v) Reduce(`+`, v, 0),
                     numeric(1L))
    incidence[as.integer(names(totals))] <- totals
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
