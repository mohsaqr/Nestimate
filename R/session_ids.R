#' The session behind each sequence
#'
#' @description
#' Names every sequence of a network or of a clustering by the columns it was
#' built from. A network built from long data with \code{actor} and
#' \code{session} has one sequence per actor-session; its rows are ordered by
#' the grouping, not by the input, and a fitted clustering reports its
#' assignments in that same order. \code{session_ids()} returns the key that
#' joins them back to the input data, so no label has to be parsed.
#'
#' @param x A \code{netobject} from \code{\link{build_network}} on long data,
#'   or a \code{net_mmm} (\code{\link{build_mmm}}) or \code{net_clustering}
#'   (\code{\link{build_clusters}}) fitted on such a network.
#' @param ... Unused.
#'
#' @return A data frame with one row per sequence, in the row order of the
#'   network's \code{$data} (and of the fit's assignments):
#'   \describe{
#'     \item{sequence}{Integer row number of the sequence.}
#'     \item{actor and session columns}{The \code{actor} and \code{session}
#'       columns given to \code{build_network()}, under their own names and
#'       with their own values.}
#'     \item{session_label}{The readable label of the sequence. With
#'       \code{time}, sessions split at time gaps carry a \code{" s<n>"}
#'       suffix, so this column separates them.}
#'     \item{cluster}{For a \code{net_mmm} or \code{net_clustering}: the
#'       assigned cluster (integer).}
#'     \item{posterior}{For a \code{net_mmm}: the posterior probability of
#'       the assigned cluster.}
#'   }
#'
#' @section Errors:
#' Raises \code{nestimate_no_session_ids} when \code{x} carries no
#' per-sequence metadata: a network built from wide data, a fit on wide data
#' or on a \code{tna} model, or a fit made before Nestimate 0.9.6 (refit it).
#' Raises \code{nestimate_session_ids_misaligned} when the metadata and the
#' sequences differ in number.
#'
#' @examples
#' events <- data.frame(
#'   student = rep(c("s1", "s2", "s3"), each = 8),
#'   step    = rep(c("a", "b"), each = 4, times = 3),
#'   action  = sample(c("read", "write", "test"), 24, replace = TRUE)
#' )
#' net <- build_network(events, actor = "student", session = "step",
#'                      action = "action", method = "relative")
#' session_ids(net)
#'
#' \donttest{
#' fit <- build_mmm(net, k = 2, n_starts = 2, seed = 1)
#' session_ids(fit)
#' }
#' @seealso \code{\link{build_network}}, \code{\link{build_mmm}},
#'   \code{\link{build_clusters}}
#' @export
session_ids <- function(x, ...) {
  UseMethod("session_ids")
}

#' @rdname session_ids
#' @export
session_ids.default <- function(x, ...) {
  stop(errorCondition(
    sprintf("session_ids() needs a netobject, net_mmm or net_clustering, not a %s.",
            class(x)[1L]),
    class = "nestimate_no_session_ids", call = NULL))
}

#' @rdname session_ids
#' @export
session_ids.netobject <- function(x, ...) {
  .session_id_table(x$metadata, x$build_args, nrow(x$data))
}

#' @rdname session_ids
#' @export
session_ids.net_mmm <- function(x, ...) {
  assigned <- x$assignments
  out <- .session_id_table(x$metadata, x$build_args, length(assigned))
  out$cluster <- as.integer(assigned)
  out$posterior <- x$posterior[cbind(seq_along(assigned), assigned)]
  out
}

#' @rdname session_ids
#' @export
session_ids.net_clustering <- function(x, ...) {
  assigned <- x$assignments
  out <- .session_id_table(x$metadata, x$build_args, length(assigned))
  out$cluster <- unname(as.integer(assigned))
  out
}

# One row per sequence: its row number, the actor and session columns it was
# built from, and its label. `metadata` must be in sequence row order, which
# prepare() guarantees.
.session_id_table <- function(metadata, build_args, n) {
  if (is.null(metadata) || !".session_label" %in% names(metadata)) {
    stop(errorCondition(
      paste0("No session ids: the sequences were not built from long data with ",
             "`actor`/`session`, or the fit predates Nestimate 0.9.6 (refit it)."),
      class = "nestimate_no_session_ids", call = NULL))
  }
  if (nrow(metadata) != n) {
    stop(errorCondition(
      sprintf("The metadata has %d rows but there are %d sequences.", nrow(metadata), n),
      class = "nestimate_session_ids_misaligned", call = NULL))
  }
  id_cols <- intersect(c(build_args$actor, build_args$session), names(metadata))
  out <- data.frame(sequence = seq_len(n))
  if (length(id_cols) > 0L) {
    ids <- metadata[id_cols]
    # columns aggregated per session by tapply() are 1-d named arrays
    ids[] <- lapply(ids, function(v) {
      dim(v) <- NULL
      unname(v)
    })
    out <- cbind(out, ids)
  }
  out$session_label <- unname(metadata$.session_label)
  rownames(out) <- NULL
  out
}
