# ==============================================================================
# extract_pathways() -- cut an event log into pathways, one row per pathway.
#
# Two shapes of pathway, one verb:
#
#   * anchored -- every occurrence of `anchor` opens a pathway that closes at
#     the next `terminal` state (a failure and its consequence);
#   * whole-unit -- each grouping unit is one pathway (an episode run to
#     whatever ends it).
#
# The anchored form is the expensive one to write by hand: the naive version
# loops over every anchor position. Here the match is a findInterval() over
# the terminal positions, so the whole log is cut in one pass.
# ==============================================================================

utils::globalVariables(c("..pid", "..pos"))

#' Cut an Event Log into Pathways
#'
#' Splits a long event log into pathways and returns them as a tidy
#' \code{data.frame}, one row per pathway, ready to tabulate against an
#' outcome or to feed to a sequence model.
#'
#' @param data A long-format \code{data.frame}, one row per event, already in
#'   chronological order within each group (or with \code{order} supplied).
#' @param action Name of the column holding the state / event label.
#' @param group Character vector of column names defining the unit a pathway
#'   may not cross -- typically the actor, or actor and session.
#' @param order Optional column giving the within-group event order. When
#'   \code{NULL} (default) the existing row order is used.
#' @param type The cut to make:
#'   \describe{
#'     \item{\code{"unit"}}{(default) one pathway per group. With
#'       \code{terminal}, the group is truncated at its last terminal state
#'       so the pathway ends on an outcome rather than on whatever followed.}
#'     \item{\code{"segments"}}{consecutive, non-overlapping pathways: the
#'       group is cut after every terminal state, so each pathway ends on one.
#'       Every event belongs to exactly one pathway.}
#'     \item{\code{"anchored"}}{one pathway per occurrence of \code{anchor},
#'       closing at the next terminal state. Pathways may overlap and events
#'       between a terminal and the next anchor belong to none.}
#'   }
#' @param anchor State that opens a pathway when \code{type = "anchored"}.
#'   Ignored by the other cuts. Default \code{NULL}.
#' @param terminal States that close a pathway, the closing state included.
#'   Required for \code{"segments"} and \code{"anchored"}; optional for
#'   \code{"unit"}, where it truncates.
#' @param resolve Optional named list giving a resolution label to append as
#'   the pathway's final state. Each element is a character vector of states;
#'   a pathway takes the label of the \strong{first} element whose states all
#'   occur in it, so order the list from most specific to least. Pathways
#'   matching nothing are labelled \code{"Unresolved"}. The label becomes the
#'   last element of \code{path} and the value of \code{closes}, and the raw
#'   final state is kept in \code{ends}.
#' @param sep Separator used when pasting the path. Default \code{" -> "}.
#'
#' @return A \code{data.frame} with one row per pathway: the \code{group}
#'   columns, \code{pathway} (an integer id within group), \code{path} (the
#'   states pasted with \code{sep}), \code{length}, \code{opens} (first state)
#'   and \code{closes} (last state). With \code{resolve}, \code{closes} holds
#'   the resolution label instead and a further \code{ends} column carries the
#'   raw final state. Pathways are returned in the order they occur.
#' @seealso \code{\link{sequence_compare}}, \code{\link{outcome_model}}
#' @examples
#' log <- data.frame(
#'   id  = c(1, 1, 1, 1, 1, 1, 2, 2, 2),
#'   act = c("Try", "Wrong", "Hint", "Retry", "Right", "Praise",
#'           "Try", "Wrong", "Retry"),
#'   stringsAsFactors = FALSE
#' )
#' # Each failure and its consequence:
#' extract_pathways(log, action = "act", group = "id", type = "anchored",
#'                  anchor = "Wrong", terminal = c("Right", "Wrong"))
#'
#' # Each id as one episode:
#' extract_pathways(log, action = "act", group = "id")
#' @export
extract_pathways <- function(data, action, group, order = NULL,
                             type = c("unit", "segments", "anchored"),
                             anchor = NULL, terminal = NULL, resolve = NULL,
                             sep = " -> ") {
  type <- match.arg(type)
  stopifnot(
    "`data` must be a data.frame" = is.data.frame(data),
    "`action` must be a single column name" =
      is.character(action) && length(action) == 1L,
    "`group` must be a character vector of column names" =
      is.character(group) && length(group) >= 1L,
    "`anchor` must be a single state or NULL" =
      is.null(anchor) || length(anchor) == 1L,
    "`terminal` must be a character vector or NULL" =
      is.null(terminal) || is.character(terminal),
    "`sep` must be a single string" = is.character(sep) && length(sep) == 1L,
    "`resolve` must be a named list of character vectors or NULL" =
      is.null(resolve) ||
      (is.list(resolve) && length(resolve) > 0L &&
       !is.null(names(resolve)) && all(nzchar(names(resolve))) &&
       all(vapply(resolve, is.character, logical(1L))))
  )
  cols <- c(action, group, order)
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Column(s) not found in `data`: ",
         paste(utils::head(missing_cols, 5L), collapse = ", "), call. = FALSE)
  }
  if (identical(type, "anchored") && (is.null(anchor) || is.null(terminal))) {
    stop("type = \"anchored\" needs both `anchor` (what opens a pathway) and ",
         "`terminal` (what closes it).", call. = FALSE)
  }
  if (identical(type, "segments") && is.null(terminal)) {
    stop("type = \"segments\" needs `terminal`: it is what each segment ends on.",
         call. = FALSE)
  }

  d <- data.table::as.data.table(data[, cols, drop = FALSE])
  data.table::setnames(d, action, "..act")
  if (!is.null(order)) data.table::setorderv(d, c(group, order))
  d[, "..act" := as.character(get("..act"))]

  # Every cut below is a span between two row positions, which is only the
  # group's own events if the group's rows are contiguous. They need not be:
  # a learner can leave a step and come back to it later, so the log
  # interleaves. Rank each group by where it FIRST appears and sort on that
  # rank -- stable, so the event order inside a group is untouched and the
  # pathways still come back in the order they occur, but a group's scattered
  # runs are now adjacent. Without it a pathway silently swallows every event
  # that sat between its first and last row, including other groups' events.
  gkey <- do.call(paste, c(d[, group, with = FALSE], sep = "\r"))
  d <- d[order(match(gkey, gkey[!duplicated(gkey)]), method = "radix")]

  spans <- switch(type,
    unit     = .xp_whole(d, group, terminal),
    segments = .xp_segments(d, group, terminal),
    anchored = .xp_anchored(d, group, anchor, terminal)
  )
  if (nrow(spans) == 0L) {
    stop(errorCondition(
      "No pathway matched: no anchor/terminal state occurs in the data.",
      class = "nestimate_no_pathway", call = NULL))
  }
  .xp_assemble(d, spans, group, sep, resolve)
}

# ---- span finders -----------------------------------------------------------

# One span per group, optionally truncated at its last terminal state.
.xp_whole <- function(d, group, terminal) {
  d[, "..row" := .I]
  spans <- d[, list(a = min(get("..row")), z = max(get("..row"))), by = group]
  if (!is.null(terminal)) {
    ends <- d[get("..act") %in% terminal,
              list(z_term = max(get("..row"))), by = group]
    spans <- merge(spans, ends, by = group, all.x = TRUE)
    spans <- spans[!is.na(spans$z_term), ]
    spans$z <- spans$z_term
    spans$z_term <- NULL
  }
  spans
}

# Consecutive, non-overlapping spans: cut after every terminal state.
.xp_segments <- function(d, group, terminal) {
  d[, "..row" := .I]
  d[, "..seg" := {
    hit <- get("..act") %in% terminal
    cumsum(c(FALSE, utils::head(hit, -1L)))
  }, by = group]
  # keep only segments that actually reach a terminal state
  spans <- d[, list(a = min(get("..row")), z = max(get("..row")),
                    ok = any(get("..act") %in% terminal)),
             by = c(group, "..seg")]
  spans <- spans[spans$ok, ]
  spans$ok <- NULL
  spans[["..seg"]] <- NULL
  spans
}

# One span per anchor occurrence, closing at the next terminal state.
# findInterval() locates that state for every anchor at once.
.xp_anchored <- function(d, group, anchor, terminal) {
  d[, "..row" := .I]
  out <- d[, {
    act <- get("..act")
    a_pos <- which(act == anchor)
    t_pos <- which(act %in% terminal)
    if (length(a_pos) == 0L || length(t_pos) == 0L) {
      list(a = integer(0), z = integer(0))
    } else {
      # first terminal strictly after each anchor
      nxt <- t_pos[findInterval(a_pos, t_pos) + 1L]
      keep <- !is.na(nxt)
      base <- get("..row")[1L] - 1L
      list(a = base + a_pos[keep], z = base + nxt[keep])
    }
  }, by = group]
  out[!is.na(out$z), ]
}

# ---- assembly ---------------------------------------------------------------

# Expand each span to its rows once, then collapse to one path per pathway.
.xp_assemble <- function(d, spans, group, sep, resolve = NULL) {
  len <- spans$z - spans$a + 1L
  idx <- unlist(Map(seq.int, spans$a, spans$z), use.names = FALSE)
  pid <- rep.int(seq_len(nrow(spans)), len)

  flat <- data.table::data.table(..pid = pid, ..act = d[["..act"]][idx])
  paths <- flat[, list(
    path   = paste(get("..act"), collapse = sep),
    length = .N,
    opens  = get("..act")[1L],
    closes = get("..act")[.N]
  ), by = "..pid"]

  if (!is.null(resolve)) {
    # states present in each pathway, as one padded string per pathway
    present <- flat[, list(bag = paste0(sep, paste(unique(get("..act")),
                                                   collapse = sep), sep)),
                    by = "..pid"]
    lab <- rep("Unresolved", nrow(present))
    for (nm in rev(names(resolve))) {
      hit <- Reduce(`&`, lapply(resolve[[nm]], function(st)
        grepl(paste0(sep, st, sep), present$bag, fixed = TRUE)))
      lab[hit] <- nm
    }
    paths$ends   <- paths$closes
    paths$path   <- paste(paths$path, lab, sep = sep)
    paths$length <- paths$length + 1L
    paths$closes <- lab
  }

  keys <- spans[, group, with = FALSE]
  out <- cbind(as.data.frame(keys, stringsAsFactors = FALSE),
               as.data.frame(paths[, -1L], stringsAsFactors = FALSE))
  out$pathway <- stats::ave(seq_len(nrow(out)),
                            do.call(paste, c(keys, sep = "\r")),
                            FUN = seq_along)
  keep <- c(group, "pathway", "path", "length", "opens", "closes")
  if (!is.null(resolve)) keep <- c(keep, "ends")
  out[, keep]
}
