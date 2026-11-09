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
#' @description
#' `extract_pathways()` cuts a long event log into pathways, runs of
#' consecutive events within one unit (an actor, or an actor's session), and
#' returns them as a data frame with one row per pathway, ready to tabulate
#' against an outcome or to pass to a sequence model. Three cuts are
#' available through `type`: each unit as one pathway, consecutive segments
#' that each end on a terminal state, or one pathway from every occurrence of
#' an `anchor` state to the next terminal state.
#'
#' With `terminal`, the cuts keep only what reaches a terminal state: a
#' unit without one, events after a unit's last terminal state, and anchors
#' with no later terminal state are left out without a message, so count
#' the units and events kept before comparing pathway frequencies across
#' groups.
#'
#' @param data A long data frame with one row per event, in time order
#'   within each unit unless `order` is given. A unit's rows need not be
#'   contiguous: rows of other units in between are not part of its
#'   pathways.
#' @param action The name of the column holding each event's state.
#' @param group The names of the columns defining the unit a pathway may not
#'   cross, typically the actor, or the actor and the session.
#' @param order The name of a column giving the order of events within a
#'   unit; the rows are sorted by `group` and then by it. Default `NULL`,
#'   which keeps the row order.
#' @param type The cut to make:
#'   \describe{
#'     \item{`"unit"`}{(default) one pathway per unit. With `terminal`, each
#'       unit is cut after its last terminal state, so the pathway ends on an
#'       outcome rather than on whatever followed; units without a terminal
#'       state are left out.}
#'     \item{`"segments"`}{consecutive, non-overlapping pathways: each unit
#'       is cut after every terminal state, so each pathway ends on one.
#'       Events after a unit's last terminal state belong to no pathway.}
#'     \item{`"anchored"`}{one pathway per occurrence of `anchor`, closing at
#'       the first terminal state after it. Pathways may overlap, and events
#'       between a terminal state and the next anchor belong to none.}
#'   }
#' @param anchor The single state that opens a pathway when `type =
#'   "anchored"`, where it is required; the other cuts ignore it. It may
#'   also be one of the `terminal` states: a pathway then closes at the next
#'   terminal state after the anchor, not at the anchor itself. Default
#'   `NULL`.
#' @param terminal The states that close a pathway, the closing state
#'   included. Required for `"segments"` and `"anchored"`; optional for
#'   `"unit"`, where it cuts each unit after its last terminal state.
#'   Default `NULL`.
#' @param resolve A named list that appends a resolution label as each
#'   pathway's final state. Each element is a character vector of states; a
#'   pathway takes the name of the first element whose states all occur in
#'   it, so order the list from most to least specific, and pathways
#'   matching none are labelled `"Unresolved"`. The label becomes the last
#'   state of `path`, counts in `length`, and is the value of `closes`; the
#'   observed last state is kept in `ends`. Default `NULL` (no label).
#' @param sep The separator placed between states in `path`. Default
#'   `" -> "`.
#'
#' @return A data frame with one row per pathway, in the order the pathways
#'   occur within each unit. Units appear in the order they first occur in
#'   `data`, or sorted by `group` when `order` is given and when `type =
#'   "unit"` is combined with `terminal`. Columns: the `group` columns,
#'   `pathway` (the pathway's number within its unit, from 1), `path` (the
#'   states joined with `sep`), `length` (the number of states), `opens`
#'   (the first state) and `closes` (the last state). With `resolve`,
#'   `closes` holds the resolution label and an `ends` column the observed
#'   last state. When no pathway is found
#'   (for example, the anchor or the terminal states never occur), an error
#'   of class `nestimate_no_pathway` is raised.
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
  type <- .ho_match_arg(type, "type")
  .assert(
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
    stop(errorCondition(paste0("Column(s) not found in `data`: ",
         paste(utils::head(missing_cols, 5L), collapse = ", ")),
                        class = "nestimate_bad_input", call = NULL))
  }
  if (identical(type, "anchored") && (is.null(anchor) || is.null(terminal))) {
    stop(errorCondition(
      paste0(
        "type = \"anchored\" needs both `anchor` (what opens a pathway) and ",
        "`terminal` (what closes it)."),
      class = "nestimate_bad_input", call = NULL))
  }
  if (identical(type, "segments") && is.null(terminal)) {
    stop(errorCondition(
      "type = \"segments\" needs `terminal`: it is what each segment ends on.",
      class = "nestimate_bad_input", call = NULL))
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
    # A loop on purpose: the labels are written in reverse so that, where a
    # pathway satisfies several rules, the first rule listed wins.
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
