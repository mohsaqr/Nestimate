# ==============================================================================
# Object-level state colours.
#
# One palette, set once on the object, honoured by every figure drawn from it -
# Nestimate's own plot family and cograph's `splot()` alike.
#
# The cograph half works through the documented `meta$splot` producer contract:
# `defaults$node_fill` is merged into `splot()`'s arguments before its TNA /
# psychometric styling branch runs, so a stamped palette wins over the internal
# defaults without cograph needing to know about this function. Verified
# byte-identical to an explicit `splot(net, node_fill = ...)` on cograph 2.4.5.
#
# The Nestimate half is a fallback: every plot verb whose `state_colors` (or
# `colors`) argument is left NULL reads the stored palette, so an explicit
# argument still wins for one figure without disturbing the object.
# ==============================================================================


#' Set the state colours carried by a network object
#'
#' @description
#' `set_state_colors()` attaches a palette to a network object, so every
#' figure drawn from it uses the same colours: [sequence_plot()],
#' [distribution_plot()], [plot_state_frequencies()] and `cograph::splot()`.
#' A colour argument given to one of those plot verbs still wins for that
#' figure, and the object keeps its palette. `state_colors(x) <- value` is
#' the replacement form, and [state_colors()] reads back the colours in use.
#'
#' @param x A `netobject`, `netobject_group`, `mcml` or `htna`; any other
#'   object raises an error.
#' @param colors The palette: a named character vector of colours, such as
#'   `c(plan = "#0072B2", monitor = "#D55E00")`. Names are states, and for an
#'   `mcml` may also be cluster names. A name the object does not carry is
#'   reported in a message and draws nothing, so one project-wide palette can
#'   be attached to every object. States left unnamed are dealt the
#'   Okabe-Ito colours the palette does not already use, in the object's
#'   state order (recycled when there are more such states than colours).
#'   The palette replaces any palette set earlier, and `NULL` removes it. An
#'   unnamed vector raises an error.
#' @param value The palette, as for `colors`, in the replacement form
#'   `state_colors(x) <- value`.
#' @return `x` with the palette stored on it; the class is unchanged.
#'   [state_colors()] shows the colours it resolves to. For a `netobject` or
#'   `htna`, the palette is also written, in node order, into the
#'   `meta$splot` producer contract (`defaults$node_fill`), so
#'   `cograph::splot()` draws the same colours without being passed any. A
#'   `netobject_group` stores the palette on the group and on each member
#'   network, so a member drawn on its own carries it too.
#' @seealso \code{\link{state_colors}} to read the resolved palette back.
#' @examples
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' net <- set_state_colors(net, c(plan = "#0072B2", monitor = "#D55E00"))
#' state_colors(net)
#' @export
set_state_colors <- function(x, colors) UseMethod("set_state_colors")

#' @rdname set_state_colors
#' @export
set_state_colors.default <- function(x, colors) {
  stop(errorCondition(
    paste0(
      "`set_state_colors()` needs a netobject, netobject_group, mcml or ",
      "htna; got '",
      class(x)[1L],
      "'."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname set_state_colors
#' @export
set_state_colors.netobject <- function(x, colors) {
  x$state_colors <- .validate_state_colors(colors, .object_states(x))
  .stamp_node_fill(x)
}

#' @rdname set_state_colors
#' @export
set_state_colors.htna <- function(x, colors) set_state_colors.netobject(x, colors)

#' @rdname set_state_colors
#' @export
set_state_colors.mcml <- function(x, colors) {
  x$state_colors <- .validate_state_colors(colors, .object_states(x))
  x
}

#' @rdname set_state_colors
#' @export
set_state_colors.netobject_group <- function(x, colors) {
  # A group is a list of networks: the palette goes on the group (where the
  # plot verbs read it) and on each member (so a member drawn on its own, or
  # handed to splot(), carries it too).
  pal <- .validate_state_colors(colors, .object_states(x))
  out <- lapply(x, function(net) {
    if (!inherits(net, "netobject")) return(net)
    net$state_colors <- pal
    .stamp_node_fill(net)
  })
  attributes(out) <- attributes(x)
  # An attribute, not an element: a netobject_group must stay a list of
  # networks, or every verb that iterates over it meets a character vector.
  attr(out, "state_colors") <- pal
  out
}

#' @rdname set_state_colors
#' @export
`state_colors<-` <- function(x, value) set_state_colors(x, value)


#' The state colours an object will draw with
#'
#' @description
#' `state_colors()` returns the palette an object draws with: the colours set
#' with [set_state_colors()] plus the defaults filled in for every other
#' state, so the table shows what the figures use.
#'
#' @param x A `netobject`, `netobject_group`, `mcml` or `htna`; any other
#'   object raises an error.
#' @param ... Not used.
#' @return A data frame with one row per colour key the object carries
#'   (for a `netobject_group`, the states of all its networks): `state` (the
#'   key), `color` (the colour it draws with) and `source` (`"set"` when the
#'   palette named it, `"default"` when it fell back to Okabe-Ito). For an
#'   `mcml`, the cluster names follow the states.
#' @seealso \code{\link{set_state_colors}}.
#' @examples
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' state_colors(set_state_colors(net, c(plan = "#0072B2")))
#' @export
state_colors <- function(x, ...) UseMethod("state_colors")

#' @rdname state_colors
#' @export
state_colors.default <- function(x, ...) {
  stop(errorCondition(
    paste0(
      "`state_colors()` needs a netobject, netobject_group, mcml or htna; ",
      "got '",
      class(x)[1L],
      "'."),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname state_colors
#' @export
state_colors.netobject <- function(x, ...) .state_colors_table(x)

#' @rdname state_colors
#' @export
state_colors.htna <- function(x, ...) .state_colors_table(x)

#' @rdname state_colors
#' @export
state_colors.mcml <- function(x, ...) .state_colors_table(x)

#' @rdname state_colors
#' @export
state_colors.netobject_group <- function(x, ...) .state_colors_table(x)


# ---- internal ---------------------------------------------------------------

# The colour keys an object carries, in the order they should be shown: states
# first, then (for an mcml) the cluster names, which are keys of the
# multichannel Summary panel.
.object_states <- function(x) {
  if (inherits(x, "mcml")) {
    states <- unlist(lapply(x$clusters, function(z) z$labels), use.names = FALSE)
    return(unique(c(states, names(x$clusters))))
  }
  if (inherits(x, "netobject_group")) {
    members <- Filter(function(z) inherits(z, "netobject"), x)
    return(unique(unlist(lapply(members, .object_states), use.names = FALSE)))
  }
  unique(as.character(x$nodes$name))
}

# A stored palette is always a named character vector. Names the object does
# not carry are kept (one project-wide palette can sit on every object) and
# are reported here, once, at set time.
.validate_state_colors <- function(colors, keys) {
  if (is.null(colors)) {
    return(NULL)
  }
  .assert(
    "`colors` must be a character vector of colours" = is.character(colors),
    "`colors` must be named: c(state = \"#0072B2\", ...)" =
      .is_named_colors(colors),
    "every colour in `colors` must be named" = all(nzchar(names(colors))))
  .report_unmatched_colors(colors, keys)
  colors
}

# The set-time counterpart of .report_unused_colors(): nothing is drawn yet
# and nothing is dropped, so the message says what holds - the names match no
# state of the object and stay in its palette.
.report_unmatched_colors <- function(colors, keys) {
  nm        <- names(colors)
  unmatched <- setdiff(nm, keys)
  if (length(unmatched) == 0L) {
    return(invisible(NULL))
  }
  shown <- paste0(paste(utils::head(unmatched, 6L), collapse = ", "),
                  if (length(unmatched) > 6L) ", ..." else "")
  if (length(unmatched) == length(nm)) {
    message("`colors`: no name matches a state of this object, so every ",
            "state keeps its default colour. Unmatched: ", shown, ".")
  } else {
    message("`colors`: ", length(unmatched), " of ", length(nm),
            " names are not drawn, because the object has no such state: ",
            shown, ". They stay in its palette.")
  }
  invisible(NULL)
}

# Mirror the palette into the cograph producer contract, in node order, so
# `cograph::splot(x)` draws the same colours without being passed anything.
# Nodes the palette does not name keep cograph's own fill, which is what
# leaving a NULL entry in `node_fill` would do - so they are filled from the
# same Okabe-Ito default the Nestimate plots use, and the two agree.
.stamp_node_fill <- function(x) {
  keys <- .object_states(x)
  if (is.null(x$state_colors)) {
    x$meta$splot$defaults$node_fill <- NULL
    if (length(x$meta$splot$defaults) == 0L) x$meta$splot$defaults <- NULL
    if (length(x$meta$splot) == 0L) x$meta$splot <- NULL
    return(x)
  }
  x$meta$splot$defaults$node_fill <- unname(.resolve_state_colors(x)[keys])
  x
}

# The full palette an object draws with: the stored colours over an Okabe-Ito
# base, named by key. This is the single resolution point - the table, the
# cograph stamp and the plot fallback all read it.
# The palette set with set_state_colors(): an attribute on a
# netobject_group, a field on a single object.
.palette_of <- function(x) {
  # `[[`, not `$`: a tibble warns "Unknown or uninitialised column" when
  # `$` probes a column it does not have; `[[` returns NULL silently.
  if (inherits(x, "netobject_group")) attr(x, "state_colors") else
    x[["state_colors"]]
}

.resolve_state_colors <- function(x) {
  keys <- .object_states(x)
  stats::setNames(.fill_state_colors(keys, .palette_of(x) %||% character(0)),
                  keys)
}

.state_colors_table <- function(x) {
  keys <- .object_states(x)
  pal  <- .resolve_state_colors(x)
  data.frame(state  = keys,
             color  = unname(pal[keys]),
             source = ifelse(keys %in% names(.palette_of(x) %||% character(0)),
                             "set", "default"),
             stringsAsFactors = FALSE)
}

# The palette a plot verb falls back to when its `state_colors` argument is
# NULL. The FULLY RESOLVED palette is returned - the keys the user set plus the
# defaults filled in - not just the keys they set. Returning only the set keys
# leaves each figure to deal the remaining colours itself, in whatever order it
# sorts its states: `plot_state_frequencies()` sorts by frequency and
# `sequence_plot()` alphabetically, so the same state came out amber in one and
# black in the other. Resolving once here, in the object's canonical key order,
# makes every downstream verb a pure name lookup, and sort order stops
# mattering. The `state_colors(x)` table is the same resolution, so the table
# and the figures cannot disagree.
.stored_state_colors <- function(x) {
  if (!is.list(x) || !length(.palette_of(x))) return(NULL)
  .resolve_state_colors(x)
}
