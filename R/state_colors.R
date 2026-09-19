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
#' Attaches a palette to the object so every figure drawn from it uses the same
#' colours: \code{\link{sequence_plot}}, \code{\link{distribution_plot}},
#' \code{\link{plot_state_frequencies}} and \code{cograph::splot()}.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml} or
#'   \code{htna}.
#' @param colors A named character vector of colours, e.g.
#'   \code{c(plan = "#0072B2", monitor = "#D55E00")}. Names are states, and for
#'   an \code{mcml} may also be cluster names. Names the object does not carry
#'   are dropped with a message, so one project-wide palette can be attached to
#'   every object. States you do not name keep the default Okabe-Ito colour.
#'   \code{NULL} removes a palette set earlier.
#' @return \code{x}, with the palette stored in \code{x$state_colors} and, for
#'   an object carrying \code{$nodes}, mirrored into
#'   \code{x$meta$splot$defaults$node_fill} in node order so
#'   \code{cograph::splot()} honours it. The class is unchanged.
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
  stop("`set_state_colors()` needs a netobject, netobject_group, mcml or ",
       "htna; got '", class(x)[1L], "'.", call. = FALSE)
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
  out$state_colors <- pal
  out
}

#' @rdname set_state_colors
#' @export
`state_colors<-` <- function(x, value) set_state_colors(x, value)


#' The state colours an object will draw with
#'
#' Reads back the palette an object resolves to: the colours set with
#' \code{\link{set_state_colors}} plus the defaults filled in for everything
#' else, so the table is what the figures actually use.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml} or
#'   \code{htna}.
#' @param ... Ignored, for method consistency.
#' @return A \code{data.frame}, one row per colour key the object carries, with
#'   columns \code{state} (the key), \code{color} (the hex colour it draws
#'   with) and \code{source} (\code{"set"} when the palette named it,
#'   \code{"default"} when it fell back to Okabe-Ito). For an \code{mcml} the
#'   cluster names appear after the states.
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
  stop("`state_colors()` needs a netobject, netobject_group, mcml or htna; ",
       "got '", class(x)[1L], "'.", call. = FALSE)
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

# A stored palette is always a named character vector restricted to nothing -
# extra names are dropped here, once, at set time, rather than at every plot.
.validate_state_colors <- function(colors, keys) {
  if (is.null(colors)) {
    return(NULL)
  }
  stopifnot(
    "`colors` must be a character vector of colours" = is.character(colors),
    "`colors` must be named: c(state = \"#0072B2\", ...)" =
      .is_named_colors(colors),
    "every colour in `colors` must be named" = all(nzchar(names(colors))))
  .report_unused_colors(colors, keys, arg = "colors")
  colors
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
.resolve_state_colors <- function(x) {
  keys <- .object_states(x)
  stats::setNames(.fill_state_colors(keys, x$state_colors %||% character(0)),
                  keys)
}

.state_colors_table <- function(x) {
  keys <- .object_states(x)
  pal  <- .resolve_state_colors(x)
  data.frame(state  = keys,
             color  = unname(pal[keys]),
             source = ifelse(keys %in% names(x$state_colors %||% character(0)),
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
  if (!is.list(x) || is.null(x$state_colors) || !length(x$state_colors)) {
    return(NULL)
  }
  .resolve_state_colors(x)
}
