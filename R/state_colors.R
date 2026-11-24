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


#' Set the state colours of a network object
#'
#' @description
#' `set_state_colors()` attaches a palette to a network object, so every
#' figure drawn from it uses the same colours. `state_colors(x) <- value` is
#' the replacement form, and [state_colors()] reads back the colours in use.
#'
#' One rule holds for colours across the package:
#' * With no colours set, the states are drawn in the Okabe-Ito colours,
#'   dealt in the alphabetical order of the states, so a state has the same
#'   colour in every figure whatever the order of the network's nodes.
#' * A palette set here replaces those colours for the states it names.
#' * Objects computed from the network carry its colours: groups and their
#'   members, [build_clusters()], [build_mmm()] and the networks
#'   [build_network()] and [cluster_network()] build from them,
#'   [bootstrap_network()], [certainty()], [permutation()],
#'   [bayes_compare()], [net_edge_betweenness()], [entropy_network()],
#'   [net_prune()], [stability_landscape()], [landscape_trend()] and
#'   [landscape_surface()]. A cluster or group that lacks a state keeps the
#'   colours of the others.
#' * Every figure that colours states follows them: [sequence_plot()],
#'   [distribution_plot()], [plot_state_frequencies()], the `plot()`
#'   methods of [stability_landscape()] and [landscape_trend()], and
#'   `cograph::splot()` (through the network's `meta$splot` contract).
#' * `state_colors =` on one of those figures changes that figure only.
#'
#' @param x A `netobject`, `netobject_group`, `mcml` or `htna`; any other
#'   object raises an error.
#' @param colors The palette: a named character vector of colours, such as
#'   `c(plan = "#0072B2", monitor = "#D55E00")`. Names are states (for a
#'   network built with `labels`, the labels it draws), and for an
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
  # a tna model, a group_tna or a bare cograph_network: convert, dispatch again
  if (.is_foreign_network(x)) {
    return(set_state_colors(.native_network(x, "set_state_colors()"), colors))
  }
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
  # An attribute, not an element: a netobject_group must stay a list of
  # networks, or every verb that iterates over it meets a character vector.
  attr(x, "state_colors") <- pal
  .stamp_group_palette(x)
}

#' @rdname set_state_colors
#' @export
`state_colors<-` <- function(x, value) set_state_colors(x, value)


#' Get the state colours of an object
#'
#' @description
#' `state_colors()` returns the palette an object draws with: the colours set
#' with [set_state_colors()] plus the defaults filled in for every other
#' state, so the table shows what the figures use.
#'
#' @param x A `netobject`, `netobject_group`, `mcml` or `htna`. A `tna` model
#'   or a `group_tna` (package tna) is first rebuilt as a network from its own
#'   sequences (see [as_netobject()]). Any other object raises an error.
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
  # a tna model, a group_tna or a bare cograph_network: convert, dispatch again
  if (.is_foreign_network(x)) {
    return(state_colors(.native_network(x, "state_colors()"), ...))
  }
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
  # an object without a node table (a clustering, a mixture model) carries
  # the resolved palette of the network it was derived from: its keys
  nodes <- x[["nodes"]]
  if (is.null(nodes)) {
    return(union(names(.palette_of(x)), names(.base_of(x))) %||% character(0))
  }
  # The names the figures draw, in cograph's order (`labels`, then `label`):
  # a network built with `labels` draws (and holds in $weights and $data) its
  # labels, so a palette keyed by the original names would colour states that
  # are never drawn. `[[` because `$` partially matches `label` to `labels`.
  unique(as.character(nodes[["labels"]] %||% nodes[["label"]] %||%
                        nodes[["name"]]))
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

# Carry the palette of `from` onto an object derived from it (a clustering,
# a mixture model, the networks built from them, an edge-betweenness or
# entropy network), so every figure of the derived object draws each state
# in the colour the source drew it. The fully resolved palette is carried -
# the colours set plus the defaults dealt to the other states - so a state
# the user did not name keeps its colour too. Nothing happens when `from`
# has no palette.
.carry_state_colors <- function(to, from) {
  if (!is.list(from) || identical(to, from)) return(to)
  explicit <- .palette_of(from)
  # an object that is not a network resolves the same Okabe-Ito defaults by
  # the same rule, so it carries a palette only when one was set upstream
  is_network <- inherits(to, c("netobject", "netobject_group"))
  if (!is_network && !length(explicit) && !length(.base_of(from))) {
    return(to)
  }
  base <- .stored_state_colors(from)
  if (!length(base)) return(to)
  if (inherits(to, "netobject_group")) {
    attr(to, "state_colors") <- explicit
    attr(to, "palette_base") <- base
    return(.stamp_group_palette(to))
  }
  to$state_colors <- explicit
  to$palette_base <- base
  if (inherits(to, "netobject")) to <- .stamp_node_fill(to)
  to
}

# The palette a network gets when build_network() returns it: a network is
# stamped for cograph with its resolved colours (Okabe-Ito when none are set),
# a group gets one palette resolved over all its states, which every member
# takes as its base, so a member lacking a state keeps the others' colours.
.finalize_palette <- function(x) {
  if (inherits(x, "netobject_group")) return(.stamp_group_palette(x))
  if (inherits(x, "netobject")) return(.stamp_node_fill(x))
  if (inherits(x, "netobject_ml")) {
    members <- vapply(x, inherits, logical(1L), "netobject")
    x[members] <- lapply(x[members], .stamp_node_fill)
  }
  x
}

.stamp_group_palette <- function(x) {
  # members take the group's resolved palette as their base, so a member
  # lacking a state keeps the others' colours; the group itself keeps only
  # the base it inherited, so a later palette on it starts from that, not
  # from the colours an earlier palette set
  resolved <- .resolve_state_colors(x)
  explicit <- attr(x, "state_colors", exact = TRUE)
  members <- vapply(x, inherits, logical(1L), "netobject")
  out <- x
  out[members] <- lapply(x[members], function(net) {
    net$state_colors <- explicit
    # the group's colours for this member's own states, in its node order
    net$palette_base <- resolved[intersect(.object_states(net), names(resolved))]
    .stamp_node_fill(net)
  })
  attributes(out) <- attributes(x)
  out
}

# A result drawn by cograph (a bootstrap, a certainty, a permutation or a
# Bayesian comparison) declares the node colours of the network it came from
# in its own meta$splot contract: cograph applies a result's contract before
# its renderer fills its default palette, so the figure draws the source's
# colours. The pruned model network inside a bootstrap carries them too.
.declare_source_fill <- function(result, source) {
  if (!inherits(source, "netobject")) return(result)
  # a network built before networks were stamped is resolved here
  fill <- source[["meta"]][["splot"]][["defaults"]][["node_fill"]] %||%
    .stamp_node_fill(source)[["meta"]][["splot"]][["defaults"]][["node_fill"]]
  if (is.null(fill)) return(result)
  result$meta$splot$defaults$node_fill <- fill
  result$meta$splot$defaults$label_color <- .label_contrast(fill)
  # the pruned model network carries the palette itself, so state_colors()
  # of it and its figures agree with the source
  if (inherits(result[["model"]], "netobject")) {
    result$model <- .carry_state_colors(result$model, source)
  }
  result
}

# Mirror the palette into the cograph producer contract, in node order, so
# `cograph::splot(x)` draws the same colours without being passed anything.
# Nodes the palette does not name keep cograph's own fill, which is what
# leaving a NULL entry in `node_fill` would do - so they are filled from the
# same Okabe-Ito default the Nestimate plots use, and the two agree.
.stamp_node_fill <- function(x) {
  keys <- .object_states(x)
  if (!length(keys)) return(x)
  # one colour per node, in node order, even when two nodes share a label
  nodes <- x[["nodes"]]
  per_node <- as.character(nodes[["labels"]] %||% nodes[["label"]] %||%
                             nodes[["name"]] %||% keys)
  fill <- unname(.resolve_state_colors(x)[per_node])
  x$meta$splot$defaults$node_fill <- fill
  x$meta$splot$defaults$label_color <- .label_contrast(fill)
  x
}

# Label colour for text drawn on a node: black, or white where black would
# fall below a 3:1 contrast ratio (WCAG 2 relative luminance), which among
# the Okabe-Ito colours is only black itself.
.label_contrast <- function(fill) {
  rgb <- grDevices::col2rgb(fill) / 255
  lin <- ifelse(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055)^2.4)
  lum <- colSums(lin * c(0.2126, 0.7152, 0.0722))
  ifelse((lum + 0.05) / 0.05 < 3, "white", "black")
}

# The full palette an object draws with: the stored colours over an Okabe-Ito
# base, named by key. This is the single resolution point - the table, the
# cograph stamp and the plot fallback all read it.
# The palette set with set_state_colors(): an attribute on a
# netobject_group, a field on a single object.
.palette_of <- function(x) {
  # `[[`, not `$`: a tibble warns "Unknown or uninitialised column" when
  # `$` probes a column it does not have; `[[` returns NULL silently.
  if (inherits(x, "netobject_group")) attr(x, "state_colors", exact = TRUE) else
    x[["state_colors"]]
}

# The base palette carried from the network an object was derived from.
.base_of <- function(x) {
  if (inherits(x, "netobject_group")) attr(x, "palette_base", exact = TRUE) else
    x[["palette_base"]]
}

# The colours set, then the carried base, then Okabe-Ito for any other key,
# dealt in the sorted order of the keys, so the default colour of a state
# does not depend on the order of a network's nodes.
.resolve_state_colors <- function(x) {
  keys <- as.character(.object_states(x))
  explicit <- .palette_of(x) %||% character(0)
  explicit <- explicit[names(explicit) %in% keys]
  base <- .base_of(x) %||% character(0)
  base <- base[setdiff(names(base), names(explicit))]
  # a carried default that a set colour now uses is dealt again, so two
  # states never share a colour the user did not give them both
  base <- base[!(unname(base) %in% unname(explicit))]
  dealt <- .order_keys(keys)
  if (inherits(x, "mcml")) {
    # states first, in the order a network of them deals them, then the
    # cluster names from the colours the states leave
    states <- unique(as.character(unlist(lapply(x$clusters, function(z)
      z$labels), use.names = FALSE)))
    dealt <- c(.order_keys(states), .order_keys(setdiff(keys, states)))
  }
  stats::setNames(.fill_state_colors(dealt, c(explicit, base)), dealt)[keys]
}

# The palette a figure of `x` draws with when the call gives `override`:
# a named override is added to the colours set on `x` and resolved by the
# same rule (the other states keep their colours, re-dealt only where the
# override took one); an unnamed override is used as given.
.with_override <- function(x, override = NULL) {
  stored <- .stored_state_colors(x)
  if (is.null(override)) return(stored)
  if (is.null(stored) || !.is_named_colors(override)) return(override)
  explicit <- .palette_of(x) %||% character(0)
  explicit[names(override)] <- unname(override)
  if (inherits(x, "netobject_group")) {
    attr(x, "state_colors") <- explicit
  } else {
    x$state_colors <- explicit
  }
  merged <- c(.resolve_state_colors(x),
              override[setdiff(names(override), .object_states(x))])
  # the figure reports the names the call gave, not the object's own keys
  attr(merged, "user_names") <- names(override)
  merged
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
  if (!is.list(x)) return(NULL)
  # a network or a group always draws with its resolved palette (Okabe-Ito
  # when none is set), so every figure of it agrees; any other object only
  # when it carries a palette
  is_network <- inherits(x, c("netobject", "netobject_group", "mcml"))
  if (!is_network && !length(.palette_of(x)) && !length(.base_of(x))) {
    return(NULL)
  }
  .resolve_state_colors(x)
}
