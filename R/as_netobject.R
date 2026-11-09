# Boundary between psychnet (psychometric math, emits a lean cograph_network)
# and Nestimate (owns the canonical netobject schema). psychnet is imported by
# Nestimate, never the reverse, so the converter lives here - on the side that
# owns the target schema and can see both types.

#' Coerce a network object to a Nestimate netobject
#'
#' @description
#' `as_netobject()` converts a \pkg{psychnets} result (class
#' `c("psychnet", "cograph_network")`) or any bare `cograph_network` into the
#' `c("netobject", "cograph_network")` object the rest of Nestimate reads, so
#' it can be passed to [net_centrality()] and the other verbs that read a
#' network. The resampling verbs ([bootstrap_network()],
#' [network_reliability()]) need the source data, which a psychnets result
#' does not carry, so build the network with [build_network()] for those. A
#' `netobject` is returned unchanged.
#'
#' @details
#' For a `psychnet`, the edge table is rebuilt with integer node indices (a
#' `psychnet` stores character-labelled edges), the estimator name is kept
#' in `method`, and every field specific to the `psychnet`, including the
#' graphical-lasso `kkt` optimality certificate, is kept in the `psychnet`
#' entry of `meta`, so nothing is lost in the conversion.
#'
#' For a bare `cograph_network`, the method is read from the `tna` entry of
#' its `meta`; without one, a symmetric weight matrix is taken as
#' `"co_occurrence"` and any other as `"relative"`. A missing `directed`
#' field is read as `TRUE`. Numeric state codes in the stored data of a
#' converted tna model are decoded to their labels.
#'
#' @param x The network to convert: a `psychnet` object, a `cograph_network`
#'   or a `netobject`. Any other object raises an error.
#' @return A `netobject` (class `c("netobject", "cograph_network")`), read by
#'   the same verbs as one from [build_network()]; [validate_netobject()]
#'   checks its structure.
#' @seealso \code{\link{validate_netobject}}
#' @examples
#' net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
#' identical(as_netobject(net), net) # netobjects pass through unchanged
#' @export
as_netobject <- function(x) {
  UseMethod("as_netobject")
}

#' @rdname as_netobject
#' @export
as_netobject.netobject <- function(x) {
  x
}

#' @rdname as_netobject
#' @export
as_netobject.psychnet <- function(x) {
  net <- .as_netobject(x)
  # .as_netobject() infers method from matrix symmetry; psychnet carries the
  # true estimator name, so prefer it.
  method <- x$method %||% net$method
  net$method <- method
  net$meta$tna$method <- method
  net$meta$source <- "psychnet"
  # Preserve every method-specific field psychnet attached (precision, lambda,
  # cor_matrix, ebic_path, kkt, ...) so the certificate survives the conversion.
  core <- c("weights", "nodes", "edges", "directed", "method", "n")
  net$meta$psychnet <- x[setdiff(names(x), core)]
  net
}

#' @rdname as_netobject
#' @export
as_netobject.cograph_network <- function(x) {
  .as_netobject(x)
}

#' @rdname as_netobject
#' @export
as_netobject.default <- function(x) {
  stop(errorCondition(paste0(
    "as_netobject() needs a psychnet, cograph_network, or netobject; got <",
    paste(class(x), collapse = "/"), ">."),
    class = "nestimate_bad_input", call = NULL))
}

#' Validate a netobject / cograph_network against the shared schema
#'
#' @description
#' `validate_netobject()` checks that an object meets the structure Nestimate
#' `netobject`s and \pkg{psychnets} objects share, so either can be passed
#' across the package boundary. It is the one place that says what a network
#' object is, so a change on either side (a renamed field, a mistyped edge
#' column) fails here instead of further down, for instance when the network
#' is drawn.
#'
#' The checks are the shared part only: the object inherits
#' `cograph_network`; `weights` is a square numeric matrix with node names as
#' dimnames; `nodes` is a data frame with columns `id`, `label` and `name`,
#' one row per node; `edges` is a data frame with columns `from`, `to` and a
#' numeric `weight`, whose endpoints are integer node indices (Nestimate) or
#' node labels (\pkg{psychnets}); `directed` is a single logical value; and
#' `method` is a single string. The `x`/`y` layout columns of `nodes` and
#' Nestimate's own fields (`data`, `level`, ...) are not required.
#'
#' @param x The object to check, expected to be a `cograph_network`.
#' @return `TRUE`, invisibly, when `x` passes every check. Otherwise an error
#'   that lists every check it failed.
#' @seealso \code{\link{as_netobject}}
#' @examples
#' net <- build_cor(data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50)))
#' validate_netobject(net)
#' @export
validate_netobject <- function(x) {
  problems <- character(0)
  note <- function(msg) problems <<- c(problems, msg)

  if (!inherits(x, "cograph_network")) {
    note("object does not inherit class 'cograph_network'")
  }

  w <- x$weights
  square <- FALSE
  if (is.null(w) || !is.matrix(w) || !is.numeric(w)) {
    note("$weights must be a numeric matrix")
  } else {
    square <- nrow(w) == ncol(w)
    if (!square) note("$weights must be a square matrix")
    if (is.null(dimnames(w))) note("$weights must carry node names as dimnames")
  }

  nodes <- x$nodes
  have_nodes <- is.data.frame(nodes) && all(c("id", "label", "name") %in% names(nodes))
  if (!have_nodes) {
    note("$nodes must be a data.frame with columns id, label, name")
  } else if (square && nrow(nodes) != ncol(w)) {
    note("$nodes row count must equal the network dimension")
  }

  edges <- x$edges
  if (!is.data.frame(edges) || !all(c("from", "to", "weight") %in% names(edges))) {
    note("$edges must be a data.frame with columns from, to, weight")
  } else {
    if (!is.numeric(edges$weight)) note("$edges$weight must be numeric")
    if (have_nodes && nrow(edges) > 0L) {
      refs_ok <- function(col) {
        if (is.numeric(col)) all(col >= 1L & col <= nrow(nodes))
        else all(col %in% nodes$label)
      }
      if (!refs_ok(edges$from)) note("$edges$from does not reference valid nodes")
      if (!refs_ok(edges$to)) note("$edges$to does not reference valid nodes")
    }
  }

  if (!is.logical(x$directed) || length(x$directed) != 1L) {
    note("$directed must be a single logical value")
  }
  if (!is.character(x$method) || length(x$method) != 1L) {
    note("$method must be a single character value")
  }

  if (length(problems) > 0L) {
    stop(errorCondition(paste0(
      "invalid netobject / cograph_network:\n",
      paste0("  - ", problems, collapse = "\n")),
      class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}
