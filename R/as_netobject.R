# Boundary between psychnet (psychometric math, emits a lean cograph_network)
# and Nestimate (owns the canonical netobject schema). psychnet is imported by
# Nestimate, never the reverse, so the converter lives here — on the side that
# owns the target schema and can see both types.

#' Coerce a network object to a Nestimate netobject
#'
#' Promotes a \pkg{psychnet} result (class \code{c("psychnet",
#' "cograph_network")}) or any bare \code{cograph_network} to the dual-class
#' \code{c("netobject", "cograph_network")} used throughout Nestimate, so it
#' dispatches to the package's verbs (\code{centrality()}, \code{plot()},
#' bootstrap, reliability, ...). A \code{netobject} is returned unchanged.
#'
#' The psychnet method re-derives the integer-indexed edge table that Nestimate
#' expects (psychnet stores character-labelled edges), preserves the estimator
#' name in \code{$method}, and parks every psychnet-specific field — including
#' the graphical-lasso \code{$kkt} optimality certificate — under
#' \code{$meta$psychnet} so nothing is lost in translation.
#'
#' @param x A \code{psychnet} object, a \code{cograph_network}, or a
#'   \code{netobject}.
#' @return A \code{c("netobject", "cograph_network")} object.
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
  stop("as_netobject() needs a psychnet, cograph_network, or netobject; got <",
       paste(class(x), collapse = "/"), ">.", call. = FALSE)
}

#' Validate a netobject / cograph_network against the shared schema
#'
#' Enforces the structural contract that both Nestimate netobjects and
#' \pkg{psychnet} objects must satisfy to be interchangeable across the package
#' boundary. This is the single place that says what "a network object" means,
#' so a drift on either side (a renamed field, a mistyped edge column) fails
#' loudly here rather than mis-rendering three layers downstream.
#'
#' The contract is deliberately the \emph{shared} subset: the \code{$nodes}
#' \code{x}/\code{y} layout columns and the Nestimate pipeline fields
#' (\code{$data}, \code{$level}, ...) are not required, and \code{$edges}
#' endpoints may be either integer node indices (Nestimate) or character labels
#' (psychnet).
#'
#' @param x An object expected to satisfy the \code{cograph_network} contract.
#' @return Invisibly \code{TRUE} if \code{x} conforms; otherwise stops with the
#'   full list of violations.
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
    stop("invalid netobject / cograph_network:\n",
         paste0("  - ", problems, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}
