# ---- Network-level descriptive summary ----
#
# Same metric set as `tna::summary.tna()` but computed from `$weights`
# directly so we do not depend on igraph at runtime. Distances reuse
# `.floyd_warshall_sp()` from R/centrality_measures.R.

#' @rdname build_network
#' @export
summary.netobject <- function(object, ...) {
  .assert("`object` must be a netobject or cograph_network" =
            inherits(object, "netobject") || inherits(object, "cograph_network"))
  W <- object$weights
  .assert("`object$weights` must be a square matrix" =
            is.matrix(W) && nrow(W) == ncol(W))
  directed <- isTRUE(object$directed)

  .summary_metrics_from_weights(W, directed)
}

#' @rdname build_network
#' @export
summary.netobject_group <- function(object, combined = TRUE, ...) {
  .assert("`object` must be a netobject_group" =
            inherits(object, "netobject_group"),
          "`combined` must be TRUE or FALSE" =
            is.logical(combined) && length(combined) == 1L && !is.na(combined))

  per_group <- lapply(object, summary)

  if (!combined) {
    return(structure(per_group, class = c("summary.netobject_group", "list"),
                     combined = FALSE))
  }

  metric <- per_group[[1L]]$metric
  values <- vapply(per_group, function(s) s$value,
                   numeric(length(metric)))
  out <- data.frame(metric = metric, values, check.names = FALSE,
                    stringsAsFactors = FALSE)
  names(out)[-1L] <- names(per_group)
  structure(out, class = c("summary.netobject_group", "data.frame"),
            combined = TRUE)
}


# ---- Internal: compute the 13 metrics from a weight matrix ----

.summary_metrics_from_weights <- function(W, directed) {
  # Diagonal is preserved as estimated: self-loops count as edges, as in
  # igraph's `edge_density` and tna's `summary.tna()`. An edge is any non-zero
  # weight, so signed networks keep their negative edges, and an undirected
  # edge is counted once (upper triangle), in line with `$n_edges`.
  n <- nrow(W)
  present <- W != 0
  edge_count <- if (directed) {
    sum(present)
  } else {
    sum(present[upper.tri(present, diag = TRUE)])
  }

  density <- if (n > 1L) {
    if (directed) edge_count / (n * (n - 1))
    else edge_count / (n * (n - 1) / 2)
  } else 0
  density <- min(density, 1)

  # Shortest paths need non-negative costs: a signed edge costs its magnitude.
  D <- .floyd_warshall_sp(abs(W), invert = FALSE)$D
  finite_off <- is.finite(D) & D > 0
  mean_distance <- if (any(finite_off)) mean(D[finite_off]) else NaN

  abs_W <- abs(W)
  row_strength <- rowSums(abs_W)
  col_strength <- colSums(abs_W)
  out_degree   <- rowSums(W != 0)
  in_degree    <- colSums(W != 0)

  # Centralization uses degrees without self-loops (matches
  # `igraph::centr_degree(..., loops = FALSE)`).
  W_no_loops <- W
  diag(W_no_loops) <- 0
  out_degree_nl <- rowSums(W_no_loops != 0)
  in_degree_nl  <- colSums(W_no_loops != 0)

  # Rows are sources: out-strength is the row sum, in-strength the column
  # sum (their means are equal by construction; the SDs are not). Until
  # 1.0.5 the two were swapped, as in tna before 1.5.4.
  out_strength <- row_strength
  in_strength  <- col_strength

  cent_out <- .freeman_degree_centralization(out_degree_nl, n, directed)
  cent_in  <- .freeman_degree_centralization(in_degree_nl,  n, directed)

  # igraph's default `reciprocity()` excludes self-loops, so we measure
  # mutual-edge fraction over off-diagonal positions only.
  present_off <- present
  diag(present_off) <- FALSE
  reciprocity <- if (directed && any(present_off)) {
    mutual <- present_off & t(present_off)
    sum(mutual) / sum(present_off)
  } else if (!directed) {
    1
  } else {
    NaN
  }

  metric <- c(
    "Node Count",
    "Edge Count",
    "Network Density",
    "Mean Distance",
    "Mean Out-Strength",
    "SD Out-Strength",
    "Mean In-Strength",
    "SD In-Strength",
    "Mean Out-Degree",
    "SD Out-Degree",
    "Centralization (Out-Degree)",
    "Centralization (In-Degree)",
    "Reciprocity"
  )
  value <- c(
    n,
    edge_count,
    density,
    mean_distance,
    mean(out_strength, na.rm = TRUE),
    stats::sd(out_strength, na.rm = TRUE),
    mean(in_strength, na.rm = TRUE),
    stats::sd(in_strength, na.rm = TRUE),
    mean(out_degree),
    stats::sd(out_degree),
    cent_out,
    cent_in,
    reciprocity
  )
  out <- data.frame(metric = metric, value = unname(value),
                    stringsAsFactors = FALSE)
  structure(out, class = c("summary.netobject", "data.frame"))
}

# Freeman degree centralization. For directed graphs the maximum possible
# sum of (max - d) across n nodes is (n - 1)^2 - the star graph attains it.
# For undirected graphs the corresponding bound is (n - 1)(n - 2).
.freeman_degree_centralization <- function(d, n, directed) {
  if (n < 2L) return(0)
  num <- sum(max(d) - d)
  denom <- if (directed) (n - 1)^2 else (n - 1) * (n - 2)
  if (denom <= 0) return(0)
  num / denom
}


# ---- print methods ----

#' @rdname build_network
#' @export
print.summary.netobject <- function(x, ...) {
  cat("Network metrics:\n")
  out <- x
  out$value <- .format_metric_values(out$metric, out$value)
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}

# Counts print as integers; other metrics as 4 significant digits, formatted
# value-by-value so `format()` cannot promote the whole column to scientific
# notation when one cell is near machine epsilon.
.format_metric_values <- function(metric, value) {
  is_count <- metric %in% c("Node Count", "Edge Count")
  vapply(seq_along(value), function(i) {
    if (is_count[i]) formatC(value[i], format = "d", big.mark = "")
    else formatC(value[i], digits = 4L, format = "g")
  }, character(1L))
}

#' @rdname build_network
#' @export
print.summary.netobject_group <- function(x, ...) {
  if (isFALSE(attr(x, "combined"))) {
    invisible(lapply(names(x), function(nm) {
      cat("Group:", nm, "\n")
      print(x[[nm]])
      cat("\n")
    }))
    return(invisible(x))
  }
  cat("Network metrics by group:\n")
  out <- x
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols],
                              function(v) .format_metric_values(out$metric, v))
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}
