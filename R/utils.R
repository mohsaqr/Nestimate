#' @title Internal Helper Functions for Nestimate
#' @name utils
#' @description Internal utility functions used by other Nestimate functions.
#' @keywords internal
#' @importFrom utils tail head
NULL

# Global variable declarations to avoid R CMD check notes
utils::globalVariables(c(
  # Common column names
  ".", ":=", ".I", ".SD", ".N",
  "cr_lower", "cr_upper", "effect_size", "sig",
  "from", "to", "value", "run_id", "setting_name",
  # Metrics
  "TP", "TN", "FP", "FN", "Sensitivity", "Specificity", "FPR", "FNR",
  "Accuracy", "MCC", "mcc_denom_sq", "Metric", "Value",
  # Totals
  "Total_TP", "Total_TN", "Total_FP", "Total_FN",
  # Bootstrap/simulation
  "p_value", "weight", "p_value_num", "weight_num",
  "bootstrap_significant_run", "is_significant",
  "n_significant", "recovery_rate", "avg_recovery_rate",
  "avg_p_value", "avg_weight",
  "ground_truth_stable", "gt_stable",
  # Parameters
  "num_rows", "max_seq_length", "min_na", "max_na", "num_states",
  "successful_runs",
  # Network simulation
  "category", "metric", "model_type", "comparison_type",
  "metric_category", "metric_name", "data_idx",
  # Data conversion
  "id", "Time", "Action",
  # Summary functions
  ".n", ".sd", "network_id"
))

#' Check Value in Range
#'
#' Check if a value falls within a specified range.
#'
#' @param value Numeric value to check.
#' @param range_val Numeric vector of length 2 with min and max, or NULL.
#'
#' @return Logical indicating whether value is in range.
#'
#' @keywords internal
check_val_in_range <- function(value, range_val) {
  if (is.null(range_val)) return(TRUE)
  if (is.na(value) || !is.numeric(value)) return(FALSE)
  return(value >= range_val[1] && value <= range_val[2])
}

#' Safe Median
#'
#' Calculate median with handling for empty vectors.
#'
#' @param x Numeric vector.
#'
#' @return Median value or NA if vector is empty.
#'
#' @keywords internal
safe_median <- function(x) {
  if (length(x) > 0) median(x, na.rm = TRUE) else NA_real_
}

#' Safe Mean
#'
#' Calculate mean with handling for empty vectors.
#'
#' @param x Numeric vector.
#'
#' @return Mean value or NA if vector is empty.
#'
#' @keywords internal
safe_mean <- function(x) {
  if (length(x) > 0) mean(x, na.rm = TRUE) else NA_real_
}

#' Safe Standard Deviation
#'
#' Calculate standard deviation with handling for single-value vectors.
#'
#' @param x Numeric vector.
#'
#' @return Standard deviation or NA if vector has fewer than 2 elements.
#'
#' @keywords internal
safe_sd <- function(x) {
  if (length(x) > 1) sd(x, na.rm = TRUE) else NA_real_
}


#' Convert pure cograph_network to dual-class netobject/cograph_network
#'
#' Internal converter so that downstream functions (bootstrap, permutation,
#' reliability, etc.) can accept either \code{netobject} or
#' \code{cograph_network} inputs transparently. Objects that already have
#' the \code{"netobject"} class are returned unchanged.
#'
#' @param x A \code{netobject} (returned unchanged) or \code{cograph_network}.
#' @return A dual-class \code{c("netobject", "cograph_network")} object.
#' @noRd
.as_netobject <- function(x) {
  if (inherits(x, "netobject")) return(x)
  if (!inherits(x, "cograph_network")) {
    stop("Expected a netobject or cograph_network.", call. = FALSE)
  }

  mat <- x$weights
  nodes_df <- x$nodes
  states <- nodes_df$label
  raw_data <- x$data
  directed <- x$directed %||% TRUE

  # Infer method from tna metadata or matrix symmetry
  tna_meta <- x$meta$tna
  method <- if (!is.null(tna_meta$method)) {
    tna_meta$method
  } else if (isSymmetric(mat)) {
    "co_occurrence"
  } else {
    "relative"
  }

  is_sequence_method <- method %in% c(
    "relative", "frequency", "co_occurrence", "attention"
  )

  # Decode integer-encoded tna data -> character labels
  # Only for sequence methods; association methods keep numeric data as-is
  if (!is.null(raw_data)) {
    raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
    if (is_sequence_method &&
        (is.integer(raw_data[[1]]) || is.numeric(raw_data[[1]]))) {
      raw_data[] <- lapply(raw_data, function(col) {
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(states),
               NA_character_, states[idx])
      })
    }
  }

  edges <- .extract_edges_from_matrix(mat, directed = directed)

  structure(list(
    data = raw_data, weights = mat, nodes = nodes_df,
    edges = edges, directed = directed, method = method,
    params = list(), scaling = NULL, threshold = 0,
    n_nodes = length(states), n_edges = nrow(edges),
    level = NULL,
    meta = x$meta %||% list(source = "cograph", layout = NULL,
                            tna = list(method = method)),
    node_groups = x$node_groups
  ), class = c("netobject", "cograph_network"))
}


#' Convert to cograph_network
#'
#' Since \code{netobject} now inherits from \code{cograph_network}, this
#' function returns the object as-is for netobjects. For other inputs it
#' uses \code{cograph::cograph()} if the cograph package is available.
#'
#' @param x An object to convert.
#' @param ... Additional arguments passed to \code{cograph::cograph()}.
#' @return A \code{cograph_network}.
#'
#' @export
as_cograph <- function(x, ...) {
  UseMethod("as_cograph")
}

#' @rdname as_cograph
#' @export
as_cograph.netobject <- function(x, ...) {
  # Already a cograph_network via dual class
  x
}

#' @rdname as_cograph
#' @export
as_cograph.cograph_network <- function(x, ...) {
  x
}

#' @rdname as_cograph
#' @export
as_cograph.netobject_ml <- function(x, ...) {
  list(
    between = as_cograph(x$between, ...),
    within = as_cograph(x$within, ...)
  )
}

#' @rdname as_cograph
#' @export
as_cograph.netobject_group <- function(x, ...) {
  result <- lapply(x, as_cograph, ...)
  names(result) <- names(x)
  result
}

#' @rdname as_cograph
#' @export
as_cograph.default <- function(x, ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop("Package 'cograph' is required.", call. = FALSE)
  }
  cograph::cograph(x, ...)
}
