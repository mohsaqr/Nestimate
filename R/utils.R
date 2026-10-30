# Internal helpers. No Rd page: nothing here is exported. The roxygen block
# below exists only to carry the NAMESPACE imports.

#' @importFrom utils tail head
NULL

# Polyfill `%||%` for R < 4.4. Base R 4.4 added it; before then it lived only
# in rlang/purrr. Many internal call sites already use it. Defined unconditionally
# because Nestimate's namespace looks here first; base's version on R >= 4.4 is
# functionally identical, so this is a harmless shadow.
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Coerce a labels argument to a named character lookup.
#' Accepts a 2-col data.frame `(name, label)`, a named character vector,
#' or a named list. Returns a named character vector or NULL on empty input.
#' @noRd
.coerce_label_map <- function(labels) {
  if (is.null(labels)) return(NULL)
  if (is.data.frame(labels)) {
    if (ncol(labels) < 2L) {
      stop("`labels` data.frame must have at least two columns.",
           call. = FALSE)
    }
    map <- setNames(as.character(labels[[2]]), as.character(labels[[1]]))
    return(.validate_label_map(map))
  }
  if (is.list(labels) && !is.null(names(labels))) {
    return(.validate_label_map(unlist(labels)))
  }
  if (is.character(labels) && !is.null(names(labels))) {
    return(.validate_label_map(labels))
  }
  stop("`labels` must be a named character vector, named list, or 2-column ",
       "data.frame (name, label).", call. = FALSE)
}

.validate_label_map <- function(map) {
  key <- names(map)
  map <- as.character(map)
  names(map) <- key
  if (is.null(key) || any(is.na(key)) || any(!nzchar(key))) {
    stop("`labels` must have non-empty source names.", call. = FALSE)
  }
  if (anyDuplicated(key)) {
    stop("`labels` source names must be unique.", call. = FALSE)
  }
  if (any(is.na(map)) || any(!nzchar(map))) {
    stop("`labels` values must be non-missing and non-empty.",
         call. = FALSE)
  }
  if (anyDuplicated(map)) {
    stop("`labels` values must be unique.", call. = FALSE)
  }
  map
}

#' Apply a name -> label remap to a netobject or mcml object.
#'
#' For a netobject: rewrites `$nodes$label` and the dimnames of `$weights`.
#' For an mcml object: rewrites within-cluster `$labels`, weight dimnames,
#' init names, `$cluster_members`, and `$edges$from`/`$to`. Macro layer is
#' left untouched (its labels are cluster names, not node names).
#' Unmapped names pass through unchanged.
#' @noRd
.apply_node_labels <- function(x, labels) {
  map <- .coerce_label_map(labels)
  if (is.null(map)) return(x)

  remap <- function(v) {
    v <- as.character(v)
    out <- map[v]; out[is.na(out)] <- v[is.na(out)]; unname(out)
  }
  check_unique <- function(v) {
    if (anyDuplicated(v)) {
      stop("`labels` must map nodes to unique labels.", call. = FALSE)
    }
    v
  }

  # Apply the label remap to the values inside a network's $data so
  # they line up with the relabeled $weights rownames. Without this,
  # bootstrap_network() in its precomputed fast path silently
  # produces all-zero matrices because match($data, states) returns
  # NA for every cell. Values not in `map` (e.g. macro $data carries
  # cluster names) are preserved.
  remap_data <- function(d, columns = colnames(d) %||% seq_len(NCOL(d))) {
    if (is.null(d)) return(d)
    # Numeric sequence states are literal labels, unlike association
    # measurements and one-hot indicators (which never enter this path).
    # When any column holds a mapped state, every numeric state column
    # becomes character: converting only those columns would leave the
    # others numeric, and as.matrix() then pads them (" 2"), so the
    # resampler's match() against the labels silently drops those cells.
    cols_in <- if (is.data.frame(d)) intersect(columns, names(d)) else
      intersect(columns, colnames(d) %||% seq_len(NCOL(d)))
    holds_states <- any(vapply(cols_in, \(j) {
      v <- if (is.data.frame(d)) d[[j]] else d[, j]
      any(as.character(v) %in% names(map))
    }, logical(1)))
    remap_column <- function(col) {
      if (is.factor(col)) {
        levels(col) <- remap(levels(col))
        return(col)
      }
      if (is.character(col) || holds_states) {
        return(remap(col))
      }
      col
    }
    if (is.data.frame(d)) {
      columns <- intersect(columns, names(d))
      d[columns] <- lapply(d[columns], remap_column)
      return(d)
    }
    if (is.matrix(d)) {
      columns <- if (is.character(columns)) {
        match(columns, colnames(d), nomatch = 0L)
      } else columns
      columns <- columns[columns > 0L & columns <= ncol(d)]
      frame <- as.data.frame(d, stringsAsFactors = FALSE)
      frame[columns] <- lapply(frame[columns], remap_column)
      out <- as.matrix(frame)
      dimnames(out) <- dimnames(d)
      return(out)
    }
    d
  }
  if (inherits(x, "mcml")) {
    relabel_one <- function(net) {
      if (is.null(net)) return(net)
      net$labels <- check_unique(remap(net$labels))
      if (!is.null(net$weights))
        dimnames(net$weights) <- list(net$labels, net$labels)
      if (!is.null(net$inits)) names(net$inits) <- net$labels
      net$data <- remap_data(net$data)
      net
    }
    if (!is.null(x$clusters))        x$clusters        <- lapply(x$clusters, relabel_one)
    if (!is.null(x$cluster_members)) x$cluster_members <- lapply(x$cluster_members, remap)
    if (!is.null(x$edges)) {
      x$edges$from <- remap(x$edges$from)
      x$edges$to   <- remap(x$edges$to)
    }
    return(x)
  }

  if (inherits(x, "netobject")) {
    original_states <- rownames(x$weights) %||% colnames(x$weights)
    if (!is.null(x$nodes) && "name" %in% names(x$nodes)) {
      x$nodes$label <- check_unique(remap(as.character(x$nodes$name)))
    }
    if (!is.null(x$weights)) {
      nm <- rownames(x$weights) %||% colnames(x$weights)
      if (!is.null(nm)) {
        mapped <- check_unique(remap(nm))
        dimnames(x$weights) <- list(mapped, mapped)
      }
    }
    # A wide sequence stores nodes in its cells, while association and
    # one-hot inputs store nodes in column names. Renaming positional wide
    # columns would invalidate params$cols, and remapping numeric indicator
    # or association values would change the observations themselves.
    sequence_methods <- c("relative", "frequency", "co_occurrence",
                          "attention", "ngram", "gap", "reverse")
    column_nodes <- !x$method %in% sequence_methods ||
      !is.null(.param_get(x$params, "codes")) ||
      identical(.param_get(x$params, "format"), "onehot")
    if (column_nodes) {
      old_columns <- colnames(x$data)
      if (!is.null(old_columns)) {
        new_columns <- old_columns
        is_node <- old_columns %in% original_states
        new_columns[is_node] <- remap(old_columns[is_node])
        colnames(x$data) <- check_unique(new_columns)
        rename_columns <- function(v) {
          if (is.null(v)) return(NULL)
          idx <- match(v, old_columns)
          matched <- !is.na(idx)
          v[matched] <- new_columns[idx[matched]]
          v
        }
        column_roles <- c("cols", "codes", "action", "id", "id_col",
                          "actor", "time", "session", "order", "group")
        update_roles <- function(args) {
          roles <- intersect(column_roles, names(args))
          args[roles] <- lapply(args[roles], rename_columns)
          args
        }
        x$params <- update_roles(x$params)
        if (!is.null(x$build_args)) x$build_args <- update_roles(x$build_args)
      }
    } else {
      id_columns <- unique(c(.param_get(x$params, "id"),
                             .param_get(x$params, "id_col"),
                             .param_get(x$params, "actor")))
      columns <- .param_get(x$params, "cols") %||%
        setdiff(colnames(x$data), id_columns)
      x$data <- remap_data(x$data, columns)
      state_roles <- intersect(c("alphabet", "begin_state", "end_state"),
                               names(x$params))
      x$params[state_roles] <- lapply(x$params[state_roles], remap)
    }
    if (!is.null(x$initial)) names(x$initial) <- remap(names(x$initial))
    if (!is.null(x$predictability)) {
      names(x$predictability) <- remap(names(x$predictability))
    }
    if (!is.null(x$frequency_matrix)) {
      dimnames(x$frequency_matrix) <- lapply(dimnames(x$frequency_matrix),
                                             function(d) if (is.null(d)) d else remap(d))
    }
    if (is.data.frame(x$node_groups) && "node" %in% names(x$node_groups)) {
      x$node_groups$node <- remap(as.character(x$node_groups$node))
    }
    return(x)
  }

  x
}


#' Convert a named numeric matrix to a long tidy data.frame.
#'
#' Used by `summary()` methods on class-stamped matrix returns to produce a
#' `(from, to, <value>)` edge-list view, keeping row/column names as strings.
#'
#' @param m A numeric matrix.
#' @param value_col Name of the value column (e.g. "weight", "count").
#' @param include One of `"nonzero"` (default) or `"positive"` - which entries
#'   to include.
#' @param sort_by One of `"abs_value"` (default, descending) or `"none"`.
#' @return A data.frame with columns `from`, `to`, and `<value_col>`.
#' @noRd
.matrix_to_long_df <- function(m, value_col = "weight",
                               include = c("nonzero", "positive"),
                               sort_by = c("abs_value", "none")) {
  include <- match.arg(include)
  sort_by <- match.arg(sort_by)
  empty_df <- function() {
    out <- data.frame(from = character(0L), to = character(0L),
                      x = numeric(0L), stringsAsFactors = FALSE)
    names(out)[3L] <- value_col
    out
  }
  if (!is.matrix(m) || nrow(m) == 0L || ncol(m) == 0L) return(empty_df())
  idx <- switch(include,
                nonzero  = which(m != 0, arr.ind = TRUE),
                positive = which(m > 0,  arr.ind = TRUE))
  if (nrow(idx) == 0L) return(empty_df())
  rn <- rownames(m); if (is.null(rn)) rn <- as.character(seq_len(nrow(m)))
  cn <- colnames(m); if (is.null(cn)) cn <- as.character(seq_len(ncol(m)))
  vals <- m[idx]
  df <- data.frame(
    from = rn[idx[, 1L]],
    to   = cn[idx[, 2L]],
    x    = if (value_col == "count") as.integer(vals) else as.numeric(vals),
    stringsAsFactors = FALSE,
    row.names        = NULL
  )
  names(df)[3L] <- value_col
  if (sort_by == "abs_value") df <- df[order(-abs(df[[3L]])), ]
  row.names(df) <- NULL
  df
}

# Global variable declarations to avoid R CMD check notes
utils::globalVariables(c(
  # Common column names
  ".", ":=", ".I", ".SD", ".N",
  "cr_lower", "cr_upper", "effect_size", "sig",
  "from", "to", "value", "run_id", "setting_name",
  # ggplot2 NSE: density is a computed stat from geom_histogram(after_stat(density))
  "density",
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
safe_sd <- function(x) {
  if (length(x) > 1) sd(x, na.rm = TRUE) else NA_real_
}


#' Coerce tna or netobject to labeled sequence data.frame
#'
#' When \code{data} is a \code{tna} or \code{netobject}, extracts the
#' sequence data and decodes positional IDs only for \code{tna} input.
#' A \code{netobject} stores literal states, including numeric states. This
#' allows \code{build_hon()}, \code{build_hypa()}, and other pathway
#' functions to accept model objects directly.
#'

#' @param data Input: data.frame, list, tna, or netobject.
#' @return A data.frame or list suitable for \code{.hon_parse_input()}.
#' @noRd
.coerce_sequence_input <- function(data) {
  if (inherits(data, "tna")) {
    if (is.null(data$data)) { # nocov start
      stop("tna object has no sequence data ($data). ",
           "Build the tna from sequence data, not a raw matrix.",
           call. = FALSE)
    } # nocov end
    df <- as.data.frame(data$data, stringsAsFactors = FALSE) # nocov start
    lbl <- attr(data$data, "labels") %||% data$labels
    if (!is.null(lbl) && length(lbl) > 0L) {
      df[] <- lapply(df, function(col) {
        if (!is.numeric(col)) return(col)
        idx <- as.integer(col)
        valid <- !is.na(idx) & idx >= 1L & idx <= length(lbl)
        out <- rep(NA_character_, length(idx))
        out[valid] <- lbl[idx[valid]]
        out
      })
      attr(df, "labels") <- NULL
      attr(df, "alphabet") <- NULL
    }
    return(df) # nocov end
  }
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }
  if (inherits(data, "netobject")) {
    if (is.null(data$data)) { # nocov start
      stop("netobject has no sequence data ($data). ",
           "Build the network from sequence data.",
           call. = FALSE) # nocov end
    }
    return(as.data.frame(data$data, stringsAsFactors = FALSE))
  }
  ## Bare sequence matrix (character / logical) -> wide data.frame.
  if (is.matrix(data) && !is.numeric(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  data
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
  if (is.null(mat)) {
    stop("cograph_network has no $weights matrix.", call. = FALSE)
  }
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (!is.numeric(mat)) storage.mode(mat) <- "double"
  nodes_df <- x$nodes
  states <- nodes_df$label
  raw_data <- x$data
  directed <- x$directed %||% TRUE

  # Infer method from tna metadata or matrix symmetry
  tna_meta <- x$meta$tna
  method <- if (!is.null(tna_meta$method)) {
    tna_meta$method
  } else if (is.matrix(mat) && isSymmetric(mat)) {
    "co_occurrence"
  } else {
    "relative"
  }

  is_sequence_method <- method %in% c(
    "relative", "frequency", "co_occurrence", "attention"
  )

  # Numeric cells alone do not identify encoded states: generic networks
  # may store literal numeric labels. A converted tna model retains its
  # data class/alphabet or its explicit source metadata.
  if (!is.null(raw_data)) {
    encoded_labels <- attr(raw_data, "labels") %||% attr(raw_data, "alphabet")
    encoded <- inherits(raw_data, "tna_seq_data") ||
      !is.null(encoded_labels) || identical(x$meta$source, "tna")
    raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
    if (is_sequence_method && encoded) {
      labels <- encoded_labels %||% states
      raw_data[] <- lapply(raw_data, function(col) {
        if (!is.numeric(col)) return(col)
        idx <- as.integer(col)
        valid <- !is.na(idx) & idx >= 1L & idx <= length(labels)
        out <- rep(NA_character_, length(idx))
        out[valid] <- labels[idx[valid]]
        out
      })
      attr(raw_data, "labels") <- NULL
      attr(raw_data, "alphabet") <- NULL
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


# ---------------------------------------------------------------------------
# Higher-order -> cograph_network bridge
# ---------------------------------------------------------------------------

#' Add cograph_network fields to a higher-order network object
#'
#' @param mat Square weight matrix with named rows/columns.
#' @param node_names Character vector of node names.
#' @param method Character. Method label for metadata.
#' @return Named list with \code{weights}, \code{nodes} (data.frame),
#'   \code{edges}, \code{directed}, \code{meta} fields.
#' @noRd
.ho_cograph_fields <- function(mat, node_names, method = "hon") {
  nodes_df <- data.frame(
    id = seq_along(node_names),
    label = node_names,
    name = node_names,
    stringsAsFactors = FALSE
  )
  edges <- .extract_edges_from_matrix(mat, directed = TRUE)
  list(
    weights = mat,
    nodes = nodes_df,
    edges = edges,
    directed = TRUE,
    n_nodes = length(node_names),
    n_edges = nrow(edges),
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna = list(method = method)
    ),
    node_groups = NULL
  )
}


#' Report resamples whose network could not be estimated
#'
#' Resamplers catch estimator errors per draw so one degenerate draw does
#' not abort the run; the count is reported here instead of being dropped
#' silently.
#' @noRd
.warn_failed_resamples <- function(n_failed, n_total, fn) {
  if (n_failed > 0L) {
    warning(warningCondition(sprintf(paste0(
      "%s: %d of %d resamples failed to estimate a network and are left ",
      "out%s."), fn, n_failed, n_total,
      if (n_failed == n_total) " (all of them: the results are NA)" else ""),
      class = "nestimate_resample_failed", call = NULL))
  }
  invisible(n_failed)
}

