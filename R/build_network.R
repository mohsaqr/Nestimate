#' Build a Network
#'
#' @description
#' Universal network estimation function that supports both transition
#' networks (relative, frequency, co-occurrence) and association networks
#' (correlation, partial correlation, graphical lasso). Uses the global
#' estimator registry, so custom estimators can also be used.
#'
#' @param data Data frame (sequences or per-observation frequencies) or a
#'   square symmetric matrix (correlation or covariance).
#' @param method Character. Required. Name of a registered estimator.
#'   Built-in methods: \code{"relative"}, \code{"frequency"},
#'   \code{"co_occurrence"}, \code{"cor"}, \code{"pcor"}, \code{"glasso"}.
#'   Aliases: \code{"tna"} and \code{"transition"} map to \code{"relative"};
#'   \code{"ftna"} and \code{"counts"} map to \code{"frequency"};
#'   \code{"cna"} maps to \code{"co_occurrence"};
#'   \code{"corr"} and \code{"correlation"} map to \code{"cor"};
#'   \code{"partial"} maps to \code{"pcor"};
#'   \code{"ebicglasso"} and \code{"regularized"} map to \code{"glasso"}.
#' @param params Named list. Method-specific parameters passed to the estimator
#'   function (e.g. \code{list(gamma = 0.5)} for glasso, or
#'   \code{list(format = "wide")} for transition methods). This is the key
#'   composability feature: downstream functions like bootstrap or grid search
#'   can store and replay the full params list without knowing method internals.
#' @param scaling Character vector or NULL. Post-estimation scaling to apply
#'   (in order). Options: \code{"minmax"}, \code{"max"}, \code{"rank"},
#'   \code{"normalize"}. Can combine: \code{c("rank", "minmax")}.
#'   Default: \code{NULL} (no scaling).
#' @param threshold Numeric. Absolute values below this are set to zero in the
#'   result matrix. Default: 0 (no thresholding).
#' @param level Character or NULL. Multilevel decomposition for association
#'   methods. One of \code{NULL}, \code{"between"}, \code{"within"},
#'   \code{"both"}. Requires \code{id_col}. Default: \code{NULL}.
#' @param actor Character. Name of the actor/person ID column for sequence
#'   grouping. Default: \code{NULL}.
#' @param action Character. Name of the action/state column (long format).
#'   Default: \code{NULL}.
#' @param time Character. Name of the time column (long format).
#'   Default: \code{NULL}.
#' @param session Character. Name of the session column. Default: \code{NULL}.
#' @param order Character. Name of the ordering column. Default: \code{NULL}.
#' @param codes Character vector. Column names of one-hot encoded states
#'   (for onehot format). Default: \code{NULL}.
#' @param group Character. Name of a grouping column for per-group networks.
#'   Returns a \code{netobject_group} (named list of netobjects).
#'   Default: \code{NULL}.
#' @param format Character. Input format: \code{"auto"}, \code{"wide"},
#'   \code{"long"}, or \code{"onehot"}. Default: \code{"auto"}.
#' @param window_size Integer. Window size for one-hot windowing.
#'   Default: \code{3L}.
#' @param mode Character. Windowing mode: \code{"non-overlapping"} or
#'   \code{"overlapping"}. Default: \code{"non-overlapping"}.
#' @param time_threshold Numeric. Maximum time gap (seconds) for long format
#'   session splitting. Default: \code{900}.
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return An object of class \code{c("netobject", "cograph_network")} containing:
#' \describe{
#'   \item{data}{The input data used for estimation, as a data frame.}
#'   \item{weights}{The estimated network weight matrix.}
#'   \item{nodes}{Data frame with columns \code{id}, \code{label}, \code{name},
#'     \code{x}, \code{y}. Node labels are in \code{$nodes$label}.}
#'   \item{edges}{Data frame of non-zero edges with integer \code{from}/\code{to}
#'     (node IDs) and numeric \code{weight}.}
#'   \item{directed}{Logical. Whether the network is directed.}
#'   \item{method}{The resolved method name.}
#'   \item{params}{The params list used (for reproducibility).}
#'   \item{scaling}{The scaling applied (or NULL).}
#'   \item{threshold}{The threshold applied.}
#'   \item{n_nodes}{Number of nodes.}
#'   \item{n_edges}{Number of non-zero edges.}
#'   \item{level}{Decomposition level used (or NULL).}
#'   \item{meta}{List with \code{source}, \code{layout}, and \code{tna} metadata
#'     (cograph-compatible).}
#'   \item{node_groups}{Node groupings data frame, or NULL.}
#' }
#' Method-specific extras (e.g. \code{precision_matrix}, \code{cor_matrix},
#' \code{frequency_matrix}, \code{lambda_selected}, etc.) are preserved
#' from the estimator output.
#'
#' When \code{level = "both"}, returns an object of class
#' \code{"netobject_ml"} with \code{$between} and \code{$within}
#' sub-networks and a \code{$method} field.
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Resolves method aliases to canonical names.
#'   \item Retrieves the estimator function from the global registry.
#'   \item For association methods with \code{level} specified, decomposes
#'     the data (between-person means or within-person centering).
#'   \item Calls the estimator: \code{do.call(fn, c(list(data = data), params))}.
#'   \item Applies scaling and thresholding to the result matrix.
#'   \item Extracts edges and constructs the \code{netobject}.
#' }
#'
#' @examples
#' \dontrun{
#' library(tna)
#'
#' # Transition network (relative probabilities)
#' net <- build_network(group_regulation, method = "relative")
#' print(net)
#'
#' # Aliases
#' net_tna <- build_network(group_regulation, method = "tna")
#' net_ftna <- build_network(group_regulation, method = "ftna")
#' net_cna <- build_network(group_regulation, method = "cna")
#'
#' # Association network (glasso)
#' freq_data <- convert_sequence_format(group_regulation, format = "frequency")
#' net_glasso <- build_network(freq_data, method = "glasso",
#'                              params = list(gamma = 0.5, nlambda = 50))
#'
#' # Partial correlation network
#' net_pcor <- build_network(freq_data, method = "pcor")
#'
#' # Correlation network with alias
#' net_cor <- build_network(freq_data, method = "corr")
#'
#' # With scaling
#' net_scaled <- build_network(group_regulation, method = "relative",
#'                              scaling = c("rank", "minmax"))
#'
#' # Composable: replay config on new data
#' config <- net_glasso$params
#' net2 <- build_network(new_data, method = net_glasso$method,
#'                        params = config)
#' }
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{bootstrap_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @export
build_network <- function(data,
                          method,
                          actor = NULL,
                          action = NULL,
                          time = NULL,
                          session = NULL,
                          order = NULL,
                          codes = NULL,
                          group = NULL,
                          format = "auto",
                          window_size = 3L,
                          mode = c("non-overlapping", "overlapping"),
                          scaling = NULL,
                          threshold = 0,
                          level = NULL,
                          time_threshold = 900,
                          params = list(),
                          ...) {
  # --- Early dispatch for net_clustering objects ---
  if (inherits(data, "net_clustering")) {
    if (missing(method)) method <- "relative"
    return(.build_network_clustering(data, method = method, ...))
  }

  stopifnot(is.character(method), length(method) == 1)
  stopifnot(is.list(params))
  stopifnot(is.numeric(threshold), length(threshold) == 1, threshold >= 0)
  mode <- match.arg(mode)

  # Merge ... into params (... takes precedence)
  dots <- list(...)
  if (length(dots) > 0) {
    params <- modifyList(params, dots)
  }

  # Resolve method aliases early (needed for format detection)
  method <- .resolve_method_alias(method)

  # ---- Group dispatch: per-group networks ----
  if (!is.null(group)) {
    stopifnot(is.character(group), all(group %in% names(data)))
    if (length(group) == 1L) {
      grp_key <- data[[group]]
    } else {
      grp_key <- interaction(data[, group, drop = FALSE], sep = "-",
                             drop = TRUE)
    }
    grp_levels <- unique(grp_key)
    # Drop group column(s) from sub-data so they don't become state cols
    drop_cols <- intersect(group, names(data))
    nets <- lapply(grp_levels, function(g) {
      sub <- data[grp_key == g, , drop = FALSE]
      if (length(drop_cols) > 0L) {
        sub <- sub[, setdiff(names(sub), drop_cols), drop = FALSE]
      }
      build_network(
        sub, method = method, actor = actor, action = action,
        time = time, session = session, order = order, codes = codes,
        group = NULL, format = format, window_size = window_size,
        mode = mode, scaling = scaling, threshold = threshold,
        level = level, time_threshold = time_threshold,
        params = params, ...
      )
    })
    names(nets) <- as.character(grp_levels)
    attr(nets, "group_col") <- group
    class(nets) <- "netobject_group"
    return(nets)
  }

  # ---- Auto-detect input format ----
  is_onehot <- FALSE
  if (format == "auto" && is.data.frame(data)) {
    if (!is.null(codes)) {
      # Explicit codes = one-hot
      is_onehot <- TRUE
      format <- "onehot"
    } else if (!is.null(action) && action %in% names(data)) {
      # Has action column = long event log
      format <- "long"
    } else {
      # Check if all columns (excluding actor/session) are binary
      exclude <- c(actor, session)
      check_cols <- setdiff(names(data), exclude)
      if (length(check_cols) >= 2L && all(vapply(
        data[, check_cols, drop = FALSE],
        function(x) is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)),
        logical(1)
      ))) {
        is_onehot <- TRUE
        format <- "onehot"
      } else {
        format <- "wide"
      }
    }
  }

  # ---- Long format: prepare event log data ----
  if (format == "long" && !is.null(action) && is.data.frame(data) &&
      action %in% names(data)) {
    prep_args <- list(data = data, action = action)
    if (!is.null(actor)) prep_args$actor <- actor
    if (!is.null(time)) prep_args$time <- time
    if (!is.null(session)) prep_args$session <- session
    if (!is.null(order)) prep_args$order <- order
    prep_args$time_threshold <- time_threshold

    prepared <- do.call(prepare_data, prep_args)
    data <- prepared$sequence_data
    format <- "wide"
    params$format <- "wide"
  } else {
    prepared <- NULL
  }

  # ---- One-hot format: dispatch to wtna/cna path ----
  if (format == "onehot" || is_onehot) {
    resolved_codes <- .resolve_codes(
      as.data.frame(data), codes,
      exclude = c(actor, session)
    )

    # Warning for one-hot without windowing or session
    if (window_size <= 1L && is.null(session)) {
      warning(
        "One-hot data without windowing or sessions produces sparse networks. ",
        "Consider setting window_size > 1 or providing session.",
        call. = FALSE
      )
    }

    # Build grouping from actor + session
    grp_col <- NULL
    if (!is.null(actor) || !is.null(session)) {
      grp_parts <- c(actor, session)
      grp_col <- grp_parts
    }

    params$codes <- resolved_codes
    params$window_size <- window_size
    params$mode <- mode
    if (!is.null(grp_col)) params$actor <- grp_col
    # Let the estimator handle it
  } else {
    # Wide sequence data: pass format through
    if (!"format" %in% names(params)) params$format <- format
  }

  # ---- Multilevel decomposition ----
  id_col <- params$id %||% actor

  # Validate level parameter
  if (!is.null(level)) {
    level <- match.arg(level, c("between", "within", "both"))
    if (is.null(id_col)) {
      stop("'id' is required when 'level' is specified.", call. = FALSE)
    }
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame when 'level' is specified.",
           call. = FALSE)
    }
  }

  # Validate scaling
  if (!is.null(scaling)) {
    valid_scaling <- c("minmax", "max", "rank", "normalize")
    bad <- setdiff(scaling, valid_scaling)
    if (length(bad) > 0) {
      stop("Unknown scaling method(s): ", paste(bad, collapse = ", "),
           ". Options: ", paste(valid_scaling, collapse = ", "),
           call. = FALSE)
    }
  }

  # Get estimator from registry
  estimator <- get_estimator(method)

  # level = "both": recursive dispatch
  if (identical(level, "both")) {
    between <- build_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "between"
    )
    within_net <- build_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "within"
    )
    result <- list(between = between, within = within_net, method = method)
    class(result) <- "netobject_ml"
    return(result)
  }

  # Multilevel decomposition for association methods
  if (!is.null(level) && !estimator$directed) {
    data <- .decompose_multilevel(data, id_col = id_col, level = level)
  }

  # Call estimator
  est_result <- do.call(estimator$fn, c(list(data = data), params))

  # Validate estimator output
  if (!is.list(est_result) ||
      is.null(est_result$matrix) ||
      is.null(est_result$nodes) ||
      is.null(est_result$directed)) {
    stop("Estimator '", method,
         "' must return a list with 'matrix', 'nodes', and 'directed'.",
         call. = FALSE)
  }

  net_matrix <- est_result$matrix
  nodes <- est_result$nodes
  directed <- est_result$directed

  # Apply scaling
  if (!is.null(scaling)) {
    net_matrix <- .apply_scaling(net_matrix, scaling)
  }

  # Apply threshold
  if (threshold > 0) {
    net_matrix[abs(net_matrix) < threshold] <- 0
  }

  # Extract edges
  edges <- .extract_edges_from_matrix(net_matrix, directed = directed)

  # Split cleaned_data into state-only $data and $metadata
  # A column is a state column if all its non-void/non-NA values are in nodes
  raw_data <- est_result$cleaned_data
  metadata <- NULL
  if (is.data.frame(raw_data)) {
    is_state_col <- vapply(raw_data, function(col) {
      vals <- .clean_states(as.character(col))
      vals <- vals[!is.na(vals)]
      length(vals) > 0L && all(vals %in% nodes)
    }, logical(1))
    state_cols <- names(raw_data)[is_state_col]
    extra_cols <- names(raw_data)[!is_state_col]
    if (length(extra_cols) > 0L) {
      metadata <- raw_data[, extra_cols, drop = FALSE]
      raw_data <- raw_data[, state_cols, drop = FALSE]
    }
    # Clean void/missing markers in character/factor state columns
    if (is.data.frame(raw_data)) {
      raw_data[] <- lapply(raw_data, function(col) {
        if (is.character(col) || is.factor(col)) {
          .clean_states(as.character(col))
        } else {
          col
        }
      })
    }
  }

  # Build unified netobject / cograph_network
  nodes_df <- data.frame(
    id = seq_along(nodes),
    label = nodes,
    name = nodes,
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
  )

  result <- list(
    data = raw_data,
    metadata = metadata,
    weights = net_matrix,
    nodes = nodes_df,
    edges = edges,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = length(nodes),
    n_edges = nrow(edges),
    level = level,
    prepared = prepared,
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna = list(method = method)
    ),
    node_groups = NULL
  )

  # Carry over method-specific extras
  extras <- setdiff(names(est_result), c("matrix", "nodes", "directed"))
  for (key in extras) {
    result[[key]] <- est_result[[key]]
  }

  structure(result, class = c("netobject", "cograph_network"))
}


# ---- S3 methods ----

#' Print Method for Network Object
#'
#' @param x A \code{netobject}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.netobject <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    ising         = "Ising Model Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)"
  )

  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label <- if (x$directed) " [directed]" else " [undirected]"
  level_label <- if (!is.null(x$level)) {
    sprintf(" [%s-person]", x$level)
  } else ""

  cat(label, dir_label, level_label, "\n", sep = "")

  # ---- Data overview ----
  if (!is.null(x$data)) {
    cat(sprintf("  Data: %d sequences x %d time points\n",
                nrow(x$data), ncol(x$data)))
  }
  if (!is.null(x$metadata)) {
    cat(sprintf("  Metadata: %s\n",
                paste(names(x$metadata), collapse = ", ")))
  }
  if (!is.null(x$n)) cat(sprintf("  Sample size: %d\n", x$n))

  # ---- Nodes ----
  labels <- x$nodes$label
  cat(sprintf("  Nodes (%d): %s\n", x$n_nodes,
              if (x$n_nodes <= 8) paste(labels, collapse = ", ")
              else paste(c(labels[1:6], sprintf("... +%d more", x$n_nodes - 6)),
                         collapse = ", ")))

  # ---- Network structure ----
  mat <- x$weights
  n <- x$n_nodes
  max_possible <- if (x$directed) n * (n - 1) else n * (n - 1) / 2
  density <- if (max_possible > 0) x$n_edges / max_possible else 0
  cat(sprintf("  Edges: %d / %d (density: %.1f%%)\n",
              x$n_edges, as.integer(max_possible), density * 100))

  # Edge weight summary
  if (x$directed) {
    nz <- mat[mat != 0 & row(mat) != col(mat)]
  } else {
    nz <- mat[upper.tri(mat) & mat != 0]
  }

  if (length(nz) > 0) {
    is_assoc <- x$method %in% c("cor", "pcor", "glasso", "ising")
    if (is_assoc) {
      n_pos <- sum(nz > 0)
      n_neg <- sum(nz < 0)
      cat(sprintf("  Weights: [%.3f, %.3f]  |  +%d / -%d edges\n",
                  min(nz), max(nz), n_pos, n_neg))
    } else {
      cat(sprintf("  Weights: [%.3f, %.3f]  |  mean: %.3f\n",
                  min(nz), max(nz), mean(nz)))
    }

    # Top edges
    if (x$directed) {
      idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)
    } else {
      idx <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
    }
    if (nrow(idx) > 0) {
      w <- mat[idx]
      top_k <- min(5L, nrow(idx))
      ord <- order(abs(w), decreasing = TRUE)[seq_len(top_k)]
      cat("  Strongest edges:\n")
      for (j in ord) {
        arrow <- if (x$directed) " -> " else " -- "
        cat(sprintf("    %s%s%s  %.3f\n",
                    labels[idx[j, 1]], arrow, labels[idx[j, 2]], w[j]))
      }
    }
  }

  # Self-loops
  diag_vals <- diag(mat)
  n_self <- sum(diag_vals != 0)
  if (n_self > 0) {
    cat(sprintf("  Self-loops: %d  |  range: [%.3f, %.3f]\n",
                n_self, min(diag_vals[diag_vals != 0]),
                max(diag_vals[diag_vals != 0])))
  }

  # ---- Method-specific ----
  if (x$method == "glasso" && !is.null(x$gamma)) {
    cat(sprintf("  Gamma: %.2f  |  Lambda: %.4f\n",
                x$gamma, x$lambda_selected))
  }
  if (x$method == "ising") {
    cat(sprintf("  Gamma: %.2f  |  Rule: %s\n", x$gamma, x$rule))
    if (!is.null(x$thresholds)) {
      thr_rng <- range(x$thresholds)
      cat(sprintf("  Thresholds: [%.3f, %.3f]\n", thr_rng[1], thr_rng[2]))
    }
  }
  if (!is.null(x$scaling)) {
    cat(sprintf("  Scaling: %s\n", paste(x$scaling, collapse = " -> ")))
  }
  if (x$threshold > 0) cat(sprintf("  Threshold: %g\n", x$threshold))

  invisible(x)
}


#' Print Method for Group Network Object
#'
#' @param x A \code{netobject_group}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.netobject_group <- function(x, ...) {
  grps <- names(x)
  cat(sprintf("Group Networks (%d groups)\n", length(grps)))
  for (g in grps) {
    net <- x[[g]]
    cat(sprintf("  %s: %d nodes, %d edges\n", g, net$n_nodes, net$n_edges))
  }
  invisible(x)
}


#' Print Method for Multilevel Network Object
#'
#' @param x A \code{netobject_ml}.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.netobject_ml <- function(x, ...) {
  cat(sprintf("Multilevel Network (method: %s)\n", x$method))
  cat("-- Between-person --\n")
  b <- x$between
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", b$n_nodes, b$n_edges))
  if (!is.null(b$n)) {
    cat(sprintf("  Sample size: %d (unique persons)\n", b$n))
  }
  cat("-- Within-person --\n")
  w <- x$within
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", w$n_nodes, w$n_edges))
  if (!is.null(w$n)) {
    cat(sprintf("  Sample size: %d (observations)\n", w$n))
  }
  invisible(x)
}


#' Generate distinct pastel node colors
#' @noRd
.node_colors <- function(p) {
  palette <- c(
    "#A8D8EA", # light blue
    "#FFCAB1", # peach
    "#B5EAD7", # mint
    "#E2B6CF", # mauve
    "#FFDAC1", # apricot
    "#C7CEEA", # lavender
    "#F3E8C0", # cream
    "#D4F0C0", # pistachio
    "#F5C6D0", # pink
    "#B8E0D2", # seafoam
    "#EAC8A0", # sand
    "#C8B8DB", # lilac
    "#A0D2DB", # teal
    "#F0D9A0", # gold
    "#D8A8C8"  # orchid
  )
  rep_len(palette, p)
}


#' Plot Method for Network Object
#'
#' @description
#' Plots the network using \code{cograph::splot()}.
#' Requires the \pkg{cograph} package to be installed.
#' For association methods (\code{"glasso"}, \code{"pcor"}, \code{"cor"}),
#' node predictability (R\eqn{^2}) is shown as pie charts by default.
#'
#' @param x A \code{netobject}.
#' @param predictability Logical. If \code{TRUE}, display node predictability
#'   as pie charts for association methods (default: \code{TRUE}).
#' @param pie_color Character. Color for the predictability pie segments.
#'   Default: \code{"#377EB8"} (blue, following mgm convention).
#' @param ... Additional arguments passed to \code{cograph::splot()}.
#'
#' @export
plot.netobject <- function(x, predictability = TRUE,
                           pie_color = "#377EB8", ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }

  node_cols <- .node_colors(x$n_nodes)

  dots <- list(
    x = x$weights,
    directed = x$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  if (predictability && x$method %in% c("glasso", "pcor", "cor")) {
    r2 <- predictability.netobject(x)
    dots$pie_values <- r2
    dots$pie_colors <- pie_color
  }

  do.call(cograph::splot, dots)
}


#' Plot Method for Multilevel Network Object
#'
#' @description
#' Plots the between-person and within-person networks side by side using
#' \code{cograph::splot()}.
#' Requires the \pkg{cograph} package to be installed.
#'
#' @param x A \code{netobject_ml}.
#' @param predictability Logical. If \code{TRUE}, display node predictability
#'   as pie charts for association methods (default: \code{TRUE}).
#' @param pie_color Character. Color for the predictability pie segments.
#'   Default: \code{"#377EB8"}.
#' @param ... Additional arguments passed to \code{cograph::splot()}.
#'
#' @importFrom graphics par
#' @export
plot.netobject_ml <- function(x, predictability = TRUE,
                              pie_color = "#377EB8", ...) {
  if (!requireNamespace("cograph", quietly = TRUE)) {
    stop(
      "Package 'cograph' is required for plotting. ",
      "Install it with: install.packages('cograph')"
    )
  }
  old_par <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(old_par))

  p <- x$between$n_nodes
  node_cols <- .node_colors(p)

  dots <- list(
    directed = x$between$directed,
    node_fill = node_cols,
    edge_labels = TRUE,
    edge_label_size = 0.65,
    node_size = 8,
    theme = "colorblind",
    ...
  )

  if (predictability && x$method %in% c("glasso", "pcor", "cor")) {
    r2 <- predictability.netobject_ml(x)

    dots_b <- c(list(x = x$between$weights,
                     title = "Between-person",
                     pie_values = r2$between,
                     pie_colors = pie_color),
                dots)
    dots_w <- c(list(x = x$within$weights,
                     title = "Within-person",
                     pie_values = r2$within,
                     pie_colors = pie_color),
                dots)
  } else {
    dots_b <- c(list(x = x$between$weights,
                     title = "Between-person"), dots)
    dots_w <- c(list(x = x$within$weights,
                     title = "Within-person"), dots)
  }

  do.call(cograph::splot, dots_b)
  do.call(cograph::splot, dots_w)
}


# ---- Predictability ----

#' Compute Node Predictability
#'
#' @description
#' Computes the proportion of variance explained (R\eqn{^2}) for each node in
#' the network, following Haslbeck & Waldorp (2018).
#'
#' For \code{method = "glasso"} or \code{"pcor"}, predictability is computed
#' analytically from the precision matrix:
#' \deqn{R^2_j = 1 - 1 / \Omega_{jj}}
#' where \eqn{\Omega} is the precision (inverse correlation) matrix.
#'
#' For \code{method = "cor"}, predictability is the multiple R\eqn{^2} from
#' regressing each node on its network neighbors (nodes with non-zero edges).
#'
#' @param object A \code{netobject} or \code{netobject_ml} object.
#' @param ... Additional arguments (ignored).
#'
#' @return For \code{netobject}: a named numeric vector of R\eqn{^2} values
#'   (one per node, between 0 and 1).
#'
#'   For \code{netobject_ml}: a list with elements \code{$between} and
#'   \code{$within}, each a named numeric vector.
#'
#' @references
#' Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
#' predict observations? On the importance of predictability in network models.
#' \emph{Behavior Research Methods}, 50(2), 853--861.
#' \doi{10.3758/s13428-017-0910-x}
#'
#' @examples
#' \dontrun{
#' freq <- convert_sequence_format(group_regulation, format = "frequency")
#' net <- build_network(freq, method = "glasso")
#' predictability(net)
#'
#' # Plot with predictability rings
#' plot(net, predictability = TRUE)
#' }
#'
#' @export
predictability <- function(object, ...) {
  UseMethod("predictability")
}


#' @rdname predictability
#' @export
predictability.netobject <- function(object, ...) {
  if (object$method %in% c("glasso", "pcor")) {
    # From precision matrix: R²_j = 1 - 1/Omega_jj
    omega_diag <- diag(object$precision_matrix)
    r2 <- 1 - 1 / omega_diag
  } else {
    # cor method: multiple R² from correlation matrix
    S <- object$cor_matrix
    net <- object$weights
    p <- ncol(net)
    r2 <- vapply(seq_len(p), function(j) {
      neighbors <- which(net[j, ] != 0)
      if (length(neighbors) == 0L) return(0)
      if (length(neighbors) == 1L) return(S[neighbors, j]^2)
      r_vec <- S[neighbors, j]
      R_nn <- S[neighbors, neighbors]
      tryCatch(
        as.numeric(crossprod(r_vec, solve(R_nn, r_vec))),
        error = function(e) 0
      )
    }, numeric(1))
  }
  r2 <- pmin(pmax(r2, 0), 1)
  names(r2) <- object$nodes$label
  r2
}


#' @rdname predictability
#' @export
predictability.netobject_ml <- function(object, ...) {
  list(
    between = predictability(object$between),
    within  = predictability(object$within)
  )
}
