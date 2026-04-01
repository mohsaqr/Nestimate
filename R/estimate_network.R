#' Estimate a Network (Deprecated)
#'
#' @description
#' This function is deprecated. Use \code{\link{build_network}} instead.
#'
#' @inheritParams build_network
#' @param method Character. Defaults to \code{"relative"} for backward
#'   compatibility.
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#'
#' @return A \code{netobject} (see \code{\link{build_network}}).
#'
#' @examples
#' data <- data.frame(A = c("x","y","z","x"), B = c("y","x","z","y"))
#' net <- estimate_network(data, method = "relative")
#'
#' @seealso \code{\link{build_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @export
estimate_network <- function(data,
                             method = "relative",
                             params = list(),
                             scaling = NULL,
                             threshold = 0,
                             level = NULL,
                             ...) {
  .Deprecated("build_network")
  build_network(
    data = data,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    level = level,
    ...
  )
}


# ---- Shared internal helpers ----
# These are used by build_network(), bootstrap_network(), and other functions.


# ---- Method alias resolution ----

#' Resolve method aliases to canonical names
#' @noRd
.resolve_method_alias <- function(method) {
  aliases <- c(
    ebicglasso        = "glasso",
    regularized       = "glasso",
    partial           = "pcor",
    correlation       = "cor",
    corr              = "cor",
    transition        = "relative",
    tna               = "relative",
    counts            = "frequency",
    ftna              = "frequency",
    cna               = "co_occurrence",
    wcna              = "co_occurrence",
    wtna_transition   = "wtna",
    wtna_cooccurrence = "wtna_cooccurrence",
    isingfit          = "ising",
    atna              = "attention"
  )
  if (method %in% names(aliases)) {
    aliases[[method]]
  } else {
    method
  }
}


# ---- Post-estimation scaling ----

#' Apply scaling transformations to a network matrix
#' @noRd
.apply_scaling <- function(mat, scaling) {
  for (s in scaling) {
    mat <- switch(s,
      minmax = {
        vals <- mat[mat != 0]
        if (length(vals) == 0) {
          mat
        } else {
          rng <- range(vals)
          if (rng[1] == rng[2]) mat
          else {
            mat[mat != 0] <- (mat[mat != 0] - rng[1]) / (rng[2] - rng[1])
            mat
          }
        }
      },
      max = {
        max_abs <- max(abs(mat))
        if (max_abs > 0) mat / max_abs else mat # nocov
      },
      rank = {
        nz <- mat != 0
        if (any(nz)) {
          mat[nz] <- rank(mat[nz])
          mat
        } else {
          mat
        }
      },
      normalize = {
        rs <- rowSums(abs(mat))
        nonzero_rows <- rs > 0
        mat[nonzero_rows, ] <- mat[nonzero_rows, ] / rs[nonzero_rows]
        mat
      },
      mat  # default: no change
    )
  }
  mat
}


# ---- Edge extraction ----

#' Extract non-zero edges from a network matrix
#'
#' For undirected networks, uses upper triangle only.
#' For directed networks, uses all non-diagonal non-zero entries.
#'
#' @noRd
.extract_edges_from_matrix <- function(mat, directed = FALSE) {
  if (directed) {
    idx <- which(mat != 0 & row(mat) != col(mat), arr.ind = TRUE)
  } else {
    idx <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  }

  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from   = as.integer(idx[, 1]),
    to     = as.integer(idx[, 2]),
    weight = mat[idx],
    stringsAsFactors = FALSE
  )
}


# ---- Multilevel decomposition ----

#' Decompose data for multilevel analysis
#'
#' @param data Data frame.
#' @param id_col Character. Grouping variable name.
#' @param level Character: "between" or "within".
#'
#' @return Transformed data frame.
#' @noRd
.decompose_multilevel <- function(data, id_col, level) {
  stopifnot(is.data.frame(data)) # nocov
  grp_var <- id_col[1]

  if (!grp_var %in% names(data)) {
    stop("id_col '", grp_var, "' not found in data.", call. = FALSE)
  }

  # Get numeric columns (exclude id columns and "rid")
  exclude <- c(id_col, "rid")
  numeric_cols <- vapply(data, is.numeric, logical(1))
  keep <- setdiff(names(data)[numeric_cols], exclude)

  if (length(keep) < 2) {
    stop("At least 2 numeric columns are required for multilevel decomposition.")
  }

  mat <- data[, keep, drop = FALSE]
  id_vals <- data[[grp_var]]

  if (level == "between") {
    # Aggregate to person means
    mat$.id <- id_vals
    agg <- aggregate(. ~ .id, data = mat, FUN = mean)
    result <- agg[, names(agg) != ".id", drop = FALSE]
    return(as.data.frame(result))

  } else if (level == "within") {
    # Drop persons with < 2 observations
    tab <- table(id_vals)
    multi <- names(tab[tab >= 2])
    keep_rows <- id_vals %in% multi
    if (any(!keep_rows)) {
      n_single <- sum(!keep_rows)
      message("Dropping ", n_single,
              " single-observation rows (within-person centering).")
      mat <- mat[keep_rows, , drop = FALSE]
      id_vals <- id_vals[keep_rows]
    }

    if (nrow(mat) < 3) {
      stop("Fewer than 3 rows remain after dropping ",
           "single-observation persons.")
    }

    # Person-mean center each variable
    mat_m <- as.matrix(mat)
    for (j in seq_len(ncol(mat_m))) {
      mat_m[, j] <- mat_m[, j] - ave(mat_m[, j], id_vals, FUN = mean)
    }

    return(as.data.frame(mat_m))
  }

  data
}
