#' Estimate a network (deprecated)
#'
#' @description
#' `estimate_network()` is deprecated: use [build_network()] instead. It
#' passes its arguments to [build_network()] unchanged and returns what that
#' returns, after a warning of class `deprecatedWarning`.
#'
#' @param data The data, as for [build_network()].
#' @param params Settings for the estimator, as a named list, as on the
#'   method's page (for example [build_glasso()]). Default `list()`.
#' @param scaling,threshold As in [build_network()].
#' @param level Multilevel decomposition for the association methods, as in
#'   [build_cor()]. Default `NULL`.
#' @param method The estimator, as in [build_network()]. Default
#'   `"relative"`, kept for backward compatibility; [build_network()] has no
#'   default.
#' @param ... Further arguments for [build_network()], such as `actor`,
#'   `action`, `time` or `group`.
#'
#' @return A `netobject`, exactly as [build_network()] returns it (a
#'   `netobject_group` or `netobject_ml` when `...` or `level` asks for one).
#'
#' @examples
#' data <- data.frame(T1 = c("x","y","z","x"), T2 = c("y","x","z","y"))
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


# ---- Convenience wrappers for build_network ----

# The convenience wrappers fix `method`; say so, instead of R's
# "formal argument matched by multiple actual arguments".
.refuse_method_arg <- function(wrapper, method, dot_names) {
  if ("method" %in% dot_names) {
    stop(errorCondition(sprintf(paste0(
      "%s() always uses method = \"%s\"; for another method call ",
      "build_network(method = ...)."), wrapper, method),
      class = "nestimate_bad_input", call = NULL))
  }
}

# Calls build_network() from a wrapper with only the arguments the caller set,
# each passed by name and evaluated in the wrapper's frame, plus the wrapper's
# `...` in its place, so unnamed extras bind to build_network()'s formals as
# they always did. build_network() treats a shaping argument the caller set
# (`scaling`, `threshold`, `params`, `start`, `end`) as an override when it
# rebuilds the networks of a clustering, so a wrapper must not pass its own
# defaults on. "Set" is read from the matched call, as build_network() reads
# it: an argument forwarded from a caller that left it missing is still set,
# and evaluating it raises the same error a direct call raises.
.forward_to_build_network <- function(wrapper, method, dot_names,
                                      rename = character()) {
  wrapper_env <- parent.frame()
  .refuse_method_arg(wrapper, method, dot_names)
  wrapper_call <- match.call(sys.function(sys.parent()),
                             sys.call(sys.parent()), expand.dots = FALSE,
                             envir = parent.frame(2L))
  given <- setdiff(names(wrapper_call)[-1L], c("data", "...", ""))
  # an alias's own name for a build_network() argument (build_latency()'s
  # `floor` is `latency_floor`)
  passed_as <- given
  renamed <- passed_as %in% names(rename)
  passed_as[renamed] <- rename[passed_as[renamed]]
  forwarded <- c(list(quote(build_network), data = quote(data),
                      method = method),
                 stats::setNames(lapply(given, as.name), passed_as),
                 list(quote(...)))
  eval(as.call(forwarded), wrapper_env)
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
    atna              = "attention",
    mixed_graphical   = "mgm",
    mixed             = "mgm",
    `co-occurrence`   = "co_occurrence",
    `n-gram`          = "ngram",
    n_gram            = "ngram"
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
.scaling_includes_zeros <- function(method) {
  method %in% c("relative", "frequency", "co_occurrence", "attention", "ngram",
                "gap", "reverse", "latency")
}

.apply_scaling <- function(mat, scaling, include_zeros = FALSE) {
  # Loop kept: the scalings are applied in the order given, each to the
  # result of the one before.
  for (s in scaling) {
    mat <- switch(s,
      minmax = {
        vals <- if (isTRUE(include_zeros)) as.vector(mat) else mat[mat != 0]
        if (length(vals) == 0) {
          mat
        } else {
          rng <- range(vals, na.rm = TRUE)
          if (rng[1] == rng[2]) mat
          else {
            if (isTRUE(include_zeros)) {
              mat[] <- (mat - rng[1]) / (rng[2] - rng[1])
            } else {
              # Structure-preserving min-max for association networks: the
              # smallest genuine edge must not collide with the structural-
              # zero value (which would drop a real edge and corrupt
              # n_edges). Map the non-zero block into the strictly positive
              # interval (eps, 1] so the block minimum becomes eps, not 0.
              eps <- sqrt(.Machine$double.eps)
              nz <- mat != 0
              mat[nz] <- eps + (1 - eps) *
                (mat[nz] - rng[1]) / (rng[2] - rng[1])
            }
            mat
          }
        }
      },
      max = {
        max_abs <- max(abs(mat))
        if (max_abs > 0) mat / max_abs else mat # nocov
      },
      rank = {
        # Rank scaling always preserves structural zeros: a non-edge must
        # remain a non-edge after scaling. include_zeros = TRUE only affects
        # min-max / range scalings; for ranks it would promote absent
        # transitions to positive weights and densify the network.
        nz <- mat != 0
        if (any(nz)) {
          mat[nz] <- rank(mat[nz], ties.method = "average")
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
    # Keep self-loops too: every non-zero entry is a real edge.
    idx <- which(mat != 0, arr.ind = TRUE)
  } else {
    # row <= col keeps the upper triangle PLUS the diagonal (one row
    # per self-loop, no double-count for undirected networks).
    idx <- which(mat != 0 & row(mat) <= col(mat), arr.ind = TRUE)
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
  .assert("`data` must be a data frame" = is.data.frame(data)) # nocov
  grp_var <- id_col[1]

  if (!grp_var %in% names(data)) {
    stop(errorCondition(paste0("id_col '", grp_var, "' not found in data."),
                        class = "nestimate_bad_input", call = NULL))
  }

  # Get numeric columns (exclude id columns and "rid")
  exclude <- c(id_col, "rid")
  numeric_cols <- vapply(data, is.numeric, logical(1))
  keep <- setdiff(names(data)[numeric_cols], exclude)

  if (length(keep) < 2) {
    stop(errorCondition(
      "At least 2 numeric columns are required for multilevel decomposition.",
      class = "nestimate_bad_input", call = NULL))
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
      stop(errorCondition(paste0(
        "Fewer than 3 rows remain after dropping ",
        "single-observation persons."),
        class = "nestimate_bad_input", call = NULL))
    }

    # Person-mean center each variable
    mat_m <- as.matrix(mat)
    mat_m[] <- vapply(seq_len(ncol(mat_m)), function(j) {
      mat_m[, j] - ave(mat_m[, j], id_vals, FUN = mean)
    }, numeric(nrow(mat_m)))

    return(as.data.frame(mat_m))
  }

  data
}
