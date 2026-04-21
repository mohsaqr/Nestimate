# Python reference helpers for equivalence tests.
#
# Gated behind NESTIMATE_EQUIV_TESTS=true via skip_equiv_tests() in the caller.
# Only loaded when reticulate + pathpy are both available.
#
# Reticulate's auto-conversion stumbles on pathpy's tuple-keyed dicts and
# nested numpy arrays, so every extraction routes through a Python helper
# function defined in .ensure_py_helpers() and only then crosses back to R.

.PATHPY_HELPERS_LOADED <- new.env(parent = emptyenv())
.PATHPY_HELPERS_LOADED$ok <- FALSE

.ensure_py_helpers <- function() {
  if (isTRUE(.PATHPY_HELPERS_LOADED$ok)) return(invisible(NULL))
  reticulate::py_run_string("
import pathpy as pp
import numpy as np

def _paths_from_sequences(seqs):
    paths = pp.Paths()
    for s in seqs:
        if len(s) >= 2:
            paths.add_path(','.join(list(s)), frequency=1)
    return paths

def hon_edges_df(paths, order):
    hon = pp.HigherOrderNetwork(paths, k=int(order))
    rows = []
    for (u, v), attrs in hon.edges.items():
        w = attrs.get('weight', np.zeros(2))
        # w is [subpath_count, longest_path_count]; Nestimate counts every occurrence
        rows.append((str(u), str(v), float(w[0]) + float(w[1])))
    return rows

def mogen_cumulative_likelihood(paths, max_order):
    # Total log-likelihood of MultiOrderModel with given max_order.
    # Matches Nestimate's $log_likelihood[[k+1]] at k = max_order,
    # which uses hierarchical decomposition across orders 0..max_order.
    mog = pp.MultiOrderModel(paths, max_order=int(max_order))
    return float(mog.likelihood(paths, log=True))

def mogen_hon_transition_matrix(paths, order):
    # Return (nodes, count_matrix) for a HigherOrderNetwork at given order.
    # Subpath + longest-path counts summed per cell. Node labels are pathpy's
    # comma-joined k-tuples.
    hon = pp.HigherOrderNetwork(paths, k=int(order))
    nodes = sorted(hon.nodes)
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    M = np.zeros((n, n), dtype=float)
    for (u, v), attrs in hon.edges.items():
        w = attrs.get('weight', np.zeros(2))
        M[idx[u], idx[v]] = float(w[0]) + float(w[1])
    return {'nodes': nodes, 'counts': M.tolist()}

def mogen_stats(paths, max_order):
    mog = pp.MultiOrderModel(paths, max_order=int(max_order))
    # per-layer likelihoods, dof, and optimal order via LR test.
    # layer_likelihood(log=True) returns (log_lik, weight) tuple in pathpy 2.2.
    layers = []
    for k in range(int(max_order) + 1):
        lk_tuple = mog.layer_likelihood(paths, l=k, log=True)
        lk = float(lk_tuple[0]) if isinstance(lk_tuple, tuple) else float(lk_tuple)
        dof = mog.degrees_of_freedom(k)
        layers.append((int(k), lk, int(dof)))
    try:
        opt = int(mog.estimate_order(paths))
    except Exception:
        opt = -1
    return {'layers': layers, 'optimal_order': opt,
            'total_likelihood': float(mog.likelihood(paths, log=True)),
            'total_dof': int(mog.degrees_of_freedom(mog.max_order))}
")
  .PATHPY_HELPERS_LOADED$ok <- TRUE
  invisible(NULL)
}

#' Skip the calling test unless reticulate + a named Python module are loadable.
#' @param py_package Character. pip/import name of the reference library (e.g. "pathpy").
#' @noRd
skip_if_no_python_ref <- function(py_package = "pathpy") {
  skip_if_pkg_broken("reticulate")
  venv <- Sys.getenv("RETICULATE_PYTHON_ENV",
                     unset = "~/.virtualenvs/r-reticulate")
  if (dir.exists(path.expand(venv))) {
    try(reticulate::use_virtualenv(venv, required = TRUE), silent = TRUE)
  }
  if (!isTRUE(reticulate::py_available(initialize = TRUE))) {
    skip("Python runtime not available for reticulate")
  }
  if (!isTRUE(reticulate::py_module_available(py_package))) {
    skip(sprintf("Python module %s not installed (pip install %s)",
                 py_package, py_package))
  }
  .ensure_py_helpers()
}

#' Seed numpy to match R's seed.
#' @noRd
py_seed <- function(seed) {
  np <- reticulate::import("numpy", convert = FALSE)
  np$random$seed(as.integer(seed))
  invisible(NULL)
}

#' Build a pathpy.Paths from a list of character sequences. Comma in state names
#' would clash with pathpy's separator, so we guard against it.
#' @noRd
py_paths_from_sequences <- function(seqs) {
  if (any(grepl(",", unlist(seqs), fixed = TRUE))) {
    stop("state names contain comma; pathpy.Paths separator clash", call. = FALSE)
  }
  reticulate::py$`_paths_from_sequences`(seqs)
}

#' Pull HON edge list (columns: from, to, count) from a pathpy.Paths.
#' `count` sums subpath and longest-path weights so every k-gram occurrence is
#' counted — matching Nestimate semantics.
#' @noRd
py_hon_edges <- function(paths, order) {
  rows <- reticulate::py$hon_edges_df(paths, as.integer(order))
  if (length(rows) == 0L) {
    return(data.frame(from = character(0), to = character(0),
                      count = numeric(0), stringsAsFactors = FALSE))
  }
  from <- vapply(rows, `[[`, character(1), 1L)
  to <- vapply(rows, `[[`, character(1), 2L)
  count <- vapply(rows, `[[`, numeric(1), 3L)
  data.frame(from = from, to = to, count = count, stringsAsFactors = FALSE)
}

#' Pull MOGen per-layer stats from a pathpy.Paths.
#' Returns list(layers = data.frame(k, log_lik, dof),
#'              optimal_order, total_likelihood, total_dof).
#' @noRd
py_mogen_stats <- function(paths, max_order) {
  res <- reticulate::py$mogen_stats(paths, as.integer(max_order))
  layers <- res$layers
  df <- data.frame(
    k = vapply(layers, `[[`, integer(1), 1L),
    log_lik = vapply(layers, `[[`, numeric(1), 2L),
    dof = vapply(layers, `[[`, integer(1), 3L),
    stringsAsFactors = FALSE
  )
  list(layers = df, optimal_order = res$optimal_order,
       total_likelihood = res$total_likelihood, total_dof = res$total_dof)
}

#' Total log-likelihood of pathpy's MultiOrderModel at a given max_order.
#' @noRd
py_mogen_likelihood <- function(paths, max_order) {
  as.numeric(reticulate::py$mogen_cumulative_likelihood(paths, as.integer(max_order)))
}

#' Pull per-order count matrix from pathpy's HigherOrderNetwork.
#' Returns list(nodes = character, counts = numeric matrix). pathpy uses "," as
#' the node-tuple separator; Nestimate uses \x01. The caller must translate.
#' @noRd
py_hon_count_matrix <- function(paths, order) {
  res <- reticulate::py$mogen_hon_transition_matrix(paths, as.integer(order))
  mat <- do.call(rbind, lapply(res$counts, unlist))
  rownames(mat) <- unlist(res$nodes)
  colnames(mat) <- unlist(res$nodes)
  mat
}
