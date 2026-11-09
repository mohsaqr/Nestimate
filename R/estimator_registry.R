# ---- Global Estimator Registry ----

#' @noRd
.estimator_registry <- new.env(parent = emptyenv())


#' Register a Network Estimator
#'
#' @description
#' `register_estimator()` adds an estimator function to the registry that
#' [build_network()] reads, so that `build_network(data, method = name)` calls
#' it the way it calls a built-in method; the package registers its fourteen
#' built-in methods with this function when it loads. A name already in the
#' registry is replaced, built-ins included, without a message.
#' [list_estimators()] shows what is
#' registered and [remove_estimator()] takes
#' an entry out.
#'
#' A name that [build_network()] maps to a different built-in method (an
#' alias such as `"tna"` or `"ebicglasso"`) is resolved to that method before
#' the registry is read, so an estimator registered under it is never called;
#' pick a name that is not an alias.
#'
#' @details
#' [build_network()] calls the function as `fn(data = data, ...)`, where `...`
#' holds the `params` list, the further arguments given to [build_network()]
#' and a `format` setting, which is why `fn` must accept `...`. When `actor`
#' names a column of wide or association data (long data is prepared by
#' [build_network()] itself) and `fn` has an `id_col` argument (or, failing
#' that and for data that is not one-hot, an `id` argument), the actor
#' column is passed under that name unless `params` already sets it. Scaling,
#' thresholding and the edge list are then applied to the returned `matrix`,
#' as for a built-in method.
#'
#' Besides `matrix`, `nodes` and `directed`, the returned list may hold
#' `cleaned_data`, the data frame the network was estimated from, which
#' [build_network()] splits into the network's `data` and `metadata` fields;
#' any other named element is copied onto the `netobject`, replacing a field
#' of the same name. A return value without `matrix`, `nodes` or `directed`
#' makes [build_network()] stop with an error.
#'
#' @param name The name the estimator is called by, as `method` in
#'   [build_network()]: a single non-empty string. A name already registered
#'   is replaced.
#' @param fn The estimator function. Its first argument is `data`, and it
#'   must accept `...`, because [build_network()] passes every setting by
#'   name. It returns a list with `matrix` (the square numeric weight matrix),
#'   `nodes` (the node names, in the matrix's order) and `directed` (`TRUE`
#'   or `FALSE`); see *Details* for the optional elements.
#' @param description A short description of the estimator, a single string,
#'   shown by [list_estimators()].
#' @param directed Whether the estimator produces directed networks, `TRUE`
#'   or `FALSE`. [build_network()] reads it to decide whether `level` (the
#'   multilevel decomposition) is allowed: a directed estimator refuses it.
#'   The edges of the network follow the `directed` element `fn` returns,
#'   but [bootstrap_network()] reads this registered flag to build its edge
#'   tables, so the two should agree.
#'
#' @return `NULL`, invisibly. The estimator is added to the registry as a
#'   side effect.
#'
#' @examples
#' my_fn <- function(data, ...) {
#'   m <- cor(data)
#'   diag(m) <- 0
#'   list(matrix = m, nodes = colnames(m), directed = FALSE)
#' }
#' register_estimator("my_cor", my_fn, "Custom correlation", directed = FALSE)
#' df <- data.frame(A = rnorm(20), B = rnorm(20), C = rnorm(20))
#' net <- build_network(df, method = "my_cor")
#' remove_estimator("my_cor")
#'
#' @seealso \code{\link{get_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{remove_estimator}}, \code{\link{build_network}}
#'
#' @export
register_estimator <- function(name, fn, description, directed) {
  .assert(
    is.character(name), length(name) == 1, nzchar(name),
    is.function(fn),
    is.character(description), length(description) == 1,
    is.logical(directed), length(directed) == 1
  )
  entry <- list(fn = fn, description = description, directed = directed)
  assign(name, entry, envir = .estimator_registry)
  invisible(NULL)
}


#' Retrieve a Registered Estimator
#'
#' @description
#' `get_estimator()` returns the registry entry of one estimator: its
#' function, description and directedness, as [build_network()] uses them. A
#' name that is not registered raises an error that lists the registered
#' names. Method aliases such as `"tna"` are not resolved here; give the
#' registered name, as [list_estimators()]
#' shows it.
#'
#' @param name The registered name of the estimator, a single string.
#'
#' @return A list with `fn` (the estimator function), `description` (a
#'   string) and `directed` (`TRUE` or `FALSE`).
#'
#' @examples
#' est <- get_estimator("relative")
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}}
#'
#' @export
get_estimator <- function(name) {
  .assert(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    available <- ls(envir = .estimator_registry)
    stop(errorCondition(
      sprintf("Estimator '%s' not found. Available: %s",
              name, paste(sort(available), collapse = ", ")),
      class = "nestimate_bad_input", call = NULL
    ))
  }
  get(name, envir = .estimator_registry, inherits = FALSE)
}


#' List All Registered Estimators
#'
#' @description
#' `list_estimators()` returns the estimators [build_network()] can call by
#' name: the fourteen built-ins and any added with
#' [register_estimator()]. Aliases such as
#' `"tna"` or `"ebicglasso"` are not listed; the *Estimators* section of
#' [build_network()] gives them.
#'
#' @return A data frame with one row per registered estimator, sorted by
#'   name, and columns `name`, `description` and `directed` (`TRUE` for
#'   directed networks). With an empty registry it has no rows.
#'
#' @examples
#' list_estimators()
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{get_estimator}}
#'
#' @export
list_estimators <- function() {
  nms <- ls(envir = .estimator_registry)
  if (length(nms) == 0) {
    return(data.frame(
      name = character(0),
      description = character(0),
      directed = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  nms <- sort(nms)
  entries <- lapply(nms, function(nm) {
    e <- get(nm, envir = .estimator_registry, inherits = FALSE)
    data.frame(
      name = nm,
      description = e$description,
      directed = e$directed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, entries)
}


#' Remove a Registered Estimator
#'
#' @description
#' `remove_estimator()` takes an estimator out of the registry, so
#' [build_network()] no longer finds it by name. Built-ins can be removed
#' too; they come back when the package is loaded again. A name that is not
#' registered raises an error.
#'
#' @param name The registered name of the estimator to remove, a single
#'   string.
#'
#' @return `NULL`, invisibly. The entry is removed from the registry as a
#'   side effect.
#'
#' @examples
#' register_estimator("test_est", function(data, ...) diag(3),
#'   description = "test", directed = FALSE)
#' remove_estimator("test_est")
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}}
#'
#' @export
remove_estimator <- function(name) {
  .assert(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    stop(errorCondition(
      sprintf("Estimator '%s' not found in registry.", name),
      class = "nestimate_bad_input", call = NULL
    ))
  }
  rm(list = name, envir = .estimator_registry)
  invisible(NULL)
}


#' Register Built-in Estimators
#' @noRd
.register_builtin_estimators <- function() {
  register_estimator("relative", .estimator_relative,
                     "Row-normalized transition probabilities", directed = TRUE)
  register_estimator("frequency", .estimator_frequency,
                     "Raw transition frequency counts", directed = TRUE)
  register_estimator("co_occurrence", .estimator_co_occurrence,
                     "Co-occurrence within sequences", directed = FALSE)
  register_estimator("cor", .estimator_cor,
                     "Pairwise correlation network", directed = FALSE)
  register_estimator("pcor", .estimator_pcor,
                     "Unregularized partial correlations", directed = FALSE)
  register_estimator("glasso", .estimator_glasso,
                     "EBICglasso regularized partial correlations",
                     directed = FALSE)
  register_estimator("ising", .estimator_ising,
                     "Ising model (L1-penalized logistic regression)",
                     directed = FALSE)
  register_estimator("mgm", .estimator_mgm,
                     "Mixed Graphical Model (nodewise lasso, EBIC, LW threshold)",
                     directed = FALSE)
  register_estimator("attention", .estimator_attention,
                     "Decay-weighted attention transitions", directed = TRUE)
  register_estimator("ngram", .estimator_ngram,
                     "n-gram transitions (adjacent pairs per n-gram window)",
                     directed = TRUE)
  register_estimator("gap", .estimator_gap,
                     "Gap-allowed transitions weighted by 1/distance",
                     directed = TRUE)
  register_estimator("reverse", .estimator_reverse,
                     "Reverse (reply) transitions: transpose of frequency",
                     directed = TRUE)
  register_estimator("wtna", .estimator_wtna,
                     "Window-based TNA transitions (one-hot)", directed = TRUE)
  register_estimator("wtna_cooccurrence",
                     function(data, codes = NULL, window_size = 3L,
                              mode = "non-overlapping", actor = NULL,
                              type = "frequency", ...) {
                       .estimator_wtna_core(data, codes = codes,
                         window_size = window_size, mode = mode,
                         actor = actor, wtna_method = "cooccurrence",
                         type = type, ...)
                     },
                     "Window-based TNA co-occurrence (one-hot)", directed = FALSE)
}
