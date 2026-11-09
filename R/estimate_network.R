#' Estimate a Network (Deprecated)
#'
#' @description
#' `estimate_network()` is deprecated: use [build_network()] instead. It
#' passes its arguments to [build_network()] unchanged and returns what that
#' returns, after a warning of class `deprecatedWarning`.
#'
#' @inheritParams build_network
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

#' Build a Transition Network (TNA)
#'
#' @description
#' `build_tna()` estimates a transition network (TNA): the probability that
#' each state is followed by each state, itself included. It is
#' [build_network()] with `method = "relative"`: transition counts
#' normalised to sum to 1 in each row.
#'
#' An edge from A to B is the share of A's outgoing transitions that go to
#' B, so a state's edges compare within its row, not across rows. The
#' weights are point estimates: check how far they move under resampling
#' with [bootstrap_network()], and compare groups with [permutation()].
#'
#' @inheritParams build_network
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, an event log (one row per event, with `action`
#'   set), one-hot code columns (declared with `codes` or
#'   `format = "onehot"`), or a square matrix of transition counts or
#'   weights, which is row-normalised (every row needs a positive sum).
#' @param actor The column naming who each row belongs to (a person, a
#'   team). In an event log, each actor's events form a sequence, split
#'   further by `session` or `time_threshold`. Default `NULL`.
#' @param start A start state put before every sequence, so the network
#'   shows where sequences begin. `FALSE` (default) adds none, `TRUE` adds
#'   one labelled `"Start"`, a single string adds one with that label.
#' @param end An end state put after each sequence's last observed state,
#'   so the network shows where sequences stop; it has no outgoing edges.
#'   `FALSE` (default), `TRUE` (`"End"`) or a single string, as `start`.
#' @param ... Further arguments for [build_network()], such as `group`,
#'   `params`, `scaling`, `threshold` or `state_cols`. `method` cannot be
#'   given: it is always `"relative"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`), exactly as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the weights with
#'   [extract_transition_matrix()], and draw it with `cograph::splot()`.
#' @seealso \code{\link{build_network}}
#' @examples
#' net <- build_tna(trajectories)
#' net
#' @export
build_tna <- function(data, start = FALSE, end = FALSE,
                      actor = NULL, action = NULL, time = NULL,
                      session = NULL, ...) {
  .refuse_method_arg("build_tna", "relative", ...names())
  build_network(data, method = "relative", actor = actor, action = action,
                time = time, session = session, start = start, end = end,
                ...)
}

#' Build a Frequency Transition Network (FTNA)
#'
#' @description
#' `build_ftna()` estimates a frequency transition network (FTNA): how many
#' times each state is followed by each state, itself included, as raw
#' counts. It is
#' [build_network()] with `method = "frequency"`, and takes the same input
#' as [build_tna()]; a matrix must hold non-negative whole counts.
#'
#' Counts grow with the number and length of sequences, so they compare
#' within one network; for transition probabilities use [build_tna()]. Check
#' how far the counts move under resampling with [bootstrap_network()].
#'
#' @inheritParams build_network
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, an event log (one row per event, with `action`
#'   set), one-hot code columns, or a square matrix of transition
#'   counts (non-negative whole numbers). Column roles are declared as in
#'   [build_network()].
#' @param actor The column naming who each row belongs to (a person, a
#'   team). In an event log, each actor's events form a sequence, split
#'   further by `session` or `time_threshold`. Default `NULL`.
#' @param start A start state put before every sequence, so the network
#'   shows where sequences begin. `FALSE` (default) adds none, `TRUE` adds
#'   one labelled `"Start"`, a single string adds one with that label.
#' @param end An end state put after each sequence's last observed state,
#'   so the network shows where sequences stop; it has no outgoing edges.
#'   `FALSE` (default), `TRUE` (`"End"`) or a single string, as `start`.
#' @param ... Further arguments for [build_network()], such as `group`,
#'   `params`, `scaling`, `threshold` or `state_cols`. `method` cannot be
#'   given: it is always `"frequency"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`), exactly as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the weights with
#'   [extract_transition_matrix()], and draw it with `cograph::splot()`.
#' @seealso \code{\link{build_network}}
#' @examples
#' net <- build_ftna(trajectories)
#' net
#' @export
build_ftna <- function(data, start = FALSE, end = FALSE,
                       actor = NULL, action = NULL, time = NULL,
                       session = NULL, ...) {
  .refuse_method_arg("build_ftna", "frequency", ...names())
  build_network(data, method = "frequency", actor = actor, action = action,
                time = time, session = session, start = start, end = end,
                ...)
}

#' Build an Attention-Weighted Transition Network (ATNA)
#'
#' @description
#' `build_atna()` estimates an attention-weighted transition network (ATNA):
#' every later state in a sequence, not only the next one, adds to the edge
#' from an earlier state, with a weight that decays with the number of
#' positions between them, `exp(-distance / lambda)`. It is
#' [build_network()] with `method = "attention"`. The decay is set in
#' `params`: `lambda` (default `1`; larger values decay more slowly) and
#' `direction` (`"forward"`, the default, `"backward"` or `"both"`).
#' `duration` (one value per position, accumulated into times) or
#' `time_matrix` (one row per sequence, one column per position) measure
#' the distance in time instead of positions, and `decay`, a function of
#' `(ti, tj, lambda)` returning one weight per pair, replaces the decay
#' function.
#'
#' The weights are summed decay weights, not counts or probabilities, so they
#' compare within one network; add `scaling = "normalize"` for rows that sum
#' to 1. Check how far they move under resampling with
#' [bootstrap_network()].
#'
#' @inheritParams build_network
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, or an event log (one row per event, with `action`
#'   set). Column roles are declared as in [build_network()].
#'   One-hot data and matrices are refused.
#' @param actor The column naming who each row belongs to (a person, a
#'   team). In an event log, each actor's events form a sequence, split
#'   further by `session` or `time_threshold`. Default `NULL`.
#' @param start A start state put before every sequence, so the network
#'   shows where sequences begin. `FALSE` (default) adds none, `TRUE` adds
#'   one labelled `"Start"`, a single string adds one with that label.
#' @param end An end state put after each sequence's last observed state,
#'   so the network shows where sequences stop; it has no outgoing edges.
#'   `FALSE` (default), `TRUE` (`"End"`) or a single string, as `start`.
#' @param ... Further arguments for [build_network()], such as `group`,
#'   `params`, `scaling`, `threshold` or `state_cols`. `method` cannot be
#'   given: it is always `"attention"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`), exactly as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the weights with
#'   [extract_transition_matrix()], and draw it with `cograph::splot()`.
#' @seealso \code{\link{build_network}}
#' @examples
#' net <- build_atna(trajectories)
#' net
#' @export
build_atna <- function(data, start = FALSE, end = FALSE,
                       actor = NULL, action = NULL, time = NULL,
                       session = NULL, ...) {
  .refuse_method_arg("build_atna", "attention", ...names())
  build_network(data, method = "attention", actor = actor, action = action,
                time = time, session = session, start = start, end = end,
                ...)
}

#' Build a Co-occurrence Network (CNA)
#'
#' @description
#' `build_cna()` estimates an undirected co-occurrence network (CNA): for
#' every pair of positions in a sequence where both are observed, the two
#' states count as co-occurring once, in either order. It is
#' [build_network()] with `method = "co_occurrence"`, and takes sequences in
#' a wide table, an event log, one-hot codes (declared with `codes` or
#' `format = "onehot"`, and counted in windows as [wtna()] does) or a square
#' matrix of counts.
#'
#' An edge counts how often two states occur together within sequences
#' (over pairs of positions, so a sequence can add to one edge several
#' times), not which came first; for order, use [build_tna()]. Check how
#' far the counts move under resampling with [bootstrap_network()].
#'
#' @inheritParams build_network
#' @param data The sequences: a data frame with one row per sequence and one
#'   column per position, an event log (one row per event, with `action`
#'   set), one-hot code columns, or a square matrix of
#'   co-occurrence counts (non-negative whole numbers). Column roles are
#'   declared as in [build_network()].
#' @param actor The column naming who each row belongs to (a person, a
#'   team). In an event log, each actor's events form a sequence, split
#'   further by `session` or `time_threshold`. Default `NULL`.
#' @param start A start state put before every sequence, so the network
#'   shows where sequences begin. `FALSE` (default) adds none, `TRUE` adds
#'   one labelled `"Start"`, a single string adds one with that label.
#' @param end An end state put after each sequence's last observed state,
#'   so the network shows where sequences stop; it has no outgoing edges.
#'   `FALSE` (default), `TRUE` (`"End"`) or a single string, as `start`.
#' @param ... Further arguments for [build_network()], such as `group`,
#'   `params`, `scaling`, `threshold` or `state_cols`. `method` cannot be
#'   given: it is always `"co_occurrence"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`), exactly as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()] and the weights with
#'   [extract_transition_matrix()], and draw it with `cograph::splot()`.
#' @seealso \code{\link{build_network}}, \code{\link{cooccurrence}} for
#'   delimited-field, bipartite, and other non-sequence co-occurrence formats.
#' @examples
#' net <- build_cna(trajectories)
#' net
#' @export
build_cna <- function(data, start = FALSE, end = FALSE,
                      actor = NULL, action = NULL, time = NULL,
                      session = NULL, ...) {
  .refuse_method_arg("build_cna", "co_occurrence", ...names())
  build_network(data, method = "co_occurrence", actor = actor, action = action,
                time = time, session = session, start = start, end = end,
                ...)
}

#' Build a Correlation Network
#'
#' @description
#' `build_cor()` estimates a correlation network: an undirected graph of
#' variables in which each edge is the correlation of two variables. It is
#' [build_network()] with `method = "cor"`. The correlation is Pearson's by
#' default; `params = list(cor_method = "spearman")` (or `"kendall"`) passes
#' another method to [stats::cor()].
#'
#' A correlation edge is a marginal association: it can come from a shared
#' link with a third variable. To separate direct from indirect links, use
#' [build_pcor()] or [build_glasso()]. Check the stability of the edges with
#' [bootstrap_network()].
#'
#' @param data The variables: a data frame with one row per observation and
#'   one column per variable, or a square symmetric correlation or
#'   covariance matrix, which needs the sample size in `params` (as in
#'   `params = list(n = 250)`). Non-numeric columns are set aside; a frame of
#'   two or more character or factor columns and fewer than two numeric ones
#'   is read as sequences and turned into per-row state counts first. A
#'   column named `rid` is set aside; columns with non-syntactic names,
#'   all-missing and constant columns, and rows with a missing value are
#'   dropped, each with a message; fewer than 2 variables or 3 complete rows
#'   raise an error.
#' @param actor The id column naming whose observation each row is (a
#'   person). With `level`, the decomposition groups by it; without, it is
#'   set aside, not modelled. Default `NULL`.
#' @param level Multilevel decomposition, with `actor`: `"between"`
#'   estimates the network of person means, `"within"` the network of
#'   deviations from them, and `"both"` returns the two as a
#'   `netobject_ml`. Default `NULL` (no decomposition).
#' @param ... Further arguments for [build_network()], such as `params`
#'   (`cor_method`, `n`), `threshold`, `scaling` or `group`. `method` cannot
#'   be given: it is always `"cor"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), exactly as [build_network()] returns it: print
#'   it for an overview, read the edges with [extract_edges()], and draw it
#'   with `cograph::splot()`.
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_cor(srl_strategies)
#' @export
build_cor <- function(data, actor = NULL, level = NULL, ...) {
  .refuse_method_arg("build_cor", "cor", ...names())
  build_network(data, method = "cor", actor = actor, level = level, ...)
}

#' Build a Partial Correlation Network
#'
#' @description
#' `build_pcor()` estimates a partial correlation network: an undirected
#' graph of variables in which each edge is the correlation of two variables
#' after the others are held constant, read from the inverse of the
#' correlation matrix. It is [build_network()] with `method = "pcor"`, with
#' no regularisation.
#'
#' Unregularised partial correlations are noisy when there are many
#' variables for the number of observations, and they cannot be estimated
#' when the correlation matrix is singular (more variables than
#' observations, or collinear ones): that raises an error, and a
#' near-singular matrix a warning. [build_glasso()] estimates a sparse
#' network in those cases. Check the stability of the edges with
#' [bootstrap_network()].
#'
#' @param data The variables: a data frame with one row per observation and
#'   one column per variable, or a square symmetric correlation or
#'   covariance matrix, which needs the sample size in `params` (as in
#'   `params = list(n = 250)`). Non-numeric columns are set aside; a frame of
#'   two or more character or factor columns and fewer than two numeric ones
#'   is read as sequences and turned into per-row state counts first. A
#'   column named `rid` is set aside; columns with non-syntactic names,
#'   all-missing and constant columns, and rows with a missing value are
#'   dropped, each with a message; fewer than 2 variables or 3 complete rows
#'   raise an error.
#' @param actor The id column naming whose observation each row is (a
#'   person). With `level`, the decomposition groups by it; without, it is
#'   set aside, not modelled. Default `NULL`.
#' @param level Multilevel decomposition, with `actor`: `"between"`
#'   estimates the network of person means, `"within"` the network of
#'   deviations from them, and `"both"` returns the two as a
#'   `netobject_ml`. Default `NULL` (no decomposition).
#' @param ... Further arguments for [build_network()], such as `params`
#'   (`cor_method`, `n`), `threshold`, `scaling` or `group`. `method` cannot
#'   be given: it is always `"pcor"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), exactly as [build_network()] returns it: print
#'   it for an overview, read the edges with [extract_edges()], and draw it
#'   with `cograph::splot()`.
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_pcor(srl_strategies)
#' @export
build_pcor <- function(data, actor = NULL, level = NULL, ...) {
  .refuse_method_arg("build_pcor", "pcor", ...names())
  build_network(data, method = "pcor", actor = actor, level = level, ...)
}

#' Build a Graphical Lasso Network (EBICglasso)
#'
#' @description
#' `build_glasso()` estimates a regularised partial correlation network with
#' the graphical lasso: an L1 penalty shrinks small partial correlations to
#' exactly zero, and the penalty is chosen along a path of `nlambda` values
#' by the extended Bayesian information criterion (EBIC; Foygel & Drton,
#' 2010). It is [build_network()] with `method = "glasso"`. The settings go
#' in `params`: `gamma` (the EBIC hyperparameter, default `0.5`; larger
#' values give sparser networks), `nlambda` (default `100`),
#' `lambda.min.ratio` (default `0.01`), `penalize.diagonal` (default
#' `FALSE`), `refit` (default `FALSE`; `TRUE` re-estimates the selected
#' edges without the penalty) and `cor_method`.
#'
#' An absent edge means the penalty set it to zero at the selected lambda,
#' not that the two variables are shown to be independent. Check the
#' stability of the edges and centralities with [boot_glasso()].
#'
#' @param data The variables: a data frame with one row per observation and
#'   one column per variable, or a square symmetric correlation or
#'   covariance matrix, which needs the sample size in `params` (as in
#'   `params = list(n = 250)`). Non-numeric columns are set aside; a frame of
#'   two or more character or factor columns and fewer than two numeric ones
#'   is read as sequences and turned into per-row state counts first. A
#'   column named `rid` is set aside; columns with non-syntactic names,
#'   all-missing and constant columns, and rows with a missing value are
#'   dropped, each with a message; fewer than 2 variables or 3 complete rows
#'   raise an error.
#' @param actor The id column naming whose observation each row is (a
#'   person). With `level`, the decomposition groups by it; without, it is
#'   set aside, not modelled. Default `NULL`.
#' @param level Multilevel decomposition, with `actor`: `"between"`
#'   estimates the network of person means, `"within"` the network of
#'   deviations from them, and `"both"` returns the two as a
#'   `netobject_ml`. Default `NULL` (no decomposition).
#' @param ... Further arguments for [build_network()], such as `params`
#'   (`gamma`, `nlambda`, `refit`, `n` and the others above), `threshold`,
#'   `scaling` or `group`. `method` cannot be given: it is always
#'   `"glasso"` (error of class `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`, a `netobject_ml`
#'   with `level = "both"`), exactly as [build_network()] returns it: print
#'   it for an overview, read the edges with [extract_edges()], and draw it
#'   with `cograph::splot()`.
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse
#' covariance estimation with the graphical lasso. \emph{Biostatistics},
#' 9(3), 432--441. \doi{10.1093/biostatistics/kxm045}
#'
#' Foygel, R., & Drton, M. (2010). Extended Bayesian information criteria
#' for Gaussian graphical models. \emph{Advances in Neural Information
#' Processing Systems}, 23, 604--612.
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_glasso(srl_strategies)
#' @export
build_glasso <- function(data, actor = NULL, level = NULL, ...) {
  .refuse_method_arg("build_glasso", "glasso", ...names())
  build_network(data, method = "glasso", actor = actor, level = level, ...)
}

#' Build an Ising Network
#'
#' @description
#' `build_ising()` estimates an Ising network for binary (0/1) variables:
#' each variable is regressed on all the others with L1-penalised logistic
#' regression, the penalty is chosen by EBIC, and the two estimates of each
#' edge are combined into one undirected weight (van Borkulo et al., 2014).
#' It is [build_network()] with `method = "ising"`, and needs the glmnet
#' package. The settings go in `params`: `gamma` (the EBIC hyperparameter,
#' default `0.25`), `nlambda` (default `100`) and `rule`, `"AND"` (default:
#' an edge only when both regressions select it) or `"OR"` (when either
#' does).
#'
#' An edge is a conditional association between two variables given all
#' the others; an absent edge means the penalty removed it. Check the
#' stability of the edges with [bootstrap_network()].
#'
#' @param data The binary variables: a data frame with one row per
#'   observation and one column per variable, holding only 0 and 1.
#'   Non-numeric columns and a column named `rid` are set aside; columns
#'   with non-syntactic names, all-missing and constant columns, and rows
#'   with a missing value are dropped, each with a message; fewer than 2
#'   variables or 3 complete rows, or any value other than 0 and 1, raise an
#'   error.
#' @param actor An id column naming whose observation each row is; it is
#'   set aside, not modelled. Default `NULL`.
#' @param ... Further arguments for [build_network()], such as `params`
#'   (`gamma`, `nlambda`, `rule`), `threshold`, `scaling` or `group`.
#'   `method` cannot be given: it is always `"ising"` (error of class
#'   `nestimate_bad_input`).
#' @return A `netobject` (a `netobject_group` with `group`), exactly as
#'   [build_network()] returns it: print it for an overview, read the edges
#'   with [extract_edges()], and draw it with `cograph::splot()`.
#' @references
#' van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo,
#' L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method for
#' constructing networks from binary data. \emph{Scientific Reports}, 4,
#' 5918. \doi{10.1038/srep05918}
#' @seealso \code{\link{build_network}}
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   bin_data <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
#'   net <- build_ising(bin_data)
#' }
#' @export
build_ising <- function(data, actor = NULL, ...) {
  .refuse_method_arg("build_ising", "ising", ...names())
  build_network(data, method = "ising", actor = actor, ...)
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
                "gap", "reverse")
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
