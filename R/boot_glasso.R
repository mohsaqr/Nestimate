# ---- Bootstrap for Regularized Partial Correlation Networks ----

#' Bootstrap for Regularised Partial Correlation Networks
#'
#' @description
#' `boot_glasso()` assesses how accurately an EBIC graphical lasso network
#' is estimated and how stable its centralities are, by resampling the
#' observations (Epskamp, Borsboom & Fried, 2018). One call runs a
#' nonparametric bootstrap, which gives percentile intervals and inclusion
#' probabilities for every edge, intervals for each node's centralities and
#' predictability, and bootstrap tests of the difference between every pair
#' of edges and every pair of nodes; and a case-dropping bootstrap, which
#' gives each centrality measure's correlation-stability (CS) coefficient.
#' It takes the raw data or a network built with `method = "glasso"`.
#'
#' An interval shows how precisely an edge weight is estimated: a wide one
#' means the weight changes a good deal between samples. Epskamp et al.
#' (2018) caution that, because the lasso shrinks estimates towards zero,
#' these intervals are not tests of whether an edge differs from zero; the
#' thresholded network (edges whose interval excludes zero) and the
#' "significant" count in `print()` read them that way, and the inclusion
#' probability (how often an edge is non-zero across resamples) is a
#' separate summary. The CS-coefficient is the largest proportion of cases
#' that can be dropped while more than 95% of the subsamples keep a
#' correlation above 0.7 with the original centralities; Epskamp et al.
#' (2018) suggest it should not be below 0.25 and preferably above 0.5
#' before the centrality order is interpreted, and `print()` labels the
#' values accordingly. The difference tests are not corrected for the
#' number of pairs compared.
#'
#' @details
#' The network is estimated as [build_glasso()] estimates it: with a
#' `netobject`, with its own `params` (`lambda.min.ratio`,
#' `penalize.diagonal`, `refit`); with data, with that function's defaults
#' (smallest lambda 0.01 times the largest, no penalty on the diagonal, the
#' selected penalty re-solved exactly, no unregularised refit). The lambda
#' path is set once from the original data and reused for every resample,
#' which is fitted with the same settings; the penalty is chosen on it by
#' EBIC with `gamma`. Strength is the sum of a node's absolute edge weights and
#' expected influence the sum of its signed weights; betweenness and
#' closeness use shortest paths with distance 1 / |weight|, computed
#' without igraph. Predictability is \eqn{1 - 1 / \Theta_{ii}}, from the
#' diagonal of the precision matrix, bounded to 0--1.
#'
#' The case-dropping bootstrap follows \pkg{bootnet}: each of the `cs_iter`
#' iterations draws one proportion from `cs_drop` at random, keeps that
#' share of the rows (at least `p + 1` of them) without replacement, and
#' correlates the subsample's centralities with the original ones; a
#' centrality that is constant across nodes counts as a correlation of 0.
#' The difference tests give, for each pair, the two-sided bootstrap p-value
#' `2 * min(P(d >= 0), P(d <= 0))` of the difference `d` between the two
#' edges or nodes, capped at 1, with ties counted in both tails. The
#' edge-difference test is skipped for networks with more than 500 edges. A
#' resample whose estimation fails (for example, a variable constant in the
#' resample) is kept as missing and left out of the intervals and the
#' case-dropping results; the failures are counted in `n_failed` and
#' `cs_n_failed`, shown by `print()`, and reported in a warning of class
#' `nestimate_resample_failed`. The function was written as a faster, single-call
#' alternative to running \pkg{bootnet} several times.
#'
#' @param x The data: a data frame or numeric matrix with one row per
#'   observation and one column per variable, cleaned as in
#'   [build_glasso()]: non-numeric columns and a column named `rid` are
#'   set aside; columns with non-syntactic names, all-missing and constant
#'   columns, and rows with a missing value are dropped, each with a
#'   message; fewer than 2 variables or 3 complete rows, or no more rows
#'   than variables, raise an error. A
#'   square symmetric matrix with matching row and column names is read as a
#'   correlation or covariance matrix and refused (class
#'   `nestimate_bad_input`): it has no observations to resample. Or a
#'   `netobject` (or `cograph_network`) built with `method = "glasso"`,
#'   whose data are resampled; any `gamma`, `nlambda` or `cor_method` set in
#'   its `params` replaces the arguments of the same name, and its
#'   `lambda.min.ratio`, `penalize.diagonal` and `refit` are used too, so
#'   `original_pcor` is the network's weights. Any other network raises an
#'   error. A correlation or
#'   covariance matrix is not accepted: a matrix is always read as data,
#'   its rows as observations. For the `print()`
#'   and `plot()` methods: a `boot_glasso`.
#' @param iter The number of nonparametric bootstrap resamples, at least 2.
#'   Default `1000`.
#' @param cs_iter The total number of case-dropping iterations, at
#'   least 1. Default `500`. As in \pkg{bootnet}, each iteration draws one
#'   proportion from `cs_drop` at random, so the iterations are spread
#'   across the proportions (about `cs_iter / length(cs_drop)` each) rather
#'   than repeated `cs_iter` times at each one.
#' @param cs_drop The proportions of cases to drop, each in (0, 1). Default
#'   `seq(0.1, 0.9, by = 0.1)`.
#' @param alpha The significance level: the intervals run from the
#'   `alpha / 2` to the `1 - alpha / 2` quantile, and `print()` and the
#'   `"edge_diff"` / `"centrality_diff"` plots mark the difference tests'
#'   p-values below it as differing. Default `0.05`.
#' @param gamma The EBIC hyperparameter; larger values give sparser
#'   networks. Default `0.5`. Replaced by the value stored in a `netobject`.
#' @param nlambda The number of lambda values in the regularisation path,
#'   at least 2. Default `100`. Replaced by the value stored in a
#'   `netobject`.
#' @param centrality The centrality measures to compute. Default all four
#'   built-in ones: `"strength"` (sum of absolute edge weights),
#'   `"expected_influence"` (sum of signed edge weights), `"betweenness"`
#'   and `"closeness"` (on shortest paths with distance 1 / |weight|).
#'   Without `centrality_fn`, names may be abbreviated. Any other name
#'   needs `centrality_fn`, which then computes it; the built-in measures
#'   are always computed internally.
#' @param centrality_fn A function for measures beyond the four built-in
#'   ones: it is called as `centrality_fn(mat)` on the partial correlation
#'   matrix and returns a named list of centrality vectors, such as
#'   `list(my_metric = ...)`. Default `NULL`; a measure in `centrality`
#'   that is neither built in nor returned by it raises an error.
#' @param cor_method The correlation method: `"pearson"` (default),
#'   `"spearman"` or `"kendall"`. Replaced by the value stored in a
#'   `netobject`.
#' @param ncores The number of cores for the resampling, used through
#'   [parallel::mclapply()], which forks the R process; on Windows it
#'   accepts only 1, and any larger value raises an error. Default `1`
#'   (sequential).
#' @param seed A seed for reproducible resampling; the session's random
#'   state is restored when the call returns. Default `NULL` uses the
#'   current random stream. With `ncores > 1` the resamples are
#'   reproducible only under `RNGkind("L'Ecuyer-CMRG")`.
#'
#' @return A `boot_glasso` object. `print()` reports the data size,
#'   `alpha` and `gamma`, the number of edges whose interval excludes zero,
#'   the mean inclusion probability, the failed resamples when there are
#'   any, each measure's CS-coefficient with its
#'   label, the number of edge pairs that differ, and the timing;
#'   [summary()] returns the table chosen by `type` (edges, centralities,
#'   case-dropping results or predictability; described at the end of this
#'   section); `plot()` draws one figure chosen by `type`. The object's
#'   fields:
#' \describe{
#'   \item{original_pcor}{The partial correlation matrix of the original
#'     data.}
#'   \item{original_precision}{The precision matrix of the original data.}
#'   \item{original_centrality}{A named list with one vector of node values
#'     per centrality measure.}
#'   \item{original_predictability}{A named vector of each node's R-squared.}
#'   \item{edge_ci}{A data frame with one row per edge: `edge`, `weight`,
#'     `ci_lower`, `ci_upper` and `inclusion`.}
#'   \item{edge_inclusion}{A named vector of edge inclusion probabilities.}
#'   \item{thresholded_pcor}{The partial correlation matrix with every edge
#'     whose interval includes zero set to zero.}
#'   \item{centrality_ci}{A named list with one data frame per centrality
#'     measure: `node`, `value`, `ci_lower`, `ci_upper`.}
#'   \item{cs_coefficient}{A named vector with the CS-coefficient of each
#'     centrality measure.}
#'   \item{cs_data}{A data frame of the case-dropping results, one row per
#'     drop proportion and measure: `drop_prop`, `measure`, `mean_cor`,
#'     `prop_above` (the share of that proportion's iterations correlating
#'     above 0.7) and `n_samples` (the iterations that drew that
#'     proportion).}
#'   \item{edge_diff_p}{A symmetric matrix of the edge-difference p-values;
#'     `NULL` when the network has more than 500 edges.}
#'   \item{centrality_diff_p}{A named list with one symmetric matrix of
#'     node-difference p-values per centrality measure.}
#'   \item{predictability_ci}{A data frame with one row per node: `node`,
#'     `r2`, `ci_lower`, `ci_upper`.}
#'   \item{boot_edges}{An `iter` x edges matrix of the bootstrap edge
#'     weights.}
#'   \item{boot_centrality}{A named list with one `iter` x nodes matrix of
#'     bootstrap values per centrality measure.}
#'   \item{boot_predictability}{An `iter` x nodes matrix of bootstrap
#'     R-squared values.}
#'   \item{nodes}{The node names.}
#'   \item{n}{The sample size.}
#'   \item{p}{The number of variables.}
#'   \item{iter, cs_iter, cs_drop, alpha, gamma, nlambda}{The settings used.}
#'   \item{centrality_measures}{The centrality measures computed.}
#'   \item{cor_method}{The correlation method used.}
#'   \item{lambda_path}{The lambda path used for every fit.}
#'   \item{lambda_selected}{The lambda selected on the original data.}
#'   \item{n_failed, cs_n_failed}{The numbers of bootstrap and
#'     case-dropping resamples whose estimation failed and that were left
#'     out.}
#'   \item{timing}{A named vector of elapsed seconds: `total`, `phase1`
#'     (the original fit), `bootstrap`, `case_drop` and `statistics`.}
#' }
#'
#' @examples
#' set.seed(1)
#' dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
#' net <- build_network(dat, method = "glasso")
#' bg <- boot_glasso(net, iter = 10, cs_iter = 5, centrality = "strength")
#' \donttest{
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' net <- build_network(as.data.frame(mat), method = "glasso")
#' # iter = 20 keeps the example fast; a real analysis uses 1000 or more.
#' boot <- boot_glasso(net, iter = 20, cs_iter = 10, seed = 42,
#'   centrality = c("strength", "expected_influence"))
#' print(boot)
#' summary(boot, type = "edges")
#' }
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods} 50(1), 195-212.
#' \doi{10.3758/s13428-017-0862-1}
#' (source of the CS-coefficient and of the case-dropping,
#' edge-difference and centrality-difference procedures reproduced here.)
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}}
#'
#' @importFrom stats cor quantile
#' @importFrom parallel mclapply
#' @export
boot_glasso <- function(x,
                        iter = 1000L,
                        cs_iter = 500L,
                        cs_drop = seq(0.1, 0.9, by = 0.1),
                        alpha = 0.05,
                        gamma = 0.5,
                        nlambda = 100L,
                        centrality = c("strength", "expected_influence",
                                       "betweenness", "closeness"),
                        centrality_fn = NULL,
                        cor_method = "pearson",
                        ncores = 1L,
                        seed = NULL) {

  # ---- Input validation ----
  .assert(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(cs_iter), length(cs_iter) == 1, cs_iter >= 1,
    is.numeric(cs_drop), all(cs_drop > 0), all(cs_drop < 1),
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.numeric(gamma), length(gamma) == 1, gamma >= 0,
    is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2,
    is.character(centrality), length(centrality) >= 1,
    is.character(cor_method), length(cor_method) == 1,
    is.numeric(ncores), length(ncores) == 1, ncores >= 1
  )
  iter <- as.integer(iter)
  cs_iter <- as.integer(cs_iter)
  nlambda <- as.integer(nlambda)
  ncores <- as.integer(ncores)
  cor_method <- match.arg(cor_method, c("pearson", "spearman", "kendall"))
  builtin_centrality <- c("strength", "expected_influence",
                          "closeness", "betweenness")

  if (!is.null(centrality_fn)) {
    .assert("centrality_fn must be a function" = is.function(centrality_fn))
    # With a centrality_fn, custom (non-builtin) names are valid: they are
    # resolved by the function. Only validate the builtin subset against the
    # whitelist; pass custom names through to .bg_compute_centrality, which
    # routes them via centrality_fn.
    is_builtin <- centrality %in% builtin_centrality
    if (any(is_builtin)) {
      centrality[is_builtin] <- match.arg(centrality[is_builtin],
                                          builtin_centrality,
                                          several.ok = TRUE)
    }
  } else {
    centrality <- match.arg(centrality, builtin_centrality,
                            several.ok = TRUE)
  }

  if (!is.null(seed)) {
    .assert(is.numeric(seed), length(seed) == 1)
    .local_seed(seed)
  }

  # ---- Input dispatch ----
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(x, "netobject")) {
    if (x$method != "glasso") {
      stop(errorCondition("netobject must use method = 'glasso'.",
        class = "nestimate_bad_input", call = NULL))
    }
    if (is.null(x$data)) {
      stop(errorCondition(
        "netobject does not contain $data. Rebuild with build_network().",
        class = "nestimate_bad_input", call = NULL))
    }
    data_mat <- x$data
    if (!is.matrix(data_mat)) data_mat <- as.matrix(data_mat)
    # Extract params from netobject if available
    if (!is.null(x$params$gamma)) gamma <- x$params$gamma
    if (!is.null(x$params$nlambda)) nlambda <- as.integer(x$params$nlambda)
    if (!is.null(x$params$cor_method)) cor_method <- x$params$cor_method
    glasso_set <- .bg_glasso_settings(x$params)
  } else if (is.data.frame(x) || is.matrix(x)) {
    # a correlation or covariance matrix has no observations to resample
    if (is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x) &&
        isSymmetric(unname(x)) &&
        identical(rownames(x), colnames(x))) {
      stop(errorCondition(paste0(
        "boot_glasso() resamples observations, so it needs the raw data (one ",
        "row per observation), not a correlation or covariance matrix."),
        class = "nestimate_bad_input", call = NULL))
    }
    # Use .prepare_association_input for cleaning
    # Convert matrix to data frame so it goes through the raw-data path
    input <- if (is.matrix(x)) as.data.frame(x) else x
    prepared <- .prepare_association_input(input, cor_method = cor_method)
    data_mat <- prepared$mat
    glasso_set <- .bg_glasso_settings(list())
  } else {
    stop(errorCondition("x must be a data frame, numeric matrix, or netobject.",
      class = "nestimate_bad_input", call = NULL))
  }

  n <- nrow(data_mat)
  p <- ncol(data_mat)
  nodes <- colnames(data_mat)

  if (p < 2) stop(errorCondition(
    "At least 2 variables are required.",
    class = "nestimate_bad_input", call = NULL))
  if (n < 3) stop(errorCondition("At least 3 observations are required.",
    class = "nestimate_bad_input", call = NULL)) # nocov end
  # Case-dropping keeps at least p + 1 rows of every subsample, which cannot
  # be drawn without replacement from n <= p rows (it used to stop inside
  # sample.int() with "cannot take a sample larger than the population").
  if (n <= p) stop(errorCondition(sprintf(paste0(
    "boot_glasso() needs more observations than variables: the data has ",
    "%d rows and %d variables."), n, p),
    class = "nestimate_bad_input", call = NULL))
  # parallel::mclapply() forks, which Windows cannot do.
  if (ncores > 1L && identical(.Platform$OS.type, "windows")) {
    stop(errorCondition( # nocov start
      "`ncores` > 1 is not supported on Windows. Use `ncores = 1`.",
      class = "nestimate_bad_input", call = NULL)) # nocov end
  }

  t_start <- proc.time()["elapsed"]

  # ---- Phase 1: Original estimation ----
  # Solver delegated to psychnets with the settings build_network(method =
  # "glasso") uses (a glasso netobject's own params, otherwise the
  # estimator's defaults), so original_pcor is the network the user built.
  # The lambda path is fixed once from the original S and reused across
  # every resample, which are fitted with the same settings.
  S <- cor(data_mat, method = cor_method)

  fit_orig <- tryCatch(
    psychnets::ebic_glasso(cor_matrix = S, n = n, gamma = gamma,
                           nlambda = as.integer(nlambda),
                           lambda_min_ratio = glasso_set$lambda_min_ratio,
                           penalize_diagonal = glasso_set$penalize_diagonal,
                           refit = glasso_set$refit, native = TRUE),
    error = function(e) stop(errorCondition(
      paste0("glassopath failed on original data: ", conditionMessage(e)),
      class = "nestimate_bad_input", call = NULL))
  )

  lambda_path <- fit_orig$lambda_path
  Wi_orig <- fit_orig$precision
  lambda_selected <- fit_orig$lambda

  pcor_orig <- .precision_to_pcor(Wi_orig, threshold = 0)
  pcor_orig <- (pcor_orig + t(pcor_orig)) / 2
  colnames(pcor_orig) <- rownames(pcor_orig) <- nodes

  # Original centrality
  cent_orig <- .bg_compute_centrality(pcor_orig, p, nodes, centrality,
                                       centrality_fn)


  # Original predictability
  pred_orig <- pmin(pmax(1 - 1 / diag(Wi_orig), 0), 1)
  names(pred_orig) <- nodes

  # Pre-allocate bootstrap storage
  n_upper <- p * (p - 1) / 2
  ut <- .bg_upper_tri_indices(p)

  t_phase1 <- proc.time()["elapsed"]

  # ---- Phase 2: Nonparametric bootstrap ----
  boot_fn <- function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    .bg_estimate_once(data_mat[idx, , drop = FALSE], p, gamma, nlambda,
                      glasso_set, cor_method, lambda_path,
                      centrality, centrality_fn)
  }

  if (ncores > 1L) {
    boot_results <- parallel::mclapply(seq_len(iter), boot_fn,
                                        mc.cores = ncores)
  } else {
    boot_results <- lapply(seq_len(iter), boot_fn)
  }

  # Unpack results into matrices
  boot_edges <- matrix(NA_real_, nrow = iter, ncol = n_upper)
  boot_centrality <- lapply(centrality, function(m) {
    matrix(NA_real_, nrow = iter, ncol = p,
           dimnames = list(NULL, nodes))
  })
  names(boot_centrality) <- centrality
  boot_pred <- matrix(NA_real_, nrow = iter, ncol = p,
                      dimnames = list(NULL, nodes))

  # A resample whose estimation fails stays a row of NA, left out of the
  # intervals; it is counted, stored and reported, never dropped silently.
  n_failed <- sum(vapply(boot_results, is.null, logical(1)))
  .warn_failed_resamples(n_failed, iter, "boot_glasso() bootstrap",
                         below_threshold = FALSE)

  # Loop kept on purpose: it fills rows of pre-allocated matrices in place.
  for (i in seq_len(iter)) {
    res <- boot_results[[i]]
    if (!is.null(res)) {
      boot_edges[i, ] <- res$edges
      for (m in centrality) {
        boot_centrality[[m]][i, ] <- res$centrality[[m]]
      }
      boot_pred[i, ] <- res$predictability
    }
  }

  edge_names <- .bg_build_edge_names(nodes, ut$row_idx, ut$col_idx)
  colnames(boot_edges) <- edge_names

  t_phase2 <- proc.time()["elapsed"]

  # ---- Phase 3: Case-dropping ----
  cs_result <- .bg_case_drop(
    data_mat, p, gamma, nlambda, cor_method, lambda_path, centrality,
    cs_iter, cs_drop, cent_orig, ncores, centrality_fn, glasso_set
  )
  .warn_failed_resamples(cs_result$n_failed, cs_iter,
                         "boot_glasso() case-dropping",
                         below_threshold = FALSE)

  cs_coef <- vapply(centrality, function(m) {
    .bg_cs_coefficient(cs_result$cors[[m]], cs_drop, 0.7, 0.95)
  }, numeric(1))
  names(cs_coef) <- centrality

  t_phase3 <- proc.time()["elapsed"]

  # ---- Phase 4: Statistics ----

  # Edge CIs
  orig_edge_vec <- pcor_orig[upper.tri(pcor_orig)]
  edge_ci_lower <- apply(boot_edges, 2, quantile, probs = alpha / 2,
                          na.rm = TRUE)
  edge_ci_upper <- apply(boot_edges, 2, quantile, probs = 1 - alpha / 2,
                          na.rm = TRUE)

  # Edge inclusion probability
  edge_inclusion <- colMeans(boot_edges != 0, na.rm = TRUE)
  names(edge_inclusion) <- edge_names

  # Thresholded network: zero out edges where CI includes zero
  thresholded <- .bg_threshold_network(pcor_orig, edge_ci_lower,
                                        edge_ci_upper, ut)

  # Edge CI data frame
  edge_ci_df <- data.frame(
    edge = edge_names,
    weight = orig_edge_vec,
    ci_lower = edge_ci_lower,
    ci_upper = edge_ci_upper,
    inclusion = edge_inclusion,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Centrality CIs
  centrality_ci <- lapply(centrality, function(m) {
    ci_lo <- apply(boot_centrality[[m]], 2, quantile, probs = alpha / 2,
                    na.rm = TRUE)
    ci_hi <- apply(boot_centrality[[m]], 2, quantile, probs = 1 - alpha / 2,
                    na.rm = TRUE)
    data.frame(
      node = nodes,
      value = cent_orig[[m]],
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  names(centrality_ci) <- centrality

  # Edge difference test
  edge_diff_p <- .bg_edge_diff_test(boot_edges)

  # Centrality difference test
  centrality_diff_p <- lapply(centrality, function(m) {
    .bg_centrality_diff_test(boot_centrality[[m]])
  })
  names(centrality_diff_p) <- centrality

  # Predictability CIs
  pred_ci_lower <- apply(boot_pred, 2, quantile, probs = alpha / 2,
                          na.rm = TRUE)
  pred_ci_upper <- apply(boot_pred, 2, quantile, probs = 1 - alpha / 2,
                          na.rm = TRUE)
  pred_ci_df <- data.frame(
    node = nodes,
    r2 = pred_orig,
    ci_lower = pred_ci_lower,
    ci_upper = pred_ci_upper,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # CS data frame
  cs_data <- cs_result$cs_data

  t_end <- proc.time()["elapsed"]

  timing <- c(
    total = unname(t_end - t_start),
    phase1 = unname(t_phase1 - t_start),
    bootstrap = unname(t_phase2 - t_phase1),
    case_drop = unname(t_phase3 - t_phase2),
    statistics = unname(t_end - t_phase3)
  )

  # ---- Assemble result ----
  result <- list(
    original_pcor          = pcor_orig,
    original_precision     = Wi_orig,
    original_centrality    = cent_orig,
    original_predictability = pred_orig,
    edge_ci                = edge_ci_df,
    edge_inclusion         = edge_inclusion,
    thresholded_pcor       = thresholded,
    centrality_ci          = centrality_ci,
    cs_coefficient         = cs_coef,
    cs_data                = cs_data,
    edge_diff_p            = edge_diff_p,
    centrality_diff_p      = centrality_diff_p,
    predictability_ci      = pred_ci_df,
    boot_edges             = boot_edges,
    boot_centrality        = boot_centrality,
    boot_predictability    = boot_pred,
    nodes                  = nodes,
    n                      = n,
    p                      = p,
    iter                   = iter,
    cs_iter                = cs_iter,
    cs_drop                = cs_drop,
    alpha                  = alpha,
    gamma                  = gamma,
    nlambda                = nlambda,
    centrality_measures    = centrality,
    cor_method             = cor_method,
    lambda_path            = lambda_path,
    lambda_selected        = lambda_selected,
    n_failed               = n_failed,
    cs_n_failed            = cs_result$n_failed,
    timing                 = timing
  )

  class(result) <- "boot_glasso"
  result
}


# ---- Internal helpers ----

#' The glasso settings of build_network(method = "glasso")
#'
#' Reads lambda.min.ratio, penalize.diagonal and refit from a network's
#' params with .estimator_glasso()'s defaults, and maps refit as the
#' estimator does: refit = TRUE is psychnets' unregularized refit, refit =
#' FALSE its exact re-solve at the selected penalty.
#' @noRd
.bg_glasso_settings <- function(params) {
  defaults <- formals(.estimator_glasso)
  get <- function(k) .param_get(params, k, eval(defaults[[k]]))
  list(
    lambda_min_ratio  = get("lambda.min.ratio"),
    penalize_diagonal = isTRUE(get("penalize.diagonal")),
    refit             = if (isTRUE(get("refit"))) "unregularized" else TRUE
  )
}

#' Single bootstrap iteration: cor -> EBIC-glasso (psychnets) -> pcor -> centrality
#' @noRd
.bg_estimate_once <- function(data_mat, p, gamma, nlambda, glasso_set,
                               cor_method, lambda_path, centrality,
                               centrality_fn = NULL) {
  n_boot <- nrow(data_mat)
  nodes <- colnames(data_mat)

  S_boot <- tryCatch(
    cor(data_mat, method = cor_method),
    error = function(e) NULL
  )
  if (is.null(S_boot)) return(NULL) # nocov

  # Check for NAs in correlation matrix (can happen with zero-variance
  # columns in bootstrap sample)
  if (any(is.na(S_boot))) return(NULL) # nocov

  # Solver delegated to psychnets; the fixed lambda_path from the original
  # S is reused (never recomputed from the resample), with the original
  # fit's settings. Iteration-level failures stay NULL by the bootstrap
  # contract (they are counted and reported upstream).
  fit <- tryCatch(
    psychnets::ebic_glasso(cor_matrix = S_boot, n = n_boot, gamma = gamma,
                           lambda_path = lambda_path,
                           penalize_diagonal = glasso_set$penalize_diagonal,
                           refit = glasso_set$refit, native = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL) # nocov

  Wi <- fit$precision

  pcor <- .precision_to_pcor(Wi, threshold = 0)
  pcor <- (pcor + t(pcor)) / 2

  # Extract upper triangle
  edges <- pcor[upper.tri(pcor)]

  # Centrality
  cent <- .bg_compute_centrality(pcor, p, nodes, centrality, centrality_fn)

  # Predictability
  pred <- pmin(pmax(1 - 1 / diag(Wi), 0), 1)

  list(edges = edges, centrality = cent, predictability = pred)
}


#' Compute requested centrality measures
#'
#' All four built-ins are dependency-free: strength and expected_influence
#' use rowSums; betweenness and closeness use the internal Floyd-Warshall
#' routines on the absolute partial correlations (undirected, inverted
#' weights). Only measures outside those four need a \code{centrality_fn}.
#'
#' @param pcor Partial correlation matrix.
#' @param p Number of nodes.
#' @param nodes Character vector of node names.
#' @param measures Character vector of requested measures.
#' @param centrality_fn Optional function taking a weight matrix and
#'   returning a named list of centrality vectors. Required only for
#'   measures outside strength/expected_influence/betweenness/closeness.
#' @noRd
.bg_compute_centrality <- function(pcor, p, nodes, measures,
                                    centrality_fn = NULL) {
  # Built-in measures (no dependencies)
  builtin  <- c("strength", "expected_influence",
                "betweenness", "closeness")
  external <- setdiff(measures, builtin)

  result <- vector("list", length(measures))
  names(result) <- measures

  # Strength measures (rowSums based)
  if ("strength" %in% measures) {
    v <- rowSums(abs(pcor)); names(v) <- nodes
    result[["strength"]] <- v
  }
  if ("expected_influence" %in% measures) {
    v <- rowSums(pcor); names(v) <- nodes
    result[["expected_influence"]] <- v
  }

  # Path-based measures (Floyd-Warshall, no igraph needed)
  path_mat <- pcor
  rownames(path_mat) <- colnames(path_mat) <- nodes
  if ("betweenness" %in% measures) {
    result[["betweenness"]] <- .betweenness(abs(path_mat),
                                             directed = FALSE, invert = TRUE)
  }
  if ("closeness" %in% measures) {
    cl <- .closeness(abs(path_mat), directed = FALSE, invert = TRUE)
    result[["closeness"]] <- cl[["Closeness"]]
  }

  if (length(external) > 0L) {
    if (is.null(centrality_fn)) {
      stop(errorCondition(paste0("centrality_fn is required for measures: ",
           paste(external, collapse = ", "),
           ". Provide a function that takes a weight matrix and returns ",
           "a named list of centrality vectors."),
        class = "nestimate_bad_input", call = NULL))
    }
    custom <- centrality_fn(pcor)
    if (!is.list(custom)) {
      stop(errorCondition("centrality_fn must return a named list.",
        class = "nestimate_bad_input", call = NULL)) # nocov
    }
    absent <- external[vapply(external, function(m) is.null(custom[[m]]),
                              logical(1))]
    if (length(absent) > 0L) {
      stop(errorCondition(
        paste0("centrality_fn did not return measure '", absent[1L], "'."),
        class = "nestimate_bad_input", call = NULL)) # nocov end
    }
    result[external] <- lapply(external, function(m) {
      stats::setNames(custom[[m]], nodes)
    })
  }

  result
}


#' Run case-dropping iterations and compute CS data
#'
#' Matches bootnet's approach: each of cs_iter iterations randomly picks one
#' drop proportion from cs_drop, draws a subsample of that size, estimates
#' the network, and computes centrality. The CS-coefficient is then computed
#' as the maximum drop proportion where >95% of samples have correlation
#' with the original centrality exceeding 0.7.
#'
#' @noRd
.bg_case_drop <- function(data_mat, p, gamma, nlambda, cor_method,
                           lambda_path, centrality, cs_iter, cs_drop,
                           original_centralities, ncores,
                           centrality_fn = NULL,
                           glasso_set = .bg_glasso_settings(list())) {
  n <- nrow(data_mat)
  nodes <- colnames(data_mat)
  n_drop <- length(cs_drop)

  # Pre-compute subsample sizes for each drop proportion
  sub_sizes <- vapply(cs_drop, function(dp) {
    as.integer(max(floor(n * (1 - dp)), p + 1L))
  }, integer(1))

  # Each iteration randomly picks a drop proportion (bootnet approach)
  cs_fn <- function(j) {
    d_idx <- sample.int(n_drop, 1L)
    n_sub <- sub_sizes[d_idx]
    idx <- sample.int(n, n_sub, replace = FALSE)
    res <- .bg_estimate_once(data_mat[idx, , drop = FALSE], p, gamma,
                              nlambda, glasso_set, cor_method, lambda_path,
                              centrality, centrality_fn)
    if (is.null(res)) return(list(d_idx = d_idx, cors = NULL)) # nocov

    # Compute Pearson correlation with original for each measure
    # (matches bootnet::cor0 which uses Pearson by default)
    cors_m <- vapply(centrality, function(m) {
      sub_vals <- res$centrality[[m]]
      orig_vals <- original_centralities[[m]]
      if (any(!is.finite(sub_vals)) || any(!is.finite(orig_vals)) ||
          sd(sub_vals, na.rm = TRUE) == 0 || sd(orig_vals, na.rm = TRUE) == 0) {
        return(0)
      }
      tryCatch(
        cor(sub_vals, orig_vals),
        warning = function(w) NA_real_,
        error = function(e) NA_real_
      )
    }, numeric(1))
    names(cors_m) <- centrality

    list(d_idx = d_idx, cors = cors_m)
  }

  if (ncores > 1L) {
    cs_results <- parallel::mclapply(seq_len(cs_iter), cs_fn,
                                      mc.cores = ncores)
  } else {
    cs_results <- lapply(seq_len(cs_iter), cs_fn)
  }

  # Organize results by drop proportion
  # cors: list by measure, each a list of vectors by drop proportion
  cors <- lapply(centrality, function(m) {
    mat <- matrix(NA_real_, nrow = 0, ncol = n_drop)
    # Collect per-proportion vectors
    done <- Filter(function(res) !is.null(res$cors), cs_results)
    done_prop <- vapply(done, function(res) res$d_idx, integer(1))
    done_cor <- vapply(done, function(res) res$cors[[m]], numeric(1))
    lapply(seq_len(n_drop), function(d_idx) {
      unname(done_cor[done_prop == d_idx])
    })
  })
  names(cors) <- centrality

  # Build cs_data: mean correlation and P(cor > 0.7) for each proportion
  cs_rows <- lapply(seq_len(n_drop), function(d_idx) {
    do.call(rbind, lapply(centrality, function(m) {
      vals <- cors[[m]][[d_idx]]
      mean_cor <- if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA_real_
      prop_above <- if (length(vals) > 0) {
        mean(!is.na(vals) & vals > 0.7)
      } else {
        NA_real_
      }
      data.frame(
        drop_prop = cs_drop[d_idx],
        measure = m,
        mean_cor = mean_cor,
        prop_above = prop_above,
        n_samples = length(vals),
        stringsAsFactors = FALSE
      )
    }))
  })

  cs_data <- do.call(rbind, cs_rows)
  rownames(cs_data) <- NULL

  list(cors = cors, cs_data = cs_data,
       n_failed = sum(vapply(cs_results, function(r) is.null(r$cors),
                             logical(1))))
}


#' Compute CS-coefficient from case-dropping correlations
#'
#' Matches bootnet::corStability: CS = max drop proportion where the
#' proportion of bootstrap samples with cor > cor_threshold exceeds
#' prop_threshold (default 0.95).
#'
#' @param cors_per_prop List of numeric vectors, one per drop proportion.
#'   Each vector contains the per-iteration Pearson correlations with the
#'   original centralities (matching bootnet::cor0's default).
#' @param cs_drop Numeric vector of drop proportions.
#' @param cor_threshold Numeric. Correlation threshold (default 0.7).
#' @param prop_threshold Numeric. Required proportion above threshold
#'   (default 0.95).
#' @return Numeric scalar: CS-coefficient.
#' @noRd
.bg_cs_coefficient <- function(cors_per_prop, cs_drop,
                                cor_threshold = 0.7,
                                prop_threshold = 0.95) {
  # For each drop proportion, compute proportion of samples exceeding threshold
  prop_above <- vapply(cors_per_prop, function(vals) {
    if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
    # a draw without a correlation does not reach the threshold
    mean(!is.na(vals) & vals > cor_threshold)
  }, numeric(1))

  # CS = max drop proportion where prop_above > prop_threshold
  above <- !is.na(prop_above) & prop_above > prop_threshold
  if (!any(above)) return(0)
  max(cs_drop[above])
}




#' Pairwise edge difference p-values
#'
#' For each pair of edges (i, j), the two-sided bootstrap p-value
#' 2 * min(P(edge_i >= edge_j), P(edge_i <= edge_j)), capped at 1; ties count
#' in both tails.
#'
#' @noRd
.bg_edge_diff_test <- function(boot_edges) {
  n_edges <- ncol(boot_edges)
  # For efficiency, only compute for manageable number of edges
  if (n_edges > 500) {
    # Return NULL for very large networks
    return(NULL)
  }

  p_mat <- matrix(1, nrow = n_edges, ncol = n_edges,
                   dimnames = list(colnames(boot_edges),
                                    colnames(boot_edges)))

  valid <- complete.cases(boot_edges)
  bm <- boot_edges[valid, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 2) return(p_mat)

  # Loop kept on purpose: one pass per edge holds an n_valid x n_edges block
  # of differences; vectorising over all pairs would hold n_edges of them.
  for (i in seq_len(n_edges - 1L)) {
    js <- seq(i + 1L, n_edges)
    diffs <- bm[, i] - bm[, js, drop = FALSE]
    # A tie (d == 0, e.g. both edges absent) is no evidence of a
    # difference, so it counts in both tails.
    p_greater <- colMeans(diffs >= 0)
    p_less    <- colMeans(diffs <= 0)
    p_vals    <- pmin(1, 2 * pmin(p_greater, p_less))
    p_mat[i, js] <- p_vals
    p_mat[js, i] <- p_vals
  }

  diag(p_mat) <- 0
  p_mat
}


#' Pairwise centrality difference p-values
#' @noRd
.bg_centrality_diff_test <- function(boot_cent) {
  p_nodes <- ncol(boot_cent)
  node_names <- colnames(boot_cent)

  p_mat <- matrix(1, nrow = p_nodes, ncol = p_nodes,
                   dimnames = list(node_names, node_names))

  valid <- complete.cases(boot_cent)
  bm <- boot_cent[valid, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 2) return(p_mat)

  # Loop kept on purpose, as in .bg_edge_diff_test().
  for (i in seq_len(p_nodes - 1L)) {
    js <- seq(i + 1L, p_nodes)
    diffs <- bm[, i] - bm[, js, drop = FALSE]
    # A tie (d == 0, e.g. both edges absent) is no evidence of a
    # difference, so it counts in both tails.
    p_greater <- colMeans(diffs >= 0)
    p_less    <- colMeans(diffs <= 0)
    p_vals    <- pmin(1, 2 * pmin(p_greater, p_less))
    p_mat[i, js] <- p_vals
    p_mat[js, i] <- p_vals
  }

  diag(p_mat) <- 0
  p_mat
}


#' Threshold network: zero out edges whose CI includes zero
#' @noRd
.bg_threshold_network <- function(pcor, ci_lower_vec, ci_upper_vec, ut) {
  result <- pcor
  # An edge is significant if both CI bounds have the same sign and neither is 0
  not_sig <- (sign(ci_lower_vec) != sign(ci_upper_vec)) |
    (ci_lower_vec == 0) | (ci_upper_vec == 0)

  drop_i <- ut$row_idx[which(not_sig)]
  drop_j <- ut$col_idx[which(not_sig)]
  result[cbind(drop_i, drop_j)] <- 0
  result[cbind(drop_j, drop_i)] <- 0

  result
}


#' Pre-compute upper triangle row/col indices
#' @noRd
.bg_upper_tri_indices <- function(p) {
  ut <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  list(row_idx = ut[, 1], col_idx = ut[, 2])
}


#' CS-coefficient label
#' @noRd
.bg_cs_label <- function(cs_value) {
  if (is.na(cs_value)) return("Unknown")
  if (cs_value >= 0.5) return("Stable")
  if (cs_value >= 0.25) return("Acceptable")
  "Unstable"
}


#' Build edge names from upper triangle indices
#' @noRd
.bg_build_edge_names <- function(nodes, row_idx, col_idx) {
  paste(nodes[row_idx], "--", nodes[col_idx])
}


# ---- S3 Methods ----

#' @rdname boot_glasso
#' @param ... For the `plot()` method with `type = "edge_diff"` or
#'   `"centrality_diff"`: `order`, `"sample"` (default) to sort the edges by
#'   their absolute weight or the nodes by their centrality, or `"id"` to
#'   sort them alphabetically; with `type = "edge_diff"` also `pos_color`
#'   and `neg_color`, the fills for a cell where the edge on the x-axis has a
#'   larger or a smaller mean bootstrap weight than the edge on the y-axis
#'   (defaults
#'   `"#C0392B"` and `"#2C6E8A"`). Ignored
#'   otherwise, and by the `print()` and `summary()` methods.
#' @param object For the `summary()` method: a `boot_glasso`.
#' @param type What to return or draw. For the `summary()` method:
#'   `"edges"` (default), the edge intervals sorted by decreasing absolute
#'   weight; `"centrality"`, the centrality intervals; `"cs"`, the
#'   case-dropping table; `"predictability"`, the predictability intervals;
#'   or `"all"`, the four together. For the `plot()` method: `"edges"`
#'   (default), each edge's weight and interval; `"stability"`, the mean
#'   correlation with the original centralities against the proportion
#'   dropped, with a dashed line at 0.7; `"edge_diff"`, a heatmap of the
#'   mean bootstrap difference between every pair of edges, faint where the
#'   pair does not differ at `alpha` and starred where it does, with the
#'   edge weights on the diagonal; `"centrality_diff"`, a heatmap of which
#'   node pairs differ on `measure`, with the node values on the diagonal;
#'   or `"inclusion"`, each edge's inclusion probability, with a dashed line
#'   at 0.5.
#' @param measure For the `plot()` method with `type = "centrality_diff"`:
#'   the centrality measure to compare, one of those computed (names may be
#'   abbreviated). Default `NULL`, the first measure computed.
#' @return `print()` returns its input invisibly. `summary()` returns a data
#'   frame: for `type = "edges"` the `edge_ci` table (`edge`, `weight`,
#'   `ci_lower`, `ci_upper`, `inclusion`) sorted by decreasing absolute
#'   weight, for `"cs"` the `cs_data` table and for `"predictability"` the
#'   `predictability_ci` table; `"centrality"` returns a named list with one
#'   data frame per measure (`node`, `value`, `ci_lower`, `ci_upper`), and
#'   `"all"` a named list of the four (`edges`, `centrality`, `cs`,
#'   `predictability`). `plot()` returns a `ggplot` object, which is drawn
#'   when the call is made at the top level. The edge-difference plot raises
#'   an error when the test was skipped (more than 500 edges).
#' @export
print.boot_glasso <- function(x, ...) {
  cat(sprintf(
    "GLASSO Bootstrap (%d iterations, %d case-drop iterations)\n",
    x$iter, x$cs_iter
  ))
  cat(sprintf("  Data: %d x %d  |  Alpha: %.2f  |  Gamma: %.2f\n",
              x$n, x$p, x$alpha, x$gamma))

  # Count significant edges (CI excludes zero)
  n_upper <- x$p * (x$p - 1) / 2
  n_sig <- sum(x$thresholded_pcor[upper.tri(x$thresholded_pcor)] != 0)
  cat(sprintf("  Edges: %d/%d significant (CI excludes zero)\n",
              n_sig, n_upper))
  cat(sprintf("  Mean inclusion probability: %.2f\n",
              mean(x$edge_inclusion)))
  n_failed <- x$n_failed %||% 0L
  cs_n_failed <- x$cs_n_failed %||% 0L
  if (n_failed > 0L || cs_n_failed > 0L) {
    cat(sprintf(paste0("  Failed resamples (left out): %d of %d bootstrap, ",
                       "%d of %d case-drop\n"),
                n_failed, x$iter, cs_n_failed, x$cs_iter))
  }

  cat("\n  Centrality Stability (CS-coefficient):\n")
  cs_lines <- vapply(x$centrality_measures, function(m) {
    cs_val <- x$cs_coefficient[[m]]
    sprintf("    %-22s %.2f [%s]\n", paste0(m, ":"), cs_val,
            .bg_cs_label(cs_val))
  }, character(1))
  cat(cs_lines, sep = "")

  # Edge differences
  if (!is.null(x$edge_diff_p)) {
    n_edge_pairs <- sum(upper.tri(x$edge_diff_p))
    n_diff_sig <- sum(x$edge_diff_p[upper.tri(x$edge_diff_p)] < x$alpha)
    cat(sprintf("\n  Edge differences: %d/%d pairs significantly different\n",
                n_diff_sig, n_edge_pairs))
  }

  cat(sprintf("  Timing: %.1fs (bootstrap: %.1fs, case-drop: %.1fs)\n",
              x$timing["total"], x$timing["bootstrap"],
              x$timing["case_drop"]))

  invisible(x)
}


#' @rdname boot_glasso
#' @export
summary.boot_glasso <- function(object, type = "edges", ...) {
  type <- match.arg(type, c("edges", "centrality", "cs", "predictability",
                             "all"))

  switch(type,
    edges = {
      df <- object$edge_ci
      df[order(-abs(df$weight)), ]
    },
    centrality = {
      object$centrality_ci
    },
    cs = {
      object$cs_data
    },
    predictability = {
      object$predictability_ci
    },
    all = {
      list(
        edges = {
          df <- object$edge_ci
          df[order(-abs(df$weight)), ]
        },
        centrality = object$centrality_ci,
        cs = object$cs_data,
        predictability = object$predictability_ci
      )
    }
  )
}


#' @rdname boot_glasso
#' @section Methods:
#' * `plot()` draws one figure of the bootstrap results per call, chosen with
#'   `type`.
#' @export
plot.boot_glasso <- function(x, type = "edges", measure = NULL, ...) {
  type <- match.arg(type, c("edges", "stability", "edge_diff",
                              "centrality_diff", "inclusion"))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for this plot type.", call. = FALSE) # nocov
  }
  switch(type,
    edges           = .bg_plot_edges(x, ...),
    stability       = .bg_plot_stability(x, ...),
    edge_diff       = .bg_plot_edge_diff(x, ...),
    centrality_diff = .bg_plot_centrality_diff(x, measure, ...),
    inclusion       = .bg_plot_inclusion(x, ...)
  )
}


#' Edge CI plot: sorted pointrange
#' @noRd
.bg_plot_edges <- function(x, ...) {
  df <- x$edge_ci
  df <- df[order(abs(df$weight)), ]
  df$edge <- factor(df$edge, levels = df$edge)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight, y = .data$edge)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      size = 0.3
    ) +
    ggplot2::labs(
      title = "Edge Weight Bootstrap CIs",
      x = "Edge Weight",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' CS stability plot: drop curves with threshold line
#' @noRd
.bg_plot_stability <- function(x, ...) {
  cs_data <- x$cs_data

  ggplot2::ggplot(cs_data,
    ggplot2::aes(x = .data$drop_prop, y = .data$mean_cor,
                 color = .data$measure, group = .data$measure)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed",
                         color = "red") +
    ggplot2::scale_x_continuous(breaks = x$cs_drop) +
    ggplot2::labs(
      title = "Centrality Stability (Case-Dropping Bootstrap)",
      x = "Proportion Dropped",
      y = "Mean Correlation with Original",
      color = "Measure"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Edge difference heatmap
#'
#' Upper-triangle heatmap coloured by the mean bootstrap difference between
#' each edge pair. Red = column (x-axis) edge stronger; blue = row (y-axis)
#' edge stronger.
#' Significant cells (p < alpha) are shown at full intensity with p-value
#' stars; non-significant cells are nearly invisible. Diagonal shows the
#' sample edge weight.
#'
#' @param x boot_glasso object.
#' @param order Character. \code{"sample"} (default) orders edges by
#'   absolute sample weight; \code{"id"} orders alphabetically.
#' @param pos_color Colour for positive differences (col edge > row edge).
#'   Default \code{"#C0392B"} (crimson).
#' @param neg_color Colour for negative differences (row edge > col edge).
#'   Default \code{"#2C6E8A"} (teal).
#' @param ... Ignored.
#' @noRd
.bg_plot_edge_diff <- function(x,
                                order     = c("sample", "id"),
                                pos_color = "#C0392B",
                                neg_color = "#2C6E8A",
                                ...) {
  if (is.null(x$edge_diff_p))
    stop(errorCondition(
      "Edge difference test was not computed (too many edges).",
      class = "nestimate_bad_input", call = NULL))

  order <- match.arg(order)

  p_mat      <- x$edge_diff_p
  alpha      <- x$alpha
  boot_mat   <- x$boot_edges
  edge_names <- colnames(p_mat)
  n_e        <- length(edge_names)

  # Mean bootstrap differences: E[edge_i - edge_j] = mean(edge_i) - mean(edge_j)
  edge_means    <- colMeans(boot_mat, na.rm = TRUE)
  mean_diff_mat <- outer(edge_means, edge_means, `-`)
  rownames(mean_diff_mat) <- colnames(mean_diff_mat) <- edge_names

  # Edge weight lookup (diagonal labels)
  edge_weights <- x$edge_ci$weight
  names(edge_weights) <- x$edge_ci$edge

  # Ordering
  if (order == "sample") {
    ordered_names <- edge_names[order(abs(edge_weights[edge_names]))]
  } else {
    ordered_names <- sort(edge_names)
  }

  # Build full grid
  grid <- expand.grid(edge1 = ordered_names, edge2 = ordered_names,
                      stringsAsFactors = FALSE)
  idx1 <- match(grid$edge1, ordered_names)
  idx2 <- match(grid$edge2, ordered_names)

  is_upper <- idx1 < idx2
  is_diag  <- idx1 == idx2

  grid$p_value   <- p_mat[cbind(grid$edge1, grid$edge2)]
  grid$mean_diff <- mean_diff_mat[cbind(grid$edge1, grid$edge2)]
  grid$sig       <- grid$p_value < alpha

  # Fill value: mean_diff for upper triangle, 0 (midpoint) otherwise
  grid$fill_val  <- ifelse(is_upper | is_diag, grid$mean_diff, NA_real_)
  grid$fill_val[is_diag] <- 0

  # Alpha: sig upper = 1, non-sig upper = 0.08, diagonal = 0.06, lower = 0
  grid$alpha_val <- 0
  grid$alpha_val[is_upper &  grid$sig] <- 1
  grid$alpha_val[is_upper & !grid$sig] <- 0.08
  grid$alpha_val[is_diag]              <- 0.06

  # Significance stars for upper-triangle significant cells
  grid$stars <- ""
  grid$stars[is_upper & grid$sig] <- vapply(
    grid$p_value[is_upper & grid$sig],
    function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else "*",
    character(1)
  )

  # Diagonal weight labels
  grid$label <- ifelse(is_diag,
                        sprintf("%.2f", edge_weights[grid$edge1]), "")

  grid$edge1 <- factor(grid$edge1, levels = ordered_names)
  grid$edge2 <- factor(grid$edge2, levels = rev(ordered_names))

  label_size <- if (n_e <= 10) 3.2 else if (n_e <= 20) 2.5 else 1.8
  star_size  <- label_size * 0.80

  # Subtitle
  ut     <- which(upper.tri(p_mat), arr.ind = TRUE)
  n_sig  <- sum(p_mat[ut] < alpha)
  n_pairs <- nrow(ut)
  subtitle <- sprintf(
    "%d of %d pairs significantly different (p < %s)  |  Red: column > row  |  Blue: row > column",
    n_sig, n_pairs, alpha
  )

  # Abs max for symmetric scale
  diff_max <- max(abs(grid$fill_val), na.rm = TRUE)
  if (diff_max == 0) diff_max <- 0.1

  ggplot2::ggplot(grid, ggplot2::aes(x = .data$edge1, y = .data$edge2)) +
    # White background for all cells
    ggplot2::geom_tile(fill = "white", colour = "#F0F0F0", linewidth = 0.25) +
    # Coloured tiles (significant: full; non-sig: faint; diagonal: near-white)
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$fill_val, alpha = .data$alpha_val),
      colour = "white", linewidth = 0.25
    ) +
    # Stars on significant cells
    ggplot2::geom_text(
      data = grid[is_upper & grid$sig, ],
      ggplot2::aes(label = .data$stars),
      size = star_size, colour = "white", fontface = "bold"
    ) +
    # Edge weight on diagonal
    ggplot2::geom_text(
      data = grid[is_diag, ],
      ggplot2::aes(label = .data$label),
      size = label_size * 0.88, colour = "#444444"
    ) +
    ggplot2::scale_fill_gradient2(
      low      = neg_color,
      high     = pos_color,
      mid      = "white",
      midpoint = 0,
      limits   = c(-diff_max, diff_max),
      na.value = "transparent",
      name     = "Mean\ndifference"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(
      title    = "Bootstrap Edge Difference Test",
      subtitle = subtitle,
      x = NULL, y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(
        angle  = 45, hjust = 1, vjust = 1,
        size   = label_size * 2.8, colour = "#333333"
      ),
      axis.text.y      = ggplot2::element_text(
        size   = label_size * 2.8, colour = "#333333"
      ),
      plot.title       = ggplot2::element_text(
        size = 13, face = "bold", colour = "#1A1A1A",
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle    = ggplot2::element_text(
        size = 9, colour = "#666666",
        margin = ggplot2::margin(b = 10)
      ),
      legend.title     = ggplot2::element_text(size = 8.5, colour = "#444444"),
      legend.text      = ggplot2::element_text(size = 8, colour = "#555555"),
      legend.key.width = ggplot2::unit(0.5, "cm"),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin      = ggplot2::margin(12, 12, 12, 12)
    )
}


#' Centrality difference heatmap
#'
#' Upper triangle heatmap with red (significant) / green (non-significant)
#' tiles and sample centrality values on the diagonal.
#'
#' @param x boot_glasso object.
#' @param measure Character. Centrality measure to plot.
#' @param order Character. \code{"sample"} (default) orders nodes by
#'   centrality value (ascending, highest at top-right); \code{"id"}
#'   orders alphabetically.
#' @param ... Ignored.
#' @noRd
.bg_plot_centrality_diff <- function(x, measure = NULL,
                                      order = c("sample", "id"), ...) {
  if (is.null(measure)) measure <- x$centrality_measures[1]
  measure <- match.arg(measure, x$centrality_measures)
  order <- match.arg(order)

  p_mat <- x$centrality_diff_p[[measure]]
  alpha <- x$alpha
  node_names <- colnames(p_mat)
  cent_vals <- x$original_centrality[[measure]]

  # Determine ordering
  if (order == "sample") {
    ordered_names <- node_names[order(cent_vals[node_names])]
  } else {
    ordered_names <- sort(node_names)
  }

  # Build grid: upper triangle + diagonal only
  grid <- expand.grid(
    node1 = ordered_names,
    node2 = ordered_names,
    stringsAsFactors = FALSE
  )
  idx1 <- match(grid$node1, ordered_names)
  idx2 <- match(grid$node2, ordered_names)

  grid$p_value <- vapply(seq_len(nrow(grid)), function(i) {
    p_mat[grid$node1[i], grid$node2[i]]
  }, numeric(1))

  is_diag <- idx1 == idx2
  is_upper <- idx1 < idx2
  is_lower <- idx1 > idx2

  grid$fill <- ifelse(
    is_diag, "diagonal",
    ifelse(is_lower, "blank",
      ifelse(grid$p_value < alpha, "significant", "non-significant")
    )
  )

  # Diagonal labels: sample centrality value
  grid$label <- ifelse(is_diag,
    sprintf("%.2f", cent_vals[grid$node1]), "")

  grid$node1 <- factor(grid$node1, levels = ordered_names)
  grid$node2 <- factor(grid$node2, levels = rev(ordered_names))

  n_p <- length(ordered_names)
  label_size <- if (n_p <= 8) 3.5 else if (n_p <= 15) 2.5 else 1.8

  # Count significant pairs
  ut <- which(upper.tri(p_mat), arr.ind = TRUE)
  n_sig <- sum(p_mat[ut] < alpha)
  n_pairs <- nrow(ut)
  subtitle <- sprintf("%d of %d pairs significantly different (p < %s)",
                       n_sig, n_pairs, alpha)

  ggplot2::ggplot(grid,
    ggplot2::aes(x = .data$node1, y = .data$node2, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      size = label_size, colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "significant"     = "#E74C3C",
        "non-significant" = "#2ECC71",
        "diagonal"        = "white",
        "blank"           = "white"
      ),
      labels = c(
        "significant"     = sprintf("p < %s", alpha),
        "non-significant" = sprintf("p >= %s", alpha)
      ),
      breaks = c("significant", "non-significant"),
      name = NULL
    ) +
    ggplot2::labs(
      title = sprintf("Centrality Difference Test (%s)", measure),
      subtitle = subtitle, x = NULL, y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text()
    )
}


#' Inclusion probability bar plot
#' @noRd
.bg_plot_inclusion <- function(x, ...) {
  df <- data.frame(
    edge = names(x$edge_inclusion),
    inclusion = as.numeric(x$edge_inclusion),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$inclusion), ]
  df$edge <- factor(df$edge, levels = df$edge)

  ggplot2::ggplot(df,
    ggplot2::aes(x = .data$inclusion, y = .data$edge)) +
    ggplot2::geom_col(fill = "#377EB8") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                         color = "grey50") +
    ggplot2::labs(
      title = "Edge Inclusion Probability",
      x = "Inclusion Probability",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


