# ---- Markov Order Test (within-context permutation) ----
#
# Principled test of Markov order for categorical sequence data.
#
# At each order k, test H0: "process is order (k-1)" against H1: "process
# is order k" by reframing as a conditional-independence test on
# (k+1)-grams (x, w, s) where:
#   - w is the (k-1)-gram context
#   - x is the additional (k-th-back) state
#   - s is the next state
# Under H0, s is independent of x given w.
#
# Test statistic: classical LRT for conditional independence
#   G^2 = 2 sum_{x,w,s} N(x,w,s) log[ N(x,w,s) / E(x,w,s) ]
# with E(x,w,s) = N(x,w,.) * N(.,w,s) / N(.,w,.),
# df = sum_w (r_w - 1)(c_w - 1).
#
# Null distribution: within-w (conditional) permutation. For each context w,
# the successor labels are exchangeable under H0; shuffling s within
# each w-group approximates the reference distribution (the windows of a
# sequence overlap, so it is not exact; Besag & Mondal 2013). No plug-in MLE
# bias (which inflated the parametric-bootstrap null on sparse k+1-gram
# tables). No refitting per replicate.


# ---------------------------------------------------------------------------
# Internal: Extract (x, w, s) tuples from trajectories
# ---------------------------------------------------------------------------

#' @noRd
.mot_extract_tuples <- function(trajectories, k) {
  stopifnot(k >= 1L)
  win <- k + 1L

  tuples <- lapply(trajectories, function(traj) {
    n <- length(traj)
    if (n < win) return(NULL)
    pos <- seq_len(n - k)
    x <- traj[pos]
    s <- traj[pos + k]
    w <- if (k == 1L) {
      rep("", length(pos))
    } else {
      vapply(pos, function(p) paste(traj[(p + 1L):(p + k - 1L)],
                                    collapse = .HON_SEP),
             character(1L))
    }
    data.frame(x = x, w = w, s = s, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, tuples[!vapply(tuples, is.null, logical(1L))])
  if (is.null(out)) {
    data.frame(x = character(0L), w = character(0L), s = character(0L),
               stringsAsFactors = FALSE)
  } else out
}


# ---------------------------------------------------------------------------
# Internal: Split trajectories at missing values
# ---------------------------------------------------------------------------

#' @noRd
.mot_split_na <- function(trajectories) {
  pieces <- unlist(lapply(trajectories, function(traj) {
    observed <- !is.na(traj)
    unname(split(traj[observed], cumsum(!observed)[observed]))
  }), recursive = FALSE)
  pieces[vapply(pieces, length, integer(1L)) >= 2L]
}


# ---------------------------------------------------------------------------
# Internal: G^2 statistic and df for one w-table
# ---------------------------------------------------------------------------

#' @noRd
.mot_g2_one_w <- function(x, s) {
  tab <- table(x, s)
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(c(stat = 0, df = 0))
  rs <- rowSums(tab); cs <- colSums(tab); total <- sum(tab)
  if (total == 0) return(c(stat = 0, df = 0))
  expt <- outer(rs, cs) / total
  nz <- tab > 0 & expt > 0
  g2 <- 2 * sum(tab[nz] * log(tab[nz] / expt[nz]))
  df <- (nrow(tab) - 1L) * (ncol(tab) - 1L)
  c(stat = g2, df = df)
}


#' @noRd
.mot_g2 <- function(tuples) {
  if (nrow(tuples) == 0L) return(list(stat = 0, df = 0L))
  by_w <- split(seq_len(nrow(tuples)), tuples$w)
  parts <- vapply(by_w, function(idx) {
    .mot_g2_one_w(tuples$x[idx], tuples$s[idx])
  }, numeric(2L))
  list(stat = sum(parts["stat", ]),
       df   = as.integer(sum(parts["df", ])))
}


# ---------------------------------------------------------------------------
# Internal: Within-w permutation null
# ---------------------------------------------------------------------------

#' @noRd
.mot_permute_null <- function(tuples, n_perm, parallel = FALSE, n_cores = 2L) {
  if (nrow(tuples) == 0L) return(numeric(0L))
  by_w <- split(seq_len(nrow(tuples)), tuples$w)
  s_orig <- tuples$s
  x <- tuples$x

  one_perm <- function(.) {
    # Shuffle the successors within each context, in the order of by_w
    s <- s_orig
    s[unlist(by_w, use.names = FALSE)] <- unlist(
      lapply(by_w, function(idx) {
        if (length(idx) > 1L) sample(s_orig[idx]) else s_orig[idx]
      }),
      use.names = FALSE
    )
    # Recompute G^2 directly on shuffled s (x, w unchanged)
    parts <- vapply(by_w, function(idx) {
      .mot_g2_one_w(x[idx], s[idx])
    }, numeric(2L))
    sum(parts["stat", ])
  }

  map_fn <- if (parallel && .Platform$OS.type != "windows") {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else lapply

  unlist(map_fn(seq_len(n_perm), one_perm), use.names = FALSE)
}


# ---------------------------------------------------------------------------
# Internal: Fit models + compute per-order log-likelihood (for IC panel)
# ---------------------------------------------------------------------------

#' @noRd
.mot_fit_models <- function(trajectories, max_order, dof = "paths") {
  marginal <- .mogen_marginal(trajectories)
  trans_mats <- vector("list", max_order + 1L)
  trans_mats[[1L]] <- marginal
  layer_dofs <- numeric(max_order + 1L)
  layer_dofs[1L] <- length(marginal) - 1L

  lapply(seq_len(max_order), function(k) {
    kg <- .mogen_count_kgrams(trajectories, k)
    tm <- .mogen_transition_matrix(kg$nodes, kg$edges)
    trans_mats[[k + 1L]] <<- tm
    layer_dofs[k + 1L] <<- .mogen_layer_dof(tm)
    NULL
  })
  if (dof == "paths") {
    layer_dofs[-1L] <- .mogen_path_dofs(trans_mats[[2L]], max_order)
  }

  logliks <- vapply(0L:max_order, function(k) {
    .mogen_log_likelihood(trajectories, k, trans_mats[seq_len(k + 1L)])
  }, numeric(1L))

  list(trans_mats = trans_mats, marginal = marginal,
       layer_dofs = layer_dofs, logliks = logliks)
}


# ---------------------------------------------------------------------------
# Main: markov_order_test()
# ---------------------------------------------------------------------------

#' Test the Markov order of a sequential process
#'
#' @description
#' `markov_order_test()` tests how many previous states are needed to
#' predict the next one in a set of categorical sequences: the order of the
#' Markov chain that describes them. At each order
#' \eqn{k = 1, \ldots,} `max_order` it tests the null hypothesis that order
#' \eqn{k - 1} is enough, with the likelihood-ratio statistic \eqn{G^2} for
#' the conditional independence \eqn{s \perp x \mid w}{s _||_ x | w},
#' where \eqn{w} is the \eqn{(k-1)}-gram context, \eqn{x} the extra state
#' \eqn{k} steps back and \eqn{s} the next state. The p-value comes from a
#' permutation null, with the asymptotic \eqn{\chi^2} p-value beside it for
#' reference, and the order is chosen sequentially. It reads wide sequence
#' tables, lists of sequences, sequence networks (single or grouped) and
#' [prepare()] results.
#'
#' The order is raised while the test rejects and stops at the first
#' non-rejection, so the selected order is the highest \eqn{k} whose test,
#' and every test below it, rejected at `alpha`: `0` when order 1 is
#' already not rejected, and `max_order` when every test rejects. AIC and
#' BIC of the fitted chains are reported beside it and can select a
#' different order. With many observations even small departures from the
#' lower order are rejected, so compare the log-likelihoods of neighbouring
#' orders as well. A first-order network from [build_network()] assumes
#' order 1; when a higher order is selected, [build_hon()] and
#' [build_mogen()] model the longer memory.
#'
#' @details
#' Under the null hypothesis the successors are exchangeable within each
#' context \eqn{w}, so shuffling \eqn{s} within each \eqn{w}-group keeps
#' every count the null hypothesis conditions on, with no plug-in estimate
#' and no refit per permutation. The windows of a sequence overlap, so
#' this approximates the reference distribution of \eqn{G^2} rather than
#' giving it exactly (an exact test would sample whole sequences with the
#' same counts; Besag & Mondal, 2013). The permutation p-value is
#' `(number of permutations with a G^2 at least the observed + 1) /
#' (n_perm + 1)`. An order at which no context holds both two different
#' extra states and two different successors has no degrees of freedom; it gets
#' no p-value and counts as not rejected.
#'
#' The AIC and BIC columns come from the log-likelihood of each order's
#' fitted chain and its cumulative number of parameters (see `dof`); BIC
#' takes the number of observations (observed states) as the sample size.
#' The deterministic columns agree with `tna::markov_order()` (tna 1.4.1
#' or later).
#'
#' @param data The sequences: a data frame with one row per sequence and
#'   one column per position, a list of character vectors (one per
#'   sequence), a sequence `netobject` or `netobject_group` built with its
#'   data, or a [prepare()] result (its sequence data is used). A data
#'   frame is read whole, so it must hold only positions: one column, or
#'   columns named as a numbered run (`T1`, `T2`, ...); anything else, such
#'   as an id column beside the positions, a column named `action`, or
#'   columns holding only 0 and 1, raises an error of class
#'   `nestimate_bad_input`, and the positions are then declared with
#'   [build_network()] (`state_cols`, or `action` for an event log) and the
#'   network passed instead. A one-hot or association network raises the
#'   same error. Sequences are split at missing values: no window spans an
#'   `NA`, and pieces with fewer than two states are dropped.
#' @param max_order The highest order to test, at least 1. Default `3`. It
#'   is lowered, with a message, to one less than the length of the longest
#'   sequence.
#' @param n_perm The number of within-context permutations per order, at
#'   least 1. Default `500`. The smallest attainable p-value is
#'   `1 / (n_perm + 1)`, and a test rejects only when its p-value is below
#'   `alpha`, so with `1 / alpha - 1` permutations or fewer (19 for the
#'   default `alpha`) no order above 0 can be selected.
#' @param alpha The significance level for the sequential selection,
#'   strictly between 0 and 1. Default `0.05`.
#' @param parallel Whether to run the permutations with
#'   [parallel::mclapply()]. Default `FALSE`. Ignored on Windows, where
#'   they run serially. With `TRUE` the permutations are reproducible from
#'   `seed` only under `RNGkind("L'Ecuyer-CMRG")`; under the default
#'   generator they differ from run to run.
#' @param n_cores The number of cores for `parallel = TRUE`. Default `2`.
#' @param seed An integer that makes the permutations reproducible (with
#'   `parallel = TRUE`, only under the L'Ecuyer-CMRG generator). Default
#'   `NULL`. The session's own random stream is left as it was; a group's
#'   networks are each
#'   started from the same seed.
#' @param dof How the AIC and BIC columns count each layer's parameters:
#'   `"paths"` (default) counts every transition the observed first-order
#'   topology allows, as in Scholtes (2017) and pathpy; `"observed"` counts
#'   the transitions seen in the data, the rule before Nestimate 0.9.28. The
#'   permutation test does not use it.
#' @return A `net_markov_order` object. `print()` shows the data size, the
#'   order each criterion selects and the test table; `summary()` returns
#'   the test table with the selected orders, `alpha` and `n_perm` as
#'   attributes; `plot()` draws the fit across orders and the permutation
#'   null of each test, and returns a ggplot, a gtable or a list of the two
#'   ggplots (see `panel` and `combined`). For a `netobject_group`, the
#'   result is a `net_markov_order_group`, a named list of
#'   `net_markov_order` objects, one per group, whose `print()` prints
#'   each. The `print()` methods return their input invisibly. The
#'   `net_markov_order`'s fields:
#' \describe{
#'   \item{optimal_order}{The order selected by the sequential permutation
#'     test.}
#'   \item{bic_order, aic_order}{The orders with the lowest BIC and AIC,
#'     shown by `print()` and `plot()` and not used in the selection.}
#'   \item{test_table}{The test table, one row per order from 0 to
#'     `max_order`, with columns `order`, `loglik`, `AIC`, `BIC`, `df`,
#'     `g2`, `p_permutation`, `p_asymptotic` and `significant`. The order-0
#'     row has `NA` in the test columns (`df` to `significant`).}
#'   \item{permutation_null}{The permutation null of \eqn{G^2} at each
#'     order, a list of numeric vectors named `order_1`, `order_2`, ...}
#'   \item{logliks}{The log-likelihood of each order, as a named numeric
#'     vector (for the AIC and BIC columns only; not used in the test).}
#'   \item{layer_dofs}{The parameters added at each order, as a named
#'     numeric vector, counted as set by `dof`.}
#'   \item{dof_method}{`"paths"` or `"observed"`.}
#'   \item{transition_matrices}{The fitted transition matrices, the
#'     order-0 state distribution first.}
#'   \item{states}{The observed states.}
#'   \item{n_sequences, n_observations}{The number of sequences (after
#'     splitting at missing values) and of observed states.}
#'   \item{n_perm, alpha, max_order}{The settings used; `max_order` is the
#'     highest order actually tested.}
#' }
#'
#' @references
#' Besag, J., & Mondal, D. (2013). Exact goodness-of-fit tests for Markov
#' chains. \emph{Biometrics}, 69(2), 488--496. \doi{10.1111/biom.12009}
#'
#' Scholtes, I. (2017). When is a network a network? Multi-order graphical
#' model selection in pathways and temporal networks. In \emph{Proceedings
#' of the 23rd ACM SIGKDD International Conference on Knowledge Discovery
#' and Data Mining} (pp. 1037--1046). \doi{10.1145/3097983.3098145}
#'
#' @examples
#' # Is one previous state enough to predict the next one?
#' res <- markov_order_test(as.data.frame(trajectories),
#'                          max_order = 2, n_perm = 99, seed = 1)
#' res
#' summary(res)
#' plot(res)
#' @export
markov_order_test <- function(data, max_order = 3L, n_perm = 500L, alpha = 0.05,
                               parallel = FALSE, n_cores = 2L, seed = NULL,
                               dof = c("paths", "observed")) {
  dof <- match.arg(dof)
  .refuse_onehot_network(data, "markov_order_test()")
  if (is.data.frame(data)) {
    .refuse_undeclared_roles(data, "markov_order_test()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  if (inherits(data, "netobject_group")) {
    out <- lapply(data, function(net) {
      seq_data <- net$data
      if (is.null(seq_data)) {
        stop(errorCondition(
          "netobject_group member has no $data; rebuild with build_network().",
          class = "nestimate_bad_input", call = NULL))
      }
      markov_order_test(net, max_order = max_order, n_perm = n_perm,
                        alpha = alpha, parallel = parallel,
                        n_cores = n_cores, seed = seed, dof = dof)
    })
    class(out) <- c("net_markov_order_group", "list")
    return(out)
  }
  if (inherits(data, "netobject")) {
    if (is.null(data$data)) {
      stop(errorCondition(
        "netobject has no $data; rebuild with build_network().",
        class = "nestimate_bad_input", call = NULL))
    }
    data <- data$data
  }
  max_order <- as.integer(max_order)
  n_perm    <- as.integer(n_perm)
  .assert(
    "'max_order' must be >= 1" = max_order >= 1L,
    "'n_perm' must be >= 1"    = n_perm >= 1L,
    "'alpha' must be in (0, 1)" = alpha > 0 && alpha < 1,
    is.logical(parallel)
  )
  if (!is.null(seed)) .local_seed(as.integer(seed))

  # Accept a prepare() result directly: use its $sequence_data (wide df)
  if (is.list(data) && !is.data.frame(data) && "sequence_data" %in% names(data)) {
    data <- data$sequence_data
  }

  trajectories <- .mot_split_na(
    .hon_parse_input(data, collapse_repeats = FALSE)
  )
  n_seqs <- length(trajectories)
  .assert("need >= 1 sequence after parsing" = n_seqs >= 1L)
  seq_lens <- vapply(trajectories, length, integer(1L))

  longest <- max(seq_lens)
  if (max_order >= longest) {
    max_order <- longest - 1L
    message(sprintf("max_order capped at %d (longest sequence = %d)",
                    max_order, longest))
  }

  fits <- .mot_fit_models(trajectories, max_order, dof = dof)

  # --- Per-order test ---
  per_order <- lapply(seq_len(max_order), function(k) {
    tuples <- .mot_extract_tuples(trajectories, k)
    g2 <- .mot_g2(tuples)
    if (nrow(tuples) > 0L && g2$df > 0L) {
      null_stats <- .mot_permute_null(tuples, n_perm, parallel, n_cores)
      p_perm  <- (sum(null_stats >= g2$stat) + 1) / (length(null_stats) + 1)
      p_asymp <- stats::pchisq(g2$stat, g2$df, lower.tail = FALSE)
    } else {
      null_stats <- numeric(0L)
      p_perm <- NA_real_; p_asymp <- NA_real_
    }
    list(stat = g2$stat, df = g2$df, null = null_stats,
         p_perm = p_perm, p_asymp = p_asymp, n_tuples = nrow(tuples))
  })

  # Sequential selection
  stat_vec   <- vapply(per_order, `[[`, numeric(1L), "stat")
  df_vec     <- vapply(per_order, `[[`, numeric(1L), "df")
  p_perm_vec <- vapply(per_order, `[[`, numeric(1L), "p_perm")
  p_asy_vec  <- vapply(per_order, `[[`, numeric(1L), "p_asymp")

  sig_vec <- !is.na(p_perm_vec) & p_perm_vec < alpha
  # The order is raised while the test rejects
  optimal_order <- if (all(sig_vec)) max_order else which(!sig_vec)[1L] - 1L

  orders_all <- 0L:max_order
  cum_dof    <- cumsum(fits$layer_dofs)
  aic_vec    <- 2 * cum_dof - 2 * fits$logliks
  bic_vec    <- log(sum(seq_lens)) * cum_dof - 2 * fits$logliks

  test_table <- data.frame(
    order         = orders_all,
    loglik        = as.numeric(fits$logliks),
    AIC           = as.numeric(aic_vec),
    BIC           = as.numeric(bic_vec),
    df            = c(NA_integer_, as.integer(df_vec)),
    g2            = c(NA_real_,    as.numeric(stat_vec)),
    p_permutation = c(NA_real_,    as.numeric(p_perm_vec)),
    p_asymptotic  = c(NA_real_,    as.numeric(p_asy_vec)),
    significant   = c(NA,          as.logical(sig_vec)),
    stringsAsFactors = FALSE,
    row.names     = NULL
  )

  null_list <- lapply(per_order, `[[`, "null")
  names(null_list) <- paste0("order_", seq_len(max_order))
  names(fits$logliks) <- paste0("order_", orders_all)
  names(fits$layer_dofs) <- paste0("order_", orders_all)

  bic_order <- as.integer(test_table$order[which.min(test_table$BIC)])
  aic_order <- as.integer(test_table$order[which.min(test_table$AIC)])

  structure(
    list(
      optimal_order       = as.integer(optimal_order),
      bic_order           = bic_order,
      aic_order           = aic_order,
      test_table          = test_table,
      permutation_null    = null_list,
      logliks             = fits$logliks,
      layer_dofs          = fits$layer_dofs,
      dof_method          = dof,
      transition_matrices = fits$trans_mats,
      states              = names(fits$marginal),
      n_sequences         = n_seqs,
      n_observations      = sum(seq_lens),
      n_perm              = n_perm,
      alpha               = alpha,
      max_order           = max_order
    ),
    class = "net_markov_order"
  )
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' @rdname markov_order_test
#' @param x For the `print()` methods: a `net_markov_order` or
#'   `net_markov_order_group`; for the `plot()` method: a
#'   `net_markov_order`.
#' @param ... For the `print()` method of a `net_markov_order_group`:
#'   passed on to the `print()` of each group. For the other methods: not
#'   used.
#' @param object For the `summary()` method: a `net_markov_order`.
#' @param panel For the `plot()` method: which panels to draw. `"both"`
#'   (default), `"ic"` (the fit across orders) or `"permutation"` (the
#'   permutation nulls); a single panel is returned as a ggplot.
#' @param combined For the `plot()` method with `panel = "both"`: `TRUE`
#'   (default) draws the two panels side by side and returns them as one
#'   gtable when \pkg{gridExtra} is installed, or, without it, draws them
#'   with \pkg{grid} viewports and returns the two ggplots (`ic`,
#'   `permutation`) invisibly as a named list. `FALSE` returns that list
#'   invisibly without drawing. Ignored for a single panel.
#' @export
print.net_markov_order <- function(x, ...) {
  cat(sprintf("Markov Order Test  [within-w permutation, n_perm = %d, alpha = %.3f]\n",
              x$n_perm, x$alpha))
  cat(sprintf("  %d sequences / %d observations / %d states\n\n",
              x$n_sequences, x$n_observations, length(x$states)))
  cat(sprintf("  Selected order  BIC: %d   AIC: %d   permutation-LRT: %d\n\n",
              x$bic_order, x$aic_order, x$optimal_order))
  tt <- x$test_table
  tt$loglik <- round(tt$loglik, 2)
  tt$AIC    <- round(tt$AIC,    2)
  tt$BIC    <- round(tt$BIC,    2)
  tt$g2     <- round(tt$g2,     2)
  print(tt, row.names = FALSE)
  invisible(x)
}

#' @rdname markov_order_test
#' @export
print.net_markov_order_group <- function(x, ...) {
  cat(sprintf("Markov Order Test -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}


#' @rdname markov_order_test
#' @export
summary.net_markov_order <- function(object, ...) {
  out <- object$test_table
  attr(out, "optimal_order") <- object$optimal_order
  attr(out, "bic_order")     <- object$bic_order
  attr(out, "aic_order")     <- object$aic_order
  attr(out, "alpha")         <- object$alpha
  attr(out, "n_perm")        <- object$n_perm
  out
}


#' @rdname markov_order_test
#' @section Plot panels:
#' \describe{
#'   \item{Fit across orders (`"ic"`)}{The log-likelihood, AIC and BIC of
#'     each order, one facet each, with a dashed line at the order the
#'     permutation test selects and a dotted line at the order with the
#'     lowest BIC.}
#'   \item{Permutation nulls (`"permutation"`)}{The density of the
#'     permutation null of \eqn{G^2} for each order, one facet each, with
#'     the observed \eqn{G^2} as a vertical line labelled with its p-value
#'     and coloured by whether the test rejects at `alpha`.}
#' }
#' The lines and labels use the Okabe-Ito palette; the null densities are
#' grey.
#' @export
plot.net_markov_order <- function(x,
                                   panel = c("both", "ic", "permutation"),
                                   combined = TRUE,
                                   ...) {
  panel <- match.arg(panel)
  .assert(is.logical(combined), length(combined) == 1L, !is.na(combined))
  pal <- c(reject = "#D55E00", accept = "#0072B2",
           ic_ll = "#009E73", ic_aic = "#E69F00", ic_bic = "#56B4E9",
           selected = "#D55E00")

  tt <- x$test_table
  ic_df <- rbind(
    data.frame(order = tt$order, metric = "log-likelihood", value = tt$loglik),
    data.frame(order = tt$order, metric = "AIC",            value = tt$AIC),
    data.frame(order = tt$order, metric = "BIC",            value = tt$BIC)
  )
  ic_df$metric <- factor(ic_df$metric,
                         levels = c("log-likelihood", "AIC", "BIC"))

  p_ic <- ggplot2::ggplot(ic_df,
      ggplot2::aes(x = .data$order, y = .data$value,
                   color = .data$metric)) +
    ggplot2::geom_vline(xintercept = x$optimal_order,
                        linetype = "dashed", color = pal["selected"],
                        linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = x$bic_order,
                        linetype = "dotted", color = pal["ic_bic"],
                        linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = unname(pal[c("ic_ll",
                                                       "ic_aic",
                                                       "ic_bic")]),
                                name = NULL) +
    ggplot2::scale_x_continuous(breaks = tt$order) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Markov order",
                  y = NULL,
                  title = "Model fit across orders",
                  subtitle = sprintf(
                    "Permutation LRT selects order %d (alpha = %.2f); BIC minimized at %d",
                    x$optimal_order, x$alpha, x$bic_order)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  # --- Permutation null ---
  null_df <- do.call(rbind, lapply(seq_along(x$permutation_null), function(k) {
    nn <- x$permutation_null[[k]]
    if (length(nn) == 0L) return(NULL)
    data.frame(order = k, stat = nn, stringsAsFactors = FALSE)
  }))
  obs_df <- data.frame(
    order = seq_len(x$max_order),
    g2    = tt$g2[-1L],
    p     = tt$p_permutation[-1L],
    sig   = tt$significant[-1L]
  )
  obs_df$decision <- ifelse(obs_df$sig, "reject H0", "accept H0")
  if (!is.null(null_df)) {
    null_df$order_fac <- factor(sprintf("order %d vs %d", null_df$order,
                                         null_df$order - 1L),
                                 levels = sprintf("order %d vs %d",
                                                   seq_len(x$max_order),
                                                   seq_len(x$max_order) - 1L))
  }
  obs_df$order_fac <- factor(sprintf("order %d vs %d", obs_df$order,
                                      obs_df$order - 1L),
                              levels = sprintf("order %d vs %d",
                                                seq_len(x$max_order),
                                                seq_len(x$max_order) - 1L))

  p_perm <- ggplot2::ggplot(null_df,
      ggplot2::aes(x = .data$stat)) +
    ggplot2::geom_density(fill = "#CCCCCC", color = "#888888",
                          alpha = 0.6) +
    ggplot2::geom_vline(data = obs_df,
                        ggplot2::aes(xintercept = .data$g2,
                                     color = .data$decision),
                        linewidth = 0.9) +
    ggplot2::geom_text(data = obs_df,
                       ggplot2::aes(x = .data$g2, y = 0,
                                    label = sprintf("p = %.3f", .data$p),
                                    color = .data$decision),
                       hjust = -0.1, vjust = -0.8, size = 3.2,
                       show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("reject H0" = unname(pal["reject"]),
                                            "accept H0" = unname(pal["accept"])),
                                name = NULL) +
    ggplot2::facet_wrap(~ .data$order_fac, scales = "free", ncol = 1) +
    ggplot2::labs(x = expression(G^2 ~ "statistic"),
                  y = "density",
                  title = "Permutation null distributions",
                  subtitle = sprintf(
                    "Within-context permutation, n_perm = %d; observed G^2 as vertical line",
                    x$n_perm)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  if (panel == "ic") return(p_ic)
  if (panel == "permutation") return(p_perm)

  if (!combined) return(invisible(list(ic = p_ic, permutation = p_perm)))

  # Side-by-side. When gridExtra is available, return its arranged gtable (a
  # drawable/saveable grob -- `ggsave(plot(fit))` keeps working). Otherwise
  # fall back to base `grid` viewports (always available): both panels still
  # render, but only a named list of the two ggplots can be returned.
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(p_ic, p_perm, ncol = 2))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, 2L)))
  print(p_ic,   vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  print(p_perm, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))
  grid::popViewport()
  invisible(list(ic = p_ic, permutation = p_perm))
}
