# ---- Permutation diagnostics for nested data ----

#' Does Nesting Bias a Permutation Test?
#'
#' @description
#' `permutation_diagnostics()` shows what treating nested sequences as
#' independent costs a network comparison. It runs [permutation()] twice on
#' the same networks and seed: once shuffling single sequences, once with
#' `actor`, shuffling whole actors (the persons whose sessions the
#' sequences are, or the teams of students). It returns the two results
#' side by side, one row per compared pair or per edge, with the
#' intraclass correlation (ICC) and design effect that explain any
#' difference between them. It takes two transition networks, or every
#' pair of groups of a `netobject_group`.
#'
#' A design effect near 1 and few changed edges say that the nesting
#' changes the comparison little. Where the two tests differ, the
#' actor-level one is the test that follows the nesting, as [permutation()]
#' sets out for nested data. The sections *Nested data and actor* and *ICC
#' and design effect* of [permutation()] define these quantities, and *How
#' to read the result* below says what each column shows.
#'
#' @section How to read the result:
#' \describe{
#'   \item{`deff_edges`, `deff_global`}{The design effect
#'     (Kish, 1965): actor-level over sequence-level null variance. 1 means
#'     the two shuffles give the same chance variation; above 1 the
#'     actor-level one varies more, below 1 less.}
#'   \item{`edges_changed`}{Edges significant under one test but not
#'     the other. Edges with p-values close to `alpha` can flip from
#'     Monte Carlo error alone; increase `iter` before reading much
#'     into one or two.}
#'   \item{`min_p_actor` above `alpha`}{Too few actors: the
#'     actor-level test cannot reject anything.}
#' }
#'
#' @param x The networks to compare: a `netobject_group`, whose every pair
#'   of groups is diagnosed, or a `netobject`, compared with `y`. An `mcml`,
#'   whose layers hold different states, raises an error of class
#'   `nestimate_bad_input`. Only
#'   the transition methods
#'   `"relative"`, `"frequency"` and `"co_occurrence"` are supported; any
#'   other method, association networks and `wtna` among them, raises an
#'   error of class `nestimate_actor_unsupported`.
#' @param y The `netobject` to compare with when `x` is a `netobject`.
#'   Default `NULL`, as it must stay when `x` is a group.
#' @param actor The column identifying the actor each sequence belongs to,
#'   as in [permutation()], looked up in the networks' metadata or sequence
#'   data. A column not found, or one with missing values, raises an error
#'   of class `nestimate_actor_missing`. Required.
#' @param iter The number of permutations for each of the two tests, at
#'   least 2. Default `1000`.
#' @param alpha The significance level that counts an edge as significant,
#'   between 0 and 1. The p-values are not adjusted for multiplicity.
#'   Default `0.05`.
#' @param level `"overall"` (default) returns one row per compared pair;
#'   `"edges"` one row per edge per pair.
#' @param seed A seed for the shuffles, used by both tests, so the result
#'   repeats exactly. Default `NULL`.
#'
#' @return A data frame. With `level = "overall"`, one row per compared
#' pair:
#' \describe{
#'   \item{pair}{`"<x> vs <y>"`, the names of the two groups (`"x vs y"`
#'     for two networks).}
#'   \item{n_sequences, n_actors}{The sequences and distinct actors in the
#'     pair.}
#'   \item{design}{`"between"` (every actor in one group), `"within"`
#'     (every actor in both groups) or `"mixed"`.}
#'   \item{icc, icc_ci_lower, icc_ci_upper}{How alike the sequences of one
#'     actor are: the jackknife bias-corrected ICC, as printed by
#'     [permutation()], with a 95% interval; with fewer than 3 actors the
#'     ICC is the uncorrected estimate and the interval is `NA`; the ICC is
#'     `NA` when no transition share varies.}
#'   \item{deff_edges}{The median over edges of the design effect, the
#'     ratio of the actor-level to the sequence-level null variance of the
#'     edge difference (Kish, 1965).}
#'   \item{deff_global}{The same ratio for the global `M` statistic.}
#'   \item{p_global_sequence, p_global_actor}{The permutation p-values of
#'     `M` when sequences or whole actors are reassigned.}
#'   \item{sig_edges_sequence, sig_edges_actor}{The number of edges with
#'     `p < alpha` under each test.}
#'   \item{edges_changed}{The number of edges significant under one test
#'     but not the other.}
#'   \item{min_p_actor}{The smallest p-value the actor-level test can
#'     reach: one over the number of distinct arrangements of the actors,
#'     or `1 / (iter + 1)`, whichever is larger.}
#' }
#' With `level = "edges"`, one row per edge present in either network, with
#' the columns `pair`, `from`, `to`, `diff` (the observed difference,
#' `x - y`), `icc` (the edge's ANOVA ICC, not bias-corrected; `NA` where
#' the share does not vary), `null_sd_sequence` and `null_sd_actor` (the
#' SD of the shuffled differences under each test), `deff` (`NaN` where
#' the edge difference never varies under either null), `p_sequence`,
#' `p_actor` and `changed` (whether the edge is significant under one test
#' only).
#'
#' The ICC and design effects are those of the actor-level run, as
#' [permutation()] prints them when given `actor`; the p-values and
#' significance counts compare it with a separate ordinary run.
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall. (jackknife, ch. 11)
#'
#' Kish, L. (1965). \emph{Survey Sampling}. Wiley. (design effect)
#'
#' Anderson, M. J., & ter Braak, C. J. F. (2003). Permutation tests for
#' multi-factorial analysis of variance. \emph{Journal of Statistical
#' Computation and Simulation}, 73(2), 85-113.
#'
#' @seealso \code{\link{permutation}}
#'
#' @examples
#' # Students are nested in teams; Achiever is a team-level label.
#' # iter = 50 keeps the example fast; a real analysis uses 1000 or more.
#' net <- build_network(group_regulation_long, method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time",
#'                      group = "Achiever")
#' permutation_diagnostics(net, actor = "Group", iter = 50, seed = 1)
#' head(permutation_diagnostics(net, actor = "Group", iter = 50,
#'                              level = "edges", seed = 1))
#' @export
permutation_diagnostics <- function(x, y = NULL, actor, iter = 1000L,
                                    alpha = 0.05,
                                    level = c("overall", "edges"),
                                    seed = NULL) {
  level <- match.arg(level)
  .assert(
    "`actor` must be a single column name" =
      is.character(actor) && length(actor) == 1L && !is.na(actor) && nzchar(actor),
    "`iter` must be a single integer >= 2" =
      is.numeric(iter) && length(iter) == 1L && iter >= 2,
    "`alpha` must be in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha < 1
  )
  iter <- as.integer(iter)

  if (inherits(x, "mcml") || inherits(y, "mcml")) {
    .stop_mcml_input("permutation_diagnostics()",
                     "compares networks over the same states")
  }
  pairs <- if (inherits(x, "netobject_group") && is.null(y)) {
    if (length(x) < 2L) {
      stop(errorCondition("Need at least 2 groups to diagnose a comparison.",
        class = "nestimate_bad_input", call = NULL))
    }
    idx <- utils::combn(length(x), 2L)
    lapply(seq_len(ncol(idx)), function(k) {
      list(x = x[[idx[1L, k]]], y = x[[idx[2L, k]]],
           name = paste(names(x)[idx[1L, k]], "vs", names(x)[idx[2L, k]]))
    })
  } else {
    .assert("`x` and `y` must be netobjects" =
                inherits(x, "netobject") && inherits(y, "netobject"))
    list(list(x = x, y = y, name = "x vs y"))
  }

  rows <- lapply(pairs, function(pr) {
    .permutation_diagnose_pair(pr$x, pr$y, pr$name, actor, iter, alpha,
                               level, seed)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Diagnose one x-vs-y comparison
#'
#' ICC and design effects come from the actor-level run itself (its
#' sequence-level reference is drawn in the same pass), so they match what
#' print(permutation(..., actor = )) reports.
#' @noRd
.permutation_diagnose_pair <- function(x, y, pair, actor, iter, alpha,
                                       level, seed) {
  method <- .resolve_method_alias(x$method)
  if (!method %in% c("relative", "frequency", "co_occurrence")) {
    .stop_block_unsupported(actor, method)
  }
  by_sequence <- permutation(x, y, iter = iter, alpha = alpha, seed = seed)
  by_actor <- permutation(x, y, iter = iter, alpha = alpha, seed = seed,
                          actor = actor)
  clus <- by_actor$clustering
  edges <- by_actor$clustering_edges
  at <- cbind(edges$from, edges$to)
  p_sequence <- by_sequence$p_values[at]
  p_actor <- by_actor$p_values[at]
  sig_sequence <- p_sequence < alpha
  sig_actor <- p_actor < alpha

  if (identical(level, "edges")) {
    return(data.frame(
      pair = pair, from = edges$from, to = edges$to,
      diff = by_actor$diff[at], icc = edges$icc,
      null_sd_sequence = edges$null_sd_sequence,
      null_sd_actor = edges$null_sd_actor,
      deff = edges$deff, p_sequence = p_sequence, p_actor = p_actor,
      changed = sig_sequence != sig_actor,
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    pair = pair,
    n_sequences = clus$n_sequences, n_actors = clus$n_actors,
    design = clus$design,
    icc = clus$icc, icc_ci_lower = clus$icc_ci_lower,
    icc_ci_upper = clus$icc_ci_upper,
    deff_edges = clus$deff_edges, deff_global = clus$deff_global,
    p_global_sequence = by_sequence$global$p_value[1L],
    p_global_actor = by_actor$global$p_value[1L],
    sig_edges_sequence = sum(sig_sequence),
    sig_edges_actor = sum(sig_actor),
    edges_changed = sum(sig_sequence != sig_actor),
    min_p_actor = clus$min_p,
    stringsAsFactors = FALSE
  )
}


#' ICC text for print methods: value with interval, or why it is missing
#' @noRd
.format_icc <- function(icc, lower, upper) {
  if (is.na(icc)) return("ICC not estimable (no variation between actors)")
  if (is.na(lower) || is.na(upper)) {
    return(sprintf("ICC = %.3f (interval needs at least 3 actors)", icc))
  }
  sprintf("ICC = %.3f [95%% CI %.3f, %.3f]", icc, lower, upper)
}


#' One-way ANOVA ICC per column, unbalanced clusters
#'
#' Uses the unbalanced cluster size n0 = (N - sum(n_k^2) / N) / (k - 1).
#' Columns without variance return NA (the ICC is undefined there).
#' @noRd
.anova_icc_cols <- function(feat, cluster) {
  n_k <- as.vector(table(cluster))
  k <- length(n_k)
  n_total <- nrow(feat)
  if (k < 2L || n_total <= k) return(rep(NA_real_, ncol(feat)))
  n0 <- (n_total - sum(n_k^2) / n_total) / (k - 1)
  cl_mean <- rowsum(feat, cluster) / n_k
  grand <- colMeans(feat)
  msb <- colSums(n_k * sweep(cl_mean, 2L, grand)^2) / (k - 1)
  msw <- colSums((feat - cl_mean[as.character(cluster), , drop = FALSE])^2) /
    (n_total - k)
  denom <- msb + (n0 - 1) * msw
  ifelse(denom > 0, (msb - msw) / denom, NA_real_)
}


#' Frequency-weighted ICC of transition shares within groups, with a
#' leave-one-block-out jackknife interval
#'
#' The estimate is the jackknife bias-corrected ICC and the interval is
#' estimate +/- t(k - 1) x jackknife SE (Efron & Tibshirani, 1993, ch. 11).
#' `interval = FALSE` skips the jackknife and returns the raw estimate.
#' Block-bootstrap percentile and BCa intervals were tried and rejected:
#' resampling k blocks shrinks the between-block variance, and in
#' calibration (30-80 blocks) they covered 57-89% against 90-96% for this
#' interval, which errs low when it misses
#' (local_testing_and_equivalence/calibration-block-icc-coverage.R).
#' @noRd
.block_icc <- function(counts, block_ids, stratum, interval = TRUE) {
  totals <- rowSums(counts)
  keep <- totals > 0                      # sequences with no transition
  counts <- counts[keep, , drop = FALSE]
  shares <- counts / totals[keep]
  block_ids <- block_ids[keep]
  stratum <- stratum[keep]

  per_edge <- function(rows, cl) {
    per <- lapply(split(rows, stratum[rows]), function(r) {
      list(icc = .anova_icc_cols(shares[r, , drop = FALSE],
                                 cl[match(r, rows)]),
           n = length(unique(cl[match(r, rows)])))
    })
    vals <- matrix(vapply(per, function(p) p$icc, numeric(ncol(shares))),
                   nrow = ncol(shares))
    n_cl <- vapply(per, function(p) p$n, numeric(1))
    w <- (!is.na(vals)) * rep(n_cl, each = nrow(vals))
    ifelse(rowSums(w) > 0,
           rowSums(ifelse(is.na(vals), 0, vals) * w) / rowSums(w), NA_real_)
  }
  # edge weights come from the rows in hand, so the jackknife variance
  # includes the variability of the weighting itself
  weighted_icc <- function(rows) {
    v <- per_edge(rows, block_ids[rows])
    weight <- colSums(counts[rows, , drop = FALSE])
    ok <- !is.na(v) & weight > 0
    if (!any(ok)) return(NA_real_)       # no transition varies: undefined
    sum(v[ok] * weight[ok]) / sum(weight[ok])
  }

  all_rows <- seq_len(nrow(shares))
  icc_edges <- per_edge(all_rows, block_ids)
  estimate <- weighted_icc(all_rows)

  rows_by_block <- split(all_rows, block_ids)
  k <- length(rows_by_block)
  if (!interval || k < 3L || is.na(estimate)) {
    return(list(estimate = estimate, ci = c(NA_real_, NA_real_),
                per_edge = icc_edges))
  }
  jack <- vapply(seq_len(k), function(j) {
    weighted_icc(unlist(rows_by_block[-j], use.names = FALSE))
  }, numeric(1))
  corrected <- k * estimate - (k - 1) * mean(jack)
  se <- sqrt((k - 1) / k * sum((jack - mean(jack))^2))

  list(estimate = corrected,
       ci = corrected + c(-1, 1) * stats::qt(0.975, df = k - 1L) * se,
       per_edge = icc_edges)
}
