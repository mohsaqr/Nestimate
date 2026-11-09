# ---- Network Comparison Test (NCT) ----

#' @noRd
.nct_bad_input <- function(msg) {
  stop(errorCondition(msg, class = "nestimate_bad_input", call = NULL))
}

#' @noRd
.nct_degenerate <- function(msg) {
  stop(errorCondition(msg, class = "nestimate_degenerate_sample", call = NULL))
}

#' @noRd
.nct_validate <- function(iter, gamma, paired, abs, weighted) {
  flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  ok <- c(
    "`iter` must be a single whole number >= 1" =
      is.numeric(iter) && length(iter) == 1L && is.finite(iter) && iter >= 1 &&
      base::abs(iter - round(iter)) < sqrt(.Machine$double.eps),
    "`gamma` must be a single non-negative number" =
      is.numeric(gamma) && length(gamma) == 1L && is.finite(gamma) && gamma >= 0,
    "`paired` must be TRUE or FALSE"   = flag(paired),
    "`abs` must be TRUE or FALSE"      = flag(abs),
    "`weighted` must be TRUE or FALSE" = flag(weighted)
  )
  if (!all(ok)) .nct_bad_input(paste(names(ok)[!ok], collapse = "; "))
  invisible(TRUE)
}

# Keep the numeric columns of each group and check that the two groups hold the
# same variables in the same order. A position-wise comparison of differently
# ordered columns would test the wrong pairs, so a mismatch is an error, not a
# silent reorder. Unnamed input is taken by position and labelled V1..Vp.
#' @noRd
.nct_columns <- function(data1, data2) {
  numeric_only <- function(d) {
    # an unnamed matrix is matched by position: as.data.frame() would give
    # it V1..Vp names, so remember that it had none
    unnamed <- is.matrix(d) && is.null(colnames(d))
    d <- as.data.frame(d, stringsAsFactors = FALSE)
    m <- as.matrix(d[, vapply(d, is.numeric, logical(1)), drop = FALSE])
    if (unnamed) colnames(m) <- NULL
    m
  }
  d1 <- numeric_only(data1)
  d2 <- numeric_only(data2)
  if (any(is.infinite(d1)) || any(is.infinite(d2)))
    .nct_bad_input("`data1` and `data2` must not contain Inf or -Inf; use finite observations or NA for missing data.")
  if (ncol(d1) < 2L || ncol(d2) < 2L)
    .nct_bad_input("`data1` and `data2` each need at least 2 numeric variables.")
  n1 <- colnames(d1)
  n2 <- colnames(d2)
  if (!is.null(n1) && !is.null(n2) && !identical(n1, n2)) {
    listing <- function(v) if (length(v)) paste(v, collapse = ", ") else "none"
    .nct_bad_input(if (setequal(n1, n2) && length(n1) == length(n2))
      "`data1` and `data2` must have the same columns in the same order; put the columns of `data2` in the order of `data1`."
    else sprintf(
      "`data1` and `data2` must have the same columns. Only in data1: %s. Only in data2: %s.",
      listing(setdiff(n1, n2)), listing(setdiff(n2, n1))))
  }
  if (ncol(d1) != ncol(d2))
    .nct_bad_input(sprintf(
      "`data1` has %d numeric variables but `data2` has %d.", ncol(d1), ncol(d2)))
  labels <- n1 %||% n2 %||% paste0("V", seq_len(ncol(d1)))
  colnames(d1) <- colnames(d2) <- labels
  list(data1 = d1, data2 = d2, labels = labels)
}

# Correlation matrix of one (sub)sample, or NULL when a variable pair is never
# co-observed (fewer than 2 shared rows): that correlation is undefined and is
# refused. A variable without variation (over the rows where a pair is observed)
# enters as an isolated node -- its correlations are 0 -- the convention of
# psychnets' estimators for a constant column. Missing values are never filled:
# with NAs the correlation is pairwise. stats::cor() announces a constant
# variable with a "standard deviation is zero" warning; that one warning
# (matched in the session's own translation) is muffled because the case is
# handled here, and any other warning passes through.
#' @noRd
.nct_cor <- function(x) {
  pairwise <- anyNA(x)
  if (pairwise) {
    co <- crossprod(!is.na(x))
    if (any(co[upper.tri(co)] < 2L)) return(NULL)
  }
  zero_sd <- gettext("the standard deviation is zero", domain = "R-stats")
  S <- withCallingHandlers(
    stats::cor(x, use = if (pairwise) "pairwise.complete.obs" else "everything"),
    warning = function(w) {
      if (identical(conditionMessage(w), zero_sd)) invokeRestart("muffleWarning")
    }
  )
  if (anyNA(S)) {
    S[is.na(S)] <- 0
    diag(S) <- 1
  }
  S
}

# EBIC-glasso network of one (sub)sample, or NULL when a correlation is
# undefined. Correlation -> Matrix::nearPD symmetrization (NCT's
# NCT_estimator_GGM default) -> psychnets' EBIC glasso. With missing values the
# correlation is pairwise and n is the median pairwise co-observation count,
# the convention of psychnets' estimators.
#' @noRd
.nct_estimate <- function(x, gamma) {
  cor_x <- .nct_cor(x)
  if (is.null(cor_x)) return(NULL)
  n <- if (anyNA(x)) {
    co <- crossprod(!is.na(x))
    max(round(stats::median(co[upper.tri(co)])), 2L)
  } else nrow(x)
  cor_x <- as.matrix(Matrix::nearPD(cor_x, corr = TRUE)$mat)
  cor_x <- (cor_x + t(cor_x)) / 2
  w <- psychnets::ebic_glasso(cor_matrix = cor_x, n = n, gamma = gamma,
                              nlambda = 100L, lambda_min_ratio = 0.01,
                              native = TRUE)$weights
  # psychnets < 0.6.0 can return an absent edge as a ~1e-17 residue (at
  # lambda = max|S_ij|). Against an exact zero that is a difference, so two
  # empty networks would test as unequal. Same rule as psychnets >= 0.6.0:
  # |w| <= 1e-10 is no edge (`+ 0` keeps zeros unsigned).
  w[abs(w) <= 1e-10] <- 0
  w + 0
}

# Display name of a group: the argument's name when it is a plain variable
# (as t.test() reports "data: x and y"), otherwise "group 1" / "group 2".
#' @noRd
.nct_label <- function(expr, default) {
  if (is.name(expr)) as.character(expr) else default
}

# Names of the variables with fewer than two distinct observed values.
#' @noRd
.nct_constant_vars <- function(x) {
  const <- vapply(seq_len(ncol(x)), function(j) {
    v <- x[!is.na(x[, j]), j]
    length(unique(v)) < 2L
  }, logical(1))
  colnames(x)[const]
}

#' Network Comparison Test
#'
#' @description
#' `nct()` tests whether two networks estimated from two samples of the same
#' variables differ, at the three levels of van Borkulo et al. (2022):
#' **network structure** (statistic M, the largest absolute edge
#' difference), **global strength** (statistic S, the difference in the
#' summed absolute edge weights, or with `abs = FALSE` in the summed signed
#' weights, the global expected influence) and **individual edges**
#' (statistic E, the absolute difference of each edge). Each network is an
#' EBIC graphical lasso network of the Pearson correlations, and each
#' p-value comes from re-estimating both networks on permuted group labels
#' (within-pair swaps with `paired = TRUE`).
#'
#' A p-value counts the permutations whose statistic is at least the
#' observed one, with the observed split added:
#' `(sum(perm >= observed) + 1) / (iter + 1)`, so the smallest attainable
#' value is `1 / (iter + 1)`. A non-significant test does not
#' show that the networks are equal; with small samples the test has little
#' power (van Borkulo et al. 2022). The edge tests are many tests at once;
#' set `p_adjust` before reading single edges.
#' Both networks are regularised estimates, so check how stable each one is
#' with [boot_glasso()].
#'
#' @details
#' `nct()` follows `NetworkComparisonTest::NCT()` with its defaults
#' `abs = TRUE`, `weighted = TRUE` and `paired = FALSE`. The network
#' estimator is EBIC-selected glasso (`psychnets::ebic_glasso()`) applied to
#' a Pearson correlation matrix made positive definite with
#' `Matrix::nearPD()` (NCT's `NCT_estimator_GGM` default). Under the same
#' seed both functions permute the same rows, so the permutation p-values
#' agree; the glasso solver is not the Fortran one NCT wraps, so the test
#' statistics agree to independent-solver precision (of the order of
#' `1e-4`) rather than bit for bit. The \pkg{Matrix} package (in Suggests)
#' must be installed.
#'
#' The two groups must hold the same variables in the same order;
#' non-numeric columns, such as an id or a group label, are dropped first.
#' Infinite values are refused; use finite observations or `NA` for missing
#' data.
#'
#' A variable without variation in a group, or in a permuted split (a rare
#' binary item can be constant in a split), enters that network as an
#' isolated node. The statistic is then defined on every split, so the
#' permutation test stays exact. A constant variable in an observed group
#' raises a warning of class `nestimate_constant_variable`.
#'
#' Malformed `iter`, `gamma`, `paired`, `abs`, `weighted` or data raise an
#' error of class `nestimate_bad_input`.
#'
#' @param data1 The observations of group 1: a numeric matrix or data frame,
#'   one row per observation. Non-numeric columns are dropped; at least two
#'   numeric variables and three rows are needed.
#' @param data2 The observations of group 2, with the same variables as
#'   `data1` in the same order; columns that differ in name or order raise an
#'   error of class `nestimate_bad_input`. An unnamed matrix is matched by
#'   position and takes the other group's names, or `V1`, `V2`, ... when
#'   neither group has names.
#' @param iter The number of permutations, a whole number of at least 1.
#'   Default `1000`.
#' @param gamma The EBIC tuning parameter of the glasso, a non-negative
#'   number; larger values give sparser networks. Default `0.5`.
#' @param paired Whether row `i` of both groups is the same unit (such as a
#'   person measured twice). Each permutation then swaps the two rows of a
#'   unit at random instead of shuffling all rows, and the groups need the
#'   same number of rows, otherwise an error of class `nestimate_bad_input`
#'   is raised. Default `FALSE`.
#' @param abs Whether global strength sums the absolute edge weights. With
#'   `FALSE`, it sums the signed weights (global expected influence). Default
#'   `TRUE`.
#' @param weighted Whether the tests compare the weighted networks. With
#'   `FALSE`, both networks are binarised (an edge is present or absent)
#'   before the statistics are computed. Default `TRUE`.
#' @param p_adjust The multiplicity correction for the per-edge p-values,
#'   any method of [stats::p.adjust()]. Default `"none"`.
#' @param na_method How missing values are handled. `"pairwise"` (default)
#'   keeps every row, permutes rows as usual, and estimates each sample or
#'   permuted split from its pairwise-complete correlation matrix, with the
#'   median pairwise co-observation count as the sample size; a variable
#'   pair sharing fewer than two observed rows there raises an error of class
#'   `nestimate_degenerate_sample`. `"listwise"` drops incomplete rows first
#'   (with `paired = TRUE`, a pair is dropped when either member is
#'   incomplete). Complete data are unaffected by the choice.
#' @return A `net_nct` object. `print()` reports the three tests, each
#'   group's global strength and the edges, largest difference first, and
#'   returns the object invisibly; [summary()] returns the M and S tests as
#'   a data frame with one row per test. The object's fields:
#' \describe{
#'   \item{nw1, nw2}{The estimated weighted adjacency matrices of the two
#'     groups (binarised with `weighted = FALSE`).}
#'   \item{structure}{The network-structure test (statistic M, the largest
#'     absolute edge difference): `observed`, `perm` (the statistic in each
#'     permutation) and `p_value`.}
#'   \item{strength}{The global-strength test (statistic S): the same, plus
#'     `by_group`, each network's global strength (NCT's `glstrinv.sep`).}
#'   \item{E}{The per-edge tests: `observed` and `p_value` hold one value per
#'     node pair (the upper triangle), `perm` is an `iter` by node-pair
#'     matrix, and `edge_names` is a data frame with the columns `from` and
#'     `to` naming each node pair.}
#'   \item{n_iter}{The number of permutations.}
#'   \item{groups}{The display names of the two groups: the argument names
#'     when they are plain variables, otherwise `"group 1"` and
#'     `"group 2"`.}
#'   \item{n}{The rows used in each group.}
#'   \item{paired}{Whether a paired test was used.}
#'   \item{params}{The settings used: `gamma`, `abs`, `weighted`, `p_adjust`
#'     and `na_method`.}
#' }
#' The letters follow van Borkulo et al. (2022) and
#' `NetworkComparisonTest`. Before Nestimate 0.9.27 the fields `M` and `S`
#' held the global-strength and network-structure tests the other way
#' round; reading them now raises an error of class `nestimate_renamed`.
#' @references van Borkulo, C. D., van Bork, R., Boschloo, L., Kossakowski,
#'   J. J., Tio, P., Schoevers, R. A., Borsboom, D., & Waldorp, L. J. (2022).
#'   Comparing network structures on three aspects: A permutation test.
#'   \emph{Psychological Methods}. \doi{10.1037/met0000476}
#' @examples
#' set.seed(1)
#' x1 <- matrix(rnorm(100 * 4), 100, 4)
#' x2 <- matrix(rnorm(100 * 4), 100, 4)
#' colnames(x1) <- colnames(x2) <- paste0("V", 1:4)
#' # iter = 20 keeps the example fast; a real analysis uses 1000 or more.
#' res <- nct(x1, x2, iter = 20)
#' res
#' summary(res)
#' @export
nct <- function(data1, data2, iter = 1000L, gamma = 0.5,
                paired = FALSE, abs = TRUE, weighted = TRUE,
                p_adjust = "none", na_method = c("pairwise", "listwise")) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("nct() requires the 'Matrix' package.", call. = FALSE)
  }
  na_method <- match.arg(na_method)
  groups <- c(.nct_label(substitute(data1), "group 1"),
              .nct_label(substitute(data2), "group 2"))
  .nct_validate(iter, gamma, paired, abs, weighted)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  cols <- .nct_columns(data1, data2)
  data1 <- cols$data1
  data2 <- cols$data2
  labs <- cols$labels

  if (paired && nrow(data1) != nrow(data2)) {
    .nct_bad_input("A paired comparison requires the two groups to have equal size.")
  }
  if (na_method == "listwise") {
    keep1 <- stats::complete.cases(data1)
    keep2 <- stats::complete.cases(data2)
    if (paired) keep1 <- keep2 <- keep1 & keep2
    data1 <- data1[keep1, , drop = FALSE]
    data2 <- data2[keep2, , drop = FALSE]
  }
  iter <- as.integer(round(iter))
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  if (n1 < 3L || n2 < 3L) {
    .nct_bad_input(sprintf("Each group needs at least 3 rows; got %d and %d.",
                           n1, n2))
  }
  dataall <- rbind(data1, data2)

  undefined <- function(where) .nct_degenerate(sprintf(
    "A variable pair shares fewer than 2 observed rows in %s, so its correlation is undefined. Use na_method = \"listwise\" or drop the sparsely observed variable.",
    where))
  est <- function(x, where) {
    nw <- .nct_estimate(x, gamma)
    if (is.null(nw)) undefined(where)
    nw
  }

  constant <- list(data1 = .nct_constant_vars(data1),
                   data2 = .nct_constant_vars(data2))
  constant <- constant[lengths(constant) > 0L]
  if (length(constant)) {
    warning(warningCondition(sprintf(
      "%s; a constant variable enters that network as an isolated node.",
      paste(sprintf("constant in `%s`: %s", names(constant),
                    vapply(constant, paste, character(1), collapse = ", ")),
            collapse = "; ")),
      class = "nestimate_constant_variable", call = NULL))
  }

  nw1 <- est(data1, "`data1`")
  nw2 <- est(data2, "`data2`")
  ut <- upper.tri(nw1)

  binarize <- function(m) (m != 0) * 1
  if (!weighted) {
    nw1 <- binarize(nw1)
    nw2 <- binarize(nw2)
  }

  # ---- observed statistics ----
  # NCT notation: M = network structure (largest absolute edge difference),
  # S = global strength (summed absolute weights; signed when !abs), E = edges.
  global <- function(m) if (abs) sum(base::abs(m[ut])) else sum(m[ut])
  strength_by_group <- c(global(nw1), global(nw2))
  str_obs <- base::abs(strength_by_group[1L] - strength_by_group[2L])
  diff_real <- base::abs(nw1 - nw2)
  max_obs <- max(diff_real[ut])
  E_obs <- diff_real[ut]
  n_edges <- length(E_obs)

  # ---- permutation ----
  # A for loop: each permutation refits two networks; vectorising would
  # materialise iter x n_edges intermediate networks.
  str_perm <- numeric(iter)
  max_perm <- numeric(iter)
  E_perm <- matrix(0, nrow = iter, ncol = n_edges)

  for (i in seq_len(iter)) {
    if (paired) {
      s <- sample(c(1L, 2L), n1, replace = TRUE)
      x1p <- rbind(data1[s == 1L, , drop = FALSE],
                   data2[s == 2L, , drop = FALSE])
      x2p <- rbind(data2[s == 1L, , drop = FALSE],
                   data1[s == 2L, , drop = FALSE])
    } else {
      s <- sample(seq_len(n1 + n2), n1, replace = FALSE)
      x1p <- dataall[s, , drop = FALSE]
      x2p <- dataall[-s, , drop = FALSE]
    }
    where <- sprintf("permuted split %d", i)
    r1 <- est(x1p, where)
    r2 <- est(x2p, where)
    if (!weighted) {
      r1 <- binarize(r1)
      r2 <- binarize(r2)
    }
    str_perm[i] <- base::abs(global(r1) - global(r2))
    diff_perm <- base::abs(r1 - r2)
    max_perm[i] <- max(diff_perm[ut])
    E_perm[i, ] <- diff_perm[ut]
  }

  # ---- p-values (matches NCT formula) ----
  str_pval <- (sum(str_perm >= str_obs) + 1) / (iter + 1)
  max_pval <- (sum(max_perm >= max_obs) + 1) / (iter + 1)
  E_pval_raw <- (colSums(E_perm >= matrix(E_obs, iter, n_edges, byrow = TRUE))
                 + 1) / (iter + 1)
  E_pval <- if (p_adjust != "none") {
    stats::p.adjust(E_pval_raw, method = p_adjust)
  } else {
    E_pval_raw
  }

  # node pairs in the upper-triangle (column-major) order of the E vectors
  ij <- which(ut, arr.ind = TRUE)
  edge_names <- data.frame(from = labs[ij[, 1L]], to = labs[ij[, 2L]],
                           stringsAsFactors = FALSE, row.names = NULL)

  structure(
    list(
      nw1       = nw1,
      nw2       = nw2,
      structure = list(observed = max_obs, perm = max_perm, p_value = max_pval),
      strength  = list(observed = str_obs, perm = str_perm, p_value = str_pval,
                       by_group = strength_by_group),
      E         = list(observed = E_obs, perm = E_perm,
                       p_value = E_pval, edge_names = edge_names),
      groups    = groups,
      n_iter    = iter,
      n         = c(n1, n2),
      paired    = paired,
      params    = list(gamma = gamma, abs = abs, weighted = weighted,
                       p_adjust = p_adjust, na_method = na_method)
    ),
    class = "net_nct"
  )
}


#' @rdname nct
#' @param x For the `print()` and `as.data.frame()` methods: a `net_nct`
#'   object.
#' @param ... For the `print()`, `summary()` and `as.data.frame()` methods:
#'   ignored.
#' @param object For the `summary()` method: a `net_nct` object.
#' @param n For the `print()` method: the number of edges to list, largest
#'   difference first. Default `Inf` (all).
#' @return For the `summary()` method: a data frame with one row per global
#'   test and the columns `test` (`"network structure"`, and
#'   `"global strength"`, or `"global expected influence"` when
#'   `abs = FALSE`), `statistic` (`"M"`, `"S"`), `observed`, `p_value` and
#'   `significant` (`p_value < 0.05`, a fixed level). The `as.data.frame()`
#'   method gives one row per node pair, with the columns `from`, `to`,
#'   `weight_1` and `weight_2` (the edge in each group's network),
#'   `difference` (statistic E), `p_value` and `significant`
#'   (`p_value < 0.05`).
#' @export
print.net_nct <- function(x, n = Inf, ...) {
  cat(.nct_display(x, n = n), sep = "\n")
  invisible(x)
}

# The printed report, as lines (same layout as psychnets::net_compare()).
# Fixed-width columns; numbers to 3 decimals; p < 0.001 shown as "<0.001";
# "*" marks p < 0.05.
#' @noRd
.nct_display <- function(x, n = Inf, alpha = 0.05) {
  g <- x$groups %||% c("group 1", "group 2")
  gl <- summary.net_nct(x)
  ed <- as.data.frame.net_nct(x)
  num <- function(v) formatC(v, format = "f", digits = 3)
  pval <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
  star <- function(p) ifelse(p < alpha, "*", "")
  pad <- function(v, w, right = TRUE) formatC(v, width = if (right) w else -w)
  adj <- x$params$p_adjust %||% "none"

  header <- c(
    sprintf("Network Comparison Test: %s vs %s", g[1L], g[2L]),
    sprintf("%d permutations, %s groups, n = %d and %d",
            x$n_iter, if (isTRUE(x$paired)) "paired" else "independent",
            x$n[1L], x$n[2L]),
    "")
  n_sig <- sum(ed$p_value < alpha)
  labels <- c("M  network structure", paste("S ", gl$test[2L]), "E  edges differing")
  wl <- max(nchar(labels))
  test_lines <- c(
    paste0("  ", pad("Test", wl, FALSE), "  ", pad("Statistic", 9), "  ", pad("p", 6)),
    paste0("  ", pad(labels[1:2], wl, FALSE), "  ", pad(num(gl$observed), 9), "  ",
           pad(pval(gl$p_value), 6), " ", star(gl$p_value)),
    sprintf("       %s %s vs %s %s", g[1L], num(x$strength$by_group[1L]),
            g[2L], num(x$strength$by_group[2L])),
    paste0("  ", pad(labels[3L], wl, FALSE), "  ",
           pad(sprintf("%d of %d", n_sig, nrow(ed)), 9)))

  k <- min(nrow(ed), n)
  edge_lines <- character(0)
  if (k > 0) {
    top <- ed[order(-ed$difference, ed$from, ed$to)[seq_len(k)], ]
    edge <- paste(top$from, "-", top$to)
    gn <- substr(g, 1L, 12L)
    we <- max(nchar(c("Edge", edge)))
    wg <- pmax(nchar(gn), 6L)
    edge_lines <- c(
      "",
      if (k < nrow(ed)) sprintf("  Largest edge differences (%d of %d)", k, nrow(ed))
      else "  Edges, largest difference first",
      paste0("  ", pad("Edge", we, FALSE), "  ", pad(gn[1L], wg[1L]), "  ",
             pad(gn[2L], wg[2L]), "  ", pad("|diff|", 6), "  ", pad("p", 6)),
      paste0("  ", pad(edge, we, FALSE), "  ", pad(num(top$weight_1), wg[1L]),
             "  ", pad(num(top$weight_2), wg[2L]), "  ",
             pad(num(top$difference), 6), "  ", pad(pval(top$p_value), 6), " ",
             star(top$p_value)))
  }
  footer <- c("", sprintf("  * p < %g%s", alpha,
                          if (identical(adj, "none")) ""
                          else sprintf("; edge p-values %s-adjusted", adj)))
  c(header, test_lines, edge_lines, footer)
}


#' @rdname nct
#' @export
summary.net_nct <- function(object, ...) {
  p <- c(object$structure$p_value, object$strength$p_value)
  data.frame(
    test = c("network structure",
             if (isFALSE(object$params$abs)) "global expected influence"
             else "global strength"),
    statistic = c("M", "S"),
    observed = c(object$structure$observed, object$strength$observed),
    p_value = p,
    significant = p < 0.05,
    stringsAsFactors = FALSE)
}


#' @rdname nct
#' @param row.names,optional For the `as.data.frame()` method: ignored, kept
#'   for compatibility with the base generic.
#' @export
as.data.frame.net_nct <- function(x, row.names = NULL, optional = FALSE, ...) {
  ut <- upper.tri(x$nw1)
  data.frame(
    from       = as.character(x$E$edge_names$from),
    to         = as.character(x$E$edge_names$to),
    weight_1   = x$nw1[ut],
    weight_2   = x$nw2[ut],
    difference = as.numeric(x$E$observed),
    p_value    = as.numeric(x$E$p_value),
    significant = x$E$p_value < 0.05,
    stringsAsFactors = FALSE,
    row.names  = NULL)
}


#' @noRd
#' @export
`$.net_nct` <- function(x, name) {
  if (name %in% c("M", "S")) {
    stop(errorCondition(sprintf(paste0(
      "`$%s` was removed in Nestimate 0.9.27: it held the %s test, which ",
      "van Borkulo et al. (2022) and NetworkComparisonTest call %s. Use ",
      "summary(x), or `$structure` (M) / `$strength` (S)."),
      name, if (name == "M") "global-strength" else "network-structure",
      if (name == "M") "S" else "M"),
      class = "nestimate_renamed", call = NULL))
  }
  .subset2(x, name)
}
