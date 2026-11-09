# ---- Sequence Pattern Comparison ----

utils::globalVariables(c(
  "freq", "p_label", "pattern", "resid", "resid_label",
  "side", "x_center", "xmax", "xmin", "yidx", "ymax", "ymin"
))

#' Compare Subsequence Patterns Between Groups
#'
#' @description
#' `sequence_compare()` compares how often short subsequences occur in two
#' or more groups of sequences. For each pattern length in `sub` it counts
#' every k-gram (a run of k consecutive states, such as `"A->B->C"`) in
#' each group, keeps the patterns that are frequent enough in every group,
#' and returns a `net_sequence_comparison` whose pattern table has one row
#' per pattern with its count, proportion and standardised residual in each
#' group and, optionally, a permutation or chi-square test of whether its
#' frequency differs between the groups. It takes a grouped network, a
#' single network or wide sequences with a grouping.
#'
#' A residual says how far a group's count of a pattern lies from the count
#' expected if pattern and group were independent, in approximately
#' standard-normal units (Haberman, 1973); a positive value means the
#' pattern is more frequent in that group than independence predicts. The
#' p-values are adjusted across every pattern tested. The permutation test
#' treats the sequence as the
#' unit, while the chi-square test treats every occurrence as independent,
#' which overstates the evidence when occurrences within a sequence depend
#' on each other (see *Details*). Patterns of different lengths overlap
#' (each 3-gram contains two 2-grams), so their tests are not independent of
#' each other. [sequence_patterns()] with `group` asks the related question
#' of whether a pattern is present in more sequences of one group.
#'
#' @details
#' The standardised residuals come from the 2 x G table of the pattern's
#' counts against the counts of all other k-grams of the same length in
#' each group, with the formula
#' \code{(o - e) / sqrt(e * (1 - r/N) * (1 - c/N))}. They are always
#' computed, whatever `test` is. A k-gram that spans a missing position is
#' not counted.
#'
#' `test` chooses the inference:
#' \describe{
#'   \item{\code{"permutation"}}{Shuffles the group labels across
#'     sequences `iter` times and recomputes, for each pattern, the
#'     Euclidean norm of its row of differences between observed and
#'     expected counts in the patterns-by-groups table of the patterns kept
#'     at that length (expected under independence). It asks whether the
#'     pattern's distribution is associated with group membership at the
#'     level of the sequence,
#'     which respects the sequence as the unit of analysis but has little
#'     power when there are few sequences. The p-value is
#'     `(b + 1) / (iter + 1)`, where `b` counts the permutations with a
#'     statistic at least as large as the observed one, so it cannot fall
#'     below `1 / (iter + 1)`. The effect size is the observed statistic
#'     minus the mean of the permuted ones, divided by their standard
#'     deviation (taken as 1 when the permuted statistics do not vary; with
#'     `iter = 1` the effect size is `NA`, since at least two permutations
#'     are needed). The shuffles use R's random number generator: call
#'     [base::set.seed()] first for reproducible p-values.}
#'   \item{\code{"chisq"}}{Runs [stats::chisq.test()] on each pattern's
#'     2 x G table (with Yates' continuity correction when there are two
#'     groups). It asks whether the groups' streams of k-grams produce the
#'     pattern at different rates, treating every occurrence as an
#'     independent event: it is fast and has power with few sequences, but
#'     the independence it assumes is optimistic when sequences are
#'     strongly autocorrelated.}
#'   \item{\code{"none"}}{No test; the counts, proportions and residuals
#'     only.}
#' }
#'
#' The p-values are adjusted once across all patterns of all lengths, with
#' [stats::p.adjust()].
#'
#' @param x The sequences: a `netobject_group` (from [build_network()] with
#'   `group`, or from a clustering), whose groups are compared; a
#'   `netobject`, `tna` model or `cograph_network` built from sequences; a
#'   wide data frame with one row per sequence and one column per position,
#'   whose position columns are named as a numbered run, as in
#'   [build_network()] (any column beside them other than the `group`
#'   column raises an error of class `nestimate_bad_input`); or a character
#'   matrix with one row per sequence and one column per position, whatever
#'   its column names. A numeric or logical matrix of only 0 and 1, and
#'   one-hot or association networks, raise `nestimate_bad_input`. Every
#'   input except a `netobject_group` needs `group`. For the `print()` and
#'   `plot()` methods: a `net_sequence_comparison`.
#' @param group The group of each sequence: the name of a column of a data
#'   frame or matrix `x` (that column is then not read as a position), or a
#'   vector with one label per sequence. Ignored for a `netobject_group`;
#'   required otherwise. At least two groups are needed. Sequences with a
#'   missing group are left out of the counts with a warning of class
#'   `nestimate_missing_group`. Default `NULL`.
#' @param sub The pattern lengths to count, positive whole numbers.
#'   Lengths longer than the sequences are dropped, and an error is raised
#'   when none remains. Default `3:5`.
#' @param min_freq The smallest count a pattern must reach in every group
#'   to be kept, a positive whole number; when no pattern of any length
#'   reaches it, an error is raised. Default `5`.
#' @param test The inference: `"permutation"` (default), `"chisq"` or
#'   `"none"`; see *Details*.
#' @param iter The number of permutations for `test = "permutation"`, a
#'   positive whole number; ignored otherwise. `iter = 1` gives a warning
#'   of class `nestimate_undefined_effect_size`. Default `1000`.
#' @param adjust The p-value adjustment, any method of [stats::p.adjust()]
#'   (`"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`,
#'   `"fdr"` or `"none"`). Default `"fdr"` (Benjamini-Hochberg).
#'
#' @return A `net_sequence_comparison` object. `print()` shows the numbers
#'   of patterns and groups, the settings and the first 10 rows of the
#'   pattern table, and returns its input invisibly; `summary()` returns the
#'   pattern table, a data frame with one row per pattern; `plot()` draws
#'   the residuals of the top patterns (see *Plot styles*). The object
#'   holds:
#' \describe{
#'   \item{patterns}{The pattern table, one row per k-gram kept. Always:
#'     `pattern`, `length`, and per group `freq_<group>` (the count),
#'     `prop_<group>` (the count as a share of all k-grams of that length
#'     in the group, rare ones included) and `resid_<group>` (the
#'     standardised residual). With `test = "permutation"` also
#'     `effect_size` and `p_value`; with `test = "chisq"`, `statistic` and
#'     `p_value`. The p-values are adjusted. Rows are ordered by ascending
#'     `p_value` (ties by the larger absolute effect size, or the larger
#'     statistic) when a test
#'     was run, and by descending largest absolute residual otherwise.}
#'   \item{groups}{The group names, sorted.}
#'   \item{n_patterns}{The number of rows of the pattern table: the
#'     patterns that reach `min_freq` in every group.}
#'   \item{params}{The settings used: `sub` (after dropping lengths longer
#'     than the sequences), `min_freq`, `test`, `iter` and `adjust`.}
#' }
#'
#' @references
#' Haberman, S. J. (1973). The analysis of residuals in cross-classified
#' tables. \emph{Biometrics}, 29(1), 205--220. (standardised residuals)
#'
#' Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate.
#' \emph{Journal of the Royal Statistical Society B}, 57(1), 289--300.
#' (the default \code{adjust = "fdr"})
#'
#' @examples
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 60, TRUE),
#'   V2 = sample(LETTERS[1:4], 60, TRUE),
#'   V3 = sample(LETTERS[1:4], 60, TRUE),
#'   V4 = sample(LETTERS[1:4], 60, TRUE)
#' )
#' grp <- rep(c("X", "Y"), 30)
#' net <- build_network(seqs, method = "relative")
#' res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
#'
#' @importFrom data.table data.table rbindlist setDT setorderv
#' @importFrom stats p.adjust sd chisq.test
#' @export
sequence_compare <- function(x, group = NULL, sub = 3:5,
                             min_freq = 5L,
                             test = c("permutation", "chisq", "none"),
                             iter = 1000L, adjust = "fdr") {

  test <- match.arg(test)
  min_freq <- as.integer(min_freq)
  iter <- as.integer(iter)
  sub <- as.integer(sub)
  .assert(
    "`sub` must be one or more pattern lengths of 1 or more" =
      length(sub) >= 1L && all(sub >= 1L),
    "`min_freq` must be a whole number of 1 or more" = all(min_freq >= 1L),
    "`iter` must be a whole number of 1 or more" = all(iter >= 1L)
  )
  if (identical(test, "permutation") && isTRUE(all(iter < 2L))) {
    warning(warningCondition(
      "`iter = 1` gives no effect size (NA): it needs at least two permutations.",
      class = "nestimate_undefined_effect_size", call = NULL))
  }
  adjust <- match.arg(adjust, stats::p.adjust.methods)

  # ---- Extract sequences and group labels ----
  parsed <- .cs_parse_input(x, group)
  seq_mat <- parsed$seq_mat
  group_vec <- parsed$group
  n_no_group <- sum(is.na(group_vec))
  if (n_no_group > 0L) {
    warning(warningCondition(
      sprintf(paste0("%d of %d sequences have a missing group and are left ",
                     "out of the counts."), n_no_group, length(group_vec)),
      class = "nestimate_missing_group", call = NULL))
  }
  groups <- sort(unique(group_vec))
  n_groups <- length(groups)
  .assert("Need at least 2 groups" = n_groups >= 2L)

  n <- nrow(seq_mat)
  p <- ncol(seq_mat)
  sub <- sub[sub <= p]
  .assert("No valid pattern lengths (sub > sequence length)" = length(sub) >= 1L)

  pat_list <- .cs_extract_patterns(seq_mat, sub)

  group_idx <- lapply(groups, function(g) which(group_vec == g))
  names(group_idx) <- groups

  row_list <- lapply(seq_along(sub), function(si) {
    k <- sub[si]
    pat_mat <- pat_list[[si]]
    unique_pats <- sort(unique(pat_mat[nzchar(pat_mat)]))
    if (!length(unique_pats)) return(NULL)

    counts_all <- vapply(groups, function(g) {
      vals <- pat_mat[, group_idx[[g]], drop = FALSE]
      idx <- match(vals, unique_pats, nomatch = 0L)
      tabulate(idx, nbins = length(unique_pats))
    }, integer(length(unique_pats)))
    if (!is.matrix(counts_all)) counts_all <- matrix(counts_all, ncol = n_groups)
    rownames(counts_all) <- unique_pats
    colnames(counts_all) <- groups

    total_by_group <- colSums(counts_all)
    N_grand <- sum(total_by_group)
    if (N_grand == 0) return(NULL)

    # min_freq filter: every group must meet threshold
    keep <- apply(counts_all, 1, min) >= min_freq
    if (!any(keep)) return(NULL)
    counts <- counts_all[keep, , drop = FALSE]
    kept_pats <- rownames(counts)
    n_pat <- nrow(counts)

    props <- sweep(counts, 2, total_by_group, `/`)
    props[!is.finite(props)] <- NA_real_

    # Standardized residuals: (o - e) / sqrt(e * (1 - r/N) * (1 - c/N))
    # 2xG table is (pattern, total - pattern). Group totals come from the
    # unfiltered count matrix so "rest" includes rare patterns.
    resids <- matrix(0, nrow = n_pat, ncol = n_groups,
                     dimnames = list(kept_pats, groups))
    col_props <- total_by_group / N_grand
    row_sums <- rowSums(counts)
    row_props <- row_sums / N_grand
    expected <- outer(row_sums, col_props)
    denom <- sqrt(expected * (1 - row_props) *
                    rep(1 - col_props, each = n_pat))
    resids[] <- (counts - expected) / denom
    resids[!is.finite(resids)] <- 0

    freq_df  <- as.data.frame(counts, stringsAsFactors = FALSE)
    names(freq_df)  <- paste0("freq_",  groups)
    prop_df  <- as.data.frame(props,  stringsAsFactors = FALSE)
    names(prop_df)  <- paste0("prop_",  groups)
    resid_df <- as.data.frame(resids, stringsAsFactors = FALSE)
    names(resid_df) <- paste0("resid_", groups)

    rows <- data.table(pattern = kept_pats, length = k)
    rows <- cbind(rows, freq_df, prop_df, resid_df)

    # ---- Inference ----
    if (test == "chisq") {
      chi <- vapply(seq_len(n_pat), function(i) {
        row1 <- counts[i, ]
        row2 <- total_by_group - row1
        tab <- rbind(row1, row2)
        ch <- suppressWarnings(stats::chisq.test(tab))
        c(unname(ch$statistic), ch$p.value)
      }, numeric(2L))
      rows$statistic <- chi[1L, ]
      rows$p_value <- chi[2L, ]
    } else if (test == "permutation") {
      stat_obs <- .cs_pattern_statistic(counts)
      stat_perm <- matrix(0, nrow = n_pat, ncol = iter)
      stat_ge <- integer(n_pat)

      # A loop on purpose: one permutation at a time keeps a single
      # patterns-by-groups table in memory, not iter of them.
      for (it in seq_len(iter)) {
        perm_group <- group_vec[sample(n)]
        perm_idx <- lapply(groups, function(g) which(perm_group == g))
        perm_counts_all <- vapply(seq_along(groups), function(gi) {
          vals <- pat_mat[, perm_idx[[gi]], drop = FALSE]
          idx <- match(vals, unique_pats, nomatch = 0L)
          tabulate(idx, nbins = length(unique_pats))
        }, integer(length(unique_pats)))
        if (!is.matrix(perm_counts_all)) {
          perm_counts_all <- matrix(perm_counts_all, ncol = n_groups)
        }
        rownames(perm_counts_all) <- unique_pats
        perm_counts <- perm_counts_all[kept_pats, , drop = FALSE]
        stat_perm[, it] <- .cs_pattern_statistic(perm_counts)
        stat_ge <- stat_ge + (stat_perm[, it] >= stat_obs)
      }

      perm_mean <- rowMeans(stat_perm)
      perm_sd <- apply(stat_perm, 1, sd)
      perm_sd[perm_sd == 0] <- 1
      rows$effect_size <- (stat_obs - perm_mean) / perm_sd
      rows$p_value <- (stat_ge + 1) / (iter + 1)
    }
    # test == "none" adds nothing

    rows
  })
  row_list <- Filter(Negate(is.null), row_list)

  if (!length(row_list)) {
    stop(errorCondition(
      sprintf("No patterns with min_freq >= %d in all groups.", min_freq),
      class = "nestimate_bad_input", call = NULL))
  }

  results <- rbindlist(row_list)

  # Adjust p-values ONCE across the full pattern set
  if ("p_value" %in% names(results) && adjust != "none") {
    results$p_value <- stats::p.adjust(results$p_value, method = adjust)
  }

  # Sort: by p_value if tested, else by max |residual|
  if ("p_value" %in% names(results)) {
    if ("effect_size" %in% names(results)) {
      results <- results[order(results$p_value, -abs(results$effect_size)), ]
    } else {
      results <- results[order(results$p_value, -results$statistic), ]
    }
  } else {
    resid_cols <- paste0("resid_", groups)
    max_r <- do.call(pmax, lapply(resid_cols, function(cn) abs(results[[cn]])))
    results <- results[order(-max_r), ]
  }
  rownames(results) <- NULL

  structure(list(
    patterns   = as.data.frame(results, stringsAsFactors = FALSE),
    groups     = groups,
    n_patterns = nrow(results),
    params     = list(sub = sub, min_freq = min_freq,
                      test = test, iter = iter, adjust = adjust)
  ), class = "net_sequence_comparison")
}


# ---- Input parsing ----

#' @noRd
.cs_parse_input <- function(x, group) {
  # netobject_group: each element has $data
  if (inherits(x, "netobject_group")) {
    .refuse_onehot_network(x, "sequence_compare()")
    group_names <- names(x)
    mats <- lapply(seq_along(x), function(i) {
      d <- .coerce_sequence_input(x[[i]])
      if (is.null(d)) stop(errorCondition(
        "netobject_group element has no $data.",
        class = "nestimate_bad_input", call = NULL))
      as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
    })
    # Align column count
    max_cols <- max(vapply(mats, ncol, integer(1)))
    mats <- lapply(mats, function(m) {
      if (ncol(m) < max_cols) {
        pad <- matrix(NA_character_, nrow(m), max_cols - ncol(m))
        cbind(m, pad)
      } else m
    })
    seq_mat <- do.call(rbind, mats)
    group_vec <- rep(group_names, vapply(mats, nrow, integer(1)))
    return(list(seq_mat = seq_mat, group = group_vec))
  }

  # netobject or tna: extract $data
  if (inherits(x, "tna") || inherits(x, "netobject") ||
      inherits(x, "cograph_network")) {
    .refuse_onehot_network(x, "sequence_compare()")
    d <- .coerce_sequence_input(x)
    if (is.null(d)) {
      stop(errorCondition("Object has no $data.",
                          class = "nestimate_bad_input", call = NULL))
    }
    seq_mat <- as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
  } else if (is.data.frame(x) || is.matrix(x)) {
    if (is.data.frame(x)) {
      taken <- if (is.character(group) && length(group) == 1L) group
      .refuse_undeclared_roles(x, "sequence_compare()", taken = taken,
                               via_network = TRUE, whole_frame = TRUE)
    }
    .refuse_binary_matrix(x, "sequence_compare()")
    seq_mat <- as.matrix(as.data.frame(x, stringsAsFactors = FALSE))
  } else {
    stop(errorCondition(
      "x must be a netobject_group, netobject, tna, data.frame, or matrix.",
      class = "nestimate_bad_input", call = NULL))
  }

  # Resolve group
  if (is.null(group)) {
    stop(errorCondition("'group' is required when x is not a netobject_group.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (length(group) == 1L && is.character(group) && group %in% colnames(seq_mat)) {
    group_vec <- seq_mat[, group]
    seq_mat <- seq_mat[, colnames(seq_mat) != group, drop = FALSE]
  } else {
    .assert("'group' length must match nrow(x)" = length(group) == nrow(seq_mat))
    group_vec <- as.character(group)
  }

  list(seq_mat = seq_mat, group = group_vec)
}


# ---- Pattern extraction (vectorized) ----

#' @noRd
.cs_extract_patterns <- function(m, sub) {
  n <- nrow(m)
  p <- ncol(m)
  mis <- is.na(m)

  lapply(sub, function(k) {
    n_pos <- p - k + 1L
    if (n_pos < 1L) {
      return(matrix("", nrow = 0, ncol = n))
    }
    # Build pattern strings for each position x sequence
    pats <- lapply(seq_len(n_pos), function(i) {
      idx <- i:(i + k - 1L)
      has_na <- rowSums(mis[, idx, drop = FALSE]) > 0L
      cols <- m[, idx, drop = FALSE]
      pat <- do.call(paste, c(as.data.frame(cols, stringsAsFactors = FALSE),
                              list(sep = "->")))
      pat[has_na] <- ""
      pat
    })
    matrix(unlist(pats, use.names = FALSE), nrow = n_pos, ncol = n,
           byrow = TRUE)
  })
}


# ---- Test statistic ----

#' @noRd
.cs_pattern_statistic <- function(x) {
  n <- sum(x)
  nr <- nrow(x)
  nc <- ncol(x)
  rs <- rowSums(x)
  cs <- colSums(x)
  expected <- outer(rs, cs) / n
  sqrt(rowSums((x - expected)^2))
}


# ---- S3 Methods ----

#' @rdname sequence_compare
#' @param ... Not used.
#' @param object For the `summary()` method: a `net_sequence_comparison`.
#' @param top_n For `plot()`: the number of patterns drawn, the first ones
#'   after ranking by `sort`. Default `10`.
#' @param style For `plot()`: `"auto"` (default) draws the back-to-back
#'   pyramid for exactly 2 groups and the heatmap for any other number;
#'   `"pyramid"` and `"heatmap"` force one. `"pyramid"` with other than 2
#'   groups raises an error.
#' @param sort For `plot()`: how to pick and order the patterns.
#'   `"statistic"` (default) ranks them by the absolute effect size after a
#'   permutation test, by the chi-square statistic after a chi-square test,
#'   and by the largest absolute residual without a test. `"frequency"`
#'   ranks them by their total count across the groups.
#' @param alpha For the pyramid: the significance level for the p-value
#'   labels. Patterns with `p_value < alpha` are starred and drawn in bold
#'   dark text, the rest in plain grey. Default `0.05`. The heatmap and a
#'   comparison without a test ignore it.
#' @param show_residuals For the pyramid: whether to print each bar's
#'   standardised residual inside it. Default `FALSE`. The heatmap always
#'   prints the residuals.
#' @export
print.net_sequence_comparison <- function(x, ...) {
  cat(sprintf("Sequence Comparison  [%d patterns | %d groups: %s]\n",
              x$n_patterns, length(x$groups), paste(x$groups, collapse = ", ")))
  cat(sprintf("  Lengths: %s  |  min_freq: %d",
              paste(x$params$sub, collapse = ", "), x$params$min_freq))
  test <- x$params$test
  if (identical(test, "permutation")) {
    cat(sprintf("  |  permutation: %d iter (%s)",
                x$params$iter, x$params$adjust))
  } else if (identical(test, "chisq")) {
    cat(sprintf("  |  chi-square  (%s)", x$params$adjust))
  }
  cat("\n")

  if (x$n_patterns > 0) {
    top <- utils::head(x$patterns, 10)
    cat("\n")
    print(top, row.names = FALSE)
    if (x$n_patterns > 10) {
      cat(sprintf("  ... and %d more patterns\n", x$n_patterns - 10))
    }
  }
  invisible(x)
}


#' @rdname sequence_compare
#' @export
summary.net_sequence_comparison <- function(object, ...) {
  object$patterns
}


#' @rdname sequence_compare
#' @section Plot styles:
#' `plot()` draws the standardised residuals of the top `top_n` patterns,
#' prints the plot and returns the `ggplot` invisibly (`NULL`, with a
#' message, when there are no patterns). `style = "auto"` (the default)
#' picks the style from the number of groups:
#' \describe{
#'   \item{\code{"pyramid"}}{Back-to-back bars of each pattern's proportion
#'     in the two groups, each bar shaded by that group's standardised
#'     residual (the colour saturates at plus or minus 3), with the
#'     (adjusted) p-value in the gap between the bars when a test was run.
#'     Needs exactly 2 groups.}
#'   \item{\code{"heatmap"}}{One tile per pattern and group, coloured by
#'     the standardised residual (the colour saturates at plus or minus 3)
#'     and labelled with its value. Works for any number of groups.}
#' }
#' The residuals are filled whatever `test` was chosen, so both styles work
#' after `test = "none"`.
#' @import ggplot2
#' @export
plot.net_sequence_comparison <- function(x, top_n = 10L,
                                          style = c("auto", "pyramid", "heatmap"),
                                          sort = c("statistic", "frequency"),
                                          alpha = 0.05,
                                          show_residuals = FALSE, ...) {
  style <- match.arg(style)
  sort  <- match.arg(sort)
  pat <- x$patterns
  groups <- x$groups
  if (nrow(pat) == 0) {
    message("No patterns to plot.")
    return(invisible(NULL))
  }

  # Smart default: the back-to-back pyramid only makes sense for 2 groups,
  # so anything else falls back to the heatmap automatically. An explicit
  # style = "pyramid" on != 2 groups still errors in .sc_plot_pyramid().
  if (style == "auto") style <- if (length(groups) == 2L) "pyramid" else "heatmap"

  resid_cols <- paste0("resid_", groups)
  freq_cols  <- paste0("freq_",  groups)

  if (sort == "frequency") {
    total_freq <- Reduce("+", lapply(freq_cols, function(cn) pat[[cn]]))
    pat <- pat[order(-total_freq), ]
  } else {
    if ("effect_size" %in% names(pat)) {
      pat <- pat[order(-abs(pat$effect_size)), ]
    } else if ("statistic" %in% names(pat)) {
      pat <- pat[order(-pat$statistic), ]
    } else {
      max_r <- do.call(pmax, lapply(resid_cols, function(cn) abs(pat[[cn]])))
      pat <- pat[order(-max_r), ]
    }
  }
  pat <- utils::head(pat, top_n)

  if (style == "heatmap") {
    return(.sc_plot_heatmap(pat, groups, resid_cols))
  }
  .sc_plot_pyramid(pat, groups, alpha, show_residuals)
}

#' @noRd
.sc_plot_pyramid <- function(pat, groups, alpha, show_residuals = FALSE) {
  if (length(groups) != 2L) {
    stop(errorCondition(
      paste0(
        "style = 'pyramid' requires exactly 2 groups; got ",
        length(groups),
        ". Use style = 'heatmap' instead."),
      class = "nestimate_bad_input", call = NULL))
  }
  g1 <- groups[1]; g2 <- groups[2]
  prop1 <- paste0("prop_", g1); prop2 <- paste0("prop_", g2)
  freq1 <- paste0("freq_", g1); freq2 <- paste0("freq_", g2)
  rcol1 <- paste0("resid_", g1); rcol2 <- paste0("resid_", g2)

  # Preserve the rank order passed in (top-ranked first) so the pyramid
  # agrees with the heatmap. Reverse factor levels so row 1 appears at
  # the top of the y-axis rather than the bottom.
  resid1 <- pat[[rcol1]]
  resid2 <- pat[[rcol2]]

  pyramid <- data.frame(
    pattern = factor(pat$pattern, levels = rev(pat$pattern)),
    left    = -pat[[prop1]],
    right   =  pat[[prop2]],
    freq_left   = pat[[freq1]],
    freq_right  = pat[[freq2]],
    resid_left  = resid1,
    resid_right = resid2,
    stringsAsFactors = FALSE
  )

  r_max <- max(abs(c(resid1, resid2)), na.rm = TRUE)
  x_max <- max(abs(c(pyramid$left, pyramid$right)), na.rm = TRUE)
  gutter <- x_max * 0.06
  has_pval <- "p_value" %in% names(pat)
  bar_hw <- 0.42
  y_num <- as.numeric(pyramid$pattern)

  rect_df <- rbind(
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = pyramid$left - gutter, xmax = -gutter,
               freq = pyramid$freq_left, resid = pyramid$resid_left,
               side = "left", stringsAsFactors = FALSE),
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = gutter, xmax = pyramid$right + gutter,
               freq = pyramid$freq_right, resid = pyramid$resid_right,
               side = "right", stringsAsFactors = FALSE)
  )

  if (has_pval) {
    p_num <- ifelse(pat$p_value < 0.001, "<.001",
              ifelse(pat$p_value < 0.01,
                     sprintf("%.3f", pat$p_value),
                     sprintf("%.2f", pat$p_value)))
    # `alpha` is the significance threshold for the pyramid's p-value
    # display: patterns with p < alpha are starred and drawn bold/dark,
    # the rest stay plain grey so the reader sees which patterns clear
    # the chosen level.
    p_sig <- !is.na(pat$p_value) & pat$p_value < alpha
    p_df <- data.frame(
      yidx       = y_num,
      p_label    = ifelse(p_sig, paste0(p_num, "*"), p_num),
      sig        = p_sig,
      p_face     = ifelse(p_sig, "bold.italic", "italic"),
      p_color    = ifelse(p_sig, "grey15", "grey55"),
      stringsAsFactors = FALSE
    )
  }

  # Residual labels inside each bar (optional)
  if (show_residuals) {
    rect_df$x_center <- (rect_df$xmin + rect_df$xmax) / 2
    rect_df$resid_label <- sprintf("%.2f", rect_df$resid)
    rect_df$text_color <- ifelse(abs(rect_df$resid) > r_max * 0.35,
                                  "white", "grey30")
  }

  p <- ggplot() +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = resid)) +
    geom_text(data = rect_df,
              aes(x = ifelse(side == "left", xmin, xmax),
                  y = yidx, label = freq),
              hjust = ifelse(rect_df$side == "left", 1.15, -0.15),
              size = 3, fontface = "bold",
              color = ifelse(abs(rect_df$resid) > r_max * 0.5, "grey25", "grey45")) +
    {if (show_residuals)
      geom_text(data = rect_df,
                aes(x = x_center, y = yidx, label = resid_label),
                size = 2.8, fontface = "bold",
                color = rect_df$text_color)
    } +
    {if (has_pval)
      geom_text(data = p_df,
                aes(x = 0, y = yidx, label = p_label),
                size = 2.5, color = p_df$p_color, fontface = p_df$p_face)
    } +
    scale_fill_gradient2(low = "#CB181D", mid = "white", high = "#2171B5",
                         midpoint = 0, limits = c(-3, 3),
                         oob = scales::squish,
                         name = "Standardized\nresidual") +
    scale_y_continuous(breaks = y_num,
                       labels = as.character(pat$pattern)) +
    scale_x_continuous(
      labels = function(x) format(abs(x), scientific = FALSE)
    ) +
    labs(x = "Proportion", y = NULL,
         title = "Sequence Pattern Comparison",
         subtitle = paste0(g1, "  <<--    -->>  ", g2)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y          = element_text(family = "mono", size = 8),
      axis.ticks.y         = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle        = element_text(hjust = 0.5, color = "grey40",
                                           size = 11),
      legend.position      = "right",
      plot.margin          = margin(5, 10, 5, 5)
    )

  print(p)
  invisible(p)
}

#' @noRd
.sc_plot_heatmap <- function(pat, groups, resid_cols) {
  long <- as.data.frame(data.table::rbindlist(lapply(seq_along(groups), function(gi) {
    data.frame(
      pattern = pat$pattern,
      group   = groups[gi],
      resid   = pat[[resid_cols[gi]]],
      stringsAsFactors = FALSE
    )
  })))
  long$pattern <- factor(long$pattern, levels = rev(pat$pattern))
  long$group   <- factor(long$group, levels = groups)

  r_max <- max(abs(long$resid), na.rm = TRUE)
  # Saturate at +/-3 (|z|>3 is the conventional "very strong" threshold).
  # Anything beyond clips to full color via scales::squish so mid-range
  # residuals (+/-1 to +/-2) pick up visible color instead of fading to white.
  r_lim <- 3

  p <- ggplot(long, aes(x = group, y = pattern, fill = resid)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", resid)),
              color = "white", size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "#CB181D", mid = "white", high = "#2171B5",
                         midpoint = 0, limits = c(-r_lim, r_lim),
                         oob = scales::squish,
                         name = "Standardized\nresidual") +
    labs(x = "Groups", y = NULL, title = "Sequence Pattern Comparison") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y        = element_text(family = "mono", size = 10),
      axis.text.x        = element_text(angle = 0, size = 11,
                                         face = "bold"),
      axis.title.x       = element_text(size = 14),
      panel.grid         = element_blank(),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position    = "right",
      legend.title       = element_text(size = 12),
      legend.text        = element_text(size = 11),
      legend.key.height  = unit(1.2, "cm"),
      plot.margin        = margin(10, 15, 10, 10)
    )

  print(p)
  invisible(p)
}
