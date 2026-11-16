# ---- mosaic_analysis(): two-variable chi-square + flat mosaic ----
#
# A self-contained categorical association analysis for two columns of a
# data.frame: contingency table -> min-count filtering -> chi-square (or
# Fisher) test -> Cramer's V effect size -> standardized residuals -> a flat
# ggplot2 mosaic (shared .mosaic_flat_draw renderer, ColorBrewer RdBu
# ramp). Ported from Saqrmisc::mosaic_analysis(), rewritten in base R
# (no dplyr / tibble / vcd) with tidy one-row-per-cell output.


# Cramer's V effect-size interpretation, using Cohen's df-adjusted thresholds.
#' @noRd
.interpret_cramers_v <- function(v, df) {
  if (is.na(v)) return(NA_character_)
  thr <- if (df <= 1L) {
    c(small = 0.10, medium = 0.30, large = 0.50)
  } else if (df == 2L) {
    c(small = 0.07, medium = 0.21, large = 0.35)
  } else if (df == 3L) {
    c(small = 0.06, medium = 0.17, large = 0.29)
  } else {
    c(small = 0.05, medium = 0.15, large = 0.25)
  }
  if (v < thr["small"])  return("negligible")
  if (v < thr["medium"]) return("small")
  if (v < thr["large"])  return("medium")
  "large"
}


#' Two-variable mosaic analysis (chi-square test + flat mosaic)
#'
#' @description
#' `mosaic_analysis()` tests the association between two categorical columns
#' of a data frame and draws it as a mosaic. It builds their contingency
#' table, drops sparse categories (fewer than `min_count` rows), runs a
#' Pearson chi-square test (or Fisher's exact test), computes Cramer's V
#' with an effect-size label, and draws a flat ggplot2 mosaic in which tile
#' area is the cell count and the fill is the cell's standardised
#' residual. The per-cell results come as a data frame, one row per cell.
#'
#' The p-value tests the null hypothesis that the two variables are
#' independent; with large samples almost any departure is significant, so
#' read Cramer's V for the size of the association. The label follows
#' Cohen's (1988) thresholds for the table's degrees of freedom,
#' `min(rows, columns) - 1`: 0.10 / 0.30 / 0.50 for small /
#' medium / large at 1, 0.07 / 0.21 / 0.35 at 2, 0.06 / 0.17 / 0.29 at 3,
#' and 0.05 / 0.15 / 0.25 above (Cohen's w divided by the square root of the
#' degrees of freedom); a V below the small threshold is labelled
#' `"negligible"`. The residuals are `chisq.test()`'s standardised
#' residuals, roughly standard normal under independence, so one beyond
#' about +-2 marks a cell that holds more (or fewer) rows than independence
#' predicts; with many cells, some pass that line by chance.
#'
#' @section Nested data:
#' Pearson's test treats every row as an independent observation. When rows
#' are nested in actors (events within students, students within classes),
#' rows of one actor resemble each other and the test is too liberal. Set
#' `nested = TRUE` and name the actor column(s) in `actor` (as in
#' [permutation()]). The test is then `nested_method`, one of five that
#' allow for the nesting; each gives a test and its own residual per cell,
#' which fills the mosaic.
#' \describe{
#'   \item{`"rao_scott"` (default)}{Rao-Scott second-order corrected chi-square,
#'     reported as an F statistic on fractional degrees of freedom (the
#'     default of `survey::svychisq()`, which it reproduces for an
#'     unweighted cluster sample). Works for any cluster sizes and when both
#'     variables vary within clusters. `design_effect` is the mean
#'     generalised design effect; 1 means the clustering costs nothing.
#'     Residual: the standardised residual divided by the square root of
#'     the cell's design effect (cluster-robust over simple-random-sampling
#'     variance of `p_ij - p_i. p_.j`); with one row per cluster it is
#'     `chisq.test()`'s, and it is `NA` where the cell's design effect is
#'     zero. With empty cells the design effects are those of the estimable
#'     contrasts. The test is liberal with few actors when one variable is
#'     constant within actors: in simulations of 3 x 4 tables its type I
#'     error at alpha = 0.05 was 0.22 with 10 actors, 0.10 with 30 and 0.06
#'     with 100.}
#'   \item{`"donner_klar"`}{Pearson's chi-square divided, group by group, by
#'     the inflation factor `1 + (m - 1) * rho` (its weighted form for
#'     unequal cluster sizes), where `rho` is the ANOVA intracluster
#'     correlation of `var2`, pooled within the levels of the other
#'     variable. For designs where one variable is constant within each
#'     cluster, as in a cluster-randomised trial; otherwise an error is
#'     raised. `icc` is `rho`; a negative estimate is reported and used as
#'     0. Residual: the standardised residual divided by the square root of
#'     its group's factor.}
#'   \item{`"cmh"`}{Generalised Cochran-Mantel-Haenszel test with each
#'     cluster as a stratum: is there an association within clusters? Suits
#'     a few clusters with enough rows each; both variables must vary within
#'     clusters, and clusters with fewer than two rows are left out
#'     (`clusters` counts those that inform the test). Residual: summed
#'     observed minus expected over the square root of the summed
#'     hypergeometric variances. This tests conditional association, whereas
#'     Rao-Scott and GEE test marginal association; the null hypotheses can
#'     differ. Stratification does not correct arbitrary serial dependence
#'     among observations within an actor.}
#'   \item{`"gee"`}{Marginal model for `var2` given `var1` with an identity
#'     link and independence working correlation; Wald chi-square with the
#'     cluster-robust (sandwich) covariance for equal distributions of
#'     `var2` across `var1`, on `(rows - 1) (columns - 1)` degrees of
#'     freedom. When the covariance has lower rank (fewer actors than
#'     degrees of freedom, or cells that never occur) the test would depend
#'     on the order of rows and columns, and an error of class
#'     `nestimate_bad_input` is raised. The Wald test is liberal unless the
#'     actors are many relative to the degrees of freedom: in simulations
#'     of 3 x 4 tables its type I error at alpha = 0.05 was 0.28 with 30
#'     actors and 0.11 with 100. Residual: the share of `var2` within the
#'     level of `var1` minus the pooled share, over its cluster-robust
#'     standard error.}
#'   \item{`"glmm"`}{Log-linear Poisson mixed model of the cluster-by-cell
#'     counts with random intercepts for cluster-by-`var1` and
#'     cluster-by-`var2` at every level of `actor`, so several columns
#'     give nested random effects (outermost first); Wald chi-square on the
#'     interaction. Fitted by `lme4::glmer(nAGQ = 0)`, a faster, less exact
#'     approximation than the default Laplace fit; needs lme4, and a
#'     failed fit raises an error of class `nestimate_bad_input`. With an
#'     empty cell the Wald test loses power (its interaction term runs to
#'     minus infinity). Residual: observed minus expected count under
#'     independence, over its cluster-robust standard error (the
#'     GLMM supplies the test, the residual is model-free). With several
#'     actor levels, this standard error clusters at the outermost level;
#'     it is undefined with fewer than two outer actors. With about 10
#'     actors this standard error is small and 7 to 8% of null cells
#'     pass 1.96.}
#' }
#' With several `actor` columns, the four other tests use the outermost
#' (the first), which carries the dependence of every level below it.
#' Cramer's V is the pooled-table effect
#' size whatever the test: the clustering changes the uncertainty, not the
#' size of the association.
#'
#' @param data A data frame holding the two variables. Rows missing either
#'   of them (or an `actor` column) are left out; a `var1` or `var2` that
#'   is not a column of `data`, or no complete rows, raises an error.
#' @param var1 The name of the first variable, drawn as the mosaic's columns
#'   (column width is its marginal share).
#' @param var2 The name of the second variable, stacked within each column.
#' @param min_count The smallest number of rows a category must have to be
#'   kept. Categories of either variable with fewer rows in the full table
#'   are dropped, in one pass, before the test. Default `10`. If fewer than
#'   two categories of either variable remain, an error is raised.
#' @param test The test reported: `"chisq"` (default) for Pearson's
#'   chi-square test, or `"fisher"` for Fisher's exact test: exact for a
#'   2 x 2 table, and otherwise with a p-value simulated from 10,000 tables
#'   (so it changes with the random seed and cannot fall below 1/10,001).
#'   Cramer's V and the residuals come from the chi-square statistic. Used
#'   when `nested = FALSE`; given with `nested = TRUE` it is ignored with a
#'   warning of class `nestimate_ignored_argument`.
#' @param nested Whether the rows are nested in actors and the test should
#'   allow for it. Default `FALSE`: every row is independent, as before.
#'   `TRUE` reports `nested_method` and needs `actor`; see *Nested data*.
#' @param nested_method With `nested = TRUE`, the test: `"rao_scott"`
#'   (default), `"donner_klar"`, `"cmh"`, `"gee"` or `"glmm"`. A test that
#'   cannot run on the design raises an error of class
#'   `nestimate_bad_input`; `"glmm"` without lme4 raises one of class
#'   `nestimate_missing_package`.
#' @param actor With `nested = TRUE` (required then): the name of the
#'   column whose values identify the actors, or several names for nested
#'   levels, outermost first (e.g. `c("Group", "Actor")`). `nested_method`
#'   and `actor` given with `nested = FALSE` are ignored with a warning of
#'   class `nestimate_ignored_argument`.
#' @param report Whether to report the test: `TRUE` draws it as the plot
#'   subtitle and prints it as a one-row data frame with the columns of
#'   [mosaic_test()] but `note`. Default `FALSE`. The test is always in
#'   the `stats` table and in `print()`.
#' @param percentage_base What the `pct` column and the `"percent"` tile
#'   label are a percentage of: `"total"` (default) the whole table, `"row"`
#'   each level of `var1`, `"column"` each level of `var2`.
#' @param tile_label What to print inside each tile: `"count"` (default),
#'   `"percent"`, `"residual"`, `"category"` (the `var2` level) or `"none"`.
#' @param title The plot title. Default `""` (none).
#' @param ... Further styling for the flat mosaic: `col_label_side`,
#'   `row_label_side`, `col_label_angle`, `row_label_angle`, `show_legend`,
#'   `legend_position`, `legend_size`, `legend_title`, `label_size`,
#'   `range`, `palette`, `min_label_h` and `min_label_w` (and `pct_base`,
#'   which overrides `percentage_base` in the plot only; the `pct` column
#'   always follows `percentage_base`). The fill is the
#'   ColorBrewer RdBu ramp unless `palette` is given. Column labels turn
#'   vertical when there are more than 6 columns; `col_label_angle` sets an
#'   angle. Any other name raises an error. For the `plot()` method: the
#'   same styling and `title`, overriding what the fit stored; for the
#'   `print()` and `summary()` methods: ignored.
#'
#' @return A `mosaic_analysis` object. `print()` shows the sample size, the
#'   test, Cramer's V with its label and the dropped categories, and returns
#'   the object invisibly; `summary()` returns the per-cell data frame (the
#'   `counts` table below) with the one-row test summary attached as the
#'   attribute `"stats"`; `plot()` draws the mosaic on the active device,
#'   with new styling if given, and returns the `ggplot` invisibly. The
#'   object is a list with:
#' \describe{
#'   \item{plot}{The flat mosaic, a `ggplot`. Building the object does not
#'     draw it.}
#'   \item{counts}{A data frame with one row per (`var1`, `var2`) cell: the
#'     two variables (named after them), `observed`, `expected`, `residual`
#'     (the standardised residual, or the nested test's residual, `NA`
#'     where it is undefined and drawn at 0) and `pct` (on
#'     `percentage_base`); `expected` and `residual` are rounded to 3
#'     decimals, `pct` to 2.}
#'   \item{stats}{A one-row data frame: `test`, `statistic`, `df`, `df2`,
#'     `p_value`, `cramers_v`, `effect_size`, `n`, `clusters`,
#'     `design_effect` and `icc`, rounded (`statistic` and `cramers_v` to 3
#'     decimals, `p_value` to 4 significant digits). `statistic` and `df`
#'     are `NA` for Fisher's test. `df2` (the denominator degrees of freedom
#'     of the Rao-Scott F), `clusters`, `design_effect` and `icc` are `NA`
#'     where the test has none.}
#'   \item{test}{The `htest` object of the test reported; `NULL` for a
#'     nested test.}
#'   \item{actor}{The `actor` column names, or `NULL`.}
#'   \item{cramers_v, effect_size}{Cramer's V, unrounded, and its label.}
#'   \item{table}{The filtered contingency table.}
#'   \item{removed}{The dropped categories, as a list with `var1` and
#'     `var2`.}
#'   \item{n_original, n_filtered}{The number of complete rows before and
#'     after the sparse categories are dropped.}
#'   \item{vars}{The two variable names, as `c(var1 = , var2 = )`.}
#'   \item{plot_parts, plot_args}{The residual matrix, the table, the test
#'     line drawn as the subtitle with `report = TRUE`, and the styling,
#'     kept so `plot()` can redraw without testing again.}
#' }
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#' sciences} (2nd ed.). Lawrence Erlbaum.
#'
#' Rao, J. N. K., & Scott, A. J. (1984). On chi-squared tests for multiway
#' contingency tables with cell proportions estimated from survey data.
#' \emph{The Annals of Statistics}, 12(1), 46--60.
#' \doi{10.1214/aos/1176346391}
#'
#' Donner, A., & Klar, N. (1994). Methods for comparing event rates in
#' intervention studies when the unit of allocation is a cluster.
#' \emph{American Journal of Epidemiology}, 140(3), 279--289.
#' \doi{10.1093/oxfordjournals.aje.a117247}
#'
#' Landis, J. R., Heyman, E. R., & Koch, G. G. (1978). Average partial
#' association in three-way contingency tables: A review and discussion of
#' alternative tests. \emph{International Statistical Review}, 46(3),
#' 237--254. \doi{10.2307/1402373}
#'
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, 73(1), 13--22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' Agresti, A. (2013). \emph{Categorical data analysis} (3rd ed.), chapter
#' 13. Wiley.
#'
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' res <- mosaic_analysis(group_regulation_long, "Course", "Action",
#'                        min_count = 20)
#' res
#' head(summary(res))
#' \donttest{
#' plot(res, tile_label = "percent")
#'
#' # Actions are nested in students: Rao-Scott test by student
#' nested <- mosaic_analysis(group_regulation_long, "Course", "Action",
#'                           nested = TRUE, actor = "Actor")
#' nested
#' }
#' @seealso \code{\link{mosaic_plot}} for the network/table mosaic (which also
#'   accepts \code{style = "flat"}); \code{\link{mosaic_test}} for every
#'   test, clustered or not, side by side in one table.
#' @export
mosaic_analysis <- function(data, var1, var2, min_count = 10L,
                            test = c("chisq", "fisher"),
                            percentage_base = c("total", "row", "column"),
                            tile_label = c("count", "percent", "residual",
                                           "category", "none"),
                            title = "", ..., nested = FALSE,
                            nested_method = c("rao_scott", "donner_klar",
                                              "cmh", "gee", "glmm"),
                            actor = NULL, report = FALSE) {
  .assert(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' has no rows"          = nrow(data) > 0L,
    "'var1' must be a single column name" =
      is.character(var1) && length(var1) == 1L,
    "'var2' must be a single column name" =
      is.character(var2) && length(var2) == 1L
  )
  method <- .resolve_nesting(nested, match.arg(nested_method),
                             !missing(test), !missing(nested_method), actor,
                             other = "test")
  test            <- if (is.null(method)) match.arg(test) else "chisq"
  percentage_base <- match.arg(percentage_base)
  tile_label      <- match.arg(tile_label)
  min_count       <- as.integer(min_count)
  if (!var1 %in% names(data)) stop(errorCondition(
    paste0("Variable '", var1, "' not found in data."),
    class = "nestimate_bad_input", call = NULL))
  if (!var2 %in% names(data)) stop(errorCondition(
    paste0("Variable '", var2, "' not found in data."),
    class = "nestimate_bad_input", call = NULL))
  nested <- !is.null(method)
  if (nested) {
    .assert("`actor` must name one or more columns of `data` when nested = TRUE" =
              is.character(actor) && length(actor) >= 1L &&
              !anyNA(actor) && !anyDuplicated(actor))
    missing_cols <- setdiff(actor, names(data))
    if (length(missing_cols)) stop(errorCondition(
      sprintf("`actor` column(s) not found in data: %s.",
              paste(missing_cols, collapse = ", ")),
      class = "nestimate_bad_input", call = NULL))
    if (any(c(var1, var2) %in% actor)) stop(errorCondition(
      "`actor` must not include `var1` or `var2`.",
      class = "nestimate_bad_input", call = NULL))
  } else {
    actor <- NULL   # already warned by .resolve_nesting()
  }

  # Complete cases on the two variables (and the actor columns), coerced
  # to factors.
  used <- c(var1, var2, actor)
  keep_rows <- Reduce(`&`, lapply(used, function(v) !is.na(data[[v]])))
  clean <- data[keep_rows, used, drop = FALSE]
  if (nrow(clean) == 0L) {
    stop(errorCondition("No complete cases after removing missing values.",
                        class = "nestimate_bad_input", call = NULL))
  }
  clean[[var1]] <- factor(clean[[var1]])
  clean[[var2]] <- factor(clean[[var2]])

  full_tab <- table(clean[[var1]], clean[[var2]])
  keep1 <- names(which(rowSums(full_tab) >= min_count))
  keep2 <- names(which(colSums(full_tab) >= min_count))
  removed <- list(var1 = setdiff(rownames(full_tab), keep1),
                  var2 = setdiff(colnames(full_tab), keep2))

  filtered <- clean[clean[[var1]] %in% keep1 & clean[[var2]] %in% keep2, ,
                    drop = FALSE]
  filtered[[var1]] <- droplevels(filtered[[var1]])
  filtered[[var2]] <- droplevels(filtered[[var2]])
  tab <- table(filtered[[var1]], filtered[[var2]])
  if (nrow(tab) < 2L || ncol(tab) < 2L) {
    stop(errorCondition(
      "Not enough categories remain after filtering; lower 'min_count'.",
      class = "nestimate_bad_input", call = NULL))
  }

  # Chi-square is always computed (for residuals + Cramer's V); the reported
  # test object is Fisher's when requested.
  chi <- suppressWarnings(stats::chisq.test(tab))
  if (identical(test, "fisher")) {
    htest    <- stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000L)
    p_value  <- htest$p.value
    stat_val <- NA_real_
    df_val   <- NA_real_
  } else {
    htest    <- chi
    p_value  <- chi$p.value
    stat_val <- as.numeric(chi$statistic)
    df_val   <- as.numeric(chi$parameter)
  }

  n  <- sum(tab)
  k  <- min(dim(tab))
  cramers_v   <- as.numeric(sqrt(as.numeric(chi$statistic) / (n * (k - 1L))))
  effect_size <- .interpret_cramers_v(cramers_v, k - 1L)

  # Clustered test: one table per actor, the residuals and the test from
  # the chosen method. Cramer's V stays the pooled-table effect size.
  residual_mat <- chi$stdres
  nested_test <- NULL
  if (nested) {
    lv <- filtered[, actor, drop = FALSE]
    key <- do.call(paste, c(unname(lapply(lv, as.character)), sep = "\r"))
    first <- !duplicated(key)
    arr <- unclass(table(filtered[[var1]], filtered[[var2]],
                         factor(key, levels = unique(key))))
    fit <- .nested_association(arr, method,
                               levels = lv[first, , drop = FALSE],
                               predictor = "rows")
    residual_mat <- fit$residuals
    nested_test <- fit$test
    htest    <- NULL
    p_value  <- fit$test$p_value
    stat_val <- fit$test$statistic
    df_val   <- fit$test$df
  }

  # Tidy one-row-per-cell output: observed / expected / std residual / percent.
  cells <- expand.grid(v1 = rownames(tab), v2 = colnames(tab),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  idx <- cbind(cells$v1, cells$v2)
  observed <- as.numeric(tab[idx])
  base_den <- switch(percentage_base,
    total  = n,
    row    = as.numeric(rowSums(tab)[cells$v1]),
    column = as.numeric(colSums(tab)[cells$v2]))
  counts <- data.frame(
    observed = observed,
    expected = round(as.numeric(chi$expected[idx]), 3),
    residual = round(as.numeric(residual_mat[idx]), 3),
    pct      = round(observed / base_den * 100, 2),
    stringsAsFactors = FALSE
  )
  counts <- cbind(stats::setNames(cells, c(var1, var2)), counts)

  # The default table keeps its 1.0.12 columns; a nested test adds df2,
  # clusters, design_effect and icc.
  stats_df <- data.frame(
    test        = if (nested) nested_test$test else
      if (identical(test, "fisher")) "Fisher's exact" else "Chi-square",
    statistic   = round(stat_val, 3),
    df          = if (nested) round(df_val, 3) else df_val,
    p_value     = signif(p_value, 4),
    cramers_v   = round(cramers_v, 3),
    effect_size = effect_size,
    n           = n,
    stringsAsFactors = FALSE
  )
  if (nested) {
    stats_df <- data.frame(
      stats_df[c("test", "statistic", "df")],
      df2 = round(nested_test$df2, 1),
      stats_df[c("p_value", "cramers_v", "effect_size", "n")],
      clusters = nested_test$clusters,
      design_effect = round(nested_test$design_effect, 3),
      icc = round(nested_test$icc, 3),
      stringsAsFactors = FALSE)
  }

  # Flat mosaic via the shared renderer, fed the residual matrix directly
  # (no second decomposition). An undefined clustered residual (NA) is drawn
  # at 0, as mosaic_plot() draws undefined residuals. `tile_label`/`pct_base`
  # default from the named args but the caller's `...` overrides win, so
  # passing e.g. pct_base = "row" never collides. The plot is stored, not
  # drawn -- construction is silent; call plot() (or print $plot) to draw.
  draw_res <- residual_mat
  draw_res[!is.finite(draw_res)] <- 0
  flat_args <- utils::modifyList(
    c(list(title = title, tile_label = tile_label, pct_base = percentage_base),
      if (nested) list(legend_title = "Cluster-adjusted\nresidual")),
    list(...))
  .mosaic_flat_check_args(names(flat_args))
  # report = TRUE: the test reported, as the subtitle and a printed row.
  .assert("`report` must be TRUE or FALSE" =
            is.logical(report) && length(report) == 1L && !is.na(report))
  subtitle <- NULL
  if (report) {
    reported <- if (nested) nested_test else
      data.frame(test = stats_df$test, statistic = stat_val, df = df_val,
                 df2 = NA_real_, p_value = p_value, clusters = NA_integer_,
                 design_effect = NA_real_, icc = NA_real_,
                 stringsAsFactors = FALSE)
    subtitle <- .nested_test_label(reported)
    .mosaic_report(reported, method %||% test, n)
  }
  p <- do.call(.mosaic_flat_draw,
               c(list(res = draw_res, tab = tab), flat_args))
  if (!is.null(subtitle)) p <- p + ggplot2::labs(subtitle = subtitle)

  # Default objects keep their 1.0.12 fields; `actor` and the stored
  # subtitle are added only when used.
  plot_parts <- list(res = draw_res, tab = tab)
  if (!is.null(subtitle)) plot_parts$subtitle <- subtitle
  out <- list(
    plot        = p,
    counts      = counts,
    stats       = stats_df,
    test        = htest,
    cramers_v   = cramers_v,
    effect_size = effect_size,
    table       = tab,
    removed     = removed,
    n_original  = nrow(clean),
    n_filtered  = nrow(filtered),
    vars        = c(var1 = var1, var2 = var2),
    plot_parts  = plot_parts,
    plot_args   = flat_args
  )
  if (nested) out <- append(out, list(actor = actor), after = 11L)
  structure(out, class = c("mosaic_analysis", "list"))
}


#' @rdname mosaic_analysis
#' @section Methods:
#' * `plot.mosaic_analysis()` redraws the flat mosaic from the stored
#'   contingency table and residuals, so the styling can change without the
#'   test being run again. Any flat-mosaic styling argument (`tile_label`,
#'   `pct_base`, `col_label_side`, `legend_size`, ...) can be overridden
#'   through `...`.
#' @param x For the `plot()` and `print()` methods: a `mosaic_analysis`.
#' @param object For the `summary()` method: a `mosaic_analysis`.
#' @export
plot.mosaic_analysis <- function(x, ...) {
  args <- utils::modifyList(x$plot_args, list(...))
  .mosaic_flat_check_args(names(args))
  p <- do.call(.mosaic_flat_draw,
               c(list(res = x$plot_parts$res, tab = x$plot_parts$tab), args))
  if (!is.null(x$plot_parts$subtitle)) {
    p <- p + ggplot2::labs(subtitle = x$plot_parts$subtitle)
  }
  print(p)
  invisible(p)
}


#' @rdname mosaic_analysis
#' @export
print.mosaic_analysis <- function(x, ...) {
  cat(sprintf("Mosaic analysis: %s x %s\n", x$vars[["var1"]], x$vars[["var2"]]))
  cat(sprintf("  N = %d (filtered from %d); table %d x %d\n",
              x$n_filtered, x$n_original, nrow(x$table), ncol(x$table)))
  s <- x$stats
  # A p-value below machine precision is stored as 0; it is not zero.
  p_text <- if (isTRUE(s$p_value < .Machine$double.eps)) {
    "p < 2.2e-16"
  } else {
    paste("p =", format(s$p_value))
  }
  if ("clusters" %in% names(s)) {        # a nested test
    cat(sprintf("  %s\n", .nested_test_label(x$stats)))
    cat(sprintf("  Actor: %s\n", paste(x$actor, collapse = " / ")))
  } else if (is.na(s$statistic)) {
    cat(sprintf("  %s: %s\n", s$test, p_text))
  } else {
    cat(sprintf("  %s: X2 = %.3f, df = %g, %s\n",
                s$test, s$statistic, s$df, p_text))
  }
  cat(sprintf("  Cramer's V = %.3f (%s)\n", s$cramers_v, s$effect_size))
  rm1 <- x$removed$var1; rm2 <- x$removed$var2
  if (length(rm1) || length(rm2)) {
    cat(sprintf("  Dropped (min_count): %s\n",
                paste(c(rm1, rm2), collapse = ", ")))
  }
  invisible(x)
}


#' @rdname mosaic_analysis
#' @export
summary.mosaic_analysis <- function(object, ...) {
  out <- object$counts
  attr(out, "stats") <- object$stats
  out
}
