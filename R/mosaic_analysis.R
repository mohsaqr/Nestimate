# ---- mosaic_analysis(): two-variable chi-square + flat mosaic ----
#
# A self-contained categorical association analysis for two columns of a
# data.frame: contingency table -> min-count filtering -> chi-square (or
# Fisher) test -> Cramer's V effect size -> standardized residuals -> a flat
# ggplot2 mosaic (shared .mosaic_flat_draw renderer, Nestimate diverging
# palette). Ported from Saqrmisc::mosaic_analysis(), rewritten in base R
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
#' @param data A data frame holding the two variables. Rows missing either
#'   of them are left out; a `var1` or `var2` that is not a column of `data`,
#'   or no complete rows, raises an error.
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
#'   Cramer's V and the residuals always come
#'   from the chi-square statistic.
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
#'     (the standardised residual) and `pct` (on `percentage_base`);
#'     `expected` and `residual` are rounded to 3 decimals, `pct` to 2.}
#'   \item{stats}{A one-row data frame: `test`, `statistic`, `df`,
#'     `p_value`, `cramers_v`, `effect_size` and `n`, rounded (`statistic`
#'     and `cramers_v` to 3 decimals, `p_value` to 4 significant digits).
#'     `statistic` and `df` are `NA` for Fisher's test.}
#'   \item{test}{The `htest` object of the test reported.}
#'   \item{cramers_v, effect_size}{Cramer's V, unrounded, and its label.}
#'   \item{table}{The filtered contingency table.}
#'   \item{removed}{The dropped categories, as a list with `var1` and
#'     `var2`.}
#'   \item{n_original, n_filtered}{The number of complete rows before and
#'     after the sparse categories are dropped.}
#'   \item{vars}{The two variable names, as `c(var1 = , var2 = )`.}
#'   \item{plot_parts, plot_args}{The residual matrix, the table and the
#'     styling, kept so `plot()` can redraw without testing again.}
#' }
#'
#' @references
#' Cohen, J. (1988). \emph{Statistical power analysis for the behavioral
#' sciences} (2nd ed.). Lawrence Erlbaum.
#'
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' res <- mosaic_analysis(group_regulation_long, "Course", "Action",
#'                        min_count = 20)
#' res
#' head(summary(res))
#' \donttest{
#' plot(res, tile_label = "percent")
#' }
#' @seealso \code{\link{mosaic_plot}} for the network/table mosaic (which also
#'   accepts \code{style = "flat"}).
#' @export
mosaic_analysis <- function(data, var1, var2, min_count = 10L,
                            test = c("chisq", "fisher"),
                            percentage_base = c("total", "row", "column"),
                            tile_label = c("count", "percent", "residual",
                                           "category", "none"),
                            title = "", ...) {
  .assert(
    "'data' must be a data.frame" = is.data.frame(data),
    "'data' has no rows"          = nrow(data) > 0L,
    "'var1' must be a single column name" =
      is.character(var1) && length(var1) == 1L,
    "'var2' must be a single column name" =
      is.character(var2) && length(var2) == 1L
  )
  test            <- match.arg(test)
  percentage_base <- match.arg(percentage_base)
  tile_label      <- match.arg(tile_label)
  min_count       <- as.integer(min_count)
  if (!var1 %in% names(data)) stop(errorCondition(
    paste0("Variable '", var1, "' not found in data."),
    class = "nestimate_bad_input", call = NULL))
  if (!var2 %in% names(data)) stop(errorCondition(
    paste0("Variable '", var2, "' not found in data."),
    class = "nestimate_bad_input", call = NULL))

  # Complete cases on the two variables, coerced to factors.
  keep_rows <- !is.na(data[[var1]]) & !is.na(data[[var2]])
  clean <- data[keep_rows, c(var1, var2), drop = FALSE]
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
    residual = round(as.numeric(chi$stdres[idx]), 3),
    pct      = round(observed / base_den * 100, 2),
    stringsAsFactors = FALSE
  )
  counts <- cbind(stats::setNames(cells, c(var1, var2)), counts)

  stats_df <- data.frame(
    test        = if (identical(test, "fisher")) "Fisher's exact" else "Chi-square",
    statistic   = round(stat_val, 3),
    df          = df_val,
    p_value     = signif(p_value, 4),
    cramers_v   = round(cramers_v, 3),
    effect_size = effect_size,
    n           = n,
    stringsAsFactors = FALSE
  )

  # Flat mosaic via the shared renderer, fed the chi-square residual matrix
  # directly (no second chi-square decomposition). `tile_label`/`pct_base`
  # default from the named args but the caller's `...` overrides win, so
  # passing e.g. pct_base = "row" never collides. The plot is stored, not
  # drawn -- construction is silent; call plot() (or print $plot) to draw.
  flat_args <- utils::modifyList(
    list(title = title, tile_label = tile_label, pct_base = percentage_base),
    list(...))
  .mosaic_flat_check_args(names(flat_args))
  p <- do.call(.mosaic_flat_draw,
               c(list(res = chi$stdres, tab = tab), flat_args))

  structure(
    list(
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
      plot_parts  = list(res = chi$stdres, tab = tab),
      plot_args   = flat_args
    ),
    class = c("mosaic_analysis", "list")
  )
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
  if (is.na(s$statistic)) {
    cat(sprintf("  %s: p = %s\n", s$test, format(s$p_value)))
  } else {
    cat(sprintf("  %s: X2 = %.3f, df = %g, p = %s\n",
                s$test, s$statistic, s$df, format(s$p_value)))
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
