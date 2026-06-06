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
#' Analyses the association between two categorical columns of a data.frame.
#' Builds the contingency table, drops sparse categories below
#' \code{min_count}, runs a Pearson chi-square test (or Fisher's exact test),
#' computes Cramer's V with a df-adjusted effect-size label, and draws a flat
#' ggplot2 mosaic whose tile area encodes counts and whose fill encodes the
#' standardized Pearson residual (Nestimate diverging palette). All tabular
#' output is a tidy one-row-per-cell \code{data.frame}.
#'
#' @param data A data.frame containing the two variables.
#' @param var1 Character. Name of the first variable (mosaic columns).
#' @param var2 Character. Name of the second variable (stacked within columns).
#' @param min_count Integer. Minimum marginal count for a category to be kept.
#'   Categories of either variable below this are dropped before testing.
#'   Default 10.
#' @param test Character. \code{"chisq"} (default) Pearson chi-square, or
#'   \code{"fisher"} Fisher's exact test (simulated p-value). Cramer's V and
#'   the residual fill are always derived from the chi-square statistic.
#' @param percentage_base Character. Base for the \code{"percent"} tile label
#'   and the \code{pct} column: \code{"total"} (default), \code{"row"}
#'   (within \code{var1}), or \code{"column"} (within \code{var2}).
#' @param tile_label Character. What to print inside each tile: \code{"count"}
#'   (default), \code{"percent"}, \code{"residual"}, \code{"category"}
#'   (\code{var2} level), or \code{"none"}.
#' @param title Character. Plot title. Default \code{""}.
#' @param ... Further flat-mosaic styling arguments passed to the renderer
#'   (e.g. \code{col_label_side}, \code{row_label_side}, \code{legend_position},
#'   \code{legend_size}, \code{label_size}, \code{palette}). Tile fill uses the
#'   ColorBrewer RdBu ramp by default (override with \code{palette}). Column
#'   labels auto-rotate to vertical when there are more than 6 columns; pass
#'   \code{col_label_angle} to force an angle.
#'
#' @return An object of class \code{"mosaic_analysis"}: a list with
#' \describe{
#'   \item{plot}{The flat mosaic \code{ggplot} object.}
#'   \item{counts}{Tidy data.frame, one row per (var1, var2) cell, with
#'     \code{observed}, \code{expected}, \code{residual} (standardized), and
#'     \code{pct} (on \code{percentage_base}).}
#'   \item{stats}{One-row data.frame: \code{test}, \code{statistic},
#'     \code{df}, \code{p_value}, \code{cramers_v}, \code{effect_size},
#'     \code{n}.}
#'   \item{test}{The raw \code{htest} object.}
#'   \item{cramers_v, effect_size}{Effect size value and label.}
#'   \item{table}{The filtered contingency \code{table}.}
#'   \item{removed}{List of dropped \code{var1} / \code{var2} categories.}
#'   \item{n_original, n_filtered}{Row counts before/after filtering.}
#' }
#'
#' @examples
#' df <- data.frame(
#'   gender = sample(c("F", "M"), 200, replace = TRUE),
#'   level  = sample(c("Low", "Mid", "High"), 200, replace = TRUE)
#' )
#' res <- mosaic_analysis(df, "gender", "level", min_count = 5)
#' res$stats
#' res$counts
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
  stopifnot(
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
  if (!var1 %in% names(data)) stop("Variable '", var1, "' not found in data.",
                                   call. = FALSE)
  if (!var2 %in% names(data)) stop("Variable '", var2, "' not found in data.",
                                   call. = FALSE)

  # Complete cases on the two variables, coerced to factors.
  keep_rows <- !is.na(data[[var1]]) & !is.na(data[[var2]])
  clean <- data[keep_rows, c(var1, var2), drop = FALSE]
  if (nrow(clean) == 0L) {
    stop("No complete cases after removing missing values.", call. = FALSE)
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
    stop("Not enough categories remain after filtering; lower 'min_count'.",
         call. = FALSE)
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


#' Plot method for mosaic_analysis objects
#'
#' @description
#' Re-renders the flat mosaic from the stored contingency table and residuals,
#' so styling can be changed without re-running the test. Any flat-mosaic
#' styling argument (\code{tile_label}, \code{pct_base}, \code{col_label_side},
#' \code{legend_size}, ...) may be overridden via \code{...}.
#'
#' @param x A \code{mosaic_analysis} object.
#' @param ... Styling overrides forwarded to the flat renderer.
#' @return The \code{ggplot} object, invisibly (drawn as a side effect).
#' @examples
#' df <- data.frame(
#'   a = sample(c("X", "Y", "Z"), 200, replace = TRUE),
#'   b = sample(c("P", "Q"), 200, replace = TRUE)
#' )
#' res <- mosaic_analysis(df, "a", "b", min_count = 5)
#' \donttest{
#' plot(res, tile_label = "percent", legend_position = "bottom")
#' }
#' @export
plot.mosaic_analysis <- function(x, ...) {
  args <- utils::modifyList(x$plot_args, list(...))
  .mosaic_flat_check_args(names(args))
  p <- do.call(.mosaic_flat_draw,
               c(list(res = x$plot_parts$res, tab = x$plot_parts$tab), args))
  print(p)
  invisible(p)
}


#' Print method for mosaic_analysis objects
#'
#' @param x A \code{mosaic_analysis} object.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
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


#' Summary method for mosaic_analysis objects
#'
#' @param object A \code{mosaic_analysis} object.
#' @param ... Ignored.
#' @return The tidy per-cell \code{counts} data.frame, with the \code{stats}
#'   row attached as an attribute.
#' @export
summary.mosaic_analysis <- function(object, ...) {
  out <- object$counts
  attr(out, "stats") <- object$stats
  out
}
