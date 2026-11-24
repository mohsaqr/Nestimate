# Results that are data frames (owner rules, 2026-11-24): a result is a data
# frame of class c("net_<x>", "data.frame") whose settings are attributes;
# groups are a `group` column, wide is an option; summary() returns every
# row as a plain data frame, and a subset is a plain data frame, since the
# settings describe the whole result. These helpers are shared by every
# such result so they print and subset alike.

#' A result table: a data frame with a class and its settings as attributes
#'
#' @param table A data frame.
#' @param class The result's class, put before `"data.frame"`.
#' @param ... The settings, stored as attributes.
#' @return `table` with class `c(class, "data.frame")` and the attributes.
#' @noRd
.result_df <- function(table, class, ...) {
  stopifnot(is.data.frame(table))
  rownames(table) <- NULL
  structure(table, class = c(class, "data.frame"), ...)
}

#' A result table as a plain data frame, its settings dropped
#'
#' @param x A data frame, possibly classed and carrying settings.
#' @return A `data.frame` with only its names and row names.
#' @noRd
.plain_df <- function(x) {
  kept <- attributes(x)[c("names", "row.names")]
  attributes(x) <- c(kept, list(class = "data.frame"))
  x
}

#' Print the first rows of a result table, rounded
#'
#' @param top The rows to show, a data frame.
#' @param digits Decimals for double columns other than p-values.
#' @param p_cols Columns printed at three significant digits.
#' @return `top`, invisibly; called for its output.
#' @noRd
.print_result_rows <- function(top, digits,
                               p_cols = c("p_value", "p_adj")) {
  p_cols <- intersect(p_cols, names(top))
  num <- setdiff(names(top)[vapply(top, is.double, logical(1L))], p_cols)
  top[num] <- lapply(top[num], round, digits = digits)
  # each p-value to 3 significant digits on its own, not padded to the
  # smallest one in the column
  top[p_cols] <- lapply(top[p_cols], function(p) {
    vapply(p, function(v) if (is.na(v)) "NA" else format(signif(v, 3L)),
           character(1L))
  })
  print(.plain_df(top), row.names = FALSE, right = FALSE)
  invisible(top)
}

#' Check the n and digits arguments of a print method
#' @noRd
.check_print_args <- function(n, digits) {
  whole <- function(v, low) {
    length(v) == 1L && is.numeric(v) && is.finite(v) && v >= low &&
      v == round(v)
  }
  # n = Inf prints every row
  if (!(whole(n, 1) || (length(n) == 1L && identical(as.numeric(n), Inf)))) {
    stop(errorCondition("`n` must be a single positive whole number, or Inf.",
                        class = "nestimate_bad_input", call = NULL))
  }
  if (!whole(digits, 0)) {
    stop(errorCondition("`digits` must be a single whole number of 0 or more.",
                        class = "nestimate_bad_input", call = NULL))
  }
  invisible(TRUE)
}
