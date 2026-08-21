# Observed-only grouping ------------------------------------------------------

# `interaction()` codes a combination as `code(last) + nlevels(last) * code(...)`
# over the *marginal* level space, before unused combinations are dropped. Two
# columns of 46,341 distinct values overflow 32-bit integers, and the resulting
# NAs are pasted into the literal string "NA" downstream, silently merging
# unrelated rows. This assigns an integer per *observed* combination instead, so
# cost tracks the data rather than the Cartesian product, and no separator can
# make two distinct combinations collide.
#
# Group numbering reproduces `interaction()`'s ordering exactly: the last column
# is the primary sort key and the first varies fastest, so prepared row order --
# and therefore same-seed bootstrap output -- is unchanged.
.observed_group_id <- function(cols, context = "grouping") {
  stopifnot(is.list(cols), length(cols) > 0L)

  codes <- lapply(cols, function(x) as.integer(factor(x)))

  incomplete <- names(cols)[vapply(codes, anyNA, logical(1))]
  if (length(incomplete)) {
    stop(
      "Missing values in ", context, " column(s): ",
      paste(incomplete, collapse = ", "),
      ". Every event needs an identifier; drop or relabel these rows first.",
      call. = FALSE
    )
  }

  n <- length(codes[[1L]])
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)

  # rev(): interaction() treats the last column as the primary sort key.
  keys <- rev(codes)
  ord <- do.call(order, keys)
  m <- do.call(cbind, keys)[ord, , drop = FALSE]

  changed <- c(
    TRUE,
    rowSums(m[-1L, , drop = FALSE] != m[-n, , drop = FALSE]) > 0L
  )

  ids <- integer(n)
  ids[ord] <- cumsum(changed)
  ids
}

# Row indices per group, groups ordered by first appearance.
#
# Replaces the `lapply(unique(key), function(k) x[key == k, ])` idiom, which
# rescans the whole key vector once per group and is therefore quadratic in the
# number of groups. `split()` does it in one pass; reordering by each group's
# smallest index restores the first-appearance order that `unique()` gave, so
# downstream row order -- and any same-seed resampling built on it -- is
# unchanged.
.split_first_seen <- function(ids) {
  idx <- split(seq_along(ids), ids)
  if (length(idx) <= 1L) return(idx)
  idx[order(vapply(idx, function(i) i[1L], integer(1)))]
}

# Human-readable label for a set of grouping columns. Display only -- never a
# key, because any separator can occur inside the values themselves.
.group_label <- function(cols, sep = " | ") {
  do.call(paste, c(unname(cols), list(sep = sep)))
}
