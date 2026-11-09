# ---- Path Dependence: Per-Context Order-k vs Order-1 Diagnostic ----

#' @noRd
.pd_extract_data <- function(x) {
  if (is.matrix(x))     return(as.data.frame(x, stringsAsFactors = FALSE))
  if (is.data.frame(x)) return(x)
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    if (!is.null(x$data) && (is.data.frame(x$data) || is.matrix(x$data)))
      return(as.data.frame(x$data, stringsAsFactors = FALSE))
    stop(errorCondition(
      paste0("netobject does not carry the source sequence data. Pass the wide ",
      "sequence data.frame directly to path_dependence()."),
      class = "nestimate_bad_input", call = NULL)
    )
  }
  stop(errorCondition(
    paste0("'x' must be a wide sequence data.frame, matrix, or netobject built ",
    "with raw data preserved."),
    class = "nestimate_bad_input", call = NULL)
  )
}

#' @noRd
.pd_kgram_counts <- function(data, k) {
  ## Wide-format only. Vectorised k-gram extraction across columns,
  ## results pooled across rows. NAs in any of the k positions zero
  ## out the contribution.
  m <- matrix(vapply(data, as.character, character(nrow(data))),
              nrow = nrow(data))
  ncols <- ncol(m)
  if (ncols < k) {
    stop(errorCondition(
      sprintf("Sequences too short for order %d (need >=%d columns, got %d).",
        k - 1L, k, ncols),
      class = "nestimate_bad_input", call = NULL))
  }
  windows <- vapply(seq_len(ncols - k + 1L), function(t) {
    do.call(paste, c(lapply(seq_len(k), function(j) m[, t + j - 1L]),
                     list(sep = "\x01")))
  }, character(nrow(m)))
  # Drop windows that touch a missing cell -- tested on the cells, not on
  # the pasted text, where a state such as "FINAL" contains "NA".
  touches_na <- vapply(seq_len(ncols - k + 1L), function(t) {
    rowSums(is.na(m[, t:(t + k - 1L), drop = FALSE])) > 0L
  }, logical(nrow(m)))
  pooled <- as.vector(windows)[!as.vector(touches_na)]
  table(pooled)
}

#' @noRd
.pd_split_kgram <- function(kgram_names, k) {
  parts  <- strsplit(kgram_names, "\x01", fixed = TRUE)
  ctx    <- vapply(parts,
                   function(p) paste(p[seq_len(k - 1L)], collapse = " -> "),
                   character(1))
  nxt    <- vapply(parts, function(p) p[k], character(1))
  list(context = ctx, next_state = nxt)
}

#' @noRd
.pd_kl <- function(p, q, base = 2) {
  ## Returns KL(p || q) in `base` units. Defines 0 log 0 := 0; if p[i] > 0
  ## but q[i] == 0, returns Inf.
  msk <- p > 0
  if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
  sum(p[msk] * log(p[msk] / q[msk], base = base))
}

#' @noRd
.pd_entropy <- function(p, base = 2) {
  msk <- p > 0
  -sum(p[msk] * log(p[msk], base = base))
}


# ---- path_dependence() ----

#' Per-Context Path Dependence at Order k
#'
#' @description
#' `path_dependence()` shows, context by context, where a first-order
#' Markov chain misses what longer histories say about the next state. For
#' every context of k = `order` consecutive states
#' \eqn{c = (s_1, \ldots, s_k)}{c = (s_1, ..., s_k)}, it compares the
#' observed distribution of the next state after that context,
#' \eqn{P_k(s_{k+1} \mid c)}{P_k(s_(k+1) | c)}, with the first-order
#' prediction \eqn{P_1(s_{k+1} \mid s_k)}{P_1(s_(k+1) | s_k)} from the last
#' state alone. Its table has one row per context, sorted by the
#' Kullback-Leibler divergence between the two, with both entropies and
#' whether the most likely next state changes. It reads wide sequence data
#' or a network built from sequences.
#'
#' A KL of 0 means the longer history adds nothing in that context; a
#' positive entropy drop means it sharpens the prediction; and a flip means
#' it changes which next state is most likely, not only how confident the
#' prediction is. The values are computed from the observed counts and
#' come without a test, and a context seen a few times gives a noisy KL:
#' read `n` beside it and raise `min_count` when small contexts lead the
#' table. Whether order k is needed for the chain as a whole is tested by
#' [markov_order_test()]; `path_dependence()` breaks that decision down by
#' context.
#'
#' @param x The sequences: a wide data frame or matrix (one row per
#'   sequence, one column per time step), or a `netobject` built from
#'   sequences, whose stored data are read. A data frame needs position
#'   columns named as a numbered run (`T1`, `T2`, ...); any other frame, a
#'   one-hot network and an association network raise an error of class
#'   `nestimate_bad_input`, and a network without its data stops
#'   `path_dependence()`. Paths that touch a missing cell are left out. For
#'   the `print()` methods: a `net_path_dependence` object or its
#'   `summary()`; for the `plot()` method: a `net_path_dependence` object.
#' @param order The length k of the conditioning context, a whole number of
#'   at least 2; a non-integer value is an error, not truncated. Default
#'   `2`, which compares a two-state history with the last state alone; `3`
#'   compares a three-state history, and so on. The data need at least
#'   `order + 1` time steps.
#' @param min_count The number of times a context must be observed to be
#'   reported, a whole number of at least 1. Default `5`. Contexts seen less
#'   often give noisy KL estimates; when no context reaches `min_count`,
#'   `path_dependence()` stops with an error.
#' @param base The logarithm base for the entropies and the KL divergence:
#'   `2` (default) gives bits, `exp(1)` nats and `10` hartleys. A positive
#'   number other than 1.
#'
#' @return A `net_path_dependence` object. `print()` shows the chain-level
#'   summaries and the `top` contexts by KL; `summary()` returns a
#'   `summary.net_path_dependence` holding the full sorted table, which
#'   prints every context and the chain-level summaries; `plot()` returns a
#'   ggplot lollipop chart of each context's KL, highest first, with point
#'   size proportional to the context's count and the contexts whose modal
#'   next state flips coloured and labelled with the change. The `print()`
#'   methods return their input invisibly. The object's fields:
#' \describe{
#'   \item{contexts}{Data frame with one row per context, sorted by KL
#'     (descending), with columns `context` (such as `"A -> B"`), `n` (its
#'     count), `H_order1` (entropy of \eqn{P_1(\cdot \mid s_k)}{P_1(. | s_k)}),
#'     `H_orderk` (entropy of \eqn{P_k(\cdot \mid c)}{P_k(. | c)}), `H_drop`
#'     (`H_order1 - H_orderk`), `KL`
#'     (\eqn{D_{KL}(P_k \| P_1)}{KL(P_k || P_1)}), `top_o1` (the most
#'     likely next state under order 1), `top_ok` (the most likely next
#'     state under order k) and `flips` (whether the two differ).}
#'   \item{chain}{List of chain-level summaries: `KL_weighted` (the
#'     count-weighted mean KL over contexts), `H_drop_weighted` (the
#'     count-weighted mean entropy drop), `n_contexts` and `n_flips` (the
#'     number of contexts whose most likely next state changes).}
#'   \item{order}{The context length used.}
#'   \item{base}{The logarithm base used.}
#'   \item{min_count}{The minimum context count used.}
#'   \item{states}{The states observed.}
#' }
#'
#' @details
#' For each context \eqn{c} observed at least `min_count` times, the
#' function computes:
#' \itemize{
#'   \item the empirical conditional \eqn{P_k(\cdot \mid c)}{P_k(. | c)}
#'     from the counts of the paths of k + 1 states;
#'   \item the first-order prediction \eqn{P_1(\cdot \mid s_k)}{P_1(. | s_k)}
#'     from the most recent state alone, estimated from all pairs of
#'     consecutive states;
#'   \item the entropy drop \eqn{H(P_1) - H(P_k)}{H(P_1) - H(P_k)}, the
#'     uncertainty removed by the longer memory in this context;
#'   \item the Kullback-Leibler divergence
#'     \eqn{D_{KL}(P_k \,\|\, P_1)}{KL(P_k || P_1)}, the expected extra
#'     surprise of predicting with the first-order model when the order-k
#'     distribution holds.
#' }
#'
#' A negative `H_drop` means the order-k context spreads its probability
#' over more outcomes than the last state alone does, which can be
#' small-sample noise or uncertainty the context adds; inspect `n`.
#'
#' @examples
#' data(trajectories, package = "Nestimate")
#' pd <- path_dependence(as.data.frame(trajectories), order = 2)
#' print(pd)
#' summary(pd)
#' plot(pd)
#'
#' @seealso \code{\link{markov_order_test}}, \code{\link{transition_entropy}},
#'   \code{\link{build_mogen}}
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information
#' Theory}, 2nd ed., chapters 2 and 4. Wiley. (KL divergence and
#' conditional entropy.)
#'
#' @export
path_dependence <- function(x, order = 2L, min_count = 5L, base = 2) {
  .assert(
    "'order' must be a whole number" =
      is.numeric(order) && length(order) == 1L &&
      order == round(order)
  )
  order     <- as.integer(order)
  min_count <- as.integer(min_count)
  .assert(
    "'order' must be >= 2" = order >= 2L,
    "'min_count' must be >= 1" = min_count >= 1L,
    is.numeric(base), length(base) == 1L, base > 0, base != 1
  )

  .refuse_onehot_network(x, "path_dependence()")
  .refuse_binary_matrix(x, "path_dependence()")
  if (is.data.frame(x)) {
    .refuse_undeclared_roles(x, "path_dependence()", via_network = TRUE,
                             whole_frame = TRUE)
  }
  data <- .pd_extract_data(x)

  ## Counts at order k (full context + next) and at order 1 (last state + next)
  c_k <- .pd_kgram_counts(data, k = order + 1L)
  c_1 <- .pd_kgram_counts(data, k = 2L)

  if (length(c_k) == 0L)
    stop(errorCondition(
      paste0("No complete order-", order, " contexts found in data."),
      class = "nestimate_bad_input", call = NULL))

  ## Order-1 conditional table P_1(next | last)
  parts1 <- .pd_split_kgram(names(c_1), k = 2L)
  o1_df  <- data.frame(
    last = parts1$context, next_state = parts1$next_state,
    count = as.integer(c_1), stringsAsFactors = FALSE
  )
  o1_totals <- tapply(o1_df$count, o1_df$last, sum)
  o1_df$p   <- o1_df$count / o1_totals[o1_df$last]

  states <- sort(unique(c(parts1$context, parts1$next_state)))

  build_p1 <- function(last_state) {
    out <- setNames(numeric(length(states)), states)
    sub <- o1_df[o1_df$last == last_state, ]
    out[sub$next_state] <- sub$p
    out
  }
  p1_cache <- vapply(states, build_p1, numeric(length(states)))
  ## p1_cache is states x states matrix: rows = next_state, cols = last_state

  ## Order-k conditional table
  parts_k <- .pd_split_kgram(names(c_k), k = order + 1L)
  ok_df   <- data.frame(
    context = parts_k$context, next_state = parts_k$next_state,
    count   = as.integer(c_k), stringsAsFactors = FALSE
  )
  ok_df$last_state <- vapply(strsplit(ok_df$context, " -> ", fixed = TRUE),
                             function(p) p[length(p)], character(1))
  ok_totals <- tapply(ok_df$count, ok_df$context, sum)
  ok_df$p_k <- ok_df$count / ok_totals[ok_df$context]

  ## Aggregate per context, dropping low-count contexts
  ctxs <- names(ok_totals)[ok_totals >= min_count]

  ctx_rows <- do.call(rbind, lapply(ctxs, function(c) {
    sub  <- ok_df[ok_df$context == c, ]
    last <- sub$last_state[1L]
    p_k  <- setNames(numeric(length(states)), states)
    p_k[sub$next_state] <- sub$p_k
    p_1  <- p1_cache[, last]

    H1 <- .pd_entropy(p_1, base = base)
    Hk <- .pd_entropy(p_k, base = base)
    KL <- .pd_kl(p_k, p_1, base = base)

    top_o1 <- states[which.max(p_1)]
    top_ok <- states[which.max(p_k)]

    data.frame(
      context  = c,
      n        = as.integer(ok_totals[c]),
      H_order1 = H1,
      H_orderk = Hk,
      H_drop   = H1 - Hk,
      KL       = KL,
      top_o1   = top_o1,
      top_ok   = top_ok,
      flips    = top_o1 != top_ok,
      stringsAsFactors = FALSE
    )
  }))
  ctx_rows <- ctx_rows[order(-ctx_rows$KL), ]
  rownames(ctx_rows) <- NULL

  finite_KL <- is.finite(ctx_rows$KL)
  total_n   <- sum(ctx_rows$n[finite_KL])
  KL_weighted     <- if (total_n > 0)
    sum(ctx_rows$KL[finite_KL]     * ctx_rows$n[finite_KL]) / total_n else 0
  H_drop_weighted <- if (total_n > 0)
    sum(ctx_rows$H_drop[finite_KL] * ctx_rows$n[finite_KL]) / total_n else 0

  structure(
    list(
      contexts   = ctx_rows,
      chain      = list(
        KL_weighted     = KL_weighted,
        H_drop_weighted = H_drop_weighted,
        n_contexts      = nrow(ctx_rows),
        n_flips         = sum(ctx_rows$flips)
      ),
      order      = order,
      base       = base,
      min_count  = min_count,
      states     = states
    ),
    class = "net_path_dependence"
  )
}


#' @rdname path_dependence
#' @param top For the `print()` and `plot()` methods: how many contexts to
#'   show, highest KL first. Default `10` for `print()` and `15` for
#'   `plot()`.
#' @param digits For the `print()` methods: decimal places of the
#'   entropies, the entropy drop and the KL, in the table and in the
#'   chain-level summaries. Default `3`.
#' @param ... For the `print()`, `summary()` and `plot()` methods: ignored.
#' @param object For the `summary()` method: a `net_path_dependence`
#'   object.
#' @param title For the `plot()` method: the plot title. Default `NULL`,
#'   which builds `"Path dependence: order k vs order 1"` from the fitted
#'   order.
#' @export
print.net_path_dependence <- function(x, top = 10L, digits = 3L, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Path Dependence (order %d vs order 1, %s)\n\n",
              x$order, unit))
  cat(sprintf("Contexts: %d (min_count = %d).  Modal-prediction flips: %d.\n",
              x$chain$n_contexts, x$min_count, x$chain$n_flips))
  cat(sprintf("Chain-level KL_weighted      = %.*f %s\n",
              digits, x$chain$KL_weighted, unit))
  cat(sprintf("Chain-level H_drop_weighted  = %.*f %s\n\n",
              digits, x$chain$H_drop_weighted, unit))

  n_top <- min(top, nrow(x$contexts))
  cat(sprintf("Top %d contexts by KL:\n", n_top))
  tbl <- x$contexts[seq_len(n_top), ]
  tbl$H_order1 <- round(tbl$H_order1, digits)
  tbl$H_orderk <- round(tbl$H_orderk, digits)
  tbl$H_drop   <- round(tbl$H_drop,   digits)
  tbl$KL       <- round(tbl$KL,       digits)
  print(tbl, row.names = FALSE)
  invisible(x)
}


#' @rdname path_dependence
#' @export
summary.net_path_dependence <- function(object, ...) {
  structure(
    list(
      contexts   = object$contexts,
      chain      = object$chain,
      order      = object$order,
      base       = object$base,
      min_count  = object$min_count,
      n_states   = length(object$states)
    ),
    class = "summary.net_path_dependence"
  )
}

#' @rdname path_dependence
#' @export
print.summary.net_path_dependence <- function(x, digits = 3L, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Path Dependence Summary (order %d vs order 1, %s)\n\n",
              x$order, unit))
  cat(sprintf("All %d contexts (sorted by KL):\n", nrow(x$contexts)))
  tbl <- x$contexts
  tbl$H_order1 <- round(tbl$H_order1, digits)
  tbl$H_orderk <- round(tbl$H_orderk, digits)
  tbl$H_drop   <- round(tbl$H_drop,   digits)
  tbl$KL       <- round(tbl$KL,       digits)
  print(tbl, row.names = FALSE)

  cat("\nChain-level:\n")
  ch <- data.frame(
    KL_weighted     = round(x$chain$KL_weighted,     digits),
    H_drop_weighted = round(x$chain$H_drop_weighted, digits),
    n_contexts      = x$chain$n_contexts,
    n_flips         = x$chain$n_flips
  )
  print(ch, row.names = FALSE)
  invisible(x)
}


#' @rdname path_dependence
#' @export
plot.net_path_dependence <- function(x,
                                     top = 15L,
                                     title = NULL,
                                     ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"
  if (is.null(title))
    title <- sprintf("Path dependence: order %d vs order 1", x$order)

  df <- x$contexts
  df <- df[seq_len(min(top, nrow(df))), ]
  df$context <- factor(df$context, levels = rev(df$context))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$KL, y = .data$context)) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = .data$context),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = .data$n,
                                     colour = .data$flips)) +
    ggplot2::geom_text(
      data = df[df$flips, ],
      ggplot2::aes(label = sprintf("%s -> %s", .data$top_o1, .data$top_ok)),
      hjust = -0.15, size = 3.2, colour = "#D55E00") +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
      labels = c("modal next-state same",
                 "modal next-state changes"),
      name = NULL) +
    ggplot2::scale_size_continuous(name = "context count") +
    ggplot2::labs(
      x = sprintf("KL(P_order%d || P_order1)  (%s)", x$order, unit),
      y = NULL,
      title = title,
      subtitle = sprintf(
        "Chain-level weighted KL = %.3f %s; %d / %d contexts flip the modal next state",
        x$chain$KL_weighted, unit, x$chain$n_flips, x$chain$n_contexts)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold"),
      legend.position    = "bottom"
    )
}
