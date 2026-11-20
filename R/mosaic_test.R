# ---- mosaic_test(): the association tests behind a mosaic, as a table ----

.mosaic_test_methods <- c("chisq", "rao_scott", "donner_klar", "cmh", "gee",
                          "glmm")

#' Test the association of a mosaic with and without clustering
#'
#' @description
#' `mosaic_test()` runs the tests behind [mosaic_plot()] and
#' [mosaic_analysis()] and returns their results as a data frame. Give it
#' the same object you would draw: a network, a `netobject_group`, a
#' three-way table, or a data frame with two variables. By default it runs
#' Pearson's chi-square, which treats every count as independent, next to
#' the five tests that allow for counts nested in actors (sequences,
#' students, teams), so the table shows what the nesting changes.
#'
#' The tests are those of *Nested data* in [mosaic_analysis()]:
#' `"rao_scott"` (Rao-Scott corrected F), `"donner_klar"` (Pearson divided
#' by the cluster-size inflation factor), `"cmh"` (generalised
#' Cochran-Mantel-Haenszel with each actor as a stratum), `"gee"` (Wald
#' test of a marginal model with cluster-robust covariance) and `"glmm"`
#' (Wald test of a log-linear mixed model, needs lme4). `"chisq"` is
#' Pearson's test on the pooled table. As in the mosaics, a network's
#' source state and a group are the predictors.
#'
#' Not every test fits every design: Donner-Klar needs one variable that is
#' constant within actors (the group of a `netobject_group`), and CMH needs
#' both to vary within actors (the states of a network). When several
#' methods are asked for, a test that cannot run on the design gives a row
#' of `NA` with the reason in `note`; asked for alone, it raises the error.
#'
#' @param x A `netobject` or `htna` network (relative or frequency, with its
#'   sequences), a `netobject_group`, a three-way
#'   `table` or array (rows x columns x actors), or a data frame (with `var1`,
#'   `var2` and `actor`). A co-occurrence network raises an error of class
#'   `nestimate_mosaic_undirected`: it counts each pair of states in both
#'   mirrored cells and the pairs of a sequence share its positions, so its
#'   counts are not independent and the tests would not hold their level
#'   (under independence Pearson's test rejected 200 of 200 simulated tables
#'   at 0.05). [mosaic_plot()] still draws it. A `tna` model or a
#'   `group_tna` (package tna) is first rebuilt as a network from its own
#'   sequences (see [as_netobject()]).
#'   `actor` is read from the metadata the model keeps (`attr(model$data,
#'   "meta")`, kept when the model is built from `tna::prepare_data()`); a
#'   model without it raises `nestimate_actor_missing`.
#' @param method The tests to run, any of `"chisq"`, `"rao_scott"`,
#'   `"donner_klar"`, `"cmh"`, `"gee"` and `"glmm"`. Default: all six.
#' @param actor The actor column(s), as in [permutation()]: for a network
#'   or `netobject_group`, `NULL` (default) makes every sequence its own
#'   actor, and the name of a metadata column (such as a team) pools the
#'   sequences that share its value; for a data frame, the column(s) that
#'   identify the actors (required). Several names give nested levels,
#'   outermost first: `"glmm"` models each level, the other tests use the
#'   outermost. Not used for a three-way table, whose third dimension is
#'   the actor.
#' @param level `"test"` (default) for one row per test, `"cells"` for one
#'   row per cell per test with the cell's residual.
#' @param var1,var2 For a data frame: the names of the two variables;
#'   `var1` is the predictor.
#' @param ... Not used; a name given here raises an error of class
#'   `nestimate_bad_input`.
#'
#' @return A data frame. With `level = "test"`, one row per test: `method`,
#'   `test` (its name), `statistic` (F for Rao-Scott, chi-square otherwise),
#'   `df`, `df2` (the Rao-Scott denominator degrees of freedom), `p_value`,
#'   `n` (counts in the table), `clusters` (actors the test used),
#'   `design_effect` (Rao-Scott's mean generalised design effect, or
#'   Donner-Klar's mean inflation factor), `icc` (Donner-Klar's
#'   intracluster correlation) and `note` (why a test did not run, the caution on the GEE Wald test, else
#'   `NA`). With `level = "cells"`, one row per cell and test: `method`,
#'   the two variables (`from`/`to` for a network, `group`/`state` for a
#'   group, the variable names for a data frame, `row`/`column` for a
#'   table), `observed` and `residual` (the residual the mosaic draws for
#'   that method; `NA` where undefined). An unsupported input raises an
#'   error of class `nestimate_bad_input`.
#'
#' @references
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
#' @seealso [mosaic_plot()], [mosaic_analysis()]
#' @examples
#' data(group_regulation_long, package = "Nestimate")
#' net <- build_network(group_regulation_long, method = "frequency",
#'                      format = "long", actor = "Actor", action = "Action",
#'                      order = "Time")
#' mosaic_test(net, method = c("chisq", "rao_scott"), actor = "Group")
#' \donttest{
#' # All tests, sequences nested in teams
#' mosaic_test(net, actor = "Group")
#'
#' # Cell residuals of one test
#' head(mosaic_test(net, method = "rao_scott", actor = "Group",
#'                  level = "cells"))
#' }
#' @export
mosaic_test <- function(x, ...) UseMethod("mosaic_test")

#' @rdname mosaic_test
#' @export
mosaic_test.default <- function(x, ...) {
  # a tna model, a group_tna or a bare cograph_network: convert, dispatch again
  if (.is_foreign_network(x)) {
    x <- .native_network(x, "mosaic_test()", need_data = TRUE,
                         actor = list(...)[["actor"]])
    return(mosaic_test(x, ...))
  }
  stop(errorCondition(sprintf(paste0(
    "mosaic_test() takes a network, a netobject_group, a three-way table ",
    "or a data frame, not an object of class %s."),
    paste(class(x), collapse = "/")),
    class = "nestimate_bad_input", call = NULL))
}

#' @rdname mosaic_test
#' @export
mosaic_test.netobject <- function(x,
                                  method = c("chisq", "rao_scott",
                                             "donner_klar", "cmh", "gee",
                                             "glmm"),
                                  actor = NULL,
                                  level = c("test", "cells"), ...) {
  .mosaic_test_network(x, "netobject", method, actor, match.arg(level), ...)
}

#' @rdname mosaic_test
#' @export
mosaic_test.htna <- function(x,
                             method = c("chisq", "rao_scott", "donner_klar",
                                        "cmh", "gee", "glmm"),
                             actor = NULL, level = c("test", "cells"), ...) {
  .mosaic_test_network(x, "htna", method, actor, match.arg(level), ...)
}

#' @rdname mosaic_test
#' @export
mosaic_test.netobject_group <- function(x,
                                        method = c("chisq", "rao_scott",
                                                   "donner_klar", "cmh",
                                                   "gee", "glmm"),
                                        actor = NULL,
                                        level = c("test", "cells"), ...) {
  .mosaic_test_network(x, "netobject_group", method, actor,
                       match.arg(level), ...)
}

#' @rdname mosaic_test
#' @export
mosaic_test.array <- function(x, ...) mosaic_test.table(as.table(x), ...)

#' @rdname mosaic_test
#' @export
mosaic_test.table <- function(x,
                              method = c("chisq", "rao_scott", "donner_klar",
                                         "cmh", "gee", "glmm"),
                              level = c("test", "cells"), ...) {
  .mosaic_test_no_dots(...)
  if (length(dim(x)) != 3L) {
    stop(errorCondition(
      "mosaic_test() needs a three-way table: rows x columns x actors, e.g. table(var1, var2, actor).",
      class = "nestimate_bad_input", call = NULL))
  }
  .mosaic_test_run(unclass(x), levels = NULL, predictor = "rows",
                   method = method, level = match.arg(level),
                   var_names = c("row", "column"))
}

#' @rdname mosaic_test
#' @export
mosaic_test.data.frame <- function(x, var1, var2, actor,
                                   method = c("chisq", "rao_scott",
                                              "donner_klar", "cmh", "gee",
                                              "glmm"),
                                   level = c("test", "cells"), ...) {
  .mosaic_test_no_dots(...)
  .assert(
    "`var1` must be a single column name" =
      !missing(var1) && is.character(var1) && length(var1) == 1L,
    "`var2` must be a single column name" =
      !missing(var2) && is.character(var2) && length(var2) == 1L,
    "`actor` must name one or more columns" =
      !missing(actor) && is.character(actor) && length(actor) >= 1L &&
      !anyNA(actor) && !anyDuplicated(actor)
  )
  used <- c(var1, var2, actor)
  absent <- setdiff(used, names(x))
  if (length(absent)) stop(errorCondition(
    sprintf("Column(s) not found in data: %s.", paste(absent, collapse = ", ")),
    class = "nestimate_bad_input", call = NULL))
  if (any(c(var1, var2) %in% actor)) stop(errorCondition(
    "`actor` must not include `var1` or `var2`.",
    class = "nestimate_bad_input", call = NULL))
  keep <- Reduce(`&`, lapply(used, function(v) !is.na(x[[v]])))
  clean <- x[keep, used, drop = FALSE]
  if (nrow(clean) == 0L) stop(errorCondition(
    "No complete rows in `var1`, `var2` and `actor`.",
    class = "nestimate_bad_input", call = NULL))
  lv <- clean[, actor, drop = FALSE]
  key <- do.call(paste, c(unname(lapply(lv, as.character)), sep = "\r"))
  counts <- unclass(table(factor(clean[[var1]]), factor(clean[[var2]]),
                          factor(key, levels = unique(key))))
  .mosaic_test_run(counts, levels = lv[!duplicated(key), , drop = FALSE],
                   predictor = "rows", method = method,
                   level = match.arg(level), var_names = c(var1, var2))
}


#' @noRd
.mosaic_test_no_dots <- function(...) {
  extra <- names(list(...))
  if (length(list(...))) stop(errorCondition(sprintf(
    "mosaic_test(): unused argument(s): %s.",
    paste(ifelse(nzchar(extra %||% ""), extra, "<unnamed>"), collapse = ", ")),
    class = "nestimate_bad_input", call = NULL))
}

#' @noRd
.mosaic_test_network <- function(x, source_class, method, actor, level,
                                 ...) {
  .mosaic_test_no_dots(...)
  panels <- .mosaic_panels(x, source_class)
  if (isTRUE(panels$pairs)) .mosaic_refuse_undirected("mosaic_test()")
  tab <- panels$tabs[[1L]]
  arr <- .mosaic_cluster_array(x, source_class, actor, tab)
  var_names <- if (identical(source_class, "netobject_group")) {
    c("group", "state")
  } else {
    c("to", "from")   # mosaic rows are targets, columns sources
  }
  .mosaic_test_run(arr$counts, levels = arr$levels,
                   predictor = arr$predictor, method = method, level = level,
                   var_names = var_names)
}

# Run each requested test on one count array and bind the results. With
# several methods, a design a test cannot handle (classed nestimate_bad_input
# / nestimate_missing_package) becomes a row with the reason in `note`.
#' @noRd
.mosaic_test_run <- function(counts, levels, predictor, method, level,
                             var_names) {
  .assert(
    "`counts` must be a numeric array of three dimensions" =
      is.array(counts) && is.numeric(counts) && length(dim(counts)) == 3L,
    "`counts` must hold finite, non-negative whole numbers" =
      all(is.finite(counts)) && all(counts >= 0) &&
      all(abs(counts - round(counts)) <= 1e-8)
  )
  .assert("`method` must name tests from mosaic_test()'s list" =
            is.character(method) && length(method) >= 1L &&
            all(method %in% .mosaic_test_methods) && !anyDuplicated(method))
  tab <- rowSums(counts, dims = 2L)
  one <- function(m) {
    if (identical(m, "chisq")) return(.mosaic_test_pearson(tab))
    .nested_association(counts, m, levels = levels, predictor = predictor)
  }
  fits <- lapply(method, function(m) {
    if (length(method) == 1L) return(one(m))
    # The reason, up to its first aside: the full advice is in the error a
    # single-method call raises.
    short <- function(e) sub("\\s*[(;].*$", "", conditionMessage(e))
    tryCatch(one(m), nestimate_bad_input = short,
             nestimate_missing_package = short)
  })
  if (identical(level, "cells")) {
    cells <- expand.grid(r = rownames(tab) %||% as.character(seq_len(nrow(tab))),
                         c = colnames(tab) %||% as.character(seq_len(ncol(tab))),
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    out <- do.call(rbind, lapply(seq_along(method), function(i) {
      res <- if (is.character(fits[[i]])) rep(NA_real_, length(tab)) else
        as.vector(fits[[i]]$residuals)
      d <- data.frame(method = method[i], cells, observed = as.vector(tab),
                      residual = res, stringsAsFactors = FALSE)
      names(d)[2:3] <- var_names
      d
    }))
    # Networks read from -> to: put the source first.
    if (identical(var_names, c("to", "from"))) {
      out <- out[, c("method", "from", "to", "observed", "residual")]
    }
    rownames(out) <- NULL
    return(out)
  }
  out <- do.call(rbind, lapply(seq_along(method), function(i) {
    f <- fits[[i]]
    if (is.character(f)) {
      return(data.frame(method = method[i],
                        test = c(.nested_labels, chisq = "Pearson")[[method[i]]],
                        statistic = NA_real_, df = NA_real_, df2 = NA_real_,
                        p_value = NA_real_, n = sum(tab),
                        clusters = NA_integer_, design_effect = NA_real_,
                        icc = NA_real_, note = f, stringsAsFactors = FALSE))
    }
    t <- f$test
    data.frame(method = method[i], test = t$test, statistic = t$statistic,
               df = t$df, df2 = t$df2, p_value = t$p_value, n = sum(tab),
               clusters = as.integer(t$clusters),
               design_effect = t$design_effect, icc = t$icc,
               note = if (identical(method[i], "gee")) .mosaic_gee_note else
                 NA_character_,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# Pearson's chi-square on the pooled table, empty rows and columns dropped.
#' @noRd
.mosaic_test_pearson <- function(tab) {
  keep_r <- rowSums(tab) > 0
  keep_c <- colSums(tab) > 0
  if (sum(keep_r) < 2L || sum(keep_c) < 2L) {
    stop(errorCondition(
      "Pearson's test needs at least two non-empty rows and two non-empty columns.",
      class = "nestimate_bad_input", call = NULL))
  }
  pr <- .nest_pearson(tab[keep_r, keep_c, drop = FALSE])
  df <- (sum(keep_r) - 1) * (sum(keep_c) - 1)
  residuals <- matrix(NA_real_, nrow(tab), ncol(tab), dimnames = dimnames(tab))
  residuals[keep_r, keep_c] <- pr$stdres
  list(residuals = residuals,
       test = data.frame(test = "Pearson", statistic = pr$x2, df = df,
                         df2 = NA_real_,
                         p_value = stats::pchisq(pr$x2, df, lower.tail = FALSE),
                         clusters = NA_integer_, design_effect = NA_real_,
                         icc = NA_real_, stringsAsFactors = FALSE))
}

# The GEE Wald test is liberal unless the actors are many relative to its
# degrees of freedom (simulated type I error at alpha = 0.05, 3 x 4 tables:
# 0.28 with 30 actors, 0.11 with 100); the row says so.
.mosaic_gee_note <- paste0(
  "Wald test with a cluster-robust covariance; liberal unless the actors ",
  "are many relative to the degrees of freedom (see ?mosaic_analysis)")
