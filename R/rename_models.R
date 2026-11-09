# ---- rename_models ----
#
# A `netobject_group` is a named list whose names label each constituent
# network ("Cluster 1", "High", "A", ...). Those labels show up in every
# downstream report (print methods, plot panel titles, bootstrap output,
# permutation tests, etc.), so the group built by `build_network()` may
# need post-hoc renaming once the analyst knows what each cluster
# actually represents.
#
# `rename_models()` is the equivalent of `tna::rename_groups()`. The
# Nestimate name uses "models" because the Nestimate vocabulary refers
# to each constituent of a `netobject_group` as a model (one fitted
# network), reserving "group" for the dispatch class.

#' Rename the models of a `netobject_group`
#'
#' @description
#' `rename_models()` replaces the names of the networks in a
#' `netobject_group` and returns the group with nothing else changed. The
#' names label every later report (print methods, plot panels, bootstrap
#' and permutation output), so this is the step that turns the generic
#' labels [build_network()] gives clusters (`"Cluster 1"`, `"Cluster 2"`)
#' into meaningful ones (`"High engagement"`, `"Low engagement"`) once it is
#' known what each cluster represents. It is the Nestimate counterpart of
#' `tna::rename_groups()`.
#'
#' @param x A `netobject_group`, or an object inheriting from it. Any other
#'   object raises an error. A `net_mlvar` is accepted too, but its
#'   `print()` and `summary()` find the three layers by their names
#'   (`temporal`, `contemporaneous`, `between`): after renaming, `print()`
#'   shows none of the model's dimensions and `summary()` fails.
#' @param new_names The new names, a character vector in the order of the
#'   networks in `x`: one per network, none missing or empty, and no two the
#'   same. Anything else raises an error.
#' @return An object of the same class and length as `x`, its networks and
#'   attributes (such as `group_col`) unchanged and its names replaced by
#'   `new_names`.
#' @examples
#' grp <- build_network(group_regulation_long, method = "tna",
#'                      actor = "Actor", action = "Action", time = "Time",
#'                      group = "Achiever")
#' names(grp)
#' names(rename_models(grp, c("High achievers", "Low achievers")))
#' @export
rename_models <- function(x, new_names) {
  UseMethod("rename_models")
}

#' @rdname rename_models
#' @export
rename_models.netobject_group <- function(x, new_names) {
  .assert(
    is.character(new_names),
    length(new_names) == length(x),
    !anyNA(new_names),
    all(nzchar(new_names)),
    !anyDuplicated(new_names)
  )
  names(x) <- new_names
  x
}

#' @rdname rename_models
#' @export
rename_models.default <- function(x, new_names) {
  stop(errorCondition(
    paste0(
      "rename_models() expects a netobject_group (or subclass). Got: ",
      paste(class(x), collapse = ", "),
      "."),
    class = "nestimate_bad_input", call = NULL))
}
