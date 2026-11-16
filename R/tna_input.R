# ---- tna and group_tna models as input ----
#
# One converter for every verb. A `tna` model (tna package) carries its
# weights, its integer-coded sequence data and, as attributes, the model
# type, the scaling and the parameters it was built with. The converter
# rebuilds the network from those sequences with build_network() and trusts
# the result only when the rebuilt weights equal the model's weights. A model
# it cannot reproduce becomes a network of the weights alone, without
# sequence data: verbs that only read weights go on, verbs that resample or
# read sequences stop and say why.
#
# Nestimate reads the fields of the object; it never calls the tna package.

# tna model type -> build_network() method.
.tna_types <- c(relative = "relative", frequency = "frequency",
                `co-occurrence` = "co_occurrence", attention = "attention",
                `n-gram` = "ngram", gap = "gap", reverse = "reverse")

# The types compare_networks() has rebuilt since it first took tna models;
# its other types stay weights-only there, so its results do not move.
.tna_types_basic <- c("relative", "frequency", "co-occurrence")

# TRUE for a tna model that is not already a Nestimate network.
#' @noRd
.is_tna <- function(x) inherits(x, "tna") && !inherits(x, "netobject")

# TRUE for an input the shared step converts: a tna model, a group_tna, or a
# cograph_network that is not yet a netobject.
#' @noRd
.is_foreign_network <- function(x) {
  .is_tna(x) || inherits(x, "group_tna") ||
    (inherits(x, "cograph_network") && !inherits(x, "netobject"))
}

# The sequences of a tna model as a data frame of state names (V1, V2, ...),
# or NULL when the model keeps none.
#' @noRd
.tna_sequences <- function(x) {
  seq_data <- x$data
  if (!is.matrix(seq_data) || !is.numeric(seq_data)) return(NULL)
  alphabet <- attr(seq_data, "labels") %||% attr(seq_data, "alphabet") %||%
    x$labels %||% rownames(x$weights)
  if (is.null(alphabet)) return(NULL)
  as.data.frame(matrix(alphabet[as.integer(seq_data)], nrow = nrow(seq_data)),
                stringsAsFactors = FALSE)
}

# The per-sequence covariates of a tna model: `attr(model$data, "meta")`, kept
# by tna when the model is built from prepared data. NULL when there is none.
#' @noRd
.tna_metadata <- function(x) {
  meta <- attr(x$data, "meta")
  if (!is.data.frame(meta) || !is.matrix(x$data) ||
      nrow(meta) != nrow(x$data)) return(NULL)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  rownames(meta) <- NULL
  meta
}

# Rebuild a tna model as a netobject. Returns list(net, reason): `reason` is
# NULL when the rebuilt weights equal the model's, otherwise one clause
# saying why the model could not be reproduced, and `net` then holds the
# weights without sequence data.
#' @noRd
.tna_rebuild <- function(x, types = names(.tna_types)) {
  W <- x$weights
  if (!is.matrix(W) || !is.numeric(W) || nrow(W) != ncol(W)) {
    stop(errorCondition(
      "The tna model has no square numeric weight matrix (`$weights`).",
      class = "nestimate_bad_input", call = NULL))
  }
  labels <- x$labels %||% rownames(W)
  if (is.null(rownames(W))) dimnames(W) <- list(labels, labels)
  type <- attr(x, "type") %||% ""
  method <- if (type %in% types) unname(.tna_types[type]) else NULL
  scaling <- attr(x, "scaling") %||% character(0)
  sequences <- .tna_sequences(x)

  reason <- if (is.null(method)) {
    sprintf("its type \"%s\" is not one this function rebuilds", type)
  } else if (length(scaling) > 0L) {
    sprintf("its weights are scaled (%s)", paste(scaling, collapse = ", "))
  } else if (is.null(sequences)) {
    "it keeps no sequence data"
  }
  if (is.null(reason)) {
    # tna's own parameters carry Nestimate's names (weighted, windowed,
    # window_size, n_gram, max_gap, lambda, ...); a parameter that is not a
    # plain value (a decay function, a time matrix) cannot be passed on, and
    # the weight check below then decides.
    given <- attr(x, "params")
    given <- if (is.list(given)) {
      given[vapply(given, function(p) is.atomic(p) && is.null(dim(p)),
                   logical(1L))]
    } else list()
    given <- given[setdiff(names(given), c("", "alphabet"))]
    net <- tryCatch(
      build_network(sequences, method = method,
                    params = c(list(alphabet = labels), given)),
      error = function(e) conditionMessage(e))
    if (is.character(net)) {
      reason <- paste0("build_network() could not rebuild it: ", net)
    } else if (!identical(dim(net$weights), dim(W)) ||
               !isTRUE(all.equal(unname(net$weights[labels, labels]),
                                 unname(W), tolerance = 1e-8))) {
      reason <- "the weights rebuilt from its sequences differ from its own"
    } else {
      # The covariates tna keeps beside the sequences (one row per sequence:
      # actor, group, session), where `actor` is looked up.
      meta <- .tna_metadata(x)
      if (!is.null(meta) && is.data.frame(net$data) &&
          nrow(meta) == nrow(net$data)) {
        net$metadata <- meta
      }
      return(list(net = net, reason = NULL))
    }
  }
  list(net = .wrap_netobject(W, method = method %||% "tna", directed = TRUE,
                             inits = x$inits),
       reason = reason)
}

#' @noRd
.stop_tna_unrebuilt <- function(reason, what = NULL, member = NULL) {
  stop(errorCondition(sprintf(paste0(
    "%s needs the sequences behind the network, and the tna model%s could ",
    "not be rebuilt from its own data: %s. Build the network with ",
    "build_network() from the sequences."),
    what %||% "This function",
    if (is.null(member)) "" else sprintf(" of group '%s'", member), reason),
    class = c("nestimate_tna_unrebuilt", "nestimate_bad_input"),
    call = NULL))
}

# A tna model as a netobject. `need_data = TRUE`: stop when the model could
# not be rebuilt, since the caller reads its sequences.
#' @noRd
.tna_to_netobject <- function(x, what = NULL, need_data = FALSE,
                              types = names(.tna_types)) {
  built <- .tna_rebuild(x, types = types)
  if (need_data && !is.null(built$reason)) {
    .stop_tna_unrebuilt(built$reason, what)
  }
  built$net
}

# A group_tna as a netobject_group: one network per group, named as the
# groups are, with the grouping variable's name where grouped
# build_network() keeps it.
#' @noRd
.group_tna_to_netobjects <- function(x, what = NULL, need_data = FALSE,
                                     types = names(.tna_types)) {
  members <- unclass(x)
  group_names <- names(members) %||% paste0("Group ", seq_along(members))
  nets <- Map(function(member, nm) {
    if (!.is_tna(member)) {
      stop(errorCondition(sprintf(
        "Group '%s' of the group_tna is not a tna model.", nm),
        class = "nestimate_bad_input", call = NULL))
    }
    built <- .tna_rebuild(member, types = types)
    if (need_data && !is.null(built$reason)) {
      .stop_tna_unrebuilt(built$reason, what, member = nm)
    }
    built$net
  }, members, group_names)
  attributes(nets) <- NULL
  names(nets) <- group_names
  label <- attr(x, "label")
  attr(nets, "group_col") <- if (is.character(label) && length(label) == 1L &&
                                 !is.na(label)) label else "group"
  class(nets) <- "netobject_group"
  nets
}

# The shared step, called first by every verb that takes a network. A tna
# model becomes a netobject and a group_tna a netobject_group; anything else
# is returned as it came. `need_data = TRUE`: the caller reads sequences, so
# a model that could not be rebuilt stops here. `actor`: the caller's actor
# argument, looked up in the metadata the model keeps; a model without
# metadata cannot serve it.
#' @noRd
.tna_input <- function(x, what = NULL, need_data = FALSE, actor = NULL) {
  if (!.is_tna(x) && !inherits(x, "group_tna")) return(x)
  if (!is.null(actor)) {
    members <- if (inherits(x, "group_tna")) unclass(x) else list(x)
    has_meta <- vapply(members, function(m) !is.null(.tna_metadata(m)),
                       logical(1L))
    if (!all(has_meta)) .stop_tna_actor(what %||% "This function")
  }
  if (inherits(x, "group_tna")) {
    return(.group_tna_to_netobjects(x, what = what, need_data = need_data))
  }
  .tna_to_netobject(x, what = what, need_data = need_data)
}

# The same for the default method of a generic, which also sees a bare
# cograph_network (the generics have netobject methods only).
#' @noRd
.native_network <- function(x, what = NULL, need_data = FALSE, actor = NULL) {
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    return(.as_netobject(x))
  }
  .tna_input(x, what = what, need_data = need_data, actor = actor)
}

# `actor` names a column of the model's metadata; this model keeps none (it
# was built from sequences alone, or by a tna version that did not keep it).
#' @noRd
.stop_tna_actor <- function(what) {
  stop(errorCondition(sprintf(paste0(
    "%s: `actor` needs the actor column, and this tna model keeps no ",
    "metadata beside its sequences (`attr(model$data, \"meta\")`). Build the ",
    "tna model from tna::prepare_data(), or build the network with ",
    "build_network(), declaring `actor`."), what),
    class = c("nestimate_actor_missing", "nestimate_bad_input"),
    call = NULL))
}
