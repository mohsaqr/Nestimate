#' Extract Sequences in Wide Format
#'
#' @description
#' `extract_sequences()` returns the sequences behind a Nestimate input as a
#' wide data frame, one row per sequence and one column per position. It
#' reads the output of [prepare()], a long event log (with `action` set), a
#' wide data frame or character matrix, a `netobject`, `netobject_group` or
#' `cograph_network` built from sequences, and a `tna` model; networks and
#' models give back their stored, decoded states. The sequence verbs, such
#' as [sequence_metrics()] and [sequence_patterns()], read their input the
#' same way.
#'
#' Column roles are declared, never guessed: an event log needs `action`,
#' and a wide table needs its positions named as a numbered run or listed in
#' `cols`. When the positions were listed and their names are not a run,
#' the result records them, so it can be passed to another sequence verb
#' without `cols` again (not when `meta = TRUE` adds covariate columns).
#'
#' @param x The sequence data: a `nestimate_data` object from [prepare()];
#'   a long event log (set `action`); a wide data frame or character
#'   matrix; or a `netobject`, `netobject_group`, `cograph_network` or `tna`
#'   model built from sequence data. A network estimated with an
#'   association method (`cor`, `pcor`, `glasso`, `ising`, `mgm`) or from
#'   one-hot codes, a network built from a matrix (it holds no sequences),
#'   and any other input raise an error of class `nestimate_bad_input`.
#' @param cols The state columns of a wide data frame; the other columns
#'   are the sequences' covariates. Default `NULL`: the frame must be a
#'   single column or hold its positions as a numbered run of names (a
#'   shared prefix and a number, `T1`, `T2`, ..., or `1`, `2`, ...), and the
#'   columns outside the run, and numeric columns beside it, are
#'   covariates. Any other frame, or columns holding only 0 and 1 (which
#'   could be one-hot codes), then raise an error of class
#'   `nestimate_bad_input`; listing them in `cols` reads them as positions.
#'   A name in `cols` that is not in the frame raises the same class.
#'   Ignored for prepared data, event logs (where a listed name must still
#'   exist) and model objects.
#' @param actor,action,time The actor, action and timestamp columns of a
#'   long event log. Setting `action` declares `x` an event log and passes
#'   `actor` and `time` to [prepare()]. On a wide table without `action`,
#'   the columns they name must exist, or an error of class
#'   `nestimate_bad_input` is raised. The columns they name are covariates,
#'   never positions: they are left out of the states (also from a numbered
#'   run) and kept with the covariates, and naming one of them in `cols`
#'   too raises an error of class `nestimate_bad_input`. Default `NULL`.
#' @param meta Whether to put the covariate columns of each sequence
#'   (actor, session, group, metadata) before the state columns; a grouped
#'   network adds its grouping column among them. Default `FALSE`.
#' @param ... Further arguments passed to [prepare()] for a long event log,
#'   such as `time_threshold` or `order`; ignored otherwise.
#'
#' @return A data frame with one row per sequence. The state columns are
#'   character, `NA` after a sequence ends and at missing positions (empty
#'   strings and TraMineR's void `%` and missing `*` markers are read as
#'   `NA`). Their names are the position names of the input (`T1`, `T2`,
#'   ... for prepared or long data, `V1`, `V2`, ... for an unnamed matrix).
#'   With `meta = TRUE` the covariate columns come first. Input that holds
#'   no observed state raises an error of class `nestimate_bad_input`.
#'
#' @examples
#' extract_sequences(trajectories)
#'
#' extract_sequences(group_regulation_long, actor = "Actor",
#'                   action = "Action", time = "Time", meta = TRUE)
#'
#' net <- build_tna(trajectories)
#' extract_sequences(net)
#'
#' @seealso \code{\link{prepare}}, \code{\link{sequence_patterns}}
#' @export
extract_sequences <- function(x, cols = NULL, actor = NULL, action = NULL,
                              time = NULL, meta = FALSE, ...) {
  .assert("`meta` must be TRUE or FALSE" = isTRUE(meta) || isFALSE(meta))
  seqs <- .as_sequences(x, cols = cols, actor = actor, action = action,
                        time = time, ...)
  states <- as.data.frame(seqs$states, stringsAsFactors = FALSE)
  if (!meta || ncol(seqs$units) == 0L) return(.mark_positions(states))
  out <- cbind(seqs$units, states)
  rownames(out) <- NULL
  out
}


#' Resolve any sequence input to states plus a unit table
#'
#' @return A list with \code{states} (character matrix, one row per
#'   sequence), \code{units} (data.frame of per-sequence covariates, same
#'   rows; possibly zero columns) and \code{alphabet} (sorted states).
#' @noRd
# Methods whose $data are sequences of states.
.sequence_methods <- c("relative", "frequency", "co_occurrence", "attention",
                       "ngram", "gap", "reverse", "entropy")

# Name of the column that holds a grouped network's group labels: the
# grouping column, or for several the joined names ("Achiever-Course"),
# matching the joined group labels ("High-A").
.group_col_name <- function(x) {
  paste(attr(x, "group_col") %||% "group", collapse = "-")
}

.as_sequences <- function(x, cols = NULL, actor = NULL, action = NULL,
                          time = NULL, exclude = NULL, ...) {
  if (is.data.frame(x) && !inherits(x, "nestimate_data")) {
    .check_declared_columns(x, "Reading sequences", cols = cols, actor = actor,
                            action = action, time = time)
  }
  if (is.data.frame(x) && !inherits(x, "nestimate_data") && is.null(action) &&
      is.null(cols)) {
    .refuse_undeclared_roles(x, "Reading sequences", taken = c(actor, time, exclude))
  }
  if (is.data.frame(x) && !inherits(x, "nestimate_data") && !is.null(action)) {
    x <- prepare(x, actor = actor, action = action, time = time, ...)
  }
  .refuse_binary_matrix(x, "Reading sequences", takes_cols = TRUE)

  if (inherits(x, "nestimate_data")) {
    states <- x$sequence_data
    units <- x$meta_data
  } else if (inherits(x, "netobject_group")) {
    group_col <- .group_col_name(x)
    parts <- lapply(x, .as_sequences)
    width <- max(vapply(parts, \(p) ncol(p$states), integer(1)))
    states <- do.call(rbind, lapply(parts, \(p) .pad_states(p$states, width)))
    units <- do.call(rbind, Map(\(p, g) {
      u <- p$units
      u[[group_col]] <- rep(g, nrow(p$states))
      u
    }, parts, names(x)))
  } else if (inherits(x, c("netobject", "cograph_network", "tna"))) {
    .refuse_onehot_network(x, "Reading sequences")
    # an association or window network carries numeric or one-hot data, not
    # states; reading it as sequences would give meaningless "patterns"
    if (inherits(x, "netobject") && !is.null(x$method) &&
        !x$method %in% .sequence_methods) {
      stop(errorCondition(sprintf(paste(
        "`x` is a \"%s\" network, whose data are not sequences of states;",
        "build it with a transition method (\"relative\", \"frequency\",",
        "...) or pass the sequences."), x$method),
        class = "nestimate_bad_input", call = NULL))
    }
    if (is.null(x$data)) {
      stop(errorCondition(
        "`x` holds no sequence data; build it from sequences, not a matrix.",
        class = "nestimate_bad_input", call = NULL))
    }
    states <- .coerce_sequence_input(x)
    if (!is.data.frame(states)) {
      stop(errorCondition("`x` holds no sequence data.",
                          class = "nestimate_bad_input", call = NULL))
    }
    units <- if (inherits(x, "netobject")) x$metadata else NULL
  } else if (is.data.frame(x) || is.matrix(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    missing_cols <- setdiff(c(cols, exclude), names(df))
    if (length(missing_cols) > 0L) {
      stop(errorCondition(
        sprintf("column(s) not found in `x`: %s.",
                paste(missing_cols, collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    # a column declared as actor or time is a covariate, never a position
    roles <- intersect(c(actor, time), names(df))
    if (length(intersect(cols, roles)) > 0L) {
      stop(errorCondition(
        sprintf(paste("column(s) %s are declared both as positions (`cols`)",
                      "and as `actor`/`time`; a column has one role."),
                paste(intersect(cols, roles), collapse = ", ")),
        class = "nestimate_bad_input", call = NULL))
    }
    state_cols <- .build_sequence_cols(df, cols = cols,
                                       exclude = c(exclude, roles)) %||%
      setdiff(names(df), roles)
    states <- df[state_cols]
    units <- df[setdiff(names(df), state_cols)]
  } else {
    stop(errorCondition(
      paste("`x` must be prepared data, a long or wide data frame, a",
            "character matrix, a netobject, a netobject_group or a tna",
            "model."),
      class = "nestimate_bad_input", call = NULL))
  }

  states <- .states_matrix(states)
  units <- if (is.null(units)) {
    data.frame(row.names = seq_len(nrow(states)))
  } else {
    as.data.frame(units, stringsAsFactors = FALSE)
  }
  rownames(units) <- NULL
  stopifnot("sequence and covariate rows must align" =
              nrow(units) == nrow(states))
  if (nrow(states) == 0L || all(is.na(states))) {
    stop(errorCondition("`x` contains no observed states.",
                        class = "nestimate_bad_input", call = NULL))
  }
  list(states = states, units = units,
       alphabet = sort(unique(states[!is.na(states)])))
}


#' Character state matrix with void and missing markers as NA
#'
#' Same markers as \code{build_network()}: empty strings and TraMineR's void
#' (\code{\%}) and missing (\code{*}) symbols, so an \code{stslist} reads
#' with its padding removed.
#' @noRd
.states_matrix <- function(states) {
  df <- as.data.frame(states, stringsAsFactors = FALSE)
  # as.matrix() on mixed columns pads numbers with format(); convert first
  df[] <- lapply(df, as.character)
  m <- .clean_states(as.matrix(df))
  rownames(m) <- NULL
  m
}


#' Pad a state matrix with NA columns up to `width`
#' @noRd
.pad_states <- function(m, width) {
  if (ncol(m) >= width) return(m)
  pad <- matrix(NA_character_, nrow(m), width - ncol(m))
  extra <- setdiff(paste0("T", seq_len(width)), colnames(m))
  colnames(pad) <- utils::head(extra, width - ncol(m))
  cbind(m, pad)
}
