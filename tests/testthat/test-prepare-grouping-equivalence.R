# Randomized grouping equivalence for prepare()
#
# Two independent arms:
#
#   1. Ground truth. Each dataset is generated from a known set of sessions, so
#      the expected sequences are known without consulting any Nestimate code.
#      Comparing an optimized path against Nestimate's own internals would be
#      circular; this compares against the construction instead.
#
#   2. Legacy ordering. interaction() is called directly as the reference for
#      the row order prepare() used to produce, pinning the backward
#      compatibility that a finite same-seed bootstrap depends on.
#
# The generator sweeps actor/session column counts, column types, time on and
# off, explicit order columns, tied timestamps, thresholds, and identifiers
# containing the separators that used to collide.

N_DATASETS <- 1000L

STATES <- c("A", "B", "C", "D", "E")

# Identifier vocabularies, given separately for the actor and session roles so
# the two can form genuinely colliding pairs. The "sep" pool is built so that
# ("a | b", "c") and ("a", "b | c") both paste to "a | b | c" -- the collision
# that used to merge two distinct sessions into one.
ID_POOLS <- list(
  plain = list(actor = sprintf("id%02d", 1:8),
               session = c("q1", "q2", "q3", "q4")),
  sep   = list(actor = c("a | b", "a", "x-y", "x", "p", "p | q", "q", "r"),
               session = c("c", "b | c", "y", "-y")),
  num   = list(actor = 1:8, session = 1:4),
  blank = list(actor = c("", " ", "  ", "z", "zz", "z z", "  z", "zzz"),
               session = c("", " ", "y", "  ")),
  utf   = list(actor = c("élève", "学生", "u1", "u|1", "u", "ü", "u 1", "é"),
               session = c("1", "|1", "é", "1 "))
)

.pick_ids <- function(pool, role, n) {
  vals <- ID_POOLS[[pool]][[role]]
  unique(as.character(rep_len(vals, min(n, length(vals)))))
}

# One random configuration. Seeded by index so a failure is reproducible.
.grouping_config <- function(i) {
  set.seed(20260720L + i)
  list(
    n_actor_cols   = sample(1:3, 1L),
    n_session_cols = sample(0:2, 1L),
    n_actors       = sample(2:8, 1L),
    n_sessions     = sample(1:4, 1L),
    pool           = sample(names(ID_POOLS), 1L),
    # Time is weighted on: the gap-splitting path is what the interval
    # options govern, and a coin flip left it thinly covered.
    use_time       = sample(c(TRUE, TRUE, TRUE, FALSE), 1L),
    use_order      = sample(c(TRUE, FALSE), 1L),
    tied_times     = sample(c(TRUE, FALSE, FALSE, FALSE), 1L),
    shuffle        = sample(c(TRUE, FALSE), 1L),
    events_per     = sample(2:5, 1L),
    threshold      = sample(list(60, 900, 3600, FALSE), 1L)[[1L]],
    actor_factor   = sample(c(TRUE, FALSE), 1L),
    # Time-gap structure inside each actor/session block. "none" keeps events
    # a second apart; "split" inserts gaps past the threshold; "boundary" uses
    # a gap of exactly time_threshold, which must NOT split because the test
    # is `gaps > time_threshold`.
    gap_mode       = sample(c("none", "split", "boundary"), 1L),
    n_sub          = sample(1:3, 1L)
  )
}

# Build an event log from an explicit session list, so the expected sequences
# are known by construction rather than derived from the code under test.
.grouping_fixture <- function(cfg) {
  actors <- .pick_ids(cfg$pool, "actor", cfg$n_actors)
  sessions <- .pick_ids(cfg$pool, "session", cfg$n_sessions)
  if (cfg$n_session_cols == 0L) sessions <- "single"

  grid <- expand.grid(actor = actors, session = sessions,
                      stringsAsFactors = FALSE)

  # Gap inserted between sub-blocks. Only "split" exceeds the threshold; the
  # boundary case sits exactly on it and must stay one session.
  thr <- if (isFALSE(cfg$threshold)) 900 else cfg$threshold
  gap <- switch(cfg$gap_mode,
                none = 1, split = thr + 10, boundary = thr)
  n_sub <- if (cfg$gap_mode == "none") 1L else cfg$n_sub

  # Each block is n_sub runs of events_per events, separated by `gap`.
  blocks <- lapply(seq_len(nrow(grid)), function(g) {
    within_pos <- seq_len(cfg$events_per)
    # A sub-block spans events_per - 1 seconds, so starting the next one
    # (events_per - 1 + gap) later makes the inter-block gap exactly `gap`.
    offsets <- unlist(lapply(seq_len(n_sub) - 1L, function(b) {
      b * (cfg$events_per - 1L + gap) + within_pos
    }))
    data.frame(
      .truth_actor   = grid$actor[g],
      .truth_session = grid$session[g],
      .truth_group   = g,
      .truth_sub     = rep(seq_len(n_sub), each = cfg$events_per),
      act = sample(STATES, cfg$events_per * n_sub, replace = TRUE),
      offset = offsets,
      stringsAsFactors = FALSE
    )
  })
  ev <- do.call(rbind, blocks)

  # Blocks sit far apart so they never merge, whatever the threshold.
  ev$tm <- 1700000000 + ev$.truth_group * 1e7 + ev$offset
  if (cfg$tied_times) ev$tm <- 1700000000 + ev$.truth_group * 1e7

  # Expected sessions: gaps only split when time is used, the gap strictly
  # exceeds the threshold, splitting is enabled, and times are not all tied.
  splits <- cfg$use_time && !cfg$tied_times &&
    !isFALSE(cfg$threshold) && cfg$gap_mode == "split"
  n_expected <- nrow(grid) * if (splits) n_sub else 1L

  # Spread the identifier across the requested number of columns. Splitting a
  # single logical id over several columns must not change the partition.
  actor_cols <- paste0("a", seq_len(cfg$n_actor_cols))
  ev[[actor_cols[1L]]] <- ev$.truth_actor
  if (cfg$n_actor_cols > 1L) {
    for (j in seq_len(cfg$n_actor_cols - 1L) + 1L) {
      ev[[actor_cols[j]]] <- paste0("k", nchar(ev$.truth_actor))
    }
  }
  session_cols <- character(0)
  if (cfg$n_session_cols > 0L) {
    session_cols <- paste0("s", seq_len(cfg$n_session_cols))
    ev[[session_cols[1L]]] <- ev$.truth_session
    if (cfg$n_session_cols > 1L) {
      ev[[session_cols[2L]]] <- paste0("w", nchar(ev$.truth_session))
    }
  }

  if (cfg$actor_factor) ev[[actor_cols[1L]]] <- factor(ev[[actor_cols[1L]]])
  ev$ord <- seq_len(nrow(ev))
  if (cfg$shuffle) ev <- ev[sample(nrow(ev)), , drop = FALSE]

  list(events = ev, actor_cols = actor_cols, session_cols = session_cols,
       n_groups = nrow(grid), n_expected = n_expected)
}

.call_prepare <- function(fx, cfg) {
  args <- list(
    data = fx$events, actor = fx$actor_cols, action = "act",
    session = if (length(fx$session_cols)) fx$session_cols else NULL,
    time = if (cfg$use_time) "tm" else NULL,
    order = if (cfg$use_order) "ord" else NULL
  )
  if (cfg$use_time) args$time_threshold <- cfg$threshold
  do.call(prepare, args)
}

test_that("randomized grouping matches ground truth across configurations", {
  skip_on_cran()

  for (i in seq_len(N_DATASETS)) {
    cfg <- .grouping_config(i)
    fx  <- .grouping_fixture(cfg)
    res <- .call_prepare(fx, cfg)
    info <- sprintf("dataset %d (pool=%s, actors=%d, sessions=%d, time=%s)",
                    i, cfg$pool, cfg$n_actor_cols, cfg$n_session_cols,
                    cfg$use_time)

    # One sequence per constructed (actor, session) block: identifiers that
    # merely contain a separator must not collapse, and splitting one id over
    # several columns must not multiply groups.
    expect_identical(nrow(res$sequence_data), fx$n_expected, info = info)

    # Rowwise alignment across all three returned frames.
    expect_identical(nrow(res$meta_data), fx$n_expected, info = info)
    if (!is.null(res$time_data)) {
      expect_identical(nrow(res$time_data), fx$n_expected, info = info)
    }
    expect_identical(anyDuplicated(res$meta_data$.session_id), 0L, info = info)

    # Every event survives, and no cell is invented.
    expect_identical(sum(!is.na(as.matrix(res$sequence_data))),
                     nrow(fx$events), info = info)

    # The multiset of actions per sequence must match the blocks that produced
    # it, so no session absorbed another session's events.
    got <- sort(table(unlist(res$sequence_data), useNA = "no"))
    want <- sort(table(fx$events$act))
    expect_identical(got, want, info = info)
  }
})

test_that("grouping reproduces the legacy interaction() row order", {
  skip_on_cran()

  for (i in seq_len(N_DATASETS)) {
    cfg <- .grouping_config(i)
    fx  <- .grouping_fixture(cfg)
    ev  <- fx$events
    info <- sprintf("dataset %d", i)

    # Legacy reference, called directly rather than through Nestimate.
    actor_key <- if (length(fx$actor_cols) > 1L) {
      interaction(ev[, fx$actor_cols, drop = FALSE], sep = "-", drop = TRUE)
    } else {
      ev[[fx$actor_cols]]
    }
    legacy <- if (length(fx$session_cols)) {
      session_key <- if (length(fx$session_cols) > 1L) {
        interaction(ev[, fx$session_cols, drop = FALSE], sep = "-", drop = TRUE)
      } else {
        ev[[fx$session_cols]]
      }
      interaction(actor_key, session_key, sep = " | ", drop = TRUE)
    } else {
      as.factor(actor_key)
    }

    got <- .observed_group_id(ev[, c(fx$actor_cols, fx$session_cols),
                                 drop = FALSE])

    if (nlevels(legacy) == fx$n_groups) {
      # Legacy agreed with the construction, so the partition and its ordering
      # must be reproduced exactly. Order is what a finite same-seed bootstrap
      # indexes into, so a permutation here would break stored analyses.
      expect_identical(max(got), nlevels(legacy), info = info)
      expect_identical(order(got, seq_along(got)),
                       order(legacy, seq_along(legacy)), info = info)
      expect_identical(duplicated(got), duplicated(as.integer(legacy)),
                       info = info)
    } else {
      # Legacy pasted two distinct identifier pairs into the same label and
      # merged real sessions. Divergence here is the fix, not a regression:
      # the new grouping must match the construction instead.
      expect_lt(nlevels(legacy), fx$n_groups)
      expect_identical(max(got), fx$n_groups, info = info)
    }
  }
})

test_that("grouping is exact on a deterministic core without skipping", {
  # A small always-run subset, so CRAN still exercises the path.
  for (i in 1:25) {
    cfg <- .grouping_config(i)
    fx  <- .grouping_fixture(cfg)
    res <- .call_prepare(fx, cfg)
    expect_identical(nrow(res$sequence_data), fx$n_expected,
                     info = sprintf("dataset %d", i))
  }
})

test_that("grouping stays exact across the 32-bit overflow boundary", {
  skip_on_cran()

  # interaction() codes two columns as code(b) + nlevels(b) * code(a), which
  # exceeds .Machine$integer.max once both marginals pass sqrt(2^31 - 1). The
  # helper is exercised directly here so the sweep stays cheap; test-prepare.R
  # carries the end-to-end prepare() case.
  boundary <- as.integer(floor(sqrt(.Machine$integer.max)))
  for (k in c(1000L, boundary - 1L, boundary, boundary + 1L, 50000L, 60000L)) {
    cols <- list(a = sprintf("a%06d", seq_len(k)),
                 b = sprintf("b%06d", seq_len(k)))
    info <- sprintf("k = %d (product %.3e)", k, as.numeric(k)^2)

    expect_no_warning(ids <- .observed_group_id(cols), message = info)
    expect_false(anyNA(ids), info = info)
    # One-to-one pairs: k rows, k distinct groups, no two merged.
    expect_identical(max(ids), k, info = info)
    expect_identical(anyDuplicated(ids), 0L, info = info)
  }
})

test_that("randomized missing identifiers always fail fast", {
  skip_on_cran()

  for (i in seq_len(200L)) {
    cfg <- .grouping_config(i)
    fx  <- .grouping_fixture(cfg)
    ev  <- fx$events
    set.seed(i)

    # Blank out one identifier cell in a randomly chosen grouping column.
    key_cols <- c(fx$actor_cols, fx$session_cols)
    target <- sample(key_cols, 1L)
    ev[[target]] <- as.character(ev[[target]])
    ev[sample(nrow(ev), 1L), target] <- NA

    args <- list(
      data = ev, actor = fx$actor_cols, action = "act",
      session = if (length(fx$session_cols)) fx$session_cols else NULL,
      time = if (cfg$use_time) "tm" else NULL
    )
    # Silently changing the session count is the failure mode being guarded.
    expect_error(do.call(prepare, args),
                 "Missing values in actor/session column",
                 info = sprintf("dataset %d, column %s", i, target))
  }
})
