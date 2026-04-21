# ---- Equivalence tests: Nestimate vs tna across windowing, grouping, methods --
#
# Per-edge comparison of every method × window_size × mode × grouping
# combination against the corresponding tna:: function.
#
#   build_network(method="relative")     vs tna::tna()
#   build_network(method="frequency")    vs tna::ftna()
#   build_network(method="co_occurrence") vs tna::ctna()
#   build_network(method="attention")    vs tna::atna()
#
# Crossed with:
#   window_size: 1 (default), 2, 3, 5
#   mode: "non-overlapping", "overlapping"
#   grouping: ungrouped, 2-group, 3-group
#   lambda (attention only): 0.5, 1, 2, 5
#
# Requires NESTIMATE_EQUIV_TESTS=true.  Skipped on CRAN.

# ---- Config ------------------------------------------------------------------

set.seed(9999)
TOL <- 1e-10
N_BASE <- 30L  # base configs per combination

base_configs <- lapply(seq_len(N_BASE), function(i) {
  list(
    n_actors   = sample(c(15, 20, 30), 1),
    n_states   = sample(c(3, 5, 7), 1),
    seq_length = sample(c(12, 20, 30), 1),
    seed       = sample.int(100000, 1)
  )
})

window_sizes <- c(1L, 2L, 3L, 5L)
modes <- c("non-overlapping", "overlapping")
lambdas <- c(0.5, 1, 2, 5)
n_groups_list <- c(1L, 2L, 3L)

# ---- Helpers -----------------------------------------------------------------

#' Per-edge comparison: one expect_equal per cell
.edge_check <- function(ours, ref, tol, label) {
  ref[is.na(ref)] <- 0
  ours[is.na(ours)] <- 0
  states <- sort(intersect(rownames(ours), rownames(ref)))
  if (length(states) < 2L) return(0L)
  a <- ours[states, states, drop = FALSE]
  b <- ref[states, states, drop = FALSE]
  n <- length(states)
  n_checked <- 0L
  invisible(lapply(seq_len(n), function(i) {
    lapply(seq_len(n), function(j) {
      n_checked <<- n_checked + 1L
      expect_equal(a[i, j], b[i, j], tolerance = tol,
                   label = sprintf("%s [%s -> %s]", label, states[i], states[j]))
    })
  }))
  n_checked
}

.cfg_label <- function(cfg, i, ws, mode, extra = "") {
  sprintf("cfg%d(a=%d,s=%d,l=%d) ws=%d %s%s",
          i, cfg$n_actors, cfg$n_states, cfg$seq_length, ws, mode, extra)
}

.add_groups <- function(data, n_groups) {
  if (n_groups <= 1L) return(data)
  n <- nrow(data)
  grp <- rep(paste0("G", seq_len(n_groups)), length.out = n)
  cbind(group = grp, data, stringsAsFactors = FALSE)
}

report <- equiv_report()


# ---- 1. RELATIVE: tna::tna() × window × mode × groups ----------------------

test_that("relative: per-edge vs tna::tna() across window/mode/group combos", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_BASE), function(i) {
    cfg  <- base_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    invisible(lapply(window_sizes, function(ws) {
      # skip overlapping ws=1 (identical to non-overlapping)
      ms <- if (ws == 1L) "non-overlapping" else modes
      invisible(lapply(ms, function(mode) {
        invisible(lapply(n_groups_list, function(ng) {
          gdata <- .add_groups(data, ng)
          params <- if (ws > 1L) list(window_size = ws, mode = mode) else list()
          lbl <- .cfg_label(cfg, i, ws, mode, sprintf(" g=%d", ng))

          if (ng == 1L) {
            ours <- build_network(gdata, method = "relative", params = params)$weights
            ref  <- tna::tna(gdata, params = params)$weights
            n <- .edge_check(ours, ref, TOL, sprintf("relative %s", lbl))
            ref[is.na(ref)] <- 0
            states <- sort(intersect(rownames(ours), rownames(ref)))
            delta <- abs(ours[states, states] - ref[states, states])
            report$log("relative", lbl, n, sum(delta >= TOL),
                       max(delta), mean(delta), median(delta),
                       quantile(delta, 0.95, names = FALSE), "tna::tna()")
          } else {
            ours <- build_network(gdata, method = "relative", group = "group", params = params)
            ref  <- tna::group_tna(gdata, group = "group", params = params)
            invisible(lapply(names(ref), function(g) {
              r <- ref[[g]]$weights
              o <- ours[[g]]$weights
              .edge_check(o, r, TOL, sprintf("relative %s group=%s", lbl, g))
            }))
          }
        }))
      }))
    }))
  }))
})


# ---- 2. FREQUENCY: tna::ftna() × window × mode × groups --------------------

test_that("frequency: per-edge vs tna::ftna() across window/mode/group combos", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_BASE), function(i) {
    cfg  <- base_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    invisible(lapply(window_sizes, function(ws) {
      ms <- if (ws == 1L) "non-overlapping" else modes
      invisible(lapply(ms, function(mode) {
        invisible(lapply(n_groups_list, function(ng) {
          gdata <- .add_groups(data, ng)
          params <- if (ws > 1L) list(window_size = ws, mode = mode) else list()
          lbl <- .cfg_label(cfg, i, ws, mode, sprintf(" g=%d", ng))

          if (ng == 1L) {
            ours <- build_network(gdata, method = "frequency", params = params)$weights
            ref  <- tna::ftna(gdata, params = params)$weights
            n <- .edge_check(ours, ref, TOL, sprintf("frequency %s", lbl))
            ref[is.na(ref)] <- 0
            states <- sort(intersect(rownames(ours), rownames(ref)))
            delta <- abs(ours[states, states] - ref[states, states])
            report$log("frequency", lbl, n, sum(delta >= TOL),
                       max(delta), mean(delta), median(delta),
                       quantile(delta, 0.95, names = FALSE), "tna::ftna()")
          } else {
            ours <- build_network(gdata, method = "frequency", group = "group", params = params)
            ref  <- tna::group_ftna(gdata, group = "group", params = params)
            invisible(lapply(names(ref), function(g) {
              .edge_check(ours[[g]]$weights, ref[[g]]$weights, TOL,
                          sprintf("frequency %s group=%s", lbl, g))
            }))
          }
        }))
      }))
    }))
  }))
})


# ---- 3. CO-OCCURRENCE: tna::ctna() × window × mode × groups ----------------

test_that("co_occurrence: per-edge vs tna::ctna() across window/mode/group combos", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_BASE), function(i) {
    cfg  <- base_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    invisible(lapply(window_sizes, function(ws) {
      ms <- if (ws == 1L) "non-overlapping" else modes
      invisible(lapply(ms, function(mode) {
        invisible(lapply(n_groups_list, function(ng) {
          gdata <- .add_groups(data, ng)
          params <- if (ws > 1L) list(window_size = ws, mode = mode) else list()
          lbl <- .cfg_label(cfg, i, ws, mode, sprintf(" g=%d", ng))

          if (ng == 1L) {
            ours <- build_network(gdata, method = "co_occurrence", params = params)$weights
            ref  <- tna::ctna(gdata, params = params)$weights
            n <- .edge_check(ours, ref, TOL, sprintf("co_occurrence %s", lbl))
            ref[is.na(ref)] <- 0
            states <- sort(intersect(rownames(ours), rownames(ref)))
            delta <- abs(ours[states, states] - ref[states, states])
            report$log("co_occurrence", lbl, n, sum(delta >= TOL),
                       max(delta), mean(delta), median(delta),
                       quantile(delta, 0.95, names = FALSE), "tna::ctna()")
          } else {
            ours <- build_network(gdata, method = "co_occurrence", group = "group", params = params)
            ref  <- tna::group_ctna(gdata, group = "group", params = params)
            invisible(lapply(names(ref), function(g) {
              .edge_check(ours[[g]]$weights, ref[[g]]$weights, TOL,
                          sprintf("co_occurrence %s group=%s", lbl, g))
            }))
          }
        }))
      }))
    }))
  }))
})


# ---- 4. ATTENTION: tna::atna() × lambda × groups ----------------------------

test_that("attention: per-edge vs tna::atna() across lambda/group combos", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(N_BASE), function(i) {
    cfg  <- base_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)

    invisible(lapply(lambdas, function(lam) {
      invisible(lapply(n_groups_list, function(ng) {
        gdata <- .add_groups(data, ng)
        lbl <- sprintf("cfg%d(a=%d,s=%d,l=%d) lam=%.1f g=%d",
                       i, cfg$n_actors, cfg$n_states, cfg$seq_length, lam, ng)

        if (ng == 1L) {
          ours <- build_network(gdata, method = "attention",
                                params = list(lambda = lam))$weights
          ref  <- tna::atna(gdata, params = list(lambda = lam))$weights
          n <- .edge_check(ours, ref, TOL, sprintf("attention %s", lbl))
          ref[is.na(ref)] <- 0
          states <- sort(intersect(rownames(ours), rownames(ref)))
          delta <- abs(ours[states, states] - ref[states, states])
          report$log("attention", lbl, n, sum(delta >= TOL),
                     max(delta), mean(delta), median(delta),
                     quantile(delta, 0.95, names = FALSE), "tna::atna()")
        } else {
          ours <- build_network(gdata, method = "attention", group = "group",
                                params = list(lambda = lam))
          ref  <- tna::group_atna(gdata, group = "group",
                                  params = list(lambda = lam))
          invisible(lapply(names(ref), function(g) {
            .edge_check(ours[[g]]$weights, ref[[g]]$weights, TOL,
                        sprintf("attention %s group=%s", lbl, g))
          }))
        }
      }))
    }))
  }))
})


# ---- 5. Long format: build_network(action=, actor=) vs tna(wide) -----------
#
# Convert group_regulation to long format, then verify that
# build_network(long, action=, actor=) produces the same weights as
# tna::tna(wide). Also test with generated data across 20 configs
# varying actor count, state count, and sequence length.

test_that("long format: build_network(action=,actor=) matches tna(wide) on group_regulation", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  data(group_regulation, package = "tna")

  # Convert to long
  nr <- nrow(group_regulation)
  long_df <- do.call(rbind, lapply(seq_len(nr), function(i) {
    vals <- as.character(group_regulation[i, ])
    ok <- !is.na(vals) & vals != ""
    if (!any(ok)) return(NULL)
    data.frame(actor = i, time = which(ok), action = vals[ok],
               stringsAsFactors = FALSE)
  }))

  # Reference: tna on wide
  tna_net <- tna::tna(group_regulation)

  # Nestimate long path
  nest_long <- build_network(long_df, method = "relative",
                             action = "action", actor = "actor")

  .edge_check(nest_long$weights, tna_net$weights, TOL,
              "group_regulation long vs tna wide")
})

test_that("long format: frequency counts match tna on group_regulation", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  data(group_regulation, package = "tna")

  nr <- nrow(group_regulation)
  long_df <- do.call(rbind, lapply(seq_len(nr), function(i) {
    vals <- as.character(group_regulation[i, ])
    ok <- !is.na(vals) & vals != ""
    if (!any(ok)) return(NULL)
    data.frame(actor = i, time = which(ok), action = vals[ok],
               stringsAsFactors = FALSE)
  }))

  tna_net  <- tna::ftna(group_regulation)
  nest_net <- build_network(long_df, method = "frequency",
                            action = "action", actor = "actor")

  .edge_check(nest_net$weights, tna_net$weights, TOL,
              "group_regulation freq long vs tna wide")
})

test_that("long format: actor/time across 20 configs with varying actors/states/lengths", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(20L), function(i) {
    cfg <- base_configs[[i]]
    set.seed(cfg$seed)
    n_act <- cfg$n_actors
    states <- LETTERS[seq_len(cfg$n_states)]
    sl <- cfg$seq_length

    # Generate long data
    long_df <- do.call(rbind, lapply(seq_len(n_act), function(a) {
      data.frame(actor = a, time = seq_len(sl),
                 action = sample(states, sl, replace = TRUE),
                 stringsAsFactors = FALSE)
    }))

    # Reference: tna on wide (gold standard)
    wide_df <- do.call(rbind, lapply(split(long_df, long_df$actor), function(sub) {
      sub <- sub[order(sub$time), ]
      t(sub$action)
    }))
    colnames(wide_df) <- paste0("T", seq_len(ncol(wide_df)))
    tna_net <- tna::tna(as.data.frame(wide_df, stringsAsFactors = FALSE))

    # Nestimate long path
    nest_net <- build_network(long_df, method = "relative",
                              action = "action", actor = "actor")

    lbl <- sprintf("long cfg%d(a=%d,s=%d,l=%d)", i, n_act, cfg$n_states, sl)
    .edge_check(nest_net$weights, tna_net$weights, TOL, lbl)
  }))
})


# ---- 6. Unequal sequence lengths per actor ---------------------------------

test_that("long format: unequal lengths per actor match tna across 20 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(20L), function(i) {
    cfg <- base_configs[[i]]
    set.seed(cfg$seed)
    n_act <- cfg$n_actors
    states <- LETTERS[seq_len(cfg$n_states)]

    # Each actor gets a DIFFERENT sequence length
    long_df <- do.call(rbind, lapply(seq_len(n_act), function(a) {
      sl <- sample(8:25, 1)
      data.frame(actor = a, time = seq_len(sl),
                 action = sample(states, sl, replace = TRUE),
                 stringsAsFactors = FALSE)
    }))

    # Build wide (pad short sequences with NA)
    wide_list <- lapply(split(long_df, long_df$actor), function(sub) {
      sub <- sub[order(sub$time), ]
      sub$action
    })
    max_len <- max(vapply(wide_list, length, integer(1)))
    wide_mat <- do.call(rbind, lapply(wide_list, function(s) {
      c(s, rep(NA_character_, max_len - length(s)))
    }))
    colnames(wide_mat) <- paste0("T", seq_len(max_len))
    wide_df <- as.data.frame(wide_mat, stringsAsFactors = FALSE)

    tna_net  <- tna::tna(wide_df)
    nest_net <- build_network(long_df, method = "relative",
                              action = "action", actor = "actor")

    lbl <- sprintf("unequal cfg%d(a=%d,s=%d)", i, n_act, cfg$n_states)
    .edge_check(nest_net$weights, tna_net$weights, TOL, lbl)
  }))
})


# ---- 7. Multiple sessions per actor (time gaps) ----------------------------

test_that("long format: multi-session per actor matches tna across 10 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(10L), function(i) {
    cfg <- base_configs[[i]]
    set.seed(cfg$seed)
    n_act <- min(cfg$n_actors, 15)
    states <- LETTERS[seq_len(min(cfg$n_states, 5))]

    # 2-3 sessions per actor, separated by >900s gaps
    long_df <- do.call(rbind, lapply(seq_len(n_act), function(a) {
      n_sess <- sample(2:3, 1)
      do.call(rbind, lapply(seq_len(n_sess), function(s) {
        sl <- sample(5:10, 1)
        t_offset <- (s - 1L) * 2000L  # 2000s gap between sessions
        data.frame(actor = a, time = seq_len(sl) + t_offset,
                   action = sample(states, sl, replace = TRUE),
                   stringsAsFactors = FALSE)
      }))
    }))

    # tna: prepare_data splits sessions at time_threshold=900
    tna_wide <- tna::prepare_data(long_df, actor = "actor",
                                  action = "action", time = "time")
    tna_net <- tna::tna(tna_wide)

    # Nestimate: manually split sessions to match tna's behavior
    wide_rows <- do.call(rbind, lapply(split(long_df, long_df$actor), function(sub) {
      sub <- sub[order(sub$time), ]
      diffs <- diff(sub$time)
      breaks <- which(diffs > 900)
      if (length(breaks) == 0) return(t(sub$action))
      starts <- c(1L, breaks + 1L)
      ends <- c(breaks, nrow(sub))
      do.call(rbind, lapply(seq_along(starts), function(s) {
        t(sub$action[starts[s]:ends[s]])
      }))
    }))
    colnames(wide_rows) <- paste0("T", seq_len(ncol(wide_rows)))
    nest_net <- build_network(as.data.frame(wide_rows, stringsAsFactors = FALSE),
                              method = "relative")

    lbl <- sprintf("session cfg%d(a=%d)", i, n_act)
    .edge_check(nest_net$weights, tna_net$weights, TOL, lbl)
  }))
})


# ---- 8. Windowed + grouped -------------------------------------------------

test_that("windowed + grouped matches tna across 10 configs", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(seq_len(10L), function(i) {
    cfg <- base_configs[[i]]
    data <- simulate_sequences(cfg$n_actors, cfg$n_states, cfg$seq_length, cfg$seed)
    ws <- sample(c(2L, 3L), 1)

    gdata <- cbind(group = rep(c("G1", "G2"), length.out = nrow(data)),
                   data, stringsAsFactors = FALSE)

    tna_group  <- tna::group_tna(gdata, group = "group",
                                 params = list(window_size = ws))
    nest_group <- build_network(gdata, method = "relative", group = "group",
                                params = list(window_size = ws))

    lbl <- sprintf("grouped+ws cfg%d ws=%d", i, ws)
    invisible(lapply(names(tna_group), function(g) {
      .edge_check(nest_group[[g]]$weights, tna_group[[g]]$weights, TOL,
                  sprintf("%s group=%s", lbl, g))
    }))
  }))
})


# ---- Report ------------------------------------------------------------------

test_that("tna equivalence report is written", {
  skip_on_cran()
  skip_equiv_tests()

  report$write_csv("tna_cross_validation")
  report$write_cvs("tna_cross_validation")
  msg <- report$summary()
  message(msg)
  expect_true(length(report$rows) > 0L)
})
