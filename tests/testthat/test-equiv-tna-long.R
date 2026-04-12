# ---- Equivalence: Nestimate vs tna on bundled long-format datasets ----------
#
# Same long data, same args, both packages, per-edge comparison.
# Covers: actor, session, time, action granularity, grouping.
# 151 combos total.

TOL <- 1e-10

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

#' Mirror Nestimate auto-match: combine actor + session_id for tna
.tna_net <- function(d, meth, act, actor, time) {
  eff_session <- NULL
  cl <- tolower(names(d))
  hit <- which(cl == "session_id")
  if (length(hit) == 1L) eff_session <- names(d)[hit]
  if (is.null(eff_session)) {
    hit2 <- which(cl == "session")
    if (length(hit2) == 1L) eff_session <- names(d)[hit2]
  }
  ta <- list(data = d, action = act)
  ta$actor <- c(actor, eff_session)
  if (!is.null(time)) ta$time <- time
  tp <- suppressMessages(suppressWarnings(do.call(tna::prepare_data, ta)))
  suppressMessages(if (meth == "relative") tna::tna(tp) else tna::ftna(tp))
}

.nest_net <- function(d, meth, act, actor, time) {
  a <- list(data = d, method = meth, action = act)
  if (!is.null(actor)) a$actor <- actor
  if (!is.null(time)) a$time <- time
  suppressMessages(suppressWarnings(do.call(build_network, a)))
}

.cmp <- function(d, meth, act, actor, time, lbl) {
  n <- .nest_net(d, meth, act, actor, time)
  t <- .tna_net(d, meth, act, actor, time)
  .edge_check(n$weights, t$weights, TOL, lbl)
}


# ---- 1. human_ai: 3 actions × 2 methods × 5 arg combos = 30 ---------------

test_that("human_ai: all arg combos match tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  d <- Nestimate::human_ai
  invisible(lapply(c("code", "category", "superclass"), function(act) {
    invisible(lapply(c("relative", "frequency"), function(meth) {
      .cmp(d, meth, act, "actor", "timestamp",
           sprintf("human_ai %s %s actor+session+time", meth, act))
      .cmp(d, meth, act, "actor", NULL,
           sprintf("human_ai %s %s actor+session", meth, act))
      .cmp(d, meth, act, NULL, NULL,
           sprintf("human_ai %s %s session_only", meth, act))

      # Per-actor subsets
      invisible(lapply(unique(d$actor), function(a) {
        sub <- d[d$actor == a, ]
        .cmp(sub, meth, act, "actor", "timestamp",
             sprintf("human_ai %s %s actor=%s", meth, act, a))
      }))
    }))
  }))
})


# ---- 2. human_detailed + ai_detailed: 3 actions × 2 methods × 4 combos ----

test_that("human_detailed and ai_detailed match tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(c("human_detailed", "ai_detailed"), function(ds) {
    d <- get(ds, envir = asNamespace("Nestimate"))
    invisible(lapply(c("code", "category", "superclass"), function(act) {
      invisible(lapply(c("relative", "frequency"), function(meth) {
        .cmp(d, meth, act, "actor", "timestamp",
             sprintf("%s %s %s actor+time", ds, meth, act))
        .cmp(d, meth, act, "actor", NULL,
             sprintf("%s %s %s actor", ds, meth, act))
        .cmp(d, meth, act, NULL, NULL,
             sprintf("%s %s %s session_only", ds, meth, act))
      }))
    }))
  }))
})


# ---- 3. group_regulation_long: actor+time, actor-only ----------------------

test_that("group_regulation_long matches tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  d <- Nestimate::group_regulation_long
  invisible(lapply(c("relative", "frequency"), function(meth) {
    .cmp(d, meth, "Action", "Actor", "Time",
         sprintf("group_reg %s actor+time", meth))
    .cmp(d, meth, "Action", "Actor", NULL,
         sprintf("group_reg %s actor", meth))
  }))
})


# ---- 4. Grouped: human_ai by actor, by project -----------------------------

test_that("human_ai grouped by actor matches tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  d <- Nestimate::human_ai
  invisible(lapply(c("code", "category", "superclass"), function(act) {
    invisible(lapply(unique(d$actor), function(a) {
      sub <- d[d$actor == a, ]
      invisible(lapply(c("relative", "frequency"), function(meth) {
        .cmp(sub, meth, act, "actor", "timestamp",
             sprintf("human_ai %s %s actor=%s", meth, act, a))
      }))
    }))
  }))
})

test_that("human_ai grouped by project (10) matches tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  d <- Nestimate::human_ai
  projs <- head(sort(unique(d$project)), 10)
  invisible(lapply(projs, function(p) {
    sub <- d[d$project == p, ]
    invisible(lapply(c("category", "superclass"), function(act) {
      .cmp(sub, "relative", act, "actor", "timestamp",
           sprintf("human_ai project=%s %s", p, act))
    }))
  }))
})


# ---- 5. Grouped: group_regulation by Group, Course, Achiever ---------------

test_that("group_regulation grouped by Group/Course/Achiever matches tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  d <- Nestimate::group_regulation_long

  # By Group (first 10)
  groups <- head(sort(unique(d$Group)), 10)
  invisible(lapply(groups, function(g) {
    sub <- d[d$Group == g, ]
    invisible(lapply(c("relative", "frequency"), function(meth) {
      .cmp(sub, meth, "Action", "Actor", "Time",
           sprintf("group_reg %s Group=%s", meth, g))
    }))
  }))

  # By Course
  invisible(lapply(unique(d$Course), function(cr) {
    sub <- d[d$Course == cr, ]
    .cmp(sub, "relative", "Action", "Actor", "Time",
         sprintf("group_reg Course=%s", cr))
  }))

  # By Achiever
  invisible(lapply(unique(d$Achiever), function(ach) {
    sub <- d[d$Achiever == ach, ]
    invisible(lapply(c("relative", "frequency"), function(meth) {
      .cmp(sub, meth, "Action", "Actor", "Time",
           sprintf("group_reg %s Achiever=%s", meth, ach))
    }))
  }))
})


# ---- 6. Wide datasets direct -----------------------------------------------

test_that("wide datasets match tna", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_not_installed("tna")

  invisible(lapply(c("human_wide", "ai_wide", "srl_strategies"), function(ds) {
    d <- get(ds, envir = asNamespace("Nestimate"))
    n <- build_network(d, method = "relative")
    t <- tna::tna(d)
    .edge_check(n$weights, t$weights, TOL, ds)
  }))
})
