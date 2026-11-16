# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# ---- actor = c(outer, inner): restricted permutation for nested levels ----

# Teams > students > sessions. Every student draws an idiosyncratic
# transition matrix (log-scale deviations with SD `tau`), so the sessions of
# one student are correlated. The condition is a student label, half the
# students of every team in each, so teams are crossed and students pure.
# Groups never differ: H0 is true.
.sim_teams <- function(n_team = 10L, n_student = 6L, n_session = 8L,
                       len = 12L, tau = 1.5, seed = 1) {
  set.seed(seed)
  states <- c("A", "B", "C")
  base <- matrix(c(.5, .3, .2, .2, .5, .3, .3, .2, .5), 3, byrow = TRUE)
  grid <- expand.grid(student = seq_len(n_student), team = seq_len(n_team))
  rows <- Map(function(tm, st) {
    P <- exp(log(base) + matrix(rnorm(9, sd = tau), 3))
    P <- P / rowSums(P)
    seqs <- replicate(n_session, Reduce(
      function(prev, i) sample.int(3L, 1L, prob = P[prev, ]),
      seq_len(len - 1L), accumulate = TRUE, sample.int(3L, 1L)))
    data.frame(Team = paste0("t", tm), Student = paste0("s", st),
               session = rep(seq_len(n_session), each = len),
               t = rep(seq_len(len), n_session),
               Action = states[as.vector(seqs)],
               cond = if (st <= n_student / 2) "A" else "B")
  }, grid$team, grid$student)
  do.call(rbind, rows)
}

.team_net <- function(...) {
  build_network(.sim_teams(...), method = "relative",
                actor = c("Team", "Student"), session = "session",
                order = "t", action = "Action", group = "cond")
}


# ---- one level and the default are unchanged ----

# Fixtures were written on macOS arm64; Linux differs in the last bits of doubles
# (~1e-15 relative), so doubles compare to 1e-10. Integers, labels and layout stay exact.
test_that("actor = NULL and one actor column reproduce the frozen outputs", {
  ref <- readRDS(test_path("fixtures", "actor-onelevel-reference.rds"))
  strip_perm <- function(r) r[setdiff(names(r), c("x", "y"))]
  strip_boot <- function(r) r[setdiff(names(r), c("original", "model"))]
  d <- group_regulation_long
  net <- build_network(d, method = "relative", actor = "Actor",
                       action = "Action", time = "Time", group = "Achiever")
  d$Half <- ifelse(as.integer(factor(d$Actor)) %% 2L == 0L, "A", "B")
  net_w <- build_network(d, method = "relative", actor = "Actor",
                         action = "Action", time = "Time", group = "Half")
  expect_equal(strip_perm(permutation(net, iter = 50, seed = 1)[[1]]),
                   ref$perm_null, tolerance = 1e-10)
  expect_equal(
    strip_perm(permutation(net, iter = 50, seed = 1, actor = "Group")[[1]]),
    ref$perm_group, tolerance = 1e-10)
  expect_equal(
    strip_perm(permutation(net_w, iter = 50, seed = 2, actor = "Group")[[1]]),
    ref$perm_within, tolerance = 1e-10)
  expect_equal(local({
    set.seed(3)
    sequence_compare(net, sub = 2, iter = 50)
  }), ref$seqcmp_null, tolerance = 1e-10)
  expect_equal(local({
    set.seed(3)
    sequence_compare(net, sub = 2, iter = 50, actor = "Group")
  }), ref$seqcmp_group, tolerance = 1e-10)
  expect_equal(strip_boot(bootstrap_network(net[[1]], iter = 50, seed = 4)),
                   ref$boot_null, tolerance = 1e-10)
  expect_equal(
    strip_boot(bootstrap_network(net[[1]], iter = 50, seed = 4,
                                 actor = "Group")),
    ref$boot_group, tolerance = 1e-10)
  expect_equal(
    unclass(compare_networks(net, test = "permutation", iter = 50, seed = 5,
                             actor = "Group"))[c("edges", "global")],
    ref$cmp_group, tolerance = 1e-10)
})


# ---- the scheme itself (property tests on the helper) ----

# T1 (x), T2 (y), T4 (x) are pure teams; T3 is crossed: s5 (x, 2 sessions),
# s6 (y), s7 (x) are pure students, s8 has one session in each group.
.hand_ids <- list(
  team = c("T1", "T1", "T1", "T2", "T2", "T3", "T3", "T3", "T3", "T3", "T3",
           "T4"),
  student = c("s1", "s1", "s2", "s3", "s3", "s5", "s5", "s6", "s7", "s8",
              "s8", "s9")
)
.hand_x <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE,
             FALSE, TRUE)

test_that("the arrangement count is the product of the exchange sets'", {
  # teams {x, y, x}: 3; students of T3 {x, y, x}: 3; sessions of s8: 2
  expect_warning(
    design <- .nested_design(.hand_ids, .hand_x, block = c("team", "student"),
                             alpha = 0.05),
    class = "nestimate_few_actors")
  expect_equal(design$log_arrangements, log(18))
  expect_identical(design$levels$design, c("mixed", "mixed"))
  expect_identical(design$levels$n_actors, c(4L, 8L))
  expect_identical(design$levels$n_pure, c(3L, 3L))
  expect_identical(design$levels$n_crossed, c(1L, 1L))
  # 40 arrangements clear alpha = 0.05 without a warning
  expect_no_warning(.nested_design(.hand_ids, .hand_x,
                                   block = c("team", "student"), alpha = 0.1))
})

test_that("draws keep pure units whole and move students within their team", {
  design <- suppressWarnings(  # 18 arrangements: the warning is tested above
    .nested_design(.hand_ids, .hand_x, block = c("team", "student"),
                   alpha = 0.05))
  set.seed(1)
  draws <- replicate(3600L, .nested_permute(design))
  # pure teams move whole, and two of the three are x
  expect_true(all(draws[1, ] == draws[2, ] & draws[1, ] == draws[3, ]))
  expect_true(all(draws[4, ] == draws[5, ]))
  expect_true(all(draws[1, ] + draws[4, ] + draws[12, ] == 2L))
  # the crossed team is never moved as a whole: it keeps both labels
  expect_true(all(colSums(draws[6:11, ]) %in% 2:5))
  expect_true(all(colSums(!draws[6:11, ]) >= 1L))
  # inside it, pure students move whole among themselves, two of three
  # are x; s8's sessions swap
  expect_true(all(draws[6, ] == draws[7, ]))
  expect_true(all(draws[6, ] + draws[8, ] + draws[9, ] == 2L))
  expect_true(all(draws[10, ] + draws[11, ] == 1L))
  # every one of the 18 arrangements occurs, about equally often
  keys <- apply(draws, 2, function(v) paste(as.integer(v), collapse = ""))
  freq <- table(keys) / length(keys)
  expect_length(freq, 18L)
  expect_true(all(abs(freq - 1 / 18) < 0.02))
})

test_that("one level gives the arrangement count of the one-level scheme", {
  ids <- c("p1", "p1", "p2", "p3", "p3", "p4", "p4", "p4", "p5", "p5")
  is_x <- c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE)
  one <- suppressWarnings(.block_design(ids, is_x, block = "b", alpha = 0.05))
  gen <- suppressWarnings(.nested_design(list(ids), is_x, block = "b",
                                         alpha = 0.05))
  expect_equal(gen$log_arrangements, one$log_arrangements)
  set.seed(2)
  support_one <- unique(apply(replicate(2000L, .block_permute(one)), 2,
                              function(v) paste(as.integer(v), collapse = "")))
  support_gen <- unique(apply(replicate(2000L, .nested_permute(gen)), 2,
                              function(v) paste(as.integer(v), collapse = "")))
  expect_setequal(support_gen, support_one)
})

test_that("pure top-level units leave the inner level unopened", {
  team <- rep(c("T1", "T2", "T3", "T4"), each = 3L)
  student <- paste0("s", seq_along(team))
  is_x <- rep(c(TRUE, FALSE, TRUE, FALSE), each = 3L)
  design <- .nested_design(list(team, student), is_x,
                           block = c("team", "student"), alpha = 0.2)
  expect_identical(design$levels$design, c("between", "not opened"))
  set.seed(3)
  draws <- replicate(500L, .nested_permute(design))
  by_team <- split(seq_along(team), team)
  expect_true(all(vapply(by_team, function(r) {
    all(apply(draws[r, , drop = FALSE], 2, function(v) length(unique(v)) == 1L))
  }, logical(1L))))
  expect_true(all(colSums(draws) == 6L))
})

test_that("several labels and missing groups follow the same recursion", {
  team <- c("T1", "T1", "T2", "T2", "T3", "T3", "T3", "T4", "T4")
  student <- c("a", "b", "c", "d", "e", "e", "f", "g", "h")
  grp <- c("g1", "g1", "g2", "g2", "g1", "g3", "g2", NA, "g3")
  design <- suppressWarnings(
    .nested_design(list(team, student), grp, block = c("team", "student"),
                   alpha = 0.05))
  set.seed(4)
  draws <- replicate(500L, .nested_permute(design))
  # the missing group stays missing, and nothing else becomes missing
  expect_true(all(is.na(draws[8, ])))
  expect_false(anyNA(draws[-8, ]))
  # pure teams (T1, T2) move whole; T3 keeps its multiset of labels
  expect_true(all(draws[1, ] == draws[2, ] & draws[3, ] == draws[4, ]))
  expect_true(all(apply(draws[5:7, ], 2, function(v) {
    identical(sort(v), sort(grp[5:7]))
  })))
  # T4 is pure among its labelled rows and moves with T1 and T2
  expect_setequal(unique(draws[9, ]), c("g1", "g2", "g3"))
})


# ---- the verbs ----

test_that("permutation() reports the design per level", {
  net <- .team_net(seed = 5)
  res <- permutation(net, iter = 19L, seed = 1,
                     actor = c("Team", "Student"))[[1]]
  expect_identical(res$actor, c("Team", "Student"))
  expect_identical(res$n_actors, 60L)
  levels <- res$actor_levels
  expect_identical(levels$level, c("Team", "Student"))
  expect_identical(levels$design, c("within", "between"))
  expect_identical(levels$n_actors, c(10L, 60L))
  expect_identical(res$clustering$design,
                   "Team: within (10 actors), Student: between (60 actors)")
  expect_output(print(res), "Actor: Team > Student \\(nested\\)")
  expect_output(print(res), "Nesting in Student \\(innermost level\\)")
  # the observed differences do not depend on the shuffle
  plain <- permutation(net, iter = 19L, seed = 1)[[1]]
  expect_identical(res$diff, plain$diff)
})

test_that("sequence_compare() accepts nested actor levels", {
  net <- .team_net(seed = 6)
  set.seed(1)
  res <- sequence_compare(net, sub = 2, iter = 19L,
                          actor = c("Team", "Student"))
  expect_identical(res$params$actor, c("Team", "Student"))
  expect_identical(res$params$actor_levels$design, c("within", "between"))
  expect_output(print(res), "actor: Team: within \\(10 actors\\)")
  plain <- sequence_compare(net, sub = 2, iter = 19L)
  expect_identical(res$patterns$freq_A[order(res$patterns$pattern)],
                   plain$patterns$freq_A[order(plain$patterns$pattern)])
})

test_that("bootstrap_network() resamples the outermost level only", {
  net <- build_network(.sim_teams(seed = 7), method = "relative",
                       actor = c("Team", "Student"), session = "session",
                       order = "t", action = "Action")
  res <- bootstrap_network(net, iter = 20L, seed = 1,
                           actor = c("Team", "Student"))
  one <- bootstrap_network(net, iter = 20L, seed = 1, actor = "Team")
  expect_identical(res$n_actors, 10L)
  expect_identical(res$actor, c("Team", "Student"))
  # resampling the outermost level is the one-level team bootstrap
  expect_identical(res$mean, one$mean)
  expect_identical(res$clustering, one$clustering)
  expect_output(print(res), "resampled: Team \\(10 actors, outermost level\\)")
})

test_that("compare_networks() and permutation_diagnostics() pass levels on", {
  net <- .team_net(seed = 8)
  cmp <- compare_networks(net, test = "permutation", iter = 19L, seed = 1,
                          actor = c("Team", "Student"))
  expect_identical(cmp$actor, c("Team", "Student"))
  expect_output(print(cmp), "actor = Team > Student")
  diag <- permutation_diagnostics(net, actor = c("Team", "Student"),
                                  iter = 19L, seed = 1)
  expect_identical(diag$design,
                   "Team: within (10 actors), Student: between (60 actors)")
  expect_identical(diag$n_actors, 60L)
})


# ---- error and warning paths, by class ----

test_that("bad nested actor arguments raise classed errors", {
  net <- .team_net(seed = 9)
  expect_error(permutation(net, iter = 5L, actor = c("Team", "Team")),
               class = "nestimate_bad_actor")
  expect_error(permutation(net, iter = 5L, actor = c("Team", NA)),
               class = "nestimate_bad_actor")
  expect_error(permutation(net, iter = 5L, actor = c("Team", "")),
               class = "nestimate_bad_actor")
  expect_error(permutation(net, iter = 5L, actor = c("Team", "nope")),
               class = "nestimate_actor_missing")
  expect_error(permutation(net[[1]], net[[1]], iter = 5L, paired = TRUE,
                           actor = c("Team", "Student")),
               class = "nestimate_bad_actor")
  expect_error(sequence_compare(net, sub = 2, iter = 5L,
                                actor = c("Team", "Team")),
               class = "nestimate_bad_actor")
  expect_error(sequence_compare(net, sub = 2, iter = 5L,
                                actor = c("Team", "nope")),
               class = "nestimate_actor_missing")
  expect_error(bootstrap_network(net[[1]], iter = 5L,
                                 actor = c("Team", "nope")),
               class = "nestimate_actor_missing")
  expect_error(bootstrap_network(net[[1]], iter = 5L,
                                 actor = c("Team", "Team")),
               class = "nestimate_bad_actor")
  expect_error(compare_networks(net, test = "permutation", iter = 5L,
                                actor = c("Team", "Team")),
               class = "nestimate_bad_input")
})

test_that("too few arrangements warn under nested levels", {
  net <- build_network(.sim_teams(n_team = 2L, n_student = 2L, seed = 10),
                       method = "relative", actor = c("Team", "Student"),
                       session = "session", order = "t", action = "Action",
                       group = "cond")
  # two teams, each with one A and one B student: 2 x 2 arrangements
  expect_warning(permutation(net, iter = 5L, seed = 1,
                             actor = c("Team", "Student")),
                 class = "nestimate_few_actors")
})


# ---- calibration under H0 ----

.edge_rejection <- function(net, actor, seed) {
  res <- permutation(net, iter = 100L, seed = seed, actor = actor)[[1]]
  mean(res$summary$p_value < 0.05)
}

test_that("labels set on students: nested levels hold the size, the team alone does not", {
  rates <- vapply(seq_len(30L), function(rep) {
    net <- .team_net(seed = rep)
    c(team = .edge_rejection(net, "Team", rep),
      nested = .edge_rejection(net, c("Team", "Student"), rep))
  }, numeric(2L))
  rates <- rowMeans(rates)
  # sessions of one student shuffled apart: far too liberal
  expect_gt(rates[["team"]], 0.15)
  expect_gt(rates[["nested"]], 0.02)
  expect_lt(rates[["nested"]], 0.09)
})

test_that("one sequence per student: every scheme is near nominal", {
  # 100 of the 200 teams keep the test fast; half the students of each
  # team are relabelled A, half B, in every replicate.
  base <- build_network(group_regulation_long, method = "relative",
                        actor = "Actor", action = "Action", time = "Time")
  wide <- cbind(base$metadata[c("Group", "Actor")], base$data)
  wide <- subset(wide, Group <= 100)
  rates <- vapply(seq_len(20L), function(rep) {
    set.seed(100 + rep)
    wide$lab <- ave(seq_len(nrow(wide)), wide$Group, FUN = function(i) {
      sample(rep(c("A", "B"), length.out = length(i)))
    })
    net <- build_network(wide, method = "relative", group = "lab",
                         state_cols = names(base$data))
    c(none = .edge_rejection(net, NULL, rep),
      team = .edge_rejection(net, "Group", rep),
      nested = .edge_rejection(net, c("Group", "Actor"), rep))
  }, numeric(3L))
  rates <- rowMeans(rates)
  expect_true(all(rates > 0.02 & rates < 0.09))
})

test_that("levels in the wrong order are refused; the ICC message names the cause", {
  d <- data.frame(team = rep(paste0("t", 1:6), each = 20),
                  who = rep(paste0("p", 1:30), each = 4),
                  stringsAsFactors = FALSE)
  ids_ok <- list(d$team, d$who)
  ids_bad <- list(d$who, d$team)
  expect_silent(Nestimate:::.check_actor_levels(ids_ok, c("team", "who")))
  expect_error(Nestimate:::.check_actor_levels(ids_bad, c("who", "team")),
               class = "nestimate_bad_actor")
  # A redundant level (one person per team) adds nothing and is refused too.
  expect_error(Nestimate:::.check_actor_levels(list(d$team, d$team),
                                               c("team", "copy")),
               class = "nestimate_bad_actor")
  expect_match(Nestimate:::.format_icc(NA, NA, NA, 100L, 100L),
               "one sequence per actor")
  expect_match(Nestimate:::.format_icc(NA, NA, NA, 100L, 1L),
               "fewer than two actors")
  expect_match(Nestimate:::.format_icc(NA, NA, NA, 100L, 20L),
               "no variation in the transition shares")
})
