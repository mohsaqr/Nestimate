# sequence_compare(actor =): opt-in blocked permutation over actors.

# ---- helpers ----

# Wide sequences nested in actors, `per_actor` per actor, with actor-specific
# state probabilities and a group label that is constant within an actor and
# unrelated to the states (a null with strong actor effects).
.sca_actor_null <- function(n_actors = 16L, per_actor = 5L, len = 8L) {
  states <- c("A", "B", "C")
  actor_probs <- lapply(seq_len(n_actors), function(a) {
    w <- stats::rgamma(length(states), shape = 0.6)
    w / sum(w)
  })
  rows <- lapply(seq_len(n_actors), function(a) {
    m <- matrix(sample(states, per_actor * len, TRUE, prob = actor_probs[[a]]),
                nrow = per_actor)
    colnames(m) <- paste0("T", seq_len(len))
    data.frame(m, team = sprintf("a%02d", a), stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$grp <- ifelse(seq_len(n_actors) %% 2L == 0L, "X", "Y")[
    as.integer(sub("a", "", d$team))]
  d
}

.sca_wide <- function() {
  data.frame(
    T1 = c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"),
    T2 = c("B", "C", "A", "B", "C", "A", "C", "A", "B", "B", "A", "C"),
    T3 = c("C", "A", "B", "C", "A", "B", "A", "C", "A", "C", "B", "A"),
    grp = rep(c("x", "y", "z"), each = 4L),
    who = rep(c("p1", "p2", "p3", "p4", "p5", "p6"), each = 2L),
    stringsAsFactors = FALSE
  )
}

# ---- default unchanged: identical to the released 1.0.12 ----

# Fixtures were written on macOS arm64; Linux differs in the last bits of doubles
# (~1e-15 relative), so doubles compare to 1e-10. Integers, labels and layout stay exact.
test_that("actor = NULL preserves released outputs for labelled sequences", {
  fx <- readRDS(test_path("fixtures", "sequence-compare-default-1.0.12.rds"))
  net2 <- build_network(group_regulation_long, method = "relative",
                        actor = "Actor", action = "Action", time = "Time",
                        group = "Achiever")
  net3 <- build_network(group_regulation_long, method = "relative",
                        actor = "Actor", action = "Action", time = "Time",
                        group = "Course")
  wide <- data.frame(
    T1 = c("A", "B", "C", "A", "B", "C", "A", "B", "C", "A", "B", "C"),
    T2 = c("B", "C", "A", "B", "C", "A", "C", "A", "B", "B", "A", "C"),
    T3 = c("C", "A", "B", "C", "A", "B", "A", "C", "A", "C", "B", "A"),
    T4 = c("A", "B", "C", "A", "C", "B", "B", "A", "C", "A", "B", "C"),
    grp = c("x", "y", "x", "y", "x", "y", "x", "y", "x", "y", NA, "x"),
    stringsAsFactors = FALSE)

  set.seed(42)
  r1 <- sequence_compare(net2, sub = 2:3, iter = 50L)
  set.seed(7)
  r2 <- sequence_compare(net3, sub = 3L, iter = 30L, adjust = "holm")
  set.seed(3)
  expect_warning(
    r3 <- sequence_compare(wide, group = "grp", sub = 2L, min_freq = 1L,
                           iter = 40L),
    class = "nestimate_missing_group")
  r4 <- sequence_compare(net2, sub = 2L, test = "chisq")
  r5 <- sequence_compare(net3, sub = 2L, test = "none")

  expect_equal(r1, fx$r1, tolerance = 1e-10)
  expect_equal(r2, fx$r2, tolerance = 1e-10)
  # The old fixture shuffled NA labels into observed groups. The corrected
  # null must instead equal an analysis that excludes the unlabelled row.
  set.seed(3)
  labelled <- sequence_compare(wide[!is.na(wide$grp), ], group = "grp",
                               sub = 2L, min_freq = 1L, iter = 40L)
  expect_equal(r3, labelled, tolerance = 1e-10)
  expect_equal(r4, fx$r4, tolerance = 1e-10)
  expect_identical(r5, fx$r5)
  expect_null(r1$params$actor)
})

# ---- the blocked permutation helper ----

test_that("pure actors move as a whole and keep the label multiset", {
  ids <- rep(sprintf("a%d", 1:6), times = c(3, 2, 4, 1, 2, 3))
  grp <- rep(c("x", "y", "z", "x", "y", "z"), times = c(3, 2, 4, 1, 2, 3))
  design <- .cs_actor_design(ids, grp, "who")
  expect_identical(design$design, "between")
  expect_identical(design$n_actors, 6L)

  set.seed(11)
  perms <- replicate(200L, .cs_actor_permute(design))
  # every actor holds one label in every permutation
  labels_per_actor <- apply(perms, 2L, function(p) {
    vapply(split(p, ids), function(v) length(unique(v)), integer(1L))
  })
  expect_true(all(labels_per_actor == 1L))
  # the actors' labels are a reshuffle of the original ones
  actor_label <- apply(perms, 2L, function(p) {
    sort(vapply(split(p, ids), function(v) v[[1L]], character(1L)))
  })
  expect_true(all(actor_label == sort(c("x", "y", "z", "x", "y", "z"))))
  # and the shuffle does move them
  expect_gt(sum(colSums(perms != grp) > 0L), 100L)
})

test_that("crossed actors keep their per-actor group counts", {
  ids <- rep(c("a", "b", "c"), each = 6L)
  grp <- c("x", "x", "y", "y", "z", "z",
           "x", "y", "y", "y", "z", "x",
           "z", "z", "z", "x", "y", "y")
  design <- .cs_actor_design(ids, grp, "who")
  expect_identical(design$design, "within")

  set.seed(12)
  perms <- replicate(200L, .cs_actor_permute(design))
  expected <- lapply(split(grp, ids), table)
  kept <- apply(perms, 2L, function(p) {
    identical(lapply(split(p, ids), table), expected)
  })
  expect_true(all(kept))
  expect_gt(sum(colSums(perms != grp) > 0L), 100L)
})

test_that("a mixed design shuffles each part its own way; NA groups stay put", {
  ids <- c("a", "a", "b", "b", "c", "c", "d", "d", "e")
  grp <- c("x", "x", "y", "y", "x", "y", "z", "z", NA)
  design <- .cs_actor_design(ids, grp, "who", warn = FALSE)
  expect_identical(design$design, "mixed")
  expect_identical(design$n_actors, 4L)

  set.seed(13)
  perms <- replicate(100L, .cs_actor_permute(design))
  expect_true(all(is.na(perms[9L, ])))
  # crossed actor c keeps one x and one y
  expect_true(all(apply(perms[5:6, ], 2L, function(v) {
    identical(sort(v), c("x", "y"))
  })))
  # pure actors a, b, d stay whole
  expect_true(all(perms[1L, ] == perms[2L, ]))
  expect_true(all(perms[3L, ] == perms[4L, ]))
  expect_true(all(perms[7L, ] == perms[8L, ]))
})

test_that("with two groups the arrangement count equals .block_design()'s", {
  ids <- rep(sprintf("a%d", 1:7), times = c(2, 3, 2, 4, 1, 3, 2))
  grp <- c("x", "x", "y", "y", "y", "x", "y", "x", "x", "x", "x",
           "y", "y", "x", "y", "x", "y")
  ours <- .cs_actor_design(ids, grp, "who", warn = FALSE)
  theirs <- .block_design(ids, grp == "x", "who", alpha = 0.05, warn = FALSE)
  expect_equal(ours$log_arrangements, theirs$log_arrangements,
               tolerance = 1e-12)
  expect_identical(ours$n_pure, theirs$n_pure)
  expect_identical(ours$n_crossed, theirs$n_crossed)
})

test_that("the arrangement count is a multinomial over pure actors", {
  # 6 pure actors, two per group: 6! / (2! 2! 2!) = 90 arrangements
  design <- .cs_actor_design(sprintf("a%d", 1:6),
                             rep(c("x", "y", "z"), 2L), "who")
  expect_equal(design$log_arrangements, log(90), tolerance = 1e-12)
})

# ---- sequence_compare() with actor ----

test_that("actor on a netobject_group records the design and prints it", {
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Achiever")
  set.seed(1)
  res <- sequence_compare(net, sub = 2L, iter = 20L, actor = "Group")
  expect_identical(res$params$actor, "Group")
  expect_identical(res$params$n_actors, 200L)
  expect_identical(res$params$design, "between")
  expect_output(print(res),
                "permutation: 20 iter \\(fdr\\), actor: Group \\(200 actors, between\\)")

  # the observed table does not depend on the shuffle
  set.seed(1)
  plain <- sequence_compare(net, sub = 2L, iter = 20L)
  keep <- c("pattern", "freq_High", "freq_Low", "resid_High", "resid_Low")
  expect_identical(
    summary(res)[order(summary(res)$pattern), keep],
    summary(plain)[order(summary(plain)$pattern), keep],
    ignore_attr = TRUE)
})

test_that("actor on a data frame is set aside from the positions", {
  d <- .sca_wide()
  set.seed(5)
  # 6 pure actors, two per group: 90 arrangements, enough for p < 0.05
  expect_no_warning(
    res <- sequence_compare(d, group = "grp", actor = "who", sub = 2L,
                            min_freq = 1L, iter = 20L))
  expect_identical(res$params$n_actors, 6L)
  expect_identical(res$params$design, "between")
  expect_false(any(grepl("p1|p2", summary(res)$pattern)))
  # positions are T1..T3 only: 2-grams use states A, B, C
  expect_true(all(grepl("^[ABC]->[ABC]$", summary(res)$pattern)))
})

test_that("few actors warn; enough actors do not", {
  set.seed(21)
  d <- .sca_actor_null(n_actors = 16L)
  expect_no_warning(
    sequence_compare(d, group = "grp", actor = "team", sub = 2L,
                     min_freq = 1L, iter = 10L))
  d4 <- .sca_actor_null(n_actors = 4L)
  expect_warning(
    sequence_compare(d4, group = "grp", actor = "team", sub = 2L,
                     min_freq = 1L, iter = 10L),
    class = "nestimate_few_actors")
})

# ---- error and warning paths ----

test_that("bad actor arguments raise classed conditions", {
  d <- .sca_wide()
  expect_error(
    sequence_compare(d, group = "grp", actor = "nope", sub = 2L,
                     min_freq = 1L, iter = 5L),
    class = "nestimate_actor_missing")
  d_na <- d
  d_na$who[3L] <- NA
  expect_error(
    sequence_compare(d_na, group = "grp", actor = "who", sub = 2L,
                     min_freq = 1L, iter = 5L),
    class = "nestimate_actor_missing")
  expect_error(
    sequence_compare(d, group = "grp", actor = c("who", "who"), sub = 2L,
                     min_freq = 1L, iter = 5L),
    class = "nestimate_bad_actor")

  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Achiever")
  expect_error(sequence_compare(net, sub = 2L, iter = 5L, actor = "Team"),
               class = "nestimate_actor_missing")
})

test_that("actor with a non-permutation test is ignored with a warning", {
  d <- .sca_wide()
  expect_warning(
    res <- sequence_compare(d, group = "grp", actor = "who", sub = 2L,
                            min_freq = 1L, test = "chisq"),
    class = "nestimate_ignored_argument")
  expect_null(res$params$actor)
  plain <- sequence_compare(subset(d, select = -who), group = "grp",
                            sub = 2L, min_freq = 1L, test = "chisq")
  expect_identical(res, plain)
  expect_warning(
    sequence_compare(d, group = "grp", actor = "who", sub = 2L,
                     min_freq = 1L, test = "none"),
    class = "nestimate_ignored_argument")
})

# ---- type I sanity ----

test_that("under an actor-level null, actor shuffling rejects less", {
  skip_on_cran()
  set.seed(2026)
  rates <- vapply(seq_len(40L), function(r) {
    d <- .sca_actor_null()
    by_sequence <- sequence_compare(subset(d, select = -team), group = "grp",
                                    sub = 2L, min_freq = 1L, iter = 100L,
                                    adjust = "none")
    by_actor <- sequence_compare(d, group = "grp", actor = "team", sub = 2L,
                                 min_freq = 1L, iter = 100L, adjust = "none")
    c(sequence = mean(summary(by_sequence)$p_value < 0.05),
      actor = mean(summary(by_actor)$p_value < 0.05))
  }, numeric(2L))
  rate <- rowMeans(rates)
  expect_gt(rate[["sequence"]], 2 * rate[["actor"]])
  expect_gt(rate[["sequence"]], 0.10)
  expect_lt(rate[["actor"]], 0.10)
})
