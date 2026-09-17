# session_ids(): the session behind each sequence of a network or a fit.

# Long events in which every session's first action encodes the session, so
# the right id for each sequence row is known from the sequence alone.
.session_events <- function(seed = 1) {
  set.seed(seed)
  keys <- expand.grid(grp = c("a", "b"), user = c("u1", "u2", "u3"),
                      step = c("x | 1", "y", "z"), stringsAsFactors = FALSE)
  keys$first <- sprintf("F%02d", seq_len(nrow(keys)))
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    data.frame(keys[i, c("grp", "user", "step")], t = 1:5,
               action = c(keys$first[i], sample(c("A", "B", "C"), 4, TRUE)),
               row.names = NULL, stringsAsFactors = FALSE)
  })
  ev <- do.call(rbind, rows)
  list(events = ev[sample(nrow(ev)), ], keys = keys)
}

test_that("each sequence gets the actor and session it was built from", {
  s <- .session_events()
  net <- build_network(s$events, actor = c("grp", "user"), session = "step",
                       action = "action", order = "t", method = "relative")
  ids <- session_ids(net)
  expect_s3_class(ids, "data.frame")
  expect_named(ids, c("sequence", "grp", "user", "step", "session_label"))
  expect_equal(nrow(ids), nrow(net$data))
  expect_equal(ids$sequence, seq_len(nrow(net$data)))
  # the first action of each sequence identifies its true key
  truth <- s$keys[match(net$data$T1, s$keys$first), c("grp", "user", "step")]
  expect_identical(ids$grp, truth$grp)
  expect_identical(ids$user, truth$user)
  # a session id containing the label separator is returned intact
  expect_identical(ids$step, truth$step)
  expect_true("x | 1" %in% ids$step)
  expect_null(dim(ids$grp))
})

test_that("the ids do not depend on the row order of the input", {
  s <- .session_events()
  build <- function(ev) {
    build_network(ev, actor = c("grp", "user"), session = "step",
                  action = "action", order = "t", method = "relative")
  }
  a <- build(s$events)
  b <- build(s$events[rev(seq_len(nrow(s$events))), ])
  expect_identical(session_ids(a), session_ids(b))
  expect_identical(a$data, b$data)
})

test_that("a mixture fit reports cluster and posterior per session", {
  s <- .session_events()
  net <- build_network(s$events, actor = c("grp", "user"), session = "step",
                       action = "action", order = "t", method = "relative")
  fit <- build_mmm(net, k = 2, n_starts = 2, max_iter = 50, seed = 1)
  ids <- session_ids(fit)
  expect_named(ids, c("sequence", "grp", "user", "step", "session_label",
                      "cluster", "posterior"))
  expect_identical(ids$cluster, as.integer(fit$assignments))
  expect_equal(ids$posterior, apply(fit$posterior, 1, max))
  expect_identical(ids[c("grp", "user", "step")], session_ids(net)[c("grp", "user", "step")])
  expect_true(all(ids$posterior >= 1 / 2 & ids$posterior <= 1))
})

test_that("a mixture fit that drops sequences keeps ids aligned", {
  s <- .session_events()
  s$events$score <- ifelse(s$events$user == "u2", NA_real_, 1)
  s$events$score[s$events$user == "u3"] <- 2
  net <- build_network(s$events, actor = c("grp", "user"), session = "step",
                       action = "action", order = "t", method = "relative")
  fit <- suppressWarnings(build_mmm(net, k = 2, n_starts = 2, max_iter = 50,
                                    seed = 1, covariates = "score"))
  ids <- session_ids(fit)
  expect_equal(nrow(ids), length(fit$assignments))
  expect_false("u2" %in% ids$user)
  truth <- s$keys[match(fit$data$T1, s$keys$first), "user"]
  expect_identical(ids$user, truth)
})

test_that("a sequence clustering reports the cluster per session", {
  s <- .session_events()
  net <- build_network(s$events, actor = c("grp", "user"), session = "step",
                       action = "action", order = "t", method = "relative")
  cl <- build_clusters(net, k = 2)
  ids <- session_ids(cl)
  expect_identical(ids$cluster, unname(as.integer(cl$assignments)))
  expect_identical(ids$step, session_ids(net)$step)
})

test_that("a single actor without a session column gives the actor only", {
  ev <- data.frame(student = rep(c("s2", "s1"), each = 4),
                   action = c("A", "B", "A", "C", "B", "C", "A", "B"))
  net <- build_network(ev, actor = "student", action = "action",
                       method = "relative")
  ids <- session_ids(net)
  expect_named(ids, c("sequence", "student", "session_label"))
  expect_setequal(ids$student, c("s1", "s2"))
})

test_that("objects without session metadata raise a classed error", {
  wide <- data.frame(T1 = c("A", "B", "C"), T2 = c("B", "C", "A"))
  net <- build_network(wide, method = "relative")
  expect_error(session_ids(net), class = "nestimate_no_session_ids")
  expect_error(session_ids(list()), class = "nestimate_no_session_ids")
  fit <- build_mmm(wide, k = 2, n_starts = 1, max_iter = 5, seed = 1)
  expect_error(session_ids(fit), class = "nestimate_no_session_ids")
})

test_that("metadata of the wrong length raises a classed error", {
  s <- .session_events()
  net <- build_network(s$events, actor = c("grp", "user"), session = "step",
                       action = "action", order = "t", method = "relative")
  net$metadata <- net$metadata[-1, ]
  expect_error(session_ids(net), class = "nestimate_session_ids_misaligned")
})

test_that("prepare() returns metadata in sequence order with the session column", {
  s <- .session_events()
  prep <- prepare(s$events, actor = c("grp", "user"), session = "step",
                  action = "action", order = "t")
  expect_true("step" %in% names(prep$meta_data))
  expect_equal(nrow(prep$meta_data), nrow(prep$sequence_data))
  truth <- s$keys[match(prep$sequence_data$T1, s$keys$first), "step"]
  expect_identical(prep$meta_data$step, truth)
})
