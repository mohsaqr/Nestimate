# method = "latency": reply network weighted by response speed
# (Templeton et al. 2022). Expected values are computed by hand.

.latency_toy <- function(offset = 0, scale = 1) {
  t0 <- as.POSIXct("2024-01-01 10:00:00", tz = "UTC")
  # session 1: A 0, B 10, A 70, B 100; gap of 2000 s; session 2: C, A, A
  secs <- c(0, 10, 70, 100, 2100, 2400, 2410)
  data.frame(
    chat = "c1",
    sender = c("A", "B", "A", "B", "C", "A", "A"),
    at = t0 + offset + scale * secs,
    stringsAsFactors = FALSE
  )
}

.latency_random <- function(n_chats = 12L, n_msgs = 25L, seed = 7L) {
  .local_seed(seed)
  t0 <- as.POSIXct("2024-03-01", tz = "UTC")
  do.call(rbind, lapply(seq_len(n_chats), function(i) {
    data.frame(
      team = if (i <= n_chats / 2) "red" else "blue",
      chat = paste0("chat", i),
      sender = sample(c("A", "B", "C", "D"), n_msgs, replace = TRUE),
      at = t0 + cumsum(stats::rexp(n_msgs, 1 / 60)),
      stringsAsFactors = FALSE
    )
  }))
}

test_that("latency weights match a hand computation (sum and mean)", {
  d <- .latency_toy()
  net <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                       action = "sender", time = "at",
                       latency_floor = 0.5, latency_cap = "60s")
  # pairs (responder -> previous sender, t): B->A 10, A->B 60, B->A 30,
  # A->C 300 (session 2), A->A 10; C -> B spans the session gap: no pair.
  w <- function(t) 1 - 0.5 * min(t / 60, 1)
  expected <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"),
                                              c("A", "B", "C")))
  expected["B", "A"] <- w(10) + w(30)
  expected["A", "B"] <- w(60)
  expected["A", "C"] <- w(300)
  expected["A", "A"] <- w(10)
  expect_equal(net$weights, expected, tolerance = 1e-12)
  expect_equal(net$weights["C", "B"], 0)
  expect_equal(net$latency_cap_seconds, 60)
  expect_equal(net$params$latency_cap_seconds, 60)

  edges <- extract_edges(net, include_self = TRUE, sort_by = "from")
  expect_named(edges, c("from", "to", "weight", "n_replies",
                        "median_latency"))
  expect_identical(nrow(edges), 4L)
  ba <- subset(edges, from == "B" & to == "A")
  expect_identical(ba$n_replies, 2L)
  expect_equal(ba$median_latency, 20)

  mean_net <- build_network(d, method = "latency", actor = "chat",
                            action = "sender", time = "at",
                            latency_cap = "60s", latency_aggregate = "mean")
  expected_mean <- expected
  expected_mean["B", "A"] <- (w(10) + w(30)) / 2
  expect_equal(mean_net$weights, expected_mean, tolerance = 1e-12)
})

test_that("latency_floor = 1 gives the counts of method = 'reverse'", {
  d <- .latency_random()
  lat <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                       action = "sender", time = "at", latency_floor = 1)
  rev <- build_network(d, method = "reverse", actor = "chat",
                       action = "sender", time = "at")
  expect_equal(lat$weights, rev$weights, tolerance = 1e-12)
  expect_equal(lat$n_replies, rev$weights, tolerance = 1e-12)
})

test_that("floor = 0 gives slow replies no weight but still counts them", {
  d <- .latency_toy()
  net <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                       action = "sender", time = "at",
                       latency_floor = 0, latency_cap = "60s")
  expect_equal(net$weights["A", "C"], 0)  # t = 300 >= cap
  expect_equal(net$weights["A", "B"], 0)  # t = 60 == cap
  expect_equal(net$n_replies["A", "C"], 1)
  expect_equal(net$n_replies["A", "B"], 1)
  expect_equal(net$weights["B", "A"], (1 - 10 / 60) + (1 - 30 / 60),
               tolerance = 1e-12)
})

test_that("a percentile cap is the quantile of sender-change latencies", {
  d <- .latency_toy()
  net <- build_network(d, method = "latency", actor = "chat",
                       action = "sender", time = "at", latency_cap = 0.8)
  # sender-change pairs only: A -> A (10 s) is left out
  expect_equal(net$latency_cap_seconds,
               stats::quantile(c(10, 60, 30, 300), 0.8, names = FALSE))
  default <- build_network(d, method = "latency", actor = "chat",
                           action = "sender", time = "at")
  expect_equal(default$latency_cap_seconds,
               stats::quantile(c(10, 60, 30, 300), 0.95, names = FALSE))
  expect_output(print(default), "Reply Network")
  expect_no_error(summary(default))
})

test_that("with group, every group shares the cap of the full data", {
  d <- .latency_random()
  full <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                        action = "sender", time = "at")
  grouped <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                           action = "sender", time = "at", group = "team")
  expect_s3_class(grouped, "netobject_group")
  caps <- vapply(grouped, function(n) n$latency_cap_seconds, numeric(1))
  expect_equal(unname(caps), rep(full$latency_cap_seconds, 2))
  own <- build_network(subset(d, team == "red"), method = "latency", latency_aggregate = "sum",
                       actor = "chat", action = "sender", time = "at")
  expect_false(isTRUE(all.equal(own$latency_cap_seconds,
                                full$latency_cap_seconds)))
  # with the shared cap, the group networks add up to the full network
  total <- Reduce(`+`, lapply(grouped, function(n) n$weights))
  expect_equal(total, full$weights, tolerance = 1e-12)
})

test_that("no reply crosses a group, also when group is the actor", {
  d <- .latency_random(n_chats = 4L, n_msgs = 10L)
  # one sequence per team is noted (class nestimate_single_sequence)
  grouped <- suppressMessages(
    build_network(d, method = "latency", actor = "team", action = "sender",
                  time = "at", group = "team", time_threshold = FALSE),
    classes = "nestimate_single_sequence")
  # one sequence per team: n - 1 replies each
  n_rep <- vapply(grouped, function(n) sum(n$n_replies), numeric(1))
  expect_equal(unname(n_rep[c("red", "blue")]), c(19, 19))
})

test_that("shifting all times leaves the weights unchanged", {
  base <- build_network(.latency_toy(), method = "latency", actor = "chat",
                        action = "sender", time = "at")
  shifted <- build_network(.latency_toy(offset = 86400 * 3 + 17),
                           method = "latency", actor = "chat",
                           action = "sender", time = "at")
  expect_equal(shifted$weights, base$weights, tolerance = 1e-12)
})

test_that("a percentile cap makes the weights free of the time unit", {
  # no time split, so stretching the times cannot move a session boundary;
  # the chat is then one sequence, which is noted
  build <- function(d) {
    suppressMessages(
      build_network(d, method = "latency", actor = "chat", action = "sender",
                    time = "at", time_threshold = FALSE),
      classes = "nestimate_single_sequence")
  }
  base <- build(.latency_toy())
  stretched <- build(.latency_toy(scale = 60))
  expect_equal(stretched$weights, base$weights, tolerance = 1e-12)
  expect_equal(stretched$latency_cap_seconds, 60 * base$latency_cap_seconds,
               tolerance = 1e-12)
})

test_that("self-continuations stay on the diagonal", {
  net <- build_network(.latency_toy(), method = "latency", actor = "chat",
                       action = "sender", time = "at", latency_cap = "60s")
  expect_equal(diag(net$weights), c(A = 1 - 0.5 * 10 / 60, B = 0, C = 0),
               tolerance = 1e-12)
  scaled <- build_network(.latency_toy(), method = "latency", actor = "chat",
                          action = "sender", time = "at", latency_cap = "60s",
                          scaling = "max")
  expect_gt(scaled$weights["A", "A"], 0)
})

test_that("resampling verbs other than bootstrap_network() refuse a latency network", {
  # bootstrap_network() resamples latency networks by their sequences; every
  # other verb that shuffles or drops sequences refuses one, alone or inside a
  # group
  d <- .latency_random()
  net <- build_network(d, method = "latency", actor = "chat",
                       action = "sender", time = "at")
  grouped <- build_network(d, method = "latency", actor = "chat",
                           action = "sender", time = "at", group = "team")
  refused <- "nestimate_latency_resampling"
  expect_s3_class(bootstrap_network(net, iter = 10L, seed = 1), "net_bootstrap")
  expect_s3_class(bootstrap_network(grouped, iter = 10L, seed = 1),
                  "net_bootstrap_group")
  expect_error(permutation(grouped$red, grouped$blue, iter = 10L, seed = 1),
               class = refused)
  expect_error(permutation(grouped, iter = 10L, seed = 1), class = refused)
  expect_error(permutation_diagnostics(grouped, actor = "chat", iter = 10L,
                                       seed = 1), class = refused)
  expect_error(centrality_stability(net, iter = 10L, seed = 1),
               class = refused)
  expect_error(centrality_stability(grouped, iter = 10L, seed = 1),
               class = refused)
  expect_error(casedrop_reliability(net, iter = 10L, seed = 1),
               class = refused)
  expect_error(network_reliability(net, iter = 10L, seed = 1),
               class = refused)
  expect_error(network_reliability(grouped, iter = 10L, seed = 1),
               class = refused)
  expect_error(compare_networks(grouped$red, grouped$blue,
                                test = "permutation", iter = 10L),
               class = refused)
  expect_error(compare_networks(grouped$red, grouped$blue,
                                test = "bootstrap", iter = 10L),
               class = refused)
  # the message names the verb that does resample them
  expect_error(permutation(grouped, iter = 10L, seed = 1),
               "bootstrap_network\\(\\) resamples latency networks")
  # without resampling, the networks are still compared
  expect_s3_class(compare_networks(grouped$red, grouped$blue),
                  "net_network_comparison")
})

test_that("bad input raises nestimate_bad_input", {
  d <- .latency_toy()
  lat <- function(...) {
    build_network(d, method = "latency", actor = "chat", action = "sender",
                  ...)
  }
  # no time
  expect_error(lat(), class = "nestimate_bad_input")
  # not an event log
  wide <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  expect_error(build_network(wide, method = "latency"),
               class = "nestimate_bad_input")
  # floor
  expect_error(lat(time = "at", latency_floor = 1.5),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", latency_floor = -0.1),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", latency_floor = "half"),
               class = "nestimate_bad_input")
  # cap
  expect_error(lat(time = "at", latency_cap = 0),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", latency_cap = 120),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", latency_cap = "two minutes"),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", latency_cap = "0s"),
               class = "nestimate_bad_input")
  # aggregate
  expect_error(lat(time = "at", latency_aggregate = "median"),
               class = "nestimate_bad_input")
  # weighted means inverse sequence length
  expect_error(lat(time = "at", weighted = TRUE),
               class = "nestimate_bad_input")
  expect_error(lat(time = "at", params = list(weighted = TRUE)),
               class = "nestimate_bad_input")
  # start / end markers have no times
  expect_error(lat(time = "at", start = TRUE), class = "nestimate_bad_input")
  # latency settings on another method
  expect_error(build_network(d, method = "reverse", actor = "chat",
                             action = "sender", time = "at",
                             latency_floor = 0.2),
               class = "nestimate_bad_input")
})

test_that("latency is a registered directed estimator", {
  est <- list_estimators()
  expect_true("latency" %in% est$name)
  expect_true(get_estimator("latency")$directed)
})

# ---- latency_scope ----

# W: session 1 alternates A/B every 5 s, then a 2000 s gap, session 2
# alternates C/D every 60 s. V: A, B, A 30 s apart (2 replies: too few for
# a cap of its own).
.latency_scoped <- function(prefix = "") {
  t0 <- as.POSIXct("2024-05-01", tz = "UTC")
  w <- data.frame(chat = paste0(prefix, "W"),
                  sender = c("A", "B", "A", "B", "A", "B",
                             "C", "D", "C", "D", "C", "D"),
                  at = t0 + c(0, 5, 10, 15, 20, 25,
                              2025, 2085, 2145, 2205, 2265, 2325))
  v <- data.frame(chat = paste0(prefix, "V"), sender = c("A", "B", "A"),
                  at = t0 + c(0, 30, 60))
  rbind(w, v)
}

test_that("global scope is the default and unchanged", {
  d <- .latency_scoped()
  default <- build_network(d, method = "latency", actor = "chat",
                           action = "sender", time = "at")
  global <- build_network(d, method = "latency", actor = "chat",
                          action = "sender", time = "at",
                          latency_scope = "global")
  expect_equal(default$weights, global$weights)
  # sender-change latencies: 5 x 5 s, 5 x 60 s, 2 x 30 s
  expect_equal(default$latency_cap_seconds,
               stats::quantile(c(rep(5, 5), rep(60, 5), 30, 30), 0.95,
                               names = FALSE))
  expect_null(default$latency_unit_caps)
})

test_that("per-actor and per-session caps match a hand computation", {
  d <- .latency_scoped()
  global_cap <- stats::quantile(c(rep(5, 5), rep(60, 5), 30, 30), 0.95,
                                names = FALSE)
  w_g <- function(t) 1 - 0.5 * min(t / global_cap, 1)

  expect_message(
    actor <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                           action = "sender", time = "at",
                           latency_scope = "actor"),
    class = "nestimate_latency_fallback")
  w_cap <- stats::quantile(c(rep(5, 5), rep(60, 5)), 0.95, names = FALSE)
  expect_equal(actor$latency_unit_caps, c(W = w_cap))
  # V falls back to the global cap
  expect_equal(actor$weights["B", "A"],
               3 * (1 - 0.5 * 5 / w_cap) + w_g(30), tolerance = 1e-12)

  expect_message(
    sess <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                          action = "sender", time = "at",
                          latency_scope = "session"),
    class = "nestimate_latency_fallback")
  expect_equal(sort(unname(sess$latency_unit_caps)), c(5, 60))
  # session 1 replies sit at their cap (5 s): weight 0.5 each
  expect_equal(sess$weights["B", "A"], 3 * 0.5 + w_g(30), tolerance = 1e-12)
  expect_equal(sess$weights["A", "B"], 2 * 0.5 + w_g(30), tolerance = 1e-12)
  expect_equal(sess$weights["D", "C"], 3 * 0.5, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(sess$weights, actor$weights)))
  # the global cap is still the one of all the data
  expect_equal(sess$latency_cap_seconds, global_cap)
})

test_that("responder scope gives each responding sender its own cap", {
  # B always answers in 5 s, A in 50 s
  t0 <- as.POSIXct("2024-05-01", tz = "UTC")
  gaps <- rep(c(5, 50), length.out = 11)
  d <- data.frame(chat = "c", sender = rep(c("A", "B"), 6),
                  at = t0 + cumsum(c(0, gaps)))
  resp <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                        action = "sender", time = "at",
                        latency_scope = "responder")
  expect_equal(resp$latency_unit_caps[c("A", "B")], c(A = 50, B = 5))
  expect_equal(resp$weights["B", "A"], 6 * 0.5, tolerance = 1e-12)
  expect_equal(resp$weights["A", "B"], 5 * 0.5, tolerance = 1e-12)
  global <- build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat",
                          action = "sender", time = "at")
  expect_equal(global$weights["B", "A"],
               6 * (1 - 0.5 * 5 / global$latency_cap_seconds),
               tolerance = 1e-12)
})

test_that("unit caps are resolved before a group split", {
  d <- rbind(transform(.latency_scoped("r"), team = "red"),
             transform(.latency_scoped("b"), team = "blue"))
  full <- suppressMessages(
    build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat", action = "sender",
                  time = "at", latency_scope = "actor"),
    classes = "nestimate_latency_fallback")
  grouped <- suppressMessages(
    build_network(d, method = "latency", latency_aggregate = "sum", actor = "chat", action = "sender",
                  time = "at", latency_scope = "actor", group = "team"),
    classes = "nestimate_latency_fallback")
  invisible(lapply(grouped, function(n) {
    expect_equal(n$latency_unit_caps, full$latency_unit_caps)
    expect_equal(n$latency_cap_seconds, full$latency_cap_seconds)
  }))
  total <- Reduce(`+`, lapply(grouped, function(n) n$weights))
  expect_equal(total, full$weights, tolerance = 1e-12)
})

test_that("the stored params and event times reproduce a scoped network", {
  d <- do.call(rbind, lapply(1:4, function(i) .latency_scoped(paste0("k", i))))
  net <- suppressMessages(
    build_network(d, method = "latency", actor = "chat", action = "sender",
                  time = "at", latency_scope = "session"),
    classes = "nestimate_latency_fallback")
  # the states with their times and units, re-estimated with the stored
  # params (which hold the resolved caps), give the network back
  rebuilt <- do.call(get_estimator("latency")$fn,
                     c(list(data = cbind(net$data, net$latency_time)),
                       net$params))
  expect_equal(rebuilt$matrix, net$weights, tolerance = 1e-12)
  expect_equal(net$params$latency_unit_caps, net$latency_unit_caps)
})

test_that("scope with a cap in seconds is ignored, a bad scope refused", {
  d <- .latency_scoped()
  expect_message(
    secs <- build_network(d, method = "latency", actor = "chat",
                          action = "sender", time = "at",
                          latency_cap = "60s", latency_scope = "actor"),
    class = "nestimate_ignored_argument")
  plain <- build_network(d, method = "latency", actor = "chat",
                         action = "sender", time = "at", latency_cap = "60s")
  expect_equal(secs$weights, plain$weights)
  expect_error(build_network(d, method = "latency", actor = "chat",
                             action = "sender", time = "at",
                             latency_scope = "team"),
               class = "nestimate_bad_input")
  expect_error(build_network(d, method = "latency", action = "sender",
                             time = "at", latency_scope = "actor"),
               class = "nestimate_bad_input")
})

# Latency settings given in `params` are used when the argument is not
# given, as the alias's arguments are; an explicit argument wins.
test_that("latency settings in params are not overwritten by the defaults", {
  replies <- data.frame(id = rep(1:6, each = 2), state = rep(c("A", "B"), 6),
                        t = rep(c(0, 1), 6))
  by_alias <- suppressMessages(build_latency(replies, action = "state",
                                             time = "t", actor = "id",
                                             cap = "2s", floor = 0.25,
                                             aggregate = "sum"))
  by_params <- suppressMessages(build_latency(replies, action = "state",
    time = "t", actor = "id",
    params = list(latency_cap = "2s", latency_floor = 0.25,
                  latency_aggregate = "sum")))
  expect_equal(by_params$weights, by_alias$weights)
  expect_equal(by_params$latency_cap_seconds, 2)
  explicit <- suppressMessages(build_latency(replies, action = "state",
    time = "t", actor = "id", cap = "4s", params = list(latency_cap = "2s")))
  expect_equal(explicit$latency_cap_seconds, 4)
})
