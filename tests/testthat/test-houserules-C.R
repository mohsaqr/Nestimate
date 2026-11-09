# House-rule fixes in the group C files (build_network, estimators, registry,
# converters, extraction, frequencies, wtna): condition classes, contract
# messages in place of deparsed calls, corrected printed text, new warnings.

msg_of <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

test_that("build_network() says what to pass for prepare() output", {
  p <- prepare(group_regulation_long[seq_len(500), ], actor = "Actor",
               action = "Action", time = "Time")
  expect_error(build_network(p, method = "relative"),
               class = "nestimate_bad_input")
  expect_match(msg_of(build_network(p, method = "relative")),
               "result of prepare()", fixed = TRUE)
})

test_that("build_network() argument checks state the contract", {
  expect_error(build_network(trajectories, method = "relative",
                             threshold = "LW"),
               class = "nestimate_bad_input")
  expect_match(msg_of(build_network(trajectories, method = "relative",
                                    threshold = "LW")),
               "`threshold` must be a single number of at least 0",
               fixed = TRUE)
  expect_match(msg_of(build_network(trajectories, method = 1)),
               "`method` must be a single character string", fixed = TRUE)
  expect_match(msg_of(build_network(trajectories, method = "relative",
                                    params = "a")),
               "`params` must be a list", fixed = TRUE)
})

test_that("a matrix given to the attention method is refused in words", {
  m <- matrix(c(1, 2, 1, 3, 1, 1, 1, 1, 1), 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(build_atna(m), class = "nestimate_bad_input")
  expect_match(msg_of(build_atna(m)), "got <matrix>", fixed = TRUE)
})

test_that("a decay function returning one value is a classed error", {
  wide <- data.frame(V1 = c("A", "B", "C"), V2 = c("B", "C", "A"),
                     V3 = c("A", "A", "B"))
  expect_error(
    build_atna(wide, params = list(decay = function(ti, tj, lam) 1)),
    class = "nestimate_bad_input")
  # a decay returning one weight per pair still gives the default network
  same <- build_atna(wide, params = list(
    decay = function(ti, tj, lam) exp(-abs(ti - tj) / lam)))
  expect_equal(same$weights, build_atna(wide)$weights)
})

test_that("print.netobject shows no level label for an mgm network", {
  net <- suppressWarnings(build_network(srl_strategies, method = "mgm"))
  header <- capture.output(print(net))[1L]
  expect_false(grepl("-person]", header, fixed = TRUE))
  expect_match(header, "method: mgm", fixed = TRUE)
})

test_that("predictability() warns when a neighbour block is singular", {
  a <- as.numeric(scale(srl_strategies[[1L]]))
  d <- data.frame(A = a, B = a, C = a + srl_strategies[[2L]] / 10)
  net <- build_network(d, method = "cor", predictability = FALSE)
  expect_warning(pred <- predictability(net),
                 class = "nestimate_singular_neighbours")
  expect_identical(pred$R2[pred$node == "C"], 0)
  # a regular network stays silent
  expect_no_warning(predictability(
    build_network(srl_strategies, method = "cor", predictability = FALSE)))
})

test_that("grouped predictability() warns about an argument it drops", {
  d <- srl_strategies
  d$half <- rep(c("a", "b"), length.out = nrow(d))
  grp <- build_network(d, method = "cor", group = "half")
  expect_warning(res <- predictability(grp, data = d),
                 class = "nestimate_ignored_argument")
  expect_identical(res, predictability(grp))
})

test_that("a whole-frame refusal names the column beside the run", {
  d <- as.data.frame(trajectories)
  d$grade <- seq_len(nrow(d))
  expect_error(frequencies(d), class = "nestimate_bad_input")
  expect_match(msg_of(frequencies(d)), "beside the numbered run", fixed = TRUE)
  expect_match(msg_of(frequencies(d)), "`grade`", fixed = TRUE)
})

test_that("frequencies(format = 'long') without action names the argument", {
  expect_error(frequencies(group_regulation_long, format = "long"),
               class = "nestimate_bad_input")
  expect_match(msg_of(frequencies(group_regulation_long, format = "long")),
               "`action`", fixed = TRUE)
})

test_that("extraction results print without a class attribute block", {
  net <- build_network(trajectories, method = "relative")
  tm <- extract_transition_matrix(net)
  out <- capture.output(res <- print.nest_transition_matrix(tm))
  expect_false(any(grepl("attr(,", out, fixed = TRUE)))
  expect_identical(out, capture.output(print(unclass(tm))))
  expect_identical(res, tm)
  ip <- extract_initial_probs(net)
  out <- capture.output(res <- print.nest_initial_probs(ip))
  expect_false(any(grepl("attr(,", out, fixed = TRUE)))
  expect_identical(res, ip)
})

test_that("the uniform initial-probability fallback is a classed warning", {
  net <- build_network(trajectories, method = "co_occurrence")
  expect_warning(ip <- extract_initial_probs(net),
                 class = "nestimate_uniform_initial")
  expect_equal(unname(unclass(ip)), rep(1 / 3, 3))
  expect_error(extract_initial_probs(1:3), class = "nestimate_bad_input")
})

test_that("wtna() argument errors state the contract", {
  la <- learning_activities
  expect_error(wtna(la, actor = "nope"), class = "nestimate_bad_input")
  expect_match(msg_of(wtna(la, actor = "nope")), "`actor` must name columns",
               fixed = TRUE)
  expect_match(msg_of(wtna(la, actor = "student", codes = "Reading")),
               "At least two one-hot code columns", fixed = TRUE)
  expect_match(msg_of(wtna(la, actor = "student", window_size = 0)),
               "`window_size` must be a single whole number", fixed = TRUE)
  expect_match(msg_of(wtna(la, actor = "student", window_size = 2.7)),
               "`window_size` must be a single whole number", fixed = TRUE)
})

test_that("registry errors state the contract and aliases are flagged", {
  expect_error(register_estimator("x", 1, "d", TRUE),
               class = "nestimate_bad_input")
  expect_match(msg_of(register_estimator("x", 1, "d", TRUE)),
               "`fn` must be a function", fixed = TRUE)
  expect_match(msg_of(register_estimator("", identity, "d", TRUE)),
               "`name` must be a single non-empty", fixed = TRUE)
  expect_match(msg_of(get_estimator(1)), "`name` must be a single character",
               fixed = TRUE)
  fn <- function(data, ...) NULL
  expect_warning(register_estimator("tna", fn, "shadowed", TRUE),
                 class = "nestimate_alias_shadowed")
  remove_estimator("tna")
  expect_no_warning(register_estimator("hr_c_custom", fn, "fine", TRUE))
  remove_estimator("hr_c_custom")
})

test_that("converter conditions are classed and worded", {
  long <- group_regulation_long[seq_len(200), ]
  expect_error(action_to_onehot(long, action_col = "nope"),
               class = "nestimate_bad_input")
  expect_match(msg_of(action_to_onehot(long, action_col = "nope")),
               "`action_col` must name one column", fixed = TRUE)
  expect_match(msg_of(wide_to_long(1)), "`data` must be a data frame",
               fixed = TRUE)
  wide <- as.data.frame(trajectories)
  names(wide) <- paste0("V", seq_len(ncol(wide)))
  expect_warning(prepare_for_tna(wide, state_names = c("Active", "Average")),
                 class = "nestimate_unknown_states")
})

test_that("prepare() column errors state the contract", {
  long <- group_regulation_long[seq_len(200), ]
  expect_error(prepare(long, action = "nope"), class = "nestimate_bad_input")
  expect_match(msg_of(prepare(long, action = "nope")),
               "`action` must name one column", fixed = TRUE)
  expect_match(msg_of(prepare(long, action = "Action", actor = "nope")),
               "`actor` must name columns", fixed = TRUE)
  expect_match(msg_of(prepare(long, action = "Action", time = "nope")),
               "`time` must name one column", fixed = TRUE)
})

test_that("a missing declared column is reported against the data", {
  d <- as.data.frame(trajectories)
  expect_error(extract_sequences(d, cols = "zz"),
               class = "nestimate_bad_input")
  expect_match(msg_of(extract_sequences(d, cols = "zz")),
               "not found in the data", fixed = TRUE)
})

test_that("vectorised converters return what the loops returned", {
  # hand-checked on a ragged two-sequence log
  log <- data.frame(Actor = c("a", "a", "a", "b"), Time = c(1, 2, 3, 1),
                    Action = factor(c("X", "Y", "X", "Y")))
  wide <- long_to_wide(log)
  expect_identical(names(wide), c("Actor", "V1", "V2", "V3"))
  expect_identical(as.character(unlist(wide[wide$Actor == "b", -1L])),
                   c("Y", NA, NA))
  tna_ready <- prepare_for_tna(data.frame(V1 = factor(c("A", "B")),
                                          V2 = c("B", "A")))
  expect_identical(vapply(tna_ready, class, character(1L)),
                   c(V1 = "character", V2 = "character"))
})
