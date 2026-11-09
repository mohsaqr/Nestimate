# Bad input to the build/conversion verbs raises class nestimate_bad_input.
testthat::skip_on_cran()

seqs <- data.frame(V1 = c("A", "B", "A"), V2 = c("B", "A", "B"),
                   V3 = c("A", "B", "B"), stringsAsFactors = FALSE)

test_that("build_network() input checks are classed", {
  expect_error(build_network(seqs, method = c("relative", "frequency")),
               class = "nestimate_bad_input")
  expect_error(build_network(seqs, method = "relative", threshold = -1),
               class = "nestimate_bad_input")
  expect_error(build_network(seqs, method = "relative", start = 3),
               class = "nestimate_bad_input")
  expect_error(build_network(seqs, method = "relative", scaling = "bogus"),
               class = "nestimate_bad_input")
  expect_error(build_network(seqs, method = "relative", state_cols = "nope"),
               class = "nestimate_bad_input")
  expect_error(build_network(srl_strategies, method = "cor", level = "between"),
               class = "nestimate_bad_input")
})

test_that("predictability() refuses a transition network by class", {
  net <- build_network(seqs, method = "relative")
  expect_error(predictability(net), class = "nestimate_bad_input")
})

test_that("format converters classify bad input", {
  expect_error(wide_to_long(seqs, time_prefix = "Z"),
               class = "nestimate_bad_input")
  expect_error(wide_to_long("not a frame"), class = "nestimate_bad_input")
  expect_error(long_to_wide(group_regulation_long, id_col = "nope"),
               class = "nestimate_bad_input")
  expect_error(prepare_for_tna(data.frame(x = 1:3), type = "sequences"),
               class = "nestimate_bad_input")
  expect_error(action_to_onehot(group_regulation_long, action_col = "nope"),
               class = "nestimate_bad_input")
})

test_that("prepare() input checks are classed", {
  expect_error(prepare(group_regulation_long, action = "nope"),
               class = "nestimate_bad_input")
  expect_error(prepare(group_regulation_long, action = "Action",
                       timezone = "Not/AZone"),
               class = "nestimate_bad_input")
})

test_that("estimator registry checks are classed", {
  expect_error(get_estimator("no_such_estimator"),
               class = "nestimate_bad_input")
  expect_error(remove_estimator("no_such_estimator"),
               class = "nestimate_bad_input")
  expect_error(register_estimator(1, identity, "d", TRUE),
               class = "nestimate_bad_input")
})

test_that("as_netobject() and validate_netobject() classify bad input", {
  expect_error(as_netobject(1:3), class = "nestimate_bad_input")
  expect_error(validate_netobject(list()), class = "nestimate_bad_input")
})

test_that("frequencies() and convert_sequence_format() classify bad input", {
  expect_error(frequencies("not a frame"), class = "nestimate_bad_input")
  expect_error(convert_sequence_format("not a frame"),
               class = "nestimate_bad_input")
})

test_that("wtna() input checks are classed", {
  expect_error(wtna(learning_activities, codes = c("Reading", "Video"),
                    window_size = 0),
               class = "nestimate_bad_input")
  expect_error(wtna(learning_activities, codes = c("Reading", "nope")),
               class = "nestimate_bad_input")
})

test_that("cooccurrence() input checks are classed", {
  expect_error(cooccurrence(list(c("a", "b")), threshold = -1),
               class = "nestimate_bad_input")
  expect_error(cooccurrence(list(character(0))),
               class = "nestimate_bad_input")
})
