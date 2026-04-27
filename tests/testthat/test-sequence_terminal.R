# Tests for actor_endpoints() and mark_terminal_state().

test_that("actor_endpoints returns one row per actor with the expected columns", {
  ep <- actor_endpoints(trajectories)
  expect_s3_class(ep, "data.frame")
  expect_equal(nrow(ep), nrow(trajectories))
  expect_setequal(names(ep),
                  c("actor", "first_state", "last_state",
                    "first_step", "last_step", "n_observed", "dropped_out"))
  expect_true(all(ep$first_step <= ep$last_step | is.na(ep$first_step)))
  expect_true(all(ep$n_observed >= 0L))
  expect_true(any(ep$dropped_out))
})

test_that("actor_endpoints flags terminal-NA rows but not all-NA rows", {
  M <- rbind(
    c("A", "B", "A", NA, NA),  # dropout at step 4
    c("A", "B", "A", "B", "A"), # completer
    c(NA, NA, NA, NA, NA),      # never observed
    c("B", NA, "A", NA, NA)     # NOT terminal-NA pattern (NA in middle)
  )
  ep <- actor_endpoints(as.data.frame(M, stringsAsFactors = FALSE))
  expect_equal(ep$dropped_out, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(ep$n_observed, c(3L, 5L, 0L, 2L))
  expect_equal(ep$last_state, c("A", "A", NA, "A"))
})

test_that("mark_terminal_state fills terminal NAs with the chosen label", {
  M <- mark_terminal_state(trajectories, state = "Dropout")
  expect_s3_class(M, "data.frame")
  expect_equal(dim(as.matrix(M)), dim(trajectories))
  expect_true(any(unlist(lapply(M, function(c) any(c == "Dropout")))))
  expect_equal(attr(M, "terminal_state"), "Dropout")
  # No actor-row should have NA after the first Dropout
  M_mat <- as.matrix(M)
  for (i in seq_len(nrow(M_mat))) {
    row <- M_mat[i, ]
    drop_pos <- which(row == "Dropout")
    if (length(drop_pos) > 0L) {
      after <- row[min(drop_pos):length(row)]
      expect_true(all(after == "Dropout"),
                  info = sprintf("row %d has NA after Dropout", i))
    }
  }
})

test_that("mark_terminal_state -> build_network -> chain_structure flags the absorbing state", {
  M <- mark_terminal_state(trajectories, state = "Dropout")
  net <- build_network(M, method = "relative")
  cs <- chain_structure(net)
  expect_true("Dropout" %in% cs$absorbing_states)
  expect_false(cs$is_regular)
  expect_false(cs$is_irreducible)
  expect_true(is.na(cs$is_reversible))
  # Absorption probabilities are 1 for every transient state (only one absorbing target)
  ap <- cs$absorption_probabilities
  expect_equal(unname(ap[, "Dropout"]), rep(1, nrow(ap)))
  # Mean absorption time is positive and finite for transient states
  expect_true(all(cs$mean_absorption_time > 0 &
                  is.finite(cs$mean_absorption_time)))
})

test_that("mark_first_state fills leading NAs with the chosen label", {
  M <- rbind(
    c(NA,  NA, "A", "B", "C"),  # late entry at step 3
    c("A","B", "A", "B", "C"),  # observed throughout
    c(NA, "A", "B", "C", "A")   # late entry at step 2
  )
  out <- mark_first_state(as.data.frame(M, stringsAsFactors = FALSE),
                          state = "Start")
  expect_s3_class(out, "data.frame")
  expect_equal(unname(unlist(out[1, ])),
               c("Start", "Start", "A", "B", "C"))
  expect_equal(unname(unlist(out[2, ])),
               c("A", "B", "A", "B", "C"))
  expect_equal(unname(unlist(out[3, ])),
               c("Start", "A", "B", "C", "A"))
  expect_equal(attr(out, "leading_state"), "Start")
})

test_that("mark_first_state -> build_network: Start is transient (not absorbing)", {
  M <- rbind(
    c(NA, "A", "B", "A", "B"),
    c(NA, NA,  "A", "A", "B"),
    c("A","B", "A", "B", "A")
  )
  filled <- mark_first_state(as.data.frame(M, stringsAsFactors = FALSE),
                             state = "Start")
  net <- build_network(filled, method = "relative")
  cs <- chain_structure(net)
  # Start has at least one transition out (to a real state), so it is
  # NOT a fixed point — should not appear as absorbing.
  expect_false("Start" %in% cs$absorbing_states)
})

test_that("mark_terminal_state warns and renames if state already present", {
  M <- as.data.frame(trajectories, stringsAsFactors = FALSE)
  expect_warning(out <- mark_terminal_state(M, state = "Active"),
                 "already appears")
  expect_true(grepl("^Active_\\d+$", attr(out, "terminal_state")))
})
