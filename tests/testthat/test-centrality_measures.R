# Randomized shortest-path betweenness (Kivimaki et al. 2016) on tiny
# hand-checked networks.

.rsp_w2 <- function() {
  matrix(c(0, 3, 1, 0,  2, 0, 4, 1,  1, 1, 0, 2,  0, 0, 1, 0), 4,
         byrow = TRUE, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
}

test_that("BetweennessRSP on a strongly connected network", {
  expect_equal(.rsp_betweenness(.rsp_w2()), c(A = 1, B = 4, C = 12, D = 4))
})

test_that("BetweennessRSP counts the absorption correction only for sources that reach the target", {
  # E is never entered: walks from E exist, walks to E do not. The
  # correction for walks absorbed at each target is counted once per source
  # that reaches it (4 for E, 5 for the others), not once per state.
  W3 <- rbind(cbind(.rsp_w2(), E = 0), E = c(1, 1, 0, 0, 0))
  expect_equal(.rsp_betweenness(W3), c(A = 7, B = 11, C = 17, D = 8, E = 1))
  # a never-entered state lies on no path between others: the minimum
  expect_equal(unname(.rsp_betweenness(W3)["E"]), 1)
})

test_that("BetweennessRSP is NA with a sink state and unchanged by relabelling", {
  W <- .rsp_w2()
  W["D", ] <- 0
  expect_true(all(is.na(.rsp_betweenness(W))))
  perm <- c(3, 1, 4, 2)
  Wp <- .rsp_w2()[perm, perm]
  expect_equal(.rsp_betweenness(Wp), .rsp_betweenness(.rsp_w2())[perm])
})
