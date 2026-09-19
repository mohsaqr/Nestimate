# build_mcml(<mcml>, combine =, expand =, clusters =): re-partition and
# re-estimate an existing mcml from the sequences it carries.

make_repart_mcml <- function() {
  build_mcml(group_regulation_long,
             clusters = list(P = c("plan", "consensus"),
                             R = c("coregulate", "discuss", "synthesis"),
                             S = c("cohesion", "emotion", "monitor", "adapt")),
             actor = "Actor", action = "Action", time = "Time")
}
repart_seqs <- function(mc) {
  Reduce(function(a, b) { a[is.na(a)] <- b[is.na(a)]; a },
         lapply(mc$clusters, function(z) z$data))
}

test_that("an mcml without changes comes back unchanged or identical", {
  mc <- make_repart_mcml()
  expect_identical(build_mcml(mc), mc)
  # same partition, re-estimated from the carried sequences: identical model
  expect_equal(unclass(build_mcml(mc, clusters = mc$cluster_members)), unclass(mc))
})

test_that("combine equals a fresh build with the merged partition", {
  mc <- make_repart_mcml()
  merged <- build_mcml(mc, combine = c("P", "R"))
  fresh  <- build_mcml(group_regulation_long,
                       clusters = list(`P + R` = c("plan", "consensus", "coregulate",
                                                   "discuss", "synthesis"),
                                       S = c("cohesion", "emotion", "monitor", "adapt")),
                       actor = "Actor", action = "Action", time = "Time")
  expect_equal(merged$macro$weights, fresh$macro$weights)
  expect_equal(merged$cluster_members, fresh$cluster_members)
  expect_equal(merged$clusters[["P + R"]]$weights, fresh$clusters[["P + R"]]$weights)

  named <- build_mcml(mc, combine = list(Task = c("P", "R")))
  expect_setequal(names(named$cluster_members), c("Task", "S"))
  # invariant: merging never changes the sequences
  expect_equal(repart_seqs(merged), repart_seqs(mc))
})

test_that("expand = 'all' turns the macro network into the state network", {
  mc <- make_repart_mcml()
  every <- build_mcml(mc, expand = "all")
  states <- sort(unlist(mc$cluster_members, use.names = FALSE))
  expect_setequal(names(every$cluster_members), states)
  # calibration: singleton clusters of a tna mcml = relative transitions
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time")
  expect_equal(every$macro$weights[states, states], net$weights[states, states],
               tolerance = 1e-12, ignore_attr = TRUE)

  one <- build_mcml(mc, expand = "S")
  expect_setequal(names(one$cluster_members),
                  c("P", "R", "cohesion", "emotion", "monitor", "adapt"))
})

test_that("combine then expand compose, in that order", {
  mc <- make_repart_mcml()
  both <- build_mcml(mc, combine = c("P", "R"), expand = "P + R")
  expect_setequal(names(both$cluster_members),
                  c("plan", "consensus", "coregulate", "discuss", "synthesis", "S"))
})

test_that("sequence_plot(combine =) draws what the re-estimated model draws", {
  skip_if_not_installed("ggplot2")
  mc <- make_repart_mcml()
  panel_data <- function(p) {
    do.call(rbind, lapply(attr(p, "panels"), function(q) {
      d <- q$data; d$channel <- as.character(d$channel); d$key <- as.character(d$key); d
    }))
  }
  a <- panel_data(sequence_plot(mc, type = "distribution", trim = 10,
                                combine = c("P", "R")))
  b <- panel_data(sequence_plot(build_mcml(mc, combine = c("P", "R")),
                                type = "distribution", trim = 10))
  a <- a[order(a$channel, a$key, a$time), ]; b <- b[order(b$channel, b$key, b$time), ]
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("re-partitioning rejects what it cannot do", {
  mc <- make_repart_mcml()
  expect_error(build_mcml(mc, clusters = mc$cluster_members, combine = c("P", "R")),
               "either")
  expect_error(build_mcml(mc, combine = c("P", "nope")), "Unknown")
  expect_error(build_mcml(mc, expand = "nope"), "Unknown")
  m <- matrix(runif(16), 4, dimnames = list(letters[1:4], letters[1:4]))
  agg <- build_mcml(m, clusters = list(A = c("a", "b"), B = c("c", "d")))
  expect_error(build_mcml(agg, combine = c("A", "B")),
               class = "nestimate_mcml_no_sequences")
})

test_that("combine/expand on fresh input equal the rebuild and a manual merge", {
  cl <- list(P = c("plan", "consensus"),
             R = c("coregulate", "discuss", "synthesis"),
             S = c("cohesion", "emotion", "monitor", "adapt"))
  one_pass <- build_mcml(group_regulation_long, clusters = cl, combine = c("P", "R"),
                         actor = "Actor", action = "Action", time = "Time")
  two_pass <- build_mcml(build_mcml(group_regulation_long, clusters = cl,
                                    actor = "Actor", action = "Action", time = "Time"),
                         combine = c("P", "R"))
  manual   <- build_mcml(group_regulation_long,
                         clusters = list(`P + R` = c(cl$P, cl$R), S = cl$S),
                         actor = "Actor", action = "Action", time = "Time")
  expect_equal(one_pass$macro$weights, two_pass$macro$weights)
  expect_equal(unclass(one_pass)[c("macro", "clusters", "cluster_members", "meta")],
               unclass(manual)[c("macro", "clusters", "cluster_members", "meta")])

  ex <- build_mcml(group_regulation_long, clusters = cl, expand = "S",
                   actor = "Actor", action = "Action", time = "Time")
  expect_setequal(names(ex$cluster_members),
                  c("P", "R", "cohesion", "emotion", "monitor", "adapt"))
})

test_that("fresh-input combine works for matrices and membership vectors", {
  m <- matrix(c(0, 2, 1, 0,  1, 0, 0, 3,  2, 1, 0, 1,  0, 1, 2, 0), 4, byrow = TRUE,
              dimnames = list(letters[1:4], letters[1:4]))
  cl <- list(A = "a", B = "b", C = c("c", "d"))
  expect_equal(build_mcml(m, clusters = cl, combine = c("A", "B"))$macro$weights,
               build_mcml(m, clusters = list(`A + B` = c("a", "b"), C = c("c", "d")))$macro$weights)

  seqs <- data.frame(T1 = c("a", "c", "b"), T2 = c("b", "d", "a"),
                     T3 = c("c", "c", "d"), T4 = c("d", "a", "c"))
  memb <- c(a = "A", b = "B", c = "C", d = "C")
  expect_equal(build_mcml(seqs, clusters = memb, combine = list(AB = c("A", "B")))$macro$weights,
               build_mcml(seqs, clusters = list(AB = c("a", "b"), C = c("c", "d")))$macro$weights)
})
