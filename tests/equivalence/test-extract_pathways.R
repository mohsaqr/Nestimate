# extract_pathways(): one row per actor's pathway, the actions joined in row order, with its length and its first
# and last action. Reference: the same written from the definition in base R. New: no comparison existed.

test_that("unit pathways equal the actions of each actor joined in order", {
  set.seed(1301)
  log <- data.frame(actor = rep(c("p1", "p2", "p3", "p4"), times = c(5L, 8L, 3L, 6L)),
                    action = sample(c("plan", "monitor", "adapt", "discuss"), 22L, replace = TRUE))
  actions <- split(log$action, factor(log$actor, unique(log$actor)))
  expected <- data.frame(actor = names(actions), path = vapply(actions, paste, character(1), collapse = " -> "),
                         length = lengths(actions), opens = vapply(actions, \(a) a[[1L]], character(1)),
                         closes = vapply(actions, \(a) utils::tail(a, 1L), character(1)))
  expect_reference(subset(extract_pathways(log, action = "action", group = "actor"), select = -pathway), expected,
                   reference = "independent: actions of each actor joined in row order, base R",
                   covers = "extract_pathways()", quantity = "unit pathways", tolerance = 0, key = "actor")
})

test_that("unit pathways equal TraMineR's concatenated sequences, lengths and end states", {
  skip_without_reference("TraMineR")
  invisible(lapply(c(1301L, 1302L), \(seed) {
    set.seed(seed)
    log <- data.frame(actor = rep(c("p1", "p2", "p3", "p4"), times = c(5L, 8L, 3L, 6L)),
                      action = sample(c("plan", "monitor", "adapt", "discuss"), 22L, replace = TRUE))
    sequences <- out_seqdef(out_pad(split(log$action, factor(log$actor, unique(log$actor)))))
    lengths <- as.integer(TraMineR::seqlength(sequences))
    cells <- as.matrix(sequences)
    expected <- data.frame(actor = unique(log$actor), path = as.vector(TraMineR::seqconc(sequences, sep = " -> ")),
                           length = lengths, opens = cells[, 1L], closes = cells[cbind(seq_along(lengths), lengths)])
    expect_reference(subset(extract_pathways(log, action = "action", group = "actor"), select = -pathway), expected,
                     reference = "TraMineR", covers = "extract_pathways()",
                     quantity = sprintf("unit pathways, seed %d, TraMineR", seed), tolerance = 0, key = "actor")
  }))
})
