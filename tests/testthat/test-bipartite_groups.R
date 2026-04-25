# ---- bipartite_groups() tests --------------------------------------------

# Helpers ------------------------------------------------------------------

.bg_sample_data <- function() {
  data.frame(
    player = c("Alice", "Bob", "Carol", "Alice", "Bob",
               "Dave", "Carol", "Dave", "Eve"),
    session = c("S1", "S1", "S1", "S2", "S2",
                "S3", "S3", "S3", "S3"),
    stringsAsFactors = FALSE
  )
}

# Structure ----------------------------------------------------------------

test_that("returns a net_hypergraph with required fields", {
  hg <- bipartite_groups(.bg_sample_data(), player = "player", group = "session")
  expect_s3_class(hg, "net_hypergraph")
  expect_named(hg, c("hyperedges", "incidence", "nodes", "n_nodes",
                     "n_hyperedges", "size_distribution", "params"))
  expect_equal(hg$n_nodes, 5L)
  expect_equal(hg$n_hyperedges, 3L)
  expect_setequal(hg$nodes, c("Alice", "Bob", "Carol", "Dave", "Eve"))
})

# Hyperedge content --------------------------------------------------------

test_that("each group becomes a hyperedge spanning its players", {
  hg <- bipartite_groups(.bg_sample_data(), player = "player", group = "session")
  # Map back: hyperedges are integer indices into hg$nodes
  members_by_session <- lapply(hg$hyperedges, function(idx) sort(hg$nodes[idx]))
  names(members_by_session) <- colnames(hg$incidence)
  expect_equal(members_by_session[["S1"]], c("Alice", "Bob", "Carol"))
  expect_equal(members_by_session[["S2"]], c("Alice", "Bob"))
  expect_equal(members_by_session[["S3"]], c("Carol", "Dave", "Eve"))
})

# Incidence matrix ---------------------------------------------------------

test_that("incidence is binary by default with correct dimensions and sums", {
  hg <- bipartite_groups(.bg_sample_data(), player = "player", group = "session")
  expect_equal(dim(hg$incidence), c(5L, 3L))
  expect_true(all(hg$incidence %in% c(0L, 1L)))
  # column sums = group sizes
  expect_equal(unname(colSums(hg$incidence)), c(3L, 2L, 3L))
  # row sums = player hyperdegree (num groups they appeared in)
  expect_equal(hg$incidence["Alice", , drop = TRUE], c(S1 = 1L, S2 = 1L, S3 = 0L))
  expect_equal(rowSums(hg$incidence)[["Carol"]], 2L)
})

# Size distribution --------------------------------------------------------

test_that("size_distribution counts hyperedges by size", {
  hg <- bipartite_groups(.bg_sample_data(), player = "player", group = "session")
  expect_equal(hg$size_distribution[["size_2"]], 1L)
  expect_equal(hg$size_distribution[["size_3"]], 2L)
})

# Weighted incidence -------------------------------------------------------

test_that("weight column produces weighted incidence (sum per cell)", {
  d <- data.frame(
    player = c("A", "B", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g2", "g2"),
    n      = c(2, 5, 3, 1, 4),
    stringsAsFactors = FALSE
  )
  hg <- bipartite_groups(d, player = "player", group = "grp", weight = "n")
  # A in g1: 2 + 3 = 5; B in g1: 5; A in g2: 1; B in g2: 4
  expect_equal(hg$incidence["A", "g1"], 5)
  expect_equal(hg$incidence["B", "g1"], 5)
  expect_equal(hg$incidence["A", "g2"], 1)
  expect_equal(hg$incidence["B", "g2"], 4)
  # hyperedges still based on membership (any nonzero), not weight
  expect_setequal(hg$nodes[hg$hyperedges[[1]]], c("A", "B"))
})

# NA handling --------------------------------------------------------------

test_that("rows with NA in player or group are dropped silently", {
  d <- data.frame(
    player = c("A", "B", NA,  "C", "D"),
    grp    = c("g1", "g1", "g1", NA,  "g2"),
    stringsAsFactors = FALSE
  )
  hg <- bipartite_groups(d, player = "player", group = "grp")
  expect_equal(hg$n_nodes, 3L)         # A, B, D (C dropped because group NA)
  expect_setequal(hg$nodes, c("A", "B", "D"))
  expect_equal(hg$n_hyperedges, 2L)    # g1, g2
  # n_observations records post-drop count
  expect_equal(hg$params$n_observations, 3L)
})

test_that("all-NA data raises an error", {
  d <- data.frame(player = NA, grp = NA, stringsAsFactors = FALSE)
  expect_error(bipartite_groups(d, "player", "grp"),
               "No complete observations")
})

# Repeated rows are deduplicated in binary mode ----------------------------

test_that("duplicate rows do not duplicate membership in binary mode", {
  d <- data.frame(
    player = c("A", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g1"),
    stringsAsFactors = FALSE
  )
  hg <- bipartite_groups(d, "player", "grp")
  expect_equal(hg$incidence["A", "g1"], 1L)
  expect_equal(hg$incidence["B", "g1"], 1L)
  expect_equal(hg$n_hyperedges, 1L)
  expect_equal(length(hg$hyperedges[[1]]), 2L)
})

# Single-player groups -----------------------------------------------------

test_that("singleton groups become size-1 hyperedges", {
  d <- data.frame(
    player = c("A", "B", "C"),
    grp    = c("g1", "g2", "g3"),
    stringsAsFactors = FALSE
  )
  hg <- bipartite_groups(d, "player", "grp")
  expect_equal(hg$n_hyperedges, 3L)
  expect_true(all(vapply(hg$hyperedges, length, integer(1)) == 1L))
})

# Validation ---------------------------------------------------------------

test_that("missing player column raises error", {
  d <- data.frame(player = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(bipartite_groups(d, player = "playor", group = "grp"))
})

test_that("missing group column raises error", {
  d <- data.frame(player = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(bipartite_groups(d, player = "player", group = "guppy"))
})

test_that("non-data.frame input rejected", {
  expect_error(bipartite_groups(list(player = "A", grp = "g"), "player", "grp"))
})

# print and summary work --------------------------------------------------

test_that("print and summary work via shared net_hypergraph methods", {
  hg <- bipartite_groups(.bg_sample_data(), player = "player", group = "session")
  expect_invisible(print(hg))
  # summary now returns a tidy node-degree data.frame (visible)
  s <- summary(hg)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("node", "degree"))
})

# Integration with bundled dataset ----------------------------------------

test_that("works on bundled human_long dataset (long-format event data)", {
  data("human_long", package = "Nestimate")
  hg <- bipartite_groups(human_long, player = "code", group = "session_id")
  expect_s3_class(hg, "net_hypergraph")
  expect_gt(hg$n_nodes, 0L)
  expect_gt(hg$n_hyperedges, 0L)
  # Each session is a hyperedge; players = codes appearing in that session
  expect_equal(ncol(hg$incidence), hg$n_hyperedges)
  expect_equal(nrow(hg$incidence), hg$n_nodes)
})
