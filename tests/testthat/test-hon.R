testthat::skip_on_cran()

# ===========================================================================
# Tests for build_hon() — Higher-Order Network construction
# ===========================================================================

# --- Helper: simple trajectories with known structure ---
.make_hon_data <- function() {
  data.frame(
    T1 = c("A", "B", "C", "A", "D"),
    T2 = c("B", "A", "A", "B", "A"),
    T3 = c("C", "B", "B", "C", "B"),
    T4 = c("D", "C", "C", "D", "C"),
    T5 = c("A", "D", "D", "A", NA),
    T6 = c("B", "A", "A", "B", NA),
    T7 = c("C", "B", "B", "C", NA),
    stringsAsFactors = FALSE
  )
}

# ===========================================================================
# Section 1: Input validation
# ===========================================================================
test_that("build_hon rejects non-data.frame non-list input", {
  expect_error(build_hon(42), "data.frame or list")
})

test_that("build_hon rejects empty data.frame", {
  expect_error(build_hon(data.frame()), "at least one")
})

test_that("build_hon rejects max_order < 1", {
  expect_error(build_hon(.make_hon_data(), max_order = 0), "max_order")
})

test_that("build_hon rejects min_freq < 1", {
  expect_error(build_hon(.make_hon_data(), min_freq = 0), "min_freq")
})

test_that("build_hon accepts data.frame input", {
  result <- build_hon(.make_hon_data(), min_freq = 1L)
  expect_s3_class(result, "net_hon")
})

test_that("build_hon accepts list input", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  result <- build_hon(trajs, min_freq = 1L)
  expect_s3_class(result, "net_hon")
})

test_that("build_hon strips trailing NAs from data.frame rows", {
  df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"), T3 = c("C", NA),
                   stringsAsFactors = FALSE)
  result <- build_hon(df, min_freq = 1L)
  expect_s3_class(result, "net_hon")
})

test_that("build_hon collapse_repeats removes adjacent duplicates", {
  trajs <- list(c("A", "A", "B", "B", "C"))
  result <- build_hon(trajs, min_freq = 1L, collapse_repeats = TRUE)
  expect_true(result$n_edges > 0)
})

# ===========================================================================
# Section 2: Observation counting
# ===========================================================================
test_that("counts correct for single trajectory A->B->C", {
  trajs <- list(c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 2L)
  expect_equal(count[[.hon_encode("A")]][["B"]], 1L)
  expect_equal(count[[.hon_encode("B")]][["C"]], 1L)
  expect_equal(count[[.hon_encode(c("A", "B"))]][["C"]], 1L)
})

test_that("counts accumulate across multiple trajectories", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  expect_equal(count[[.hon_encode("A")]][["B"]], 3L)
  expect_equal(count[[.hon_encode("B")]][["C"]], 2L)
  expect_equal(count[[.hon_encode("B")]][["D"]], 1L)
})

test_that("max_order limits observation depth", {
  trajs <- list(c("A", "B", "C", "D"))
  count1 <- .hon_build_observations(trajs, max_order = 1L)
  expect_false(is.null(count1[[.hon_encode("A")]]))
  expect_true(is.null(count1[[.hon_encode(c("A", "B"))]]))
  count2 <- .hon_build_observations(trajs, max_order = 2L)
  expect_false(is.null(count2[[.hon_encode(c("A", "B"))]]))
  expect_true(is.null(count2[[.hon_encode(c("A", "B", "C"))]]))
})

test_that("short trajectories contribute only possible orders", {
  trajs <- list(c("X", "Y"))
  count <- .hon_build_observations(trajs, max_order = 5L)
  expect_equal(count[[.hon_encode("X")]][["Y"]], 1L)
  all_keys <- ls(count)
  expect_true(all(.hon_key_len(all_keys) == 1L))
})

# ===========================================================================
# Section 3: Distribution building
# ===========================================================================
test_that("distributions sum to 1 for each source", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 1L)
  for (key in ls(distr)) {
    probs <- distr[[key]]
    if (length(probs) > 0L) {
      expect_equal(sum(probs), 1.0, tolerance = 1e-10)
    }
  }
})

test_that("min_freq filters low-count transitions", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 2L)
  b_distr <- distr[[.hon_encode("B")]]
  expect_true(is.na(b_distr["D"]) || is.null(b_distr["D"]))
  expect_equal(unname(b_distr["C"]), 1.0)
})

test_that("sources with all counts below min_freq get empty distribution", {
  trajs <- list(c("X", "Y"))
  count <- .hon_build_observations(trajs, max_order = 1L)
  distr <- .hon_build_distributions(count, min_freq = 5L)
  expect_equal(length(distr[[.hon_encode("X")]]), 0L)
})

# ===========================================================================
# Section 4: KL-divergence and threshold
# ===========================================================================
test_that("KLD of identical distributions is 0", {
  d <- c(A = 0.5, B = 0.5)
  expect_equal(.hon_kld(d, d), 0.0)
})

test_that("KLD of peaked vs uniform is log2(2) = 1", {
  a <- c(X = 1.0)
  b <- c(X = 0.5, Y = 0.5)
  expect_equal(.hon_kld(a, b), 1.0)
})

test_that("KLD threshold decreases with more data", {
  count_env <- new.env(hash = TRUE, parent = emptyenv())
  count_env[["k"]] <- c(A = 5L, B = 5L)
  t1 <- .hon_kld_threshold(2L, "k", count_env)
  count_env[["k"]] <- c(A = 50L, B = 50L)
  t2 <- .hon_kld_threshold(2L, "k", count_env)
  expect_true(t2 < t1)
})

test_that("KLD threshold increases with order", {
  count_env <- new.env(hash = TRUE, parent = emptyenv())
  count_env[["k"]] <- c(A = 10L, B = 10L)
  t2 <- .hon_kld_threshold(2L, "k", count_env)
  t3 <- .hon_kld_threshold(3L, "k", count_env)
  expect_true(t3 > t2)
})

# ===========================================================================
# Section 5: End-to-end pipeline
# ===========================================================================
test_that("build_hon returns correct structure", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_true(is.matrix(result$matrix))
  expect_true(is.data.frame(result$ho_edges))
  expect_true(is.data.frame(result$nodes))
  expect_true(result$directed)
  expect_true(result$n_nodes > 0L)
  expect_true(result$n_edges > 0L)
  # Check unified column structure
  expect_true(all(c("path", "from", "to", "count", "probability",
                     "from_order", "to_order") %in% names(result$ho_edges)))
})

test_that("build_hon nodes use arrow notation", {
  trajs <- list(c("A", "B", "C"))
  result <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  # First-order nodes should be simple state names (no arrows)
  expect_true(all(result$nodes$label %in% c("A", "B", "C")))
})

test_that("sequence_to_node produces readable arrow notation", {
  expect_equal(.hon_sequence_to_node("A"), "A")
  expect_equal(.hon_sequence_to_node(c("A", "B")), "A -> B")
  expect_equal(.hon_sequence_to_node(c("X", "A", "B")), "X -> A -> B")
})

test_that("print and summary work without error", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_output(print(result), "Higher-Order Network")
  expect_output(summary(result), "Summary")
})

test_that("edge probabilities are in (0, 1]", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_true(all(result$ho_edges$probability > 0))
  expect_true(all(result$ho_edges$probability <= 1))
})

test_that("encode/decode roundtrip preserves tuple", {
  tup <- c("alpha", "beta", "gamma")
  expect_equal(.hon_decode(.hon_encode(tup)), tup)
})

test_that("key_len returns correct lengths", {
  keys <- c(.hon_encode("A"), .hon_encode(c("A", "B")),
            .hon_encode(c("X", "Y", "Z")))
  expect_equal(.hon_key_len(keys), c(1L, 2L, 3L))
})

test_that("build_hon with max_order=1 produces only first-order nodes", {
  trajs <- list(c("A", "B", "C", "D"))
  result <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  # All nodes should have order 1 (no arrows in name)
  node_orders <- vapply(result$nodes$label, function(nd) {
    length(strsplit(nd, " -> ", fixed = TRUE)[[1L]])
  }, integer(1L))
  expect_true(all(node_orders == 1L))
})

test_that("adjacency matrix dimensions match n_nodes", {
  result <- build_hon(.make_hon_data(), max_order = 2L, min_freq = 1L)
  expect_equal(nrow(result$matrix), result$n_nodes)
  expect_equal(ncol(result$matrix), result$n_nodes)
})

# ===========================================================================
# Section 6: pyHON equivalence validation
# ===========================================================================

# --- Helper: convert pipe notation to arrow notation ---
.pipe_to_arrow <- function(pipe_node) {
  # "A|" -> "A", "B|A" -> "A -> B", "C|B.A" -> "A -> B -> C"
  parts <- strsplit(pipe_node, "|", fixed = TRUE)[[1L]]
  current <- parts[1L]
  if (length(parts) < 2L || parts[2L] == "") return(current)
  history <- strsplit(parts[2L], ".", fixed = TRUE)[[1L]]
  paste(c(rev(history), current), collapse = " -> ")
}

# --- Helper: run pyHON via reticulate and return edges + rules ---
.run_pyhon <- function(trajectories, max_order = 5L, min_freq = 1L) {
  skip_if_not_installed("reticulate")
  Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3")
  skip_if(!reticulate::py_available(initialize = TRUE),
          "Python not available")

  reticulate::py_run_string("
from collections import defaultdict, Counter
import math

def run_buildhon(trajectories, max_order, min_support):
    # --- Stage 1: Extract Rules ---
    Count = {}
    Distribution = defaultdict(dict)
    SourceToExtSource = {}
    Rules = defaultdict(dict)

    # BuildObservations
    for traj in trajectories:
        n = len(traj)
        for order in range(2, max_order + 2):
            if order > n:
                break
            for start in range(n - order + 1):
                subseq = tuple(traj[start:start + order])
                Target = subseq[-1]
                Source = subseq[:-1]
                if Source not in Count:
                    Count[Source] = Counter()
                Count[Source][Target] += 1

    # BuildDistributions
    for Source in Count:
        for Target in list(Count[Source].keys()):
            if Count[Source][Target] < min_support:
                Count[Source][Target] = 0
        for Target in Count[Source]:
            if Count[Source][Target] > 0:
                Distribution[Source][Target] = Count[Source][Target] / sum(Count[Source].values())

    # BuildSourceToExtSource
    for source in Distribution:
        if len(source) > 1:
            NewOrder = len(source)
            for starting in range(1, len(source)):
                curr = source[starting:]
                if curr not in SourceToExtSource:
                    SourceToExtSource[curr] = {}
                if NewOrder not in SourceToExtSource[curr]:
                    SourceToExtSource[curr][NewOrder] = set()
                SourceToExtSource[curr][NewOrder].add(source)

    def ExtendSource(Curr, NewOrder):
        if Curr in SourceToExtSource:
            if NewOrder in SourceToExtSource[Curr]:
                return SourceToExtSource[Curr][NewOrder]
        return []

    def KLD(a, b):
        divergence = 0
        for target in a:
            pa = a.get(target, 0)
            pb = b.get(target, 0)
            if pa > 0 and pb > 0:
                divergence += pa * math.log(pa / pb, 2)
            elif pa > 0 and pb == 0:
                divergence = float('inf')
        return divergence

    def KLDThreshold(NewOrder, ExtSource):
        return NewOrder / math.log(1 + sum(Count[ExtSource].values()), 2)

    def AddToRules(Source):
        if len(Source) > 0:
            Rules[Source] = dict(Distribution[Source])
            PrevSource = Source[:-1]
            AddToRules(PrevSource)

    def ExtendRule(Valid, Curr, order, MaxOrder):
        if order >= MaxOrder:
            AddToRules(Valid)
        else:
            Distr = Distribution[Valid]
            NewOrder = order + 1
            Extended = ExtendSource(Curr, NewOrder)
            if len(Extended) == 0:
                AddToRules(Valid)
            else:
                for ExtSource in Extended:
                    ExtDistr = Distribution[ExtSource]
                    if KLD(ExtDistr, Distr) > KLDThreshold(NewOrder, ExtSource):
                        ExtendRule(ExtSource, ExtSource, NewOrder, MaxOrder)
                    else:
                        ExtendRule(Valid, ExtSource, NewOrder, MaxOrder)

    # GenerateAllRules
    for Source in list(Distribution.keys()):
        if len(Source) == 1:
            AddToRules(Source)
            ExtendRule(Source, Source, 1, max_order)

    # --- Stage 2: Build Network ---
    Graph = defaultdict(dict)
    SortedSource = sorted(Rules, key=lambda x: len(x))
    for source in SortedSource:
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                PrevSource = source[:-1]
                PrevTarget = (source[-1],)
                if PrevSource not in Graph or source not in Graph[PrevSource]:
                    if PrevTarget in Graph.get(PrevSource, {}):
                        Graph[PrevSource][source] = Graph[PrevSource][PrevTarget]
                        del Graph[PrevSource][PrevTarget]

    # RewireTails
    ToAdd = []
    ToRemove = []
    for source in list(Graph.keys()):
        for target in list(Graph[source].keys()):
            if len(target) == 1:
                NewTarget = source + target
                while len(NewTarget) > 1:
                    if NewTarget in Graph:
                        ToAdd.append((source, NewTarget, Graph[source][target]))
                        ToRemove.append((source, target))
                        break
                    else:
                        NewTarget = NewTarget[1:]
    for (source, target, weight) in ToAdd:
        Graph[source][target] = weight
    for (source, target) in ToRemove:
        if target in Graph.get(source, {}):
            del Graph[source][target]

    # Convert to edge list
    def SequenceToNode(seq):
        curr = seq[-1]
        node = curr + '|'
        seq = seq[:-1]
        while len(seq) > 0:
            curr = seq[-1]
            node = node + curr + '.'
            seq = seq[:-1]
        if node[-1] == '.':
            return node[:-1]
        else:
            return node

    edges = []
    for source in Graph:
        for target in Graph[source]:
            edges.append({
                'from': SequenceToNode(source),
                'to': SequenceToNode(target),
                'weight': Graph[source][target]
            })

    return {'edges': edges, 'n_edges': len(edges)}
  ")

  py_trajs <- reticulate::r_to_py(trajectories)
  reticulate::py$run_buildhon(
    py_trajs,
    max_order = reticulate::r_to_py(as.integer(max_order)),
    min_support = reticulate::r_to_py(as.integer(min_freq))
  )
}

test_that("pyHON equivalence: simple 4-state trajectories", {
  trajs <- list(
    c("A", "B", "C", "D", "B", "A"),
    c("A", "B", "C", "D", "B", "A"),
    c("A", "B", "C", "D", "B", "A"),
    c("A", "B", "C", "D", "B", "A"),
    c("A", "B", "C", "D", "B", "A")
  )
  py <- .run_pyhon(trajs, max_order = 5L, min_freq = 1L)
  r <- build_hon(trajs, max_order = 5L, min_freq = 1L, method = "hon")

  expect_equal(r$n_edges, py$n_edges,
    info = sprintf("Edge count: R=%d, Python=%d", r$n_edges, py$n_edges))

  for (i in seq_along(py$edges)) {
    edge <- py$edges[[i]]
    py_from <- .pipe_to_arrow(edge$from)
    py_to_node <- .pipe_to_arrow(edge[["to"]])
    # R edges$to is next state (last element of to-node)
    py_next_state <- tail(strsplit(py_to_node, " -> ", fixed = TRUE)[[1L]], 1L)
    r_match <- r$ho_edges[r$ho_edges$from == py_from & r$ho_edges$to == py_next_state, ]
    expect_equal(nrow(r_match), 1L,
      info = sprintf("Missing edge: %s -> %s", py_from, py_next_state))
    if (nrow(r_match) == 1L) {
      expect_equal(r_match$probability, edge$weight, tolerance = 1e-10,
        info = sprintf("Prob mismatch: %s -> %s", py_from, py_next_state))
    }
  }
})

test_that("pyHON equivalence: higher-order dependency dataset", {
  trajs <- list(
    c("A", "B", "C", "A", "B", "C", "A"),
    c("D", "B", "A", "D", "B", "A", "D"),
    c("A", "B", "C", "A", "B", "C", "A"),
    c("D", "B", "A", "D", "B", "A", "D"),
    c("A", "B", "C", "D", "B", "A", "D"),
    c("A", "B", "C", "A", "B", "C", "A"),
    c("D", "B", "A", "D", "B", "A", "D")
  )
  py <- .run_pyhon(trajs, max_order = 5L, min_freq = 1L)
  r <- build_hon(trajs, max_order = 5L, min_freq = 1L, method = "hon")

  expect_equal(r$n_edges, py$n_edges,
    info = sprintf("Edge count: R=%d, Python=%d", r$n_edges, py$n_edges))

  for (i in seq_along(py$edges)) {
    edge <- py$edges[[i]]
    py_from <- .pipe_to_arrow(edge$from)
    py_next <- tail(strsplit(.pipe_to_arrow(edge[["to"]]), " -> ",
                             fixed = TRUE)[[1L]], 1L)
    r_match <- r$ho_edges[r$ho_edges$from == py_from & r$ho_edges$to == py_next, ]
    expect_equal(nrow(r_match), 1L,
      info = sprintf("Missing: %s -> %s", py_from, py_next))
    if (nrow(r_match) == 1L) {
      expect_equal(r_match$probability, edge$weight, tolerance = 1e-10,
        info = sprintf("Prob: %s -> %s", py_from, py_next))
    }
  }
})

test_that("pyHON equivalence: 5-state diverse trajectories", {
  trajs <- list(
    c("E", "A", "B", "C", "D", "E", "A"),
    c("B", "C", "D", "A", "E", "B", "C"),
    c("A", "B", "C", "E", "D", "A", "B"),
    c("D", "E", "A", "B", "C", "D", "E"),
    c("C", "D", "A", "B", "E", "C", "D"),
    c("A", "B", "C", "D", "E", "A", "B"),
    c("E", "D", "C", "B", "A", "E", "D"),
    c("B", "A", "E", "D", "C", "B", "A"),
    c("A", "B", "C", "D", "A", "B", "C"),
    c("D", "C", "B", "A", "D", "C", "B")
  )
  py <- .run_pyhon(trajs, max_order = 5L, min_freq = 1L)
  r <- build_hon(trajs, max_order = 5L, min_freq = 1L, method = "hon")

  expect_equal(r$n_edges, py$n_edges,
    info = sprintf("Edge count: R=%d, Python=%d", r$n_edges, py$n_edges))

  for (i in seq_along(py$edges)) {
    edge <- py$edges[[i]]
    py_from <- .pipe_to_arrow(edge$from)
    py_next <- tail(strsplit(.pipe_to_arrow(edge[["to"]]), " -> ",
                             fixed = TRUE)[[1L]], 1L)
    r_match <- r$ho_edges[r$ho_edges$from == py_from & r$ho_edges$to == py_next, ]
    expect_equal(nrow(r_match), 1L,
      info = sprintf("Missing: %s -> %s", py_from, py_next))
    if (nrow(r_match) == 1L) {
      expect_equal(r_match$probability, edge$weight, tolerance = 1e-10)
    }
  }
})

test_that("pyHON equivalence: with min_freq = 3", {
  trajs <- list(
    c("A", "B", "C", "D"),
    c("A", "B", "D", "C"),
    c("A", "B", "C", "D"),
    c("D", "C", "B", "A"),
    c("A", "B", "C", "D"),
    c("C", "B", "A", "D")
  )
  py <- .run_pyhon(trajs, max_order = 5L, min_freq = 3L)
  r <- build_hon(trajs, max_order = 5L, min_freq = 3L, method = "hon")

  expect_equal(r$n_edges, py$n_edges,
    info = sprintf("Edge count: R=%d, Python=%d", r$n_edges, py$n_edges))

  for (i in seq_along(py$edges)) {
    edge <- py$edges[[i]]
    py_from <- .pipe_to_arrow(edge$from)
    py_next <- tail(strsplit(.pipe_to_arrow(edge[["to"]]), " -> ",
                             fixed = TRUE)[[1L]], 1L)
    r_match <- r$ho_edges[r$ho_edges$from == py_from & r$ho_edges$to == py_next, ]
    expect_equal(nrow(r_match), 1L,
      info = sprintf("Missing: %s -> %s", py_from, py_next))
    if (nrow(r_match) == 1L) {
      expect_equal(r_match$probability, edge$weight, tolerance = 1e-10)
    }
  }
})

test_that("pyHON equivalence: max_order = 2", {
  trajs <- list(
    c("A", "B", "C", "A", "B", "D"),
    c("A", "B", "C", "A", "B", "C"),
    c("D", "B", "A", "D", "B", "A"),
    c("A", "B", "C", "D", "B", "A"),
    c("A", "B", "D", "A", "B", "C")
  )
  py <- .run_pyhon(trajs, max_order = 2L, min_freq = 1L)
  r <- build_hon(trajs, max_order = 2L, min_freq = 1L, method = "hon")

  expect_equal(r$n_edges, py$n_edges,
    info = sprintf("Edge count: R=%d, Python=%d", r$n_edges, py$n_edges))

  for (i in seq_along(py$edges)) {
    edge <- py$edges[[i]]
    py_from <- .pipe_to_arrow(edge$from)
    py_next <- tail(strsplit(.pipe_to_arrow(edge[["to"]]), " -> ",
                             fixed = TRUE)[[1L]], 1L)
    r_match <- r$ho_edges[r$ho_edges$from == py_from & r$ho_edges$to == py_next, ]
    expect_equal(nrow(r_match), 1L,
      info = sprintf("Missing: %s -> %s", py_from, py_next))
    if (nrow(r_match) == 1L) {
      expect_equal(r_match$probability, edge$weight, tolerance = 1e-10)
    }
  }
})

# ===========================================================================
# Section 7: HON+ internals
# ===========================================================================

# ---------------------------------------------------------------------------
# Task 1: .honp_max_divergence
# ---------------------------------------------------------------------------

test_that("honp_max_divergence puts all mass on least probable target", {
  d <- c(A = 0.6, B = 0.3, C = 0.1)
  result <- .honp_max_divergence(d)
  expect_equal(length(result), 1L)
  expect_equal(names(result), "C")
  expect_equal(unname(result), 1.0)
})

test_that("honp_max_divergence handles single-target distribution", {
  d <- c(X = 1.0)
  result <- .honp_max_divergence(d)
  expect_equal(names(result), "X")
  expect_equal(unname(result), 1.0)
})

test_that("honp_max_divergence handles ties by picking first minimum", {
  d <- c(A = 0.5, B = 0.5)
  result <- .honp_max_divergence(d)
  expect_equal(length(result), 1L)
  expect_equal(unname(result), 1.0)
})

# ---------------------------------------------------------------------------
# Task 2: .honp_build_order1
# ---------------------------------------------------------------------------

test_that("honp_build_order1 produces correct counts and starting_points", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "A"))
  result <- .honp_build_order1(trajs, min_freq = 1L)

  expect_true(is.list(result))
  expect_true(all(c("count", "distr", "starting_points") %in% names(result)))

  count_a <- result$count[[.hon_encode("A")]]
  expect_equal(unname(count_a["B"]), 2L)

  distr_a <- result$distr[[.hon_encode("A")]]
  expect_equal(sum(distr_a), 1.0)

  sp_a <- result$starting_points[[.hon_encode("A")]]
  expect_true(length(sp_a) >= 2L)
})

test_that("honp_build_order1 applies min_freq filtering", {
  trajs <- list(c("A", "B", "C", "B", "C", "B"))
  result <- .honp_build_order1(trajs, min_freq = 2L)

  distr_b <- result$distr[[.hon_encode("B")]]
  expect_true(all(distr_b > 0))

  count_b <- result$count[[.hon_encode("B")]]
  expect_true(any(count_b == 0L) || all(count_b >= 2L))
})

# ---------------------------------------------------------------------------
# Task 3: .honp_extend_observation
# ---------------------------------------------------------------------------

test_that("honp_extend_observation builds higher-order counts lazily", {
  trajs <- list(c("X", "A", "B", "C"), c("Y", "A", "B", "D"))
  o1 <- .honp_build_order1(trajs, min_freq = 1L)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())

  key_a <- .hon_encode("A")
  .honp_extend_observation(key_a, source_to_ext, o1$count, o1$distr,
                           o1$starting_points, trajs, 1L)

  ext_keys <- ls(source_to_ext[[key_a]])
  expect_true(length(ext_keys) >= 1L)

  for (ek in ext_keys) {
    expect_false(is.null(o1$distr[[ek]]))
  }
})

test_that("honp_extend_observation records starting_points for extensions", {
  trajs <- list(c("X", "A", "B"), c("Y", "A", "B"))
  o1 <- .honp_build_order1(trajs, min_freq = 1L)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())

  key_a <- .hon_encode("A")
  .honp_extend_observation(key_a, source_to_ext, o1$count, o1$distr,
                           o1$starting_points, trajs, 1L)

  key_xa <- .hon_encode(c("X", "A"))
  sp <- o1$starting_points[[key_xa]]
  expect_true(!is.null(sp))
  expect_true(length(sp) >= 1L)
})

# ---------------------------------------------------------------------------
# Task 4: .honp_extend_source_fast
# ---------------------------------------------------------------------------

test_that("honp_extend_source_fast returns extensions lazily", {
  trajs <- list(c("X", "A", "B"), c("Y", "A", "B"))
  o1 <- .honp_build_order1(trajs, min_freq = 1L)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())

  key_a <- .hon_encode("A")
  exts <- .honp_extend_source_fast(key_a, source_to_ext, o1$count,
                                    o1$distr, o1$starting_points, trajs, 1L)
  expect_true(length(exts) >= 1L)

  exts2 <- .honp_extend_source_fast(key_a, source_to_ext, o1$count,
                                     o1$distr, o1$starting_points, trajs, 1L)
  expect_equal(sort(exts), sort(exts2))
})

# ---------------------------------------------------------------------------
# Task 5: .honp_extend_rule and .honp_add_to_rules
# ---------------------------------------------------------------------------

test_that("honp_extend_rule finds rules with MaxDivergence pre-check", {
  # 4 obs per branch: threshold = 2/log2(5) ≈ 0.86 < KLD = 1.0  → extends
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  o1 <- .honp_build_order1(trajs, min_freq = 1L)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())
  rules <- new.env(hash = TRUE, parent = emptyenv())

  key_a <- .hon_encode("A")
  .honp_add_to_rules(key_a, o1$distr, rules, source_to_ext, o1$count,
                      o1$starting_points, trajs, 1L)
  .honp_extend_rule(key_a, key_a, 1L, 99L, o1$distr, o1$count,
                    source_to_ext, o1$starting_points, trajs, 1L, rules)

  rule_keys <- ls(rules)
  rule_lens <- .hon_key_len(rule_keys)
  expect_true(any(rule_lens == 2L))
})

test_that("honp_extend_rule prunes when MaxDivergence below threshold", {
  trajs <- list(
    c("X", "A", "B"), c("Y", "A", "B"),
    c("X", "A", "B"), c("Y", "A", "B")
  )
  o1 <- .honp_build_order1(trajs, min_freq = 1L)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())
  rules <- new.env(hash = TRUE, parent = emptyenv())

  key_a <- .hon_encode("A")
  .honp_add_to_rules(key_a, o1$distr, rules, source_to_ext, o1$count,
                      o1$starting_points, trajs, 1L)
  .honp_extend_rule(key_a, key_a, 1L, 99L, o1$distr, o1$count,
                    source_to_ext, o1$starting_points, trajs, 1L, rules)

  rule_keys <- ls(rules)
  rule_lens <- .hon_key_len(rule_keys)
  expect_true(all(rule_lens == 1L))
})

# ---------------------------------------------------------------------------
# Task 6: .honp_extract_rules
# ---------------------------------------------------------------------------

test_that("honp_extract_rules produces rules from trajectories", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("A", "B", "C"))
  result <- .honp_extract_rules(trajs, max_order = 99L, min_freq = 1L)
  rules <- result$rules

  rule_keys <- ls(rules)
  expect_true(length(rule_keys) > 0L)

  expect_false(is.null(rules[[.hon_encode("A")]]))
  expect_false(is.null(rules[[.hon_encode("B")]]))
})

# ===========================================================================
# Section 8: build_hon() method parameter
# ===========================================================================
test_that("build_hon accepts method parameter", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"))
  r1 <- build_hon(trajs, max_order = 2L, min_freq = 1L, method = "hon")
  expect_s3_class(r1, "net_hon")

  r2 <- build_hon(trajs, max_order = 2L, min_freq = 1L, method = "hon+")
  expect_s3_class(r2, "net_hon")
})

test_that("build_hon default method is hon+", {
  trajs <- list(c("A", "B", "C"))
  r <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  expect_s3_class(r, "net_hon")
})

test_that("build_hon rejects invalid method", {
  trajs <- list(c("A", "B", "C"))
  expect_error(build_hon(trajs, method = "invalid"), "arg")
})

# ===========================================================================
# Section 9: HON vs HON+ equivalence
# ===========================================================================
test_that("hon and hon+ produce identical networks on simple data", {
  trajs <- list(
    c("A", "B", "C", "D"), c("A", "B", "D", "C"),
    c("A", "C", "B", "D"), c("D", "C", "B", "A")
  )
  r_hon  <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon")
  r_honp <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")

  expect_equal(nrow(r_hon$ho_edges), nrow(r_honp$ho_edges))

  e1 <- r_hon$ho_edges[order(r_hon$ho_edges$from, r_hon$ho_edges$to),
                    c("from", "to", "probability")]
  e2 <- r_honp$ho_edges[order(r_honp$ho_edges$from, r_honp$ho_edges$to),
                     c("from", "to", "probability")]
  rownames(e1) <- rownames(e2) <- NULL
  expect_equal(e1$from, e2$from)
  expect_equal(e1$to, e2$to)
  expect_equal(e1$probability, e2$probability, tolerance = 1e-10)
})

test_that("hon and hon+ match on higher-order dependency data", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"),
    c("A", "B", "A"), c("B", "A", "B")
  )
  r_hon  <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon")
  r_honp <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")

  e1 <- r_hon$ho_edges[order(r_hon$ho_edges$from, r_hon$ho_edges$to), ]
  e2 <- r_honp$ho_edges[order(r_honp$ho_edges$from, r_honp$ho_edges$to), ]
  expect_equal(nrow(e1), nrow(e2))
  expect_equal(e1$probability, e2$probability, tolerance = 1e-10)
})

test_that("hon and hon+ match on 5-state data with min_freq filtering", {
  set.seed(42)
  states <- c("A", "B", "C", "D", "E")
  trajs <- lapply(1:50, function(i) sample(states, sample(5:15, 1), replace = TRUE))
  r_hon  <- build_hon(trajs, max_order = 3L, min_freq = 3L, method = "hon")
  r_honp <- build_hon(trajs, max_order = 3L, min_freq = 3L, method = "hon+")

  e1 <- r_hon$ho_edges[order(r_hon$ho_edges$from, r_hon$ho_edges$to),
                    c("from", "to", "probability")]
  e2 <- r_honp$ho_edges[order(r_honp$ho_edges$from, r_honp$ho_edges$to),
                     c("from", "to", "probability")]
  rownames(e1) <- rownames(e2) <- NULL
  expect_equal(e1, e2, tolerance = 1e-10)
})

# ===========================================================================
# Section 10: pyHON+ (BuildRulesFastParameterFree) equivalence
# ===========================================================================

.run_pyhon_plus <- function(trajectories, max_order, min_freq) {
  skip_if_not_installed("reticulate")
  Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/bin/python3")
  skip_if_not(reticulate::py_available(initialize = TRUE),
              "Python not available")

  # Write Python code for HON+ algorithm inline
  reticulate::py_run_string("
import math
from collections import defaultdict

def run_hon_plus(trajectories, max_order, min_freq):
    Count = defaultdict(lambda: defaultdict(int))
    Distribution = defaultdict(dict)
    Rules = defaultdict(dict)
    SourceToExtSource = defaultdict(set)
    StartingPoints = defaultdict(set)

    # Build order 1
    for ti in range(len(trajectories)):
        traj = trajectories[ti]
        for index in range(len(traj) - 1):
            Source = tuple(traj[index:index+1])
            Target = traj[index+1]
            Count[Source][Target] += 1
            StartingPoints[Source].add((ti, index))

    for Source in list(Count.keys()):
        if len(Source) == 1:
            for Target in list(Count[Source].keys()):
                if Count[Source][Target] < min_freq:
                    Count[Source][Target] = 0
            total = sum(Count[Source].values())
            if total > 0:
                for Target in Count[Source]:
                    if Count[Source][Target] > 0:
                        Distribution[Source][Target] = 1.0 * Count[Source][Target] / total

    def KLD(a, b):
        d = 0
        for t in a:
            pa = a.get(t, 0)
            pb = b.get(t, 0)
            if pa > 0 and pb > 0:
                d += pa * math.log(pa / pb, 2)
            elif pa > 0:
                d = float('inf')
        return d

    def KLDThreshold(no, es):
        total = sum(Count[es].values())
        if total == 0:
            return float('inf')
        return no / math.log(1 + total, 2)

    def MaxDivergence(Distr):
        mk = sorted(Distr, key=Distr.__getitem__)
        return {mk[0]: 1}

    def ExtendObservation(Source):
        if len(Source) > 1:
            if Source[1:] not in Count or len(Count[Source]) == 0:
                ExtendObservation(Source[1:])
        order = len(Source)
        C = defaultdict(lambda: defaultdict(int))
        for ti, index in StartingPoints[Source]:
            traj = trajectories[ti]
            if index - 1 >= 0 and index + order < len(traj):
                ExtSource = tuple(traj[index-1:index+order])
                Target = traj[index+order]
                C[ExtSource][Target] += 1
                StartingPoints[ExtSource].add((ti, index - 1))
        if len(C) == 0:
            return
        for s in C:
            for t in C[s]:
                if C[s][t] < min_freq:
                    C[s][t] = 0
                Count[s][t] += C[s][t]
            total = sum(C[s].values())
            if total > 0:
                for t in C[s]:
                    if C[s][t] > 0:
                        Distribution[s][t] = 1.0 * C[s][t] / total
                        SourceToExtSource[s[1:]].add(s)

    def ExtendSourceFast(Curr):
        if Curr in SourceToExtSource:
            return SourceToExtSource[Curr]
        else:
            ExtendObservation(Curr)
            if Curr in SourceToExtSource:
                return SourceToExtSource[Curr]
            else:
                return []

    def AddToRules(Source):
        for order in range(1, len(Source)+1):
            s = Source[0:order]
            if s not in Distribution or len(Distribution[s]) == 0:
                ExtendSourceFast(s[1:])
            Rules[s] = Distribution[s]

    def ExtendRule(Valid, Curr, order):
        if order >= max_order:
            AddToRules(Valid)
        else:
            Distr = Distribution[Valid]
            CurrDistr = Distribution.get(Curr, {})
            if len(CurrDistr) == 0:
                AddToRules(Valid)
                return
            if KLD(MaxDivergence(CurrDistr), Distr) < KLDThreshold(order+1, Curr):
                AddToRules(Valid)
            else:
                Extended = ExtendSourceFast(Curr)
                if len(Extended) == 0:
                    AddToRules(Valid)
                else:
                    for ext in Extended:
                        ed = Distribution.get(ext, {})
                        if len(ed) > 0 and KLD(ed, Distr) > KLDThreshold(order+1, ext):
                            ExtendRule(ext, ext, order+1)
                        else:
                            ExtendRule(Valid, ext, order+1)

    for Source in tuple(Distribution.keys()):
        AddToRules(Source)
        ExtendRule(Source, Source, 1)

    # Build network (same as BuildNetwork.py)
    Graph = defaultdict(dict)

    def SequenceToNode(seq):
        l = seq[-1]
        if len(seq) == 1:
            return l + '|'
        return l + '|' + '.'.join(reversed(seq[:-1]))

    def Rewire(source, target):
        ps = source[:-1]
        pt = (source[-1],)
        if ps in Graph and source not in Graph[ps]:
            if pt in Graph.get(ps, {}):
                Graph[ps][source] = Graph[ps][pt]
                del Graph[ps][pt]

    SortedSource = sorted(Rules, key=lambda x: len(x))
    for source in SortedSource:
        for target in Rules[source]:
            Graph[source][(target,)] = Rules[source][target]
            if len(source) > 1:
                Rewire(source, (target,))

    # Tail rewiring
    ta = []
    tr = []
    for source in Graph:
        for target in list(Graph[source].keys()):
            if len(target) == 1:
                nt = source + target
                while len(nt) > 1:
                    if nt in Graph:
                        ta.append((source, nt, Graph[source][target]))
                        tr.append((source, target))
                        break
                    nt = nt[1:]
    for s, t, w in ta:
        Graph[s][t] = w
    for s, t in tr:
        if t in Graph.get(s, {}):
            del Graph[s][t]

    return [{'from': SequenceToNode(s), 'to': SequenceToNode(t), 'weight': Graph[s][t]}
            for s in Graph for t in Graph[s]]
")

  reticulate::py_set_attr(reticulate::py, "trajectories", trajectories)
  reticulate::py_run_string(sprintf(
    "edges = run_hon_plus(trajectories, %d, %d)", max_order, min_freq))
  py_edges <- reticulate::py_to_r(
    reticulate::py_get_attr(reticulate::py, "edges"))

  if (length(py_edges) == 0L) {
    return(data.frame(from = character(0L), to = character(0L),
                      weight = numeric(0L), stringsAsFactors = FALSE))
  }

  do.call(rbind, lapply(py_edges, function(e) {
    data.frame(from = e[["from"]], to = e[["to"]], weight = e[["weight"]],
               stringsAsFactors = FALSE)
  }))
}

test_that("hon+ matches pyHON+ on 4-state data", {
  trajs <- list(
    c("A", "B", "C", "D"), c("A", "B", "D", "C"),
    c("B", "C", "A", "D"), c("D", "A", "B", "C")
  )
  r <- build_hon(trajs, max_order = 99L, min_freq = 1L, method = "hon+")
  py <- .run_pyhon_plus(trajs, 99L, 1L)

  # Convert Python pipe notation to arrow notation for comparison
  py$from <- vapply(py$from, .pipe_to_arrow, character(1L))
  py$to <- vapply(py$to, function(x) {
    arrow <- .pipe_to_arrow(x)
    tail(strsplit(arrow, " -> ", fixed = TRUE)[[1L]], 1L)
  }, character(1L))

  expect_equal(nrow(r$ho_edges), nrow(py),
    info = sprintf("Edge count: R=%d, Python=%d", nrow(r$ho_edges), nrow(py)))
  merged <- merge(
    r$ho_edges[, c("from", "to", "probability")],
    py, by = c("from", "to"), suffixes = c("_r", "_py"))
  expect_equal(nrow(merged), nrow(r$ho_edges),
    info = "Not all edges matched by from/to")
  expect_equal(merged$probability, merged$weight, tolerance = 1e-10)
})

test_that("hon+ matches pyHON+ on higher-order dependency data", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"),
    c("A", "B", "A"), c("B", "A", "B")
  )
  r <- build_hon(trajs, max_order = 99L, min_freq = 1L, method = "hon+")
  py <- .run_pyhon_plus(trajs, 99L, 1L)

  py$from <- vapply(py$from, .pipe_to_arrow, character(1L))
  py$to <- vapply(py$to, function(x) {
    arrow <- .pipe_to_arrow(x)
    tail(strsplit(arrow, " -> ", fixed = TRUE)[[1L]], 1L)
  }, character(1L))

  expect_equal(nrow(r$ho_edges), nrow(py),
    info = sprintf("Edge count: R=%d, Python=%d", nrow(r$ho_edges), nrow(py)))
  merged <- merge(
    r$ho_edges[, c("from", "to", "probability")],
    py, by = c("from", "to"), suffixes = c("_r", "_py"))
  expect_equal(nrow(merged), nrow(r$ho_edges),
    info = "Not all edges matched by from/to")
  expect_equal(merged$probability, merged$weight, tolerance = 1e-10)
})

test_that("hon+ matches pyHON+ with min_freq filtering", {
  trajs <- list(
    c("A", "B", "C", "D"), c("A", "B", "D", "C"),
    c("A", "B", "C", "D"), c("D", "C", "B", "A"),
    c("A", "B", "C", "D"), c("C", "B", "A", "D")
  )
  r <- build_hon(trajs, max_order = 99L, min_freq = 3L, method = "hon+")
  py <- .run_pyhon_plus(trajs, 99L, 3L)

  py$from <- vapply(py$from, .pipe_to_arrow, character(1L))
  py$to <- vapply(py$to, function(x) {
    arrow <- .pipe_to_arrow(x)
    tail(strsplit(arrow, " -> ", fixed = TRUE)[[1L]], 1L)
  }, character(1L))

  expect_equal(nrow(r$ho_edges), nrow(py),
    info = sprintf("Edge count: R=%d, Python=%d", nrow(r$ho_edges), nrow(py)))
  merged <- merge(
    r$ho_edges[, c("from", "to", "probability")],
    py, by = c("from", "to"), suffixes = c("_r", "_py"))
  expect_equal(nrow(merged), nrow(r$ho_edges),
    info = "Not all edges matched by from/to")
  expect_equal(merged$probability, merged$weight, tolerance = 1e-10)
})

test_that("hon+ matches pyHON+ on group_regulation subset", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  gr <- group_regulation[1:100, ]
  trajs <- lapply(seq_len(nrow(gr)), function(i) {
    row <- as.character(gr[i, ])
    row[!is.na(row)]
  })

  r <- build_hon(trajs, max_order = 3L, min_freq = 5L, method = "hon+")
  py <- .run_pyhon_plus(trajs, 3L, 5L)

  py$from <- vapply(py$from, .pipe_to_arrow, character(1L))
  py$to <- vapply(py$to, function(x) {
    arrow <- .pipe_to_arrow(x)
    tail(strsplit(arrow, " -> ", fixed = TRUE)[[1L]], 1L)
  }, character(1L))

  expect_equal(nrow(r$ho_edges), nrow(py),
    info = sprintf("Edge count: R=%d, Python=%d", nrow(r$ho_edges), nrow(py)))
  merged <- merge(
    r$ho_edges[, c("from", "to", "probability")],
    py, by = c("from", "to"), suffixes = c("_r", "_py"))
  expect_equal(nrow(merged), nrow(r$ho_edges),
    info = "Not all edges matched by from/to")
  expect_equal(merged$probability, merged$weight, tolerance = 1e-10)
})

# ===========================================================================
# Section 11: Coverage for previously uncovered paths
# ===========================================================================

# --- .hon_parse_input: all-NA row returns empty trajectory and is dropped ---
test_that(".hon_parse_input drops all-NA rows", {
  df <- data.frame(T1 = c(NA, "A"), T2 = c(NA, "B"), stringsAsFactors = FALSE)
  result <- .hon_parse_input(df)
  # The first row is all-NA: it should be dropped, leaving 1 trajectory
  expect_equal(length(result), 1L)
  expect_equal(result[[1L]], c("A", "B"))
})

# --- .hon_parse_input: invalid input type stops ---
test_that(".hon_parse_input stops on invalid input", {
  expect_error(.hon_parse_input(42), "data.frame or list")
})

# --- .hon_parse_input: collapse_repeats with length-1 trajectory ---
test_that(".hon_parse_input collapse_repeats returns length-1 traj unchanged", {
  # A length-1 traj is filtered out (needs >=2 for transitions), so test
  # that a 2-element traj after collapsing still works
  trajs <- list(c("A", "A", "B"))
  result <- .hon_parse_input(trajs, collapse_repeats = TRUE)
  expect_equal(result[[1L]], c("A", "B"))
})

# --- .hon_parse_input: collapse_repeats with single-element input ---
test_that(".hon_parse_input collapse_repeats skips length-1 element", {
  # After collapsing, length-1 elements stay as-is (branch: length(traj)<=1)
  # We pass a list with a 2-element traj where all are same => collapses to 1
  # => that traj gets filtered (< 2 states), so result is empty
  trajs <- list(c("A", "A"))
  result <- .hon_parse_input(trajs, collapse_repeats = TRUE)
  expect_equal(length(result), 0L)
})

# --- .hon_kld: empty distribution a returns 0 ---
test_that(".hon_kld returns 0 for empty distribution a", {
  a <- numeric(0L)
  b <- c(X = 0.5, Y = 0.5)
  expect_equal(.hon_kld(a, b), 0.0)
})

# --- .hon_kld: p_b = 0 for a target with p_a > 0 returns Inf ---
test_that(".hon_kld returns Inf when p_b = 0 for supported target", {
  a <- c(X = 1.0)
  b <- c(Y = 1.0)  # X not in b => p_b(X) = 0
  expect_equal(.hon_kld(a, b), Inf)
})

# --- .hon_get_extensions: order key missing returns character(0) ---
test_that(".hon_get_extensions returns character(0) when order key missing", {
  cache <- new.env(hash = TRUE, parent = emptyenv())
  # Put an entry in cache but without the queried order
  inner <- new.env(hash = TRUE, parent = emptyenv())
  inner[["3"]] <- c("some_key")
  cache[["key_a"]] <- inner
  result <- .honp_extend_source_fast  # just checking .hon_get_extensions indirectly
  ext <- .hon_get_extensions("key_a", 99L, cache)
  expect_equal(ext, character(0L))
})

# --- .hon_sequence_to_node: empty sequence returns "" ---
test_that(".hon_sequence_to_node returns empty string for empty input", {
  expect_equal(.hon_sequence_to_node(character(0L)), "")
})

# --- .hon_graph_to_edgelist: empty graph returns empty data.frame ---
test_that(".hon_graph_to_edgelist returns empty data.frame for empty graph", {
  graph <- new.env(hash = TRUE, parent = emptyenv())
  result <- .hon_graph_to_edgelist(graph)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0L)
  expect_true(all(c("path", "from", "to", "count", "probability",
                     "from_order", "to_order") %in% names(result)))
})

# --- .hon_assemble_output: empty rules gives observed_max_order = 1 ---
test_that(".hon_assemble_output uses max_order_observed=1 when no rules", {
  graph <- new.env(hash = TRUE, parent = emptyenv())
  rules <- new.env(hash = TRUE, parent = emptyenv())
  trajs <- list(c("A", "B"))
  result <- .hon_assemble_output(graph, rules, trajs, 3L, 1L)
  expect_equal(result$max_order_observed, 1L)
})

# --- build_hon: all-NA trajectories error ---
test_that("build_hon stops when all trajectories are too short after parsing", {
  # All rows are single-state (no transitions possible)
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(build_hon(df, max_order = 2L, min_freq = 1L),
               "No valid trajectories")
})


# ===========================================================================
# Section 12: pathways() tests for HON
# ===========================================================================

test_that("pathways.net_hon returns character vector of paths", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  hon <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")
  pw <- pathways(hon)
  expect_true(is.character(pw))
})

test_that("pathways.net_hon returns empty for first-order-only HON", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  hon <- build_hon(trajs, max_order = 1L, min_freq = 1L)
  pw <- pathways(hon)
  expect_equal(pw, character(0))
})

test_that("pathways.net_hon min_count filters rare paths", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  hon <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")
  pw_low  <- pathways(hon, min_count = 1L)
  pw_high <- pathways(hon, min_count = 999L)
  expect_true(length(pw_low) >= length(pw_high))
  expect_equal(length(pw_high), 0L)
})

test_that("pathways.net_hon top parameter limits results", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  hon <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")
  pw_all <- pathways(hon)
  pw_top <- pathways(hon, top = 1L)
  expect_true(length(pw_all) >= length(pw_top))
  if (length(pw_all) > 0L) expect_true(length(pw_top) <= 1L)
})

test_that("pathways.net_hon min_prob filters low-probability paths", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  hon <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")
  pw_all  <- pathways(hon)
  pw_filt <- pathways(hon, min_prob = 0.99)
  expect_true(length(pw_all) >= length(pw_filt))
})

test_that("pathways.net_hon order parameter selects specific order", {
  trajs <- list(
    c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"), c("X", "A", "C"),
    c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D"), c("Y", "A", "D")
  )
  hon <- build_hon(trajs, max_order = 3L, min_freq = 1L, method = "hon+")
  pw2 <- pathways(hon, order = 2L)
  expect_true(is.character(pw2))
})