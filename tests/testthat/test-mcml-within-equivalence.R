# Verify the within-cluster networks produced by build_mcml(type = "tna") are
# bit-identical to the conditional Markov chain obtained by:
#   1. masking out non-focal-cluster states with NA in the original sequences
#   2. counting consecutive (from, to) pairs where both endpoints are non-NA
#   3. row-normalising
# i.e. P(next = j | current = i, both in cluster k).
#
# Tested across actor / time / session-split scenarios and cross-validated
# against tna::tna() on the masked sequences.

mask_to_cluster <- function(wide_df, cluster_nodes) {
  as.data.frame(lapply(wide_df, function(col) {
    v <- as.character(col); v[!v %in% cluster_nodes] <- NA_character_; v
  }), stringsAsFactors = FALSE)
}

ref_within_tna <- function(wide_df, cluster_nodes) {
  pairs <- do.call(rbind, lapply(seq_len(ncol(wide_df) - 1), function(t) {
    data.frame(from = as.character(wide_df[[t]]),
               to   = as.character(wide_df[[t + 1]]),
               stringsAsFactors = FALSE)
  }))
  pairs <- pairs[!is.na(pairs$from) & !is.na(pairs$to), ]
  pairs <- pairs[pairs$from %in% cluster_nodes & pairs$to %in% cluster_nodes, ]
  fc <- factor(pairs$from, levels = cluster_nodes)
  tc <- factor(pairs$to,   levels = cluster_nodes)
  cnt <- table(fc, tc)
  m <- matrix(as.numeric(cnt), length(cluster_nodes), length(cluster_nodes),
              dimnames = list(cluster_nodes, cluster_nodes))
  rs <- rowSums(m)
  m / ifelse(rs == 0, 1, rs)
}

# Synthetic dataset shared across scenarios A-D
make_fixture <- function() {
  set.seed(1)
  states <- c("a1","a2","a3","b1","b2","c1","c2","c3")
  clu <- list(A = c("a1","a2","a3"),
              B = c("b1","b2"),
              C = c("c1","c2","c3"))
  n_actors <- 80; n_steps <- 25
  wide_mat <- matrix(sample(states, n_actors * n_steps, TRUE),
                     n_actors, n_steps)
  wide_df <- as.data.frame(wide_mat, stringsAsFactors = FALSE)
  names(wide_df) <- paste0("T", seq_len(n_steps))
  list(wide = wide_df, mat = wide_mat, clu = clu,
       n_actors = n_actors, n_steps = n_steps)
}

test_that("within-cluster network = NA-masked tna on wide-sequence input", {
  fx <- make_fixture()
  mc <- build_mcml(fx$wide, clusters = fx$clu, type = "tna")
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(fx$wide, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster network = NA-masked tna on long-format actor input", {
  fx <- make_fixture()
  long_df <- data.frame(
    actor  = rep(paste0("u", seq_len(fx$n_actors)), each = fx$n_steps),
    action = as.character(t(fx$mat)),
    stringsAsFactors = FALSE
  )
  mc <- suppressMessages(build_mcml(long_df, clusters = fx$clu,
                                    actor = "actor", action = "action",
                                    type = "tna"))
  prep <- suppressMessages(prepare(long_df, actor = "actor", action = "action"))
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(prep$sequence_data, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster network = NA-masked tna with timestamps (no split)", {
  fx <- make_fixture()
  ts_seq <- as.POSIXct("2024-01-01") +
    (seq_len(fx$n_actors * fx$n_steps) - 1) * 30
  long_ts <- data.frame(
    actor  = rep(paste0("u", seq_len(fx$n_actors)), each = fx$n_steps),
    action = as.character(t(fx$mat)),
    ts     = ts_seq,
    stringsAsFactors = FALSE
  )
  mc <- suppressMessages(build_mcml(long_ts, clusters = fx$clu,
                                    actor = "actor", action = "action",
                                    time = "ts", time_threshold = 1e6,
                                    type = "tna"))
  prep <- suppressMessages(prepare(long_ts, actor = "actor", action = "action",
                                   time = "ts", time_threshold = 1e6))
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(prep$sequence_data, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster network = NA-masked tna with session splits", {
  fx <- make_fixture()
  ts_seq <- as.POSIXct("2024-01-01") +
    (seq_len(fx$n_actors * fx$n_steps) - 1) * 30
  long_split <- data.frame(
    actor  = rep(paste0("u", seq_len(fx$n_actors)), each = fx$n_steps),
    action = as.character(t(fx$mat)),
    ts     = ts_seq,
    stringsAsFactors = FALSE
  )
  idx_gap <- which(rep(seq_len(fx$n_steps), fx$n_actors) > 12)
  long_split$ts[idx_gap] <- long_split$ts[idx_gap] + 60 * 60   # +1h gap
  mc <- suppressMessages(build_mcml(long_split, clusters = fx$clu,
                                    actor = "actor", action = "action",
                                    time = "ts", type = "tna"))
  prep <- suppressMessages(prepare(long_split, actor = "actor",
                                   action = "action", time = "ts"))
  # The reference now uses prep$sequence_data, which has split sessions —
  # equivalence still holds because both routes drop the same cross-session
  # transitions.
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(prep$sequence_data, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster network = NA-masked tna with explicit session column (no time)", {
  fx <- make_fixture()
  F_df <- data.frame(
    actor   = rep(paste0("u", seq_len(fx$n_actors)), each = fx$n_steps),
    action  = as.character(t(fx$mat)),
    session = rep(c(rep("morning",   floor(fx$n_steps / 2)),
                    rep("afternoon", ceiling(fx$n_steps / 2))),
                  fx$n_actors),
    stringsAsFactors = FALSE
  )
  mc <- suppressMessages(build_mcml(F_df, clusters = fx$clu,
                                    actor = "actor", action = "action",
                                    session = "session", type = "tna"))
  prep <- suppressMessages(prepare(F_df, actor = "actor", action = "action",
                                   session = "session"))
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(prep$sequence_data, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster network = NA-masked tna with actor + time + session", {
  fx <- make_fixture()
  ts_seq <- as.POSIXct("2024-01-01") +
    (seq_len(fx$n_actors * fx$n_steps) - 1) * 30
  G_df <- data.frame(
    actor   = rep(paste0("u", seq_len(fx$n_actors)), each = fx$n_steps),
    action  = as.character(t(fx$mat)),
    ts      = ts_seq,
    session = rep(c(rep("morning",   floor(fx$n_steps / 2)),
                    rep("afternoon", ceiling(fx$n_steps / 2))),
                  fx$n_actors),
    stringsAsFactors = FALSE
  )
  # Inject a +1h gap inside the morning session to force an extra split
  gap_idx <- which(rep(seq_len(fx$n_steps), fx$n_actors) %in%
                     seq(floor(fx$n_steps / 2) - 4, floor(fx$n_steps / 2)))
  G_df$ts[gap_idx] <- G_df$ts[gap_idx] + 60 * 60

  mc <- suppressMessages(build_mcml(G_df, clusters = fx$clu,
                                    actor = "actor", action = "action",
                                    time = "ts", session = "session",
                                    type = "tna"))
  prep <- suppressMessages(prepare(G_df, actor = "actor", action = "action",
                                   time = "ts", session = "session"))
  expect_gt(nrow(prep$sequence_data), fx$n_actors)   # session-splitting active
  for (cl_name in names(fx$clu)) {
    ref <- ref_within_tna(prep$sequence_data, fx$clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("within-cluster equivalence on bundled ai_long (actor + time + session)", {
  skip_if_not(exists("ai_long"))
  clu <- list(work = c("Plan", "Execute", "Investigate", "Repair"),
              comm = c("Delegate", "Ask", "Explain", "Report"))
  mc <- suppressMessages(build_mcml(ai_long, clusters = clu,
                                    actor = "project", action = "code",
                                    time = "timestamp", session = "session_id",
                                    type = "tna"))
  prep <- suppressMessages(prepare(ai_long, actor = "project", action = "code",
                                   time = "timestamp", session = "session_id"))
  for (cl_name in names(clu)) {
    ref <- ref_within_tna(prep$sequence_data, clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("ai_long cluster", cl_name))
  }
})

test_that("within-cluster equivalence on group_regulation_long incl. singleton cluster", {
  skip_if_not(exists("group_regulation_long"))
  clu <- list(cog = c("plan", "monitor", "adapt", "synthesis"),
              soc = c("cohesion", "consensus", "discuss", "coregulate"),
              emo = c("emotion"))
  mc <- suppressMessages(build_mcml(group_regulation_long, clusters = clu,
                                    actor = "Actor", action = "Action",
                                    time = "Time", session = "Course",
                                    type = "tna"))
  prep <- suppressMessages(prepare(group_regulation_long, actor = "Actor",
                                   action = "Action", time = "Time",
                                   session = "Course"))
  for (cl_name in names(clu)) {
    ref <- ref_within_tna(prep$sequence_data, clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("group_regulation_long cluster", cl_name))
  }
})

test_that("within-cluster equivalence on trajectories matrix incl. singleton cluster", {
  skip_if_not(exists("trajectories"))
  traj_df <- as.data.frame(trajectories, stringsAsFactors = FALSE)
  clu <- list(present = c("Active", "Average"),
              absent  = c("Disengaged"))
  mc <- build_mcml(traj_df, clusters = clu, type = "tna")
  for (cl_name in names(clu)) {
    ref <- ref_within_tna(traj_df, clu[[cl_name]])
    got <- mc$clusters[[cl_name]]$weights
    got <- got[rownames(ref), colnames(ref), drop = FALSE]
    expect_equal(unname(got), unname(ref), tolerance = 0,
                 info = paste("trajectories cluster", cl_name))
  }
})

test_that("macro and within inits match tna::tna() (column-1 first state)", {
  skip_if_not_installed("tna")
  fx <- make_fixture()
  mc_tna <- as_tna(build_mcml(fx$wide, clusters = fx$clu, type = "tna"))

  # Macro inits: tna::tna() on cluster-recoded sequence
  node2clu <- unlist(lapply(names(fx$clu), function(k) {
    setNames(rep(k, length(fx$clu[[k]])), fx$clu[[k]])
  }))
  recoded <- as.data.frame(lapply(fx$wide, function(c) node2clu[c]),
                           stringsAsFactors = FALSE)
  tna_macro <- suppressWarnings(tna::tna(recoded))
  expect_equal(unname(mc_tna$macro$inits[names(tna_macro$inits)]),
               unname(tna_macro$inits), tolerance = 0)

  # Within inits per cluster: tna::tna() on NA-masked sequence
  for (cl_name in names(fx$clu)) {
    masked <- mask_to_cluster(fx$wide, fx$clu[[cl_name]])
    tna_net <- suppressWarnings(tna::tna(masked))
    common <- intersect(names(tna_net$inits), names(mc_tna[[cl_name]]$inits))
    expect_equal(unname(mc_tna[[cl_name]]$inits[common]),
                 unname(tna_net$inits[common]), tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})

test_that("build_mcml on attention netobject preserves model-derived weights", {
  skip_if_not_installed("tna")
  set.seed(2)
  states <- c("a","b","c","d","e","f")
  wide <- as.data.frame(matrix(sample(states, 50 * 20, TRUE), 50, 20),
                        stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(20))
  clu <- list(left = c("a","b","c"), right = c("d","e","f"))

  # Attention method: weights are NOT raw transition counts.
  net_attn <- build_network(wide, method = "attention")
  expect_false(net_attn$method %in% c("relative", "frequency", "co_occurrence"))

  mc_attn <- build_mcml(net_attn, clusters = clu)
  # Macro must come from cluster-aggregating $weights (matrix path), not from
  # re-counting $data. Reference: cluster_summary() on the same matrix.
  ref <- cluster_summary(net_attn$weights, clusters = clu, method = "sum")
  expect_equal(unname(mc_attn$macro$weights),
               unname(ref$macro$weights), tolerance = 0)
})

test_that("within-cluster network matches tna::tna() on NA-masked sequences", {
  skip_if_not_installed("tna")
  fx <- make_fixture()
  mc <- build_mcml(fx$wide, clusters = fx$clu, type = "tna")
  for (cl_name in names(fx$clu)) {
    masked <- mask_to_cluster(fx$wide, fx$clu[[cl_name]])
    tna_net <- suppressWarnings(tna::tna(masked))
    tna_w <- tna_net$weights
    got <- mc$clusters[[cl_name]]$weights
    common <- intersect(rownames(tna_w), rownames(got))
    expect_gt(length(common), 0L)
    expect_equal(unname(tna_w[common, common]),
                 unname(got[common, common]),
                 tolerance = 0,
                 info = paste("cluster", cl_name))
  }
})


# ---------------------------------------------------------------------------
# Comprehensive cluster-by-cluster verification (mirrors the verification Rmd
# at Tutorial_docs/mcml_within_equivalence_verification.Rmd, but on a fully
# synthetic dataset so it runs in CI without external files).
#
# The fixture exercises ALL FOUR long-format arguments at once:
# actor + action + time + session, plus session-splitting via a time gap.
# Coverage:
#   * 4 multi-state clusters of varying size + 1 singleton cluster
#   * Within-cluster transition matrix vs NA-mask reference
#   * Within-cluster transition matrix vs tna::tna(masked) [if installed]
#   * Macro inits vs first-state distribution of the cluster-recoded sequence
#   * Within inits vs tna::tna(masked)$inits [if installed]
# ---------------------------------------------------------------------------

.make_full_fixture <- function() {
  set.seed(99)
  states <- c("a1","a2","a3", "b1","b2","b3","b4",
              "c1","c2", "d1","d2","d3", "z")    # z = singleton
  clu <- list(A = c("a1","a2","a3"),
              B = c("b1","b2","b3","b4"),
              C = c("c1","c2"),
              D = c("d1","d2","d3"),
              Z = c("z"))
  n_actors <- 90; n_steps <- 24
  wide_mat <- matrix(sample(states, n_actors * n_steps, TRUE),
                     n_actors, n_steps)
  ts <- as.POSIXct("2024-01-01") +
    (seq_len(n_actors * n_steps) - 1) * 30
  long_df <- data.frame(
    actor   = rep(paste0("u", seq_len(n_actors)), each = n_steps),
    action  = as.character(t(wide_mat)),
    ts      = ts,
    session = rep(c(rep("morning",   floor(n_steps / 2)),
                    rep("afternoon", ceiling(n_steps / 2))), n_actors),
    stringsAsFactors = FALSE
  )
  # Inject a 1-h gap right before the morning/afternoon boundary in some actors
  gap_idx <- which(rep(seq_len(n_steps), n_actors) %in%
                     seq(floor(n_steps / 2) - 3, floor(n_steps / 2)))
  long_df$ts[gap_idx] <- long_df$ts[gap_idx] + 60 * 60

  list(states = states, clu = clu, long = long_df,
       n_actors = n_actors, n_steps = n_steps)
}

test_that("comprehensive verification: weights & inits across all clusters (synthetic)", {
  fx <- .make_full_fixture()

  mc <- suppressMessages(build_mcml(
    fx$long, clusters = fx$clu,
    actor = "actor", action = "action", time = "ts", session = "session",
    type = "tna"
  ))
  mc_tna <- as_tna(mc)
  prep   <- suppressMessages(prepare(
    fx$long, actor = "actor", action = "action",
    time = "ts", session = "session"
  ))
  sequences <- prep$sequence_data

  # Macro inits: first state of column 1 of the cluster-recoded sequence
  node2clu <- unlist(lapply(names(fx$clu), function(k) {
    setNames(rep(k, length(fx$clu[[k]])), fx$clu[[k]])
  }))
  first_clu <- node2clu[as.character(sequences[[1]])]
  first_clu <- first_clu[!is.na(first_clu)]
  macro_ref <- as.numeric(table(factor(first_clu, levels = names(fx$clu))))
  macro_ref <- macro_ref / sum(macro_ref)
  expect_equal(unname(mc_tna$macro$inits[names(fx$clu)]),
               macro_ref, tolerance = 0)

  # Per-cluster: weights and inits
  for (cl_name in names(fx$clu)) {
    cl_nodes <- fx$clu[[cl_name]]
    masked   <- mask_to_cluster(sequences, cl_nodes)

    # Weights: NA-mask reference
    ref_W <- ref_within_tna(sequences, cl_nodes)
    got_W <- mc$clusters[[cl_name]]$weights
    got_W <- got_W[rownames(ref_W), colnames(ref_W), drop = FALSE]
    expect_equal(unname(got_W), unname(ref_W), tolerance = 0,
                 info = paste("weights:", cl_name))

    # Inits: column-1 first-state distribution within the cluster
    first_in_cl <- as.character(sequences[[1]])
    first_in_cl[!first_in_cl %in% cl_nodes] <- NA_character_
    first_in_cl <- first_in_cl[!is.na(first_in_cl)]
    if (length(first_in_cl) > 0L) {
      ref_I <- as.numeric(table(factor(first_in_cl, levels = cl_nodes)))
      ref_I <- ref_I / sum(ref_I)
    } else {
      ref_I <- rep(0, length(cl_nodes))
    }
    names(ref_I) <- cl_nodes
    got_I <- mc_tna[[cl_name]]$inits
    expect_equal(unname(got_I[cl_nodes]), unname(ref_I), tolerance = 0,
                 info = paste("inits:", cl_name))
  }
})

test_that("comprehensive verification: cross-validate every cluster vs tna::tna() (synthetic)", {
  skip_if_not_installed("tna")
  fx <- .make_full_fixture()

  mc <- suppressMessages(build_mcml(
    fx$long, clusters = fx$clu,
    actor = "actor", action = "action", time = "ts", session = "session",
    type = "tna"
  ))
  mc_tna <- as_tna(mc)
  prep   <- suppressMessages(prepare(
    fx$long, actor = "actor", action = "action",
    time = "ts", session = "session"
  ))
  sequences <- prep$sequence_data

  for (cl_name in names(fx$clu)) {
    cl_nodes <- fx$clu[[cl_name]]
    masked   <- mask_to_cluster(sequences, cl_nodes)
    tna_net  <- suppressWarnings(tna::tna(masked))

    # Weights
    common <- intersect(rownames(tna_net$weights),
                        rownames(mc$clusters[[cl_name]]$weights))
    expect_gt(length(common), 0L)
    expect_equal(
      unname(mc$clusters[[cl_name]]$weights[common, common]),
      unname(tna_net$weights[common, common]),
      tolerance = 0, info = paste("weights vs tna::", cl_name)
    )

    # Inits
    common_i <- intersect(names(tna_net$inits),
                          names(mc_tna[[cl_name]]$inits))
    expect_equal(
      unname(mc_tna[[cl_name]]$inits[common_i]),
      unname(tna_net$inits[common_i]),
      tolerance = 0, info = paste("inits vs tna::", cl_name)
    )
  }
})


# ---------------------------------------------------------------------------
# Structural invariants for mcml objects: every construction path must
# produce the same field schema, and within a path, label/rowname/init-name
# alignment must hold across macro and per-cluster slots.
# ---------------------------------------------------------------------------

test_that("mcml object has uniform schema across all construction paths", {
  set.seed(7)
  states <- c("a","b","c","d","e","f")
  clu    <- list(L = c("a","b","c"), R = c("d","e","f"))
  wide   <- as.data.frame(matrix(sample(states, 50 * 15, TRUE), 50, 15),
                          stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(15))
  edges  <- data.frame(from = sample(states, 200, TRUE),
                       to   = sample(states, 200, TRUE),
                       weight = sample(1:5, 200, TRUE),
                       stringsAsFactors = FALSE)
  W <- matrix(runif(36), 6, 6); rownames(W) <- colnames(W) <- states

  paths <- list(
    seq  = build_mcml(wide,  clusters = clu, type = "tna"),
    edge = build_mcml(edges, clusters = clu, type = "tna"),
    mat  = cluster_summary(W, clusters = clu, method = "sum")
  )

  expected_top <- c("macro", "clusters", "cluster_members",
                    "edges", "meta")
  expected_macro <- c("weights", "inits", "labels", "data")
  expected_clu   <- c("weights", "inits", "labels", "data")
  expected_meta  <- c("type", "method", "directed", "n_nodes",
                      "n_clusters", "cluster_sizes", "source")

  for (p in names(paths)) {
    o <- paths[[p]]
    expect_setequal(names(o), expected_top)
    expect_setequal(names(o$macro), expected_macro)
    for (cl in names(o$clusters)) {
      expect_setequal(names(o$clusters[[cl]]), expected_clu)
    }
    # meta is a strict superset (some paths add extra fields)
    expect_true(all(expected_meta %in% names(o$meta)),
                info = paste("meta missing fields on path", p))
  }
})

test_that("mcml object label / rowname / init-name alignment", {
  set.seed(7)
  states <- c("a","b","c","d","e","f")
  clu    <- list(L = c("a","b","c"), R = c("d","e","f"))
  wide   <- as.data.frame(matrix(sample(states, 50 * 15, TRUE), 50, 15),
                          stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(15))
  W <- matrix(runif(36), 6, 6); rownames(W) <- colnames(W) <- states

  for (mc in list(build_mcml(wide, clusters = clu, type = "tna"),
                  cluster_summary(W, clusters = clu, method = "sum"))) {
    # Macro
    expect_identical(mc$macro$labels, rownames(mc$macro$weights))
    expect_identical(mc$macro$labels, colnames(mc$macro$weights))
    expect_setequal(names(mc$macro$inits), mc$macro$labels)
    expect_equal(sum(mc$macro$inits), 1, tolerance = 1e-12)

    # Per-cluster
    for (cl in names(mc$clusters)) {
      n <- mc$clusters[[cl]]
      expect_identical(n$labels, rownames(n$weights))
      expect_identical(n$labels, colnames(n$weights))
      expect_setequal(names(n$inits), n$labels)
    }

    # cluster_members keys equal macro labels
    expect_setequal(names(mc$cluster_members), mc$macro$labels)
  }
})

test_that("as_tna(mcml) returns netobject_group with inits propagated", {
  set.seed(7)
  states <- c("a","b","c","d","e","f")
  clu    <- list(L = c("a","b","c"), R = c("d","e","f"))
  wide   <- as.data.frame(matrix(sample(states, 50 * 15, TRUE), 50, 15),
                          stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(15))

  res <- as_tna(build_mcml(wide, clusters = clu, type = "tna"))
  expect_s3_class(res, "netobject_group")

  # Macro netobject
  expect_true(inherits(res$macro, "netobject"))
  expect_false(is.null(res$macro$inits))
  expect_equal(sum(res$macro$inits), 1, tolerance = 1e-12)
  expect_setequal(names(res$macro$inits), rownames(res$macro$weights))

  # Each cluster netobject
  for (cl in c("L", "R")) {
    expect_true(inherits(res[[cl]], "netobject"))
    expect_false(is.null(res[[cl]]$inits))
    expect_setequal(names(res[[cl]]$inits), rownames(res[[cl]]$weights))
  }
})
