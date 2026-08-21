# ---- Tests for transition_entropy() ----

.P3 <- matrix(
  c(0.7, 0.2, 0.1,
    0.3, 0.5, 0.2,
    0.2, 0.3, 0.5),
  nrow = 3, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

# ---- structure ----

test_that("transition_entropy returns net_transition_entropy with all fields", {
  te <- transition_entropy(.P3, base = 2)
  expect_s3_class(te, "net_transition_entropy")
  expect_named(te, c("row_entropy", "row_entropy_norm",
                     "stationary",
                     "stationary_entropy", "stationary_entropy_norm",
                     "entropy_rate", "entropy_rate_norm",
                     "redundancy", "redundancy_norm",
                     "max_entropy", "base", "states"))
  expect_length(te$row_entropy, 3L)
  expect_length(te$stationary,  3L)
  expect_equal(te$states, c("A", "B", "C"))
  expect_equal(te$base, 2)
})

test_that("normalised fields equal raw / log_b(n) and lie in [0, 1]", {
  te <- transition_entropy(.P3, base = 2)
  H_max <- log2(3)
  expect_equal(te$max_entropy,             H_max,                       tolerance = 1e-12)
  expect_equal(te$entropy_rate_norm,       te$entropy_rate       / H_max, tolerance = 1e-12)
  expect_equal(te$stationary_entropy_norm, te$stationary_entropy / H_max, tolerance = 1e-12)
  expect_equal(unname(te$row_entropy_norm),
               unname(te$row_entropy) / H_max, tolerance = 1e-12)
  expect_true(te$entropy_rate_norm       >= 0 && te$entropy_rate_norm       <= 1 + 1e-10)
  expect_true(te$stationary_entropy_norm >= 0 && te$stationary_entropy_norm <= 1 + 1e-10)
  expect_true(all(te$row_entropy_norm >= 0 - 1e-10) &&
              all(te$row_entropy_norm <= 1 + 1e-10))
})

test_that("redundancy_norm is fraction of marginal uncertainty captured by memory", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$redundancy_norm,
               te$redundancy / te$stationary_entropy,
               tolerance = 1e-12)
  expect_true(te$redundancy_norm >= 0 - 1e-10 &&
              te$redundancy_norm <= 1 + 1e-10)
})

# ---- numerical identities ----

test_that("entropy_rate equals sum(pi * row_entropy)", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$entropy_rate,
               sum(te$stationary * te$row_entropy),
               tolerance = 1e-12)
})

test_that("redundancy equals stationary_entropy minus entropy_rate", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(te$redundancy,
               te$stationary_entropy - te$entropy_rate,
               tolerance = 1e-12)
})

test_that("entropy_rate <= stationary_entropy (information inequality)", {
  te <- transition_entropy(.P3, base = 2)
  expect_lte(te$entropy_rate, te$stationary_entropy + 1e-10)
  expect_gte(te$redundancy, -1e-10)
})

test_that("uniform rows give zero redundancy and h = log_b(n)", {
  P <- matrix(1 / 3, 3, 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  te <- transition_entropy(P, base = 2)
  expect_equal(te$redundancy, 0, tolerance = 1e-10)
  expect_equal(te$entropy_rate,       log2(3), tolerance = 1e-10)
  expect_equal(te$stationary_entropy, log2(3), tolerance = 1e-10)
})

test_that("deterministic permutation chain has zero entropy rate", {
  P <- matrix(c(0, 1, 0,
                0, 0, 1,
                1, 0, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  te <- transition_entropy(P, base = 2)
  expect_equal(te$entropy_rate, 0, tolerance = 1e-12)
  expect_equal(te$row_entropy, c(A = 0, B = 0, C = 0), tolerance = 1e-12)
})

test_that("base 2 vs base e differ by ln(2) factor", {
  te2 <- transition_entropy(.P3, base = 2)
  tee <- transition_entropy(.P3, base = exp(1))
  expect_equal(te2$entropy_rate * log(2), tee$entropy_rate,
               tolerance = 1e-12)
  expect_equal(te2$row_entropy * log(2), tee$row_entropy,
               tolerance = 1e-12)
})

test_that("row entropies bounded by log_b(n_nonzero_columns)", {
  te <- transition_entropy(.P3, base = 2)
  expect_true(all(te$row_entropy <= log2(3) + 1e-10))
  expect_true(all(te$row_entropy >= 0 - 1e-10))
})

# ---- stationary distribution checks ----

test_that("stationary distribution sums to 1 and satisfies pi P = pi", {
  te <- transition_entropy(.P3, base = 2)
  expect_equal(sum(te$stationary), 1, tolerance = 1e-10)
  expect_equal(as.vector(te$stationary %*% .P3),
               unname(te$stationary), tolerance = 1e-8)
})

# ---- input dispatch ----

test_that("transition_entropy accepts a netobject", {
  data(trajectories, package = "Nestimate")
  net <- build_network(as.data.frame(trajectories), method = "relative")
  te  <- transition_entropy(net, base = 2)
  expect_s3_class(te, "net_transition_entropy")
  expect_length(te$states, ncol(net$weights))
})

test_that("transition_entropy accepts a wide data.frame", {
  data(trajectories, package = "Nestimate")
  te <- transition_entropy(as.data.frame(trajectories), base = 2)
  expect_s3_class(te, "net_transition_entropy")
})

test_that("transition_entropy normalises rows that don't sum to 1", {
  P <- .P3 * 2
  expect_warning(te <- transition_entropy(P, base = 2),
                 regexp = "normalizing")
  te_ref <- transition_entropy(.P3, base = 2)
  expect_equal(te$entropy_rate, te_ref$entropy_rate, tolerance = 1e-10)
})

test_that("base argument is validated", {
  expect_error(transition_entropy(.P3, base = 1),
               regexp = "base != 1|stopifnot")
  expect_error(transition_entropy(.P3, base = -2))
  expect_error(transition_entropy(.P3, base = c(2, 10)))
})

# ---- group dispatch ----

test_that("transition_entropy dispatches on netobject_group", {
  data(trajectories, package = "Nestimate")
  df <- as.data.frame(trajectories)
  df$grp <- rep(c("g1", "g2"), length.out = nrow(df))
  grp <- build_network(df, method = "relative", group = "grp")
  te_grp <- transition_entropy(grp, base = 2)
  expect_s3_class(te_grp, "net_transition_entropy_group")
  expect_length(te_grp, 2L)
  expect_true(all(vapply(te_grp, inherits, logical(1),
                         "net_transition_entropy")))
})

# ---- S3 methods dispatch ----

test_that("print/summary/plot dispatch correctly", {
  te <- transition_entropy(.P3, base = 2)
  expect_output(print(te), "Transition Entropy")
  s <- summary(te)
  expect_s3_class(s, "summary.net_transition_entropy")
  expect_true(is.data.frame(s$table))
  expect_true(is.data.frame(s$chain))
  ## Tidy table has 6 columns including normalised + percent
  expect_named(s$table, c("state", "stationary", "row_entropy",
                          "row_entropy_norm",
                          "contribution", "contribution_pct"))
  ## Sorted by contribution_pct descending
  expect_equal(s$table$contribution_pct,
               sort(s$table$contribution_pct, decreasing = TRUE))
  ## Contributions sum to entropy_rate; pct sums to 100
  expect_equal(sum(s$table$contribution), te$entropy_rate, tolerance = 1e-10)
  expect_equal(sum(s$table$contribution_pct), 100, tolerance = 1e-10)
  ## Chain table has 4 rows: h(P), H(pi), redundancy, ceiling
  expect_equal(nrow(s$chain), 4L)
  g <- plot(te)
  expect_true(inherits(g, "gg"))
})

# ---- entropy_network() ----

test_that("entropy_network returns a netobject whose weights sum to h(P)", {
  ent <- entropy_network(.P3, base = 2)
  te  <- transition_entropy(.P3, base = 2)
  expect_s3_class(ent, "netobject")
  expect_s3_class(ent, "cograph_network")
  expect_equal(ent$method, "entropy")
  expect_true(ent$directed)
  expect_equal(dimnames(ent$weights), list(c("A", "B", "C"), c("A", "B", "C")))
  ## Exact edge-level decomposition of the entropy rate
  expect_equal(sum(ent$weights), te$entropy_rate, tolerance = 1e-12)
  expect_equal(ent$params$entropy_rate, te$entropy_rate, tolerance = 1e-12)
  ## Row sums of contributions equal pi_i * H(row_i)
  expect_equal(rowSums(ent$weights),
               te$stationary * te$row_entropy, tolerance = 1e-12)
  ## Stationary distribution carried in params
  expect_equal(ent$params$stationary, te$stationary, tolerance = 1e-12)
})

test_that("entropy_network surprisal weights are -log2(P) on the support", {
  ent <- entropy_network(.P3, weight = "surprisal", base = 2)
  expect_equal(ent$weights, -log2(.P3), tolerance = 1e-12,
               ignore_attr = FALSE)
  expect_equal(ent$params$weight, "surprisal")
})

test_that("entropy_network zeroes impossible and deterministic transitions", {
  P <- matrix(
    c(0,   1,   0,
      0.5, 0,   0.5,
      1,   0,   0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  ent <- entropy_network(P, base = 2)
  ## P == 0 and P == 1 both carry zero weight
  expect_equal(ent$weights[P == 0], rep(0, sum(P == 0)))
  expect_equal(ent$weights["A", "B"], 0)  # deterministic A -> B
  expect_equal(ent$weights["C", "A"], 0)  # deterministic C -> A
  ## Only B's genuine 50/50 branch contributes
  expect_true(all(ent$weights[2, c(1, 3)] > 0))
  expect_equal(sum(ent$weights),
               transition_entropy(P, base = 2)$entropy_rate,
               tolerance = 1e-12)
})

test_that("entropy_network accepts a fitted netobject and preserves states", {
  seqs <- data.frame(
    T1 = c("plan", "do", "check", "plan"),
    T2 = c("do", "check", "plan", "do"),
    T3 = c("check", "plan", "do", "check"),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  ent <- entropy_network(net)
  expect_setequal(ent$nodes$name, net$nodes$name)
  expect_equal(sum(ent$weights),
               transition_entropy(net)$entropy_rate, tolerance = 1e-12)
})

test_that("entropy_network validates arguments", {
  expect_error(entropy_network(.P3, base = 1), "base")
  expect_error(entropy_network(.P3, weight = "nope"))
})

test_that("entropy_network inherits styling metadata from the source network", {
  seqs <- data.frame(
    T1 = c("plan", "do", "check", "plan"),
    T2 = c("do", "check", "plan", "do"),
    T3 = c("check", "plan", "do", "check"),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  ent <- entropy_network(net)
  ## Clone of the source: same nodes, inits, meta, groups; new weights/method
  expect_identical(ent$nodes, net$nodes)
  expect_identical(ent$inits, net$inits)
  ## meta inherited except for the added $splot styling contract
  expect_identical(ent$meta[setdiff(names(ent$meta), "splot")], net$meta)
  expect_identical(ent$node_groups, net$node_groups)
  expect_equal(ent$method, "entropy")
  expect_false(identical(ent$weights, net$weights))
  ## Edge table rebuilt from the entropy weights
  expect_equal(nrow(ent$edges), sum(ent$weights != 0))
})

test_that("entropy_network declares the entropy house style via meta$splot", {
  ent <- entropy_network(.P3, base = 2)
  spec <- ent$meta$splot
  expect_true(is.list(spec))
  expect_equal(spec$defaults$minimum, 0)
  expect_equal(spec$defaults$edge_label_digits, 2)
  expect_equal(spec$defaults$edge_color, "#D55E00")
  ## Node rings carry the stationary distribution, in node order
  expect_equal(spec$defaults$donut_fill,
               unname(ent$params$stationary), tolerance = 1e-12)
  expect_false(spec$defaults$donut_empty)
})

# ---- entropy_network(): production variant ----

test_that("production weights decompose the entropy production rate", {
  ent <- entropy_network(.P3, weight = "production", base = 2)
  W   <- ent$weights
  pi  <- ent$params$stationary
  Fm  <- pi * .P3
  ## Independent restricted-sum computation
  ok  <- Fm > 0 & t(Fm) > 0
  ref <- sum(Fm[ok] * log2(Fm[ok] / t(Fm)[ok]))
  expect_equal(sum(W), ref, tolerance = 1e-12)
  expect_equal(ent$params$production_rate, ref, tolerance = 1e-12)
  expect_true(ref >= 0)
  ## Pair sums are non-negative: (F_ij - F_ji) log(F_ij / F_ji) >= 0
  expect_true(all(W + t(W) >= -1e-12))
  expect_equal(ent$params$n_oneway_pairs, 0L)
  ## Signed weights: production keeps sign-based colouring (no flat edge_color)
  expect_null(ent$meta$splot$defaults$edge_color)
})

test_that("production excludes one-way pairs and counts them", {
  P <- matrix(
    c(0.5, 0.5, 0,
      0.2, 0.3, 0.5,
      0.5, 0,   0.5),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )
  ent <- entropy_network(P, weight = "production")
  ## A->C flow is zero but C->A positive, and B->C positive but C->B zero
  expect_equal(ent$params$n_oneway_pairs, 2L)
  expect_equal(ent$weights["C", "A"], 0)  # one-way: excluded
  expect_equal(ent$weights["B", "C"], 0)  # one-way: excluded
})

# ---- entropy_bayes() ----

test_that("entropy_bayes recovers plug-in entropy on dense counts", {
  counts <- round(.P3 * 20000)  # dense: posterior concentrates on truth
  eb <- entropy_bayes(counts, draws = 500, seed = 42)
  te <- transition_entropy(counts / rowSums(counts), base = 2)
  h  <- subset(eb$summary, quantity == "entropy_rate")
  expect_equal(h$mean, te$entropy_rate, tolerance = 0.01)
  expect_true(h$ci_lower < te$entropy_rate && te$entropy_rate < h$ci_upper)
  expect_s3_class(eb, "net_entropy_bayes")
})

test_that("entropy_bayes output structure is tidy and consistent", {
  counts <- round(.P3 * 500)
  eb <- entropy_bayes(counts, draws = 300, seed = 7)
  expect_named(eb, c("summary", "states", "edges", "network", "model",
                     "draws_entropy_rate", "prior", "draws", "ci",
                     "min_share", "base", "states_names"))
  expect_equal(nrow(eb$summary), 3L)
  expect_equal(nrow(eb$states), 3L)
  expect_equal(nrow(eb$edges), sum(counts > 0))
  expect_true(all(eb$edges$ci_lower <= eb$edges$contribution))
  expect_true(all(eb$edges$contribution <= eb$edges$ci_upper))
  expect_length(eb$draws_entropy_rate, 300L)
  ## Networks: posterior-mean full and credible-pruned model
  expect_s3_class(eb$network, "netobject")
  expect_s3_class(eb$model, "netobject")
  pruned <- eb$model$weights == 0 & eb$network$weights > 0
  flagged <- subset(eb$edges, !credible)
  expect_equal(sum(pruned), nrow(flagged))
})

test_that("entropy_bayes is reproducible under seed and flags sparse edges", {
  counts <- round(.P3 * 500)
  counts["A", "C"] <- 1  # razor-thin edge: must not be credible
  eb1 <- entropy_bayes(counts, draws = 300, seed = 11)
  eb2 <- entropy_bayes(counts, draws = 300, seed = 11)
  expect_identical(eb1$summary, eb2$summary)
  expect_identical(eb1$edges, eb2$edges)
  thin <- subset(eb1$edges, from == "A" & to == "C")
  expect_false(thin$credible)
})

test_that("entropy_bayes recovers counts from a relative netobject via $data", {
  seqs <- data.frame(
    T1 = c("plan", "do", "check", "plan", "do", "check"),
    T2 = c("do", "check", "plan", "do", "check", "plan"),
    T3 = c("check", "plan", "do", "check", "plan", "do"),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  eb  <- entropy_bayes(net, draws = 200, seed = 3)
  expect_s3_class(eb, "net_entropy_bayes")
  expect_setequal(eb$states_names, c("plan", "do", "check"))
})

test_that("entropy_bayes validates input", {
  expect_error(entropy_bayes(.P3), "COUNTS")           # probabilities rejected
  expect_error(entropy_bayes(round(.P3 * 100), prior = 0), "prior")
  expect_error(entropy_bayes(round(.P3 * 100), draws = 10), "draws")
})

test_that("entropy_bayes S3 methods dispatch", {
  counts <- round(.P3 * 500)
  eb <- entropy_bayes(counts, draws = 200, seed = 5)
  expect_output(print(eb), "Bayesian Transition Entropy")
  expect_output(summary(eb), "Per-edge posterior")
  g <- plot(eb)
  expect_true(inherits(g, "gg"))
})

# ---- entropy_trajectory() ----

.make_stream <- function(P, n, states = rownames(P), seed = 1) {
  set.seed(seed)
  s <- integer(n)
  s[1] <- 1L
  for (k in 2:n) s[k] <- sample.int(nrow(P), 1L, prob = P[s[k - 1L], ])
  data.frame(who = "a", t = seq_len(n), act = states[s],
             stringsAsFactors = FALSE)
}

test_that("trajectory of a stationary chain hovers at the entropy rate", {
  d  <- .make_stream(.P3, 6000, seed = 2)
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           window = 500, step = 100)
  te <- transition_entropy(.P3, base = 2)
  expect_s3_class(tr, "net_entropy_trajectory")
  expect_lt(abs(mean(tr$trajectory$entropy) - te$entropy_rate), 0.05)
  expect_lt(stats::sd(tr$trajectory$entropy), 0.1)
  ## Window bookkeeping: every window has exactly `window` transitions
  expect_true(all(tr$trajectory$n_transitions == 500L))
  expect_equal(nrow(tr$trajectory),
               length(seq.int(1L, (6000 - 1) - 500 + 1, by = 100)))
})

test_that("trajectory detects a regime switch", {
  states <- c("A", "B", "C")
  det <- states[rep(1:3, length.out = 2000)]                # deterministic cycle
  set.seed(4)
  rnd <- sample(states, 2000, replace = TRUE)               # iid uniform
  d <- data.frame(who = "a", t = seq_len(4000), act = c(det, rnd),
                  stringsAsFactors = FALSE)
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           window = 400, step = 100)$trajectory
  early <- subset(tr, time_end <= 2000)
  late  <- subset(tr, time_start > 2000)
  expect_lt(max(early$entropy), 0.01)                       # h = 0 regime
  expect_gt(min(late$entropy), log2(3) - 0.15)              # h = log2(3) regime
})

test_that("trajectory respects actor boundaries and group split", {
  d <- .make_stream(.P3, 3000, seed = 5)
  d$who <- rep(c("a", "b"), each = 1500)
  d$grp <- rep(c("G1", "G2"), each = 1500)
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           group = "grp", window = 300, step = 300)
  expect_setequal(unique(tr$trajectory$group), c("G1", "G2"))
  ## 1499 transitions per group (actor boundary breaks one pair)
  expect_true(all(table(tr$trajectory$group) ==
                    length(seq.int(1L, 1499 - 300 + 1, by = 300))))
})

test_that("trajectory falls back to one window on short streams", {
  d <- .make_stream(.P3, 120, seed = 6)
  expect_warning(
    tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                             window = 500),
    regexp = "shorter than window")
  expect_equal(nrow(tr$trajectory), 1L)
  expect_equal(tr$trajectory$n_transitions, 119L)
})

test_that("trajectory validates input and S3 methods dispatch", {
  d <- .make_stream(.P3, 1500, seed = 7)
  expect_error(entropy_trajectory(d, action = "nope"))
  expect_error(entropy_trajectory(d, action = "act", window = 5))
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           window = 300, step = 150)
  expect_output(print(tr), "Entropy Trajectory")
  s <- summary(tr)
  expect_true(is.data.frame(s))
  expect_named(s, c("group", "windows", "mean", "sd", "min", "max",
                    "first", "last", "change"))
  g <- plot(tr)
  expect_true(inherits(g, "gg"))
  ## entropy_norm within [0, 1]
  expect_true(all(tr$trajectory$entropy_norm >= 0 - 1e-12))
  expect_true(all(tr$trajectory$entropy_norm <= 1 + 1e-12))
})

# ---- surprisal chance scaling ----

test_that("chance-scaled surprisal is -log(P)/log(n), 1 at chance level", {
  ent <- entropy_network(.P3, weight = "surprisal", scaling = "chance")
  expect_equal(ent$weights, -log2(.P3) / log2(3), tolerance = 1e-12)
  ## An exactly chance-level transition (P = 1/n) scores exactly 1
  P <- matrix(1 / 3, 3, 3, dimnames = dimnames(.P3))
  entu <- entropy_network(P, weight = "surprisal", scaling = "chance")
  expect_true(all(abs(entu$weights - 1) < 1e-12))
  ## Provenance recorded
  expect_equal(ent$params$scaling, "chance")
  expect_equal(ent$scaling, "chance")
})

test_that("chance scaling is rejected for non-surprisal weights", {
  expect_error(entropy_network(.P3, weight = "contribution",
                               scaling = "chance"), "surprisal")
  expect_error(entropy_network(.P3, weight = "production",
                               scaling = "chance"), "surprisal")
})

test_that("chance-scaled surprisal ships the colored rendering contract", {
  ent <- entropy_network(.P3, weight = "surprisal", scaling = "chance")
  spec <- ent$meta$splot
  ## Drawn edge set/widths come from the probabilities via the weight field
  expect_equal(spec$weight, "prob")
  expect_equal(ent$prob, .P3, tolerance = 1e-12)
  ## Below-chance edges dashed grey, above-chance solid blues
  W <- ent$weights
  expect_true(all(spec$defaults$edge_style[W >= 1] == 2))
  expect_true(all(spec$defaults$edge_style[W < 1] == 1))
  expect_true(all(spec$defaults$edge_color[W >= 1] == "grey72"))
  expect_true(all(grepl("^#", spec$defaults$edge_color[W < 1])))
  ## Stronger transitions get darker color (monotone ramp)
  ramp_pos <- function(col) which(
    grDevices::colorRampPalette(c("#C6DBEF", "#08306B"))(100) == col)[1]
  ## chance for 3 states is 1/3: compare P = 0.7 vs P = 0.5 (both above)
  strong <- spec$defaults$edge_color[which.max(.P3)]
  weak   <- spec$defaults$edge_color[.P3 == 0.5][1]
  expect_gt(ramp_pos(strong), ramp_pos(weak))
})

test_that("share scaling gives percentages of H summing to 100", {
  ent <- entropy_network(.P3, scaling = "share")
  expect_equal(sum(ent$weights), 100, tolerance = 1e-10)
  ## Proportions identical to the bits decomposition
  bits <- entropy_network(.P3)
  expect_equal(ent$weights / 100, bits$weights / sum(bits$weights),
               tolerance = 1e-12)
  expect_equal(ent$scaling, "share")
  expect_equal(ent$meta$splot$defaults$edge_label_digits, 1)
  ## share is contribution-only
  expect_error(entropy_network(.P3, weight = "surprisal",
                               scaling = "share"), "contribution")
})
