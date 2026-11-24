# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

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
               class = "nestimate_bad_input")
  expect_error(transition_entropy(.P3, base = -2),
               class = "nestimate_bad_input")
  expect_error(transition_entropy(.P3, base = c(2, 10)),
               class = "nestimate_bad_input")
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
  expect_identical(class(s), "data.frame")
  ## Tidy table has 6 columns including normalised + percent
  expect_named(s, c("state", "stationary", "row_entropy",
                    "row_entropy_norm",
                    "contribution", "contribution_pct"))
  ## Sorted by contribution_pct descending
  expect_equal(s$contribution_pct,
               sort(s$contribution_pct, decreasing = TRUE))
  ## Contributions sum to entropy_rate; pct sums to 100
  expect_equal(sum(s$contribution), te$entropy_rate, tolerance = 1e-10)
  expect_equal(sum(s$contribution_pct), 100, tolerance = 1e-10)
  ## Chain table has 4 rows: h(P), H(pi), redundancy, ceiling
  chain <- summary(te, table = "chain")
  expect_identical(class(chain), "data.frame")
  expect_equal(nrow(chain), 4L)
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
  expect_identical(ent$meta[setdiff(names(ent$meta), "splot")],
                   net$meta[setdiff(names(net$meta), "splot")])
  # the source's colours survive the entropy styling
  expect_identical(ent$meta$splot$defaults$node_fill,
                   net$meta$splot$defaults$node_fill)
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

# The positions of a plot's layers of one geom class.
.te_layer <- function(plot, geom) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, geom),
               logical(1)))
}

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

test_that("trajectory groups take the colours the states leave", {
  # three states (A, B, C) take the first three Okabe-Ito colours, so the
  # groups are dealt the six that remain
  spare <- setdiff(.okabe_ito, .okabe_ito[1:3])
  d <- .make_stream(.P3, 4200, seed = 5)
  d$who <- rep(sprintf("a%02d", 1:7), each = 600)
  d$grp <- rep(sprintf("G%02d", 1:7), each = 600)
  six <- subset(d, grp != "G07")
  tr6 <- entropy_trajectory(six, action = "act", actor = "who", time = "t",
                            group = "grp", window = 100, step = 50)
  expect_identical(unname(tr6$state_palette), .okabe_ito[1:3])
  drawn6 <- plot(tr6, iter = 0)
  built6 <- ggplot2::ggplot_build(drawn6)
  colours6 <- unique(built6$data[[.te_layer(drawn6, "GeomLine")[1L]]]$colour)
  expect_equal(nrow(built6$layout$layout), 1L)
  expect_setequal(colours6, spare)
  expect_length(intersect(colours6, .okabe_ito[1:3]), 0L)
  # one more group than those colours: one panel each, one colour
  tr7 <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                            group = "grp", window = 100, step = 50)
  drawn7 <- plot(tr7, iter = 0)
  built7 <- ggplot2::ggplot_build(drawn7)
  expect_equal(nrow(built7$layout$layout), 7L)
  expect_equal(unique(built7$data[[.te_layer(drawn7, "GeomLine")[1L]]]$colour),
               spare[1L])
  # a network's own state colours are the ones the groups avoid
  net <- build_network(subset(d, grp %in% c("G01", "G02")), method = "relative",
                       actor = "who", action = "act", time = "t",
                       group = "grp")
  net <- set_state_colors(net, c(A = "#0072B2", B = "#D55E00", C = "#CC79A7"))
  from_net <- entropy_trajectory(net, window = 100, step = 50)
  drawn <- plot(from_net, iter = 0)
  colours <- unique(ggplot2::ggplot_build(drawn)$data[[
    .te_layer(drawn, "GeomLine")[1L]]]$colour)
  expect_length(intersect(colours, c("#0072B2", "#D55E00", "#CC79A7")), 0L)
})

test_that("no group of a trajectory plot disappears, and a sparse one is named", {
  # Four groups of very different lengths: 2000, 1200, 800 and 260 events, so
  # with window 250 and step 100 they get 18, 10, 6 and 1 windows.
  sizes <- c(G1 = 2000, G2 = 1200, G3 = 800, G4 = 260)
  d <- .make_stream(.P3, sum(sizes), seed = 9)
  d$grp <- rep(names(sizes), sizes)
  d$who <- d$grp
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           group = "grp", window = 250, step = 100)
  windows <- table(tr$trajectory$group)
  expect_equal(as.vector(windows[c("G1", "G2", "G3", "G4")]), c(18, 10, 6, 1))

  # one actor per group, so the band is left out: this test is about which
  # groups get a trend
  expect_warning(plotted <- plot(tr, iter = 0),
                 class = "nestimate_sparse_trajectory")
  warned <- tryCatch(plot(tr, iter = 0), warning = conditionMessage)
  expect_match(warned, "G3 (6 windows) and G4 (1 window)", fixed = TRUE)
  expect_match(warned, "it needs 10", fixed = TRUE)

  built <- ggplot2::ggplot_build(plotted)
  colour_of <- stats::setNames(setdiff(.okabe_ito, .okabe_ito[1:3])[1:4],
                               names(sizes))
  lines <- .te_layer(plotted, "GeomLine")
  points <- .te_layer(plotted, "GeomPoint")
  layer_groups <- function(layer) {
    drawn <- built$data[[layer]]
    colours <- if ("fill" %in% names(drawn) && all(drawn$colour == "grey20"))
      drawn$fill else drawn$colour
    names(colour_of)[match(unique(colours), colour_of)]
  }
  # trend lines only where span x windows >= 4; the other groups keep every
  # window as a point (G3 also joined by a line), so every group is drawn
  expect_setequal(layer_groups(lines[1L]), c("G1", "G2"))
  expect_setequal(layer_groups(points), c("G3", "G4"))
  expect_equal(nrow(built$data[[points]]),
               as.numeric(windows["G3"] + windows["G4"]))
  expect_setequal(layer_groups(lines[2L]), "G3")
  expect_setequal(unique(c(layer_groups(points), layer_groups(lines[1L]))),
                  names(sizes))
  # every group is named on the plot at the end of its trend or windows
  labels <- built$data[[.te_layer(plotted, "GeomText")]]
  expect_setequal(labels$label, names(sizes))
  # no other condition while the plot is drawn
  expect_no_warning(ggplot2::ggplotGrob(plotted))

  # a larger span gives G3 its curve; G4's one window stays a point
  expect_warning(wider <- plot(tr, span = 0.7, iter = 0),
                 class = "nestimate_sparse_trajectory")
  built <- ggplot2::ggplot_build(wider)
  expect_setequal(layer_groups(.te_layer(wider, "GeomLine")[1L]),
                  c("G1", "G2", "G3"))
  # enough windows everywhere: no warning
  even <- subset(d, grp != "G4" & grp != "G3")
  tr2 <- entropy_trajectory(even, action = "act", actor = "who", time = "t",
                            group = "grp", window = 250, step = 100)
  expect_no_warning(ggplot2::ggplotGrob(plot(tr2, iter = 0)))
})

# The stream of `n_actors` actors, each `length` events of the chain P, in
# `groups` groups; actors' events interleave in time.
.actor_streams <- function(P, actors, length, seed) {
  set.seed(seed)
  states <- rownames(P)
  stationary <- Re(eigen(t(P))$vectors[, 1L])
  stationary <- stationary / sum(stationary)
  do.call(rbind, lapply(names(actors), function(g) {
    do.call(rbind, lapply(seq_len(actors[[g]]), function(a) {
      s <- integer(length)
      s[1L] <- sample.int(nrow(P), 1L, prob = stationary)
      for (k in 2:length) s[k] <- sample.int(nrow(P), 1L, prob = P[s[k - 1L], ])
      data.frame(grp = g, who = sprintf("%s-%03d", g, a),
                 t = seq_len(length) + stats::runif(1), act = states[s],
                 stringsAsFactors = FALSE)
    }))
  }))
}

test_that("the trajectory band is a bootstrap over actors and covers the rate", {
  d <- .actor_streams(.P3, c(G1 = 80), length = 60, seed = 21)
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           group = "grp", window = 1000, step = 100)
  drawn <- plot(tr, seed = 1)
  built <- ggplot2::ggplot_build(drawn)
  band <- built$data[[.te_layer(drawn, "GeomRibbon")]]
  trend <- built$data[[.te_layer(drawn, "GeomLine")[1L]]]
  expect_true(all(band$ymin <= band$ymax))
  expect_true(all(band$ymax - band$ymin > 0))
  # the band covers at least the middle of the grid, and the trend lies in
  # it at the times it covers
  inside <- merge(trend[, c("x", "y")], band[, c("x", "ymin", "ymax")],
                  by = "x")
  expect_gte(nrow(inside), 40L)
  expect_gte(mean(inside$y >= inside$ymin - 1e-12 &
                    inside$y <= inside$ymax + 1e-12), 0.9)
  # the band holds the chain's true entropy rate in its middle
  stationary <- Re(eigen(t(.P3))$vectors[, 1L])
  stationary <- stationary / sum(stationary)
  rate <- -sum(stationary * rowSums(.P3 * log2(.P3)))
  middle <- band[ceiling(nrow(band) / 2), ]
  expect_true(middle$ymin <= rate && rate <= middle$ymax)
  expect_match(drawn$labels$subtitle, "band: 95% of 200 resamples",
               fixed = TRUE)
})

test_that("the band widens with fewer actors and is reproducible", {
  d <- .actor_streams(.P3, c(many = 80, few = 20), length = 60, seed = 22)
  tr <- entropy_trajectory(d, action = "act", actor = "who", time = "t",
                           group = "grp", window = 100, step = 20)
  drawn <- plot(tr, seed = 3)
  band <- ggplot2::ggplot_build(drawn)$data[[.te_layer(drawn, "GeomRibbon")]]
  width <- tapply(band$ymax - band$ymin, band$fill, mean)
  colour_of <- stats::setNames(setdiff(.okabe_ito, .okabe_ito[1:3])[1:2],
                               c("few", "many"))
  expect_gt(width[[colour_of[["few"]]]], width[[colour_of[["many"]]]])

  # the same seed gives the same band, and the caller's stream is untouched
  again <- plot(tr, seed = 3)
  expect_identical(
    ggplot2::ggplot_build(again)$data[[.te_layer(again, "GeomRibbon")]],
    band)
  set.seed(99); expected <- stats::runif(1)
  set.seed(99); invisible(plot(tr, seed = 3)); observed <- stats::runif(1)
  expect_identical(observed, expected)

  # iter = 0: no band
  bare <- plot(tr, iter = 0)
  expect_length(.te_layer(bare, "GeomRibbon"), 0L)
  expect_false(grepl("band", bare$labels$subtitle))
  expect_error(plot(tr, ci_level = 1.2), class = "nestimate_bad_input")
})

# Actors of two lengths (10 or 80 actions), so a resample of the actors holds
# a different number of transitions, and of windows, than the data.
.uneven_actors <- function(times = NULL) {
  set.seed(5)
  lengths <- sample(c(10L, 80L), 30L, replace = TRUE)
  do.call(rbind, lapply(seq_along(lengths), function(a) {
    n <- lengths[[a]]
    data.frame(grp = "G", who = sprintf("a%02d", a),
               t = if (is.null(times)) seq_len(n) + stats::runif(1) else
                 times(n),
               act = sample(c("A", "B", "C"), n, replace = TRUE),
               stringsAsFactors = FALSE)
  }))
}

test_that("resamples left out of the band are counted and reported", {
  tr <- entropy_trajectory(.uneven_actors(), action = "act", actor = "who",
                           time = "t", group = "grp", window = 200,
                           step = 100)
  # 10 windows: the data has a trend at span 0.4, most resamples do not
  expect_identical(nrow(tr$trajectory), 10L)
  expect_warning(drawn <- plot(tr, iter = 50, seed = 1),
                 class = "nestimate_dropped_resamples")
  expect_no_error(ggplot2::ggplot_build(drawn))
  # with no resamples there is nothing to report
  expect_no_warning(plot(tr, iter = 0))
})

test_that("a trend loess cannot fit is named with its reason, not an error", {
  # Three distinct times: loess warns, and its prediction would then fail.
  tr <- entropy_trajectory(
    .uneven_actors(times = function(n) rep(c(1, 2), length.out = n)),
    action = "act", actor = "who", time = "t", group = "grp", window = 50,
    step = 25)
  expect_gte(nrow(tr$trajectory) * 0.4, 4)
  expect_warning(drawn <- plot(tr, iter = 0),
                 class = "nestimate_trend_not_fitted")
  built <- ggplot2::ggplot_build(drawn)
  points <- built$data[[.te_layer(drawn, "GeomPoint")]]
  expect_identical(nrow(points), nrow(tr$trajectory))
})

# Actors observed at the same times, so the stream has ties between actors.
.synchronised_actors <- function(states = c("A", "B", "C"), groups = NULL,
                                 n_actors = 20L, length = 40L) {
  set.seed(1)
  do.call(rbind, lapply(seq_len(n_actors), function(a) {
    data.frame(who = sprintf("a%02d", a), t = seq_len(length),
               g = if (is.null(groups)) "all" else
                 groups[[1L + (a > n_actors / 2)]],
               act = sample(states, length, replace = TRUE),
               stringsAsFactors = FALSE)
  }))
}

test_that("a resample drawing every actor once rebuilds the data's windows", {
  # windows that cut through times shared by several actors
  tr <- entropy_trajectory(.synchronised_actors(), action = "act",
                           actor = "who", time = "t", window = 50, step = 7)
  stream <- tr$stream
  by_actor <- split(seq_along(stream$actor), stream$actor)
  every_once <- rev(seq_along(by_actor))
  windows <- .te_windows(stream, .te_draw_rows(by_actor, every_once),
                         tr$window, tr$step, tr$base, warn = FALSE)
  expect_equal(windows$entropy, tr$trajectory$entropy, tolerance = 1e-12)
  # an actor drawn twice: its copies sit together, in time order
  twice <- .te_draw_rows(by_actor, c(1L, 1L))
  expect_identical(twice, rep(by_actor[[1L]], each = 2L))
})

test_that("the band is the basic interval, reflected about the trend", {
  grid <- seq_len(5L)
  # resampled trends shifted 0.1 below an observed trend of 1
  draws <- lapply(seq(0.85, 0.95, length.out = 101L), rep, times = 5L)
  band <- .te_band("G", draws, 101L, grid, rep(1, 5L), 0.95, "G")
  expect_equal(band$lower, rep(2 - 0.9475, 5L), tolerance = 1e-12)
  expect_equal(band$upper, rep(2 - 0.8525, 5L), tolerance = 1e-12)
  # a time reached by fewer than half of all the draws has no band
  draws[seq_len(60L)] <- lapply(draws[seq_len(60L)], replace, 1L, NA)
  expect_identical(.te_band("G", draws, 101L, grid, rep(1, 5L), 0.95,
                            "G")$time, 2:5)
})

test_that("a group with one actor gets no band, and says so", {
  tr <- entropy_trajectory(.make_stream(.P3, 1500, seed = 8), action = "act",
                           actor = "who", time = "t", window = 100,
                           step = 50)
  expect_warning(drawn <- plot(tr, seed = 1),
                 class = "nestimate_single_actor_band")
  expect_length(.te_layer(drawn, "GeomRibbon"), 0L)
  expect_false(grepl("band", drawn$labels$subtitle))
})

test_that("a missing time is refused", {
  d <- .synchronised_actors()
  d$t[[5L]] <- NA
  expect_error(entropy_trajectory(d, action = "act", actor = "who",
                                  time = "t"),
               class = "nestimate_bad_input")
})

test_that("a group labelled \"\" gets its trend", {
  tr <- entropy_trajectory(.synchronised_actors(groups = c("", "x")),
                           action = "act", actor = "who", time = "t",
                           group = "g", window = 50, step = 10)
  drawn <- plot(tr, iter = 0)
  trend <- ggplot2::ggplot_build(drawn)$data[[
    .te_layer(drawn, "GeomLine")[1L]]]
  expect_identical(length(unique(trend$colour)), 2L)
})

test_that("groups are faceted when the states take every colour", {
  tr <- entropy_trajectory(.synchronised_actors(LETTERS[1:9], c("x", "y")),
                           action = "act", actor = "who", time = "t",
                           group = "g", window = 50, step = 10)
  expect_setequal(unname(tr$state_palette), .okabe_ito)
  expect_no_warning(drawn <- plot(tr, iter = 0))
  expect_s3_class(drawn$facet, "FacetWrap")
  colours <- unique(unlist(lapply(ggplot2::ggplot_build(drawn)$data,
                                  function(layer) layer$colour)))
  expect_identical(intersect(colours, setdiff(.okabe_ito, "#000000")),
                   character(0))
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
  # 8 windows at span 0.4 are too few for a trend line: said, not silent.
  expect_warning(g <- plot(tr), class = "nestimate_sparse_trajectory")
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

# ---- entropy_trajectory() on sequences ----

test_that("entropy_trajectory reads a network the same as its event log", {
  # asymmetric sequences: read backwards, their transitions and windows change
  seqs <- data.frame(T1 = c("A", "A", "B", "A", "C", "B"),
                     T2 = c("A", "B", "C", "A", "A", "B"),
                     T3 = c("B", "B", "C", "A", "B", "C"),
                     T4 = c("C", "C", "A", "B", "B", "A"),
                     T5 = c("C", "A", "A", "C", "B", "B"),
                     T6 = c("A", "A", "B", "C", "C", "C"),
                     stringsAsFactors = FALSE)
  net <- build_network(seqs, method = "relative")
  log <- data.frame(id = rep(1:6, times = 6), pos = rep(1:6, each = 6),
                    st = unlist(seqs, use.names = FALSE),
                    stringsAsFactors = FALSE)
  from_net <- entropy_trajectory(net, window = 10, step = 5)
  from_log <- entropy_trajectory(log, action = "st", actor = "id",
                                 time = "pos", window = 10, step = 5)
  expect_equal(nrow(from_net$trajectory), 5L)
  expect_equal(from_net$trajectory, from_log$trajectory)
  # a wide table and a character matrix give the same windows as the network
  expect_equal(entropy_trajectory(seqs, window = 10, step = 5)$trajectory,
               from_net$trajectory)
  expect_equal(entropy_trajectory(as.matrix(seqs), window = 10,
                                  step = 5)$trajectory,
               from_net$trajectory)
})

test_that("entropy_trajectory on a window of all transitions equals the plug-in rate", {
  # hand value: one window holding every transition is the empirical
  # conditional entropy of the pooled counts
  seqs <- data.frame(T1 = c("A", "A", "C"), T2 = c("B", "B", "A"),
                     T3 = c("C", "B", "B"), T4 = c("A", "C", "A"),
                     stringsAsFactors = FALSE)
  # 9 transitions, fewer than the smallest window: one window, with a warning
  expect_warning(tr <- entropy_trajectory(seqs, window = 10, step = 1),
                 class = "nestimate_short_stream")
  # by hand: A->B 3, B->C 2, C->A 2, B->B 1, B->A 1 (read backwards the
  # counts differ, so the order of the positions is tested too)
  n_ij <- c(AB = 3, BC = 2, CA = 2, BB = 1, BA = 1)
  n_i <- c(AB = 3, BC = 4, CA = 2, BB = 4, BA = 4)
  h <- -sum(n_ij / 9 * log2(n_ij / n_i))
  expect_equal(tr$trajectory$n_transitions, 9L)
  expect_equal(tr$trajectory$entropy, h)
})

test_that("entropy_trajectory takes the groups of a grouped network", {
  seqs <- data.frame(T1 = rep(c("A", "B", "C"), 4),
                     T2 = rep(c("B", "C", "A"), 4),
                     T3 = rep(c("C", "A", "B"), 4),
                     cond = rep(c("x", "y"), each = 6),
                     stringsAsFactors = FALSE)
  nets <- build_network(seqs, method = "relative", group = "cond")
  tr <- entropy_trajectory(nets, window = 10, step = 2)
  expect_setequal(unique(tr$trajectory$group), c("x", "y"))
  # each group has 6 sequences x 2 transitions = 12; windows of 10, step 2
  expect_equal(as.vector(table(tr$trajectory$group)), c(2L, 2L))
  expect_true(all(tr$trajectory$n_transitions == 10L))
})

test_that("entropy_trajectory refuses a network that holds no sequences", {
  num <- data.frame(a = sin(1:40), b = cos(1:40), c = (1:40) %% 7)
  cor_net <- build_network(num, method = "cor")
  expect_error(entropy_trajectory(cor_net), class = "nestimate_bad_input")
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"),
                     stringsAsFactors = FALSE)
  expect_error(entropy_trajectory(seqs, group = "nope"),
               class = "nestimate_bad_input")
})
