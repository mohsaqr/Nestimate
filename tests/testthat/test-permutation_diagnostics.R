# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

# ---- permutation_diagnostics(): does nesting bias a permutation test? ----

# Sessions nested in persons; person-specific transition matrices with
# log-scale SD `tau`. H0 true: groups never differ.
.sim_diag_long <- function(n_person = 20L, n_session = 6L, len = 10L,
                           tau = 1, design = "between", seed = 1) {
  set.seed(seed)
  states <- c("A", "B", "C")
  base <- matrix(c(.5, .3, .2, .2, .5, .3, .3, .2, .5), 3, byrow = TRUE)
  sessions <- lapply(seq_len(n_person), function(p) {
    P <- exp(log(base) + matrix(rnorm(9, sd = tau), 3))
    P <- P / rowSums(P)
    lapply(seq_len(n_session), function(s) {
      seq_int <- Reduce(function(prev, i) sample.int(3L, 1L, prob = P[prev, ]),
                        seq_len(len - 1L), accumulate = TRUE,
                        sample.int(3L, 1L))
      cond <- if (design == "between") p %% 2L else s %% 2L
      data.frame(person = p, session = s, t = seq_len(len),
                 Action = states[seq_int], cond = paste0("g", cond))
    })
  })
  build_network(do.call(rbind, unlist(sessions, recursive = FALSE)),
                method = "relative", actor = "person", session = "session",
                order = "t", action = "Action", group = "cond")
}


# ---- the ICC estimator against a hand computation ----

test_that("ANOVA ICC matches the textbook formula on balanced data", {
  set.seed(1)
  cluster <- rep(1:6, each = 4)
  feat <- cbind(a = rnorm(6)[cluster] + rnorm(24), b = rnorm(24))
  # balanced: ICC = (MSB - MSW) / (MSB + (m - 1) MSW), m = 4
  icc_by_hand <- function(v) {
    cm <- tapply(v, cluster, mean)
    msb <- 4 * sum((cm - mean(v))^2) / 5
    msw <- sum((v - cm[cluster])^2) / (24 - 6)
    (msb - msw) / (msb + 3 * msw)
  }
  expect_equal(.anova_icc_cols(feat, cluster),
               c(a = icc_by_hand(feat[, "a"]), b = icc_by_hand(feat[, "b"])),
               tolerance = 1e-12)
})

test_that("ANOVA ICC is NA where a column has no variance", {
  feat <- cbind(const = rep(1, 12), x = c(1, 2, 3, 2, 5, 6, 5, 4, 1, 1, 2, 3))
  icc <- .anova_icc_cols(feat, rep(1:3, each = 4))
  expect_true(is.na(icc[1]))
  expect_true(is.finite(icc[2]))
})

# per-sequence transition counts + person ids, as permutation() sees them
.diag_counts <- function(net) {
  nodes <- net[[1]]$nodes$label
  parts <- lapply(net, function(n) {
    list(counts = .precompute_per_sequence(n$data, "relative", n$params, nodes),
         ids = as.character(n$metadata$person))
  })
  list(counts = do.call(rbind, lapply(parts, function(p) p$counts)),
       ids = unlist(lapply(parts, function(p) p$ids), use.names = FALSE),
       st = rep(names(parts), vapply(parts, function(p) nrow(p$counts), 1L)))
}

test_that("ICC is invariant to relabelling and reordering blocks", {
  d <- .diag_counts(.sim_diag_long(seed = 2))
  counts <- d$counts
  ids <- d$ids
  st <- d$st
  ref <- .block_icc(counts, ids, st)
  perm <- rev(seq_along(ids))
  relabelled <- .block_icc(counts[perm, ], paste0("id_", ids[perm]), st[perm])
  expect_equal(relabelled$estimate, ref$estimate, tolerance = 1e-12)
  expect_equal(relabelled$ci, ref$ci, tolerance = 1e-12)
  expect_lte(ref$ci[1], ref$estimate)
  expect_gte(ref$ci[2], ref$estimate)
})


# ---- the verb ----

test_that("overall level returns one tidy row per pair", {
  net <- .sim_diag_long(tau = 1, seed = 3)
  d <- permutation_diagnostics(net, actor = "person", iter = 99L, seed = 1)
  expect_s3_class(d, "data.frame")
  expect_identical(nrow(d), 1L)
  expect_identical(names(d), c(
    "pair", "n_sequences", "n_actors", "design", "icc", "icc_ci_lower",
    "icc_ci_upper", "deff_edges", "deff_global", "p_global_sequence",
    "p_global_actor", "sig_edges_sequence", "sig_edges_actor",
    "edges_changed", "min_p_actor"))
  expect_identical(d$pair, "g1 vs g0")
  expect_identical(d$n_actors, 20L)
  expect_identical(d$design, "between")
  expect_equal(d$min_p_actor, 1 / 100)
})

test_that("diagnostics agree with the two permutation() runs they wrap", {
  net <- .sim_diag_long(tau = 1, seed = 4)
  d <- permutation_diagnostics(net, actor = "person", iter = 99L, seed = 7)
  plain <- permutation(net, iter = 99L, seed = 7)[[1]]
  blocked <- permutation(net, iter = 99L, seed = 7, actor = "person")[[1]]
  expect_equal(d$p_global_sequence, plain$global$p_value[1])
  expect_equal(d$p_global_actor, blocked$global$p_value[1])
  expect_equal(d$deff_global, blocked$clustering$deff_global)
  expect_equal(d$icc, blocked$clustering$icc)
  expect_identical(d$sig_edges_sequence, sum(plain$summary$sig))
})

test_that("clustering shows up as design effect > 1 between, < 1 within", {
  between <- permutation_diagnostics(.sim_diag_long(tau = 2, seed = 5),
                                     actor = "person", iter = 199L, seed = 1)
  within <- permutation_diagnostics(
    .sim_diag_long(tau = 2, design = "within", seed = 5),
    actor = "person", iter = 199L, seed = 1)
  expect_gt(between$icc, 0.1)
  expect_gt(between$deff_global, 1.5)
  expect_gt(between$deff_edges, 1.5)
  expect_identical(within$design, "within")
  expect_lt(within$deff_global, 1)
})

test_that("edges level returns one row per edge with the null SDs", {
  net <- .sim_diag_long(seed = 6)
  e <- permutation_diagnostics(net, actor = "person", iter = 99L,
                               level = "edges", seed = 1)
  plain <- permutation(net, iter = 99L, seed = 1)[[1]]
  expect_identical(nrow(e), nrow(plain$summary))
  expect_identical(e$from, plain$summary$from)
  expect_identical(names(e), c("pair", "from", "to", "diff", "icc",
                               "null_sd_sequence", "null_sd_actor", "deff",
                               "p_sequence", "p_actor", "changed"))
  expect_equal(e$deff, (e$null_sd_actor / e$null_sd_sequence)^2)
  expect_identical(e$changed, (e$p_sequence < 0.05) != (e$p_actor < 0.05))
})

test_that("two netobjects can be diagnosed directly", {
  net <- .sim_diag_long(seed = 7)
  d <- permutation_diagnostics(net[[1]], net[[2]], actor = "person",
                               iter = 50L, seed = 1)
  expect_identical(d$pair, "x vs y")
})


# ---- error paths, by class ----

test_that("unsupported or missing block inputs raise classed errors", {
  net <- .sim_diag_long(seed = 8)
  expect_error(permutation_diagnostics(net, actor = "nope", iter = 10L),
               class = "nestimate_actor_missing")
  set.seed(9)
  d <- data.frame(a = rnorm(40), b = rnorm(40), c = rnorm(40))
  n1 <- build_network(d, method = "cor")
  expect_error(permutation_diagnostics(n1, n1, actor = "id", iter = 10L),
               class = "nestimate_actor_unsupported")
  expect_error(permutation_diagnostics(net, actor = c("a", "b")),
               "single column name")
})
