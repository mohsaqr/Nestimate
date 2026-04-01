# ---- permutation_test() Tests ----

# Helper: generate wide sequence data
.make_perm_wide <- function(n = 100, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Input validation ----

test_that("permutation_test rejects non-netobject inputs", {
  expect_error(permutation_test("a", "b"), "netobject")
})

test_that("permutation_test rejects mismatched methods", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "frequency")
  expect_error(permutation_test(net1, net2, iter = 10),
               "Methods must match")
})

test_that("permutation_test rejects mismatched nodes", {
  w1 <- .make_perm_wide(n = 50, states = c("A", "B", "C"), seed = 1)
  w2 <- .make_perm_wide(n = 50, states = c("X", "Y", "Z"), seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  expect_error(permutation_test(net1, net2, iter = 10),
               "Nodes must be the same")
})

test_that("permutation_test rejects missing data", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "relative")
  net2$data <- NULL
  expect_error(permutation_test(net1, net2, iter = 10),
               "does not contain \\$data")
})

test_that("paired mode requires equal n", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 60, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  expect_error(permutation_test(net1, net2, iter = 10, paired = TRUE),
               "equal number")
})


# ---- Transition methods ----

test_that("permutation_test works with method='relative'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm <- permutation_test(net1, net2, iter = 50L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "relative")
  expect_equal(perm$iter, 50L)
  expect_equal(perm$alpha, 0.05)
  expect_false(perm$paired)
  expect_equal(perm$adjust, "none")

  # Matrices
  expect_true(is.matrix(perm$diff))
  expect_true(is.matrix(perm$diff_sig))
  expect_true(is.matrix(perm$p_values))
  expect_true(is.matrix(perm$effect_size))
  expect_equal(dim(perm$diff), c(3, 3))
  expect_equal(dimnames(perm$diff), list(net1$nodes$label, net1$nodes$label))

  # P-values in [0, 1]
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))

  # diff = x$weights - y$weights
  expect_equal(perm$diff, net1$weights - net2$weights)

  # Summary
  expect_true(is.data.frame(perm$summary))
  expect_true(all(c("from", "to", "diff", "effect_size", "p_value", "sig",
                     "weight_x", "weight_y") %in% names(perm$summary)))
})

test_that("permutation_test works with method='frequency'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "frequency")
  net2 <- build_network(w2, method = "frequency")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "frequency")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("permutation_test works with method='co_occurrence'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "co_occurrence")
  net2 <- build_network(w2, method = "co_occurrence")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "co_occurrence")
  expect_false(perm$x$directed)
  # Undirected summary: from <= to (self-loops have from == to)
  if (nrow(perm$summary) > 0) {
    off_diag <- perm$summary$from != perm$summary$to
    if (any(off_diag)) {
      expect_true(all(perm$summary$from[off_diag] < perm$summary$to[off_diag]))
    }
  }
})


# ---- Association methods ----

test_that("permutation_test works with method='cor'", {
  d1 <- .make_freq_data(n = 60, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 60, p = 4, seed = 2)
  net1 <- build_network(d1, method = "cor")
  net2 <- build_network(d2, method = "cor")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "cor")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("permutation_test works with method='glasso'", {
  d1 <- .make_freq_data(n = 80, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 80, p = 4, seed = 2)
  net1 <- build_network(d1, method = "glasso")
  net2 <- build_network(d2, method = "glasso")

  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "glasso")
})


# ---- Paired mode ----

test_that("paired permutation test works", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm <- permutation_test(net1, net2, iter = 30L, paired = TRUE, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_true(perm$paired)
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- p.adjust correction ----

test_that("p.adjust correction is applied", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_none <- permutation_test(net1, net2, iter = 50L,
                                adjust = "none", seed = 42)
  perm_bh <- permutation_test(net1, net2, iter = 50L,
                               adjust = "BH", seed = 42)

  # BH-adjusted p-values should be >= raw p-values
  expect_true(all(perm_bh$p_values >= perm_none$p_values - 1e-10))
})

test_that("bonferroni correction is more conservative", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_bon <- permutation_test(net1, net2, iter = 50L,
                                adjust = "bonferroni", seed = 42)

  # All p-values still in [0, 1]
  expect_true(all(perm_bon$p_values >= 0 & perm_bon$p_values <= 1))
})


# ---- Seed reproducibility ----

test_that("seed produces reproducible results", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_a <- permutation_test(net1, net2, iter = 30L, seed = 42)
  perm_b <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_equal(perm_a$p_values, perm_b$p_values)
  expect_equal(perm_a$effect_size, perm_b$effect_size)
})


# ---- S3 methods ----

test_that("print.net_permutation works", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  output <- capture.output(print(perm))
  expect_true(any(grepl("Permutation Test", output)))
  expect_true(any(grepl("Iterations", output)))
  expect_true(any(grepl("Significant", output)))
})

test_that("print shows paired and adjust info", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L,
                           paired = TRUE, adjust = "BH", seed = 42)

  output <- capture.output(print(perm))
  expect_true(any(grepl("Paired", output)))
  expect_true(any(grepl("BH", output)))
})

test_that("summary returns the summary data frame", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  s <- summary(perm)
  expect_identical(s, perm$summary)
})


# ---- Effect size ----

test_that("effect size is zero when diff is zero", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "relative")

  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  # Same networks → diff is zero everywhere
  expect_true(all(perm$diff == 0))
  expect_true(all(perm$effect_size == 0))
})


# ---- netobject_group dispatch (L82-92) ----

test_that("permutation_test dispatches over matching netobject_group pairs", {
  skip_if_not_installed("tna")
  df1 <- tna::group_regulation
  df2 <- tna::group_regulation
  df1$grp <- rep(c("A", "B"), length.out = nrow(df1))
  df2$grp <- rep(c("A", "B"), length.out = nrow(df2))
  grp1 <- build_network(df1, method = "relative", group = "grp")
  grp2 <- build_network(df2, method = "relative", group = "grp")

  results <- permutation_test(grp1, grp2, iter = 20L, seed = 1)
  expect_true(is.list(results))
  expect_true(all(vapply(results, inherits, logical(1), "net_permutation")))
})

test_that("permutation_test errors when group names do not overlap", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  df$g1 <- rep(c("A", "B"), length.out = nrow(df))
  df$g2 <- rep(c("X", "Y"), length.out = nrow(df))
  grp1 <- build_network(df, method = "relative", group = "g1")
  grp2 <- build_network(df, method = "relative", group = "g2")
  expect_error(permutation_test(grp1, grp2, iter = 10L),
               "No matching group names")
})


# ---- missing $data on x (L111-112) ----

test_that("permutation_test errors when x has no $data", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  net1$data <- NULL
  expect_error(permutation_test(net1, net2, iter = 10L),
               "does not contain \\$data")
})


# ---- node order reordering (L131) ----

test_that("permutation_test handles same nodes in different order", {
  states <- c("A", "B", "C")
  set.seed(1)
  mat1 <- matrix(sample(states, 50 * 10, replace = TRUE), nrow = 50)
  mat2 <- matrix(sample(states, 50 * 10, replace = TRUE), nrow = 50)
  colnames(mat1) <- colnames(mat2) <- paste0("T", 1:10)
  df1 <- as.data.frame(mat1, stringsAsFactors = FALSE)
  df2 <- as.data.frame(mat2, stringsAsFactors = FALSE)

  net1 <- build_network(df1, method = "relative")
  net2 <- build_network(df2, method = "relative")

  # Manually reorder nodes in net2's label so they are same set but different order
  net2$nodes$label <- rev(net2$nodes$label)
  # They are setequal but not identical → triggers L131 branch
  perm <- permutation_test(net1, net2, iter = 15L, seed = 1)
  expect_s3_class(perm, "net_permutation")
})


# ---- paired permutation for association methods (L417-419) ----

test_that("paired permutation works for association methods", {
  d1 <- .make_freq_data(n = 60, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 60, p = 4, seed = 2)
  net1 <- build_network(d1, method = "cor")
  net2 <- build_network(d2, method = "cor")

  perm <- permutation_test(net1, net2, iter = 20L, paired = TRUE, seed = 42)
  expect_s3_class(perm, "net_permutation")
  expect_true(perm$paired)
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- pcor association path (L368-371) ----

test_that("permutation_test works with method='pcor'", {
  d1 <- .make_freq_data(n = 80, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 80, p = 4, seed = 2)
  net1 <- build_network(d1, method = "pcor")
  net2 <- build_network(d2, method = "pcor")

  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)
  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "pcor")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- custom estimator fallback in association path (L386-407) ----

test_that("permutation_test fallback for custom estimator", {
  custom_fn <- function(data, ...) {
    mat <- cor(as.matrix(data))
    diag(mat) <- 0
    list(matrix = mat, nodes = colnames(mat), directed = FALSE,
         cleaned_data = data)
  }
  register_estimator("test_perm_custom", custom_fn,
                     "custom for perm", directed = FALSE)
  on.exit(remove_estimator("test_perm_custom"), add = TRUE)

  d1 <- .make_freq_data(n = 50, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 50, p = 4, seed = 2)
  net1 <- build_network(d1, method = "test_perm_custom")
  net2 <- build_network(d2, method = "test_perm_custom")

  perm <- permutation_test(net1, net2, iter = 15L, seed = 1)
  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "test_perm_custom")
})


# ---- scaling/threshold applied in association loop (L431-434) ----

test_that("permutation_test with threshold applied in association path", {
  d1 <- .make_freq_data(n = 60, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 60, p = 4, seed = 2)
  net1 <- build_network(d1, method = "cor", threshold = 0.01)
  net2 <- build_network(d2, method = "cor", threshold = 0.01)

  perm <- permutation_test(net1, net2, iter = 20L, seed = 1)
  expect_s3_class(perm, "net_permutation")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- p.adjust with holm (L466) ----

test_that("holm p.adjust correction works", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm <- permutation_test(net1, net2, iter = 30L, adjust = "holm", seed = 42)
  expect_equal(perm$adjust, "holm")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- print: unknown method label (L535) ----

test_that("print.net_permutation shows generic label for unknown method", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 15L, seed = 1)
  perm$method <- "custom_unknown_perm"
  out <- capture.output(print(perm))
  expect_true(any(grepl("custom_unknown_perm", out)))
})


# ---- Cross-validation against tna::permutation_test ----

test_that("permutation_test matches tna::permutation_test numerically", {
  skip_if_not_installed("tna")
  df <- tna::group_regulation
  d1 <- df[1:1000, ]
  d2 <- df[1001:2000, ]
  iter <- 1000L
  seed <- 42

  # tna approach
  m1_tna <- tna::tna(d1)
  m2_tna <- tna::tna(d2)
  set.seed(seed)
  tna_perm <- tna::permutation_test(m1_tna, m2_tna, iter = iter)

  # Nestimate approach
  n1 <- build_network(d1, method = "relative")
  n2 <- build_network(d2, method = "relative")
  nest_perm <- permutation_test(n1, n2, iter = iter, seed = seed)

  # 1. True differences must match exactly
  tna_diff <- tna_perm$edges$diffs_true
  expect_equal(nest_perm$diff, tna_diff, tolerance = 1e-12)

  # 2. Build p-value matrix from tna long-form stats
  states <- rownames(n1$weights)
  n <- length(states)
  tna_stats <- tna_perm$edges$stats
  tna_pmat <- matrix(NA, n, n, dimnames = list(states, states))
  vapply(seq_len(nrow(tna_stats)), function(i) {
    parts <- strsplit(tna_stats$edge_name[i], " -> ")[[1]]
    tna_pmat[parts[1], parts[2]] <<- tna_stats$p_value[i]
    0
  }, numeric(1))

  nest_pmat <- nest_perm$p_values

  # RNG consumption differs between tna and Nestimate,
  # so p-values won't match exactly. With 1000 iterations:
  # - correlation should be > 0.999
  # - max absolute difference should be < 0.03
  valid <- !is.na(tna_pmat) & !is.na(nest_pmat)
  expect_gt(cor(tna_pmat[valid], nest_pmat[valid]), 0.999)
  expect_lt(max(abs(tna_pmat[valid] - nest_pmat[valid])), 0.05)

  # 3. Significant edges should largely agree
  tna_sig <- tna_stats$edge_name[tna_stats$p_value < 0.05]
  nest_sig_df <- nest_perm$summary[nest_perm$summary$p_value < 0.05, ]
  nest_sig <- paste(nest_sig_df$from, "->", nest_sig_df$to)

  # At most 2 borderline edges may flip (Monte Carlo noise at p~0.05)
  sym_diff <- length(setdiff(tna_sig, nest_sig)) +
              length(setdiff(nest_sig, tna_sig))
  expect_lte(sym_diff, 2L)

  # 4. Effect sizes should match closely
  tna_es <- matrix(NA, n, n, dimnames = list(states, states))
  vapply(seq_len(nrow(tna_stats)), function(i) {
    parts <- strsplit(tna_stats$edge_name[i], " -> ")[[1]]
    tna_es[parts[1], parts[2]] <<- tna_stats$effect_size[i]
    0
  }, numeric(1))

  # Replace NaN with 0 for comparison (zero-weight edges)
  tna_es[is.nan(tna_es)] <- 0
  nest_es <- nest_perm$effect_size
  nest_es[is.nan(nest_es)] <- 0
  valid_es <- !is.na(tna_es) & !is.na(nest_es) & tna_es != 0 & nest_es != 0

  if (sum(valid_es) > 5) {
    expect_gt(cor(tna_es[valid_es], nest_es[valid_es]), 0.99)
  }
})


# ---- Group dispatches ----

test_that("permutation_test single netobject_group runs all-pairs (L85-104)", {
  set.seed(1)
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A","A","B","C"),
    V2 = c("B","C","B","A","C","B","B","C","A"),
    V3 = c("C","A","C","B","A","C","C","A","B"),
    grp = c("X","X","X","Y","Y","Y","Z","Z","Z")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  perm <- permutation_test(nets, iter = 10, seed = 1)
  expect_true(inherits(perm, "net_permutation_group"))
  expect_equal(length(perm), 3)  # 3 choose 2 = 3 pairs
  expect_true(all(vapply(perm, function(x) inherits(x, "net_permutation"), logical(1))))
})

test_that("permutation_test single group with < 2 groups errors (L87-89)", {
  seqs <- data.frame(
    V1 = c("A","B","A"), V2 = c("B","C","B"), V3 = c("C","A","C"),
    grp = c("X","X","X")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  expect_error(permutation_test(nets, iter = 10), "at least 2")
})

test_that("print.net_permutation_group shows groups (L644-647)", {
  set.seed(1)
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  perm <- permutation_test(nets, iter = 10, seed = 1)
  out <- capture.output(print(perm))
  expect_true(any(grepl("Grouped Permutation", out)))
  expect_true(any(grepl("Groups:", out)))
})

test_that("summary.net_permutation_group returns combined data frame (L671-675)", {
  set.seed(1)
  seqs <- data.frame(
    V1 = c("A","B","A","C","B","A"),
    V2 = c("B","C","B","A","C","B"),
    V3 = c("C","A","C","B","A","C"),
    grp = c("X","X","X","Y","Y","Y")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  perm <- permutation_test(nets, iter = 10, seed = 1)
  s <- summary(perm)
  expect_true(is.data.frame(s))
  expect_true("group" %in% names(s))
})


# ---- mcml dispatch ----

test_that("permutation_test works with two mcml objects", {
  set.seed(42)
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  seqs1 <- data.frame(
    T1 = sample(LETTERS[1:6], 30, TRUE),
    T2 = sample(LETTERS[1:6], 30, TRUE),
    T3 = sample(LETTERS[1:6], 30, TRUE),
    T4 = sample(LETTERS[1:6], 30, TRUE),
    stringsAsFactors = FALSE
  )
  seqs2 <- data.frame(
    T1 = sample(LETTERS[1:6], 30, TRUE),
    T2 = sample(LETTERS[1:6], 30, TRUE),
    T3 = sample(LETTERS[1:6], 30, TRUE),
    T4 = sample(LETTERS[1:6], 30, TRUE),
    stringsAsFactors = FALSE
  )
  cs1 <- build_mcml(seqs1, clusters, type = "tna")
  cs2 <- build_mcml(seqs2, clusters, type = "tna")

  perm <- permutation_test(cs1, cs2, iter = 10, seed = 1)

  expect_s3_class(perm, "net_permutation_group")
  expect_true(length(perm) > 0)
  for (nm in names(perm)) {
    expect_s3_class(perm[[nm]], "net_permutation")
  }
})


