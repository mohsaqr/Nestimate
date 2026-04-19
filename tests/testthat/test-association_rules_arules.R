testthat::skip_on_cran()

test_that("arules is available for cross-validation", {
  skip_if_not_installed("arules")
  expect_true(requireNamespace("arules", quietly = TRUE))
})


.ar_make_key <- function(ante_str, cons_str) {
  ante <- paste(sort(strsplit(ante_str, ", ", fixed = TRUE)[[1]]), collapse = "\t")
  cons <- paste(sort(strsplit(cons_str, ", ", fixed = TRUE)[[1]]), collapse = "\t")
  paste(ante, "=>", cons)
}


.ar_compare <- function(trans_list, min_support, min_confidence,
                         min_lift = 0, max_length = 5L) {
  skip_if_not_installed("arules")

  ours <- association_rules(trans_list, min_support = min_support,
                            min_confidence = min_confidence,
                            min_lift = min_lift, max_length = max_length)

  t_obj <- as(trans_list, "transactions")
  theirs <- arules::apriori(
    t_obj,
    parameter = list(supp = min_support, conf = min_confidence,
                     minlen = 2L, maxlen = max_length, target = "rules"),
    control = list(verbose = FALSE))
  if (min_lift > 0) theirs <- subset(theirs, arules::quality(theirs)$lift >= min_lift)

  q <- arules::quality(theirs)
  q$conviction <- arules::interestMeasure(theirs, "conviction", transactions = t_obj)
  lhs_list <- as(arules::lhs(theirs), "list")
  rhs_list <- as(arules::rhs(theirs), "list")

  arules_key <- vapply(seq_along(lhs_list), function(i) {
    paste(paste(sort(lhs_list[[i]]), collapse = "\t"), "=>",
          paste(sort(rhs_list[[i]]), collapse = "\t"))
  }, character(1))

  r <- ours$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  single_mask <- n_cons == 1L
  our_single <- r[single_mask, , drop = FALSE]

  if (nrow(our_single) == 0 && length(theirs) == 0) {
    return(list(ours = ours, n_arules = 0L, n_single = 0L,
                n_multi = sum(!single_mask), matched = TRUE,
                max_support_diff = 0, max_confidence_diff = 0,
                max_lift_diff = 0, max_conviction_diff = 0,
                missing_from_ours = character(0), extra_in_ours = character(0)))
  }

  our_key <- vapply(seq_len(nrow(our_single)), function(i) {
    .ar_make_key(our_single$antecedent[i], our_single$consequent[i])
  }, character(1))

  missing_from_ours <- setdiff(arules_key, our_key)
  extra_in_ours <- setdiff(our_key, arules_key)
  common <- intersect(our_key, arules_key)

  sup_d <- conf_d <- lift_d <- conv_d <- numeric(length(common))
  for (i in seq_along(common)) {
    key <- common[i]
    oi <- which(our_key == key)[1]; ai <- which(arules_key == key)[1]
    sup_d[i] <- abs(our_single$support[oi] - q$support[ai])
    conf_d[i] <- abs(our_single$confidence[oi] - q$confidence[ai])
    lift_d[i] <- abs(our_single$lift[oi] - q$lift[ai])
    oc <- our_single$conviction[oi]; ac <- q$conviction[ai]
    conv_d[i] <- if (is.infinite(oc) && is.infinite(ac)) 0 else abs(oc - ac)
  }

  list(ours = ours, n_arules = length(theirs), n_single = nrow(our_single),
       n_multi = sum(!single_mask), missing_from_ours = missing_from_ours,
       extra_in_ours = extra_in_ours, n_common = length(common),
       max_support_diff = if (length(sup_d) > 0) max(sup_d) else 0,
       max_confidence_diff = if (length(conf_d) > 0) max(conf_d) else 0,
       max_lift_diff = if (length(lift_d) > 0) max(lift_d) else 0,
       max_conviction_diff = if (length(conv_d) > 0) max(conv_d) else 0,
       matched = length(missing_from_ours) == 0 && length(extra_in_ours) == 0)
}


.ar_filter_single <- function(rules_obj) {
  r <- rules_obj$rules
  if (nrow(r) == 0) return(r)
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  r[n_cons == 1L, , drop = FALSE]
}


test_that("bread/milk/eggs/butter: single-cons match arules exactly", {
  skip_if_not_installed("arules")
  trans <- list(c("bread","milk","eggs"), c("bread","butter"),
                c("milk","eggs","butter"), c("bread","milk","eggs","butter"),
                c("bread","milk"))
  res <- .ar_compare(trans, min_support = 0.3, min_confidence = 0.0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
  expect_true(res$max_confidence_diff < 1e-10)
  expect_true(res$max_lift_diff < 1e-10)
  expect_true(res$max_conviction_diff < 1e-10)
})

test_that("bread/milk: confidence >= 0.5 matches", {
  skip_if_not_installed("arules")
  trans <- list(c("bread","milk","eggs"), c("bread","butter"),
                c("milk","eggs","butter"), c("bread","milk","eggs","butter"),
                c("bread","milk"))
  res <- .ar_compare(trans, min_support = 0.3, min_confidence = 0.5)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
})

test_that("bread/milk: lift >= 1 matches", {
  skip_if_not_installed("arules")
  trans <- list(c("bread","milk","eggs"), c("bread","butter"),
                c("milk","eggs","butter"), c("bread","milk","eggs","butter"),
                c("bread","milk"))
  res <- .ar_compare(trans, min_support = 0.3, min_confidence = 0.0, min_lift = 1.0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("max_length=2: no multi-cons, exact match", {
  skip_if_not_installed("arules")
  trans <- list(c("bread","milk","eggs"), c("bread","butter"),
                c("milk","eggs","butter"), c("bread","milk","eggs","butter"),
                c("bread","milk"))
  res <- .ar_compare(trans, min_support = 0.3, min_confidence = 0.0, max_length = 2L)
  expect_equal(res$n_multi, 0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("dense 4-item: single-cons match", {
  skip_if_not_installed("arules")
  trans <- list(c("A","B","C","D"), c("A","B","C"), c("A","B","D"),
                c("A","C","D"), c("B","C","D"), c("A","B","C","D"))
  res <- .ar_compare(trans, min_support = 0.3, min_confidence = 0.0, max_length = 4L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
  expect_true(res$max_conviction_diff < 1e-10)
})

test_that("10-item planted patterns: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(101)
  items <- LETTERS[1:10]
  trans <- lapply(seq_len(200), function(i) {
    base <- sample(items, sample(2:5, 1))
    if (runif(1) < 0.7) base <- union(base, c("A", "B"))
    if (runif(1) < 0.6) base <- union(base, c("C", "D"))
    if (runif(1) < 0.4) base <- union(base, c("A", "B", "E"))
    sort(unique(base))
  })
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.3, max_length = 4L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
  expect_true(res$max_conviction_diff < 1e-10)
})

test_that("15-item random (500 trans): single-cons match", {
  skip_if_not_installed("arules")
  set.seed(202)
  trans <- lapply(seq_len(500), function(i) {
    sort(sample(paste0("item_", sprintf("%02d", 1:15)), sample(2:7, 1)))
  })
  res <- .ar_compare(trans, min_support = 0.05, min_confidence = 0.3, max_length = 4L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
})

test_that("20-item random (1000 trans), lift>=1: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(303)
  trans <- lapply(seq_len(1000), function(i) sort(sample(paste0("X", 1:20), sample(3:8, 1))))
  res <- .ar_compare(trans, min_support = 0.05, min_confidence = 0.4, min_lift = 1.0, max_length = 3L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_support_diff < 1e-10)
})

test_that("high min_support filters identically", {
  skip_if_not_installed("arules")
  set.seed(404)
  trans <- lapply(seq_len(100), function(i) sort(sample(LETTERS[1:8], sample(2:5, 1))))
  res <- .ar_compare(trans, min_support = 0.4, min_confidence = 0.0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("high min_confidence filters identically", {
  skip_if_not_installed("arules")
  set.seed(505)
  trans <- lapply(seq_len(150), function(i) sort(sample(LETTERS[1:6], sample(2:4, 1))))
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.8)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("max_length=3 caps identically", {
  skip_if_not_installed("arules")
  set.seed(606)
  trans <- lapply(seq_len(200), function(i) sort(sample(LETTERS[1:8], sample(3:6, 1))))
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.3, max_length = 3L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("max_length=5 deep itemsets: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(707)
  trans <- lapply(seq_len(300), function(i) sort(sample(LETTERS[1:6], sample(3:6, 1))))
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.0, max_length = 5L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("sparse data (40% single-item): single-cons match", {
  skip_if_not_installed("arules")
  set.seed(808)
  trans <- lapply(seq_len(200), function(i) {
    if (runif(1) < 0.4) return(sample(LETTERS[1:10], 1))
    sort(sample(LETTERS[1:10], sample(2:4, 1)))
  })
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.3)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("perfect co-occurrence: single-cons match", {
  skip_if_not_installed("arules")
  trans <- c(replicate(50, c("X","Y","Z"), simplify = FALSE),
             replicate(30, c("X","Y"), simplify = FALSE),
             replicate(20, c("A","B"), simplify = FALSE))
  res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.5)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_conviction_diff < 1e-10)
})

test_that("no rules found: both return empty", {
  skip_if_not_installed("arules")
  trans <- lapply(LETTERS[1:10], function(x) x)
  res <- .ar_compare(trans, min_support = 0.5, min_confidence = 0.5)
  expect_equal(res$n_single, 0)
  expect_equal(res$n_arules, 0)
})

test_that("all-identical transactions: single-cons match", {
  skip_if_not_installed("arules")
  trans <- replicate(50, c("A","B","C"), simplify = FALSE)
  res <- .ar_compare(trans, min_support = 0.5, min_confidence = 0.0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("binary matrix: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(909)
  mat <- matrix(rbinom(800, 1, 0.4), nrow = 100, ncol = 8)
  colnames(mat) <- paste0("V", 1:8)
  ours <- association_rules(mat, min_support = 0.1, min_confidence = 0.3, min_lift = 0)
  our_single <- .ar_filter_single(ours)
  trans_list <- lapply(seq_len(100), function(i) colnames(mat)[mat[i, ] == 1])
  trans_list <- trans_list[lengths(trans_list) > 0]
  theirs <- arules::apriori(as(trans_list, "transactions"),
    parameter = list(supp = 0.1, conf = 0.3, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  expect_equal(nrow(our_single), length(theirs))
})

test_that("netobject input: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(111)
  seqs <- data.frame(V1 = sample(LETTERS[1:5], 80, TRUE),
                     V2 = sample(LETTERS[1:5], 80, TRUE),
                     V3 = sample(LETTERS[1:5], 80, TRUE),
                     V4 = sample(LETTERS[1:5], 80, TRUE), stringsAsFactors = FALSE)
  net <- build_network(seqs, method = "relative")
  ours <- association_rules(net, min_support = 0.1, min_confidence = 0.3, min_lift = 0)
  our_single <- .ar_filter_single(ours)
  trans_list <- lapply(seq_len(nrow(seqs)), function(i) {
    vals <- as.character(unlist(seqs[i, ], use.names = FALSE))
    unique(vals[!is.na(vals) & vals != ""])
  })
  theirs <- arules::apriori(as(trans_list, "transactions"),
    parameter = list(supp = 0.1, conf = 0.3, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  expect_equal(nrow(our_single), length(theirs))
})

test_that("conviction across 5 seeds", {
  skip_if_not_installed("arules")
  for (seed in c(1001, 1002, 1003, 1004, 1005)) {
    set.seed(seed)
    trans <- lapply(seq_len(150), function(i) sort(sample(LETTERS[1:8], sample(2:5, 1))))
    res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.0)
    expect_equal(res$n_single, res$n_arules, info = sprintf("seed=%d", seed))
    expect_true(res$matched, info = sprintf("seed=%d", seed))
    expect_true(res$max_conviction_diff < 1e-10, info = sprintf("seed=%d", seed))
  }
})

test_that("count values match exactly", {
  skip_if_not_installed("arules")
  set.seed(1234)
  trans <- lapply(seq_len(200), function(i) sort(sample(LETTERS[1:7], sample(2:5, 1))))
  ours <- association_rules(trans, min_support = 0.1, min_confidence = 0.3, min_lift = 0)
  our_single <- .ar_filter_single(ours)
  t_obj <- as(trans, "transactions")
  theirs <- arules::apriori(t_obj,
    parameter = list(supp = 0.1, conf = 0.3, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  q <- arules::quality(theirs)
  lhs_l <- as(arules::lhs(theirs), "list"); rhs_l <- as(arules::rhs(theirs), "list")
  ar_key <- vapply(seq_along(lhs_l), function(i) {
    paste(paste(sort(lhs_l[[i]]), collapse = "\t"), "=>", paste(sort(rhs_l[[i]]), collapse = "\t"))
  }, character(1))
  our_key <- vapply(seq_len(nrow(our_single)), function(i) {
    .ar_make_key(our_single$antecedent[i], our_single$consequent[i])
  }, character(1))
  common <- intersect(our_key, ar_key)
  expect_true(length(common) > 0)
  diffs <- vapply(common, function(k) {
    abs(our_single$count[which(our_key == k)[1]] - q$count[which(ar_key == k)[1]])
  }, numeric(1))
  expect_true(all(diffs == 0))
})

test_that("coverage matches arules", {
  skip_if_not_installed("arules")
  set.seed(5555)
  trans <- lapply(seq_len(100), function(i) sort(sample(LETTERS[1:6], sample(2:4, 1))))
  ours <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  our_single <- .ar_filter_single(ours)
  t_obj <- as(trans, "transactions")
  theirs <- arules::apriori(t_obj,
    parameter = list(supp = 0.1, conf = 0.0, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  q <- arules::quality(theirs)
  lhs_l <- as(arules::lhs(theirs), "list"); rhs_l <- as(arules::rhs(theirs), "list")
  ar_key <- vapply(seq_along(lhs_l), function(i) {
    paste(paste(sort(lhs_l[[i]]), collapse = "\t"), "=>", paste(sort(rhs_l[[i]]), collapse = "\t"))
  }, character(1))
  our_key <- vapply(seq_len(nrow(our_single)), function(i) {
    .ar_make_key(our_single$antecedent[i], our_single$consequent[i])
  }, character(1))
  common <- intersect(our_key, ar_key)
  expect_true(length(common) > 0)
  diffs <- vapply(common, function(k) {
    oi <- which(our_key == k)[1]; ai <- which(ar_key == k)[1]
    abs(our_single$support[oi] / our_single$confidence[oi] - q$coverage[ai])
  }, numeric(1))
  expect_true(max(diffs) < 1e-10)
})

test_that("Groceries-style 2000 trans: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(7777)
  items <- c("bread","milk","eggs","butter","cheese","yogurt","juice","cereal")
  trans <- lapply(seq_len(2000), function(i) {
    b <- character(0)
    if (runif(1) < 0.6) b <- c(b, "bread"); if (runif(1) < 0.5) b <- c(b, "milk")
    if (runif(1) < 0.3) b <- c(b, "eggs"); if (runif(1) < 0.2) b <- c(b, "butter")
    if (runif(1) < 0.25) b <- c(b, "cheese"); if (runif(1) < 0.15) b <- c(b, "yogurt")
    if (runif(1) < 0.1) b <- c(b, "juice"); if (runif(1) < 0.35) b <- c(b, "cereal")
    if (length(b) == 0) b <- sample(items, 1)
    sort(unique(b))
  })
  res <- .ar_compare(trans, min_support = 0.05, min_confidence = 0.3, min_lift = 1.0, max_length = 4L)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
  expect_true(res$max_conviction_diff < 1e-10)
})

test_that("correlated A->B: single-cons match", {
  skip_if_not_installed("arules")
  set.seed(8888)
  trans <- lapply(seq_len(500), function(i) {
    b <- character(0); has_a <- runif(1) < 0.6
    if (has_a) { b <- c(b, "A"); if (runif(1) < 0.9) b <- c(b, "B") }
    if (runif(1) < 0.3) b <- c(b, "C"); if (runif(1) < 0.2) b <- c(b, "D")
    if (runif(1) < 0.4) b <- c(b, "E"); if (length(b) == 0) b <- "E"
    sort(unique(b))
  })
  res <- .ar_compare(trans, min_support = 0.05, min_confidence = 0.5, min_lift = 1.0)
  expect_equal(res$n_single, res$n_arules)
  expect_true(res$matched)
})

test_that("2-item transactions: exact match", {
  skip_if_not_installed("arules")
  set.seed(9999)
  trans <- lapply(seq_len(300), function(i) sort(sample(LETTERS[1:6], 2)))
  res <- .ar_compare(trans, min_support = 0.05, min_confidence = 0.0, max_length = 2L)
  expect_equal(res$n_multi, 0)
  expect_equal(res$n_single, res$n_arules)
})

test_that("human_long bundled data: single-cons match", {
  skip_if_not_installed("arules")
  data(human_long)
  net <- build_network(human_long, method = "relative",
                       actor = "session_id", action = "cluster", time = "timestamp")
  ours <- association_rules(net, min_support = 0.3, min_confidence = 0.5, min_lift = 1.0)
  our_single <- .ar_filter_single(ours)
  seqs_df <- as.data.frame(net$data, stringsAsFactors = FALSE)
  trans_list <- lapply(seq_len(nrow(seqs_df)), function(i) {
    vals <- as.character(unlist(seqs_df[i, ], use.names = FALSE))
    unique(vals[!is.na(vals) & vals != ""])
  })
  theirs <- arules::apriori(as(trans_list, "transactions"),
    parameter = list(supp = 0.3, conf = 0.5, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  theirs <- subset(theirs, arules::quality(theirs)$lift >= 1.0)
  expect_equal(nrow(our_single), length(theirs))
})

test_that("srl_strategies binary: single-cons match", {
  skip_if_not_installed("arules")
  data(srl_strategies)
  binary <- as.data.frame(lapply(srl_strategies, function(col) as.integer(col > median(col, na.rm = TRUE))))
  ours <- association_rules(binary, min_support = 0.2, min_confidence = 0.5, min_lift = 1.0)
  our_single <- .ar_filter_single(ours)
  trans_list <- lapply(seq_len(nrow(binary)), function(i) colnames(binary)[binary[i, ] == 1])
  theirs <- arules::apriori(as(trans_list, "transactions"),
    parameter = list(supp = 0.2, conf = 0.5, minlen = 2, maxlen = 5),
    control = list(verbose = FALSE))
  theirs <- subset(theirs, arules::quality(theirs)$lift >= 1.0)
  expect_equal(nrow(our_single), length(theirs))
})

test_that("frequent itemset counts match arules", {
  skip_if_not_installed("arules")
  set.seed(4321)
  trans <- lapply(seq_len(200), function(i) sort(sample(LETTERS[1:6], sample(2:5, 1))))
  ours <- association_rules(trans, min_support = 0.15, min_confidence = 0.0, min_lift = 0)
  theirs_fi <- arules::apriori(as(trans, "transactions"),
    parameter = list(supp = 0.15, target = "frequent itemsets", minlen = 1, maxlen = 5),
    control = list(verbose = FALSE))
  expect_equal(nrow(ours$frequent), length(theirs_fi))
})

test_that("10 random seeds: single-cons all match", {
  skip_if_not_installed("arules")
  for (seed in seq(2001, 2010)) {
    set.seed(seed)
    trans <- lapply(seq_len(100), function(i) sort(sample(LETTERS[1:7], sample(2:5, 1))))
    res <- .ar_compare(trans, min_support = 0.1, min_confidence = 0.3)
    expect_equal(res$n_single, res$n_arules, info = sprintf("seed=%d", seed))
    expect_true(res$matched, info = sprintf("seed=%d", seed))
    expect_true(res$max_support_diff < 1e-10, info = sprintf("seed=%d", seed))
    expect_true(res$max_lift_diff < 1e-10, info = sprintf("seed=%d", seed))
  }
})

test_that("multi-consequent support is correct by manual count", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    items <- unique(c(strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]],
                      strsplit(multi$consequent[i], ", ", fixed = TRUE)[[1]]))
    mask <- Reduce(`&`, lapply(items, function(it) trans_mat[, it]))
    expect_equal(multi$support[i], sum(mask) / 5)
  }
})

test_that("multi-consequent confidence = support / ante_support", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    ante_items <- strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]]
    ante_mask <- Reduce(`&`, lapply(ante_items, function(it) trans_mat[, it]))
    expect_equal(multi$confidence[i], multi$support[i] / (sum(ante_mask) / 5), tolerance = 1e-12)
  }
})

test_that("multi-consequent lift is correct", {
  trans <- list(c("A","B","C"), c("A","B","C"), c("A","B"), c("B","C"), c("A","C"))
  rules <- association_rules(trans, min_support = 0.1, min_confidence = 0.0, min_lift = 0)
  r <- rules$rules
  n_cons <- lengths(strsplit(r$consequent, ", ", fixed = TRUE))
  multi <- r[n_cons > 1L, , drop = FALSE]
  if (nrow(multi) == 0) skip("No multi-consequent rules")
  trans_mat <- matrix(FALSE, 5, 3, dimnames = list(NULL, c("A","B","C")))
  for (i in seq_along(trans)) trans_mat[i, trans[[i]]] <- TRUE
  for (i in seq_len(nrow(multi))) {
    ante_items <- strsplit(multi$antecedent[i], ", ", fixed = TRUE)[[1]]
    cons_items <- strsplit(multi$consequent[i], ", ", fixed = TRUE)[[1]]
    ante_sup <- sum(Reduce(`&`, lapply(ante_items, function(it) trans_mat[, it]))) / 5
    cons_sup <- sum(Reduce(`&`, lapply(cons_items, function(it) trans_mat[, it]))) / 5
    expect_equal(multi$lift[i], multi$support[i] / (ante_sup * cons_sup), tolerance = 1e-12)
  }
})

test_that("subscript out of bounds fixed: high min_support with no 2-itemsets", {
  trans <- lapply(seq_len(20), function(i) sample(LETTERS[1:15], sample(2:4, 1)))
  expect_no_error(association_rules(trans, min_support = 0.8, min_confidence = 0.0, min_lift = 0))
})

test_that("$frequent is a tidy data frame", {
  trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
  rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.0, min_lift = 0)
  expect_true(is.data.frame(rules$frequent))
  expect_equal(names(rules$frequent), c("itemset", "size", "support", "count"))
  expect_type(rules$frequent$itemset, "character")
  expect_type(rules$frequent$size, "integer")
  expect_true(all(rules$frequent$support > 0))
  expect_true(all(rules$frequent$size[rules$frequent$size == 1] == 1))
})

test_that("$rules has character antecedent/consequent columns", {
  trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
  rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.0, min_lift = 0)
  expect_type(rules$rules$antecedent, "character")
  expect_type(rules$rules$consequent, "character")
  expect_true(all(nchar(rules$rules$antecedent) > 0))
})