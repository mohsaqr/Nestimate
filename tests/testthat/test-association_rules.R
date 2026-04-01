# ---- Association Rules Tests ----

# Helper: standard test transactions
.make_ar_trans <- function() {
  list(
    c("bread", "milk", "eggs"),
    c("bread", "butter"),
    c("milk", "eggs", "butter"),
    c("bread", "milk", "eggs", "butter"),
    c("bread", "milk")
  )
}

# ---- 1. Basic functionality ----

test_that("association_rules returns correct class and structure", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  expect_s3_class(rules, "net_association_rules")
  expect_true(is.data.frame(rules$rules))
  expect_true(rules$n_rules > 0)
  expect_equal(rules$n_transactions, 5)
  expect_true(all(c("antecedent", "consequent", "support", "confidence",
                     "lift", "conviction") %in% names(rules$rules)))
})


# ---- 2. Support correctness ----

test_that("support values are correct", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  r <- rules$rules
  # {bread} => {milk}: both in transactions 1, 4, 5 → support = 3/5 = 0.6
  bm <- r[vapply(r$antecedent, function(a) identical(a, "bread"), logical(1)) &
           vapply(r$consequent, function(c) identical(c, "milk"), logical(1)), ]
  expect_equal(nrow(bm), 1)
  expect_equal(bm$support, 0.6)
  expect_equal(bm$count, 3)
})


# ---- 3. Confidence correctness ----

test_that("confidence = support / antecedent_support", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  r <- rules$rules
  # {eggs} => {milk}: support = 3/5, P(eggs) = 3/5, conf = (3/5) / (3/5) = 1
  em <- r[vapply(r$antecedent, function(a) identical(a, "eggs"), logical(1)) &
           vapply(r$consequent, function(c) identical(c, "milk"), logical(1)), ]
  expect_equal(em$confidence, 1.0)
})


# ---- 4. Lift correctness ----

test_that("lift = support / (P(A) * P(B))", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  r <- rules$rules
  # {eggs} => {milk}: support = 0.6, P(eggs) = 0.6, P(milk) = 0.8
  # lift = 0.6 / (0.6 * 0.8) = 1.25
  em <- r[vapply(r$antecedent, function(a) identical(a, "eggs"), logical(1)) &
           vapply(r$consequent, function(c) identical(c, "milk"), logical(1)), ]
  expect_equal(em$lift, 1.25)
})


# ---- 5. Conviction correctness ----

test_that("conviction computed correctly", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  r <- rules$rules
  # {eggs} => {milk}: conf = 1.0 → conviction = Inf
  em <- r[vapply(r$antecedent, function(a) identical(a, "eggs"), logical(1)) &
           vapply(r$consequent, function(c) identical(c, "milk"), logical(1)), ]
  expect_equal(em$conviction, Inf)

  # {bread} => {butter}: sup=0.4, P(bread)=0.8, conf=0.5, P(butter)=0.6
  # conviction = (1-0.6)/(1-0.5) = 0.4/0.5 = 0.8
  bb <- r[vapply(r$antecedent, function(a) identical(a, "bread"), logical(1)) &
           vapply(r$consequent, function(c) identical(c, "butter"), logical(1)), ]
  expect_equal(bb$conviction, 0.8)
})


# ---- 6. min_support filters correctly ----

test_that("min_support filters low-support rules", {
  rules_low <- association_rules(.make_ar_trans(), min_support = 0.1,
                                 min_confidence = 0, min_lift = 0)
  rules_high <- association_rules(.make_ar_trans(), min_support = 0.6,
                                  min_confidence = 0, min_lift = 0)
  expect_true(rules_low$n_rules >= rules_high$n_rules)
  expect_true(all(rules_high$rules$support >= 0.6))
})


# ---- 7. min_confidence filters correctly ----

test_that("min_confidence filters low-confidence rules", {
  rules_low <- association_rules(.make_ar_trans(), min_support = 0.3,
                                 min_confidence = 0.3, min_lift = 0)
  rules_high <- association_rules(.make_ar_trans(), min_support = 0.3,
                                  min_confidence = 0.9, min_lift = 0)
  expect_true(rules_low$n_rules >= rules_high$n_rules)
  expect_true(all(rules_high$rules$confidence >= 0.9))
})


# ---- 8. min_lift filters correctly ----

test_that("min_lift filters low-lift rules", {
  rules_all <- association_rules(.make_ar_trans(), min_support = 0.3,
                                 min_confidence = 0, min_lift = 0)
  rules_pos <- association_rules(.make_ar_trans(), min_support = 0.3,
                                 min_confidence = 0, min_lift = 1.0)
  expect_true(rules_all$n_rules >= rules_pos$n_rules)
  expect_true(all(rules_pos$rules$lift >= 1.0))
})


# ---- 9. max_length limits itemset size ----

test_that("max_length limits rule complexity", {
  rules2 <- association_rules(.make_ar_trans(), min_support = 0.3,
                              min_confidence = 0, min_lift = 0, max_length = 2)
  rules5 <- association_rules(.make_ar_trans(), min_support = 0.3,
                              min_confidence = 0, min_lift = 0, max_length = 5)
  max_size2 <- max(vapply(rules2$rules$antecedent, length, integer(1)) +
                     vapply(rules2$rules$consequent, length, integer(1)))
  expect_true(max_size2 <= 2)
  expect_true(rules5$n_rules >= rules2$n_rules)
})


# ---- 10. Netobject input ----

test_that("association_rules works on netobject", {
  set.seed(42)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:5], 50, TRUE),
    V2 = sample(LETTERS[1:5], 50, TRUE),
    V3 = sample(LETTERS[1:5], 50, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  rules <- association_rules(net, min_support = 0.05,
                             min_confidence = 0.3, min_lift = 0)
  expect_s3_class(rules, "net_association_rules")
  expect_true(rules$n_transactions == 50)
  # Items should be network states
  expect_true(all(rules$items %in% net$nodes$label))
})


# ---- 11. Data frame input (character columns) ----

test_that("association_rules works on character data frame", {
  df <- data.frame(
    V1 = c("A", "B", "A", "C"),
    V2 = c("B", "C", "B", "A"),
    V3 = c("C", "A", "C", "B"),
    stringsAsFactors = FALSE
  )
  rules <- association_rules(df, min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  expect_s3_class(rules, "net_association_rules")
  expect_equal(rules$n_transactions, 4)
})


# ---- 12. Binary matrix input ----

test_that("association_rules works on binary matrix", {
  mat <- matrix(c(1,1,0, 1,0,1, 0,1,1, 1,1,1), nrow = 4, byrow = TRUE)
  colnames(mat) <- c("X", "Y", "Z")
  rules <- association_rules(mat, min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  expect_s3_class(rules, "net_association_rules")
  expect_equal(rules$items, c("X", "Y", "Z"))
})


# ---- 13. Empty result ----

test_that("no rules found returns empty result gracefully", {
  trans <- list(c("A"), c("B"), c("C"))
  rules <- association_rules(trans, min_support = 0.9, min_confidence = 0.9)
  expect_s3_class(rules, "net_association_rules")
  expect_equal(rules$n_rules, 0)
  expect_equal(nrow(rules$rules), 0)
})


# ---- 14. Itemsets stored as lists, not strings ----

test_that("antecedent and consequent are list columns", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  expect_true(is.list(rules$rules$antecedent))
  expect_true(is.list(rules$rules$consequent))
  # Each element is a character vector
  expect_true(all(vapply(rules$rules$antecedent, is.character, logical(1))))
  expect_true(all(vapply(rules$rules$consequent, is.character, logical(1))))
})


# ---- 15. Rules are sorted by lift (descending) ----

test_that("rules sorted by lift descending", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  lifts <- rules$rules$lift
  expect_true(all(diff(lifts) <= 1e-10))
})


# ---- 16. Pruning step works at k >= 3 ----

test_that("3-itemset rules exist when data supports them", {
  trans <- list(
    c("A", "B", "C", "D"),
    c("A", "B", "C"),
    c("A", "B", "D"),
    c("A", "C", "D"),
    c("B", "C", "D"),
    c("A", "B", "C", "D")
  )
  rules <- association_rules(trans, min_support = 0.3,
                             min_confidence = 0, min_lift = 0, max_length = 4)
  # Should have rules from 3- and 4-itemsets
  sizes <- vapply(rules$rules$antecedent, length, integer(1)) +
           vapply(rules$rules$consequent, length, integer(1))
  expect_true(any(sizes >= 3))
})


# ---- 17. tna::group_regulation ----

test_that("association_rules works on tna::group_regulation via netobject", {
  skip_if_not_installed("tna")
  data(group_regulation, package = "tna")
  net <- build_network(group_regulation, method = "relative")
  rules <- association_rules(net, min_support = 0.3,
                             min_confidence = 0.5, min_lift = 1)
  expect_s3_class(rules, "net_association_rules")
  expect_true(rules$n_rules > 0)
})


# ---- 18. print and summary work ----

test_that("print and summary methods work", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  expect_output(print(rules), "Association Rules")
  s <- summary(rules)
  expect_true(is.data.frame(s))
})


# ---- 19. Frequent itemsets stored correctly ----

test_that("frequent_itemsets has correct structure", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  fi <- rules$frequent_itemsets
  expect_true(is.list(fi))
  # Level 1: single items
  expect_true(all(vapply(fi[[1]], function(x) length(x$items) == 1, logical(1))))
  # Level 2: pairs
  if (length(fi) >= 2) {
    expect_true(all(vapply(fi[[2]], function(x) length(x$items) == 2,
                           logical(1))))
  }
})


# ---- 20. Symmetry: {A}=>{B} and {B}=>{A} both generated ----

test_that("both directions of rules are generated", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  r <- rules$rules
  has_ab <- any(vapply(seq_len(nrow(r)), function(i) {
    identical(r$antecedent[[i]], "bread") && identical(r$consequent[[i]], "milk")
  }, logical(1)))
  has_ba <- any(vapply(seq_len(nrow(r)), function(i) {
    identical(r$antecedent[[i]], "milk") && identical(r$consequent[[i]], "bread")
  }, logical(1)))
  expect_true(has_ab)
  expect_true(has_ba)
})


# ---- 21. Input validation ----

test_that("invalid parameters are rejected", {
  expect_error(association_rules(.make_ar_trans(), min_support = -0.1))
  expect_error(association_rules(.make_ar_trans(), min_support = 1.5))
  expect_error(association_rules(.make_ar_trans(), min_confidence = 2))
  expect_error(association_rules(.make_ar_trans(), max_length = 1))
})


# ---- 22. cograph_network input ----

test_that("association_rules works on cograph_network", {
  set.seed(42)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:4], 30, TRUE),
    V2 = sample(LETTERS[1:4], 30, TRUE),
    V3 = sample(LETTERS[1:4], 30, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(seqs, method = "relative")
  cg <- structure(list(
    weights = net$weights, nodes = net$nodes, edges = net$edges,
    directed = net$directed, data = net$data,
    meta = list(source = "test", tna = list(method = "relative"))
  ), class = c("cograph_network", "list"))
  rules <- association_rules(cg, min_support = 0.1,
                             min_confidence = 0.3, min_lift = 0)
  expect_s3_class(rules, "net_association_rules")
})


# ---- 23. Large dataset performance ----

test_that("runs in reasonable time on 1000 transactions", {
  set.seed(42)
  items <- LETTERS[1:10]
  trans <- lapply(seq_len(1000), function(i) {
    sample(items, sample(3:6, 1))
  })
  t1 <- system.time(
    rules <- association_rules(trans, min_support = 0.05,
                               min_confidence = 0.3, min_lift = 1)
  )
  expect_true(t1["elapsed"] < 30)
  expect_true(rules$n_rules > 0)
})


# ---- 24. Bundled data: human_cat ----

test_that("association_rules works on human_cat data", {
  data(human_cat)
  net <- build_network(human_cat, method = "relative",
                       actor = "session_id", action = "category",
                       time = "timestamp")
  rules <- association_rules(net, min_support = 0.3,
                             min_confidence = 0.5, min_lift = 1)
  expect_s3_class(rules, "net_association_rules")
})


# ---- 25. pathways.net_association_rules ----

test_that("pathways returns arrow-notation strings from rules", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  pw <- pathways(rules)
  expect_type(pw, "character")
  expect_true(length(pw) > 0)
  expect_true(all(grepl("->", pw, fixed = TRUE)))
})

test_that("pathways top parameter limits output", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  pw_all <- pathways(rules)
  pw_top <- pathways(rules, top = 3)
  expect_true(length(pw_all) >= length(pw_top))
  expect_equal(length(pw_top), min(3, length(pw_all)))
})

test_that("pathways min_lift filter works", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  pw_all <- pathways(rules)
  pw_high <- pathways(rules, min_lift = 1.5)
  expect_true(length(pw_all) >= length(pw_high))
})

test_that("pathways min_confidence filter works", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0, min_lift = 0)
  pw_all <- pathways(rules)
  pw_high <- pathways(rules, min_confidence = 0.9)
  expect_true(length(pw_all) >= length(pw_high))
})

test_that("pathways returns empty for no-rules object", {
  trans <- list(c("A"), c("B"), c("C"))
  rules <- association_rules(trans, min_support = 0.9, min_confidence = 0.9)
  pw <- pathways(rules)
  expect_equal(length(pw), 0)
})

test_that("pathways format is compatible with plot_simplicial", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  pw <- pathways(rules)
  # Each pathway: "source1 source2 -> target"
  # Split on " -> " → exactly 2 parts
  parts <- strsplit(pw, " -> ", fixed = TRUE)
  expect_true(all(vapply(parts, length, integer(1)) == 2L))
  # Sources and target are non-empty
  expect_true(all(vapply(parts, function(p) nchar(p[1]) > 0, logical(1))))
  expect_true(all(vapply(parts, function(p) nchar(p[2]) > 0, logical(1))))
})


# ---- 26. Improved print uses arrow notation ----

test_that("print uses arrow notation", {
  rules <- association_rules(.make_ar_trans(), min_support = 0.3,
                             min_confidence = 0.5, min_lift = 0)
  out <- capture.output(print(rules))
  # Should have numbered rules with "->"
  rule_lines <- out[grepl("->", out, fixed = TRUE)]
  expect_true(length(rule_lines) > 0)
  # Should have numbered format "  1. ..."
  expect_true(any(grepl("^\\s+\\d+\\.", rule_lines)))
})
