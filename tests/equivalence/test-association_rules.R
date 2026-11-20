# association_rules(), summary.net_association_rules(): Apriori association rules (Agrawal and Srikant 1994) with
# support, confidence, lift and conviction (Brin et al. 1997). References: arules::apriori() and
# arules::interestMeasure() for the rules with one consequent item (arules mines only those), and every rule's
# measures written from their definitions on the transactions. Re-expressed from local_testing_and_equivalence/
# test-association_rules_arules-equiv-extracted.R and test-equiv-association-rules.R (108 random configurations by
# their largest difference).

ho_basket <- list(c("bread", "milk", "eggs"), c("bread", "butter"), c("milk", "eggs", "butter"),
                  c("bread", "milk", "eggs", "butter"), c("bread", "milk"))

ho_random_baskets <- function(n, items, density, seed) {
  set.seed(seed)
  present <- matrix(stats::runif(n * items) < density, n, items)
  Filter(length, lapply(seq_len(n), \(i) paste0("i", seq_len(items))[present[i, ]]))
}

inputs <- list(list(name = "basket", transactions = ho_basket, support = 0.3, confidence = 0),
               list(name = "basket, confidence 0.5", transactions = ho_basket, support = 0.2, confidence = 0.5),
               list(name = "random 60 x 6", transactions = ho_random_baskets(60L, 6L, 0.35, 1101L), support = 0.1,
                    confidence = 0.3),
               list(name = "random 120 x 8", transactions = ho_random_baskets(120L, 8L, 0.3, 1102L), support = 0.08,
                    confidence = 0.4))

ho_items_key <- \(items) vapply(items, \(set) paste(sort(set), collapse = ", "), character(1))

# One key per rule: the helper joins several key columns with ", ", which would make "a, b => c" and
# "a => b, c" the same key.
ho_rules <- \(rules) data.frame(rule = paste(rules$antecedent, "=>", rules$consequent), subset(rules, select = -c(antecedent, consequent)))

test_that("rules with one consequent and their measures equal arules::apriori()", {
  skip_without_reference("arules")
  invisible(lapply(inputs, \(input) {
    transactions <- methods::as(input$transactions, "transactions")
    rules <- arules::apriori(transactions, parameter = list(supp = input$support, conf = input$confidence, minlen = 2L,
                                                           maxlen = 5L, target = "rules"), control = list(verbose = FALSE))
    quality <- arules::quality(rules)
    expected <- data.frame(antecedent = ho_items_key(methods::as(arules::lhs(rules), "list")),
                           consequent = ho_items_key(methods::as(arules::rhs(rules), "list")),
                           support = quality$support, confidence = quality$confidence, lift = quality$lift,
                           conviction = arules::interestMeasure(rules, "conviction", transactions = transactions))
    found <- association_rules(input$transactions, min_support = input$support, min_confidence = input$confidence,
                               min_lift = 0, max_length = 5L)
    expect_reference(ho_rules(subset(summary(found), !grepl(",", consequent))), ho_rules(expected), reference = "arules",
                     covers = c("association_rules()", "summary.net_association_rules()"),
                     quantity = sprintf("%s: rules with one consequent", input$name), tolerance = 1e-10,
                     key = "rule")
  }))
})

test_that("support, confidence, lift and conviction of every rule equal their definitions", {
  invisible(lapply(inputs, \(input) {
    found <- summary(association_rules(input$transactions, min_support = input$support,
                                       min_confidence = input$confidence, min_lift = 0, max_length = 5L))
    share <- \(items) vapply(items, \(set) mean(vapply(input$transactions, \(t) all(set %in% t), logical(1))), numeric(1))
    antecedent <- strsplit(found$antecedent, ", ", fixed = TRUE)
    consequent <- strsplit(found$consequent, ", ", fixed = TRUE)
    both <- share(Map(c, antecedent, consequent))
    confidence <- both / share(antecedent)
    expected <- data.frame(antecedent = found$antecedent, consequent = found$consequent, support = both,
                           confidence = confidence, lift = confidence / share(consequent),
                           conviction = (1 - share(consequent)) / (1 - confidence))
    expect_reference(ho_rules(found), ho_rules(expected), reference = "independent: support, confidence, lift and conviction from the transactions",
                     covers = c("association_rules()", "summary.net_association_rules()"),
                     quantity = sprintf("%s: every rule", input$name), tolerance = 1e-10,
                     key = "rule")
  }))
})
