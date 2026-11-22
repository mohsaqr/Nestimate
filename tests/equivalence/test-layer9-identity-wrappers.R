# Layer 9, composite identity: the fourteen aliases fix `method` and forward every other argument to
# build_network(). Each wrapper fit must equal the direct build_network(method = ) fit exactly (tolerance 0):
# weights, edges, initial probabilities and settings, over wide, long, one-hot, grouped and multilevel input and
# with start/end markers, thresholds, scaling, windows and estimator arguments given in `...` or in `params`.

wide <- tn_simulate(n_actors = 40L, n_states = 5L, seq_length = 10L, seed = 911L, ragged = TRUE)
long <- tn_long(wide)
codes <- c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")
onehot <- subset(learning_activities, student <= 30)
continuous <- psy_continuous(n = 200L, p = 6L, rho = 0.6, seed = 912L)
binary <- psy_binary(n = 300L, p = 5L, rho = 0.4, seed = 913L)
# a message log with times, for the latency network
messages <- subset(group_regulation_long, Course == "A")

# One case per row: the wrapper call and the direct call it must equal.
wrapper_cases <- list(
  list(covers = "build_tna()", quantity = "wide sequences",
       wrapper = \() build_tna(wide), direct = \() build_network(wide, method = "relative")),
  list(covers = "build_tna()", quantity = "start and a labelled end",
       wrapper = \() build_tna(wide, start = TRUE, end = "Fin"),
       direct = \() build_network(wide, method = "relative", start = TRUE, end = "Fin")),
  list(covers = "build_tna()", quantity = "threshold and minmax scaling",
       wrapper = \() build_tna(wide, threshold = 0.1, scaling = "minmax"),
       direct = \() build_network(wide, method = "relative", threshold = 0.1, scaling = "minmax")),
  list(covers = "build_tna()", quantity = "long event log with actor and order",
       wrapper = \() build_tna(long, actor = "actor", action = "action", order = "time"),
       direct = \() build_network(long, method = "relative", actor = "actor", action = "action", order = "time")),
  list(covers = "build_tna()", quantity = "one-hot data with actor",
       wrapper = \() build_tna(onehot, codes = codes, actor = "student"),
       direct = \() build_network(onehot, method = "relative", codes = codes, actor = "student")),
  list(covers = "build_tna()", quantity = "declared state columns",
       wrapper = \() build_tna(wide, state_cols = paste0("T", 1:6)),
       direct = \() build_network(wide, method = "relative", state_cols = paste0("T", 1:6))),
  list(covers = "build_ftna()", quantity = "sequence windows in params",
       wrapper = \() build_ftna(wide, params = list(window_size = 2L)),
       direct = \() build_network(wide, method = "frequency", params = list(window_size = 2L))),
  list(covers = "build_ftna()", quantity = "max scaling, then threshold 0.2",
       wrapper = \() build_ftna(wide, threshold = 0.2, scaling = "max"),
       direct = \() build_network(wide, method = "frequency", threshold = 0.2, scaling = "max")),
  list(covers = "build_ftna()", quantity = "one-hot data, window_size 2 in `...`",
       wrapper = \() build_ftna(onehot, codes = codes, actor = "student", window_size = 2L),
       direct = \() build_network(onehot, method = "frequency", codes = codes, actor = "student", window_size = 2L)),
  list(covers = "build_atna()", quantity = "lambda in `...`",
       wrapper = \() build_atna(wide, lambda = 2), direct = \() build_network(wide, method = "attention", lambda = 2)),
  list(covers = "build_atna()", quantity = "lambda in params, with start",
       wrapper = \() build_atna(wide, start = TRUE, params = list(lambda = 3)),
       direct = \() build_network(wide, method = "attention", start = TRUE, params = list(lambda = 3))),
  list(covers = "build_cna()", quantity = "start marker",
       wrapper = \() build_cna(wide, start = TRUE), direct = \() build_network(wide, method = "co_occurrence", start = TRUE)),
  list(covers = "build_cna()", quantity = "one-hot data, window_size 2",
       wrapper = \() build_cna(onehot, codes = codes, actor = "student", window_size = 2L),
       direct = \() build_network(onehot, method = "co_occurrence", codes = codes, actor = "student", window_size = 2L)),
  list(covers = "build_cor()", quantity = "continuous variables",
       wrapper = \() build_cor(continuous), direct = \() build_network(continuous, method = "cor")),
  list(covers = "build_pcor()", quantity = "threshold 0.05",
       wrapper = \() build_pcor(continuous, threshold = 0.05),
       direct = \() build_network(continuous, method = "pcor", threshold = 0.05)),
  list(covers = "build_glasso()", quantity = "gamma in `...`",
       wrapper = \() build_glasso(continuous, gamma = 0.25), direct = \() build_network(continuous, method = "glasso", gamma = 0.25)),
  list(covers = "build_glasso()", quantity = "gamma in params",
       wrapper = \() build_glasso(continuous, params = list(gamma = 0.25)),
       direct = \() build_network(continuous, method = "glasso", params = list(gamma = 0.25))),
  list(covers = "build_ising()", quantity = "binary variables",
       wrapper = \() build_ising(binary), direct = \() build_network(binary, method = "ising")),
  list(covers = "build_ising()", quantity = "gamma 0.1 in `...`",
       wrapper = \() build_ising(binary, gamma = 0.1), direct = \() build_network(binary, method = "ising", gamma = 0.1)),
  # the method settings named on the alias go to the estimator as build_network()'s `...` and `params` send them
  list(covers = "build_atna()", quantity = "lambda and direction named on the alias",
       wrapper = \() build_atna(wide, lambda = 2, direction = "both"),
       direct = \() build_network(wide, method = "attention", params = list(lambda = 2, direction = "both"))),
  list(covers = "build_cna()", quantity = "type named on the alias",
       wrapper = \() build_cna(wide, type = "relative"),
       direct = \() build_network(wide, method = "co_occurrence", params = list(type = "relative"))),
  list(covers = "build_cor()", quantity = "cor_method named on the alias",
       wrapper = \() build_cor(continuous, cor_method = "spearman"),
       direct = \() build_network(continuous, method = "cor", params = list(cor_method = "spearman"))),
  list(covers = "build_glasso()", quantity = "gamma and nlambda named on the alias",
       wrapper = \() build_glasso(continuous, gamma = 0.25, nlambda = 50L),
       direct = \() build_network(continuous, method = "glasso", params = list(gamma = 0.25, nlambda = 50L))),
  list(covers = "build_ising()", quantity = "rule named on the alias",
       wrapper = \() build_ising(binary, rule = "OR"),
       direct = \() build_network(binary, method = "ising", params = list(rule = "OR"))),
  list(covers = "build_ngram()", quantity = "n_gram 3",
       wrapper = \() build_ngram(wide, n_gram = 3L), direct = \() build_network(wide, method = "ngram", params = list(n_gram = 3L))),
  list(covers = "build_ngram()", quantity = "start and end",
       wrapper = \() build_ngram(wide, start = TRUE, end = TRUE),
       direct = \() build_network(wide, method = "ngram", start = TRUE, end = TRUE)),
  list(covers = "build_gap()", quantity = "max_gap 2",
       wrapper = \() build_gap(wide, max_gap = 2L), direct = \() build_network(wide, method = "gap", params = list(max_gap = 2L))),
  list(covers = "build_gap()", quantity = "long event log",
       wrapper = \() build_gap(long, actor = "actor", action = "action", order = "time"),
       direct = \() build_network(long, method = "gap", actor = "actor", action = "action", order = "time")),
  list(covers = "build_reverse()", quantity = "weighted by sequence length",
       wrapper = \() build_reverse(wide, weighted = TRUE),
       direct = \() build_network(wide, method = "reverse", params = list(weighted = TRUE))),
  list(covers = "build_reverse()", quantity = "threshold 1",
       wrapper = \() build_reverse(wide, threshold = 1), direct = \() build_network(wide, method = "reverse", threshold = 1)),
  list(covers = "build_latency()", quantity = "floor and cap named on the alias",
       wrapper = \() build_latency(messages, action = "Action", time = "Time", actor = "Actor", floor = 0.3, cap = "60s"),
       direct = \() build_network(messages, method = "latency", action = "Action", time = "Time", actor = "Actor",
                                  latency_floor = 0.3, latency_cap = "60s")),
  list(covers = "build_latency()", quantity = "mean aggregate",
       wrapper = \() build_latency(messages, action = "Action", time = "Time", actor = "Actor", aggregate = "mean"),
       direct = \() build_network(messages, method = "latency", action = "Action", time = "Time", actor = "Actor",
                                  latency_aggregate = "mean")),
  list(covers = "build_wtna()", quantity = "one-hot data with actor",
       wrapper = \() build_wtna(onehot, codes = codes, actor = "student"),
       direct = \() build_network(onehot, method = "wtna", codes = codes, actor = "student")),
  list(covers = "build_wtna()", quantity = "relative type, window_size 2",
       wrapper = \() build_wtna(onehot, codes = codes, actor = "student", window_size = 2L, type = "relative"),
       direct = \() build_network(onehot, method = "wtna", codes = codes, actor = "student", window_size = 2L,
                                  params = list(type = "relative"))),
  list(covers = "build_mgm()", quantity = "continuous variables",
       wrapper = \() build_mgm(continuous), direct = \() build_network(continuous, method = "mgm")),
  list(covers = "build_mgm()", quantity = "ruleReg OR named on the alias",
       wrapper = \() build_mgm(continuous, ruleReg = "OR"),
       direct = \() build_network(continuous, method = "mgm", params = list(ruleReg = "OR")))
)

test_that("every wrapper fit equals the direct build_network() fit", {
  lapply(wrapper_cases, \(case) {
    li_expect_same_network(case$wrapper(), case$direct(), route = "build_network(method = ) on the same input",
                           covers = c(case$covers, "build_network()"), quantity = case$quantity)
  })
})

test_that("the association wrappers forward actor and level to the multilevel decomposition", {
  panel <- transform(continuous, id = rep(seq_len(20L), each = 10L))
  between_cor <- build_cor(panel, actor = "id", level = "between")
  li_expect_same_network(between_cor, build_network(panel, method = "cor", actor = "id", level = "between"),
                         route = "build_network(method = 'cor', level = 'between')", covers = c("build_cor()", "build_network()"),
                         quantity = "between-person correlations")
  within_glasso <- build_glasso(panel, actor = "id", level = "within")
  li_expect_same_network(within_glasso, build_network(panel, method = "glasso", actor = "id", level = "within"),
                         route = "build_network(method = 'glasso', level = 'within')", covers = c("build_glasso()", "build_network()"),
                         quantity = "within-person glasso")
})

test_that("the wrappers forward `group` and every group network equals the direct grouped fit", {
  grouped <- transform(wide, g = rep(c("x", "y"), 20L))
  via_wrapper <- build_tna(grouped, group = "g", start = TRUE)
  direct <- build_network(grouped, method = "relative", group = "g", start = TRUE)
  expect_reference(names(via_wrapper), names(direct), reference = li_ref("build_network(group = )"),
                   covers = c("build_tna()", "build_network()"), quantity = "group names and order", tolerance = 0,
                   authority = "internal")
  Map(\(a, e, name) li_expect_same_network(a, e, route = "build_network(group = )", covers = c("build_tna()", "build_network()"),
                                         quantity = sprintf("group %s", name)),
      via_wrapper, direct, names(direct))
  grouped_cor <- transform(continuous, g = rep(c("a", "b"), 100L))
  Map(\(a, e, name) li_expect_same_network(a, e, route = "build_network(group = ) association", covers = c("build_cor()", "build_network()"),
                                         quantity = sprintf("correlation group %s", name)),
      build_cor(grouped_cor, group = "g"), build_network(grouped_cor, method = "cor", group = "g"), c("a", "b"))
})
