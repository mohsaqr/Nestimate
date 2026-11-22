# Part 2 of test-layer5-invariance-networks.R, split so that one Prova case holds at most 7 references;
# the shared set-up is repeated.

# Layer 5, invariance under transformation (metamorphic tests): a network is refitted after a transformation of
# its input whose effect on every output is known in advance, and the refit must equal the original carried
# through the same transformation. Rows permuted, actors' blocks of a long log shuffled, states, actors and groups
# renamed, factor levels reversed: nothing changes. Every sequence duplicated: relative weights and initial
# probabilities unchanged, counts doubled. Association data with every variable multiplied by a positive constant:
# correlations, partial correlations and the graphical lasso unchanged.
# Transition counts are integers, so those checks are exact (tolerance 0). A correlation of permuted or rescaled
# rows sums the same terms in another order; those checks allow 1e-12.

wide <- tn_simulate(n_actors = 50L, n_states = 5L, seq_length = 10L, seed = 951L, ragged = TRUE)
set.seed(952L)
row_order <- sample(nrow(wide))
state_map <- c(A = "zeta", B = "alpha", C = "mu", D = "beta", E = "omega")
renamed <- li_rename_states(wide, state_map)
back <- li_invert(state_map)
doubled <- rbind(wide, wide)
sequence_methods <- c("relative", "frequency", "co_occurrence", "attention", "ngram", "gap", "reverse")
count_methods <- setdiff(sequence_methods, "relative")

test_that("net_centrality() is invariant to row order and follows renamed states", {
  measures <- c("InStrength", "OutStrength", "Betweenness", "Closeness", "Diffusion")
  quiet <- function(net) {
    withCallingHandlers(net_centrality(net, measures = measures), message = \(condition) invokeRestart("muffleMessage"))
  }
  as_table <- function(central, map = NULL) {
    state <- as.character(central$state)
    if (!is.null(map)) state <- unname(map[state])
    do.call(rbind, lapply(measures, \(measure) data.frame(state = state, measure = measure, value = central[[measure]])))
  }
  lapply(c("relative", "frequency", "co_occurrence"), \(method) {
    original <- as_table(quiet(build_network(wide, method = method)))
    expect_reference(as_table(quiet(build_network(wide[row_order, ], method = method))), original,
                     reference = li_ref("centrality of the fit on the original rows"), covers = "net_centrality()",
                     quantity = sprintf("%s centralities, rows permuted", method), key = c("state", "measure"),
                     tolerance = 0, authority = "internal")
    expect_reference(as_table(quiet(build_network(renamed, method = method)), back), original,
                     reference = li_ref("centrality of the fit on the original state names"), covers = "net_centrality()",
                     quantity = sprintf("%s centralities, states renamed", method), key = c("state", "measure"),
                     tolerance = 1e-12, authority = "internal",
                     note = "Renaming reorders the nodes, so shortest paths are relaxed in another order; 1e-12 allows for it.")
  })
})

test_that("association networks are invariant to row order, positive rescaling and duplication", {
  continuous <- psy_continuous(n = 200L, p = 6L, rho = 0.8, seed = 954L)
  rescaled <- continuous
  rescaled[] <- Map(\(column, factor) column * factor, continuous, c(1e-5, 3, 1e3, 0.2, 1e6, 7))
  set.seed(955L)
  permuted <- continuous[sample(nrow(continuous)), ]
  reordered <- continuous[, c("V6", "V2", "V4", "V1", "V3", "V5")]
  # The refit sums the same terms in another order: 1e-12 allows for it. The graphical lasso is solved by
  # coordinate descent over the variables in their order, which converges to that order's last digits: 1e-10.
  lapply(c("cor", "pcor", "glasso"), \(method) {
    original <- build_network(continuous, method = method)
    cases <- list(rows_permuted = permuted, variables_rescaled = rescaled, columns_reordered = reordered)
    Map(\(data, label) {
      li_expect_same_network(build_network(data, method = method), original, route = "fit on the original variables",
                             covers = "build_network()", quantity = sprintf("%s, %s", method, gsub("_", " ", label)),
                             parts = c("weights", "edges"), tolerance = if (method == "glasso") 1e-10 else 1e-12)
    }, cases, names(cases))
  })
  # glasso's EBIC depends on n, so duplication changes it by design; correlations do not
  lapply(c("cor", "pcor"), \(method) {
    li_expect_same_network(build_network(rbind(continuous, continuous), method = method),
                           build_network(continuous, method = method), route = "fit on the rows once",
                           covers = "build_network()", quantity = sprintf("%s, every row twice", method),
                           parts = c("weights", "edges"), tolerance = 1e-12)
  })
  binary <- psy_binary(n = 300L, p = 5L, rho = 0.4, seed = 956L)
  set.seed(957L)
  original <- build_network(binary, method = "ising")
  li_expect_same_network(build_network(binary[sample(nrow(binary)), ], method = "ising"), original,
                         route = "fit on the original rows", covers = "build_network()", quantity = "ising, rows permuted",
                         parts = c("weights", "edges"), tolerance = 1e-12)
  expect_reference(build_network(binary[, c("V3", "V1", "V5", "V2", "V4")], method = "ising")$weights, original$weights,
                   reference = li_ref("fit on the variables in their original order"), covers = "build_network()",
                   quantity = "ising weights, columns reordered", tolerance = 1e-5, authority = "internal",
                   note = "Each node's lasso path is solved by coordinate descent (glmnet), whose result depends on the order of the predictors within its convergence threshold; 1e-5 is that precision.")
})
