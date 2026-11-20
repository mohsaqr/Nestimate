# chain_structure(): the structure of a Markov chain from its transition matrix (communicating, recurrent and
# absorbing classes, hitting probabilities, absorption probabilities and mean absorption times), against
# markovchain 1.2 on three regimes of random matrices: irreducible, two absorbing states, and two closed classes.
# Re-expressed from local_testing_and_equivalence/test-equiv-chain_structure.R (a sweep of 100 matrices that
# compared maxima of differences); here the whole table of each matrix is compared number by number.

# A row-stochastic matrix with strictly positive entries, states s1..sn.
ref_stochastic <- function(n_states, seed) {
  set.seed(seed)
  weights <- matrix(stats::runif(n_states^2, 0.05, 1), n_states, n_states)
  states <- paste0("s", seq_len(n_states))
  matrix(weights / rowSums(weights), n_states, n_states, dimnames = list(states, states))
}

# The same matrix with the named states made absorbing.
ref_absorbing <- function(n_states, absorbing, seed) {
  matrix_p <- ref_stochastic(n_states, seed)
  matrix_p[absorbing, ] <- 0
  matrix_p[cbind(absorbing, absorbing)] <- 1
  matrix_p
}

# Two closed classes of sizes k1 and k2 that never communicate.
ref_two_classes <- function(k1, k2, seed) {
  n_states <- k1 + k2
  states <- paste0("s", seq_len(n_states))
  matrix_p <- matrix(0, n_states, n_states, dimnames = list(states, states))
  matrix_p[seq_len(k1), seq_len(k1)] <- ref_stochastic(k1, seed)
  matrix_p[k1 + seq_len(k2), k1 + seq_len(k2)] <- ref_stochastic(k2, seed + 1L)
  matrix_p
}

# markovchain's view of a transition matrix.
ref_chain <- function(matrix_p) {
  methods::new(methods::getClassDef("markovchain", package = "markovchain"),
               states = rownames(matrix_p), byrow = TRUE, transitionMatrix = matrix_p)
}

# One row per state: markovchain's classification (absorbing / recurrent / transient).
ref_states <- function(matrix_p) {
  chain <- ref_chain(matrix_p)
  states <- rownames(matrix_p)
  absorbing <- markovchain::absorbingStates(chain)
  transient <- markovchain::transientStates(chain)
  data.frame(state = states,
             classification = ifelse(states %in% absorbing, "absorbing", ifelse(states %in% transient, "transient", "recurrent")))
}

# summary() rounds its numbers for display (to 3 or 4 digits), so the unrounded hitting probabilities, absorption
# probabilities and absorption times are read from the result's fields: the package has no unrounded tidy
# accessor for them (API gap).

test_that("irreducible chains: classification and return probabilities equal markovchain", {
  skip_without_reference("markovchain")
  inputs <- list(ref_stochastic(3L, 11L), ref_stochastic(4L, 12L), ref_stochastic(5L, 13L))
  invisible(lapply(seq_along(inputs), \(i) {
    structure <- chain_structure(inputs[[i]])
    expect_reference(summary(structure), ref_states(inputs[[i]]), reference = "markovchain",
                     covers = c("chain_structure()", "summary.chain_structure()"),
                     quantity = sprintf("irreducible chain %d, classification", i), tolerance = 0, key = "state")
    expect_reference(structure$hitting_probabilities,
                     as.matrix(markovchain::hittingProbabilities(ref_chain(inputs[[i]]))), reference = "markovchain",
                     covers = "chain_structure()", quantity = sprintf("irreducible chain %d, hitting probabilities", i),
                     tolerance = 1e-10)
  }))
})

test_that("irreducible chains: period equals markovchain::period()", {
  skip_without_reference("markovchain")
  matrix_p <- ref_stochastic(4L, 12L)
  # A two-cycle has period 2: the case where the period is not trivially 1.
  cycle <- matrix(c(0, 1, 1, 0), 2L, 2L, dimnames = list(c("a", "b"), c("a", "b")))
  expect_reference(summary(chain_structure(matrix_p)),
                   data.frame(state = rownames(matrix_p), period = markovchain::period(ref_chain(matrix_p))),
                   reference = "markovchain", covers = c("chain_structure()", "summary.chain_structure()"),
                   quantity = "period, positive chain", tolerance = 0, key = "state")
  expect_reference(summary(chain_structure(cycle)),
                   data.frame(state = c("a", "b"), period = markovchain::period(ref_chain(cycle))),
                   reference = "markovchain", covers = c("chain_structure()", "summary.chain_structure()"),
                   quantity = "period, two-cycle", tolerance = 0, key = "state")
})

test_that("absorbing chains: absorption probabilities and mean absorption times equal markovchain", {
  skip_without_reference("markovchain")
  inputs <- list(ref_absorbing(5L, c(2L, 5L), 21L), ref_absorbing(4L, 4L, 22L), ref_absorbing(6L, c(1L, 3L), 23L))
  invisible(lapply(seq_along(inputs), \(i) {
    matrix_p <- inputs[[i]]
    chain <- ref_chain(matrix_p)
    structure <- chain_structure(matrix_p)
    expect_reference(structure$absorption_probabilities, as.matrix(markovchain::absorptionProbabilities(chain)),
                     reference = "markovchain", covers = "chain_structure()",
                     quantity = sprintf("absorbing chain %d, absorption probabilities", i), tolerance = 1e-10)
    expect_reference(structure$mean_absorption_time, markovchain::meanAbsorptionTime(chain),
                     reference = "markovchain", covers = "chain_structure()",
                     quantity = sprintf("absorbing chain %d, mean absorption time", i), tolerance = 1e-10)
    expect_reference(structure$hitting_probabilities, as.matrix(markovchain::hittingProbabilities(chain)),
                     reference = "markovchain", covers = "chain_structure()",
                     quantity = sprintf("absorbing chain %d, hitting probabilities", i), tolerance = 1e-10)
    expect_reference(summary(structure), ref_states(matrix_p), reference = "markovchain",
                     covers = c("chain_structure()", "summary.chain_structure()"),
                     quantity = sprintf("absorbing chain %d, classification", i), tolerance = 0, key = "state")
  }))
})

test_that("two closed classes: classification and return probabilities equal markovchain", {
  skip_without_reference("markovchain")
  inputs <- list(ref_two_classes(2L, 3L, 31L), ref_two_classes(3L, 3L, 32L))
  invisible(lapply(seq_along(inputs), \(i) {
    structure <- chain_structure(inputs[[i]])
    expect_reference(summary(structure), ref_states(inputs[[i]]), reference = "markovchain",
                     covers = c("chain_structure()", "summary.chain_structure()"),
                     quantity = sprintf("two closed classes %d, classification", i), tolerance = 0, key = "state")
    expect_reference(structure$hitting_probabilities,
                     as.matrix(markovchain::hittingProbabilities(ref_chain(inputs[[i]]))), reference = "markovchain",
                     covers = "chain_structure()", quantity = sprintf("two closed classes %d, hitting probabilities", i),
                     tolerance = 1e-10)
  }))
})
