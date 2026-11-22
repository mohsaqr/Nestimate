# as_tna() and as_netobject(): conversions between tna models, cograph networks and netobjects. They carry the
# weights over, so the weights after the conversion are compared with tna's networks of the same sequences (and the
# macro layer of an mcml with tna's network of the sequences recoded to their clusters). Moved here from
# not_numeric: they return networks, whose numbers can be checked.

wide <- tn_simulate(n_actors = 40L, n_states = 6L, seq_length = 12L, seed = 501L, ragged = TRUE)
clusters <- list(g1 = c("A", "B"), g2 = c("C", "D"), g3 = c("E", "F"))

test_that("as_tna() of an mcml holds tna's network of the sequences recoded to their clusters", {
  skip_without_reference("tna")
  lookup <- stats::setNames(rep(names(clusters), lengths(clusters)), unlist(clusters, use.names = FALSE))
  recoded <- out_pad(lapply(tn_sequences(wide), \(s) unname(lookup[s])))
  converted <- as_tna(build_mcml(wide, clusters = clusters))
  expect_reference(extract_transition_matrix(converted$macro), unclass(tna::tna(recoded)$weights), reference = "tna",
                   covers = c("as_tna()", "as_tna.mcml()"), quantity = "macro layer, tna of the recoded sequences",
                   tolerance = 1e-10)
  expect_reference(extract_initial_probs(converted$macro), tna::tna(recoded)$inits, reference = "tna",
                   covers = c("as_tna()", "as_tna.mcml()"),
                   quantity = "macro layer initial probabilities, tna of the recoded sequences", tolerance = 1e-10)
})

test_that("as_tna() of a tna model and as_netobject() of tna models, cograph networks and netobjects keep tna's weights", {
  skip_without_reference(c("tna", "cograph"))
  model <- tna::tna(wide)
  expected <- unclass(model$weights)
  expect_reference(unclass(as_tna(model)$weights), expected, reference = "tna", covers = c("as_tna()", "as_tna.default()"),
                   quantity = "weights of a tna model passed through", tolerance = 0)
  expect_reference(extract_transition_matrix(as_netobject(model)), expected, reference = "tna",
                   covers = c("as_netobject()", "as_netobject.default()"), quantity = "netobject of a tna model",
                   tolerance = 0)
  expect_reference(extract_initial_probs(as_netobject(model)), model$inits, reference = "tna",
                   covers = c("as_netobject()", "as_netobject.default()"),
                   quantity = "initial probabilities of the netobject of a tna model", tolerance = 0)
  expect_reference(extract_transition_matrix(as_netobject(cograph::as_cograph(expected))), expected, reference = "tna",
                   covers = c("as_netobject()", "as_netobject.cograph_network()"),
                   quantity = "netobject of a cograph network of tna's weights", tolerance = 0)
  expect_reference(extract_transition_matrix(as_netobject(build_tna(wide))), expected, reference = "tna",
                   covers = c("as_netobject()", "as_netobject.netobject()"), quantity = "netobject passed through",
                   tolerance = 1e-10)
  grouped <- tna::group_tna(cbind(wide, grp = rep(c("a", "b"), length.out = nrow(wide))), group = "grp")
  converted <- as_netobject(grouped)
  invisible(lapply(c("a", "b"), \(g) {
    expect_reference(extract_transition_matrix(converted[[g]]), unclass(grouped[[g]]$weights), reference = "tna",
                     covers = c("as_netobject()", "as_netobject.default()"),
                     quantity = sprintf("netobject group of a group_tna, group %s", g), tolerance = 0)
  }))
})
