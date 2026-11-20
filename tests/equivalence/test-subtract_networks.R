# subtract_networks(): the edge-wise difference of two networks' weights, against the weights of tna::tna() of
# each set of sequences, for Nestimate networks and for tna models given directly. as_netdifference() of a
# difference returns it unchanged.

halves <- list(x = tna::group_regulation[seq_len(500L), ], y = tna::group_regulation[1000L + seq_len(500L), ])

test_that("the difference equals tna's weights of x minus those of y", {
  skip_without_reference("tna")
  expected <- tna::tna(halves$x)$weights - tna::tna(halves$y)$weights
  difference <- subtract_networks(build_network(halves$x, method = "relative"), build_network(halves$y, method = "relative"))
  expect_reference(extract_transition_matrix(difference), expected, reference = "tna", covers = "subtract_networks()",
                   quantity = "group_regulation halves", tolerance = 1e-12)
  expect_reference(extract_transition_matrix(as_netdifference(difference)), expected, reference = "tna",
                   covers = "as_netdifference.netdifference()", quantity = "group_regulation halves, as_netdifference()",
                   tolerance = 1e-12)
  expect_reference(extract_transition_matrix(subtract_networks(tna::tna(halves$x), tna::tna(halves$y))), expected,
                   reference = "tna", covers = "subtract_networks()", quantity = "tna models as input", tolerance = 1e-12)
})
