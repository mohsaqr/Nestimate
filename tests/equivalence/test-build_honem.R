# build_honem(): HONEM, higher-order network embedding (Saebi, Ciampaglia, Kaplan and Chawla 2020). The proximity
# matrix S = sum_{k=0..L} exp(-k) D^k / sum exp(-k), D the row-normalised higher-order adjacency, is factorised by
# a truncated SVD and the embedding is U sqrt(Sigma). References: the singular values of S from base::svd() written
# from that definition, and from RSpectra::svds(); the embedding compared through E E' = U Sigma U', which does not
# depend on the signs of the singular vectors. Re-expressed from local_testing_and_equivalence/test-equiv-honem.R
# (40 random inputs by their largest difference). The HON adjacency, the embedding matrix and the singular values
# are read from the results' fields (API gap).

ho_honem_proximity <- function(adjacency, max_power) {
  strength <- rowSums(adjacency)
  transition <- adjacency / ifelse(strength > 0, strength, 1)
  weights <- exp(-(0:max_power))
  powers <- Reduce(\(power, k) power %*% transition, seq_len(max_power), accumulate = TRUE, init = transition)
  Reduce(`+`, Map(`*`, weights, powers)) / sum(weights)
}

inputs <- list(list(seed = 801L, states = 4L, max_order = 2L, dim = 3L, max_power = 5L),
               list(seed = 802L, states = 5L, max_order = 1L, dim = 4L, max_power = 10L),
               list(seed = 803L, states = 4L, max_order = 2L, dim = 2L, max_power = 3L))

honem_case <- function(input) {
  wide <- ho_sequences(20L, 30L, ho_stochastic(input$states, input$seed), input$seed)
  hon <- build_hon(lapply(seq_len(nrow(wide)), \(i) unlist(wide[i, ], use.names = FALSE)), max_order = input$max_order,
                   min_freq = 1L, method = "hon")
  proximity <- ho_honem_proximity(hon$matrix, input$max_power)
  list(found = build_honem(hon, dim = input$dim, max_power = input$max_power), proximity = proximity,
       svd = svd(proximity))
}

test_that("singular values and E E' equal the SVD of the HONEM proximity matrix", {
  invisible(lapply(inputs, \(input) {
    case <- honem_case(input)
    dim <- input$dim
    sigma <- case$svd$d[seq_len(dim)]
    vectors <- case$svd$u[, seq_len(dim), drop = FALSE]
    gram <- vectors %*% diag(sigma, dim) %*% t(vectors)
    dimnames(gram) <- dimnames(case$proximity)
    label <- sprintf("seed %d, order %d, dim %d", input$seed, input$max_order, dim)
    expect_reference(case$found$singular_values, sigma,
                     reference = "independent: HONEM proximity matrix (Saebi et al. 2020) and base::svd()",
                     covers = "build_honem()", quantity = sprintf("singular values, %s", label), tolerance = 1e-10)
    expect_reference(tcrossprod(case$found$embeddings), gram,
                     reference = "independent: HONEM proximity matrix (Saebi et al. 2020) and base::svd()",
                     covers = "build_honem()", quantity = sprintf("embedding Gram matrix, %s", label), tolerance = 1e-8)
  }))
})

test_that("singular values equal RSpectra::svds() of the proximity matrix", {
  skip_without_reference("RSpectra")
  invisible(lapply(inputs, \(input) {
    case <- honem_case(input)
    expect_reference(case$found$singular_values, RSpectra::svds(case$proximity, k = input$dim)$d, reference = "RSpectra",
                     covers = "build_honem()",
                     quantity = sprintf("singular values, seed %d, dim %d", input$seed, input$dim), tolerance = 1e-8)
  }))
})

test_that("the embedding table of summary() gives the same Gram matrix U Sigma U'", {
  input <- inputs[[2L]]
  case <- honem_case(input)
  table <- summary(case$found)
  embedding <- as.matrix(subset(table, select = -node))
  rownames(embedding) <- table$node
  vectors <- case$svd$u[, seq_len(input$dim), drop = FALSE]
  gram <- vectors %*% diag(case$svd$d[seq_len(input$dim)], input$dim) %*% t(vectors)
  dimnames(gram) <- dimnames(case$proximity)
  expect_reference(tcrossprod(embedding), gram,
                   reference = "independent: HONEM proximity matrix (Saebi et al. 2020) and base::svd()",
                   covers = c("build_honem()", "summary.net_honem()"), quantity = "embedding Gram matrix from summary()",
                   tolerance = 1e-8)
})

test_that("the embedding table of summary() gives the Gram matrix of RSpectra::svds()", {
  skip_without_reference("RSpectra")
  invisible(lapply(inputs, \(input) {
    case <- honem_case(input)
    table <- summary(case$found)
    embedding <- as.matrix(subset(table, select = -node))
    rownames(embedding) <- table$node
    decomposition <- RSpectra::svds(case$proximity, k = input$dim)
    gram <- decomposition$u %*% diag(decomposition$d, input$dim) %*% t(decomposition$u)
    dimnames(gram) <- dimnames(case$proximity)
    expect_reference(tcrossprod(embedding), gram, reference = "RSpectra", covers = "summary.net_honem()",
                     quantity = sprintf("embedding Gram matrix from summary(), seed %d, RSpectra", input$seed),
                     tolerance = 1e-8)
  }))
})
