# Layer 8, the argument sweep of build_network() on the association methods: cor, pcor and glasso at every
# correlation method, every scaling and threshold, the glasso settings (gamma, nlambda, lambda.min.ratio,
# penalize.diagonal, refit), a correlation or covariance matrix with `n` as input, the multilevel decomposition
# (level between / within / both), the threshold given in `params` instead of the formal, predictability on and
# off, and the Ising rule and gamma. References: stats::cor() and solve() (closed forms), glasso::glasso() at the
# selected lambda, IsingFit::IsingFit(), and the documented scaling written in base R (helper-layer-sweep.R).
# test-build_cor.R, test-build_glasso.R and test-build_ising.R check the defaults; this file sweeps the options.

as_data <- psy_continuous(n = 160L, p = 6L, rho = 1, seed = 811L)
sw_cor <- function(data, method = "pearson") {
  r <- stats::cor(as.matrix(data), method = method)
  diag(r) <- 0
  r
}
sw_pcor <- function(data, method = "pearson") {
  precision <- solve(stats::cor(as.matrix(data), method = method))
  partial <- -precision / sqrt(outer(diag(precision), diag(precision)))
  diag(partial) <- 0
  # solve() leaves the inverse asymmetric in the last bits; a partial correlation is symmetric.
  (partial + t(partial)) / 2
}

test_that("cor and pcor at every correlation method, scaling and threshold equal their closed forms", {
  scalings <- list(NULL, "minmax", "max", "rank", "normalize", c("rank", "minmax"), c("max", "normalize"))
  cases <- expand.grid(method = c("cor", "pcor"), cor_method = c("pearson", "spearman", "kendall"),
                       scaling = seq_along(scalings), threshold = c(0, 0.08, 0.25), stringsAsFactors = FALSE)
  invisible(Map(\(method, cor_method, scaling_index, threshold) {
    scaling <- scalings[[scaling_index]]
    raw <- if (identical(method, "cor")) sw_cor(as_data, cor_method) else sw_pcor(as_data, cor_method)
    expected <- sw_threshold(sw_scale(raw, scaling %||% character(), include_zeros = FALSE), threshold)
    fit <- build_network(as_data, method = method, params = list(cor_method = cor_method), scaling = scaling,
                         threshold = threshold)
    expect_reference(extract_transition_matrix(fit), expected,
                     reference = "independent: stats::cor() / -D solve(cor) D, then the documented scaling and threshold",
                     covers = "build_network()",
                     quantity = sprintf("%s, %s, scaling %s, threshold %s", method, cor_method,
                                        paste(scaling %||% "none", collapse = "+"), threshold),
                     tolerance = 1e-10)
  }, cases$method, cases$cor_method, cases$scaling, cases$threshold))
})

test_that("a threshold in params equals the threshold formal for cor, pcor and glasso", {
  invisible(lapply(c("cor", "pcor", "glasso"), \(method) {
    lapply(c(0.05, 0.2), \(threshold) {
      expect_reference(extract_transition_matrix(build_network(as_data, method = method, params = list(threshold = threshold))),
                       extract_transition_matrix(build_network(as_data, method = method, threshold = threshold)),
                       reference = "independent: ?build_network, weights below the threshold set to zero",
                       covers = "build_network()", quantity = sprintf("%s, threshold %s in params", method, threshold),
                       tolerance = 1e-12)
    })
  }))
})

test_that("a correlation or covariance matrix with n gives the network of the data it came from", {
  r <- stats::cor(as_data)
  s <- stats::cov(as_data) * 4
  cases <- expand.grid(method = c("cor", "pcor", "glasso"), input = c("cor", "cov auto", "cov declared"),
                       stringsAsFactors = FALSE)
  invisible(Map(\(method, input) {
    fit <- switch(input,
      cor = build_network(r, method = method, params = list(n = nrow(as_data))),
      `cov auto` = build_network(s, method = method, params = list(n = nrow(as_data))),
      `cov declared` = build_network(s, method = method, params = list(n = nrow(as_data), input_type = "cov")))
    expect_reference(extract_transition_matrix(fit), extract_transition_matrix(build_network(as_data, method = method)),
                     reference = "independent: the network of the raw data the matrix summarises",
                     covers = "build_network()", quantity = sprintf("%s from a %s matrix", method, input),
                     tolerance = 1e-10)
  }, cases$method, cases$input))
  expect_error(build_network(r, method = "cor"), class = "nestimate_bad_input")
  expect_error(build_network(r + upper.tri(r) * 0.1, method = "cor", params = list(n = 100L)),
               class = "nestimate_bad_input")
})

test_that("every glasso setting gives the glasso::glasso() solution at its selected lambda, on its lambda path", {
  skip_without_reference("glasso")
  settings <- list(list(gamma = 0), list(gamma = 0.25), list(gamma = 1), list(nlambda = 30L),
                   list(lambda.min.ratio = 0.1), list(penalize.diagonal = TRUE),
                   list(cor_method = "spearman"), list(gamma = 0.25, nlambda = 50L, lambda.min.ratio = 0.05))
  invisible(lapply(settings, \(params) {
    fit <- build_network(as_data, method = "glasso", params = params)
    s <- stats::cor(as_data, method = params$cor_method %||% "pearson")
    ratio <- params$lambda.min.ratio %||% 0.01
    lambda_max <- max(abs(s[upper.tri(s)]))
    path <- exp(seq(log(lambda_max), log(lambda_max * ratio), length.out = params$nlambda %||% 100L))
    label <- paste(names(params), unlist(params), sep = " ", collapse = ", ")
    expect_reference(sort(fit$lambda_path), sort(path),
                     reference = "independent: log-spaced path from max |r| to max |r| * lambda.min.ratio",
                     covers = "build_network()", quantity = sprintf("lambda path, %s", label), tolerance = 1e-12,
                     measure = "relative")
    reference <- glasso::glasso(s, rho = fit$lambda_selected, penalize.diagonal = isTRUE(params$penalize.diagonal),
                                thr = 1e-12, maxit = 1e5)
    expected <- psy_partial(reference$wi)
    dimnames(expected) <- dimnames(s)
    expect_reference(extract_transition_matrix(fit), expected, reference = "glasso", covers = "build_network()",
                     quantity = sprintf("partial correlations at the selected lambda, %s", label), tolerance = 1e-8)
  }))
})

test_that("glasso with refit = TRUE is the unpenalised fit of the selected edges", {
  skip_without_reference("glasso")
  selected <- extract_transition_matrix(build_network(as_data, method = "glasso"))
  absent <- which(selected == 0 & upper.tri(selected), arr.ind = TRUE)
  s <- stats::cor(as_data)
  # glasso warns that rho = 0 can converge poorly on a matrix not of full rank; this correlation matrix is of full
  # rank (160 rows, 6 variables), so that one warning is muffled.
  reference <- withCallingHandlers(
    glasso::glasso(s, rho = 0, zero = if (nrow(absent)) absent else NULL, thr = 1e-12, maxit = 1e5),
    warning = \(w) if (grepl("rho=0", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
  expected <- psy_partial(reference$wi)
  dimnames(expected) <- dimnames(s)
  expect_reference(extract_transition_matrix(build_network(as_data, method = "glasso", params = list(refit = TRUE))),
                   expected, reference = "glasso", covers = "build_network()",
                   quantity = "refit = TRUE: unpenalised fit on the selected zero pattern", tolerance = 1e-8)
})

test_that("level between, within and both equal the correlations of person means and of deviations from them", {
  set.seed(812)
  ids <- rep(sprintf("p%02d", seq_len(25L)), times = sample(2:6, 25L, replace = TRUE))
  ml <- psy_continuous(n = length(ids), p = 4L, rho = 1, seed = 813L)
  person <- psy_continuous(n = 25L, p = 4L, rho = 0.5, seed = 814L)
  ml[] <- Map(\(column, shift) column + shift[match(ids, sprintf("p%02d", seq_len(25L)))], ml, person)
  ml <- cbind(id = ids, ml, stringsAsFactors = FALSE)
  variables <- setdiff(names(ml), "id")
  means <- stats::aggregate(ml[variables], by = list(id = ml$id), FUN = mean)
  deviations <- ml[variables] - means[match(ml$id, means$id), variables]
  references <- list(cor = sw_cor, pcor = sw_pcor)
  invisible(Map(\(method, ref) {
    between <- build_network(ml, method = method, actor = "id", level = "between")
    within_net <- build_network(ml, method = method, actor = "id", level = "within")
    both <- build_network(ml, method = method, actor = "id", level = "both")
    expect_reference(extract_transition_matrix(between), ref(means[variables]),
                     reference = "independent: correlations of the person means, base R", covers = "build_network()",
                     quantity = sprintf("%s, level between", method), tolerance = 1e-10)
    expect_reference(extract_transition_matrix(within_net), ref(deviations),
                     reference = "independent: correlations of the deviations from the person means, base R",
                     covers = "build_network()", quantity = sprintf("%s, level within", method), tolerance = 1e-10)
    # level = "both" holds the two (a netobject_ml, read with `$`: no accessor returns one level).
    expect_reference(extract_transition_matrix(both$between), ref(means[variables]),
                     reference = "independent: correlations of the person means, base R", covers = "build_network()",
                     quantity = sprintf("%s, level both, between", method), tolerance = 1e-10)
    expect_reference(extract_transition_matrix(both$within), ref(deviations),
                     reference = "independent: correlations of the deviations from the person means, base R",
                     covers = "build_network()", quantity = sprintf("%s, level both, within", method), tolerance = 1e-10)
    # The id given in params instead of actor is the same decomposition.
    expect_reference(extract_transition_matrix(build_network(ml, method = method, params = list(id_col = "id"),
                                                             level = "between")),
                     ref(means[variables]), reference = "independent: correlations of the person means, base R",
                     covers = "build_network()", quantity = sprintf("%s, level between, id in params", method),
                     tolerance = 1e-10)
  }, names(references), references))
  expect_error(build_network(ml[variables], method = "cor", level = "between"), class = "nestimate_bad_input")
})

test_that("predictability is stored when asked for and absent when not", {
  invisible(lapply(c("cor", "pcor", "glasso"), \(method) {
    with_r2 <- build_network(as_data, method = method)
    without <- build_network(as_data, method = method, predictability = FALSE)
    table <- predictability(with_r2)
    expect_reference(unname(with_r2$predictability[table$node]), table$R2,
                     reference = "independent: predictability() of the same network", covers = "build_network()",
                     quantity = sprintf("stored predictability, %s", method), tolerance = 1e-12)
    expect_reference(as.numeric(is.null(without$predictability)), 1,
                     reference = "independent: ?build_network, predictability = FALSE stores none",
                     covers = "build_network()", quantity = sprintf("no predictability, %s", method), tolerance = 0)
    expect_reference(extract_transition_matrix(without), extract_transition_matrix(with_r2),
                     reference = "independent: predictability does not change the weights", covers = "build_network()",
                     quantity = sprintf("weights with and without predictability, %s", method), tolerance = 0)
  }))
})

test_that("Ising weights at every rule and gamma equal IsingFit::IsingFit()", {
  skip_without_reference("IsingFit")
  binary <- psy_binary(n = 300L, p = 6L, rho = 0.5, seed = 815L)
  cases <- expand.grid(rule = c("AND", "OR"), gamma = c(0, 0.25, 0.5), stringsAsFactors = FALSE)
  invisible(Map(\(rule, gamma) {
    reference <- IsingFit::IsingFit(as.matrix(binary), family = "binomial", AND = identical(rule, "AND"), gamma = gamma,
                                    plot = FALSE, progressbar = FALSE)
    fit <- build_network(binary, method = "ising", params = list(rule = rule, gamma = gamma))
    expected <- reference$weiadj
    dimnames(expected) <- list(names(binary), names(binary))
    expect_reference(extract_transition_matrix(fit), expected, reference = "IsingFit", covers = "build_network()",
                     quantity = sprintf("Ising, rule %s, gamma %s", rule, gamma), tolerance = 1e-10)
    expect_reference(fit$thresholds, stats::setNames(reference$thresholds, names(binary)), reference = "IsingFit",
                     covers = "build_network()", quantity = sprintf("Ising thresholds, rule %s, gamma %s", rule, gamma),
                     tolerance = 1e-10)
  }, cases$rule, cases$gamma))
})

test_that("every default of the association methods equals the value written out", {
  binary <- psy_binary(n = 200L, p = 5L, rho = 0.5, seed = 816L)
  checks <- list(
    list(data = as_data, method = "cor", params = list(cor_method = "pearson", input_type = "auto", threshold = 0)),
    list(data = as_data, method = "pcor", params = list(cor_method = "pearson", input_type = "auto", threshold = 0)),
    list(data = as_data, method = "glasso",
         params = list(gamma = 0.5, nlambda = 100L, lambda.min.ratio = 0.01, penalize.diagonal = FALSE, refit = FALSE,
                       cor_method = "pearson", threshold = 0)),
    list(data = binary, method = "ising", params = list(gamma = 0.25, nlambda = 100L, rule = "AND"))
  )
  invisible(lapply(checks, \(check) {
    omitted <- build_network(check$data, method = check$method)
    written <- build_network(check$data, method = check$method, params = check$params, scaling = NULL, threshold = 0,
                             predictability = TRUE, level = NULL)
    expect_reference(extract_transition_matrix(written), extract_transition_matrix(omitted),
                     reference = "independent: the documented defaults written out", covers = "build_network()",
                     quantity = sprintf("defaults, %s", check$method), tolerance = 0)
  }))
  invisible(Map(\(fn, method, data) {
    expect_reference(extract_transition_matrix(fn(data)), extract_transition_matrix(build_network(data, method = method)),
                     reference = "independent: the wrapper is build_network() with its method",
                     covers = sprintf("build_%s()", method), quantity = sprintf("wrapper defaults, %s", method),
                     tolerance = 0)
  }, list(build_cor, build_pcor, build_glasso, build_ising), c("cor", "pcor", "glasso", "ising"),
     list(as_data, as_data, as_data, binary)))
})
