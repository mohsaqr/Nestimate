# tna models and group_tna objects as input (R/tna_input.R).
#
# A tna model is rebuilt as a Nestimate network from its own sequences, so
# every verb must give what it gives on the network built from the same
# sequences. Most blocks use a hand-built object with the layout of a tna
# model (weights, inits, labels, integer-coded data; attributes type,
# scaling, params), so they need no tna package.

# The sequences: 300 students of group_regulation_long, wide, with a group.
.tna_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    net <- build_network(group_regulation_long, method = "relative",
                         actor = "Actor", action = "Action", time = "Time")
    keep <- seq_len(300L)
    wide <- as.data.frame(net$data, stringsAsFactors = FALSE)[keep, seq_len(12L)]
    names(wide) <- paste0("V", seq_along(wide))
    rownames(wide) <- NULL
    half <- rep(c("High", "Low"), length.out = nrow(wide))
    cache <<- list(wide = wide, half = half)
    cache
  }
})

# A list with the fields and attributes of a tna model, built by hand.
.fake_tna <- function(wide, method = "relative", type = method,
                      scaling = character(0), params = list()) {
  net <- build_network(wide, method = method)
  labels <- rownames(net$weights)
  codes <- matrix(match(as.matrix(wide), labels), nrow = nrow(wide),
                  dimnames = list(NULL, paste0("T", seq_len(ncol(wide)))))
  structure(
    list(weights = net$weights, inits = net$initial, labels = labels,
         data = structure(codes, class = c("tna_seq_data", "matrix", "array"),
                          alphabet = labels, labels = labels)),
    type = type, scaling = scaling, params = params, class = "tna")
}

.fake_group_tna <- function(wide, group) {
  members <- lapply(split(wide, group), .fake_tna)
  structure(members, label = "Achiever", levels = names(members),
            type = "relative", scaling = character(0), class = "group_tna")
}

# The Nestimate group built from the same sequences.
.native_group <- function(wide, group) {
  build_network(cbind(wide, Achiever = group), method = "relative",
                group = "Achiever")
}


test_that("a tna model converts to the network built from its sequences", {
  f <- .tna_fixture()
  net <- build_network(f$wide, method = "relative")
  model <- .fake_tna(f$wide)

  built <- .tna_rebuild(model)
  expect_null(built$reason)
  converted <- as_netobject(model)
  expect_s3_class(converted, "netobject")
  expect_equal(converted$weights, net$weights, tolerance = 1e-12)
  expect_identical(converted$method, "relative")
  expect_identical(dim(converted$data), dim(net$data))
  expect_identical(unname(as.matrix(converted$data)),
                   unname(as.matrix(net$data)))
  expect_true(inherits(converted, "cograph_network"))
})

test_that("the converted weights equal the model's weights for every type", {
  f <- .tna_fixture()
  types <- c(relative = "relative", frequency = "frequency",
             co_occurrence = "co-occurrence", attention = "attention")
  same <- vapply(names(types), function(method) {
    model <- .fake_tna(f$wide, method = method, type = types[[method]])
    built <- .tna_rebuild(model)
    is.null(built$reason) && identical(built$net$method, method) &&
      isTRUE(all.equal(unname(built$net$weights), unname(model$weights),
                       tolerance = 1e-12))
  }, logical(1L))
  expect_true(all(same), info = paste(names(types)[!same], collapse = ", "))
})

test_that("a group_tna converts to a netobject_group named by its groups", {
  f <- .tna_fixture()
  group <- as_netobject(.fake_group_tna(f$wide, f$half))
  native <- .native_group(f$wide, f$half)
  expect_s3_class(group, "netobject_group")
  expect_identical(names(group), c("High", "Low"))
  expect_identical(attr(group, "group_col"), "Achiever")
  expect_equal(group$High$weights, native$High$weights, tolerance = 1e-12)
  expect_equal(group$Low$weights, native$Low$weights, tolerance = 1e-12)
})

test_that("resampling verbs give the netobject result on a tna model", {
  f <- .tna_fixture()
  net <- build_network(f$wide, method = "relative")
  model <- .fake_tna(f$wide)

  boot_net <- bootstrap_network(net, iter = 30, seed = 1)
  boot_tna <- bootstrap_network(model, iter = 30, seed = 1)
  expect_identical(boot_tna$ci_lower, boot_net$ci_lower)
  expect_identical(boot_tna$p_values, boot_net$p_values)

  expect_identical(
    summary(centrality_stability(model, iter = 10, seed = 1)),
    summary(centrality_stability(net, iter = 10, seed = 1)))
  strengths <- c("InStrength", "OutStrength")
  expect_identical(net_centrality(model, measures = strengths, loops = TRUE),
                   net_centrality(net, measures = strengths, loops = TRUE))
})

test_that("comparisons give the netobject result on tna models and groups", {
  skip_on_cran()
  f <- .tna_fixture()
  native <- .native_group(f$wide, f$half)
  models <- .fake_group_tna(f$wide, f$half)

  ref <- permutation(native$High, native$Low, iter = 30, seed = 1)
  expect_identical(
    permutation(models$High, models$Low, iter = 30, seed = 1)$p_values,
    ref$p_values)
  # the group is compared pairwise, as a netobject_group is
  from_group <- permutation(models, iter = 30, seed = 1)
  expect_s3_class(from_group, "net_permutation_group")
  expect_identical(names(from_group), "High vs Low")
  expect_identical(from_group[[1L]]$p_values, ref$p_values)

  cmp <- compare_networks(models, test = "permutation", iter = 30, seed = 1)
  expect_identical(cmp$permutations[[1L]]$p_values, ref$p_values)

  expect_identical(
    {
      set.seed(2)
      sequence_compare(models, sub = 2, iter = 20)$patterns
    },
    {
      set.seed(2)
      sequence_compare(native, sub = 2, iter = 20)$patterns
    })
  expect_identical(sequence_patterns(models, len = 2),
                   sequence_patterns(native, len = 2))
  expect_identical(mosaic_test(models), mosaic_test(native))
  expect_identical(state_distribution(models), state_distribution(native))
})

test_that("markov_order_test() reads the sequences of a tna model", {
  f <- .tna_fixture()
  net <- build_network(f$wide, method = "relative")
  model <- .fake_tna(f$wide)
  from_net <- markov_order_test(net, max_order = 2, n_perm = 20, seed = 1)
  from_tna <- markov_order_test(model, max_order = 2, n_perm = 20, seed = 1)
  # it used to read the fields of the object as sequences (other counts of
  # sequences and states, no warning)
  expect_identical(summary(from_tna), summary(from_net))
  # a scaled model has the same sequences
  scaled <- .fake_tna(f$wide, scaling = "minmax")
  expect_identical(
    summary(markov_order_test(scaled, max_order = 2, n_perm = 20, seed = 1)),
    summary(from_net))
  expect_identical(summary(path_dependence(model)),
                   summary(path_dependence(net)))
})

test_that("a model that cannot be rebuilt is refused where sequences are read", {
  f <- .tna_fixture()
  scaled <- .fake_tna(f$wide, scaling = "minmax")
  expect_match(.tna_rebuild(scaled)$reason, "scaled", fixed = TRUE)
  expect_error(bootstrap_network(scaled, iter = 10),
               class = "nestimate_tna_unrebuilt")
  expect_error(bootstrap_network(scaled, iter = 10),
               class = "nestimate_bad_input")
  expect_error(permutation(scaled, .fake_tna(f$wide), iter = 10),
               class = "nestimate_tna_unrebuilt")
  expect_error(mosaic_test(scaled), class = "nestimate_tna_unrebuilt")
  # verbs that read the weights only go on, on the model's own weights
  converted <- as_netobject(scaled)
  expect_null(converted$data)
  expect_identical(unname(converted$weights), unname(scaled$weights))
  expect_s3_class(net_centrality(scaled, measures = "InStrength",
                                 loops = TRUE), "data.frame")

  # weights the sequences do not reproduce
  altered <- .fake_tna(f$wide)
  altered$weights <- altered$weights * 0.5
  expect_match(.tna_rebuild(altered)$reason, "differ", fixed = TRUE)
  expect_error(centrality_stability(altered, iter = 5),
               class = "nestimate_tna_unrebuilt")

  # a type without an estimator, and a model without sequence data
  expect_match(.tna_rebuild(.fake_tna(f$wide, type = "betweenness"))$reason,
               "type", fixed = TRUE)
  no_data <- .fake_tna(f$wide)
  no_data$data <- NULL
  expect_match(.tna_rebuild(no_data)$reason, "no sequence data", fixed = TRUE)

  # one unreproducible member refuses the group, naming it
  models <- .fake_group_tna(f$wide, f$half)
  attr(models$Low, "scaling") <- "max"
  expect_error(permutation(models, iter = 10),
               regexp = "group 'Low'", class = "nestimate_tna_unrebuilt")
})

test_that("actor is read from the metadata a tna model keeps", {
  f <- .tna_fixture()
  n <- nrow(f$wide)
  meta <- data.frame(team = sprintf("t%02d", rep(seq_len(ceiling(n / 5)),
                                                 each = 5, length.out = n)),
                     stringsAsFactors = FALSE)
  with_meta <- function(model, rows) {
    attr(model$data, "meta") <- meta[rows, , drop = FALSE]
    model
  }
  model <- with_meta(.fake_tna(f$wide), seq_len(n))
  net <- build_network(f$wide, method = "relative")
  net$metadata <- meta
  expect_identical(.tna_input(model)$metadata, meta)
  from_tna <- bootstrap_network(model, iter = 30, seed = 2, actor = "team")
  from_net <- bootstrap_network(net, iter = 30, seed = 2, actor = "team")
  expect_equal(from_tna$ci_lower, from_net$ci_lower, tolerance = 1e-12)
  expect_equal(from_tna$actor, from_net$actor, tolerance = 1e-12)
  expect_equal(mosaic_test(model, actor = "team", method = "rao_scott"),
               mosaic_test(net, actor = "team", method = "rao_scott"),
               tolerance = 1e-12)

  # a group: every member carries its own rows of the metadata
  models <- .fake_group_tna(f$wide, f$half)
  rows <- split(seq_len(n), f$half)
  models[] <- Map(with_meta, unclass(models), rows[names(models)])
  perm <- permutation(models, iter = 30, seed = 3, actor = "team")
  expect_identical(perm[[1L]]$actor, "team")
  expect_s3_class(sequence_compare(models, sub = 2, iter = 20, min_freq = 1,
                                   actor = "team"),
                  "net_sequence_comparison")
  expect_s3_class(permutation_diagnostics(models, actor = "team", iter = 20,
                                          seed = 1),
                  "data.frame")
  # a column the metadata does not hold
  expect_error(bootstrap_network(model, iter = 10, actor = "absent"),
               class = "nestimate_actor_missing")
})

test_that("actor is refused on a tna model that keeps no metadata", {
  f <- .tna_fixture()
  model <- .fake_tna(f$wide)
  models <- .fake_group_tna(f$wide, f$half)
  expect_error(bootstrap_network(model, iter = 10, actor = "Group"),
               class = "nestimate_actor_missing")
  expect_error(permutation(models, iter = 10, actor = "Group"),
               class = "nestimate_actor_missing")
  expect_error(permutation_diagnostics(models, actor = "Group", iter = 10),
               class = "nestimate_actor_missing")
  expect_error(mosaic_test(model, actor = "Group"),
               class = "nestimate_actor_missing")
  expect_error(sequence_compare(models, sub = 2, iter = 10, actor = "Group"),
               class = "nestimate_actor_missing")
})

test_that("objects that are not tna models pass through the shared step", {
  f <- .tna_fixture()
  net <- build_network(f$wide, method = "relative")
  expect_identical(.tna_input(net), net)
  expect_identical(.tna_input(f$wide), f$wide)
  expect_null(.tna_input(NULL))
  expect_error(as_netobject(list(1)), class = "nestimate_bad_input")
})

test_that("models from the tna package give the netobject results", {
  skip_on_cran()
  skip_if_not_installed("tna")
  f <- .tna_fixture()
  net <- build_network(f$wide, method = "relative")
  native <- .native_group(f$wide, f$half)
  model <- tna::tna(f$wide)
  models <- tna::group_tna(f$wide, group = f$half)

  expect_null(.tna_rebuild(model)$reason)
  expect_identical(bootstrap_network(model, iter = 20, seed = 1)$ci_lower,
                   bootstrap_network(net, iter = 20, seed = 1)$ci_lower)
  expect_identical(
    summary(markov_order_test(model, max_order = 2, n_perm = 20, seed = 1)),
    summary(markov_order_test(net, max_order = 2, n_perm = 20, seed = 1)))

  group <- as_netobject(models)
  expect_identical(names(group), names(native))
  expect_identical(
    permutation(models, iter = 20, seed = 1)[[1L]]$p_values,
    permutation(native, iter = 20, seed = 1)[[1L]]$p_values)
  expect_identical(mosaic_test(models), mosaic_test(native))

  scaled <- tna::tna(f$wide, scaling = "minmax")
  expect_error(bootstrap_network(scaled, iter = 10),
               class = "nestimate_tna_unrebuilt")
})

test_that("a pruned tna model resamples as the full model; weights stay pruned", {
  f <- .tna_fixture()
  full <- .fake_tna(f$wide)
  pruned <- full
  pruned$weights[pruned$weights < 0.05] <- 0
  attr(pruned, "pruning") <- list(weights = pruned$weights,
                                  original = full$weights, active = TRUE)
  expect_false(isTRUE(all.equal(pruned$weights, full$weights)))
  # verbs that resample or read sequences use the model's data
  expect_equal(bootstrap_network(pruned, iter = 30, seed = 1)$ci_lower,
               bootstrap_network(full, iter = 30, seed = 1)$ci_lower,
               tolerance = 1e-12)
  # verbs that read weights see the model as pruned
  expect_equal(unname(as_netobject(pruned)$weights), unname(pruned$weights),
               tolerance = 1e-12)
  # a pruning record that does not match the sequences is not trusted
  wrong <- pruned
  attr(wrong, "pruning")$original[1L, 2L] <- 0.9
  expect_error(bootstrap_network(wrong, iter = 10),
               class = "nestimate_tna_unrebuilt")
})
