# House-rule fixes in mcml.R, macro_network.R, mlvar.R, gimme.R,
# outcome_model.R, mogen.R, hon.R, honem.R, hypa.R, extract_pathways.R and
# pathways.R: classed conditions, corrected printed text and labels, warnings
# for arguments that were dropped without a word. Bundled data only.

testthat::skip_on_cran()

wide_traj <- as.data.frame(trajectories)
traj_net <- build_network(wide_traj, method = "relative")
traj_clusters <- list(Engaged = c("Active", "Average"),
                      Out = "Disengaged")

# ---- classed match.arg() failures ------------------------------------------

test_that("a wrong choice raises nestimate_bad_input naming the argument", {
  expect_error(build_mogen(traj_net, criterion = "zz"),
               class = "nestimate_bad_input")
  expect_error(build_mogen(traj_net, criterion = "zz"), "'criterion'")
  expect_error(build_hon(traj_net, method = "zz"),
               class = "nestimate_bad_input")
  expect_error(pathways(traj_net, ho_method = "zz"),
               class = "nestimate_bad_input")
  expect_error(build_mcml(wide_traj, traj_clusters, type = "zz"),
               class = "nestimate_bad_input")
  expect_error(build_mcml(wide_traj, traj_clusters, method = "zz"),
               "'method' should be one of")
  expect_error(cluster_summary(traj_net$weights, traj_clusters, method = "zz"),
               class = "nestimate_bad_input")
  expect_error(
    extract_pathways(group_regulation_long, action = "Action",
                     group = "Actor", type = "zz"),
    class = "nestimate_bad_input")
  expect_error(
    outcome_model(chatgpt_srl, outcome = "CSU", predictors = "IV",
                  family = "zz"),
    class = "nestimate_bad_input")
})

test_that("a valid choice is matched as match.arg() matches it", {
  expect_identical(.ho_match_arg("b", "x", c("a", "b")), "b")
  expect_identical(.ho_match_arg(c("a", "b"), "x", c("a", "b")), "a")
  f <- function(kind = c("first", "second")) .ho_match_arg(kind, "kind")
  expect_identical(f(), "first")
  expect_identical(f("sec"), "second")
  expect_error(f("zz"), class = "nestimate_bad_input")
})

# ---- mogen / hon / honem / hypa --------------------------------------------

test_that("mogen_transitions() names the object it needs", {
  expect_error(mogen_transitions(list(), order = 1),
               class = "nestimate_bad_input")
  expect_error(mogen_transitions(list(), order = 1), "net_mogen")
  mg <- build_mogen(traj_net, max_order = 2)
  expect_error(plot(mg, type = "zz"), class = "nestimate_bad_input")
})

test_that("plot.net_honem() refuses dimensions the embedding lacks", {
  emb <- build_honem(build_hon(traj_net, max_order = 2), dim = 2)
  expect_error(plot(emb, dims = c(1, 5)), class = "nestimate_bad_input")
  expect_error(plot(emb, dims = 1), class = "nestimate_bad_input")
})

test_that("summary.net_hon() prints the node-order lines unchanged", {
  hon <- build_hon(traj_net, max_order = 2)
  out <- utils::capture.output(res <- summary(hon))
  orders <- table(lengths(strsplit(hon$nodes$label, " -> ", fixed = TRUE)))
  expect_true(all(sprintf("    Order %s: %d nodes", names(orders),
                          as.integer(orders)) %in% out))
  expect_identical(res, hon$edges)
})

test_that("build_hypa() warns when a requested order yields no edges", {
  # the trajectories hold 15 positions, so order 20 has no path to count
  expect_warning(
    hy <- build_hypa(traj_net, order = c(1, 20)),
    class = "nestimate_ignored_argument")
  expect_identical(hy$order, 1L)
  expect_no_warning(build_hypa(traj_net, order = 1))
  expect_error(summary(hy, type = "zz"), class = "nestimate_bad_input")
  expect_error(pathways(hy, type = "zz"), class = "nestimate_bad_input")
})

# ---- outcome_model ----------------------------------------------------------

test_that("a non-numeric outcome with 3+ values is a classed input error", {
  d <- transform(group_regulation_long, n = seq_along(Actor))
  expect_error(outcome_model(d, outcome = "Course", predictors = "n"),
               class = "nestimate_bad_input")
  expect_error(outcome_model(d, outcome = "Course", predictors = "n"),
               "numeric outcome")
})

test_that("the dropped-predictor text and the CI label say what happened", {
  d <- transform(chatgpt_srl, flag = IV > 4)
  fit <- outcome_model(d, outcome = "CSU", predictors = c("IV", "flag"),
                       ci_level = 0.9)
  expect_identical(fit$dropped, "flag")
  out <- utils::capture.output(print(fit))
  expect_true(any(grepl("flag (not numeric or no variance)", out,
                        fixed = TRUE)))
  expect_identical(plot(fit)$labels$x, "Estimate (90% CI)")
  expect_identical(
    plot(outcome_model(d, outcome = "CSU", predictors = "IV"))$labels$x,
    "Estimate (95% CI)")
  expect_error(outcome_model(d, outcome = "CSU", predictors = "flag"),
               class = "nestimate_no_variance")
  expect_error(outcome_model(d, outcome = "CSU", predictors = "flag"),
               "non-numeric or constant")
})

# ---- mcml / macro_network ---------------------------------------------------

test_that("internal mcml validators raise classed errors", {
  expect_error(.mcml_layer("a", NULL, "a"), class = "nestimate_bad_input")
  expect_error(.build_cluster_lookup(list(G = "A"), c("A", "B")),
               class = "nestimate_bad_input")
  expect_error(.build_cluster_lookup(c("x", "y"), c("A", "B", "C")),
               class = "nestimate_bad_input")
  expect_error(net_aggregate_weights(1:3, method = "zz"),
               class = "nestimate_bad_input")
})

test_that("unnamed extra arguments are named by position", {
  mc <- build_mcml(wide_traj, traj_clusters)
  expect_error(summary(mc, 1), class = "nestimate_bad_input")
  expect_error(summary(mc, 1), "unsupported argument: ..1", fixed = TRUE)
})

test_that("print.mcml_layer() describes the weights it holds", {
  prob <- utils::capture.output(print(build_mcml(wide_traj, traj_clusters)$macro))
  expect_identical(prob[1L],
                   "MCML layer (transition probabilities)  [directed]")
  raw <- utils::capture.output(
    print(build_mcml(wide_traj, traj_clusters, type = "raw")$macro))
  expect_identical(raw[1L], "MCML layer (weights)  [directed]")
  co <- utils::capture.output(
    print(build_mcml(wide_traj, traj_clusters, type = "cooccurrence")$macro))
  expect_identical(co[1L], "MCML layer (weights)  [symmetric]")
})

test_that("print.mcml() lists every cluster with its members", {
  mc <- build_mcml(wide_traj, traj_clusters)
  out <- utils::capture.output(print(mc))
  expect_true("  Engaged (2): Active, Average" %in% out)
  expect_true("  Out (1): Disengaged" %in% out)
})

test_that("cluster_summary() aggregates each cluster pair as before", {
  w <- traj_net$weights
  for_pair <- function(i, j, method) {
    net_aggregate_weights(
      as.vector(w[traj_clusters[[i]], traj_clusters[[j]]]), method,
      length(traj_clusters[[i]]) * length(traj_clusters[[j]]))
  }
  check <- function(method) {
    got <- cluster_summary(w, traj_clusters, method = method)$macro$weights
    want <- outer(1:2, 1:2, Vectorize(function(i, j) for_pair(i, j, method)))
    expect_equal(unname(got), want, tolerance = 1e-12)
  }
  invisible(lapply(c("sum", "mean", "max", "density"), check))
})

test_that("shaping arguments on an unchanged mcml warn instead of vanishing", {
  mc <- build_mcml(wide_traj, traj_clusters)
  expect_warning(same <- build_mcml(mc, exclude = "Active"),
                 class = "nestimate_ignored_argument")
  expect_identical(same, mc)
  expect_no_warning(build_mcml(mc))
})

test_that("`type` on a non-count network warns that it is ignored", {
  att <- build_network(wide_traj, method = "attention")
  expect_warning(build_mcml(att, traj_clusters, type = "raw"),
                 class = "nestimate_ignored_argument")
  expect_no_warning(build_mcml(att, traj_clusters))
  expect_warning(build_mcml(traj_net$weights, traj_clusters, type = "tna"),
                 class = "nestimate_ignored_argument")
  expect_no_warning(build_mcml(traj_net, traj_clusters, type = "raw"))
})

test_that("an `end` label without a cluster names the argument", {
  expect_error(build_mcml(wide_traj, traj_clusters, end = TRUE),
               class = "nestimate_bad_input")
  expect_error(build_mcml(wide_traj, traj_clusters, end = TRUE),
               "Unmapped nodes: End. `end` adds", fixed = TRUE)
  with_end <- c(traj_clusters, list(Stop = "End"))
  expect_s3_class(build_mcml(wide_traj, with_end, end = TRUE), "mcml")
})

test_that("the no-source errors describe every route that lacks sequences", {
  from_net <- build_mcml(traj_net$weights, traj_clusters)
  expect_error(macro_network(from_net, expand = "Engaged"),
               class = "nestimate_no_expand_source")
  expect_error(macro_network(from_net, expand = "Engaged"),
               "an edge list or a network")
  expect_error(as_htna(from_net), "no expandable node-level source")
  expect_error(as_htna(from_net), "an edge list or a network")
})

# ---- mlvar / gimme ----------------------------------------------------------

mlvar_panel <- transform(srl_strategies,
                         id = rep(seq_len(10), length.out = nrow(srl_strategies)))
mlvar_vars <- c("Planning", "Monitoring", "Evaluating")

test_that("build_mlvar() validation errors are classed and say what is wrong", {
  expect_error(build_mlvar(mlvar_panel, vars = "Planning", id = "id"),
               class = "nestimate_bad_input")
  expect_error(build_mlvar(mlvar_panel, vars = "Planning", id = "id"),
               "`vars` must be a character vector of at least 2")
  expect_error(build_mlvar(mlvar_panel, vars = mlvar_vars, id = "id", lag = 0),
               "`lag` must be a single number >= 1")
  expect_error(coefs(1), class = "nestimate_bad_input")
})

test_that("print/summary/plot of a net_mlvar name verbs and survive renaming", {
  fit <- suppressWarnings(build_mlvar(mlvar_panel, vars = mlvar_vars, id = "id"))
  out <- utils::capture.output(print(fit))
  expect_false(any(grepl("$", out, fixed = TRUE)))
  expect_true(any(grepl("coefs(fit)", out, fixed = TRUE)))
  expect_true(any(grepl("^  temporal +3 x 3 directed", out)))

  renamed <- rename_models(fit, c("T", "C", "B"))
  expect_identical(utils::capture.output(print(renamed)), out)
  utils::capture.output(tab <- summary(renamed))
  expect_identical(tab, coefs(fit))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(plot(fit, layer = "zz"), class = "nestimate_bad_input")
})

test_that("build_gimme() validation errors are classed and named", {
  skip_if_not_installed("lavaan")
  expect_error(build_gimme(mlvar_panel, vars = "Planning", id = "id"),
               "`vars` must be a character vector of at least 2")
  expect_error(build_gimme(mlvar_panel, vars = mlvar_vars, id = "nope"),
               class = "nestimate_bad_input")
  expect_error(
    build_gimme(mlvar_panel, vars = mlvar_vars, id = "id", groupcutoff = 0),
    "`groupcutoff` must be a single number in (0, 1]", fixed = TRUE)
})

test_that("build_gimme(seed =) leaves the session's random stream alone, and coefs() works", {
  skip_if_not_installed("lavaan")
  panel <- transform(chatgpt_srl,
                     id = rep(seq_len(3), length.out = nrow(chatgpt_srl)))
  set.seed(1)
  expected <- stats::runif(1)
  set.seed(1)
  fit <- suppressWarnings(
    build_gimme(panel, vars = c("CSU", "IV"), id = "id", seed = 5))
  expect_identical(stats::runif(1), expected)
  # Nestimate's coefs() used to have no method for a net_gimme
  tab <- coefs.net_gimme(fit)
  expect_s3_class(tab, "data.frame")
  expect_true(all(c("subject", "network", "from", "to", "weight") %in%
                    names(tab)))
})
