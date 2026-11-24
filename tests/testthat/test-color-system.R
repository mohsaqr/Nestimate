testthat::skip_on_cran()

# One colour system: Okabe-Ito by default, dealt over the sorted states, the
# same colour for a state in every figure; set_state_colors() changes them
# everywhere; `state_colors =` changes one figure.

.okabe3 <- c(Active = "#E69F00", Average = "#56B4E9", Disengaged = "#009E73")
.fills_of <- function(p, keys) {
  built <- ggplot2::ggplot_build(p)
  sc <- built$plot$scales$get_scales("fill")
  stats::setNames(sc$map(keys), keys)
}

test_that("with no colours set, a network uses Okabe-Ito in sorted state order", {
  net <- build_network(trajectories, method = "relative")
  expect_equal(stats::setNames(state_colors(net)$color, state_colors(net)$state),
               .okabe3)
  # the cograph contract carries the same colours, in node order
  expect_equal(net$meta$splot$defaults$node_fill,
               unname(.okabe3[net$nodes$label]))
})

test_that("the default colour of a state does not depend on node order", {
  codes <- data.frame(T1 = c("1", "2", "3", "4"), T2 = c("2", "3", "4", "1"),
                      T3 = c("3", "4", "1", "2"), stringsAsFactors = FALSE)
  net <- build_network(codes, method = "relative",
                       labels = c("1" = "Mixed", "2" = "Positive",
                                  "3" = "Low", "4" = "Distressed"))
  tb <- state_colors(net)
  expected <- stats::setNames(c("#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                              c("Distressed", "Low", "Mixed", "Positive"))
  expect_equal(stats::setNames(tb$color, tb$state)[names(expected)], expected)
})

test_that("a group and its members share one palette, also for a missing state", {
  d <- data.frame(T1 = c("A", "B", "A", "C"), T2 = c("B", "A", "B", "A"),
                  T3 = c("A", "B", "A", "C"), g = c("x", "x", "y", "y"),
                  stringsAsFactors = FALSE)
  grp <- build_network(d, method = "relative", group = "g")
  # group x has no state C, so on its own it would deal A, B the same colours
  # anyway; group y has A, B, C - all three must keep the group's colours
  pal <- stats::setNames(state_colors(grp)$color, state_colors(grp)$state)
  member <- grp[["y"]]
  expect_equal(member$meta$splot$defaults$node_fill,
               unname(pal[member$nodes$label]))
  set <- set_state_colors(grp, c(C = "#000000"))
  expect_equal(set[["y"]]$meta$splot$defaults$node_fill[set[["y"]]$nodes$label == "C"],
               "#000000")
})

test_that("a landscape draws each state in the network's colour, by name", {
  skip_if_not_installed("ggplot2")
  states <- c("Disengaged", "Average", "Active")   # not alphabetical
  net <- build_network(trajectories, method = "relative")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  land <- stability_landscape(net, states = states, iter = 20, seed = 1)
  expect_equal(.fills_of(plot(land), states), .okabe3[states])
  coloured <- set_state_colors(net, c(Average = "#000000"))
  land2 <- stability_landscape(coloured, states = states, iter = 20, seed = 1)
  expect_equal(unname(.fills_of(plot(land2), "Average")), "#000000")
  # a one-figure override
  expect_equal(unname(.fills_of(plot(land2, state_colors = c(Active = "#CC79A7")),
                                "Active")), "#CC79A7")
  trend <- landscape_trend(coloured, states = states, iter = 20, seed = 1)
  expect_equal(unname(.fills_of(plot(trend), "Average")), "#000000")
})

test_that("plot_state_frequencies() takes state_colors; colors is deprecated", {
  skip_if_not_installed("ggplot2")
  net <- build_network(trajectories, method = "relative")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_no_warning(plot_state_frequencies(net, state_colors = c(Active = "#000000")))
  expect_warning(plot_state_frequencies(net, colors = c(Active = "#000000")),
                 class = "nestimate_deprecated")
})

# ---- audit 2026-11-23 (Codex) ------------------------------------------------

test_that("numeric state codes are dealt 2 before 10 and never miss a colour", {
  d <- data.frame(t1 = c(2L, 10L, 2L), t2 = c(10L, 2L, 10L),
                  t3 = c(2L, 2L, 10L))
  net <- build_network(d, method = "relative")
  pal <- .stored_state_colors(net)
  expect_equal(unname(pal[c("2", "10")]), c("#E69F00", "#56B4E9"))
  expect_false(anyNA(.fill_state_colors(c(2L, 10L), pal)))
})

test_that("clearing or re-setting a group palette starts from the defaults", {
  g <- build_network(data.frame(T1 = c("A", "B", "A", "B"),
                                T2 = c("B", "A", "B", "A"),
                                g = c("x", "x", "y", "y")),
                     method = "relative", group = "g")
  cleared <- set_state_colors(set_state_colors(g, c(A = "red")), NULL)
  expect_equal(state_colors(cleared)$color, c("#E69F00", "#56B4E9"))
  taken <- set_state_colors(g, c(A = "#56B4E9"))
  expect_equal(state_colors(taken)$color, c("#56B4E9", "#E69F00"))
})

test_that("a palette entry for an absent state does not shift the defaults", {
  net <- build_network(data.frame(T1 = c("A", "B"), T2 = c("B", "A")),
                       method = "relative")
  shifted <- suppressMessages(set_state_colors(net, c(absent = "#E69F00")))
  expect_equal(state_colors(shifted)$color, c("#E69F00", "#56B4E9"))
})

test_that("a figure override keeps the object's colours for the other states", {
  net <- build_network(data.frame(T1 = c("A", "B"), T2 = c("B", "A")),
                       method = "relative")
  red <- set_state_colors(net, c(B = "#D55E00"))
  expect_equal(unname(.with_override(red, c(A = "#CC79A7"))[c("A", "B")]),
               c("#CC79A7", "#D55E00"))
})

test_that("plot_state_frequencies() draws the palette, by name and abbreviated", {
  skip_if_not_installed("ggplot2")
  net <- set_state_colors(build_network(trajectories, method = "relative"),
                          c(Disengaged = "#000000"))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  p <- plot_state_frequencies(net, style = "bar")$plot
  expect_equal(unname(.fills_of(p, "Disengaged")), "#000000")
  ab <- plot_state_frequencies(net, style = "bar", abbreviate = TRUE)$plot
  short <- unname(abbreviate(c("Active", "Average", "Disengaged"), 3L)[3])
  expect_equal(unname(.fills_of(ab, short)), "#000000")
  # an old positional call still binds label_size, not the new argument
  expect_identical(names(formals(getS3method("plot_state_frequencies",
                                             "netobject")))[10], "label_size")
})

test_that("a group member missing the first state keeps the group's colours", {
  d <- data.frame(T1 = c("A", "B", "B", "C"), T2 = c("B", "A", "C", "B"),
                  g = c("x", "x", "y", "y"), stringsAsFactors = FALSE)
  grp <- build_network(d, method = "relative", group = "g")
  pal <- stats::setNames(state_colors(grp)$color, state_colors(grp)$state)
  y <- grp[["y"]]                      # states B and C, no A
  expect_equal(y$meta$splot$defaults$node_fill, unname(pal[y$nodes$label]))
  expect_false(identical(unname(pal["B"]), "#E69F00"))
})

test_that("pruned groups, differences and entropy networks keep the colours", {
  net <- set_state_colors(build_network(trajectories, method = "relative"),
                          c(Active = "#000000"))
  grp <- set_state_colors(
    build_network(cbind(trajectories,
                        g = rep(c("a", "b"), length.out = nrow(trajectories))),
                  method = "relative", group = "g"),
    c(Active = "#000000"))
  pruned <- net_prune(grp, threshold = 0.1)
  expect_identical(attr(pruned, "state_colors"), attr(grp, "state_colors"))
  expect_identical(attr(pruned, "group_col"), attr(grp, "group_col"))
  diff_net <- subtract_networks(net, net)
  expect_identical(diff_net$meta$splot$renderer, "difference")
  expect_equal(diff_net$meta$splot$defaults$node_fill,
               net$meta$splot$defaults$node_fill)
  eb <- entropy_bayes(net, draws = 100, seed = 1)
  expect_equal(eb$network$meta$splot$defaults$node_fill,
               net$meta$splot$defaults$node_fill)
})

test_that("wtna, mcml and their conversions carry the colours", {
  w <- wtna(data.frame(B = c(1, 0, 1), A = c(0, 1, 0)), codes = c("B", "A"))
  expect_equal(w$meta$splot$defaults$node_fill,
               unname(c(A = "#E69F00", B = "#56B4E9")[w$nodes$label]))
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time")
  net <- set_state_colors(net, c(plan = "#000000"))
  clusters <- list(R = c("monitor", "plan", "coregulate", "adapt"),
                   S = c("cohesion", "consensus", "discuss", "emotion",
                         "synthesis"))
  mc <- build_mcml(net, clusters = clusters)
  states <- unlist(clusters, use.names = FALSE)
  from_net <- stats::setNames(state_colors(net)$color, state_colors(net)$state)
  from_mc <- stats::setNames(state_colors(mc)$color, state_colors(mc)$state)
  expect_equal(from_mc[states], from_net[states])
  layers <- as_tna(mc)
  expect_equal(layers$R$meta$splot$defaults$node_fill,
               unname(from_net[layers$R$nodes$label]))
})

test_that("two nodes with one label still get one colour each", {
  net <- build_network(trajectories, method = "relative")
  net$nodes$labels <- c("X", "X", "Y")
  stamped <- .stamp_node_fill(net)$meta$splot$defaults$node_fill
  expect_length(stamped, 3L)
})

test_that("the magnitude plot colours states as every figure does", {
  skip_if_not_installed("ggplot2")
  md <- magnitude_difference(as.data.frame(trajectories))
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  p <- plot(md, state_colors = c(Active = "#000000"))
  expect_s3_class(p, "ggplot")
  # the source states colour their labels
  drawn <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) d$colour))
  expect_true("#000000" %in% drawn)
})

test_that("labels are white only on fills where black text cannot be read", {
  expect_equal(.label_contrast(c("#000000", "#0072B2", "#E69F00", "#999999")),
               c("white", "black", "black", "black"))
  net <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time")
  stamped <- net$meta$splot$defaults
  # nine states: the ninth Okabe-Ito colour is black, its label white
  expect_equal(stamped$label_color[stamped$node_fill == "#000000"], "white")
  expect_true(all(stamped$label_color[stamped$node_fill != "#000000"] == "black"))
})
