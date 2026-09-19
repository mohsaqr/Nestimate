# Object-level state colours: set once, honoured by every figure.

make_net <- function() {
  build_network(group_regulation_long, method = "relative",
                actor = "Actor", action = "Action", time = "Time")
}

make_mc <- function() {
  build_mcml(group_regulation_long,
             clusters = list(Cognitive  = c("discuss", "synthesis",
                                            "consensus", "cohesion"),
                             Regulation = c("plan", "monitor", "adapt",
                                            "coregulate"),
                             Affective  = "emotion"),
             actor = "Actor", action = "Action", time = "Time")
}

test_that("state_colors() reports what the object will draw with", {
  net <- set_state_colors(make_net(), c(plan = "#0072B2", emotion = "#CC79A7"))
  tb  <- state_colors(net)

  expect_s3_class(tb, "data.frame")
  expect_named(tb, c("state", "color", "source"))
  expect_equal(nrow(tb), length(net$nodes$name))
  expect_equal(tb$color[tb$state == "plan"], "#0072B2")
  expect_equal(tb$source[tb$state == "plan"], "set")
  expect_equal(tb$source[tb$state == "adapt"], "default")
  # a default never repeats a colour the palette pinned
  expect_equal(anyDuplicated(tb$color), 0L)
})

test_that("the figures agree with the table", {
  skip_if_not_installed("ggplot2")
  net <- set_state_colors(make_net(), c(plan = "#0072B2", emotion = "#CC79A7"))
  tb  <- state_colors(net)

  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d <- sequence_plot(net, type = "distribution")
  expect_equal(unname(d$palette[seq_len(nrow(tb))]), tb$color)
  # an explicit argument still wins for one figure
  one <- sequence_plot(net, type = "distribution",
                       state_colors = c(plan = "#000000"))
  expect_equal(unname(one$palette[tb$state == "plan"]), "#000000")
  # ... and does not disturb the object
  expect_equal(state_colors(net)$color[tb$state == "plan"], "#0072B2")
})

test_that("the palette is mirrored into the cograph producer contract", {
  net <- set_state_colors(make_net(), c(plan = "#0072B2", emotion = "#CC79A7"))
  stamped <- net$meta$splot$defaults$node_fill

  expect_equal(length(stamped), nrow(net$nodes))
  # in node order, so cograph needs no name matching of its own
  expect_equal(stamped, unname(state_colors(net)$color[
    match(net$nodes$name, state_colors(net)$state)]))
  # removing the palette removes the stamp
  expect_null(set_state_colors(net, NULL)$meta$splot$defaults$node_fill)
})

test_that("set_state_colors works on every data-bearing class", {
  mc  <- set_state_colors(make_mc(), c(plan = "#0072B2", Affective = "#CC79A7"))
  tb  <- state_colors(mc)
  # an mcml carries its cluster names as keys too
  expect_true(all(c("plan", "Affective") %in% tb$state))
  expect_equal(tb$color[tb$state == "Affective"], "#CC79A7")

  grp <- build_network(group_regulation_long, method = "relative",
                       actor = "Actor", action = "Action", time = "Time",
                       group = "Course")
  skip_if_not(inherits(grp, "netobject_group"))
  grp <- set_state_colors(grp, c(plan = "#0072B2"))
  expect_equal(state_colors(grp)$color[state_colors(grp)$state == "plan"],
               "#0072B2")
  # each member carries it, so a member drawn on its own matches the group
  expect_equal(grp[[1L]]$state_colors[["plan"]], "#0072B2")

  expect_error(set_state_colors(1:3, c(a = "red")), "needs a netobject")
  expect_error(state_colors(1:3), "needs a netobject")
})

test_that("the replacement form and a reused palette behave", {
  net <- make_net()
  state_colors(net) <- c(plan = "#0072B2")
  expect_equal(state_colors(net)$color[state_colors(net)$state == "plan"],
               "#0072B2")

  # a project-wide palette naming states this object does not carry
  expect_message(set_state_colors(make_net(), c(plan = "#0072B2",
                                                Approve = "#2CA02C")),
                 "1 of 2 names are not drawn")
  expect_error(set_state_colors(make_net(), c("#0072B2", "#D55E00")),
               "must be named")
})

test_that("every figure of one object draws the same colour for a state", {
  skip_if_not_installed("ggplot2")
  # The regression: plot_state_frequencies() orders states by frequency and
  # sequence_plot() alphabetically, so dealing the default colours per figure
  # gave the same state two different colours.
  net <- set_state_colors(make_net(),
                          c(plan = "#0072B2", monitor = "#D55E00",
                            emotion = "#CC79A7", discuss = "#009E73"))
  tb <- state_colors(net)

  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  d  <- sequence_plot(net, type = "distribution")
  sf <- plot_state_frequencies(net)

  seq_pal <- stats::setNames(d$palette[seq_len(nrow(tb))], sort(tb$state))
  g       <- ggplot2::get_guide_data(sf$plot, "fill")
  freq_pal <- stats::setNames(g$fill, g$.label)

  expect_equal(unname(seq_pal[tb$state]), tb$color)
  expect_equal(unname(freq_pal[tb$state]), tb$color)
})

test_that("in-tile labels take the ink with more contrast than the fill", {
  # whichever of the two inks wins the WCAG contrast ratio
  expect_equal(.contrast_label_color(c("#000000", "#FFFFFF")),
               c("white", "grey15"))
  expect_equal(.contrast_label_color("#0072B2"), "white")
  # a mid grey: the dark ink has 5.4:1 against it, white only 2.9:1, so a
  # fixed 0.4-luminance cutoff would pick white and be wrong
  expect_equal(.contrast_label_color("#999999"), "grey15")
  # the dark colours a user might pin
  expect_equal(.contrast_label_color(c("#4B1D3F", "#1F3A93", "#005F73")),
               rep("white", 3L))
  # luminance is the WCAG definition, not a mean of the channels
  expect_equal(.rel_luminance(c("#000000", "#FFFFFF")), c(0, 1))
})

test_that("a dark tile gets a light label in the drawn figure", {
  skip_if_not_installed("ggplot2")
  net <- set_state_colors(make_net(),
                          c(plan = "#0072B2", monitor = "#D55E00",
                            emotion = "#CC79A7", discuss = "#009E73"))
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  sf <- plot_state_frequencies(net)
  layers <- ggplot2::ggplot_build(sf$plot)$data
  inks <- unique(unlist(lapply(layers, function(d) d$colour[!is.na(d$colour)])))

  # the object puts one state on #000000, so both inks must be in use
  expect_true(any(state_colors(net)$color == "#000000"))
  expect_true(all(c("white", "grey15") %in% inks))
})
