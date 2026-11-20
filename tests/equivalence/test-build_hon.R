# build_hon(): higher-order networks by BuildHON (Xu, Wickramarathne and Chawla 2016) and BuildHON+ (Saebi et al.
# 2020). The reference is a Python reimplementation of both algorithms written for these tests (data/pyhon.py,
# data/pyhon_plus.py, standard library only), run through reticulate with the interpreter of RETICULATE_PYTHON.
# Re-expressed from local_testing_and_equivalence/test-hon-pyhon-equiv.R (which hard-coded
# /opt/homebrew/bin/python3 and so skipped off the Mac). Every edge must exist on both sides (from, to keys) with
# the same transition probability. The edge table is read from the result's field `ho_edges` (API gap:
# summary.net_hon() returns integer node indices, and there is no edge accessor).

# Python's node notation "C|B.A" (state C after A then B) as Nestimate's "A -> B -> C".
ho_pipe_to_arrow <- function(node) {
  parts <- strsplit(node, "|", fixed = TRUE)[[1L]]
  if (length(parts) < 2L || !nzchar(parts[[2L]])) return(parts[[1L]])
  paste(c(rev(strsplit(parts[[2L]], ".", fixed = TRUE)[[1L]]), parts[[1L]]), collapse = " -> ")
}

# The Python edges as a table: from (history node), to (the next state), probability.
ho_python_edges <- function(edges) {
  data.frame(from = vapply(edges, \(edge) ho_pipe_to_arrow(edge[["from"]]), character(1)),
             to = vapply(edges, \(edge) utils::tail(strsplit(ho_pipe_to_arrow(edge[["to"]]), " -> ", fixed = TRUE)[[1L]], 1L),
                         character(1)),
             probability = vapply(edges, \(edge) edge[["weight"]], numeric(1)))
}

ho_run_pyhon <- function(trajectories, max_order, min_freq) {
  reticulate::py_run_file(reference_data("data/pyhon.py"))
  result <- reticulate::py$run_buildhon(reticulate::r_to_py(trajectories), max_order = as.integer(max_order),
                                        min_support = as.integer(min_freq))
  ho_python_edges(result$edges)
}

ho_run_pyhon_plus <- function(trajectories, max_order, min_freq) {
  reticulate::py_run_file(reference_data("data/pyhon_plus.py"))
  ho_python_edges(reticulate::py$run_hon_plus(reticulate::r_to_py(trajectories), as.integer(max_order),
                                              as.integer(min_freq)))
}

hon_inputs <- list(
  list(name = "4 states, repeated path", max_order = 5L, min_freq = 1L,
       trajectories = rep(list(c("A", "B", "C", "D", "B", "A")), 5L)),
  list(name = "higher-order dependency", max_order = 5L, min_freq = 1L,
       trajectories = list(c("A", "B", "C", "A", "B", "C", "A"), c("D", "B", "A", "D", "B", "A", "D"),
                           c("A", "B", "C", "A", "B", "C", "A"), c("D", "B", "A", "D", "B", "A", "D"),
                           c("A", "B", "C", "D", "B", "A", "D"), c("A", "B", "C", "A", "B", "C", "A"),
                           c("D", "B", "A", "D", "B", "A", "D"))),
  list(name = "5 states", max_order = 5L, min_freq = 1L,
       trajectories = list(c("E", "A", "B", "C", "D", "E", "A"), c("B", "C", "D", "A", "E", "B", "C"),
                           c("A", "B", "C", "E", "D", "A", "B"), c("D", "E", "A", "B", "C", "D", "E"),
                           c("C", "D", "A", "B", "E", "C", "D"), c("A", "B", "C", "D", "E", "A", "B"),
                           c("E", "D", "C", "B", "A", "E", "D"), c("B", "A", "E", "D", "C", "B", "A"),
                           c("A", "B", "C", "D", "A", "B", "C"), c("D", "C", "B", "A", "D", "C", "B"))),
  list(name = "min_freq 3", max_order = 5L, min_freq = 3L,
       trajectories = list(c("A", "B", "C", "D"), c("A", "B", "D", "C"), c("A", "B", "C", "D"),
                           c("D", "C", "B", "A"), c("A", "B", "C", "D"), c("C", "B", "A", "D"))),
  list(name = "max_order 2", max_order = 2L, min_freq = 1L,
       trajectories = list(c("A", "B", "C", "A", "B", "D"), c("A", "B", "C", "A", "B", "C"),
                           c("D", "B", "A", "D", "B", "A"), c("A", "B", "C", "D", "B", "A"),
                           c("A", "B", "D", "A", "B", "C")))
)

hon_plus_inputs <- list(
  list(name = "4 states", min_freq = 1L,
       trajectories = list(c("A", "B", "C", "D"), c("A", "B", "D", "C"), c("B", "C", "A", "D"), c("D", "A", "B", "C"))),
  list(name = "higher-order dependency", min_freq = 1L,
       trajectories = c(rep(list(c("X", "A", "C")), 4L), rep(list(c("Y", "A", "D")), 4L),
                        list(c("A", "B", "A"), c("B", "A", "B")))),
  list(name = "min_freq 3", min_freq = 3L,
       trajectories = list(c("A", "B", "C", "D"), c("A", "B", "D", "C"), c("A", "B", "C", "D"),
                           c("D", "C", "B", "A"), c("A", "B", "C", "D"), c("C", "B", "A", "D")))
)

test_that("BuildHON edges and probabilities equal the Python reimplementation", {
  skip_without_reference("reticulate")
  invisible(lapply(hon_inputs, \(input) {
    found <- build_hon(input$trajectories, max_order = input$max_order, min_freq = input$min_freq, method = "hon")
    expect_reference(subset(found$ho_edges, select = c(from, to, probability)),
                     ho_run_pyhon(input$trajectories, input$max_order, input$min_freq),
                     reference = "independent: BuildHON reimplemented in Python (data/pyhon.py)", covers = "build_hon()",
                     quantity = sprintf("HON, %s", input$name), tolerance = 1e-10, key = c("from", "to"))
  }))
})

test_that("BuildHON+ edges and probabilities equal the Python reimplementation", {
  skip_without_reference("reticulate")
  invisible(lapply(hon_plus_inputs, \(input) {
    found <- build_hon(input$trajectories, max_order = 99L, min_freq = input$min_freq, method = "hon+")
    expect_reference(subset(found$ho_edges, select = c(from, to, probability)),
                     ho_run_pyhon_plus(input$trajectories, 99L, input$min_freq),
                     reference = "independent: BuildHON+ reimplemented in Python (data/pyhon_plus.py)",
                     covers = "build_hon()", quantity = sprintf("HON+, %s", input$name), tolerance = 1e-10,
                     key = c("from", "to"))
  }))
})

test_that("summary() returns the rewired network's edges, as documented with integer node indices, equal to BuildHON's", {
  skip_without_reference("reticulate")
  invisible(lapply(hon_inputs, \(input) {
    found <- build_hon(input$trajectories, max_order = input$max_order, min_freq = input$min_freq, method = "hon")
    edges <- summary(found)
    # ?build_hon: summary() returns integer `from`/`to` indices into the node table, whose `label` names the node.
    named <- data.frame(from = found$nodes$label[edges$from], to = found$nodes$label[edges$to], weight = edges$weight)
    reticulate::py_run_file(reference_data("data/pyhon.py"))
    python <- reticulate::py$run_buildhon(reticulate::r_to_py(input$trajectories), max_order = input$max_order,
                                          min_support = input$min_freq)$edges
    reference <- data.frame(from = vapply(python, \(edge) ho_pipe_to_arrow(edge[["from"]]), character(1)),
                            to = vapply(python, \(edge) ho_pipe_to_arrow(edge[["to"]]), character(1)),
                            weight = vapply(python, \(edge) edge[["weight"]], numeric(1)))
    expect_reference(named, reference, reference = "independent: BuildHON reimplemented in Python (data/pyhon.py)",
                     covers = "summary.net_hon()", quantity = sprintf("rewired edges of summary(), %s", input$name),
                     tolerance = 1e-10, key = c("from", "to"))
  }))
})
