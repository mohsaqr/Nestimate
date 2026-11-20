# build_gimme(), coefs.net_gimme(): Group Iterative Multiple Model Estimation (Gates and Molenaar 2012). Nestimate
# forwards the fit to idiographic::fit_gimme() (internal); the reference is gimme::gimme() 10.0 run on the same
# person series. Re-expressed from local_testing_and_equivalence/test-equiv-gimme.R (path counts and coefficients
# against gimme run in a separate R process); the comparison shape follows idiographic's own declared gimme test.
# gimme returns its path estimates rounded to four decimals, as fit_gimme() does, hence the tolerance 0.

# One person's matrix of a gimme result as a long table of non-zero cells, in coefs()' shape.
gimme_paths <- function(reference) {
  rows <- Map(\(cells, subject) {
    data.frame(subject = subject, row = rownames(cells)[as.vector(row(cells))],
               column = colnames(cells)[as.vector(col(cells))], weight = as.vector(cells))
  }, reference$path_est_mats, as.character(seq_along(reference$path_est_mats)))
  cells <- subset(do.call(rbind, rows), weight != 0)
  data.frame(subject = cells$subject, network = ifelse(grepl("lag$", cells$column), "temporal", "contemporaneous"),
             from = sub("lag$", "", cells$column), to = cells$row, weight = cells$weight)
}

test_that("gimme path estimates equal gimme::gimme() on simulated person series", {
  skip_without_reference("gimme")
  panel <- psy_panel(n_id = 10L, days = 1L, beeps = 100L, p = 3L, seed = 501L)
  vars <- c("V1", "V2", "V3")
  series <- lapply(split(panel, panel$id), \(person) as.matrix(person[vars]))
  printed <- utils::capture.output(reference <- gimme::gimme(data = series, ar = TRUE, plot = FALSE, subgroup = FALSE))
  found <- build_gimme(panel, vars = vars, id = "id", time = "beep", ar = TRUE)
  expect_reference(coefs(found), gimme_paths(reference), reference = "gimme",
                   covers = c("build_gimme()", "coefs()", "coefs.net_gimme()"), quantity = "path estimates, 10 people",
                   tolerance = 0, key = c("subject", "network", "from", "to"))
})
