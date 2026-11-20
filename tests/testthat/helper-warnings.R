# A verb warns once for each network it reads (a group, an mcml layer, a
# cluster), so expect_warning(), which catches only the first, lets the rest
# through. This expects at least one warning of `class`, absorbs every one of
# them, and returns the value of `expr`.
expect_warnings_of_class <- function(expr, class) {
  caught <- 0L
  value <- withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, class)) {
      caught <<- caught + 1L
      invokeRestart("muffleWarning")
    }
  })
  testthat::expect(caught > 0L,
                   sprintf("Expected a warning of class `%s`; none was raised.",
                           class))
  invisible(value)
}
