#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @importFrom stats cor var sd median quantile setNames aggregate
#'   approx spline cmdscale coef confint cov fitted lm na.omit
#'   p.adjust pnorm predict pt qnorm qt residuals rnorm runif
#'   rbinom rpois reshape symnum ecdf complete.cases cov2cor
#' @importFrom utils head tail combn modifyList
#' @importFrom glasso glasso
#' @importFrom data.table data.table setDT rbindlist
## usethis namespace: end
NULL

.onLoad <- function(libname, pkgname) {
  .register_builtin_estimators() # nocov
}

# Global variable declarations for R CMD check
utils::globalVariables(c(
  "from", "to", "weight", "value", "name", "variable",
  "edge", "node", "state", "time", "onset", "terminus",
  "lower", "upper", "significant", "p_value", "ci_lower", "ci_upper",
  "metric", "centrality", "vertex", "group",
  # data.table NSE column names
  ".orig_row", ".seq_grp", ".grp_key",
  # permutation_test data.table columns
  "weight_x", "weight_y",
  # mmm ggplot2 aes variables
  "max_posterior", "cluster", "k", "criterion",
  # simplicial plot aes variables
  "threshold", "betti", "dim_label", "birth", "death",
  "persistence", "components", "count", "total"
))
