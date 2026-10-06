# CRAN's incoming checks flag CPU time far above elapsed time; cap the
# threads that BLAS/OpenMP and data.table may spawn during the tests.
Sys.setenv(OMP_THREAD_LIMIT = 2)

library(testthat)
library(Nestimate)

test_check("Nestimate")
