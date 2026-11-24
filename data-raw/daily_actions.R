# daily_actions: a simulated log of learning actions over seven days.
# Same generator and seed as Tutorial_docs/variable-length-sequence-plots.Rmd,
# so the bundled data reproduce that tutorial. Run from the package root:
#   Rscript data-raw/daily_actions.R
set.seed(2026)
states <- c("plan", "read", "write", "discuss", "review")
students <- sprintf("s%02d", seq_len(30))
daily_actions <- do.call(rbind, lapply(students, \(s) {
  n <- sample(50:100, 1)
  # every day once, the rest spread unevenly over the days
  day <- sort(c(1:7, sample(1:7, n - 7, replace = TRUE,
                            prob = stats::rgamma(7, shape = 2))))
  data.frame(student = s, day = paste("day", day),
             action = sample(states, n, replace = TRUE,
                             prob = stats::rgamma(5, shape = 0.7)),
             stringsAsFactors = FALSE)
}))
rownames(daily_actions) <- NULL

stopifnot(
  "30 students" = length(unique(daily_actions$student)) == 30L,
  "7 days" = length(unique(daily_actions$day)) == 7L,
  "every student on every day" = all(table(daily_actions$student,
                                           daily_actions$day) >= 1L),
  "no missing values" = !anyNA(daily_actions)
)
save(daily_actions, file = file.path("data", "daily_actions.rda"),
     compress = "xz")
