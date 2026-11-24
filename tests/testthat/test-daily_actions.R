test_that("daily_actions holds 30 students over seven days with uneven volume", {
  expect_s3_class(daily_actions, "data.frame")
  expect_identical(names(daily_actions), c("student", "day", "action"))
  expect_equal(nrow(daily_actions), 2280L)
  expect_equal(length(unique(daily_actions$student)), 30L)
  expect_setequal(unique(daily_actions$day), paste("day", 1:7))
  expect_setequal(unique(daily_actions$action),
                  c("plan", "read", "write", "discuss", "review"))
  # every student acts on every day, between 1 and 41 times
  per_day <- table(daily_actions$student, daily_actions$day)
  expect_equal(range(per_day), c(1, 41))
  expect_false(anyNA(daily_actions))
})
