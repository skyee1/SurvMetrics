test_that("Gt", {
  library(survival)
  time <- c(1, 1, 2, 2, 2, 2, 2, 2)
  status <- c(0, 1, 1, 0, 1, 1, 0, 1)
  timepoint <- 2

  expect_type(Gt(Surv(time, status), timepoint), "double")
  expect_error(
    Gt(NA, timepoint),
    "object is not of class Surv"
  )
  expect_error(
    Gt(Surv(time, status), -1),
    "The timepoint must be positive"
  )
  expect_error(
    Gt(Surv(time, status), c(1, 2)),
    "Gt can only be calculated at a single time point"
  )
  time[1] <- NA
  expect_error(
    Gt(Surv(time, status), timepoint),
    "The input vector cannot have NA"
  )
})
