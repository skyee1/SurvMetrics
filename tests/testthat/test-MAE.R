test_that("MAE", {
  library(survival)
  time <- c(1, 1, 2, 2, 2, 2, 2, 2)
  status <- c(0, 1, 1, 0, 1, 1, 0, 1)
  predicted <- c(2, 3, 3, 3, 4, 2, 4, 3)

  expect_type(MAE(Surv(time, status), predicted), "double")
  expect_equal(as.numeric(MAE(Surv(time, status), predicted)), 1.2)
  expect_error(
    MAE(Surv(time[-1], status), predicted),
    "Time and status are different lengths"
  )
  expect_error(
    MAE(Surv(time, status), predicted[-1]),
    "The lengths of time and pre_time are not equal"
  )
  time[1] <- NA
  expect_error(
    MAE(Surv(time, status), predicted),
    "The input vector cannot have NA"
  )
})
