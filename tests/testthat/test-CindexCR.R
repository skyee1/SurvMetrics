test_that("CindexCR", {
  time <- c(4, 7, 5, 8)
  status <- rep(1, 4)
  predicted <- c(3, 5, 7, 10)
  Cause_int <- 1
  CindexCR(time, status, predicted, Cause_int)
  expect_type(CindexCR(time, status, predicted, Cause_int), "double")
  status[1] <- 4
  expect_error(
    CindexCR(time, status, predicted, Cause_int),
    "The status must be 0 or 1 or 2"
  )
  status[1] <- 1
  expect_error(
    CindexCR(time, status, predicted, 2),
    "Invalid input of Cause_int"
  )
  time[1] <- NA
  expect_error(
    CindexCR(time, status, predicted, Cause_int),
    "The input vector cannot have NA"
  )
  time[1] <- -1
  expect_error(
    CindexCR(time, status, predicted, Cause_int),
    "Survival time must be positive"
  )
})
