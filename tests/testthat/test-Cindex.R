test_that("Cindex", {
  # case 1 Surv object
  library(survival)
  time <- c(1, 1, 2, 2, 2, 2, 2, 2)
  status <- c(0, 1, 1, 0, 1, 1, 0, 1)
  predicted <- c(2, 3, 3, 3, 4, 2, 4, 3)
  expect_type(Cindex(Surv(time, status), predicted), "double")
  expect_equal(as.numeric(Cindex(Surv(time, status), predicted)), 0.642857)
  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 1)), c(0.5, 0.6))), 1)
  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 1)), c(0.6, 0.5))), 0)
  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 1)), c(0.5, 0.5))), 0.5)
  expect_equal(as.numeric(Cindex(Surv(c(1, 1), c(1, 1)), c(0.5, 0.5))), 1)
  expect_equal(as.numeric(Cindex(Surv(c(1, 1), c(1, 1)), c(0.5, 0.6))), 0.5)

  expect_equal(as.numeric(Cindex(Surv(c(1, 1), c(1, 0)), c(0.5, 0.5))), 0.5)
  expect_equal(as.numeric(Cindex(Surv(c(1, 1), c(1, 0)), c(0.5, 0.6))), 1)
  expect_equal(as.numeric(Cindex(Surv(c(1, 1), c(1, 0)), c(0.6, 0.5))), 0.5)

  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 0)), c(0.5, 0.6))), 1)
  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 0)), c(0.5, 0.5))), 0.5)
  expect_equal(as.numeric(Cindex(Surv(c(1, 2), c(1, 0)), c(0.6, 0.5))), 0)


  expect_error(
    Cindex(Surv(time[-1], status), predicted),
    "Time and status are different lengths"
  )
  expect_error(
    Cindex(Surv(time, status), predicted[-1]),
    "The lengths of time and predicted are not equal"
  )
  time[1] <- NA
  expect_error(
    Cindex(Surv(time, status), predicted),
    "The input vector cannot have NA"
  )

  # case2 fit object
  library(randomForestSRC)
  library(pec)
  set.seed(1234)
  mydata <- kidney[, -1]
  train_index <- sample(1:nrow(mydata), 0.7 * nrow(mydata))
  train_data <- mydata[train_index, ]
  test_data <- mydata[-train_index, ]
  # test coxph
  coxfit <- coxph(Surv(time, status) ~ ., data = train_data, x = TRUE)
  expect_type(Cindex(coxfit, test_data), "double")

  # test RSF
  rsffit <- rfsrc(Surv(time, status) ~ ., data = train_data)
  expect_type(Cindex(rsffit, test_data), "double")

  # test survreg
  for (dist in c(
    "weibull",
    "exponential",
    "gaussian",
    "logistic",
    "lognormal",
    "loglogistic"
  )) {
    survregfit <- survreg(Surv(time, status) ~ .,
      dist = dist, data = train_data
    )
    expect_type(Cindex(survregfit, test_data), "double")
  }
})
