test_that("Brier", {
  # case 1 Surv object
  library(survival)
  time <- c(1, 1, 2, 2, 2, 2, 2, 2)
  status <- c(0, 1, 1, 0, 1, 1, 0, 1)
  pre_sp <- c(0.3, 0.2, 0.4, 0.5, 0.9, 0.1, 0.2, 0.7)
  expect_type(Brier(Surv(time, status), pre_sp), "double")
  expect_equal(as.numeric(Brier(Surv(time, status), pre_sp)), 0.468571)

  expect_error(
    Brier(Surv(time[-1], status), pre_sp),
    "Time and status are different lengths"
  )
  expect_error(
    Brier(Surv(time, status), pre_sp, NA),
    "Cannot calculate Brier Score at NA"
  )
  expect_error(
    Brier(Surv(time, status), NA),
    "The input probability vector cannot have NA"
  )
  expect_error(
    Brier(Surv(time, status), pre_sp[-1]),
    "The prediction survival probability and the survival object have different lengths"
  )
  time[1] <- NA
  expect_error(
    Brier(Surv(time, status), pre_sp),
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
  expect_type(Brier(coxfit, test_data), "double")

  # test RSF
  rsffit <- rfsrc(Surv(time, status) ~ ., data = train_data)
  expect_type(Brier(rsffit, test_data), "double")

  # test survreg
  for (dist in c(
    "weibull", "exponential", "gaussian",
    "logistic", "lognormal", "loglogistic"
  )) {
    survregfit <- survreg(Surv(time, status) ~ .,
      dist = dist, data = train_data
    )
    expect_type(Brier(survregfit, test_data), "double")
  }
})
