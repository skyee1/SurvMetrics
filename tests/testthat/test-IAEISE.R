test_that("IAEISE", {
  # case 1 Surv object
  library(survival)
  library(randomForestSRC)
  library(pec)
  set.seed(1234)
  mydata <- kidney[, -1]
  train_index <- sample(1:nrow(mydata), 0.7 * nrow(mydata))
  train_data <- mydata[train_index, ]
  test_data <- mydata[-train_index, ]
  coxfit <- coxph(Surv(time, status) ~ ., data = train_data, x = TRUE)
  distime <- sort(unique(as.vector(coxfit$y[coxfit$y[, 2] == 1])))
  sp_matrix <- predictSurvProb(coxfit, test_data, distime)
  time <- test_data$time
  status <- test_data$status

  expect_type(IAEISE(Surv(time, status), sp_matrix), "double")
  expect_equal(as.numeric(IAEISE(Surv(time, status), sp_matrix)), c(132.0101, 44.3022))

  expect_error(
    IAEISE(Surv(time[-1], status), sp_matrix),
    "Time and status are different lengths"
  )
  expect_error(
    IAEISE(Surv(time, status), sp_matrix, c(2, 1)),
    "The interval must increase"
  )
  expect_error(
    IAEISE(Surv(time, status), sp_matrix, c(1, NA)),
    "Cannot calculate IAE or ISE in the interval containing NA"
  )
  expect_error(
    IAEISE(Surv(time, status), NA),
    "The input probability matrix cannot have NA"
  )
  expect_error(
    IAEISE(Surv(time, status), sp_matrix, 1),
    "Can not calculate the integration at a single point"
  )
  time[1] <- NA
  expect_error(
    IAEISE(Surv(time, status), sp_matrix),
    "The input vector cannot have NA"
  )
  # case2 fit object
  # test coxph
  coxfit <- coxph(Surv(time, status) ~ ., data = train_data, x = TRUE)
  expect_type(IAEISE(coxfit, test_data), "double")
  # test RSF
  rsffit <- rfsrc(Surv(time, status) ~ ., data = train_data)
  expect_type(IAEISE(rsffit, test_data), "double")
  # test survreg
  for (dist in c(
    "weibull", "exponential", "gaussian",
    "logistic", "lognormal", "loglogistic"
  )) {
    survregfit <- survreg(Surv(time, status) ~ .,
      dist = dist, data = train_data
    )
    expect_type(IAEISE(survregfit, test_data), "double")
  }
})
