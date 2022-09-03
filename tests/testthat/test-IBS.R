test_that("IBS", {
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

  expect_type(IBS(Surv(time, status), sp_matrix), "double")
  expect_equal(as.numeric(IBS(Surv(time, status), sp_matrix)), 0.177957)

  expect_error(
    IBS(Surv(time[-1], status), sp_matrix),
    "Time and status are different lengths"
  )
  expect_error(
    IBS(Surv(time, status), sp_matrix, NA),
    "The class of the IBSrange must be numeric! or use the default setting"
  )
  expect_error(
    IBS(Surv(time, status), sp_matrix, c(1, NA)),
    "Cannot calculate IBS in the interval containing NA"
  )
  expect_error(
    IBS(Surv(time, status), sp_matrix, c(1, 1)),
    "The integral interval value must increase"
  )
  expect_error(
    IBS(Surv(time, status), sp_matrix, c(-1, 1)),
    "The integration interval must be positive"
  )
  expect_error(
    IBS(Surv(time, status), NA),
    "The input probability matrix cannot have NA"
  )
  expect_error(
    IBS(Surv(time, status), sp_matrix[-1, ]),
    "number of rows of the sp_matrix and the survival object have different lengths"
  )
  time[1] <- NA
  expect_error(
    IBS(Surv(time, status), sp_matrix),
    "The input vector cannot have NA"
  )
  # case2 fit object
  # test coxph
  coxfit <- coxph(Surv(time, status) ~ ., data = train_data, x = TRUE)
  expect_type(IBS(coxfit, test_data), "double")

  # test RSF
  rsffit <- rfsrc(Surv(time, status) ~ ., data = train_data)
  expect_type(IBS(rsffit, test_data), "double")

  # test survreg
  for (dist in c(
    "weibull", "exponential", "gaussian",
    "logistic", "lognormal", "loglogistic"
  )) {
    survregfit <- survreg(Surv(time, status) ~ .,
      dist = dist, data = train_data
    )
    expect_type(IBS(survregfit, test_data), "double")
  }
})
