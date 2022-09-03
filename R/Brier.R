#' The Brier Score
#'
#' The Brier Score was proposed by Glenn W. Brier in 1950 which is a proper score function that measures the accuracy of probabilistic predictions, usually used to measure the accuracy of a model fit for survival data.
#' Brier can calculate the value of Brier Score at any timepoint,
#' regardless of whether it is the event time.
#'
#' The Brier Score is the mean square difference between the true classes and the predicted probabilities. So the Brier Score can be thought of as a cost function.
#' Therefore, the lower the Brier Score is for a set of predictions, the better the predictions are calibrated.
#' The Brier Score takes on a value between zero and one, since this is the square of the largest possible difference between a predicted probability and the actual outcome.
#' As we all know, for the cencoring samples, we do not know the real time of death, so the residual cannot be directly calculated when making the prediction.
#' So the Brier Score is widely used in survival analysis.
#'
#' The Brier Score is a strictly proper score (Gneiting and Raftery, 2007), which means that it takes its minimal value only when the predicted probabilities match the empirical probabilities.
#'
#' Judging from the sparse empirical evidence, predictions of duration of survival tend to be rather inaccurate. More precision is achieved by using patient-specific survival probabilities and the Brier score as predictions to discriminate future survivors from failures.
#' @param object object of class \code{Surv} created by Surv function
#' or a fitted survival model, including the survival model
#' fitted by \code{coxph}, \code{survreg}, and \code{rfsrc}.
#' @param pre_sp
#' If the input of the parameter \code{object} is a fitted survival model,
#' this parameter should be a survival dataset on which
#' you want to caculate the Brier Score of the fitted model.
#' Or with an object of class \code{Surv} as the parameter \code{object},
#' it should be a vector of predicted values of survival probabilities
#' of each observation in testing set at time t_star.
#' @param t_star the timepoint at which the Brier Score you want to calculate.
#' If the input of the parameter \code{object} is a fitted survival model, the timepoint is necessary to be specified
#' at which the survival probability is predicted,
#' and this function will calculate the Brier Score at that moment.
#' If the input object is a survival object, this parameter can be ignored
#' and the value of this parameter will not have any effect on the result of this fuction.
#'
#' @return the Brier Score at time t_star between the true classes and the predicted probabilities.
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#' Brier, G. W. (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.
#'
#' Gneiting, T. , &  Raftery, A. E. . (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.
#' @examples
#' library(survival)
#' time <- rexp(50)
#' status <- sample(c(0, 1), 50, replace = TRUE)
#' pre_sp <- runif(50)
#' t_star <- runif(1)
#' Brier(Surv(time, status), pre_sp, t_star)
#'
#' @importFrom survival Surv
#' @importFrom pec predictSurvProb
#' @importFrom randomForestSRC predict.rfsrc
#' @importFrom stats median
#'
#' @export
Brier <- function(object, pre_sp, t_star = -1) {
  # case1、coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point

    vec_coxph <-
      predictSurvProb(obj, test_data, t_star0) # get the survival probability vector
    object_coxph <- Surv(test_data$time, test_data$status)

    object <- object_coxph
    pre_sp <- vec_coxph
    t_star <- t_star0
  }


  # case2、RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star

    # the interesting times of training set
    distime <- obj$time.interest

    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point

    med_index <- order(abs(distime - t_star0))[1]

    mat_rsf <-
      predict(obj, test_data)$survival # get the survival probability matrix
    vec_rsf <-
      mat_rsf[, med_index] # get the survival probability vector
    object_rsf <- Surv(test_data$time, test_data$status)

    object <- object_rsf
    pre_sp <- vec_rsf
    t_star <- t_star0
  }

  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- pre_sp
    t_star0 <- t_star

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point

    pre_sp <- predictSurvProb2survreg(obj, test_data, t_star0)
    object <- Surv(test_data$time, test_data$status)
    t_star <- t_star0
  }

  if (is.na(t_star)) {
    stop("Cannot calculate Brier Score at NA")
  }

  # default time point for Brier(t_star)
  if (t_star <= 0) {
    t_star <- median(object[, 1][object[, 2] == 1])
  }

  if (length(t_star) != 1) {
    stop("Brier Score can only be calculated at a single time point")
  }

  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }

  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (any(is.na(pre_sp))) {
    stop("The input probability vector cannot have NA")
  }

  if (length(object) != length(pre_sp)) {
    stop("The prediction survival probability and the survival object have different lengths")
  }

  time <- object[, 1]
  status <- object[, 2]

  t_order <- order(time)
  time <- time[t_order]
  status <- status[t_order]
  pre_sp <- pre_sp[t_order]

  # Initialization
  sum_before_t <- 0
  sum_after_t <- 0

  Gtstar <- Gt(object, t_star)
  for (i in c(1:length(time))) {
    # survival time is less than t_star and sample died
    if (time[i] < t_star & (status[i] == 1)) {
      Gti <- Gt(Surv(time, status), time[i])
      if (is.na(Gti)) {
        next
      }
      sum_before_t <- sum_before_t + 1 / Gti * (pre_sp[i]) ^ 2 # IPCW
      next
    }
    # survival time is greater than t_star
    if (time[i] >= t_star) {
      if (is.na(Gtstar)) {
        next
      }
      sum_after_t <-
        sum_after_t + 1 / Gt(Surv(time, status), t_star) * (1 -
                                                              pre_sp[i]) ^ 2
    } # IPCW
  }

  BSvalue <- (sum_before_t + sum_after_t) / length(time)
  names(BSvalue) <- "Brier Score"

  return(round(BSvalue, 6))
}
