#' Predicting Survival Probabilities for a 'survreg' Object
#'
#' Function to extract survival probability predictions from \code{survreg} modeling approach.
#'
#'
#' @param object A model fitted by \code{survreg} from which to extract predicted survival probabilities
#' @param newdata
#' A data frame containing predictor variable combinations
#' for which to compute predicted survival probabilities.
#' @param time_days A vector of times in the range of the response variable,
#' We.g. times when the response is a survival object, at which to return the survival probabilities.
#'
#' @return A matrix with as many rows as NROW(newdata) and as many columns as length(time_days).
#' Each entry should be a probability and in rows the values should be decreasing.
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#'
#' @examples
#' library(survival)
#' set.seed(1234)
#' mydata <- kidney[, -1]
#' train_index <- sample(1:nrow(mydata), 0.7 * nrow(mydata))
#' train_data <- mydata[train_index, ]
#' test_data <- mydata[-train_index, ]
#' survregfit <- survreg(Surv(time, status) ~ ., dist = 'weibull', data = train_data)
#' pre_sb <- predictSurvProb2survreg(survregfit, test_data, c(10, 20))
#'
#' @importFrom stats pnorm predict
#'
#' @export
predictSurvProb2survreg <- function(object, newdata, time_days) {
  lp <- predict(object, newdata = newdata, type = "link") # vector
  B <- object$scale # fixed value
  dist <- object$dist
  surv_fun <- function(lp_i) {
    if (dist %in% c("weibull", "exponential")) {
      exp(-exp((log(time_days) - lp_i) / B))
    } else if (dist == "lognormal") {
      1 - pnorm((log(time_days) - lp_i) / B, mean = 0, sd = 1)
    } else if (dist == "gaussian") {
      1 - pnorm((time_days - lp_i) / B, mean = 0, sd = 1)
    } else if (dist == "logistic") {
      1 / (1 + exp((time_days - lp_i) / B))
    } else if (dist == "loglogistic") {
      1 / (1 + exp((log(time_days) - lp_i) / B))
    } else {
      stop("This distribution is not supported")
    }
  }

  sp <- t(sapply(lp, surv_fun))
}
