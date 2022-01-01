#' Mean Absolute Error
#'
#' A somewhat naive criterion that is sometimes used consists of simply omitting all censored cases from the data set.
#' For survival analysis problems, the mean absolute error (MAE) can be defined as an average of the differences between the predicted time values and the actual observation time values.
#' Only the samples for which the event occurs are being considered in this metric.
#'
#' Condition: MAE can only be used for the evaluation of survival models which can provide the event time as the predicted target value.
#'
#' @param object object of class \code{Surv} created by Surv function.
#' @param pre_time a vector of predicted values of survival time of each observation.
#'
#' @return the value of Mean Absolute Error
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' Matsuo, K. ,  Purushotham, S. ,  Jiang, B. ,  Mandelbaum, R. S. ,  Takiuchi, T. , &  Liu, Y. , et al. (2018). Survival outcome prediction in cervical cancer: cox models vs deep-learning model. American Journal of Obstetrics & Gynecology.
#' Coyle, E. J. , &  Lin, J. H. . (1988). Stack filters and the mean absolute error criterion. IEEE Trans Acoustics Speech Signal Processing, 36(8), 1244-1254.
#' @examples
#' library(survival)
#' time = rexp(50)
#' status = sample(c(0,1),50,replace = TRUE)
#' pre_time = rexp(50)
#' MAE(Surv(time,status),pre_time)
#'
#' @importFrom survival Surv
#'
#' @export
MAE <- function(object,pre_time){

  if (TRUE %in% (is.na(object)|is.na(pre_time))) {
    stop("The input vector cannot have NA")
  }

  time = object[,1]
  status = object[,2]
  if (missing(time) | missing(status) | missing(pre_time)) {
    stop("Input value is missing")
  }

  if (length(time) != length(status)) {
    stop("The lengths of time and status are not equal")
  }

  if (length(time) != length(pre_time)) {
    stop("The lengths of time and pre_time are not equal")
  }

  status = ifelse(status == min(status), 0, 1)
  t_mae = abs(time - pre_time)*status
  mae = sum(t_mae)/sum(status)
  names(mae) = 'MAE'
  return(mae)
}
