#' Gt
#'
#' G(t)=P(C>t) denote the Kaplan-Meier estimate of the censoring distribution which is used to adjust for censoring.
#' Gt is used to calculate G(t) at any timepoint you want.
#'
#' @param object object of class \code{Surv} created by Surv function.
#' @param timepoint any point in time you want to get the Kaplan–Meier estimate of the censoring.
#'
#' @return The Kaplan–Meier estimate of the censoring in (0,1).
#'
#' @author Zhou HanPu \email{zhouhanpu@csu.edu.cn}
#' @references
#' Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#' Kaplan, E. L. , &  Meier, P. . (1958). Nonparametric estimation from incomplete observations. Journal of the American Statistical Association, 53, 457-481.
#'
#' @examples
#' library(survival)
#' time = rexp(50)
#' status = sample(c(0,1),50,replace = TRUE)
#' pre_sp = runif(50)
#' timepoint = runif(1)
#' Gt(Surv(time,status),timepoint)
#'
#' @importFrom stats na.omit
#' @importFrom survival Surv
#' @importFrom survival survfit
#'
#' @export
Gt <- function(object,timepoint){

  if(!inherits(object, "Surv"))
    stop("object is not of class Surv")

  if (missing(object)){
    stop("The survival object is missing")
  }

  if(missing(timepoint)){
    stop("The time is missing with no default")
  }

  if(!(timepoint > 0)){
    stop("The timepoint must be positive")
  }

  if (TRUE%in%(is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (length(timepoint) != 1) {
    stop("Gt can only be calculated at a single time point")
  }

  if (is.na(timepoint)) {
    stop("Cannot calculate Gt at NA")
  }

  time = object[,1]
  status = object[,2]
  status0 = ifelse(status == min(status), 1, 0)
  fit = survfit(Surv(time,status0)~1)
  res.sum = survminer::surv_summary(fit)
  res.sum = na.omit(res.sum) #The last point time may have NA

  # The observed survival times include this timepoint
  if(timepoint %in% res.sum$time){
    Gvalue = res.sum[which(res.sum$time == timepoint),]$surv
  }else{
    index1=ifelse(length(which(res.sum$time < timepoint)),
                  max(which(res.sum$time < timepoint)),1)

    index2=ifelse(length(which(res.sum$time > timepoint)),
                  min(which(res.sum$time > timepoint)),length(res.sum$time))
    Gtemp = res.sum$surv
    if (index1 == index2) {
      Gvalue = Gtemp[index1]
    }else{
      Gvalue = ((time[index2]-timepoint)*Gtemp[index1]+(timepoint-time[index1])*Gtemp[index2])/(time[index2]-time[index1])

    }
  }
  #if No dead sample, the denominator is zero
  #The minimum survival probability is used instead
  if(is.na(Gvalue) | Gvalue <= 0){
    Gvalue = min(res.sum$surv[res.sum$surv != 0])
  }
  return(Gvalue)
  }
