#' The Integrate Absolute Error and The Integrate Square Error
#'
#' Two ways of the continuous-time approach to continuous-time identification based on least-squares and least-absolute errors are proposed.
#' Integrate Absolute Error and Integrate Square Error.To evaluate the performance of survival models methods
#' Lower values of IAE or ISE indicate better performances.
#'
#' @aliases IAEISE
#' @param object object of class \code{Surv} created by Surv function
#' or a fitted survival model, including the survival model
#' fitted by \code{coxph}, \code{survreg}, and \code{rfsrc}.
#' @param sp_matrix a matrix or data.frame of predicted values of survival probabilities for the testing set.
#' rows denote different samples, columns denote different time points,
#' and the values in row i and column j of the matrix denote the predicted survival probability
#' of the ith sample at the time point corresponding to the jth column.
#' @param IRange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#' Or the scale you want to get the IAE and ISE;  .
#'
#' @return Estimates of the Integrate Absolute Error and the Integrate Square Error of the predicted values of survival probabilities.
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' Marron, J. S. , &  Wand, M. P. . (1992). Exact mean integrated squared error. Annals of Statistics, 20(2), 712-736.
#'
#' HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). L1 splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#'
#' Kowalczuk, & Z. (1998). Integrated squared error and integrated absolute error in recursive identification of continuous-time plants. Control 98 Ukacc International Conference on (Vol.1998, pp.693-698). IET.
#'
#' @examples
#'
#' library(survival)
#' library(SurvMetrics)
#' set.seed(123)
#' N <- 100
#' mydata <- SDGM4(N, p = 20, c_step = -0.5)
#' index.train <- sample(1:N, 2 / 3 * N)
#' data.train <- mydata[index.train, ]
#' data.test <- mydata[-index.train, ]
#'
#' time_interest <- sort(data.train$time[data.train$status == 1])
#' sp_matrix <- matrix(sort(runif(nrow(data.test) * length(time_interest)),
#'   decreasing = TRUE
#' ), nrow = nrow(data.test))
#' object <- Surv(data.test$time, data.test$status)
#'
#' # a vector for all the distinct time
#' IAEISE(object, sp_matrix, time_interest)
#' # a range
#' IAEISE(object, sp_matrix, c(12, 350))
#'
#' @importFrom survival Surv
#' @importFrom survival survfit
#' @importFrom pec predictSurvProb
#' @importFrom randomForestSRC predict.rfsrc
#'
#' @export

IAEISE <- function(object, sp_matrix, IRange = c(-2, -1)) {
  # case1、coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- sp_matrix

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    mat_coxph <-
      predictSurvProb(obj, test_data, distime) # get the survival probability matrix
    object_coxph <- Surv(test_data$time, test_data$status)

    object <- object_coxph
    sp_matrix <- mat_coxph
    if (max(IRange) <= 0) {
      IRange <- range(object[, 1])
    } # the fixed time point
  }


  # case2、RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- sp_matrix

    mat_rsf <-
      predict(obj, test_data)$survival # get the survival probability matrix
    object_rsf <- Surv(test_data$time, test_data$status)

    object <- object_rsf
    sp_matrix <- mat_rsf
    if (max(IRange) <= 0) {
      IRange <- range(object[, 1])
    } # the fixed time point
  }

  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- sp_matrix

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    object <- Surv(test_data$time, test_data$status)
    sp_matrix <- predictSurvProb2survreg(obj, test_data, distime)
    if (max(IRange) <= 0) {
      IRange <- range(object[, 1])
    } # the fixed time point
  }

  if (any(is.na(IRange))) {
    stop("Cannot calculate IAE or ISE in the interval containing NA")
  }

  if (!is.numeric(IRange)) {
    stop("The class of the IRange must be numeric! or use the default setting")
  }

  # default time range for IAEISE()
  if (max(IRange) <= 0) {
    IRange <- range(object[, 1])
  }

  if (!inherits(object, "Surv")) {
    stop("object is not of class Surv")
  }

  if (any(is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (any(is.na(sp_matrix))) {
    stop("The input probability matrix cannot have NA")
  }

  if (length(IRange) <= 1) {
    stop("Can not calculate the integration at a single point")
  }

  if (any(IRange <= 0)) {
    stop("The interval value must be positive")
  }

  if (any(diff(IRange) <= 0)) {
    stop("The interval must increase")
  }

  # Simulation to generate discrete time series
  sp_matrix2sp <- apply(sp_matrix, 2, mean)
  if (length(IRange) != ncol(sp_matrix)) {
    t_IRange <- range(IRange)
    p <- ncol(sp_matrix)
    IRange <- seq(t_IRange[1], t_IRange[2], length = p)
  }


  fit0 <- survfit(object ~ 1)
  res_sum0 <- survminer::surv_summary(fit0)
  KM_frame0 <- data.frame(
    "time" = res_sum0$time,
    "sp0" = res_sum0$surv,
    "sp1" = 0
  )[res_sum0$n.event != 0, ]


  KM_frame1 <- data.frame(
    "time" = IRange,
    "sp0" = sp_matrix2sp,
    "sp1" = 0
  )

  for (i in 1:nrow(KM_frame0)) {
    t <- KM_frame0$time[i]
    if (t < min(KM_frame1$time)) {
      KM_frame0$sp1[i] <- 1
    } else {
      index <- max(which(KM_frame1$time <= t))
      KM_frame0$sp1[i] <- KM_frame1$sp0[index]
    }
  }

  KM_frame0$t_IAE <- abs(KM_frame0$sp0 - KM_frame0$sp1)
  KM_frame0$t_ISE <- (KM_frame0$sp0 - KM_frame0$sp1)^2
  IAE <- 0
  ISE <- 0
  for (i in 1:(nrow(KM_frame0) - 1)) {
    IAE <- IAE + (KM_frame0$time[i + 1] - KM_frame0$time[i]) * KM_frame0$t_IAE[i]
    ISE <- ISE + (KM_frame0$time[i + 1] - KM_frame0$time[i]) * KM_frame0$t_ISE[i]
  }
  IAEISE_value <- c(IAE, ISE)
  names(IAEISE_value) <- c("IAE", "ISE")
  return(round(IAEISE_value, 4))
}
