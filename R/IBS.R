#' The Integration of the Brier Score
#'
#' IBS is  an integrated version of the Brier which is used to calculate the integration of the Brier Score.
#' The Brier Score is the mean square difference between the true classes and the predicted probabilities.
#' Basically, the IBS is an integrated
#' weighted squared distance between the estimated survival function and the empirical
#' survival function. The inverse probability censoring weighting(IPCW) is used to adjust for censoring.
#'
#' The percentage of censored observations increases in time, and this will surely affect the dispersion of the empirical Brier Score.
#' The question of how censoring in finite samples acts on the distribution of our measures of inaccuracy is an interesting subject.
#' Our recommendation is to choose t* in a way that censoring is not too heavy (for example, the median follow-up time).
#' We also prefer measures with integrated loss functions since they will reflect inaccuracy over an interval rather than just at one point in time.
#' In addition, the corresponding empirical measures are likely to have lower dispersion, because censored observations contribute their estimated event-free probabilities to the integrand until the censoring occurs.
#'
#' @param object object of class \code{Surv} created by Surv function
#' or a fitted survival model, including the survival model
#' fitted by \code{coxph}, \code{survreg}, and \code{rfsrc}.
#' @param sp_matrix a matrix or data.frame of predicted values of survival probabilities for the testing set.
#' Rows denote different samples, columns denote different time points,
#' and the values in entry (i,j) of the matrix denote the predicted survival probability
#' of the ith sample at the time point corresponding to the jth column.
#' @param IBSrange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#' Or the scale you want to get the IBS; and if it is a single point the return value will be the Brier Score at the timepoint.
#'
#' @return The integration of the Brier score of the predicted values of survival probabilities
#' on the discrete time points or the time scale of interest to users.
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). \\(l_1\\) splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#'
#' Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#' Brier, G. W. . (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.
#'
#' Gneiting, T. , &  Raftery, A. E. . (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.

#' @examples
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
#' # the default time points
#' IBS(object, sp_matrix, time_interest)
#' # a time range
#' IBS(object, sp_matrix, c(18:100))
#'
#' @importFrom survival Surv
#' @importFrom pec predictSurvProb
#' @importFrom randomForestSRC predict.rfsrc
#'
#' @export
IBS <- function(object, sp_matrix, IBSrange = c(-2, -1)) {
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
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
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
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
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
    if (max(IBSrange) <= 0) {
      IBSrange <- range(object[, 1])
    } # the fixed time point
  }

  if (!is.numeric(IBSrange)) {
    stop("The class of the IBSrange must be numeric! or use the default setting")
  }

  if (any(is.na(IBSrange))) {
    stop("Cannot calculate IBS in the interval containing NA")
  }

  # default time range for IBS()
  if (max(IBSrange) <= 0) {
    IBSrange <- range(object[, 1])
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

  if (length(object) != nrow(sp_matrix)) {
    stop(
      "number of rows of the sp_matrix and the survival object have different lengths"
    )
  }

  if (any(IBSrange <= 0)) {
    stop("The integration interval must be positive")
  }

  if (any(diff(IBSrange) <= 0)) {
    stop("The integral interval value must increase")
  }

  # Determine the input dimension
  if (ncol(data.frame("sp" = sp_matrix)) == 1 &
    length(IBSrange) == 1) {
    bs <- Brier(object, sp_matrix, IBSrange)
    names(bs) <- "Brier Score"
    return(round(bs, 6))
  } else if ((ncol(data.frame("sp" = sp_matrix)) - 1) * (length(IBSrange) -
    1) == 0) {
    stop("The length is illegal")
  } else if (ncol(sp_matrix) == length(IBSrange)) {
    IBSrange <- sort(IBSrange)
    t_brier <- rep(0)
    for (i in 1:length(IBSrange)) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      t_brier[i] <- Brier(object, pre_sp, t_star)
    }
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    return(round(t_IBS, 6))
  } else {
    t_brier <- rep(0)
    t_IBSrange <- range(IBSrange)
    p <- ncol(sp_matrix)
    # Simulation to generate discrete time series
    IBSrange <- seq(t_IBSrange[1], t_IBSrange[2], length = p)
    for (i in 1:length(IBSrange)) {
      pre_sp <- sp_matrix[, i]
      t_star <- IBSrange[i]
      t_brier[i] <- Brier(object, pre_sp, t_star)
    }
    t_IBS <- 0
    for (i in 1:(length(IBSrange) - 1)) {
      t_IBS <- t_IBS + (IBSrange[i + 1] - IBSrange[i]) * t_brier[i]
    }
    t_IBS <- t_IBS / (range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) <- "IBS"
    return(round(t_IBS, 6))
  }
}
