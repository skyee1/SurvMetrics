#' The Concordance Index for Right-Censored Survival Time Data
#'
#' Concordance index is a rank correlation measures between a variable X and a possibly
#' censored variable Y, with event/censoring indicator. In survival analysis, a pair of patients is called concordant if the risk of
#' the event predicted by a model is lower for the patient who experiences the
#' event at a later timepoint. The concordance probability (C-index) is the
#' frequency of concordant pairs among all pairs of subjects. It can be used to
#' measure and compare the discriminative power of a risk prediction models.
#'
#' Pairs with identical observed times, where one is uncensored and one is
#' censored, are always considered usuable (independent of the value of
#' \code{tiedOutcomeIn}), as it can be assumed that the event occurs at a later
#' timepoint for the censored observation.
#'
#' For uncensored response the result equals the one obtained with the
#' functions \code{rcorr.cens} and \code{rcorrcens} from the \code{Hmisc}
#' package (see examples).
#'
#' @param object object of class \code{Surv} created by Surv function
#' or a fitted survival model, including the survival model
#' fitted by \code{coxph}, \code{survreg}, and \code{rfsrc}.
#' @param predicted
#' If the input of the parameter \code{object} is a fitted survival model,
#' this parameter should be a survival dataset on which you want to caculate the C-index of the fitted model.
#' Or with an object of class \code{Surv} as the parameter \code{object},
#' it should be a vector containing the predicted survival time or probability of each observation.
#' @param t_star the timepoint at which the C-index you want to calculate.
#' If the input of the parameter \code{object} is a fitted survival model, the timepoint is necessary to be specified
#' at which the survival probability is predicted,
#' and this function will calculate the C-index between the survival time and the predicted result at that moment.
#' If the input object is a survival object, this parameter can be ignored
#' and the value of this parameter will not have any effect on the result of this fuction.
#'
#' @return Estimates of the C-index between the survival time and the predicted result.
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#'
#' Ishwaran, H. ,  Kogalur, U. B. ,  Blackstone, E. H. , &  Lauer, M. S. . (2008). Random survival forests. Journal of Thoracic Oncology Official Publication of the International Association for the Study of Lung Cancer, 2(12), 841-860.
#'
#' Kang, L. ,  Chen, W. ,  Petrick, N. A. , &  Gallas, B. D. . (2015). Comparing two correlated c indices with right-censored survival outcome: a one-shot nonparametric approach. Statistics in Medicine, 34(4).
#'
#' TA Gerds, MW Kattan, M Schumacher, and C Yu. Estimating a time-dependent
#' concordance index for survival prediction models with covariate dependent
#' censoring. Statistics in Medicine, Ahead of print:to appear, 2013. DOI =
#' 10.1002/sim.5681
#'
#' Wolbers, M and Koller, MT and Witteman, JCM and Gerds, TA (2013) Concordance
#' for prognostic models with competing risks Research report 13/3. Department
#' of Biostatistics, University of Copenhagen
#'
#' Andersen, PK (2012) A note on the decomposition of number of life years lost
#' according to causes of death Research report 12/2. Department of
#' Biostatistics, University of Copenhagen
#'
#' Paul Blanche, Michael W Kattan, and Thomas A Gerds. The c-index is not
#' proper for the evaluation of-year predicted risks. Biostatistics, 20(2):
#' 347--357, 2018.
#'
#' @examples
#' library(survival)
#' time <- c(1, 1, 2, 2, 2, 2, 2, 2)
#' status <- c(0, 1, 1, 0, 1, 1, 0, 1)
#' predicted <- c(2, 3, 3, 3, 4, 2, 4, 3)
#' Cindex(Surv(time, status), predicted)
#'
#' @importFrom survival Surv
#' @importFrom pec predictSurvProb
#' @importFrom randomForestSRC predict.rfsrc
#' @importFrom stats median
#'
#' @export
Cindex <- function(object, predicted, t_star = -1) {
  # case1、coxph AND testing set
  if (inherits(object, "coxph")) {
    obj <- object
    test_data <- predicted
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
    predicted <- vec_coxph
  }


  # case2、RSF AND testing set
  if (inherits(object, c("rfsrc"))) {
    obj <- object
    test_data <- predicted
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
    predicted <- vec_rsf
  }


  # case3 survreg AND testing set
  if (inherits(object, c("survreg"))) {
    obj <- object
    test_data <- predicted
    t_star0 <- t_star

    # the interesting times of training set
    distime <- sort(unique(as.vector(obj$y[obj$y[, 2] == 1])))

    if (t_star0 <= 0) {
      t_star0 <- median(distime)
    } # the fixed time point

    predicted <- predictSurvProb2survreg(obj, test_data, t_star0)

    object <- Surv(test_data$time, test_data$status)
  }

  time <- object[, 1]
  status <- object[, 2]

  if (length(time) != length(status)) {
    stop("The lengths of time and status are not equal")
  }

  if (length(time) != length(predicted)) {
    stop("The lengths of time and predicted are not equal")
  }

  if (any(is.na(time) | is.na(status) | is.na(predicted))) {
    stop("The input vector cannot have NA")
  }

  # initialization
  permissible <- 0 # comparable pairs
  concord <- 0 # completely concordance
  par_concord <- 0 # partial concordance

  n <- length(time)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Exclude incomparable pairs
      if ((time[i] < time[j] &
        status[i] == 0) | (time[j] < time[i] & status[j] == 0)) {
        next
      }

      if (time[i] == time[j] & status[i] == 0 & status[j] == 0) {
        next
      }

      permissible <- permissible + 1

      # survival times unequal
      if (time[i] != time[j]) {
        if ((time[i] < time[j] &
          predicted[i] < predicted[j]) |
          (time[j] < time[i] & predicted[j] < predicted[i])) {
          concord <- concord + 1
        } else if (predicted[i] == predicted[j]) {
          par_concord <- par_concord + 0.5
        }
      }

      # survival times equal
      if (time[i] == time[j] & status[i] == 1 & status[j] == 1) {
        if (predicted[i] == predicted[j]) {
          concord <- concord + 1
        } else {
          par_concord <- par_concord + 0.5
        }
      }
      # one censored one died
      if (time[i] == time[j] &
        ((status[i] == 1 &
          status[j] == 0) | (status[i] == 0 & status[j] == 1))) {
        if ((status[i] == 1 &
          predicted[i] < predicted[j]) |
          (status[j] == 1 & predicted[j] < predicted[i])) {
          concord <- concord + 1
        } else {
          par_concord <- par_concord + 0.5
        }
      }
    }
  }
  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}
