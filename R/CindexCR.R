#' Concordance index in the Presence of Competing Risks
#'
#' The C-index (Concordance index) of the prognostic model in the presence of competing risks according to Marcel, W et al.(2014).
#'
#' @param time minimum value of deletion time and survival time.
#' @param status the status indicator, for models with competing risks, the status indicator is 0=censored, 1=event at \code{time}, 2= competing risks at \code{time}.
#' @param predicted a vector of predicted values or the survival time of survival probabilities of each observation.
#' @param Cause_int event type of interest, the default value is 1.
#'
#' @return Estimates of the C-index in the presence of competing risks.
#' @author HanPu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#'
#' Marcel, W. ,  Paul, B. ,  Koller, M. T. ,  Witteman, J. , &  Gerds, T. A. . (2014).Concordance for prognostic models with competing risks. Biostatistics(3), 526.
#'
#' Ishwaran, H. ,  Kogalur, U. B. ,  Blackstone, E. H. , &  Lauer, M. S. . (2008). Random survival forests. Journal of Thoracic Oncology Official Publication of the International Association for the Study of Lung Cancer, 2(12), 841-860.
#'
#' @examples
#' time <- c(4, 7, 5, 8)
#' status <- rep(1, 4)
#' predicted <- c(3, 5, 7, 10)
#' Cause_int <- 1
#' CindexCR(time, status, predicted, Cause_int)
#'
#' @export

CindexCR <- function(time, status, predicted, Cause_int = 1) {
  if (any(is.na(time))) {
    stop("The input vector cannot have NA")
  }
  if (any(is.na(status))) {
    stop("The input vector cannot have NA")
  }
  if (any(!(status %in% c(0, 1, 2)))) {
    stop("The status must be 0 or 1 or 2")
  }
  if (any(is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  if (!(Cause_int %in% status)) {
    stop("Invalid input of Cause_int")
  }
  if (min(time) <= 0) {
    stop("Survival time must be positive")
  }

  Time_survival <- time
  Censoring <- ifelse(status == 0, 0, 1)
  Cause <- ifelse(status == 2, 2, 1)
  Prediction <- -log(predicted)
  Time <- max(Time_survival) + 1

  n <- length(Prediction)
  A <- matrix(0, nrow = n, ncol = n)
  B <- matrix(0, nrow = n, ncol = n)
  Q <- matrix(0, nrow = n, ncol = n)
  N_t <- matrix(0, nrow = n, ncol = n)
  Num_mat <- matrix(0, nrow = n, ncol = n)
  Den_mat <- matrix(0, nrow = n, ncol = n)
  Num <- 0
  Den <- 0
  for (i in 1:n) {
    A[i, which(Time_survival[i] < Time_survival)] <- 1
    B[i, intersect(intersect(which((
      Time_survival[i] >= Time_survival
    )), which(Cause != Cause_int)), which(Censoring == 1))] <- 1
    Q[i, which(Prediction[i] > Prediction)] <- 1
  }
  for (i in 1:n) {
    if (Time_survival[i] <= Time &&
      Cause[i] == Cause_int && Censoring[i] == 1) {
      N_t[i, ] <- 1
    }
  }
  Num_mat <- (A + B) * Q * N_t
  Den_mat <- (A + B) * N_t
  Num <- sum(Num_mat)
  Den <- sum(Den_mat)
  return(Num / Den)
}
