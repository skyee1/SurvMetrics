#' Survival Data Generation Method 1
#'
#' Survival data generation method. An example of the proportional hazards model where in the Cox model is expected to perform best.
#'
#' @param N  The sample size of the simulated dataset.
#' @param p  The covariate dimension of the simulated dataset.
#' @param c_mean  The parameter which is used to control the censoring rate.
#'
#' @return the simulated dataset
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#'
#' Steingrimsson, J. A. ,  Diao, L. , &  Strawderman, R. L. . (2019). Censoring unbiased regression trees and ensembles. Journal of the American Statistical Association, 114.
#'
#' Zhu, R. , &  Kosorok, M. R. . (2012). Recursively imputed survival trees. Journal of the American Statistical Association, 107(497), 331-340.
#'
#' Ishwaran, H. ,  Kogalur, U. B. ,  Gorodeski, E. Z. ,  Minn, A. J. , &  Lauer, M. S. . (2010). High-dimensional variable selection for survival data. Journal of the American Statistical Association, 105(489), 205-217.
#'
#' @examples
#' SDGM1(N = 200, p = 15, c_mean = 0.4)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rexp
#'
#' @export
#'
SDGM1 <- function(N = 200, p = 15, c_mean = 0.4) {
  finalcenrate <- rep(0, 20)
  for (repnum in 1:20) {
    mu <- rep(0, p)
    Si <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        Si[i, j] <- 0.9^abs(i - j)
      }
    }
    W <- mvrnorm(N, mu, Si)
    Ti <- rep(0, N)
    for (i in 1:N) {
      t.mu <- exp(0.1 * sum(W[i, (floor(p / 2) + 1):p]))
      Ti[i] <- rexp(1, 1 / t.mu)
    }

    c.time <- rexp(N, c_mean)
    mydata.x <- as.data.frame(W)
    mydata.time <- data.frame("time" = Ti, "c.time" = c.time, "status" = 0)
    mydata.time$status <- ifelse(mydata.time$time < mydata.time$c.time, 1, 0)
    mydata.time$time <- apply(mydata.time[, 1:2], 1, min)

    mydata0 <- data.frame("time" = mydata.time$time, "status" = mydata.time$status, W)
    cen.rate <- 1 - sum(mydata0$status) / N
    finalcenrate[repnum] <- cen.rate
  }
  print(paste("censoring rate is:", round(mean(finalcenrate), 3)))
  return(mydata0)
}
