#' SDGM2
#'
#' Survival data generation method. The dataset represents mild violations of the proportional hazards assumption.
#'
#' @param N  The sample size of the simulated dataset.
#' @param p  The covariate dimension of the simulated dataset.
#' @param u_max  The parameter which is used to control the censoring rate.
#'
#' @return the simulated dataset
#'
#' @author Zhou HanPu \email{zhouhanpu@csu.edu.cn}
#' @references
#' Steingrimsson, J. A. ,  Diao, L. , &  Strawderman, R. L. . (2019). Censoring unbiased regression trees and ensembles. Journal of the American Statistical Association, 114.
#'
#' Zhu, R. , &  Kosorok, M. R. . (2012). Recursively imputed survival trees. Journal of the American Statistical Association, 107(497), 331-340.
#'
#' Ishwaran, H. ,  Kogalur, U. B. ,  Gorodeski, E. Z. ,  Minn, A. J. , &  Lauer, M. S. . (2010). High-dimensional variable selection for survival data. Journal of the American Statistical Association, 105(489), 205-217.
#'
#' @examples
#' SDGM2(N=200,p = 15,u_max = 4)
#'
#'
#' @importFrom stats rexp
#' @importFrom stats runif
#'
#' @export
#'
SDGM2 <- function(N=200,p = 15,u_max = 4){
  finalcenrate = rep(0,20)
  for (repnum in 1:20) {
  W = matrix(0,N,p)
  for (i in 1:p) {
    W[,i] <- runif(N)
  }
  Ti = rep(0,N)
  for (i in 1:N) {
    t.mu = sin(W[i,1]*pi)+2*abs(W[i,2]-0.5)+W[i,3]^3
    Ti[i] = rexp(1,1/t.mu)
  }

  c.time = runif(N,0,u_max)
  mydata.x = as.data.frame(W)
  mydata.time = data.frame('time'=Ti,'c.time' = c.time,'status' = 0)
  mydata.time$status = ifelse(mydata.time$time < mydata.time$c.time,1,0)
  mydata.time$time = apply(mydata.time[,1:2],1,min)

  mydata0 = data.frame('time' = mydata.time$time,'status' = mydata.time$status,W)
  cen.rate = 1 - sum(mydata0$status)/N

  finalcenrate[repnum] = cen.rate
  }

  print(paste("censoring rate is:",round(mean(finalcenrate),3)))
  return(mydata0)
}
