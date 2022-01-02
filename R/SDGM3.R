#' SDGM3
#'
#' Survival data generation method. The proportional hazards assumption is strongly violated in this dataset.
#'
#' @param N  The sample size of the simulated dataset.
#' @param p  The covariate dimension of the simulated dataset.
#' @param u_max  The parameter which is used to control the censoring rate.
#'
#' @return the simulated dataset
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#' Steingrimsson, J. A. ,  Diao, L. , &  Strawderman, R. L. . (2019). Censoring unbiased regression trees and ensembles. Journal of the American Statistical Association, 114.
#'
#' Zhu, R. , &  Kosorok, M. R. . (2012). Recursively imputed survival trees. Journal of the American Statistical Association, 107(497), 331-340.
#'
#' Ishwaran, H. ,  Kogalur, U. B. ,  Gorodeski, E. Z. ,  Minn, A. J. , &  Lauer, M. S. . (2010). High-dimensional variable selection for survival data. Journal of the American Statistical Association, 105(489), 205-217.
#'
#' @examples
#' SDGM3(N=200,p = 15,u_max = 7)
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats quantile
#' @importFrom stats rgamma
#' @importFrom stats runif
#'
#' @export
#'
SDGM3 <- function(N=200,p = 15,u_max = 7)
{
  finalcenrate = rep(0,20)
  for (repnum in 1:20) {
  mu = rep(0,p)
  Si = matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      Si[i,j] = 0.75^abs(i-j)
    }
  }
  W = mvrnorm(N,mu,Si)
  Ti = rep(0,N)
  q <- floor(quantile(1:p,c(2/5,3/5)))
  for (i in 1:N) {
    shape = 0.5 + 0.3 * abs(sum(W[i,(q[1]+1):q[2]]))
    scale = 2
    Ti[i] = rgamma(1,shape = shape,scale = scale)
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

