#' IAEISE
#'
#' Two ways of the continuous-time approach to continuous-time identification based on least-squares and least-absolute errors are proposed.
#' Integrate Absolute Error and Integrate Square Error.To evaluate the performance of survival models methods
#' Lower values of IAE or ISE indicate better performances.
#'
#' @aliases IAEISE
#' @param object object of class \code{Surv} on the testing set created by Surv function.
#' @param sp_matrix a matrix of predicted values of survival probabilities for the testing set.
#' @param IRange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#'Or the scale you want to get the IAE and ISE;  .
#'
#' @return Estimates of the IAE and ISE
#' @author Zhou HanPu \email{zhouhanpu@csu.edu.cn}
#' @references
#' Marron, J. S. , &  Wand, M. P. . (1992). Exact mean integrated squared error. Annals of Statistics, 20(2), 712-736.
#' HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). \\(l_1\\) splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#' Kowalczuk, & Z. (1998). Integrated squared error and integrated absolute error in recursive identification of continuous-time plants. Control 98 Ukacc International Conference on (Vol.1998, pp.693-698). IET.
#' @examples
#'
#' library(randomForestSRC)
#' library(MASS)
#' library(survival)
#' set.seed(123)
#' N=100
#' p = 25
#' c.mean = 0.4
#' mu = rep(0,p)
#' Si = matrix(0,p,p)
#' for (i in 1:p) {
#'   for (j in 1:p) {
#'     Si[i,j] = 0.9^abs(i-j)
#'   }
#' }
#' W = mvrnorm(N,mu,Si)
#' Ti = rep(0,N)
#' for (i in 1:N) {
#'   t.mu = exp(0.1*sum(W[i,1:round(p/2)]))
#'   Ti[i] = rexp(1,1/t.mu)
#' }
#' c.time = rexp(N,c.mean)
#' mydata.x = as.data.frame(W)
#' mydata.time = data.frame('time'=Ti,'c.time' = c.time,'status' = 0)
#' mydata.time$status = ifelse(mydata.time$time < mydata.time$c.time,1,0)
#' mydata.time$time = apply(mydata.time[,1:2],1,min)
#' mydata = data.frame('time' = mydata.time$time,'status' = mydata.time$status,W)
#'
#' index.train = sample(1:N,2/3*N)
#' data.train = mydata[index.train,]
#' data.test = mydata[-index.train,]
#' fit.RSF = rfsrc(Surv(time,status)~.,data.train,nsplit=3,ntree=500)
#' predicted = predict(fit.RSF,data.test)
#' sp_matrix = predicted$survival
#' object = Surv(data.test$time,data.test$status)
#' #a vector for all the distinct time
#' IAEISE(object,sp_matrix,with(data.train,time[status==1]))
#' #a range
#' IAEISE(object,sp_matrix,c(12,350))
#'
#' @importFrom survival Surv
#' @importFrom survival survfit
#'
#' @export
IAEISE <- function(object,sp_matrix,IRange = range(object[,1])){

  if(!inherits(object, "Surv")){
    stop("object is not of class Surv")
  }

  if(missing(object)){
    stop("The survival object of the testing set is missing")
  }

  if(missing(sp_matrix)) {
    stop("The prediction of the survival probability matrix is missing")
  }

  if (TRUE %in% (is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (TRUE %in% (is.na(sp_matrix))) {
    stop("The input probability matrix cannot have NA")
  }

  if (TRUE %in% is.na(IRange)) {
    stop("Cannot calculate IAE or ISE in the interval containing NA")
  }

  if(!is.numeric(IRange)){
    stop("The class of the IRange must be numeric! or use the default setting")
  }

  if (length(IRange) <= 1) {
    stop("Can not calculate the integration at a single point")
  }

  # Simulation to generate discrete time series
  sp_matrix2sp = apply(sp_matrix, 2, mean)
  if(length(IRange)!=ncol(sp_matrix)){
    t_IRange = range(IRange)
    p = ncol(sp_matrix)
    IRange = seq(t_IRange[1],t_IRange[2],length = p)
  }


  fit0 = survfit(object~1)
  res_sum0 = survminer::surv_summary(fit0)
  KM_frame0 = data.frame('time' = res_sum0$time,
                         'sp0' = res_sum0$surv,'sp1' = 0)[res_sum0$n.event!=0,]


  KM_frame1 = data.frame('time' = IRange,
                         'sp0' = sp_matrix2sp,'sp1' = 0)

  for (i in 1:nrow(KM_frame0)) {
    t = KM_frame0$time[i]
    if (t < min(KM_frame1$time) ) {
      KM_frame0$sp1[i] = 1
    }else{
      index = max(which(KM_frame1$time <= t))
      KM_frame0$sp1[i] = KM_frame1$sp0[index]
    }
  }

  KM_frame0$t_IAE = abs(KM_frame0$sp0 - KM_frame0$sp1)
  KM_frame0$t_ISE = (KM_frame0$sp0 - KM_frame0$sp1)^2
  IAE = 0
  ISE = 0
  for (i in 1:(nrow(KM_frame0)-1)) {
    IAE = IAE + (KM_frame0$time[i+1] - KM_frame0$time[i])*KM_frame0$t_IAE[i]
    ISE = ISE + (KM_frame0$time[i+1] - KM_frame0$time[i])*KM_frame0$t_ISE[i]
  }
  IAEISE_value = c(IAE,ISE)
  names(IAEISE_value) = c('IAE','ISE')
  return(IAEISE_value)
}

