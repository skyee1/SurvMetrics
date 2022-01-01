#' IBS
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
#' @param object object of class \code{Surv} in the testing set created by Surv function.
#' @param sp_matrix a matrix or data.frame of predicted values of survival probabilities for the testing set.
#' @param IBSrange a vector contains all discrete time points corresponding to the predicted probability in sp_matrix.
#'Or the scale you want to get the IBS; and if it is a single point the return value will be the Brier Score at the timepoint.
#'
#' @return The integration of brierscore
#'
#' @author Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
#' @references
#'HooraMoradian, DenisLarocque, & FranoisBellavance. (2017). \\(l_1\\) splitting rules in survival forests. Lifetime Data Analysis, 23(4), 671–691.
#'
#'Graf, Erika, Schmoor, Claudia, Sauerbrei, & Willi, et al. (1999). Assessment and comparison of prognostic classification schemes for survival data. Statist. Med., 18(1718), 2529-2545.
#'
#'Brier, G. W. . (1950). Verification of forecasts expressed in terms of probability. Monthly Weather Review, 78.
#'
#'Gneiting, T. , &  Raftery, A. E. . (2007). Strictly Proper Scoring Rules, Prediction, and Estimation.

#' @examples
#' library(randomForestSRC)
#' library(survival)
#' library(SurvMetrics)
#' set.seed(123)
#' N = 100
#' mydata = SDGM4(N, p = 20, c_step = -0.5)
#' index.train = sample(1:N,2/3*N)
#' data.train = mydata[index.train,]
#' data.test = mydata[-index.train,]
#'
#' fit.RSF = rfsrc(Surv(time,status)~.,data.train,nsplit=3,ntree=500)
#' predicted = predict(fit.RSF,data.test)
#' sp_matrix = predicted$survival
#'
#' object = Surv(data.test$time,data.test$status)
#'
#' #the default time points
#' IBS(object, sp_matrix, predicted$time.interest)
#' #a time range
#' IBS(object,sp_matrix,c(18:100))
#'
#'
#' @importFrom survival Surv
#'
#'
#' @export
IBS <- function(object,sp_matrix,IBSrange = range(object[,1])){

  if(!inherits(object, "Surv"))
    stop("object is not of class Surv")

  if (missing(object)) {
    stop("The survival object of the testing set is missing")
  }

  if(missing(sp_matrix)){
    stop("The prediction of the survival probability matrix is missing")
  }

  if(!is.numeric(IBSrange)){
    stop("The class of the IBSrange must be numeric! or use the default setting")
  }

  if (TRUE %in% (is.na(object))) {
    stop("The input vector cannot have NA")
  }

  if (TRUE %in% (is.na(sp_matrix))) {
    stop("The input probability matrix cannot have NA")
  }

  if (TRUE %in% is.na(IBSrange)) {
    stop("Cannot calculate IBS in the interval containing NA")
  }

  if(FALSE %in% (IBSrange > 0)){
    stop("The integration interval must be positive")
  }

  if(FALSE %in% (diff(IBSrange) > 0)){
    stop("The integral interval value must increase")
  }

  if(length(object)!=nrow(sp_matrix)){
    stop("The instances be predicted (number of rows of the sp_matrix)and the survival object have different lengths")
    }

  #Determine the input dimension
  if (ncol(data.frame('sp' = sp_matrix))==1 & length(IBSrange)==1) {
    bs = Brier(object,sp_matrix,IBSrange)
    names(bs) = 'Brier Score'
    return(bs)

  }else if((ncol(data.frame('sp' = sp_matrix))-1)*(length(IBSrange)-1)==0){
    stop("The length is illegal")

  }else if(ncol(sp_matrix) == length(IBSrange)){
    IBSrange = sort(IBSrange)
    t_brier = rep(0)
    for (i in 1:length(IBSrange)) {
      pre_sp = sp_matrix[,i]
      t_star = IBSrange[i]
      t_brier[i] = Brier(object,pre_sp,t_star)
    }
    t_IBS = 0
    for (i in 1:(length(IBSrange)-1)) {
      t_IBS = t_IBS + (IBSrange[i+1] - IBSrange[i])*t_brier[i]
    }
    t_IBS = t_IBS/(range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) = 'IBS'
    return(t_IBS)
  }else {
    t_brier = rep(0)
    t_IBSrange = range(IBSrange)
    p = ncol(sp_matrix)
    #Simulation to generate discrete time series
    IBSrange = seq(t_IBSrange[1],t_IBSrange[2],length = p)

    for (i in 1:length(IBSrange)) {
      pre_sp = sp_matrix[,i]
      t_star = IBSrange[i]
      t_brier[i] = Brier(object,pre_sp,t_star)
    }
    t_IBS = 0
    for (i in 1:(length(IBSrange)-1)) {
      t_IBS = t_IBS + (IBSrange[i+1] - IBSrange[i])*t_brier[i]
    }
    t_IBS = t_IBS/(range(IBSrange)[2] - range(IBSrange)[1])
    names(t_IBS) = 'IBS'
    return(t_IBS)
  }
}


