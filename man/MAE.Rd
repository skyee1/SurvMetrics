% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAE.R
\name{MAE}
\alias{MAE}
\title{The Mean Absolute Error}
\usage{
MAE(object, pre_time)
}
\arguments{
\item{object}{object of class \code{Surv} created by Surv function.}

\item{pre_time}{a vector of predicted values of survival time of each observation.}
}
\value{
the value of the Mean Absolute Error between the survival time and the predicted result.
}
\description{
A somewhat naive criterion that is sometimes used consists of simply omitting all censored cases from the data set.
For survival analysis problems, the mean absolute error (MAE) can be defined as an average of the differences between the predicted time values and the actual observation time values.
Only the samples for which the event occurs are being considered in this metric.
}
\details{
Condition: MAE can only be used for the evaluation of survival models which can provide the event time as the predicted target value.
}
\examples{
library(survival)
time <- rexp(50)
status <- sample(c(0, 1), 50, replace = TRUE)
pre_time <- rexp(50)
MAE(Surv(time, status), pre_time)

}
\references{
Matsuo, K. ,  Purushotham, S. ,  Jiang, B. ,  Mandelbaum, R. S. ,  Takiuchi, T. , &  Liu, Y. , et al. (2018). Survival outcome prediction in cervical cancer: cox models vs deep-learning model. American Journal of Obstetrics & Gynecology.
Coyle, E. J. , &  Lin, J. H. . (1988). Stack filters and the mean absolute error criterion. IEEE Trans Acoustics Speech Signal Processing, 36(8), 1244-1254.
}
\author{
Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
}
