% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDGM2.R
\name{SDGM2}
\alias{SDGM2}
\title{Survival Data Generation Method 2}
\usage{
SDGM2(N = 200, p = 15, u_max = 4)
}
\arguments{
\item{N}{The sample size of the simulated dataset.}

\item{p}{The covariate dimension of the simulated dataset.}

\item{u_max}{The parameter which is used to control the censoring rate.}
}
\value{
the simulated dataset
}
\description{
Survival data generation method. The dataset represents mild violations of the proportional hazards assumption.
}
\examples{
SDGM2(N = 200, p = 15, u_max = 4)

}
\references{
Steingrimsson, J. A. ,  Diao, L. , &  Strawderman, R. L. . (2019). Censoring unbiased regression trees and ensembles. Journal of the American Statistical Association, 114.

Zhu, R. , &  Kosorok, M. R. . (2012). Recursively imputed survival trees. Journal of the American Statistical Association, 107(497), 331-340.

Ishwaran, H. ,  Kogalur, U. B. ,  Gorodeski, E. Z. ,  Minn, A. J. , &  Lauer, M. S. . (2010). High-dimensional variable selection for survival data. Journal of the American Statistical Association, 105(489), 205-217.
}
\author{
Hanpu Zhou \email{zhouhanpu@csu.edu.cn}
}
