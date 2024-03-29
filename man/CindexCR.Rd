% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CindexCR.R
\name{CindexCR}
\alias{CindexCR}
\title{Concordance index in the Presence of Competing Risks}
\usage{
CindexCR(time, status, predicted, Cause_int = 1)
}
\arguments{
\item{time}{minimum value of deletion time and survival time.}

\item{status}{the status indicator, for models with competing risks, the status indicator is 0=censored, 1=event at \code{time}, 2= competing risks at \code{time}.}

\item{predicted}{a vector of predicted values or the survival time of survival probabilities of each observation.}

\item{Cause_int}{event type of interest, the default value is 1.}
}
\value{
Estimates of the C-index in the presence of competing risks.
}
\description{
The C-index (Concordance index) of the prognostic model in the presence of competing risks according to Marcel, W et al.(2014).
}
\examples{
time <- c(4, 7, 5, 8)
status <- rep(1, 4)
predicted <- c(3, 5, 7, 10)
Cause_int <- 1
CindexCR(time, status, predicted, Cause_int)

}
\references{
Marcel, W. ,  Paul, B. ,  Koller, M. T. ,  Witteman, J. , &  Gerds, T. A. . (2014).Concordance for prognostic models with competing risks. Biostatistics(3), 526.

Ishwaran, H. ,  Kogalur, U. B. ,  Blackstone, E. H. , &  Lauer, M. S. . (2008). Random survival forests. Journal of Thoracic Oncology Official Publication of the International Association for the Study of Lung Cancer, 2(12), 841-860.
}
\author{
HanPu Zhou \email{zhouhanpu@csu.edu.cn}
}
