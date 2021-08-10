
# SurvMetrics

An implementation of popular evaluation metrics that are commonly used in survival prediction 
  including Concordance Index, Brier Score, Integrated Brier Score, 
  Integrated Square Error, Integrated Absolute Error and Mean Absolute Error.
  For a detailed information, see (Ishwaran H, Kogalur UB, Blackstone EH and Lauer MS (2008) [doi:10.1214/08-AOAS169](https://doi.org/10.1214/08-AOAS169)) and
  (Moradian H, Larocque D and Bellavance F (2017) [doi:10.1007/s10985-016-9372-1](https://doi.org/10.1007/s10985-016-9372-1)) for different evaluation metrics.

## Installation

You can install the released version of SurvMetrics from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SurvMetrics")
```
or Install devtools and build the development version by:
``` r
install.packages("devtools", repos = "https://cloud.r-project.org/")
devtools::install_github("skyee1/SurvMetrics")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SurvMetrics)
library(survival)
time = c(1,1,2,2,2,2,2,2)
status = c(0,1,1,0,1,1,0,1)
predicted = c(2,3,3,3,4,2,4,3)
Cindex(Surv(time,status),predicted)
```

