# cbce

Consider two sets of high-dimensional measurements on the same set of samples. CBCE (Correlation Bi-Community Extraction method) finds sets of variables from the first measurement which and sets of variables from the second measurement which are correlated to each other.

## Installation

You can install the latest version of cbce directly from the github repo:

``` r
if("devtools" %in% rownames(installed.packages()) == FALSE) {
  install.packages("devtools")
}
devtools::install_github("miheerdew/cbce")
```

## Example

``` r
library(cbce)

#Sample size
n <- 40
#Dimension of measurement 1
dx <- 20
#Dimension of measurement 2
dy <- 50

#Correlation strength
rho <- 0.5

set.seed(1245)

# Assume first measurement is gaussian
X <- matrix(rnorm(dx*n), nrow=n, ncol=dx)

# Measurements 3:6 in set 2 are correlated to 4:7 in set 1
Y <- matrix(rnorm(dy*n), nrow=n, ncol=dy)
Y[, 3:6] <- sqrt(1-rho)*Y[, 3:6] + sqrt(rho)*rowSums(X[, 4:5])

res <- cbce(X, Y)

#Recovers the indices 4:5 for X and 3:6 for Y
res$comms[[1]]
```

