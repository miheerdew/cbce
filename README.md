# CBCE

This repository provides an R package for Multi-view data analysis. Consider two types of high-dimensional measurements on the same samples. CBCE (Correlation Bi-Community Extraction method) finds a set of features `A`, from the first measurement type, and set of features `B`, from the second measurement type, such that features in `A` and `B` are correlated to each other in aggregate. 

Formally the pair `(A,B)` is called a bimodule and the algorithm called the Bimodule Search Procedure (BSP) is introduced in [[1]](#1). We have used this method for the analysis of multi-view data in areas like genomics and climate science. 


## Features of CBCE

-  [RCpp](http://www.rcpp.org/) implementation of the iterative testing framework; multicore if using [ROpen](https://mran.microsoft.com/rro).  
- Multiple backends to calculate p-values. It is also easy to use your own backend. 
- Code tested using [testthat](https://github.com/r-lib/testthat/).
- A simple GUI interface to monitor progress and terminate early. 
- Documented using [Roxygen](https://roxygen2.r-lib.org/) and [pkgdown](https://pkgdown.r-lib.org/).

## How to install CBCE

You can install the latest version of cbce directly from the github repo by first installing [devtools](https://github.com/r-lib/devtools).

``` r
if("devtools" %in% rownames(installed.packages()) == FALSE) {
  install.packages("devtools")
}
devtools::install_github("miheerdew/cbce")
```

### MacOS compilation issues 

(July 2023, update)

If you are facing compilation issues on MacOS due to the inability to load the `gfortran` library, you may need to [check your setup](https://mac.r-project.org/tools/). I was able to address this issue on my computer by using the `macrtools::gfortran_install()` function from the [macrtools package](https://github.com/coatless-mac/macrtools).

## Example usage

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

# Recovers the indices 4:5 for X and 3:6 for Y
# If the strength of the correlation was higher
# all the indices could be recovered.
res$comms
```
## Documentation
More information is available on the [software webpage](https://miheerdew.github.io/cbce/reference/index.html).

## Acknowledgement

This project has been funded by NIH R01 HG009125-01 grant.


## References
<a id="1">[1]</a> 
Dewaskar, Miheer, John Palowitch, Mark He, Michael I. Love, and Andrew Nobel. "Finding Stable Groups of Cross-Correlated Features in Multi-View data." arXiv preprint arXiv:2009.05079 (2020).
