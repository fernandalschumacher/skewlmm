
<!-- README.md is generated from README.Rmd. Please edit that file -->

# skewlmm

<!-- badges: start -->

<!-- badges: end -->

The goal of skewlmm is to fit skew robust linear mixed models, using
scale mixture of skew-normal linear mixed models with possible
within-subject dependence structure, using an EM-type algorithm.

## Installation

<!-- You can install the released version of lmmsmsn from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("lmmsmsn") -->

<!-- ``` -->

You can install skewlmm from GitHub with:

``` r
library(devtools)
install_github("fernandalschumacher/skewlmm")
```

## Example

This is a basic example which shows you how to fit a SMSN-LMM:

``` r
library(skewlmm)
    dat1 <- as.data.frame(nlme::Orthodont)
    fm1 <- smsn.lmm(dat1,formFixed=distance ~ age,groupVar="Subject",quiet=T)
    fm1
#> Linear mixed models with distribution sn and dependency structure CI 
#> Call:
#> smsn.lmm(data = dat1, formFixed = distance ~ age, groupVar = "Subject", 
#>     quiet = T)
#> 
#> Fixed:distance ~ age
#> Random:~1
#> <environment: 0x0000000017927488>
#>   Estimated variance (D):
#>             (Intercept)
#> (Intercept)    6.599775
#> 
#> Estimated parameters:
#>      (Intercept)    age sigma2 Dsqrt1 lambda1
#>          16.7630 0.6602 2.0245 2.5690  1.1062
#> s.e.      1.0067 0.0699 0.1914 1.1731  1.7696
#> 
#> Model selection criteria:
#>    logLik     AIC     BIC
#>  -221.658 453.316 466.726
#> 
#> Number of observations: 108 
#> Number of groups: 27
    summary(fm1)
#> Linear mixed models with distribution sn and dependency structure CI 
#> Call:
#> smsn.lmm(data = dat1, formFixed = distance ~ age, groupVar = "Subject", 
#>     quiet = T)
#> 
#> Distribution sn
#> Random effects: ~1
#> <environment: 0x0000000017927488>
#>   Estimated variance (D):
#>             (Intercept)
#> (Intercept)    6.599775
#> 
#> Fixed effects: distance ~ age
#> with approximate confidence intervals
#>                  Value  Std.error IC 95% lower IC 95% upper
#> (Intercept) 16.7629611 1.00673455    14.789798   18.7361245
#> age          0.6601852 0.06987075     0.523241    0.7971293
#> 
#> Dependency structure: CI
#>   Estimate(s):
#>  sigma2 
#> 2.02447 
#> 
#> Skewness parameter estimate: 1.10616
#> 
#> Model selection criteria:
#>    logLik     AIC     BIC
#>  -221.658 453.316 466.726
#> 
#> Number of observations: 108 
#> Number of groups: 27
```

For more examples, see help(smsn.lmm).
