
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lmmsmsn

<!-- badges: start -->

<!-- badges: end -->

The goal of lmmsmsn is to fit scale mixture of skew-normal linear mixed
models with possible within-subject dependence structure, using an
EM-type algorithm.

## Installation

<!-- You can install the released version of lmmsmsn from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("lmmsmsn") -->

<!-- ``` -->

You can install lmmsmsn from GitHub with:

``` r
library(devtools)
install_github("fernandalschumacher/lmmsmsn")
```

## Example

This is a basic example which shows you how to fit a SMSN-LMM:

``` r
library(lmmsmsn)
    dat1 <- as.data.frame(nlme::Orthodont)
    fm1 <- smsn.lmm(dat1,formFixed=distance ~ age,groupVar="Subject")
#> Iteration  1   of  200 
Iteration  2   of  200 
Iteration  3   of  200 
Iteration  4   of  200 
Iteration  5   of  200 
Iteration  6   of  200 
Iteration  7   of  200 
Iteration  8   of  200 
Iteration  9   of  200 
Iteration  10   of  200 
Iteration  11   of  200 
Iteration  12   of  200 
Iteration  13   of  200 
Iteration  14   of  200 
Iteration  15   of  200 
Iteration  16   of  200 
Iteration  17   of  200 
Iteration  18   of  200 
Iteration  19   of  200 
Iteration  20   of  200 
Iteration  21   of  200 
Iteration  22   of  200 
Iteration  23   of  200 
Iteration  24   of  200 
Iteration  25   of  200 
Iteration  26   of  200 
Iteration  27   of  200 
Iteration  28   of  200 
Iteration  29   of  200 
Iteration  30   of  200 
Iteration  31   of  200 
Iteration  32   of  200 
Iteration  33   of  200 
Iteration  34   of  200 
Iteration  35   of  200 
Iteration  36   of  200 
Iteration  37   of  200 
Iteration  38   of  200 
Iteration  39   of  200 
Iteration  40   of  200 
Iteration  41   of  200 
Iteration  42   of  200 
Iteration  43   of  200 
Iteration  44   of  200 
Iteration  45   of  200 
Iteration  46   of  200 
Iteration  47   of  200 
Iteration  48   of  200 
Iteration  49   of  200 
Iteration  50   of  200 
Iteration  51   of  200 
    fm1
#> Linear mixed models with distribution sn and dependency structure CI 
#> Call:
#> smsn.lmm(data = dat1, formFixed = distance ~ age, groupVar = "Subject")
#> 
#> Fixed:distance ~ age
#> Random:~1
#> <environment: 0x0000000012430348>
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
#> smsn.lmm(data = dat1, formFixed = distance ~ age, groupVar = "Subject")
#> 
#> Distribution sn
#> Random effects: ~1
#> <environment: 0x0000000012430348>
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
