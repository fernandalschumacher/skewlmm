## skewlmm 1.1.0 _(2022-09-27)_

* `rsmsn.clmm` was added to simulate data from SMSN-CLMM with censored responses.
* `smn.clmm` was added to fit linear mixed models with censored responses.
* Generic functions such as `print`, `summary`, `plot`, `residuals`, `fitted`, `predict`, and `update` work for objects of class `SMNclmm`.
* `mahalDist.SMNclmm` computes the squared Mahalanobis distance from a fitted SMN-CLMM. 

## skewlmm 1.0.0 _(2021-09-17)_

* The main functions `smn.lmm()` and `smsn.lmm()` were changed to pass control options through the `lmmControl` function.
* The default estimation method was changed to the DAAREM method, to improve the general performance.
* A parallel optimization was added using optimParallel.
* Some options of the diagnostic tools were adjusted.

## skewlmm 0.2.3 _(2021-02-03)_

* A bug in the `acfresid()` plot was fixed.
* `lr.test()` now allows testing DEC versus AR1 or CAR1.
* `smn.lmm()` and `smnn.lmm()` were adjusted to work when `data` is a tibble.
* The abbreviation for the conditionally uncorrelated model was changed to "UNC" to match the notation in the theoretical paper, but "CI" is still accepted.

## skewlmm 0.2.2 _(2020-07-08)_

* The package rlist was removed from imports.

## skewlmm 0.2.1 _(2020-07-06)_

* The default initial values were changed.
* Some bugs were fixed (for R version 4.0.1).

## skewlmm 0.2.0 _(2020-05-15)_

* Some functions for model diagnostic were added: `acfresid()`, `mahalDist()`, `healy.plot()`.
* The argument `showCriterium` was added to functions `smn.lmm()` and `smnn.lmm()`to allow the criterium to be shown at each iteration.
* The argument `calc.bi` was removed from functions `smn.lmm()` and `smnn.lmm()`.

## skewlmm 0.1.0 _(2020-03-02)_

* Initial release.
