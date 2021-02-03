
## skewlmm 0.2.3 _(2021-02-03)_

* A bug in the `acfresid()` plot was fixed.
* `lr.test()` now allows testing DEC versus AR1 or CAR1.
* `smn.lmm()` and `smnn.lmm()` were ajusted to work when `data` is a tibble.
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
