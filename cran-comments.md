## Submission summary

This is an update from skewlmm 1.1.3 (current CRAN version) to 1.1.5.

Main changes:

* Added variance-type estimators for standard errors in `confint()`
  (sandwich, bootstrap, and asymptotic intervals).
* Removed the placeholder censored-response functions (`smn.clmm()`,
  `rsmsn.clmm()`, and related methods). In 1.1.3 these were already
  unavailable and only returned an error, because they depended on
  archived packages. Their current dependency, 'relliptical', was also
  archived on CRAN on 2026-07-10. They will be restored once a suitable
  dependency is available.
* Added unit tests with testthat.

## Test environments

* local Windows 11 x64, R 4.6.0
* [win-builder (devel and release)]
* [R-hub / GitHub Actions, if used]

## R CMD check results

0 errors | 0 warnings | [n] notes

[Paste and briefly explain any NOTEs here.]

## Reverse dependencies

[There are currently no reverse dependencies. / We checked N reverse
dependencies and found no new problems.]
