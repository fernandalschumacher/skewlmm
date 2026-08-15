fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar = "Subject")

test_that("confint default (asymptotic) returns intervals for all parameters", {
  ci <- confint(fm1)

  expect_equal(nrow(ci), length(fm1$theta))
  expect_equal(colnames(ci), c("Estimate", "Std Error", "CI 95% lower", "CI 95% upper"))
})

test_that("confint with parm = 'beta' returns intervals for fixed effects only", {
  ci <- confint(fm1, parm = "beta")

  expect_equal(nrow(ci), length(fm1$estimates$beta))
})

test_that("confint defaults parm to 'beta' when method = 'sandwich'", {
  ci <- confint(fm1, method = "sandwich", parallel = FALSE, MCiter = 5, seed = 1)

  expect_equal(nrow(ci), length(fm1$estimates$beta))
})

test_that("confint warns and resets parm to 'beta' when method = 'sandwich' and parm = 'all'", {
  expect_warning(
    ci <- confint(fm1, parm = "all", method = "sandwich", parallel = FALSE, MCiter = 5, seed = 1),
    "parm has been changed back to 'beta'"
  )
  expect_equal(nrow(ci), length(fm1$estimates$beta))
})

test_that("confint errors for an invalid level", {
  expect_error(confint(fm1, level = 1.5), "level must be a number between 0 and 1")
})
