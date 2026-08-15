test_that("smn.lmm fits a basic model and returns expected structure", {
  fm1 <- smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar = "Subject")

  expect_s3_class(fm1, "SMN")
  expect_length(fm1$estimates$beta, 3)
  expect_false(is.null(fm1$std.error))
  expect_true(all(fm1$std.error > 0))
})

test_that("smn.lmm validates its inputs", {
  expect_error(
    smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar = "Subject",
            formRandom = "not a formula"),
    "formRandom must be a formula"
  )
  expect_error(
    smn.lmm(distance ~ age + Sex, data = as.matrix(nlme::Orthodont), groupVar = "Subject"),
    "data must be a data.frame"
  )
  expect_error(
    smn.lmm(distance ~ age + Sex, data = nlme::Orthodont, groupVar = "not_a_column"),
    "not found in data"
  )
})
