context("fast marginal lm solver (mfastLmCpp)")


test_that("estimation", {
  expect_equal(mfastLmCpp(1:10, as.matrix(1:10))[1,1], 1)
})
