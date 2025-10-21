test_that("myncurve returns correct mu", {
  result <- MATH4753F25fabian::myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$mu, 10)
})

test_that("myncurve returns correct sigma", {
  result <- MATH4753F25fabian::myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$sigma, 5)
})

test_that("myncurve returns correct probability", {
  result <- MATH4753F25fabian::myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$prob ,pnorm(6,10,5) )
})
