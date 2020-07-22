context("Differencing multiple time series")

test_that("d only contains non-negative integers", {
  y <- matrix(rnorm(10), nrow = 5, ncol = 2)
  d <- c(1, -1)
  expect_error(diff_mult(y, d))
  d <- c(1, 1.1)
  expect_error(diff_mult(y, d))
})

test_that("dimensions work", {
  y <- matrix(rnorm(15), nrow = 5, ncol = 3)
  d <- rep(1, 4)
  expect_error(diff_mult(y, d))
  expect_error(diff_mult(y, d[1:2]))
})

test_that("differencing works", {
  x <- rep(1, 10)
  y <- cbind(x, cumsum(x), cumsum(cumsum(x)))
  colnames(y) <- c("I0", "I1", "I2")
  d <- 0:2
  xx <- matrix(x, nrow = length(x), ncol = 3)
  xx[1, 2] <- NA
  xx[1:2, 3] <- NA
  colnames(xx) <- colnames(y)
  expect_identical(diff_mult(y, d), xx)
  expect_identical(diff_mult(y, d, keep_NAs = FALSE), xx[-(1:2), ])
})

test_that("class is kept", {
  x <- rep(1, 10)
  y <- cbind(x, cumsum(x), cumsum(cumsum(x)))
  colnames(y) <- c("I0", "I1", "I2")
  y <- ts(y, start = 2011)
  y2 <- data.frame(y)
  rownames(y2) <- 2011:2020
  colnames(y2) <- colnames(y)
  d <- 0:2
  xx <- matrix(x, nrow = length(x), ncol = 3)
  xx[1, 2] <- NA
  xx[1:2, 3] <- NA
  colnames(xx) <- colnames(y)
  xx <- ts(xx, start = 2011)
  xx2 <- data.frame(xx)
  rownames(xx2) <- 2011:2020
  colnames(xx2) <- colnames(y)
  expect_identical(diff_mult(y, d), xx)
  expect_identical(diff_mult(y, d, keep_NAs = FALSE), xx[-(1:2), ])
  expect_identical(diff_mult(y2, d), xx2)
  expect_identical(diff_mult(y2, d, keep_NAs = FALSE), xx2[-(1:2), ])
})
