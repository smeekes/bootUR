context("Test missing value functions")

test_that("matrix without internal missings passes checks", {
  simple_matrix <- matrix(1, nrow = 10, ncol = 5)
  t1 <- check_missing_insample_values(simple_matrix)
  expect_equal(t1, rep(FALSE, 5))

  simple_matrix[1, 3] <- NA
  simple_matrix[10, 5] <- NA
  t1 <- check_missing_insample_values(simple_matrix)
  expect_equal(t1, rep(FALSE, 5))
})

test_that("matrix with internal missings fails checks", {
  simple_matrix <- matrix(1, nrow = 10, ncol = 5)
  simple_matrix[4, 2] <- NA
  t1 <- check_missing_insample_values(simple_matrix)
  expect_equal(t1[-2], rep(FALSE, 4))
  expect_equal(t1[2], TRUE)
})

test_that("find_missing_values gives correct output when no missing", {
  simple_matrix <- matrix(1, nrow = 10, ncol = 5)
  t1 <- find_nonmissing_subsample(simple_matrix)
  expect_equal(t1$all_equal, TRUE)
  for (i in 1:5) {
    expect_equivalent(t1$range[1, i], 1)
    expect_equivalent(t1$range[2, i], 10)
  }
})

test_that("find_missing_values gives correct output when missing", {
  simple_matrix <- matrix(1, nrow = 10, ncol = 5)
  simple_matrix[1, 3] <- NA
  simple_matrix[10, 3] <- NA
  t1 <- find_nonmissing_subsample(simple_matrix)
  expect_equal(t1$all_equal, FALSE)
  expect_equivalent(t1$range[1, 3], 2)
  expect_equivalent(t1$range[2, 3], 9)
  for (i in c(1, 2, 4, 5)) {
    expect_equivalent(t1$range[1, i], 1)
    expect_equivalent(t1$range[2, i], 10)
  }
})

