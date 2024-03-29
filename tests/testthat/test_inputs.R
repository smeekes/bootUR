context("Test inputs of boot functions")

test_that("univ boot refuse matrix", {
  y <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  expect_error(boot_adf(y, B = 19, do_parallel = FALSE, show_progress = FALSE))
  expect_error(boot_union(y, B = 19, do_parallel = FALSE, show_progress = FALSE))
})

test_that("Missing values in sample", {
  y <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  y[5, 1] <- NA
  expect_error(boot_adf(y[, 1], B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
  expect_error(boot_union(y[, 1], union_quantile = 0.1, B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
  expect_error(boot_ur(y, level = 0.1, B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
  expect_error(boot_fdr(y, FDR_level = 0.1, B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
  expect_error(boot_sqt(y, SQT_level = 0.1, B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
  expect_error(boot_panel(y, B = 19, do_parallel = FALSE, show_progress = FALSE), "Missing values detected inside sample.")
})

test_that("Unbalanced Panels", {
  y <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
  y[1, 1] <- NA
  expect_warning(boot_ur(y, level = 0.1, B = 19, bootstrap = "MBB", do_parallel = FALSE, show_progress = FALSE), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_warning(boot_ur(y, level = 0.1, B = 19, bootstrap = "MBB", do_parallel = FALSE, show_progress = FALSE), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_error(boot_fdr(y, FDR_level = 0.1, B = 19, bootstrap = "MBB", do_parallel = FALSE, show_progress = FALSE), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(boot_sqt(y, SQT_level = 0.1, B = 19, bootstrap = "MBB", do_parallel = FALSE, show_progress = FALSE), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(boot_panel(y, B = 19, bootstrap = "MBB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_warning(boot_ur(y, level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_warning(boot_ur(y, level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_error(boot_fdr(y, FDR_level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(boot_sqt(y, SQT_level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(boot_panel(y, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
})

test_that("Panel Sieve Bootstrap", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_warning(boot_fdr(y, FDR_level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
  expect_warning(boot_sqt(y, SQT_level = 0.1, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
  expect_warning(boot_panel(y, B = 19, bootstrap = "SB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
  expect_warning(boot_fdr(y, FDR_level = 0.1, B = 19, bootstrap = "SWB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
  expect_warning(boot_sqt(y, SQT_level = 0.1, B = 19, bootstrap = "SWB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
  expect_warning(boot_panel(y, B = 19, bootstrap = "SWB", do_parallel = FALSE, show_progress = FALSE), "SB and SWB bootstrap only recommended for boot_ur; see help for details.")
})

test_that("Union specs", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_warning(boot_ur(y, level = 0.1, B = 19, deterministics = "intercept", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Deterministic specification in argument deterministics is ignored, as union test is applied.")
  expect_warning(boot_fdr(y, FDR_level = 0.1, B = 19, deterministics = "trend", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Deterministic specification in argument deterministics is ignored, as union test is applied.")
  expect_warning(boot_sqt(y, SQT_level = 0.1, B = 19, deterministics = "intercept", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Deterministic specification in argument deterministics is ignored, as union test is applied.")
  expect_warning(boot_panel(y, B = 19, deterministics = "trend", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Deterministic specification in argument deterministics is ignored, as union test is applied.")
  expect_warning(boot_ur(y, level = 0.1, B = 19, detrend = "OLS", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Detrending method in argument detrend is ignored, as union test is applied.")
  expect_warning(boot_fdr(y, FDR_level = 0.1, B = 19, detrend = "QD", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Detrending method in argument detrend is ignored, as union test is applied.")
  expect_warning(boot_sqt(y, SQT_level = 0.1, B = 19, detrend = "OLS", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Detrending method in argument detrend is ignored, as union test is applied.")
  expect_warning(boot_panel(y, B = 19, detrend = "QD", bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Detrending method in argument detrend is ignored, as union test is applied.")
})

test_that("BSQT q spec", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_error(boot_sqt(y, steps = c(0, 0.5, 1, 2), SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Invalid input values for steps: must be quantiles or positive integers.")
  expect_error(boot_sqt(y, steps = c(-1, 0, 1, 2), SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Invalid input values for steps: must be quantiles or positive integers.")
  expect_error(boot_sqt(y, steps = c(0, 1, NA, 3), SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Invalid input values for steps: must be quantiles or positive integers.")
  expect_warning(boot_sqt(y, steps = 1:5, SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Input to argument steps transformed to fit sequential test:")
  expect_warning(boot_sqt(y, steps = 0:4, SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Input to argument steps transformed to fit sequential test:")
  expect_warning(boot_sqt(y, steps = 1:5/5, SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Input to argument steps transformed to fit sequential test:")
  expect_warning(boot_sqt(y, steps = 0:4/5, SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Input to argument steps transformed to fit sequential test:")
  expect_warning(boot_sqt(y, steps = 0:9/10, SQT_level = 0.1, B = 19, bootstrap = "AWB", do_parallel = FALSE, show_progress = FALSE), "Input to argument steps transformed to remove duplicate groups after transformation to integers:")
})
