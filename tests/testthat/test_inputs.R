context("Test inputs of boot functions")

test_that("univ boot refuse matrix", {
  y <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  expect_error(boot_df(y, level = 0.1, B = 19))
  expect_error(boot_union(y, level = 0.1, B = 19))
})

test_that("Missing values in sample", {
  y <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  y[5, 1] <- NA
  expect_error(boot_df(y[, 1], level = 0.1, B = 19), "Missing values detected inside sample.")
  expect_error(boot_union(y[, 1], level = 0.1, B = 19), "Missing values detected inside sample.")
  expect_error(iADFtest(y, level = 0.1, B = 19), "Missing values detected inside sample.")
  expect_error(bFDRtest(y, level = 0.1, B = 19), "Missing values detected inside sample.")
  expect_error(BSQTtest(y, level = 0.1, B = 19), "Missing values detected inside sample.")
  expect_error(paneltest(y, level = 0.1, B = 19), "Missing values detected inside sample.")
})

test_that("Unbalanced Panels", {
  y <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)
  y[1, 1] <- NA
  expect_warning(iADFtest(y, level = 0.1, B = 19, boot = "MBB"), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_warning(iADFtest(y, level = 0.1, B = 19, boot = "MBB"), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_error(bFDRtest(y, level = 0.1, B = 19, boot = "MBB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(BSQTtest(y, level = 0.1, B = 19, boot = "MBB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(paneltest(y, level = 0.1, B = 19, boot = "MBB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_warning(iADFtest(y, level = 0.1, B = 19, boot = "SB"), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_warning(iADFtest(y, level = 0.1, B = 19, boot = "SB"), "Missing values cause resampling bootstrap to be executed for each time series individually.")
  expect_error(bFDRtest(y, level = 0.1, B = 19, boot = "SB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(BSQTtest(y, level = 0.1, B = 19, boot = "SB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
  expect_error(paneltest(y, level = 0.1, B = 19, boot = "SB"), "Resampling-based bootstraps MBB and SB cannot handle unbalanced series.")
})

test_that("Panel Sieve Bootstrap", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_warning(bFDRtest(y, level = 0.1, B = 19, boot = "SB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  expect_warning(BSQTtest(y, level = 0.1, B = 19, boot = "SB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  expect_warning(paneltest(y, level = 0.1, B = 19, boot = "SB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  expect_warning(bFDRtest(y, level = 0.1, B = 19, boot = "SWB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  expect_warning(BSQTtest(y, level = 0.1, B = 19, boot = "SWB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
  expect_warning(paneltest(y, level = 0.1, B = 19, boot = "SWB"), "SB and SWB bootstrap only recommended for iADFtest; see help for details.")
})

test_that("Union specs", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_warning(iADFtest(y, level = 0.1, B = 19, dc = 1, boot = "AWB"), "Deterministic specification in argument dc is ignored, as union test is applied.")
  expect_warning(bFDRtest(y, level = 0.1, B = 19, dc = 2, boot = "AWB"), "Deterministic specification in argument dc is ignored, as union test is applied.")
  expect_warning(BSQTtest(y, level = 0.1, B = 19, dc = 1, boot = "AWB"), "Deterministic specification in argument dc is ignored, as union test is applied.")
  expect_warning(paneltest(y, level = 0.1, B = 19, dc = 2, boot = "AWB"), "Deterministic specification in argument dc is ignored, as union test is applied.")
  expect_warning(iADFtest(y, level = 0.1, B = 19, detr = "OLS", boot = "AWB"), "Detrending method in argument detr is ignored, as union test is applied.")
  expect_warning(bFDRtest(y, level = 0.1, B = 19, detr = "QD", boot = "AWB"), "Detrending method in argument detr is ignored, as union test is applied.")
  expect_warning(BSQTtest(y, level = 0.1, B = 19, detr = "OLS", boot = "AWB"), "Detrending method in argument detr is ignored, as union test is applied.")
  expect_warning(paneltest(y, level = 0.1, B = 19, detr = "QD", boot = "AWB"), "Detrending method in argument detr is ignored, as union test is applied.")
})

test_that("BSQT q spec", {
  y <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
  expect_error(BSQTtest(y, q = c(0, 0.5, 1, 2), level = 0.1, B = 19, boot = "AWB"), "Invalid input values for q: must be quantiles or positive integers.")
  expect_error(BSQTtest(y, q = c(-1, 0, 1, 2), level = 0.1, B = 19, boot = "AWB"), "Invalid input values for q: must be quantiles or positive integers.")
  expect_error(BSQTtest(y, q = c(0, 1, NA, 3), level = 0.1, B = 19, boot = "AWB"), "Invalid input values for q: must be quantiles or positive integers.")
  expect_warning(BSQTtest(y, q = 1:5, level = 0.1, B = 19, boot = "AWB"), "Input to argument q transformed to fit sequential test:")
  expect_warning(BSQTtest(y, q = 0:4, level = 0.1, B = 19, boot = "AWB"), "Input to argument q transformed to fit sequential test:")
  expect_warning(BSQTtest(y, q = 1:5/5, level = 0.1, B = 19, boot = "AWB"), "Input to argument q transformed to fit sequential test:")
  expect_warning(BSQTtest(y, q = 0:4/5, level = 0.1, B = 19, boot = "AWB"), "Input to argument q transformed to fit sequential test:")
  expect_warning(BSQTtest(y, q = 0:9/10, level = 0.1, B = 19, boot = "AWB"), "Input to argument q transformed to remove duplicate groups after transformation to integers:")
})
