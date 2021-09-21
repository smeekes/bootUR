context("Test replicability of RNG sequences")

test_that("Replicability ADF", {
  set.seed(1111)
  out1a <- boot_adf(MacroTS[, 1], B = 199, bootstrap = "AWB", dc = 2)
  out2a <- boot_adf(MacroTS[, 6], B = 199, bootstrap = "MBB", dc = 2)
  set.seed(1111)
  out1b <- boot_adf(MacroTS[, 1], B = 199, bootstrap = "AWB", dc = 2)
  out2b <- boot_adf(MacroTS[, 6], B = 199, bootstrap = "MBB", dc = 2)
  expect_identical(out1a, out1b)
  expect_identical(out2a, out2b)
})

test_that("Replicability union", {
  set.seed(1111)
  out1a <- boot_union(MacroTS[, 12], B = 199, bootstrap = "SWB")
  out2a <- boot_union(MacroTS[, 16], B = 199, bootstrap = "BWB")
  set.seed(1111)
  out1b <- boot_union(MacroTS[, 12], B = 199, bootstrap = "SWB")
  out2b <- boot_union(MacroTS[, 16], B = 199, bootstrap = "BWB")
  expect_identical(out1a, out1b)
  expect_identical(out2a, out2b)
})
