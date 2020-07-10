context("Test replicability of RNG sequences")

test_that("Replicability DF", {
  set.seed(1111)
  out1a <- boot_df(MacroTS[, 1], B = 199, boot = "AWB", dc = 2)
  out2a <- boot_df(MacroTS[, 6], B = 199, boot = "MBB", dc = 2)
  set.seed(1111)
  out1b <- boot_df(MacroTS[, 1], B = 199, boot = "AWB", dc = 2)
  out2b <- boot_df(MacroTS[, 6], B = 199, boot = "MBB", dc = 2)
  expect_identical(out1a, out1b)
  expect_identical(out2a, out2b)
})

test_that("Replicability union", {
  set.seed(1111)
  out1a <- boot_union(MacroTS[, 12], B = 199, boot = "SWB")
  out2a <- boot_union(MacroTS[, 16], B = 199, boot = "BWB")
  set.seed(1111)
  out1b <- boot_union(MacroTS[, 12], B = 199, boot = "SWB")
  out2b <- boot_union(MacroTS[, 16], B = 199, boot = "BWB")
  expect_identical(out1a, out1b)
  expect_identical(out2a, out2b)
})
