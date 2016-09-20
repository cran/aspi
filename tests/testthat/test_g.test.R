test_that("g.test throws expected error when fed data containing zero counts",{
  expect_error(g.test(diplostomum_lenses),"Data contains zero counts!")
})

result <- g.test(diplostomum_eyes_excl_lenses)

test_that("g.test summary reports results of three tests",{
  expect_equal(length(result$summary$Test), 3)
  expect_equal(length(names(result$summary)), 4)
  expect_match(as.character(result$summary$Test[1]), "Pooled")
  expect_match(as.character(result$summary$Test[2]), "Heterogeneity")
  expect_match(as.character(result$summary$Test[3]), "Total")
})

test_that("Total df equals sum of Pooled and Heterogeneity df",{
  expect_equal(result$summary$df[1] + result$summary$df[2], result$summary$df[3])
})

test_that("Total G equals sum of Pooled and Hetereogeneity G",{
  expect_equal(result$summary$G[1] + result$summary$G[2], result$summary$G[3])
})

test_that("Results are correct for host 3", {
  expect_match(as.character(result$hosts$Host[3]), "3")
  expect_equal(result$hosts$Left[3], 195)
  expect_equal(result$hosts$Right[3], 133)
  expect_equal(result$hosts$G[3], 11.790, tolerance=0.001)
  expect_equal(result$hosts$p[3], 0.0006, tolerance=0.0001)
  expect_equal(result$hosts$BH[3], 0.0033, tolerance=0.0001)
  expect_equal(result$hosts$Holm[3], 0.0250, tolerance=0.0001)
})