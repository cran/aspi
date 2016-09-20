test_that("plotHistogram throws expected error when fed data containing zero counts",{
  expect_error(plotHistogram(diplostomum_lenses),"Data contains zero counts!")
})

test_that("plotVolcano throws expected error when fed data containing zero counts",{
  expect_error(plotVolcano(diplostomum_lenses),"Data contains zero counts!")
})
