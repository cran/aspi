result <- eb.test(diplostomum_lenses)
host28 <- result$hosts[3,]

test_that(".prepareData removes uninfected hosts",{
  expect_equal(length(result$hosts[,1]), 31)
})

test_that("binomial exact test on pooled data returns correct value",{
  expect_equal(result$pooled, 0.8439, tolerance=0.0001)
})

test_that("binomial exact test on an individual host returns correct values",{
  expect_match(as.character(host28$Host), "28")
  expect_equal(host28$Left, 0)
  expect_equal(host28$Right, 3)
  expect_equal(host28$p, 0.25, tolerance=0.001)
  expect_equal(host28$BH, 0.9688, tolerance=0.0001)
  expect_equal(host28$Holm, 1, tolerance=0.1)
})