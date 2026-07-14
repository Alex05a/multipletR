test_that("detect_multiplets returns the expected structure", {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")
  res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                           plotPercent = FALSE, plotTotalReads = FALSE)

  # output is a data frame with the added columns
  expect_s3_class(res, "data.frame")
  expect_true(all(c("our_classification", "pct_human", "pct_mouse") %in%
                    names(res)))
})

test_that("detect_multiplets classifies cells as Multiplet or Singlet", {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")
  res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                           plotPercent = FALSE, plotTotalReads = FALSE)

  # only the two expected labels appear
  expect_setequal(unique(res$our_classification), c("Multiplet", "Singlet"))

  # finds a sensible number of multiplets on the PC65 example
  n_mult <- sum(res$our_classification == "Multiplet")
  expect_gt(n_mult, 0)
  expect_lt(n_mult, nrow(res))
})

test_that("detect_multiplets percentages are valid", {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")
  res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                           plotPercent = FALSE, plotTotalReads = FALSE)

  # percentages are in [0, 100] and roughly sum to 100
  expect_true(all(res$pct_human >= 0 & res$pct_human <= 100))
  expect_true(all(res$pct_mouse >= 0 & res$pct_mouse <= 100))
  expect_true(all(abs(res$pct_human + res$pct_mouse - 100) < 0.5))
})

test_that("detect_multiplets writes an output file", {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")
  out_file <- tempfile(fileext = ".csv")
  detect_multiplets(gem_file, out_file,
                    plotPercent = FALSE, plotTotalReads = FALSE)

  expect_true(file.exists(out_file))
})

test_that("detect_multiplets errors on a missing input file", {
  expect_error(
    detect_multiplets("does_not_exist.csv", tempfile(fileext = ".csv"),
                      plotPercent = FALSE, plotTotalReads = FALSE),
    "not found"
  )
})

test_that("adjusting thresholds changes the number of multiplets", {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")

  res_default <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                                   plotPercent = FALSE, plotTotalReads = FALSE)
  res_strict  <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                                   plotPercent = FALSE, plotTotalReads = FALSE,
                                   overlapDrop = 5, modeDiff = 0.5)

  n_default <- sum(res_default$our_classification == "Multiplet")
  n_strict  <- sum(res_strict$our_classification == "Multiplet")

  # stricter stopping should not find more multiplets than the default
  expect_lte(n_strict, n_default)
})
