# Tests for the internal plotting helpers .plot_percent and .plot_total_reads.
# These build a small valid data frame and check that each helper returns a
# ggplot/patchwork object without error, covering plotting_functions.R.

# Helper: a minimal valid `cc` data frame with all columns both plots need.
make_cc <- function(n = 30) {
  set.seed(1)
  data.frame(
    total_reads        = sample(500:5000, n, replace = TRUE),
    human_10x          = sample(100:3000, n, replace = TRUE),
    mouse_10x          = sample(100:3000, n, replace = TRUE),
    pct_mouse_10x      = runif(n, 0, 100),
    call_10x           = sample(c("Human", "Mouse", "Multiplet"), n, replace = TRUE),
    our_classification = sample(c("Singlet", "Multiplet"), n, replace = TRUE),
    stringsAsFactors   = FALSE
  )
}

test_that(".plot_percent returns a ggplot/patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  p <- multipletR:::.plot_percent(make_cc())
  expect_s3_class(p, "ggplot")
})

test_that(".plot_total_reads returns a ggplot/patchwork object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  p <- multipletR:::.plot_total_reads(make_cc())
  expect_s3_class(p, "ggplot")
})

test_that("the plotting helpers run on the bundled example data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR"
  )
  res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                           plotPercent = FALSE, plotTotalReads = FALSE
  )
  raw <- read.csv(gem_file, check.names = FALSE)
  cc <- data.frame(
    human_10x          = as.numeric(raw$GRCh38),
    mouse_10x          = as.numeric(raw$GRCm39),
    call_10x           = raw$call,
    our_classification = res$our_classification,
    stringsAsFactors   = FALSE
  )
  cc$total_reads   <- cc$human_10x + cc$mouse_10x
  cc$pct_mouse_10x <- cc$mouse_10x / cc$total_reads * 100
  expect_s3_class(multipletR:::.plot_percent(cc), "ggplot")
  expect_s3_class(multipletR:::.plot_total_reads(cc), "ggplot")
})
