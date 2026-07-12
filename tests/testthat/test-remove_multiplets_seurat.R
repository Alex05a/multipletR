# Helper: build a small Seurat object whose cells match the detect_multiplets
# output barcodes. Skipped automatically when Seurat is not installed.
make_test_objects <- function() {
  gem_file <- system.file("extdata", "12G_gem_classification.csv",
                          package = "multipletR")
  res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                           plotPercent = FALSE, plotTotalReads = FALSE)
  counts <- matrix(rpois(20 * nrow(res), 5), nrow = 20,
                   dimnames = list(paste0("gene", 1:20), res$barcode))
  seu <- Seurat::CreateSeuratObject(counts)
  list(res = res, seu = seu)
}

test_that("remove_multiplets_seurat adds the metadata columns", {
  skip_if_not_installed("Seurat")
  obj <- make_test_objects()

  seu <- remove_multiplets_seurat(obj$seu, obj$res, remove = FALSE)

  expect_true(all(c("multipletR_class", "multipletR_pct_human",
                    "multipletR_pct_mouse") %in% names(seu[[]])))
})

test_that("remove = FALSE keeps all cells and annotates them", {
  skip_if_not_installed("Seurat")
  obj <- make_test_objects()

  n_before <- ncol(obj$seu)
  seu <- remove_multiplets_seurat(obj$seu, obj$res, remove = FALSE)

  # no cells removed
  expect_equal(ncol(seu), n_before)

  # classes are the three expected labels
  expect_true(all(seu$multipletR_class %in% c("Human", "Mouse", "Multiplet")))
})

test_that("remove = TRUE drops the multiplet cells", {
  skip_if_not_installed("Seurat")
  obj <- make_test_objects()

  n_mult <- sum(obj$res$our_classification == "Multiplet")
  n_before <- ncol(obj$seu)

  seu_clean <- remove_multiplets_seurat(obj$seu, obj$res, remove = TRUE)

  # the multiplets are gone
  expect_equal(ncol(seu_clean), n_before - n_mult)
  expect_false("Multiplet" %in% seu_clean$multipletR_class)
})

test_that("remove_multiplets_seurat errors on a non-data-frame input", {
  skip_if_not_installed("Seurat")
  obj <- make_test_objects()

  expect_error(
    remove_multiplets_seurat(obj$seu, "not a data frame"),
    "must be the data frame"
  )
})

test_that("remove_multiplets_seurat warns when no barcodes match", {
  skip_if_not_installed("Seurat")
  obj <- make_test_objects()

  # rename the multiplets' barcodes so nothing matches the Seurat object
  bad_res <- obj$res
  bad_res$barcode <- paste0("nomatch_", seq_len(nrow(bad_res)))

  expect_warning(
    remove_multiplets_seurat(obj$seu, bad_res, remove = FALSE),
    "None of the Seurat cell names matched"
  )
})
