# Helper: run detect_multiplets on the bundled example and return the result.
get_multiplets <- function() {
  gem_file <- system.file("extdata", "PC65_gem_classification.csv",
                          package = "multipletR")
  detect_multiplets(gem_file, tempfile(fileext = ".csv"),
                    plotPercent = FALSE, plotTotalReads = FALSE)
}

# Helper: build a minimal counts matrix whose colnames are the barcodes.
make_counts <- function(res) {
  matrix(rpois(20 * nrow(res), 5), nrow = 20,
         dimnames = list(paste0("gene", 1:20), res$barcode))
}

# ---- Seurat branch ---------------------------------------------------------

test_that("remove_multiplets (seurat) adds the metadata columns", {
  skip_if_not_installed("Seurat")
  res <- get_multiplets()
  seu <- Seurat::CreateSeuratObject(make_counts(res))

  seu <- remove_multiplets(seu, res, object = "seurat", remove = FALSE)

  expect_true(all(c("multipletR_class", "multipletR_pct_human",
                    "multipletR_pct_mouse") %in% names(seu[[]])))
})

test_that("remove_multiplets (seurat) remove = FALSE keeps all cells", {
  skip_if_not_installed("Seurat")
  res <- get_multiplets()
  seu <- Seurat::CreateSeuratObject(make_counts(res))
  n_before <- ncol(seu)

  seu <- remove_multiplets(seu, res, object = "seurat", remove = FALSE)

  expect_equal(ncol(seu), n_before)
  expect_true(all(seu$multipletR_class %in% c("Human", "Mouse", "Multiplet")))
})

test_that("remove_multiplets (seurat) remove = TRUE drops the multiplets", {
  skip_if_not_installed("Seurat")
  res <- get_multiplets()
  seu <- Seurat::CreateSeuratObject(make_counts(res))
  n_mult <- sum(res$our_classification == "Multiplet")
  n_before <- ncol(seu)

  seu_clean <- remove_multiplets(seu, res, object = "seurat", remove = TRUE)

  expect_equal(ncol(seu_clean), n_before - n_mult)
  expect_false("Multiplet" %in% seu_clean$multipletR_class)
})

# ---- SingleCellExperiment branch -------------------------------------------

test_that("remove_multiplets (sce) adds the metadata columns", {
  skip_if_not_installed("SingleCellExperiment")
  res <- get_multiplets()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = make_counts(res)))

  sce <- remove_multiplets(sce, res, object = "sce", remove = FALSE)

  cd <- SummarizedExperiment::colData(sce)
  expect_true(all(c("multipletR_class", "multipletR_pct_human",
                    "multipletR_pct_mouse") %in% names(cd)))
})

test_that("remove_multiplets (sce) remove = FALSE keeps all cells", {
  skip_if_not_installed("SingleCellExperiment")
  res <- get_multiplets()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = make_counts(res)))
  n_before <- ncol(sce)

  sce <- remove_multiplets(sce, res, object = "sce", remove = FALSE)

  expect_equal(ncol(sce), n_before)
  cls <- SummarizedExperiment::colData(sce)$multipletR_class
  expect_true(all(cls %in% c("Human", "Mouse", "Multiplet")))
})

test_that("remove_multiplets (sce) remove = TRUE drops the multiplets", {
  skip_if_not_installed("SingleCellExperiment")
  res <- get_multiplets()
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = make_counts(res)))
  n_mult <- sum(res$our_classification == "Multiplet")
  n_before <- ncol(sce)

  sce_clean <- remove_multiplets(sce, res, object = "sce", remove = TRUE)

  expect_equal(ncol(sce_clean), n_before - n_mult)
  cls <- SummarizedExperiment::colData(sce_clean)$multipletR_class
  expect_false("Multiplet" %in% cls)
})

# ---- shared validation -----------------------------------------------------

test_that("remove_multiplets errors on a non-data-frame multiplets arg", {
  skip_if_not_installed("Seurat")
  res <- get_multiplets()
  seu <- Seurat::CreateSeuratObject(make_counts(res))

  expect_error(
    remove_multiplets(seu, "not a data frame", object = "seurat"),
    "must be the data frame"
  )
})

test_that("remove_multiplets warns when no barcodes match", {
  skip_if_not_installed("Seurat")
  res <- get_multiplets()
  seu <- Seurat::CreateSeuratObject(make_counts(res))
  bad_res <- res
  bad_res$barcode <- paste0("nomatch_", seq_len(nrow(bad_res)))

  expect_warning(
    remove_multiplets(seu, bad_res, object = "seurat", remove = FALSE),
    "None of the cell names matched"
  )
})
