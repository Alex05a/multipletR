#' Annotate (and optionally remove) multiplets in a Seurat or SingleCellExperiment object
#'
#' Adds multipletR's per-cell classification to a single-cell object's cell
#' metadata. By default all cells are kept and only annotated (useful for
#' visualizing where the multiplets fall, e.g. on a UMAP colored by
#' \code{multipletR_class} before deciding whether to filter). Set
#' \code{remove = TRUE} to drop the detected multiplets so the object is ready
#' for downstream analysis. The object type (Seurat or
#' \code{SingleCellExperiment}) is detected automatically.
#'
#' The classification follows the same idea as doublet-detection tools such as
#' DoubletFinder: each cell is labeled Human, Mouse, or Multiplet, and the
#' percent human / percent mouse are added per cell.
#'
#' @param x A single-cell object whose cell names (\code{colnames}) are
#'   barcodes: either a Seurat object or a \code{SingleCellExperiment}. The
#'   object type is detected automatically.
#' @param multiplets The data frame returned by \code{\link{detect_multiplets}}
#'   (must contain \code{barcode}, \code{our_classification}, \code{pct_human},
#'   and \code{pct_mouse}).
#' @param remove Logical. If \code{TRUE}, the detected multiplets are removed
#'   from the returned object. If \code{FALSE} (the default), all cells are kept
#'   and only the metadata is added.
#' @param barcode_col Name of the barcode column in \code{multiplets}. Default
#'   "barcode".
#' @param verbose Logical. If \code{TRUE} (the default), print a summary of the
#'   annotation and removal. Set \code{FALSE} to suppress these messages.
#'
#' @return The input object with three added cell-metadata columns
#'   (\code{multipletR_class}, \code{multipletR_pct_human},
#'   \code{multipletR_pct_mouse}). If \code{remove = TRUE}, the multiplet cells
#'   are also dropped.
#'
#' @examplesIf requireNamespace("Seurat", quietly = TRUE)
#' # Detect multiplets in the bundled example dataset
#' gem_file <- system.file("extdata", "PC65_gem_classification.csv",
#'   package = "multipletR"
#' )
#' res <- detect_multiplets(gem_file, tempfile(fileext = ".csv"),
#'   plotPercent = FALSE, plotTotalReads = FALSE
#' )
#'
#' # Build a minimal Seurat object whose cell names are the same barcodes
#' counts <- matrix(rpois(20 * nrow(res), 5),
#'   nrow = 20,
#'   dimnames = list(paste0("gene", 1:20), res$barcode)
#' )
#' seu <- Seurat::CreateSeuratObject(counts)
#'
#' # Annotate only (keep all cells), e.g. to visualize before filtering
#' seu <- remove_multiplets(seu, res)
#' table(seu$multipletR_class)
#'
#' # Or annotate and remove the multiplets in one step
#' seu_clean <- remove_multiplets(seu, res, remove = TRUE)
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
remove_multiplets <- function(x,
                              multiplets,
                              remove = FALSE,
                              barcode_col = "barcode",
                              verbose = TRUE) {
  # -- auto-detect the object type -------------------------------------------
  sce <- methods::is(x, "SingleCellExperiment")
  seu <- methods::is(x, "Seurat")
  if (!sce && !seu) {
    stop("'x' must be a Seurat or SingleCellExperiment object.")
  }
  # -- validate arguments ----------------------------------------------------
  if (!is.data.frame(multiplets)) {
    stop("'multiplets' must be the data frame returned by detect_multiplets().")
  }
  if (!is.logical(remove) || length(remove) != 1L) {
    stop("'remove' must be a single logical value (TRUE or FALSE).")
  }
  if (!is.character(barcode_col) || length(barcode_col) != 1L) {
    stop("'barcode_col' must be a single character string.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a single logical value (TRUE or FALSE).")
  }
  needed <- c(barcode_col, "our_classification", "pct_human", "pct_mouse")
  miss <- needed[!needed %in% names(multiplets)]
  if (length(miss) > 0) {
    stop(
      "Missing column(s) in 'multiplets': ", paste(miss, collapse = ", "),
      ". Did you pass the output of detect_multiplets()?"
    )
  }
  # -- derive the 3-way class (Human / Mouse / Multiplet) --------------------
  cls <- ifelse(
    multiplets$our_classification == "Multiplet", "Multiplet",
    ifelse(multiplets$pct_human >= multiplets$pct_mouse, "Human", "Mouse")
  )
  # -- match the object's cells to the multiplets table by barcode -----------
  cells <- colnames(x)
  idx <- match(cells, multiplets[[barcode_col]]) # NA where no match
  n_match <- sum(!is.na(idx))
  if (n_match == 0) {
    warning(
      "None of the cell names matched the multiplet barcodes. ",
      "Check that barcodes are in the same format (e.g. '-1' suffixes)."
    )
  }
  class_vec <- cls[idx]
  pct_human_vec <- multiplets$pct_human[idx]
  pct_mouse_vec <- multiplets$pct_mouse[idx]
  # -- add the metadata: separate branch per object type ---------------------
  if (seu) {
    x$multipletR_class <- class_vec
    x$multipletR_pct_human <- pct_human_vec
    x$multipletR_pct_mouse <- pct_mouse_vec
    class_in_obj <- x$multipletR_class
  } else {
    SummarizedExperiment::colData(x)$multipletR_class <- class_vec
    SummarizedExperiment::colData(x)$multipletR_pct_human <- pct_human_vec
    SummarizedExperiment::colData(x)$multipletR_pct_mouse <- pct_mouse_vec
    class_in_obj <- SummarizedExperiment::colData(x)$multipletR_class
  }
  if (verbose) {
    message(
      "Annotated ", n_match, " of ", length(cells), " cells. ",
      "multipletR_class counts: ",
      paste(names(table(class_in_obj)),
        table(class_in_obj),
        sep = "=", collapse = ", "
      )
    )
  }
  # -- optionally remove the multiplets --------------------------------------
  if (isTRUE(remove)) {
    keep <- is.na(class_in_obj) | class_in_obj != "Multiplet"
    n_removed <- sum(!keep)
    x <- x[, keep]
    if (verbose) {
      message(
        "Removed ", n_removed, " multiplet cells; ",
        sum(keep), " cells remain."
      )
    }
  }
  x
}
