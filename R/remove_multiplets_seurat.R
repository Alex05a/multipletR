#' Annotate (and optionally remove) multiplets in a Seurat object
#'
#' Adds multipletR's per-cell classification to a Seurat object's metadata and,
#' by default, removes the detected multiplets so the object is ready for
#' downstream analysis. Set \code{remove = FALSE} to keep all cells and only
#' annotate them (useful for visualizing where the multiplets fall, e.g. on a
#' UMAP colored by \code{multipletR_class} before deciding whether to filter).
#'
#' The classification follows the same idea as doublet-detection tools such as
#' DoubletFinder: each cell is labeled Human, Mouse, or Multiplet, and the
#' percent human / percent mouse are added per cell.
#'
#' @param seurat_obj A Seurat object whose cell names (colnames) are barcodes.
#' @param multiplets The data frame returned by \code{\link{detect_multiplets}}
#'   (must contain \code{barcode}, \code{our_classification}, \code{pct_human},
#'   and \code{pct_mouse}).
#' @param remove Logical. If \code{TRUE} (the default), the detected multiplets
#'   are removed from the returned object. If \code{FALSE}, all cells are kept and
#'   only the metadata is added.
#' @param barcode_col Name of the barcode column in \code{multiplets}. Default
#'   "barcode".
#'
#' @return The Seurat object with three added metadata columns
#'   (\code{multipletR_class}, \code{multipletR_pct_human},
#'   \code{multipletR_pct_mouse}). If \code{remove = TRUE}, the multiplet cells
#'   are also dropped.
#'
#' @examples
#' \dontrun{
#'   res <- detect_multiplets("gem_classification.csv", "out.csv")
#'
#'   # Default: annotate and remove multiplets, ready for downstream analysis
#'   seu_clean <- remove_multiplets_seurat(seu, res)
#'
#'   # Annotate only, keeping all cells (e.g. to visualize on a UMAP first)
#'   seu <- remove_multiplets_seurat(seu, res, remove = FALSE)
#'   DimPlot(seu, group.by = "multipletR_class")
#' }
#' @export
remove_multiplets_seurat <- function(seurat_obj,
                                     multiplets,
                                     remove      = TRUE,
                                     barcode_col = "barcode") {
  # -- validate inputs -------------------------------------------------------
  if (!is.data.frame(multiplets))
    stop("'multiplets' must be the data frame returned by detect_multiplets().")
  needed <- c(barcode_col, "our_classification", "pct_human", "pct_mouse")
  miss   <- needed[!needed %in% names(multiplets)]
  if (length(miss) > 0)
    stop("Missing column(s) in 'multiplets': ", paste(miss, collapse = ", "),
         ". Did you pass the output of detect_multiplets()?")

  # -- derive the 3-way class (Human / Mouse / Multiplet) --------------------
  cls <- ifelse(
    multiplets$our_classification == "Multiplet", "Multiplet",
    ifelse(multiplets$pct_human >= multiplets$pct_mouse, "Human", "Mouse")
  )

  # -- match to the Seurat object's cells by barcode -------------------------
  cells <- colnames(seurat_obj)
  idx   <- match(cells, multiplets[[barcode_col]])   # NA where no match

  n_match <- sum(!is.na(idx))
  if (n_match == 0)
    warning("None of the Seurat cell names matched the multiplet barcodes. ",
            "Check that barcodes are in the same format (e.g. '-1' suffixes).")

  # -- add metadata (NA for cells not found in the multiplets table) ---------
  seurat_obj$multipletR_class     <- cls[idx]
  seurat_obj$multipletR_pct_human <- multiplets$pct_human[idx]
  seurat_obj$multipletR_pct_mouse <- multiplets$pct_mouse[idx]

  message("Annotated ", n_match, " of ", length(cells), " cells. ",
          "multipletR_class counts: ",
          paste(names(table(seurat_obj$multipletR_class)),
                table(seurat_obj$multipletR_class), sep = "=", collapse = ", "))

  # -- optionally remove the multiplets --------------------------------------
  if (isTRUE(remove)) {
    keep_cells <- colnames(seurat_obj)[
      is.na(seurat_obj$multipletR_class) |
        seurat_obj$multipletR_class != "Multiplet"]
    n_removed <- length(cells) - length(keep_cells)
    seurat_obj <- seurat_obj[, keep_cells]
    message("Removed ", n_removed, " multiplet cells; ",
            length(keep_cells), " cells remain.")
  }

  seurat_obj
}
