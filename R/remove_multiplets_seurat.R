#' Remove detected multiplets from a Seurat object
#'
#' Convenience function to subset a Seurat object to the singlet barcodes,
#' removing the barcodes that \code{\link{detect_multiplets}} flagged as
#' multiplets. This is a thin wrapper around Seurat's subsetting so users do not
#' have to match barcodes by hand.
#'
#' @param seurat_obj A Seurat object whose cell names (colnames) are barcodes.
#' @param multiplets Either (a) the data frame returned by
#'   \code{detect_multiplets} (which has an \code{our_classification} column and
#'   a barcode column), or (b) a character vector of multiplet barcodes to remove.
#' @param barcode_col Name of the barcode column when \code{multiplets} is a data
#'   frame. Default "barcode".
#' @param classification_col Name of the classification column when
#'   \code{multiplets} is a data frame. Default "our_classification".
#'
#' @return The Seurat object subset to the cells that are NOT multiplets.
#'
#' @examples
#' \dontrun{
#'   res <- detect_multiplets("gem_classification.csv", "out.csv")
#'   seu <- remove_multiplets_seurat(seu, res)
#' }
#' @export
remove_multiplets_seurat <- function(seurat_obj,
                                     multiplets,
                                     barcode_col        = "barcode",
                                     classification_col = "our_classification") {

  # -- work out the list of multiplet barcodes -------------------------------
  if (is.data.frame(multiplets)) {
    if (!barcode_col %in% names(multiplets))
      stop("Column '", barcode_col, "' not found in the multiplets data frame.")
    if (classification_col %in% names(multiplets)) {
      multiplet_barcodes <- multiplets[[barcode_col]][
        multiplets[[classification_col]] == "Multiplet"]
    } else {
      # no classification column: assume every row is a multiplet barcode
      multiplet_barcodes <- multiplets[[barcode_col]]
    }
  } else if (is.character(multiplets)) {
    multiplet_barcodes <- multiplets
  } else {
    stop("'multiplets' must be a data frame (from detect_multiplets) or a ",
         "character vector of barcodes.")
  }

  # -- keep the cells that are NOT multiplets ---------------------------------
  all_cells   <- colnames(seurat_obj)
  keep_cells  <- setdiff(all_cells, multiplet_barcodes)

  n_removed <- length(all_cells) - length(keep_cells)
  if (n_removed == 0)
    warning("No cells were removed - none of the multiplet barcodes matched ",
            "the Seurat object's cell names. Check that the barcodes are in the ",
            "same format (e.g. suffixes like '-1').")

  seurat_obj[, keep_cells]
}
