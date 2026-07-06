#' Annotate a Seurat object with multipletR classifications
#'
#' Adds multipletR's per-cell classification and species percentages to a Seurat
#' object's metadata, following the pattern of doublet-detection tools such as
#' DoubletFinder (which annotate rather than remove). Cells are labeled Human,
#' Mouse, or Multiplet, and percent human / percent mouse are added per cell. The
#' object is returned unchanged except for the added metadata, so the user
#' decides whether and how to filter.
#'
#' @param seurat_obj A Seurat object whose cell names (colnames) are barcodes.
#' @param multiplets The data frame returned by \code{\link{detect_multiplets}}
#'   (must contain \code{barcode}, \code{our_classification}, \code{pct_human},
#'   and \code{pct_mouse}).
#' @param barcode_col Name of the barcode column. Default "barcode".
#'
#' @return The Seurat object with three added metadata columns:
#'   \code{multipletR_class} (Human/Mouse/Multiplet),
#'   \code{multipletR_pct_human}, and \code{multipletR_pct_mouse}. Cells present
#'   in the object but absent from \code{multiplets} receive NA.
#'
#' @examples
#' \dontrun{
#'   res <- detect_multiplets("gem_classification.csv", "out.csv")
#'   seu <- annotate_multiplets_seurat(seu, res)
#'   table(seu$multipletR_class)
#'   # then the user may filter if they wish:
#'   seu_clean <- seu[, seu$multipletR_class != "Multiplet"]
#' }
#' @export
annotate_multiplets_seurat <- function(seurat_obj,
                                       multiplets,
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
  seurat_obj$multipletR_class      <- cls[idx]
  seurat_obj$multipletR_pct_human  <- multiplets$pct_human[idx]
  seurat_obj$multipletR_pct_mouse  <- multiplets$pct_mouse[idx]

  message("Annotated ", n_match, " of ", length(cells), " cells. ",
          "multipletR_class counts: ",
          paste(names(table(seurat_obj$multipletR_class)),
                table(seurat_obj$multipletR_class), sep = "=", collapse = ", "))

  seurat_obj
}
