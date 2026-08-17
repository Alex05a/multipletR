#' multipletR: Adaptive Detection of Human-Mouse Multiplets in PDX Single-Cell Data
#'
#' Detects human-mouse multiplets in patient-derived xenograft (PDX) single-cell
#' RNA-seq data using an adaptive threshold method that does not assume a fixed
#' species proportion. It reads a 10x CellRanger GEM classification file and
#' returns the data with added multiplet classifications, with optional
#' diagnostic plots and helpers to annotate or filter a Seurat or
#' SingleCellExperiment object.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{detect_multiplets}}: read a 10x CellRanger GEM
#'     classification file, detect human-mouse multiplets with the adaptive
#'     threshold method, and write out the classified data, with optional
#'     diagnostic plots.
#'   \item \code{\link{remove_multiplets}}: add multipletR's per-cell
#'     classification to a Seurat or SingleCellExperiment object and, optionally,
#'     remove the detected multiplets so the object is ready for downstream
#'     analysis.
#' }
#'
#' @keywords internal
"_PACKAGE"
