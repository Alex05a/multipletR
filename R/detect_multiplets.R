#' Detect human-mouse multiplets in PDX single-cell data
#'
#' Reads a 10x CellRanger GEM classification file, applies the adaptive
#' threshold method to identify human-mouse multiplets without assuming a fixed
#' species proportion, and writes the input data back out with added
#' classification columns. Optionally produces diagnostic plots.
#'
#' The three starting thresholds T1, T2 and T3 follow the notation used in the
#' manuscript. They define the conservative starting region that is then
#' expanded adaptively; the defaults are the values we recommend and they
#' rarely need to be changed.
#'
#' @param fileIn Path to the input GEM classification CSV (from 10x CellRanger).
#'   Expected columns: a barcode column, a human read-count column (GRCh38),
#'   a mouse read-count column (e.g. GRCm39 / mm10 / mm39 / mouse_reads), and a
#'   10x call column.
#' @param fileOut Path to write the output CSV. The output is the input data
#'   with added columns: our_classification (Multiplet / Singlet), pct_human,
#'   and pct_mouse.
#' @param T1 Upper percent-mouse threshold: the starting upper bound on percent
#'   mouse content. Cells above this are too mouse-heavy. Default 70.
#' @param T2 Lower percent-mouse threshold: the starting lower bound on percent
#'   mouse content (equivalently, the percent-human side). Cells below this are
#'   too human-heavy. Default 30.
#' @param T3 Total-reads threshold: the starting lower bound on total reads,
#'   given as a percentile of the total-reads range (e.g. 25 = 25th percentile).
#'   Default 25.
#' @param plotPercent Logical, whether to draw the two percent plots (one for
#'   the original 10x classification, one for our classification), total reads
#'   vs percent mouse. Default TRUE.
#' @param plotTotalReads Logical, whether to draw the total-reads plots
#'   (mouse reads vs human reads, colored by 10x and our classification).
#'   Default FALSE.
#' @param overlapDrop Stop expanding a threshold when the overlap between the
#'   human and mouse read distributions drops by more than this percentage from
#'   its running peak. Lower values are more conservative. Default 10.
#' @param modeDiff Stop expanding when the difference between the human and
#'   mouse distribution modes increases by more than this amount from the
#'   previous step. Default 0.9.
#' @param verbose Logical, print the full step-by-step algorithm trace. The
#'   final thresholds and multiplet count are always printed; verbose adds the
#'   detailed expansion log. Default FALSE.
#'
#' @return Invisibly, a data frame: the input data with added columns
#'   our_classification, pct_human, and pct_mouse. Also written to \code{fileOut}.
#'
#' @export
detect_multiplets <- function(fileIn,
                              fileOut,
                              T1 = 70,
                              T2 = 30,
                              T3 = 25,
                              plotPercent    = TRUE,
                              plotTotalReads = FALSE,
                              overlapDrop    = 10,
                              modeDiff       = 0.9,
                              verbose        = FALSE) {

  # -- 1. Read and validate the input file ----------------------------------
  if (!file.exists(fileIn)) stop("Input file not found: ", fileIn)
  raw <- utils::read.csv(fileIn, stringsAsFactors = FALSE, check.names = FALSE)

  # locate the barcode, human, mouse, and call columns
  barcode_col <- intersect(c("barcode", "Barcode", "barcodes"), names(raw))[1]
  human_col   <- intersect(c("GRCh38", "hg38", "human_reads", "human"), names(raw))[1]
  mouse_col   <- intersect(c("mouse_reads", "GRCm39", "mm10", "mm39", "mouse"), names(raw))[1]
  call_col    <- intersect(c("call", "Call", "classification"), names(raw))[1]

  if (is.na(human_col) || is.na(mouse_col))
    stop("Could not find human (GRCh38) and mouse read-count columns in ", fileIn)
  if (is.na(barcode_col)) {
    raw$barcode <- seq_len(nrow(raw))   # fall back to row index if no barcode
    barcode_col <- "barcode"
  }

  # -- 2. Build the working data frame the engine expects -------------------
  cc <- data.frame(
    barcode   = raw[[barcode_col]],
    human_10x = as.numeric(raw[[human_col]]),
    mouse_10x = as.numeric(raw[[mouse_col]]),
    stringsAsFactors = FALSE
  )
  cc$call_10x      <- if (!is.na(call_col)) raw[[call_col]] else NA_character_
  cc$total_reads   <- cc$human_10x + cc$mouse_10x
  cc$pct_mouse_10x <- cc$mouse_10x / cc$total_reads * 100
  cc$pct_human_10x <- cc$human_10x / cc$total_reads * 100

  # normalize the 10x call labels (used only for the comparison plot)
  if (!all(is.na(cc$call_10x))) {
    cc$call_10x[cc$call_10x == human_col] <- "Human"
    cc$call_10x[!cc$call_10x %in% c("Human", "Multiplet")] <- "Mouse"
  }

  # -- 3. Run the adaptive threshold engine ---------------------------------
  # T1 = upper percent mouse, T2 = lower percent mouse, T3 = lower reads percentile.
  res <- find_adaptive_thresholds_3t(
    data                         = cc,
    percent_lower                = T2,
    percent_upper                = T1,
    reads_percentile             = T3,
    overlap_drop_threshold       = overlapDrop,
    mode_diff_increase_threshold = modeDiff,
    verbose                      = verbose
  )

  # -- 4. Add our classification columns -------------------------------------
  multiplet_barcodes <- res$detected_multiplets$barcode
  out <- raw
  out$our_classification <- ifelse(cc$barcode %in% multiplet_barcodes,
                                   "Multiplet", "Singlet")
  out$pct_human <- round(cc$pct_human_10x, 2)
  out$pct_mouse <- round(cc$pct_mouse_10x, 2)

  # -- 5. Write the output file ---------------------------------------------
  utils::write.csv(out, fileOut, row.names = FALSE)

  # -- 6. Always print the final thresholds and multiplet count -------------
  th <- res$thresholds
  n_mult <- sum(out$our_classification == "Multiplet")
  message("Final thresholds: T1 (upper % mouse) = ", th$t1_mouse_pct, "%, ",
          "T2 (lower % mouse) = ", th$t2_human_pct, "%, ",
          "T3 (lower reads) = ", round(th$t3_lower_reads), ". ",
          "Detected ", n_mult, " multiplets. Wrote ", fileOut, ".")

  # -- 7. Optional plots -----------------------------------------------------
  cc$our_classification <- out$our_classification
  if (isTRUE(plotPercent))    print(.plot_percent(cc))
  if (isTRUE(plotTotalReads)) print(.plot_total_reads(cc))

  invisible(out)
}
