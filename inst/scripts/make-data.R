# make-data.R
#
# Provenance of inst/extdata/PC65_gem_classification.csv
#
# This script documents how the bundled example dataset was created, so the
# data in inst/extdata/ can be reproduced. The example is a real
# patient-derived xenograft (PDX) prostate sample (PC65), the same sample shown
# in Figure 1D and 1E of the accompanying manuscript.
#
# The bundled file is a Cell Ranger-style GEM classification table with four
# columns (barcode, GRCh38, GRCm39, call): the human read count, the mouse read
# count, and Cell Ranger's own per-barcode call. It is derived from the
# per-barcode summary produced during the analysis of the PC65 sample.
#
# The starting file is a per-barcode summary with, among other columns:
#   barcode      - the 10x cell barcode (with a -1 suffix)
#   GRCh38       - human read count
#   mouse_reads  - mouse read count
#   call         - Cell Ranger's classification (GRCh38 / GRCm39 / Multiplet)
# together with additional analysis columns that the package does not need.
#
# The steps below keep only the four columns detect_multiplets() uses and
# rename the mouse column to the reference-build name (GRCm39) so the file
# matches a standard Cell Ranger gem_classification.csv.

# ---- 1. Read the per-barcode summary for PC65 -------------------------------
# (path on the analysis machine; adjust as needed to reproduce)
src <- file.path(
  "scrublet_results_CR10_seed42_all18",
  "results_scrublet_CR10_seed42",
  "VCU-PC-065_GRCm39",
  "scrublet_percent_plot_data.csv"
)
raw <- read.csv(src, stringsAsFactors = FALSE)

# ---- 2. Keep the four columns the package needs -----------------------------
# rename mouse_reads -> GRCm39 so the columns match a Cell Ranger
# gem_classification.csv (barcode, GRCh38, GRCm39, call)
pc65 <- data.frame(
  barcode = raw$barcode,
  GRCh38  = raw$GRCh38,
  GRCm39  = raw$mouse_reads,
  call    = raw$call,
  stringsAsFactors = FALSE
)

# ---- 3. Write the bundled example file --------------------------------------
write.csv(
  pc65,
  file.path("inst", "extdata", "PC65_gem_classification.csv"),
  row.names = FALSE
)
