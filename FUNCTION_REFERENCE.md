# multipletR — function reference

Documentation for every function in the package, grouped by role. Two functions
are **exported** (part of the public API); the rest are **internal** — they run
the algorithm but are not called directly by users. Signatures and defaults match
the source as of package version 0.1.0.

Dependencies: `diptest`, `ggplot2`, `patchwork`, `stats`, `utils` (Imports);
`Seurat` (Suggests, only for the Seurat helper).

---

## Public API

### `detect_multiplets()`
**The main entry point.** Reads a Cell Ranger GEM classification file, runs the
adaptive threshold engine, writes an annotated CSV, and optionally draws
diagnostic plots.

```r
detect_multiplets(fileIn, fileOut,
                  T1 = 70, T2 = 30, T3 = 25,
                  plotPercent = TRUE, plotTotalReads = FALSE,
                  overlapDrop = 10, modeDiff = 0.9, verbose = FALSE)
```

**Arguments**

- `fileIn` — path to the input GEM classification CSV. Must contain a human
  read-count column and a mouse read-count column (see column-name flexibility in
  the input docs); a barcode column and a `call` column are used if present.
- `fileOut` — path where the annotated CSV is written.
- `T1` — starting **upper** bound on percent mouse (default 70). Barcodes above
  this are too mouse-dominated to start as multiplets.
- `T2` — starting **lower** bound on percent mouse / upper bound on percent human
  (default 30). Barcodes below this are too human-dominated.
- `T3` — starting **lower total-reads** bound, given as a **percentile** of the
  total-reads range (default 25 = 25th percentile). Excludes low-read droplets.
- `plotPercent` — draw the two "total reads vs percent mouse" panels (10x vs
  ours). Default `TRUE`.
- `plotTotalReads` — draw the two "mouse reads vs human reads" panels. Default
  `FALSE`.
- `overlapDrop` — stop expanding a threshold when the human/mouse distribution
  overlap falls more than this many percentage points below its running peak.
  Lower = more conservative. Default 10.
- `modeDiff` — stop expanding when the difference between the human and mouse
  distribution modes grows by more than this from the previous step. Default 0.9.
- `verbose` — print the full step-by-step expansion trace. The final thresholds
  and multiplet count are always printed regardless. Default `FALSE`.

**Returns** — invisibly, a data frame: the original input plus three added
columns — `our_classification` (`"Multiplet"` / `"Singlet"`), `pct_human`, and
`pct_mouse` (each rounded to 2 dp). The same data frame is written to `fileOut`.

**Behavior notes**

- Column detection is automatic; if no barcode column exists, the row index is
  used as a fallback barcode.
- The `call` labels are normalized to `Human` / `Mouse` / `Multiplet` for the
  comparison plot only — they do not affect detection.
- A one-line summary (`Final thresholds… Detected N multiplets. Wrote <file>.`)
  is always emitted via `message()`.

---

### `remove_multiplets_seurat()`
**Seurat convenience wrapper.** Subsets a Seurat object to its singlet barcodes,
dropping the barcodes flagged as multiplets — so you don't match barcodes by hand.

```r
remove_multiplets_seurat(seurat_obj, multiplets,
                         barcode_col = "barcode",
                         classification_col = "our_classification")
```

**Arguments**

- `seurat_obj` — a Seurat object whose cell names (`colnames`) are barcodes.
- `multiplets` — either **(a)** the data frame returned by `detect_multiplets()`
  (must have a barcode column and, ideally, `our_classification`), or **(b)** a
  plain character vector of multiplet barcodes.
- `barcode_col` — barcode column name when `multiplets` is a data frame. Default
  `"barcode"`.
- `classification_col` — classification column name when `multiplets` is a data
  frame. Default `"our_classification"`. If this column is absent, **every row**
  in `multiplets` is treated as a multiplet barcode.

**Returns** — the Seurat object subset to the non-multiplet cells.

**Behavior notes**

- If no barcodes match (e.g. a `-1` suffix mismatch), the object is returned
  unchanged with a **warning** — check barcode formatting.
- The README describes an annotate-only mode (`remove = FALSE`) that writes a
  `multipletR_class` metadata column and keeps all cells; if your installed copy
  does not expose `remove`/metadata annotation, annotate manually by joining
  `our_classification` onto `seurat_obj@meta.data` by barcode. *(Flag for
  Alexandra: confirm the installed signature matches the README, which advertises
  `remove =` and a `multipletR_class` column.)*

---

## Internal engine

### `find_adaptive_thresholds_3t()`  *(internal)*
The algorithm itself — a round-robin expansion of three thresholds from the
conservative starting region until the human and mouse read distributions begin
to separate. Called by `detect_multiplets()`.

```r
find_adaptive_thresholds_3t(data,
    percent_lower = 30, percent_upper = 70, reads_percentile = 25,
    step_size_reads = 500, step_size_pct = 1,
    overlap_drop_threshold = 10, mode_diff_increase_threshold = 0.9,
    overlap_patience = 3, min_cells_for_bimodality = 30, verbose = FALSE)
```

- `data` — working data frame with `total_reads`, `pct_mouse_10x`, `human_10x`,
  `mouse_10x`.
- `percent_lower` / `percent_upper` — T2 / T1 starting percent-mouse bounds.
- `reads_percentile` — T3 starting lower-reads percentile.
- `step_size_reads` (500) / `step_size_pct` (1) — per-step expansion increments
  for the reads and percent thresholds.
- `overlap_drop_threshold` (10) / `mode_diff_increase_threshold` (0.9) — stopping
  sensitivities (see `should_stop_expanding`).
- `overlap_patience` (3) — how many consecutive below-peak steps are tolerated
  before stopping.
- `min_cells_for_bimodality` (30) — minimum cells before the dip test is run.
- `verbose` — print a detailed trace.

**Returns** a list: `thresholds`, `initial_thresholds`, `stop_reasons`,
`iterations_count`, `initial_metrics`, `final_metrics`, `detected_multiplets`
(the selected barcodes), `history` (per-step log), and `total_rounds`.

**Logic in brief.** The upper reads boundary is fixed at the maximum. Each round
advances T1 (+1% mouse), T2 (−1% mouse), and T3 (−500 reads) one step each,
recomputing the distribution metrics on the combined selection after every step.
A threshold stops when the selection would become bimodal, when overlap drops past
the per-step or peak-drop (patience) limits, when the mode difference jumps, or
when it hits a data boundary. Expansion ends when all three thresholds have
stopped. *(A safety cap of 1000 rounds prevents infinite loops.)*

---

## Internal helpers (`internal_helpers.R`)

### `initialize_thresholds(data, percent_lower = 30, percent_upper = 70, reads_percentile = 25)`
Builds the conservative starting thresholds: `t1_mouse_pct` (upper % mouse),
`t2_human_pct` (lower % mouse), `t3_lower_reads` (the reads percentile, via
`quantile`), and `upper_reads_fixed` (the max). Returns a thresholds list (with a
few compatibility-alias fields).

### `select_cells_in_rectangle(data, thresholds)`
Returns the subset of barcodes inside the current rectangle: total reads within
`[t3_lower_reads, upper_reads_fixed]` **and** percent mouse within
`[t2_human_pct, t1_mouse_pct]`. This selection is the candidate multiplet set at
each step.

### `expand_threshold(thresholds, which_threshold, step_size, min_reads, max_reads)`
Advances a single threshold by one step: `1` = T1 up (capped at 100%), `2` = T2
down (floored at 0%), `3` = T3 down (floored at `min_reads`). Returns the updated
thresholds list. The upper reads boundary is never moved.

### `calculate_adaptive_metrics(data, min_cells_for_bimodality = 30)`
The quality metrics for a candidate selection. Log-transforms reads
(`log(reads + 1)`), builds kernel-density estimates of the human and mouse
distributions on a common grid (n = 2000), and returns:

- `n_cells` — selection size.
- `pct_overlap` — **symmetric percent overlap** of the two densities (mean of
  overlap/area for each), via the trapezoidal rule. Higher = both species well
  represented.
- `area_overlap` — the raw overlap area.
- `diff_modes` — absolute difference between the two density modes (lower = more
  similar). Modes are taken from the same ranged density (matches the manuscript).
- `human_bimodal` / `mouse_bimodal` — dip-test results per distribution.
- `human_dip_pval` / `mouse_dip_pval`, `human_mode` / `mouse_mode`.

Returns a safe zero/`NA` result when the selection is too small (< 10 cells).

### `should_stop_expanding(current_metrics, previous_metrics, overlap_drop_threshold = 10, mode_diff_increase_threshold = 0.9)`
The stopping rule applied after each step. Returns `list(stop, reason)`. Stops if:
the selection is unimodal-eligible (≥ 30 cells) **and** either distribution is
bimodal (dip p < 0.05); or overlap dropped more than `overlap_drop_threshold` from
the previous step; or the mode difference increased more than
`mode_diff_increase_threshold`. Below 30 cells it never stops ("need more cells").

### `test_bimodality(x, alpha = 0.05, min_cells_for_test = 30)`
Wraps **Hartigan's dip test** (`diptest::dip.test`). Returns `is_bimodal`,
`p_value`, `dip_statistic`, and a human-readable `reason`. Returns not-bimodal if
there are fewer than `min_cells_for_test` finite values, and degrades gracefully
on error.

### `find_mode(x)`
Returns the x-position of the kernel-density peak of `x` (n = 512), or `NA` for
fewer than 10 values. *(Utility; the main metric path computes modes from the
ranged density inside `calculate_adaptive_metrics`.)*

---

## Plotting (`plotting_functions.R`, internal)

### `.plot_percent(cc)`
Two side-by-side scatter panels of **total reads (x) vs percent mouse (y)** — one
colored by the 10x `call`, one by our classification. Returns a `patchwork`
object. Drawn by `detect_multiplets()` when `plotPercent = TRUE`. This is the view
where multiplets sit in the central percent-mouse band.

### `.plot_total_reads(cc)`
Two side-by-side scatter panels of **mouse reads (x) vs human reads (y)** — one by
10x `call`, one by ours. Returns a `patchwork` object. Drawn when
`plotTotalReads = TRUE`.

`global.R` / the header of `plotting_functions.R` register the ggplot aesthetic
column names with `utils::globalVariables()` so `R CMD check` does not flag them.

---

## Minimal worked example

```r
library(multipletR)

res <- detect_multiplets(
  fileIn  = "SAMPLE/outs/analysis/gem_classification.csv",
  fileOut = "SAMPLE_multiplets.csv",
  verbose = TRUE            # print the expansion trace
)

table(res$our_classification)
# Singlet Multiplet
#   ...       ...

# Then, with Seurat:
library(Seurat)
seu_clean <- remove_multiplets_seurat(seu, res)   # drop multiplets
```

---

## Verification

The package was verified by running it end-to-end on a real PDX sample. Using
`detect_multiplets()` (defaults, `verbose = TRUE`) on the BC37 Cell Ranger v10 /
GRCm39 GEM classification file (`VCU-BC-037_GRCm39_gem_classification.csv`; 9,735
barcodes), the function loaded the file, ran the adaptive three-threshold engine,
and wrote the annotated output without error. Starting from the conservative
region (T1 = 70%, T2 = 30%, T3 = 6,529 reads), which held 28 cells at 77.9%
human–mouse distribution overlap, the algorithm expanded the thresholds to a final
selection of **72 multiplets** (T1 = 75%, T2 = 25%, T3 = 1,029 reads) with the
distribution overlap **rising to 84.5%** and both the human and mouse read
distributions remaining **unimodal** (Hartigan dip test, p ≥ 0.05). Expansion
stopped exactly as designed: the two percent thresholds (T1, T2) halted on the
overlap peak-drop criterion, and the total-reads threshold (T3) halted when
lowering it further would have made the mouse distribution bimodal. The run
produced the expected outputs — the `our_classification`, `pct_human`, and
`pct_mouse` columns in `BC37_multiplets_out.csv`, the printed
`table(our_classification)` (72 Multiplet / 9,663 Singlet), and the 10X-vs-Ours
diagnostic percent plot. These values match the BC37 result reported in the
manuscript (72 multiplets at 84.5% overlap), confirming the packaged code
reproduces the published behavior on real Cell Ranger output.
