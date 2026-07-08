# multipletR

Adaptive detection of human–mouse multiplets in patient-derived xenograft (PDX)
single-cell RNA-seq data.

In mixed-species experiments (such as PDX models, where a human tumor grows in a
mouse host), some droplets capture both a human and a mouse cell. These
**multiplets** must be removed before downstream analysis. `multipletR` detects
them using an adaptive threshold method that, unlike fixed cutoffs, does **not**
assume a fixed human/mouse proportion — so it handles the imbalanced species
mixtures typical of real PDX samples.

## Installation

```r
# install.packages("devtools")
devtools::install_github("Alex05a/multipletR")
```

The Seurat helper functions additionally require the `Seurat` package.

## Quick start

The core function, `detect_multiplets()`, reads a 10x Cell Ranger GEM
classification file and returns the same data with our multiplet classification
added.

```r
library(multipletR)

res <- detect_multiplets(
  fileIn  = "gem_classification.csv",   # 10x Cell Ranger output
  fileOut = "multiplets_out.csv"        # where to write the annotated data
)

# The result has the original data plus our classification and percentages:
head(res)
#>   barcode  GRCh38  GRCm39  call       our_classification  pct_human  pct_mouse
#>   ...      ...     ...     GRCh38     Singlet             99.5       0.5
#>   ...      ...     ...     Multiplet  Multiplet           54.8       45.2
```

By default this also draws diagnostic plots (total reads vs. percent mouse),
colored by the 10x classification and by our method. The multiplets our method
finds appear in the central band where human and mouse reads are balanced.

## Using the results with Seurat

`remove_multiplets_seurat()` works on a Seurat object. It adds our per-cell
classification to the object's metadata (Human / Mouse / Multiplet, plus percent
human and percent mouse per cell) and, by default, removes the detected
multiplets so the object is ready for downstream analysis. Set `remove = FALSE`
to keep all cells and only annotate them — useful for visualizing where the
multiplets fall before deciding whether to filter, following the same idea as
tools like DoubletFinder.

```r
library(Seurat)

# Default: annotate and remove the multiplets
seu_clean <- remove_multiplets_seurat(seu, res)

# Or annotate only (keep all cells) to inspect them first:
seu <- remove_multiplets_seurat(seu, res, remove = FALSE)
table(seu$multipletR_class)
#>     Human     Mouse  Multiplet
#>     10937      2432        280

# Visualize where the multiplets fall on a UMAP:
DimPlot(seu, group.by = "multipletR_class")
```

## Adjusting the parameters

The defaults are the values we recommend and rarely need changing, but every
parameter is adjustable — just pass it as a named argument. Anything you do not
specify keeps its default.

```r
res <- detect_multiplets(
  fileIn  = "gem_classification.csv",
  fileOut = "multiplets_out.csv",
  T1          = 70,    # starting upper bound on percent mouse (default 70)
  T2          = 30,    # starting lower bound on percent mouse (default 30)
  T3          = 25,    # starting lower reads bound, as a percentile (default 25)
  overlapDrop = 10,    # stop when overlap drops more than this from its peak (default 10)
  modeDiff    = 0.9    # stop when the mode difference grows more than this (default 0.9)
)
```

Lower `overlapDrop` and `modeDiff` make the method stop expanding sooner (more
conservative, fewer multiplets). See `?detect_multiplets` for full details on
every argument.

## Functions

| Function | Purpose |
|---|---|
| `detect_multiplets()` | Detect multiplets from a Cell Ranger GEM classification file; add our classification and draw diagnostic plots. |
| `remove_multiplets_seurat()` | Annotate a Seurat object with our classification (Human / Mouse / Multiplet and percent human/mouse) and, by default, remove the multiplets. Set `remove = FALSE` to annotate only. |

## How it works

The method defines a conservative starting region using three thresholds — T1
(upper percent mouse), T2 (lower percent mouse), and T3 (lower total reads) —
then expands them step by step to capture additional multiplets. It stops
expanding when the human and mouse read distributions start to look like
singlets: when a distribution becomes bimodal, when their overlap drops, or when
their modes diverge. This lets the multiplet region adapt to each dataset rather
than relying on a fixed cutoff.
