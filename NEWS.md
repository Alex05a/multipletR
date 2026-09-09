# multipletR 0.99.5

* Vignette: the Seurat and `SingleCellExperiment` examples now run (removed
  `eval = FALSE`), using a small self-contained toy object, and show that the
  two object types give an identical classification. Also noted
  `Seurat::as.SingleCellExperiment()` for conversion.
* `remove_multiplets()`: simplified the metadata-writing code using the
  `x$<name>` accessor, which works for both Seurat and SingleCellExperiment.
* Started tracking changes in this NEWS file.

# multipletR 0.99.4

* Removed the placeholder `inst/CITATION` file, since Bioconductor guidance is
  to include a CITATION only when there is an associated preprint or
  publication (will be restored with a DOI once published).
* Added the `Classification` biocView and expanded the package Description.

# multipletR 0.99.1

* `remove_multiplets()`: renamed from `remove_multiplets_seurat()`; now
  auto-detects whether the input is a `Seurat` or `SingleCellExperiment` object
  instead of requiring an `object` argument.
* `remove_multiplets()`: `remove = FALSE` is now the default (annotate only);
  added a `verbose` argument and argument validation.
* `detect_multiplets()`: added a `verbose` argument gating progress messages.
* Added a package-level help page (`?multipletR`).
* Vignette: switched to `BiocStyle::html_document`, added package hyperlinks
  (`CRANpkg()`/`Biocpkg()`/`Githubpkg()`), and added author/date.
* Updated installation instructions to use `BiocManager::install()`.
* Added the `fnd` role and funding information to `Authors@R`.
* Added unit tests for the plotting functions.

# multipletR 0.99.0

* Initial Bioconductor submission.
* `detect_multiplets()`: adaptive threshold detection of human-mouse multiplets
  from a 10x Cell Ranger GEM classification file, with diagnostic plots.
* `remove_multiplets_seurat()`: annotate a Seurat object with the multiplet
  classification and optionally remove the detected multiplets.
* Vignette illustrating the workflow on an example PDX dataset.
