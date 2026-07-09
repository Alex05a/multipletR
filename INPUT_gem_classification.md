# The input: Cell Ranger's GEM classification file

`multipletR` takes a single file as its primary input — the **GEM classification
CSV** produced by 10x Genomics Cell Ranger when reads are aligned to a combined
two-species reference. This section explains where that file comes from, where to
find it, and what it contains, so that new users understand the one input the
package needs.

## How the data is generated

In a patient-derived xenograft (PDX) experiment, a human tumor grows inside a
mouse host, so a dissociated sample contains **both human tumor cells and mouse
stromal/immune cells**. When this sample is run through droplet-based single-cell
RNA-sequencing (10x Genomics Chromium), each droplet — a *GEM* (Gel Bead-in-
EMulsion) — captures one or more cells together with a barcoded gel bead, so that
every cell's transcripts inherit a shared cell barcode and each molecule gets a
unique molecular identifier (UMI). Most droplets contain a single cell (a
*singlet*), but some capture two or more cells; when those cells come from
different species, the droplet is a **human–mouse multiplet**.

To tell the two species apart, the sequencing reads are aligned with
**`cellranger count` against a combined ("barnyard") reference** that contains
both genomes at once — for example GRCh38 (human) and GRCm39 (mouse). In this
combined reference every gene is prefixed by its genome of origin (e.g.
`GRCh38-EPCAM`, `GRCm39-Col1a1`), so each read can be attributed to human or
mouse. For every cell barcode, Cell Ranger then counts how many reads map to each
genome and, using fixed count thresholds, labels the barcode as human, mouse, or
a multiplet. This per-barcode, two-genome summary is written to the GEM
classification file.

## Where to find it in the Cell Ranger output

The file is written only when the reference is a **multi-genome (barnyard)**
reference. For a run with `--id=<SAMPLE>`, it is located at:

```
<SAMPLE>/outs/analysis/gem_classification.csv
```

(The full expression matrix used elsewhere lives separately at
`<SAMPLE>/outs/filtered_feature_bc_matrix/`. `multipletR` does **not** need the
matrix — only `gem_classification.csv` — although the Seurat helper later joins
the calls back onto a matrix-derived object by barcode.)

## What it contains

One row per cell barcode, with these columns:

| Column | Meaning |
|---|---|
| `barcode` | The 10x cell barcode (often carries a `-1` suffix, e.g. `AAACCTGAG….-1`). |
| `<human genome>` | Per-barcode count of reads/UMIs assigned to the human genome. The column is **named after the reference build**, e.g. `GRCh38`. |
| `<mouse genome>` | Per-barcode count of reads/UMIs assigned to the mouse genome, named after the build: `GRCm39` (2024-A references), or `mm10` / `mm39` for older builds. |
| `call` | Cell Ranger's own classification of the barcode: the human genome name (e.g. `GRCh38`), the mouse genome name (e.g. `GRCm39`), or `Multiplet`. |

A typical file looks like:

```
barcode,GRCh38,GRCm39,call
AAACCTGAGAAACCAT-1,10432,58,GRCh38
AAACCTGAGATCCGAG-1,71,9903,GRCm39
AAACCTGCACGTGTGA-1,4821,5210,Multiplet
```

## Why this is the input `multipletR` needs

The package's whole premise is that a **true human–mouse multiplet has
substantial reads from both genomes**, whereas a singlet is dominated by one. The
GEM classification file provides exactly the two numbers needed to test that —
the human and mouse read counts per barcode — from which `multipletR` derives, for
each barcode, the **total reads** (human + mouse) and the **percent mouse**
(mouse / total × 100). It then applies its adaptive thresholds in that
total-reads × percent-mouse space, rather than trusting Cell Ranger's fixed
`call`. Cell Ranger's `call` column is kept only for comparison plots.

### Column-name flexibility

Because the genome columns are named after whichever reference build was used,
`multipletR` recognizes several common names automatically so you do not have to
rename anything:

- **Human:** `GRCh38`, `hg38`, `human_reads`, `human`
- **Mouse:** `GRCm39`, `mm10`, `mm39`, `mouse_reads`, `mouse`
- **Barcode:** `barcode`, `Barcode`, `barcodes` (if none is present, the row index
  is used as a fallback)
- **Call:** `call`, `Call`, `classification`

As long as the file has a human read-count column and a mouse read-count column,
`detect_multiplets()` will run; the barcode and `call` columns are optional but
recommended (the barcode is what lets you match results back to a Seurat object).

> **Tip — one caveat that bites people:** barcode suffixes must match between the
> GEM file and any Seurat object you later filter. Cell Ranger writes barcodes
> with a trailing `-1`; if your Seurat object's cell names lack it (or vice
> versa), `remove_multiplets_seurat()` will match nothing. Make the two consistent
> before joining.
