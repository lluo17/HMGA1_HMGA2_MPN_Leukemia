# GSE122198 / 3'-TARGET-seq data provenance for Figure 2A

Figure 2A uses the public **GSE122198** 3'-TARGET-seq dataset. The historical workflow combined public single-cell expression values with mutation/diagnosis annotations, organized the resulting table in Excel, separated cells into the groups used for the figure, and then plotted/analyzed the grouped values in GraphPad Prism.

This repository keeps the compact Excel intermediates needed to document that workflow rather than redistributing the large public count matrix.

## Reconstructed provenance

The file history supplied with this repository supports the following workflow:

```text
GSE122198 public data
  |
  |-- GEO sample metadata -> mutation and diagnosis annotations
  |     (documented by R/03_GSE122198_metadata.R)
  |
  |-- 3'-TARGET-seq expression matrix
        |
        |-- retain HMGA1, HMGA2 and JAK2 expression
        |-- combine with mutation + diagnosis metadata
        v
GSE122198_extracted_HMGA1_HMGA2_JAK2.xlsx
        |
        |-- Excel organization / transformation
        |-- group cells as Control, MF JAK2 WT, or MF JAK2 mutant
        v
GSE122198_grouped_HMGA1_HMGA2.xlsx
        |
        v
GraphPad Prism -> Figure 2A
```

The original working files show **2,798 cells**: **678 healthy-control cells** and **2,120 myelofibrosis cells**. The MF cells separate into **413 JAK2-WT cells** and **1,707 JAK2-mutant cells**, matching the three groups represented in the grouped workbook.

## Files deposited here

### `GSE122198_extracted_HMGA1_HMGA2_JAK2.xlsx`

This workbook preserves the expression/metadata extraction stage. Its first sheet contains the cell identifiers together with mutation, diagnosis, HMGA1, HMGA2, JAK2 and CALR rows. A second representation places one cell per row with columns for mutation, diagnosis, HMGA1, HMGA2 and JAK2 and includes historical Excel log2 calculations. A working sheet separates HMGA1/HMGA2 values into normal, MF JAK2-WT and MF JAK2-mutant groups.

This file is retained as a provenance/intermediate workbook; it is not presented as a newly regenerated analysis.

### `GSE122198_grouped_HMGA1_HMGA2.xlsx`

This is the compact analysis-ready workbook used for the final Excel-to-Prism step. Separate `HMGA1` and `HMGA2` sheets contain three columns corresponding to:

- `WT` (healthy-control cells; n = 678)
- `MF JAK2 WT` (MF cells without JAK2 mutation; n = 413)
- `MF JAK2 mut` / `MF JAK2 mutant` (JAK2-mutant MF cells; n = 1,707)

A small summary sheet contains historical cell-count/expression calculations. The workbook is deposited as the analysis input rather than being recreated or reformatted.

## Metadata helper

`R/03_GSE122198_metadata.R` documents retrieval of mutation and diagnosis information from the GEO record using `GEOquery`. It is included to make the metadata provenance explicit. It should not be interpreted as a claim that the complete historical Figure 2A workflow was executed end-to-end in R; the final figure used Excel and GraphPad Prism.

## Public source data

The full 3'-TARGET-seq count matrix is public through **GSE122198** and is therefore not duplicated in this repository. Historical source/intermediate files supplied during repository reconstruction included the public expression matrix and tab-delimited mutation/diagnosis tables; these were used to verify the provenance above but are omitted from the repository to avoid unnecessary duplication of large public data and redundant working files.

## Reproducibility note

This repository documents the analysis materials that can be traced confidently to the accepted manuscript. Historical intermediate filenames and spreadsheet manipulations were not all recorded contemporaneously, so the provenance above distinguishes directly preserved files from reconstructed file-to-file relationships. No attempt was made to replace the historical Excel/GraphPad workflow with a new computational pipeline.
