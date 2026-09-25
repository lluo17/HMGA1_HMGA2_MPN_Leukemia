# HMGA1 and HMGA2 in myeloproliferative neoplasms

Analysis materials corresponding to the accepted *Leukemia* manuscript:

**HMGA1 and HMGA2: Epigenetic regulators in myeloproliferative neoplasms and novel platelet transcriptomic biomarkers of myelofibrosis**

This repository reflects the software actually used during manuscript preparation: R/Seurat for the single-cell analyses, the Broad Institute GSEA Java desktop application for pathway enrichment, and Excel/GraphPad Prism for selected expression/statistical panels.

## Repository structure

```text
R/
  00_human_cell_markers.R
  01_GSE144568_all21_obj2.R
  02_GSE144568_untreated_obj3.R
  03_GSE122198_metadata.R
  04_platelet_HMGA_biomarker.R
  05_verify_Fig2D_biomarker_statistics.R

data/
  GSEA_ranked_lists/
  GSE122198/
  platelet/

results/
  Fig2A_GraphPad/
  GSEA/

METHODS_AND_PROVENANCE.md
```

## Figure-to-file map

| Panel | Dataset / workflow | Repository material |
|---|---|---|
| Fig. 1A | GSE144568, 6 controls + 15 MF | `R/01_GSE144568_all21_obj2.R` |
| Fig. 1B-F | GSE144568 untreated subset, 4 controls + 6 MF | `R/02_GSE144568_untreated_obj3.R` |
| Fig. 1G | GSEA Java desktop application | exact `.rnk` inputs in `data/GSEA_ranked_lists/` |
| Fig. 2A | GSE122198/Target-seq; Excel + GraphPad Prism | `R/03_GSE122198_metadata.R`, `data/GSE122198/README.md`, `data/GSE122198/`, `results/Fig2A_GraphPad/` |
| Fig. 2B | platelet transcriptomes; Excel + GraphPad Prism | source + log1p workbooks in `data/platelet/` |
| Fig. 2C-D | platelet biomarker analysis | `data/platelet/synergy.xlsx`, `R/04_platelet_HMGA_biomarker.R` |
| Fig. 2D check | independent calculation verification | `R/05_verify_Fig2D_biomarker_statistics.R` |

## Public datasets

The manuscript analyzes public datasets including **GSE144568** and **GSE122198**. Large raw count matrices and serialized Seurat objects are not included here; obtain the public data from GEO and place the GSE144568 count matrices under `scMatrices/` as described in the R scripts.

## R packages

The cleaned scripts use packages including `Seurat`, `dplyr`, `stringr`, `tidyr`, `ggplot2`, `patchwork`, `scales`, `readxl`, and `GEOquery`.

## GSEA

Figure 1G was run in the Broad Institute GSEA Java application. The repository includes the exact ranked lists supplied for the analysis rather than replacing the historical workflow with a newly written R implementation. Complete GSEA output folders may be deposited under `results/GSEA/` without modifying their internal structure.

## Excel and GraphPad Prism

Figure 2A was prepared from GSE122198/Target-seq data through an Excel-to-Prism workflow. The reconstructed file provenance and the three analysis groups (healthy control, MF JAK2-WT, and MF JAK2-mutant) are documented in `data/GSE122198/README.md`. Figure 2B platelet expression/statistics were also performed in GraphPad Prism after Excel preparation. The compact Excel intermediates that could be traced confidently to the accepted figures are included; superseded working spreadsheets were omitted.

## Reproducibility and cleanup

The original analysis notebook accumulated exploratory work from multiple projects. This public version intentionally excludes unrelated mouse analyses, `obj1`/`obj4`/`obj5`, Harmony tests not used in the manuscript, and earlier versions of figures. See `METHODS_AND_PROVENANCE.md` for details.
