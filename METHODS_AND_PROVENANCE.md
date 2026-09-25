# Methods and file provenance

This repository was pruned against the accepted manuscript rather than against the full historical analysis notebook. Exploratory objects (`obj1`, `obj4`, `obj5`), mouse-oriented analyses, and superseded plotting attempts were intentionally excluded.

## Figure 1

### Figure 1A
GSE144568, all 21 samples (6 healthy controls and 15 MF samples). R/Seurat workflow in `R/01_GSE144568_all21_obj2.R`.

### Figure 1B-F
Untreated/age-matched GSE144568 subset (4 controls and 6 untreated MF samples; 4 JAK2V617F and 2 CALR). R/Seurat workflow in `R/02_GSE144568_untreated_obj3.R`. Clustering uses resolution 0.25.

### Figure 1G
GSEA was run in the Broad Institute GSEA Java desktop application. The four exact ranked-list inputs supplied for the manuscript analysis are deposited under `data/GSEA_ranked_lists/`. The repository does not substitute an R GSEA implementation for the software actually used.

## Figure 2

### Figure 2A
GSE122198/Target-seq data were extracted, organized in Excel, and plotted/analyzed in GraphPad Prism. The supplied historical files support a three-group workflow comprising 678 healthy-control cells, 413 MF JAK2-WT cells, and 1,707 JAK2-mutant MF cells (2,798 cells total). The repository includes the extraction/grouping workbooks and the most recent supplied Prism file whose group labels match this workflow (`WT`, `MF JAK2 WT`, and `MF JAK2 mutant`). Earlier supplied Prism working versions were omitted. A detailed reconstruction of the file lineage, including the distinction between preserved files and reconstructed relationships, is provided in `data/GSE122198/README.md`.

### Figure 2B
Platelet expression values were organized/transformed in Excel and plotted/analyzed in GraphPad Prism. The source and log1p-transformed workbooks are included. A separate final platelet Prism project was not among the supplied files, so no Prism file is claimed for this panel.

### Figure 2C-D
The compact `data/platelet/synergy.xlsx` table and `R/04_platelet_HMGA_biomarker.R` reproduce the median-based expression groups and Fisher-exact biomarker calculations. `R/05_verify_Fig2D_biomarker_statistics.R` is an independent spot-check using the reported contingency-table counts.

## Median convention
The manuscript defines high expression as greater than the cohort median. In `synergy.xlsx` (120 samples), the HMGA1 median is 1187.7435 and the HMGA2 median is 7.078; neither median is an observed value, so `>` and `>=` give the same 60/60 split for each gene.

## Files intentionally not included
- `obj1`, `obj4`, `obj5` analyses and related statistics.
- mouse gene-marker, ImmGen/SingleR, and mouse GSEA code.
- superseded figure versions and exploratory plotting attempts.
- unrelated historical working scripts.
- the accepted manuscript PDF itself.
- external supplementary workbooks from prior publications when the public source can be cited instead.
