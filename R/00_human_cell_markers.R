# Human cell-marker panel used for annotation of the GSE144568 CD34+ dataset.
# Reconstructed from the human marker groups retained in Hoffman_violinplot_A1A2.R.
# This is intentionally not the mouse marker file that was present in older working code.

human_cell_markers <- list(
  HSC_MPP = c("AVP", "MLLT3", "MEIS1", "HOPX", "GATA2", "KIT", "PROM1"),
  LMPP_CLP = c("IL7R", "DNTT", "RAG1", "RAG2", "CD79A", "PAX5", "TCF3", "BCL11A"),
  Mono_DC = c("CSF1R", "LST1", "S100A8", "S100A9", "LGALS3", "SPI1", "CTSS", "LYZ"),
  Granulocyte = c("CSF3R", "AZU1", "LTF", "PRTN3", "CXCR2", "S100A8", "S100A9"),
  Megakaryocyte = c("PPBP", "GP9", "GP1BB", "ITGB3", "THBS1", "NRGN"),
  Erythroid = c("HBB", "HBA1", "HBA2", "SLC4A1", "ALAS2", "TFRC"),
  Mast_Basophil = c("TPSAB1", "TPSB2", "KIT", "HDC", "MS4A2", "HPGDS"),
  NK_T = c("NKG7", "GNLY", "PRF1", "GZMH", "KLRD1", "KLRB1", "TRAC"),
  ISG_Progenitor = c("IFI6", "ISG15", "MX1", "IFIT1", "IFIT2", "IFIT3")
)

cellGW <- unique(unlist(human_cell_markers, use.names = FALSE))
