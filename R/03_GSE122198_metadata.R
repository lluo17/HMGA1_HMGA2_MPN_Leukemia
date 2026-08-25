# Leukemia paper - GSE122198 metadata helper for Figure 2A
#
# Extract sample title, mutation, and diagnosis metadata from GEO.
# The accepted paper compares JAK2V617F-mutant and matched JAK2-WT CD34+ cells.
# This helper documents cohort metadata; the expression/plotting code used to
# generate the final Figure 2A was not present in the supplied working scripts.

suppressPackageStartupMessages(library(GEOquery))

gse_id <- "GSE122198"
gse <- getGEO(GEO = gse_id, GSEMatrix = FALSE)
gsm_list <- GSMList(gse)

sample_titles <- sapply(gsm_list, function(gsm) Meta(gsm)$title)

extract_characteristic <- function(gsm, key) {
  characs <- Meta(gsm)$characteristics_ch1
  hit <- characs[grepl(paste0("^", key, ":"), characs, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  sub(paste0("^", key, ":\\s*"), "", hit[1], ignore.case = TRUE)
}

metadata <- data.frame(
  GSM = names(gsm_list),
  sample_title = sample_titles,
  mutation = sapply(gsm_list, extract_characteristic, key = "mutation"),
  diagnosis = sapply(gsm_list, extract_characteristic, key = "diagnosis"),
  stringsAsFactors = FALSE
)

write.table(
  metadata,
  "GSE122198_mutation_diagnosis.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
