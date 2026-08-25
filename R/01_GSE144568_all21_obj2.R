# Leukemia paper - GSE144568: all 21 samples (obj2)
#
# Purpose:
#   Build the full 21-sample Seurat object used for Figure 1A and for the
#   downstream untreated-subset analysis in 02_GSE144568_Hoffman_obj3.R.
#
# Cohort in accepted manuscript:
#   healthy controls n = 6; MF n = 15
#   MF driver mutations: JAK2V617F n = 12; CALR n = 3

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

set.seed(89234)

# -----------------------------------------------------------------------------
# Input
# -----------------------------------------------------------------------------
# Expected directory structure:
#   scMatrices/<sample>/filtered_feature_bc_matrix/
# Adjust data_dir for your local copy of GSE144568 count matrices.

data_dir <- "scMatrices"
output_dir <- "results/Fig1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- list.files(data_dir, full.names = TRUE)
if (length(dirs) != 21) {
  warning("Expected 21 sample directories; found ", length(dirs), ".")
}

# -----------------------------------------------------------------------------
# Read and merge 21 samples
# -----------------------------------------------------------------------------
objs <- lapply(seq_along(dirs), function(i) {
  counts <- Read10X(file.path(dirs[i], "filtered_feature_bc_matrix"))
  CreateSeuratObject(
    counts = counts,
    project = "MPN_CD34",
    min.cells = 50,
    min.features = 700,
    max.features = 7000
  )
})

sample_ids <- paste0("Sample_", seq_along(objs))
obj2 <- merge(
  x = objs[[1]],
  y = objs[-1],
  add.cell.ids = sample_ids
)

# Recover numeric sample ID from merged cell names.
obj2$sample <- vapply(
  strsplit(colnames(obj2), "_"),
  function(x) x[2],
  character(1)
)
obj2$sample <- str_trim(as.character(obj2$sample))

# Sample-level annotation used in the paper.
cond <- c(
  "NM","MF","MF","MF","MF","NM","MF","MF","NM","MF","MF",
  "MF","NM","MF","MF","MF","NM","NM","MF","MF","MF"
)
names(cond) <- as.character(1:21)

cond2 <- c(
  "NM","JAK2","JAK2","CALR","JAK2","NM","JAK2","JAK2","NM",
  "JAK2","JAK2","JAK2","NM","JAK2","CALR","CALR","NM","NM",
  "JAK2","JAK2","JAK2"
)
names(cond2) <- as.character(1:21)

obj2$cond <- unname(cond[obj2$sample])
obj2$cond2 <- unname(cond2[obj2$sample])

stopifnot(!anyNA(obj2$cond), !anyNA(obj2$cond2))
print(table(obj2$cond))
print(table(obj2$cond2))

# -----------------------------------------------------------------------------
# Standard Seurat workflow
# -----------------------------------------------------------------------------
obj2 <- NormalizeData(obj2)
obj2 <- FindVariableFeatures(obj2, selection.method = "vst", nfeatures = 3000)
obj2 <- ScaleData(obj2, verbose = FALSE)
obj2 <- RunPCA(obj2, verbose = FALSE)
obj2 <- RunUMAP(obj2, dims = 1:20)

saveRDS(obj2, file = "obj2_GSE144568_all21.rds")

# -----------------------------------------------------------------------------
# Figure 1A: HMGA1/HMGA2 expression, all controls versus all MF samples
# -----------------------------------------------------------------------------
obj2$cond <- factor(obj2$cond, levels = c("NM", "MF"))
genes <- c("HMGA1", "HMGA2")
figure_cols <- c("NM" = "#00096B", "MF" = "purple")

make_violin <- function(gene) {
  plot_df <- FetchData(obj2, vars = c(gene, "cond"))
  names(plot_df)[1] <- "expr"

  mean_df <- plot_df %>%
    group_by(cond) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

  VlnPlot(
    obj2,
    features = gene,
    group.by = "cond",
    cols = figure_cols,
    pt.size = 0,
    raster = FALSE
  ) +
    geom_jitter(
      data = plot_df,
      aes(x = cond, y = expr),
      width = 0.15,
      size = 0.3,
      alpha = 0.2,
      color = "grey60",
      inherit.aes = FALSE
    ) +
    geom_point(
      data = mean_df,
      aes(x = cond, y = mean_expr),
      color = "black",
      size = 2,
      inherit.aes = FALSE
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0, 4))
}

fig1a <- lapply(genes, make_violin)
pdf(file.path(output_dir, "Fig1A_HMGA1_HMGA2_all21.pdf"), width = 6, height = 4)
for (p in fig1a) print(p)
dev.off()

# Summary values shown in the Figure 1A table.
fig1a_summary <- do.call(rbind, lapply(levels(obj2$cond), function(group) {
  cells <- WhichCells(obj2, expression = cond == group)
  x <- FetchData(obj2, vars = genes, cells = cells)
  data.frame(
    condition = group,
    total_cells = nrow(x),
    HMGA1_positive_n = sum(x$HMGA1 > 0),
    HMGA1_positive_pct = 100 * mean(x$HMGA1 > 0),
    HMGA1_mean = mean(x$HMGA1),
    HMGA2_positive_n = sum(x$HMGA2 > 0),
    HMGA2_positive_pct = 100 * mean(x$HMGA2 > 0),
    HMGA2_mean = mean(x$HMGA2)
  )
}))
write.csv(fig1a_summary, file.path(output_dir, "Fig1A_summary.csv"), row.names = FALSE)
print(fig1a_summary)
