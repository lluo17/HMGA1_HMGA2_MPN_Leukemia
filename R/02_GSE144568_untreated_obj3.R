# Leukemia paper - GSE144568 untreated subset (obj3)
#
# Purpose:
#   Reproduce the untreated/age-matched subset used for Figure 1B-F.
#
# Accepted-manuscript subset:
#   healthy controls n = 4
#   untreated MF n = 6 (JAK2V617F n = 4; CALR n = 2)
#
# Run after 01_GSE144568_all21_obj2.R, or provide obj2_GSE144568_all21.rds.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

source("R/00_human_cell_markers.R")
set.seed(89234)
output_dir <- "results/Fig1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("obj2")) {
  obj2 <- readRDS("obj2_GSE144568_all21.rds")
}

# -----------------------------------------------------------------------------
# Define obj3: 4 controls + 6 untreated MF samples
# -----------------------------------------------------------------------------
obj3 <- subset(obj2, subset = sample %in% c(1, 6, 9, 13, 4, 5, 7, 8, 12, 15))
obj3$Cond2 <- ifelse(
  obj3$sample %in% c(1, 6, 9, 13), "NM",
  ifelse(obj3$sample %in% c(4, 15), "CALR", "JAK2")
)
obj3$cond <- ifelse(obj3$Cond2 == "NM", "NM", "MF")
obj3$Cond2 <- factor(obj3$Cond2, levels = c("NM", "JAK2", "CALR"))
obj3$cond <- factor(obj3$cond, levels = c("NM", "MF"))

print(table(obj3$Cond2))
print(table(obj3$cond))

# Re-run dimensional reduction and clustering on the subset, as used in paper.
obj3 <- NormalizeData(obj3)
obj3 <- FindVariableFeatures(obj3, selection.method = "vst", nfeatures = 3000)
obj3 <- ScaleData(obj3, verbose = FALSE)
obj3 <- RunPCA(obj3, verbose = FALSE)
obj3 <- RunUMAP(obj3, dims = 1:20)
obj3 <- FindNeighbors(obj3, dims = 1:20)
obj3 <- FindClusters(obj3, resolution = 0.25)
obj3$seurat_clusters <- Idents(obj3)

saveRDS(obj3, "obj3_GSE144568_untreated_subset.rds")

# -----------------------------------------------------------------------------
# Figure 1B: HMGA1/HMGA2 in 4 controls versus 6 untreated MF samples
# -----------------------------------------------------------------------------
genes <- c("HMGA1", "HMGA2")
figure_cols <- c("NM" = "#00096B", "MF" = "purple")

make_violin <- function(gene) {
  df <- FetchData(obj3, vars = c(gene, "cond"))
  names(df)[1] <- "expr"
  mean_df <- df %>%
    group_by(cond) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

  VlnPlot(
    obj3,
    features = gene,
    group.by = "cond",
    cols = figure_cols,
    pt.size = 0.001,
    raster = FALSE
  ) +
    geom_point(
      data = mean_df,
      aes(x = cond, y = mean_expr),
      color = "black",
      size = 2,
      inherit.aes = FALSE
    ) +
    theme_classic() +
    theme(axis.title.x = element_blank(), legend.position = "none") +
    coord_cartesian(ylim = c(0, 4))
}

fig1b <- lapply(genes, make_violin)
pdf(file.path(output_dir, "Fig1B_HMGA1_HMGA2_untreated_subset.pdf"), width = 6, height = 4)
for (p in fig1b) print(p)
dev.off()

fig1b_pvalues <- data.frame(
  gene = genes,
  wilcoxon_p = vapply(genes, function(g) {
    x <- FetchData(obj3, vars = c(g, "cond"))
    names(x)[1] <- "expr"
    wilcox.test(expr ~ cond, data = x)$p.value
  }, numeric(1))
)
write.csv(fig1b_pvalues, file.path(output_dir, "Fig1B_wilcoxon_pvalues.csv"), row.names = FALSE)

fig1b_summary <- do.call(rbind, lapply(levels(obj3$cond), function(group) {
  cells <- WhichCells(obj3, expression = cond == group)
  x <- FetchData(obj3, vars = genes, cells = cells)
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
write.csv(fig1b_summary, file.path(output_dir, "Fig1B_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Figure 1C: UMAP by driver group and by Seurat cluster
# -----------------------------------------------------------------------------
condition_cols <- c("NM" = "#00008B", "JAK2" = "red", "CALR" = "green")

p_cond <- UMAPPlot(obj3, group.by = "Cond2", cols = condition_cols, raster = FALSE) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

p_cluster <- UMAPPlot(obj3, group.by = "seurat_clusters", label = TRUE, raster = FALSE)

pdf(file.path(output_dir, "Fig1C_UMAP_condition_and_clusters.pdf"), width = 10, height = 5)
print(p_cond + p_cluster)
dev.off()

# Marker plot used to support manual cluster annotation.
p_marker <- DotPlot(obj3, features = cellGW) + RotatedAxis()
pdf(file.path(output_dir, "Fig1C_human_cell_markers.pdf"), width = 16, height = 7)
print(p_marker)
dev.off()

# -----------------------------------------------------------------------------
# Figure 1D: HMGA1/HMGA2 feature plots using the same real expression scale
# -----------------------------------------------------------------------------
expr_df <- FetchData(obj3, vars = genes)
global_max <- max(expr_df, na.rm = TRUE)

feature_plots <- FeaturePlot(
  obj3,
  features = genes,
  pt.size = 0.1,
  order = TRUE,
  raster = FALSE,
  combine = FALSE
)
feature_plots <- lapply(feature_plots, function(p) {
  p + scale_color_gradientn(
    colours = c("lightgrey", "#6BAED6", "#08519C"),
    limits = c(0, global_max),
    oob = scales::squish
  )
})

pdf(file.path(output_dir, "Fig1D_HMGA1_HMGA2_featureplots.pdf"), width = 10, height = 5)
print(wrap_plots(feature_plots, ncol = 2, guides = "collect"))
dev.off()

# -----------------------------------------------------------------------------
# Figure 1E: HMGA1/HMGA2 dot plot with shared true expression scale
# -----------------------------------------------------------------------------
expr <- FetchData(obj3, vars = c(genes, "seurat_clusters"))
names(expr)[ncol(expr)] <- "Cluster"

plot_df <- expr %>%
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expr") %>%
  group_by(Cluster, Gene) %>%
  summarise(
    AvgExpr = mean(Expr, na.rm = TRUE),
    PctExpr = 100 * mean(Expr > 0, na.rm = TRUE),
    .groups = "drop"
  )

p_dot <- ggplot(plot_df, aes(x = Gene, y = Cluster)) +
  geom_point(aes(color = AvgExpr, size = PctExpr)) +
  scale_color_gradientn(
    colours = c("lightgrey", "#6BAED6", "#08519C"),
    limits = range(plot_df$AvgExpr, na.rm = TRUE),
    oob = scales::squish,
    name = "Average Expression\n(Log-Normalized)"
  ) +
  scale_size_continuous(name = "% Cells expressing") +
  theme_classic() +
  labs(x = NULL, y = NULL)

pdf(file.path(output_dir, "Fig1E_HMGA1_HMGA2_dotplot.pdf"), width = 4, height = 6)
print(p_dot)
dev.off()

# -----------------------------------------------------------------------------
# Figure 1F: paired HMGA1/HMGA2 violin distributions by cluster
# -----------------------------------------------------------------------------
violin_df <- FetchData(obj3, vars = c(genes, "seurat_clusters")) %>%
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expression")

p_violin <- ggplot(violin_df, aes(x = seurat_clusters, y = Expression, fill = Gene)) +
  geom_violin(position = position_dodge(width = 0.9), scale = "width", alpha = 0.7) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9),
    size = 0.25,
    alpha = 0.4,
    color = "grey40"
  ) +
  stat_summary(
    fun.data = function(x) data.frame(
      y = median(x),
      ymin = quantile(x, 0.25),
      ymax = quantile(x, 0.75)
    ),
    geom = "errorbar",
    width = 0.3,
    position = position_dodge(width = 0.9)
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.9),
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(values = c("HMGA1" = "magenta", "HMGA2" = "#00FF30")) +
  theme_classic() +
  labs(x = NULL, y = "Normalized Expression")

pdf(file.path(output_dir, "Fig1F_HMGA1_HMGA2_cluster_violins.pdf"), width = 14, height = 6)
print(p_violin)
dev.off()
