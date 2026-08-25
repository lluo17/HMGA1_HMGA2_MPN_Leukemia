# Leukemia paper - platelet HMGA1/HMGA2 biomarker analysis (Figure 2C-D)
#
# Input: data/platelet/synergy.xlsx
# Expected columns:
#   category : CTRL, ET, PV, MF
#   HMGA1    : platelet HMGA1 expression
#   HMGA2    : platelet HMGA2 expression
#
# The accepted manuscript dichotomized HMGA1/HMGA2 around the cohort median,
# displayed the four expression groups as a stacked bar plot (Fig. 2C), and
# compared MF with ET, PV, and all non-MF samples using Fisher's exact tests
# (Fig. 2D).

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

output_dir <- "results/Fig2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_excel("data/platelet/synergy.xlsx") %>%
  mutate(category = factor(category, levels = c("CTRL", "ET", "PV", "MF")))

stopifnot(all(c("category", "HMGA1", "HMGA2") %in% names(df)))

hmga1_median <- median(df$HMGA1, na.rm = TRUE)
hmga2_median <- median(df$HMGA2, na.rm = TRUE)

# The manuscript defines high expression as > median. In the supplied
# synergy.xlsx data, neither HMGA1 nor HMGA2 has a value exactly equal to the
# cohort median, so > and >= produce identical classifications.
df_q <- df %>%
  mutate(
    HMGA1_high = HMGA1 > hmga1_median,
    HMGA2_high = HMGA2 > hmga2_median,
    expression_group = case_when(
      !HMGA1_high & !HMGA2_high ~ "Low HMGA1/Low HMGA2",
      !HMGA1_high &  HMGA2_high ~ "Low HMGA1/High HMGA2",
       HMGA1_high & !HMGA2_high ~ "High HMGA1/Low HMGA2",
       HMGA1_high &  HMGA2_high ~ "High HMGA1/High HMGA2"
    ),
    expression_group = factor(
      expression_group,
      levels = c(
        "Low HMGA1/Low HMGA2",
        "Low HMGA1/High HMGA2",
        "High HMGA1/Low HMGA2",
        "High HMGA1/High HMGA2"
      )
    )
  )

# -----------------------------------------------------------------------------
# Figure 2C: four HMGA1/HMGA2 expression groups by diagnosis
# -----------------------------------------------------------------------------
stacked <- df_q %>%
  count(category, expression_group, name = "count") %>%
  group_by(category) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()

p2c <- ggplot(stacked, aes(x = category, y = percentage, fill = expression_group)) +
  geom_col(color = "white", width = 0.7) +
  geom_text(
    aes(label = percent(percentage, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      "Low HMGA1/Low HMGA2" = "#d7e1ee",
      "Low HMGA1/High HMGA2" = "#00FF30",
      "High HMGA1/Low HMGA2" = "#FF00FF",
      "High HMGA1/High HMGA2" = "#336b87"
    ),
    name = "Expression Group"
  ) +
  theme_classic() +
  labs(x = "Clinical Condition", y = "Percentage of Samples")

ggsave(file.path(output_dir, "Fig2C_HMGA1_HMGA2_expression_groups.pdf"), p2c, width = 7, height = 5)

# -----------------------------------------------------------------------------
# Figure 2D: biomarker comparison
# -----------------------------------------------------------------------------
# Definitions:
#   both_high = HMGA1 high AND HMGA2 high
#   HMGA1_high = HMGA1 high regardless of HMGA2
#   HMGA2_high = HMGA2 high regardless of HMGA1

df_q <- df_q %>% mutate(both_high = HMGA1_high & HMGA2_high)

comparison_sets <- list(
  `MF vs ET` = c("MF", "ET"),
  `MF vs PV` = c("MF", "PV"),
  `MF vs non-MF` = c("MF", "CTRL", "ET", "PV")
)

analyse_marker <- function(marker_col, comparison_name, categories) {
  x <- df_q %>% filter(category %in% categories)

  if (comparison_name == "MF vs non-MF") {
    x <- x %>% mutate(case = category == "MF")
  } else {
    comparator <- setdiff(categories, "MF")
    x <- x %>% mutate(case = category == "MF")
    stopifnot(length(comparator) == 1)
  }

  positive <- x[[marker_col]]
  tab <- table(case = x$case, positive = positive)
  ft <- fisher.test(tab)

  # Sensitivity/specificity are reported in the paper for MF vs all non-MF.
  sensitivity <- specificity <- NA_real_
  if (comparison_name == "MF vs non-MF") {
    sensitivity <- mean(positive[x$case], na.rm = TRUE)
    specificity <- mean(!positive[!x$case], na.rm = TRUE)
  }

  data.frame(
    marker = marker_col,
    comparison = comparison_name,
    odds_ratio = unname(ft$estimate),
    fisher_p = ft$p.value,
    sensitivity = sensitivity,
    specificity = specificity
  )
}

results <- bind_rows(lapply(c("both_high", "HMGA1_high", "HMGA2_high"), function(marker) {
  bind_rows(lapply(names(comparison_sets), function(comp) {
    analyse_marker(marker, comp, comparison_sets[[comp]])
  }))
}))

write.csv(results, file.path(output_dir, "Fig2D_biomarker_Fisher_tests.csv"), row.names = FALSE)
print(results)
