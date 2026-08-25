# Independent verification of selected Figure 2D biomarker calculations.
#
# This script reproduces the manuscript spot-checks using the reported 2x2
# contingency-table counts. It is intentionally independent of the plotting code.

odds_pvalue <- function(tab) {
  odds_ratio <- tab[1, 1] * tab[2, 2] / (tab[1, 2] * tab[2, 1])
  p_value <- fisher.test(tab)$p.value
  c(odds_ratio = odds_ratio, fisher_p = p_value)
}

sensitivity_specificity <- function(tab) {
  sensitivity <- tab[1, 1] / (tab[1, 1] + tab[1, 2])
  specificity <- tab[2, 2] / (tab[2, 2] + tab[2, 1])
  c(sensitivity = sensitivity, specificity = specificity)
}

# HMGA1 high: MF vs ET
hmga1_mf_et <- rbind(
  MF = c(high = 36, low = 6),
  ET = c(high = 7, low = 17)
)
print(hmga1_mf_et)
print(odds_pvalue(hmga1_mf_et))

# HMGA1 high: MF vs PV
hmga1_mf_pv <- rbind(
  MF = c(high = 36, low = 6),
  PV = c(high = 17, low = 16)
)
print(hmga1_mf_pv)
print(odds_pvalue(hmga1_mf_pv))

# HMGA1 high: MF vs all non-MF
hmga1_mf_nonmf <- rbind(
  MF = c(high = 36, low = 6),
  nonMF = c(high = 24, low = 54)
)
print(hmga1_mf_nonmf)
print(odds_pvalue(hmga1_mf_nonmf))
print(sensitivity_specificity(hmga1_mf_nonmf))
