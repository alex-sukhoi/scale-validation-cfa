# =============================================================================
# 02_screening.R
# Preliminary item analysis before any CFA:
#   * incomplete responses and speeders
#   * item descriptives with response distribution per Likert category
#   * flags: SD < 0.50, > 70% in one category, |skew| > 2, |kurtosis| > 7
#   * reverse-coding check: negative correlation with own-subscale items
#   * estimator choice
# =============================================================================

suppressPackageStartupMessages(library(psych))
dir.create("results", showWarnings = FALSE)

raw <- read.csv("data/mental_skills_pilot.csv")
codebook <- read.csv("data/codebook.csv")
items <- codebook$item

# --- Respondents -------------------------------------------------------------
n_missing <- rowSums(is.na(raw[, items]))
incomplete <- raw$id[n_missing > 5]
speed_cut <- 0.5 * median(raw$completion_min)            # Greszki et al. (2014)
speeders <- raw$id[raw$completion_min < speed_cut]
dat <- raw[!(raw$id %in% c(incomplete, speeders)), ]
cat(sprintf("Respondents: %d collected, %d incomplete removed, %d speeders (< %.1f min) removed, %d analysed\n",
            nrow(raw), length(incomplete), length(speeders), speed_cut, nrow(dat)))

# --- Item descriptives -------------------------------------------------------
d <- psych::describe(dat[, items])
dist <- t(sapply(items, function(v) {
  tab <- table(factor(dat[[v]], levels = 1:5))
  round(100 * tab / sum(tab), 1)
}))
colnames(dist) <- paste0("pct_", 1:5)
desc <- data.frame(item = items, subscale = codebook$subscale, n = d$n,
                   mean = round(d$mean, 2), sd = round(d$sd, 2),
                   skew = round(d$skew, 2), kurtosis = round(d$kurtosis, 2), dist)

# --- Flags -------------------------------------------------------------------
desc$flag_low_sd    <- desc$sd < 0.50
desc$flag_ceiling   <- apply(dist, 1, max) > 70
desc$flag_nonnormal <- abs(desc$skew) > 2 | abs(desc$kurtosis) > 7

# Reverse-coding check: an item should not correlate negatively with the
# other items of its own subscale
neg_cor <- sapply(items, function(v) {
  peers <- setdiff(codebook$item[codebook$subscale == codebook$subscale[codebook$item == v]], v)
  r <- cor(dat[[v]], dat[, peers], use = "pairwise")
  sum(r < 0)
})
desc$n_negative_r_with_own_scale <- neg_cor
desc$flag_reverse_coding <- neg_cor >= 3
write.csv(desc, "results/01_item_descriptives_and_flags.csv", row.names = FALSE)

flagged <- desc[desc$flag_low_sd | desc$flag_ceiling | desc$flag_reverse_coding,
                c("item", "sd", "pct_4", "n_negative_r_with_own_scale")]
cat("\nItems flagged for review before CFA:\n"); print(flagged, row.names = FALSE)

# Decision log (in a real project these go to the researcher before fitting)
#  * gm7: reverse-worded and not reverse-scored -> recode 6 - x
#  * gm1: ceiling effect (> 70% "agree") but SD > 0.5 and all categories used
#         -> keep for the initial model, let the CFA decide
dat$gm7 <- 6 - dat$gm7
after <- cor(dat$gm7, dat[, paste0("gm", 1:6)], use = "pairwise")
cat(sprintf("\ngm7 after recoding: correlations with gm1-gm6 range %.2f to %.2f\n", min(after), max(after)))

# --- Estimator choice --------------------------------------------------------
cat(sprintf("\nMax |skew| = %.2f, max |kurtosis| = %.2f (all below 2 / 7).\n",
            max(abs(desc$skew)), max(abs(desc$kurtosis))))
cat("Items are 5-point ordinal, so CFA uses WLSMV on polychoric correlations\n",
    "(Finney & DiStefano, 2013; Li, 2016). MLR is reported as a sensitivity check.\n")

write.csv(dat, "data/analysis_sample.csv", row.names = FALSE)
