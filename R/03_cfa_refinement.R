# =============================================================================
# 03_cfa_refinement.R
# Ordinal CFA (WLSMV) with transparent, rule-based item refinement.
#
# Refinement rules (one item at a time, full re-estimation after each step;
# Costello & Osborne, 2005):
#   Phase 1  - remove the item with the lowest standardised loading < .50
#   Phase 2a - in subscales with AVE < .50, remove the lowest loading < .70
#   Phase 2b - if fit is still below the "good" thresholds, remove the item
#              with the largest modification indices (cross-loadings and
#              correlated residuals), i.e. the main local source of misfit
#   Guard    - a subscale never drops below 3 items
# Every step is logged (audit trail) and exported as an Mplus input file, so
# the same sequence can be re-run in Mplus.
# =============================================================================

suppressPackageStartupMessages({
  library(lavaan); library(semTools); library(psych); library(ggplot2)
})
dat <- read.csv("data/analysis_sample.csv")
codebook <- read.csv("data/codebook.csv")
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
dir.create("mplus", showWarnings = FALSE)

short <- c("Growth Mindset" = "GM", "Resilience" = "RES", "Stress Mindset" = "SM",
           "Self-Confidence" = "SC", "Performance Orientation" = "PO")
blocks <- split(codebook$item, short[codebook$subscale])[short]
names(blocks) <- short

thr <- list(cfi = .95, tli = .95, rmsea = .06, srmr = .08)  # "good" (Hu & Bentler, 1999)
good_fit <- function(f) f["cfi"] >= thr$cfi & f["tli"] >= thr$tli & f["rmsea"] <= thr$rmsea & f["srmr"] <= thr$srmr
# total modification index attached to each item (cross-loadings + residual covariances)
item_mi <- function(fit, items) {
  mi <- modindices(fit, sort. = FALSE)
  mi <- mi[mi$op %in% c("=~", "~~") & mi$lhs != mi$rhs, ]
  sapply(items, function(i) sum(mi$mi[(mi$op == "=~" & mi$rhs == i) |
                                      (mi$op == "~~" & (mi$lhs == i | mi$rhs == i))]))
}

model_syntax <- function(b) paste(sprintf("%s =~ %s", names(b), sapply(b, paste, collapse = " + ")),
                                  collapse = "\n")
fit_cfa <- function(b) cfa(model_syntax(b), data = dat, ordered = unlist(b),
                           estimator = "WLSMV", std.lv = TRUE)
get_fit <- function(fit) {
  f <- fitMeasures(fit, c("chisq.scaled", "df.scaled", "cfi.scaled", "tli.scaled",
                          "rmsea.scaled", "srmr", "wrmr"))
  c(chisq_df = unname(f[1] / f[2]), cfi = unname(f[3]), tli = unname(f[4]),
    rmsea = unname(f[5]), srmr = unname(f[6]), wrmr = unname(f[7]))
}
get_loadings <- function(fit) {
  s <- subset(standardizedSolution(fit), op == "=~")
  setNames(s$est.std, s$rhs)
}
vec <- function(x) {  # semTools >= 0.5-7 may return a data frame or a list
  if (is.data.frame(x)) x <- x[1, , drop = TRUE]
  if (is.list(x)) x <- sapply(x, as.numeric)
  x
}
factor_of <- function(item, b) names(b)[sapply(b, function(x) item %in% x)]

write_mplus <- function(b, file, title) {
  used <- unlist(b)
  wrap <- function(v) paste(strwrap(paste(v, collapse = " "), 70, prefix = "    "), collapse = "\n")
  txt <- c(sprintf("TITLE: %s - WLSMV", title), "",
           "DATA:", '  FILE = "../data/analysis_sample.dat";', "",
           "VARIABLE:", "  NAMES =", wrap(codebook$item), "  ;",
           "  USEVARIABLES =", wrap(used), "  ;",
           "  CATEGORICAL =", wrap(used), "  ;",
           "  MISSING = ALL (-999);", "",
           "ANALYSIS:", "  ESTIMATOR = WLSMV;", "  PARAMETERIZATION = DELTA;", "",
           "MODEL:", sprintf("  %s BY %s;", names(b), sapply(b, paste, collapse = " ")), "",
           "OUTPUT:", "  STDYX; MODINDICES(3.84);")
  writeLines(txt, file)
}
# Mplus-ready data file (no header, -999 for missing)
mp <- dat[, codebook$item]; mp[is.na(mp)] <- -999
write.table(mp, "data/analysis_sample.dat", row.names = FALSE, col.names = FALSE)

# --- Iterative refinement ----------------------------------------------------
b <- blocks; step <- 0; log <- list(); history <- list()
repeat {
  fit <- fit_cfa(b); f <- get_fit(fit); lam <- get_loadings(fit)
  ave <- vec(semTools::AVE(fit))[names(b)]
  write_mplus(b, sprintf("mplus/cfa_step_%02d.inp", step),
              if (step == 0) "Initial 5-factor model" else sprintf("Refinement step %d", step))
  history[[step + 1]] <- data.frame(step = step, t(round(f, 3)),
                                    ave_summary = paste(sprintf("%s=%.2f", names(ave), ave), collapse = "; "),
                                    n_items = length(unlist(b)))
  if (step == 0) { initial_fit <- fit; initial_b <- b }

  removable <- names(lam)[sapply(names(lam), function(i) length(b[[factor_of(i, b)]]) > 3)]
  lam_r <- lam[removable]
  low_ave <- names(ave)[ave < .50]

  cand <- NULL; reason <- NULL
  if (any(lam_r < .50)) {
    cand <- names(which.min(lam_r)); reason <- "Phase 1: loading < .50"
  } else {
    in_low <- lam_r[sapply(names(lam_r), factor_of, b = b) %in% low_ave & lam_r < .70]
    if (length(in_low)) {
      cand <- names(which.min(in_low)); reason <- sprintf("Phase 2a: loading < .70 in %s (AVE = %.2f)",
                                                         factor_of(cand, b), ave[factor_of(cand, b)])
    } else if (!good_fit(f) && length(lam_r)) {
      mis <- item_mi(fit, names(lam_r))
      cand <- names(which.max(mis))
      reason <- sprintf("Phase 2b: fit below 'good'; largest total MI (%.0f)", max(mis))
    }
  }
  if (is.null(cand)) break
  step <- step + 1
  log[[step]] <- data.frame(step = step, item_removed = cand, subscale = factor_of(cand, b),
                            loading = round(lam[cand], 3), reason = reason)
  b[[factor_of(cand, b)]] <- setdiff(b[[factor_of(cand, b)]], cand)
}
final_fit <- fit; final_b <- b

audit <- merge(do.call(rbind, history), rbind(data.frame(step = 0, item_removed = "-", subscale = "-",
                                                        loading = NA, reason = "Initial model"),
                                             do.call(rbind, log)), by = "step")
audit <- audit[, c("step", "item_removed", "subscale", "loading", "reason", "n_items",
                   "chisq_df", "cfi", "tli", "rmsea", "srmr", "wrmr", "ave_summary")]
write.csv(audit, "results/02_refinement_audit_trail.csv", row.names = FALSE)
cat("\nRefinement audit trail:\n")
print(audit[, c("step", "item_removed", "loading", "cfi", "tli", "rmsea", "srmr")], row.names = FALSE)

# --- Competing model: drop a subscale that cannot reach AVE >= .50 -----------
ave_final <- vec(semTools::AVE(final_fit))[names(final_b)]
weak <- names(ave_final)[ave_final < .50]
models <- list("5-factor initial" = initial_fit, "5-factor refined" = final_fit)
if (length(weak)) {
  b4 <- final_b[setdiff(names(final_b), weak)]
  fit4 <- fit_cfa(b4)
  models[[sprintf("4-factor refined (without %s)", paste(weak, collapse = ", "))]] <- fit4
  write_mplus(b4, "mplus/cfa_final_4factor.inp", "Final 4-factor model")
  final_fit <- fit4; final_b <- b4
}
comp <- do.call(rbind, lapply(names(models), function(m)
  data.frame(model = m, items = length(lavNames(models[[m]], "ov")), t(round(get_fit(models[[m]]), 3)))))
write.csv(comp, "results/03_model_comparison.csv", row.names = FALSE)
cat("\nModel comparison:\n"); print(comp, row.names = FALSE)

# --- Validity of the final model ---------------------------------------------
lam <- subset(standardizedSolution(final_fit), op == "=~")[, c("lhs", "rhs", "est.std", "se", "pvalue")]
names(lam) <- c("subscale", "item", "std_loading", "se", "p")
lam[, 3:5] <- lapply(lam[, 3:5], round, 3)
write.csv(lam, "results/04_final_loadings.csv", row.names = FALSE)

alpha_v <- sapply(final_b, function(v) psych::alpha(dat[, v], warnings = FALSE)$total$raw_alpha)
ord_alpha <- sapply(final_b, function(v) psych::alpha(psych::polychoric(dat[, v])$rho, warnings = FALSE)$total$raw_alpha)
omega_v <- vec(semTools::compRelSEM(final_fit))[names(final_b)]
ave_v <- vec(semTools::AVE(final_fit))[names(final_b)]
rel <- data.frame(subscale = names(final_b), items = lengths(final_b),
                  alpha = round(alpha_v, 3), ordinal_alpha = round(ord_alpha, 3),
                  omega = round(omega_v, 3), AVE = round(ave_v, 3))
write.csv(rel, "results/05_final_reliability_ave.csv", row.names = FALSE)

phi <- lavInspect(final_fit, "cor.lv")
fl <- round(phi, 3); diag(fl) <- round(sqrt(ave_v), 3); fl[upper.tri(fl)] <- NA
write.csv(cbind(subscale = rownames(fl), as.data.frame(fl)), "results/06_fornell_larcker.csv", row.names = FALSE)
htmt <- semTools::htmt(model_syntax(final_b), data = dat, ordered = unlist(final_b))
ht <- round(htmt, 3); ht[upper.tri(ht, diag = TRUE)] <- NA
write.csv(cbind(subscale = rownames(ht), as.data.frame(ht)), "results/07_htmt.csv", row.names = FALSE)

# Sensitivity: same final structure with robust ML
mlr <- cfa(model_syntax(final_b), data = dat, estimator = "MLR", missing = "fiml")
sens <- fitMeasures(mlr, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
write.csv(data.frame(index = names(sens), value = round(sens, 3)), "results/08_sensitivity_mlr.csv", row.names = FALSE)

cat("\nFinal reliability and AVE:\n"); print(rel, row.names = FALSE)
cat("\nMax HTMT:", max(ht, na.rm = TRUE), "\n")
cat("\nMLR sensitivity:\n"); print(round(sens, 3))

# --- Figures -----------------------------------------------------------------
h <- audit
long <- rbind(data.frame(step = h$step, index = "CFI", value = h$cfi),
              data.frame(step = h$step, index = "TLI", value = h$tli),
              data.frame(step = h$step, index = "RMSEA", value = h$rmsea),
              data.frame(step = h$step, index = "SRMR", value = h$srmr))
long$panel <- ifelse(long$index %in% c("CFI", "TLI"), "CFI / TLI (higher is better)", "RMSEA / SRMR (lower is better)")
cut_lines <- data.frame(panel = c("CFI / TLI (higher is better)", "RMSEA / SRMR (lower is better)"),
                        y = c(.95, .06))
g1 <- ggplot(long, aes(step, value, colour = index)) +
  geom_hline(data = cut_lines, aes(yintercept = y), linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = h$step, labels = ifelse(h$step == 0, "start", h$item_removed)) +
  scale_colour_manual(values = c(CFI = "#2b6cb0", TLI = "#63b3ed", RMSEA = "#dd6b20", SRMR = "#f6ad55")) +
  labs(title = "Fit indices across refinement steps (WLSMV)",
       subtitle = "Dashed lines: 'good fit' thresholds, CFI/TLI >= .95 and RMSEA <= .06",
       x = "Item removed at step", y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/fit_trajectory.png", g1, width = 8, height = 4.8, dpi = 110)

g2 <- ggplot(lam, aes(reorder(item, std_loading), std_loading, fill = subscale)) +
  geom_col(width = 0.7) + coord_flip() +
  geom_hline(yintercept = .70, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = sprintf("%.2f", std_loading)), hjust = 1.15, colour = "white", size = 3.4) +
  facet_grid(subscale ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("#2b6cb0", "#319795", "#dd6b20", "#805ad5", "#4a5568")) +
  labs(title = "Final model: standardised loadings", x = NULL, y = "Standardised loading") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave("figures/final_loadings.png", g2, width = 7, height = 6, dpi = 110)
