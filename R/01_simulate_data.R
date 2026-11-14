# =============================================================================
# 01_simulate_data.R
# Synthetic data for a new 35-item "Workplace Mental Skills" questionnaire:
# 5 hypothesised subscales x 7 items, 5-point Likert, n = 300 professionals.
#
# The population model deliberately contains the problems a real pilot
# usually has, so the validation pipeline has something to find:
#   * some weak items (loadings 0.30-0.50)
#   * one subscale (Growth Mindset) that is weakly defined overall
#   * one reverse-worded item that was NOT reverse-scored
#   * a ceiling effect (one item with > 70% in a single category)
#   * one respondent who stopped half-way
# =============================================================================

suppressPackageStartupMessages(library(lavaan))
set.seed(4242)
RNGkind(sample.kind = "Rejection")

load <- list(
  GM  = c(0.62, 0.55, 0.48, 0.58, 0.40, 0.52, 0.33),   # Growth Mindset (weak)
  RES = c(0.78, 0.80, 0.74, 0.76, 0.52, 0.49, 0.58),   # Resilience
  SM  = c(0.80, 0.74, 0.80, 0.55, 0.74, 0.47, 0.78),   # Stress Mindset
  SC  = c(0.79, 0.75, 0.51, 0.74, 0.46, 0.55, 0.77),   # Self-Confidence
  PO  = c(0.48, 0.77, 0.75, 0.80, 0.57, 0.50, 0.74)    # Performance Orientation
)
prefix <- c(GM = "gm", RES = "res", SM = "sm", SC = "sc", PO = "po")

meas <- unlist(lapply(names(load), function(f) {
  sprintf("%s =~ %s", f, paste0(load[[f]], "*", prefix[f], 1:7, collapse = " + "))
}))
latent_cor <- '
  GM ~~ 0.35*RES + 0.30*SM + 0.30*SC + 0.15*PO
  RES ~~ 0.40*SM + 0.50*SC + 0.20*PO
  SM ~~ 0.40*SC + 0.25*PO
  SC ~~ 0.35*PO
'
# Real questionnaires are never perfectly simple-structured: some items
# cross-load and some share wording (nuisance "doublet" factors)
misfit <- '
  SM  =~ 0.50*res5 + 0.40*res7
  RES =~ 0.50*sc5 + 0.40*sm6
  SC  =~ 0.45*po6 + 0.40*gm3 + 0.40*sm4
  PO  =~ 0.40*sc3
  W1  =~ 0.55*sm2 + 0.55*sm3
  W2  =~ 0.55*res1 + 0.55*res2
  W3  =~ 0.50*po2 + 0.50*po3
  W1 ~~ 0*GM + 0*RES + 0*SM + 0*SC + 0*PO + 0*W2 + 0*W3
  W2 ~~ 0*GM + 0*RES + 0*SM + 0*SC + 0*PO + 0*W3
  W3 ~~ 0*GM + 0*RES + 0*SM + 0*SC + 0*PO
'
pop_model <- paste(c(meas, latent_cor, misfit), collapse = "\n")

n <- 300
cont <- simulateData(pop_model, sample.nobs = n, standardized = TRUE)

# Agreement-skewed 5-point thresholds (typical for self-report skill scales)
to_likert <- function(x, shift = 0) cut(x + shift, c(-Inf, -2.0, -1.1, -0.2, 1.0, Inf), labels = FALSE)
dat <- as.data.frame(lapply(cont, to_likert))
dat$gm1 <- to_likert(cont$gm1, shift = 0.75)       # ceiling effect: most answer "4"
dat$gm1[dat$gm1 == 5 & runif(n) < 0.7] <- 4

# gm7 is reverse-worded ("Talent is something you either have or you don't")
# and was left un-reversed in the export
dat$gm7 <- 6 - dat$gm7

dat <- dat[, unlist(lapply(names(prefix), function(f) paste0(prefix[f], 1:7)))]

# One respondent abandoned the survey after item 20
dat[137, 21:35] <- NA

dat <- cbind(id = sprintf("P%03d", 1:n),
             completion_min = round(rlnorm(n, log(10), 0.35), 1),
             dat)
dir.create("data", showWarnings = FALSE)
write.csv(dat, "data/mental_skills_pilot.csv", row.names = FALSE)

codebook <- data.frame(
  item = names(dat)[-(1:2)],
  subscale = rep(c("Growth Mindset", "Resilience", "Stress Mindset",
                   "Self-Confidence", "Performance Orientation"), each = 7),
  reverse_worded = names(dat)[-(1:2)] %in% c("gm7")
)
write.csv(codebook, "data/codebook.csv", row.names = FALSE)
cat("Saved data/mental_skills_pilot.csv (", n, "respondents, 35 items )\n")
