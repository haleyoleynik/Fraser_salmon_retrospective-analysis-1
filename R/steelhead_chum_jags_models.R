# ============================================================
# Stock-recruit JAGS models: Chilcotin steelhead, Thompson steelhead,
# and Fraser chum -- all in one script, per the updated dredge results
# (Fraser_steelhead_dredge-results_v2.csv, Fraser_chum_dredge-results.csv).
#
# COVARIATE CHANGES vs. the previous versions:
#   Chilcotin steelhead: DROPS max_flow entirely. The updated dredge
#     top model is Spawners + SL + SST only -- identical in structure
#     to Thompson now (previously Chilcotin had max_flow, Thompson didn't).
#     steelhead_data.csv no longer even has a max_flow column, consistent
#     with this.
#   Thompson steelhead: unchanged (Spawners + SL + SST).
#   Chum: NEW model, not previously implemented as JAGS. The original
#     WinBUGS model (model_cow_m1_2-v4.odc) used SSL + PDOA + PDOJ + NPGO;
#     the updated dredge top model keeps only PDO + SL (renamed from
#     PDOA/SSL to match this pipeline's naming convention elsewhere,
#     e.g. coho/sockeye's "SL"), dropping PDOJ and NPGO. chum_data.csv's
#     columns (Year, spawners, recruits, lnrs, PDO, SL) already match
#     this exactly.
#
# WinBUGS -> JAGS mapping notes (same as the steelhead scripts):
#   check/data/compile(2)  -> jags.model()
#   inits(1,...)/inits(2,...) -> inits list passed to jags.model()
#   update(2000)            -> update(jm, 2000)                 [burn-in]
#   dic.set()/dic.stats()   -> dic.samples(jm, ...)              [post-hoc DIC]
#   update(20000/30000)     -> coda.samples(jm, ..., n.iter=...)
#   beg(501)/beg(2001)      -> window(samples, start=...)
#   set(<param>)            -> variable.names in coda.samples()
#
# Requires: JAGS itself installed (https://mcmc-jags.sourceforge.io/)
# plus the R packages below.
# ============================================================

library(rjags)
library(coda)
library(dplyr)
library(readr)
library(ggplot2)

set.seed(1)

# ============================================================
# SHARED HELPER -- extracts year-indexed reference-point series
# (ayst, smsyst, rmsyst, msyst, umsyst) from a coda samples object,
# used identically for all three models below
# ============================================================

extract_year_series <- function(combined, prefix, years) {
  cols <- grep(paste0("^", prefix, "\\["), colnames(combined), value = TRUE)
  cols <- cols[order(as.numeric(gsub(paste0(prefix, "\\[|\\]"), "", cols)))]
  data.frame(
    Year   = years,
    metric = prefix,
    mean   = colMeans(combined[, cols, drop = FALSE]),
    q2.5   = apply(combined[, cols, drop = FALSE], 2, quantile, 0.025),
    q97.5  = apply(combined[, cols, drop = FALSE], 2, quantile, 0.975)
  )
}

# ============================================================
# 1. CHILCOTIN STEELHEAD -- Spawners, SL, SST (max_flow DROPPED)
# ============================================================

steelhead_data <- read_csv("data/steelhead_data.csv")

chilcotin_data <- steelhead_data %>%
  filter(Stock == "Chilcotin", !is.na(Spawners), !is.na(lnRS), !is.na(SST), !is.na(SL)) %>%
  arrange(Year)

stopifnot(nrow(chilcotin_data) == 46)   # matches steelhead_data.csv's current Chilcotin row count

chilcotin_model_string <- "
model{
  intercept ~ dnorm(0, 0.000001)
  LSE       ~ dunif(LminSE, LmaxSE)
  LminSE    <- log(0.01)
  LmaxSE    <- log(10)
  SE        <- exp(LSE)
  tau       <- 1 / (SE * SE)

  # Bounded Ricker b prior (Skeena methodology, same as sockeye): lognormal
  # centered on the frequentist MLE (median = 1/prSmax), truncated below at
  # 1/Smaxmax so Smax can't exceed 5x the MLE's own Smax. Replaces the
  # original WinBUGS slope ~ dunif(0.0001, 1000) prior.
  ms    <- log(1 / prSmax)
  msmin <- 1 / Smaxmax
  taus  <- 1 / (prCV * prCV)
  b     ~ dlnorm(ms, taus) T(msmin, )

  s ~ dnorm(0, 0.0001)   # SL coefficient
  t ~ dnorm(0, 0.0001)   # SST coefficient

  a <- intercept

  PPslope <- step(intercept)
  pvalues <- step(-s)
  pvaluet <- step(-t)

  smax <- 1 / b

  for (i in 1:ndata) {
    lnrs_pred[i] <- a - b * sp[i] + s * SL[i] + t * SST[i]
    lnrs[i]      ~ dnorm(lnrs_pred[i], tau)
    lnrs_rep[i]  ~ dnorm(lnrs_pred[i], tau)
    pvalue[i]    <- step(lnrs_rep[i] - lnrs[i])

    ayst[i]    <- a + s * SL[i] + t * SST[i]
    cstx[i]    <- ayst[i] / b
    cst[i]     <- max(cstx[i], 0.0001)
    smsystx[i] <- cst[i] * (0.5 - 0.07 * ayst[i])
    smsyst[i]  <- max(smsystx[i], 0.0001)
    rmsyst[i]  <- smsyst[i] * exp(ayst[i] - b * smsyst[i])
    msystx[i]  <- rmsyst[i] - smsyst[i]
    msyst[i]   <- max(msystx[i], 0.00000001)
    umsystx[i] <- msyst[i] / rmsyst[i]
    umsysty[i] <- max(umsystx[i], 0.0001)
    umsyst[i]  <- min(umsysty[i], 0.999)
  }
}
"

chilcotin_ricker_priors <- list(prSmax = 0.8950, Smaxmax = 4.4750, prCV = 2)

chilcotin_jags_data <- list(
  ndata = nrow(chilcotin_data), sp = chilcotin_data$Spawners, lnrs = chilcotin_data$lnRS,
  SL = chilcotin_data$SL, SST = chilcotin_data$SST,
  prSmax = chilcotin_ricker_priors$prSmax, Smaxmax = chilcotin_ricker_priors$Smaxmax,
  prCV = chilcotin_ricker_priors$prCV
)

# b starting values: the prior median (1/prSmax, i.e. the MLE itself) for
# chain 1, a nearby distinct value for chain 2 -- both comfortably above
# msmin since Smaxmax is 5x Smax_mle
chilcotin_b_init_1 <- 1 / chilcotin_ricker_priors$prSmax
chilcotin_b_init_2 <- 1.3 / chilcotin_ricker_priors$prSmax

chilcotin_inits <- list(
  list(LSE = 0.9, intercept = 1.1, s = -0.75, t = -0.127, b = chilcotin_b_init_1),
  list(LSE = 0.8, intercept = 1.2, s = -0.61, t = -0.137, b = chilcotin_b_init_2)
)

chilcotin_params <- c("intercept", "a", "b", "s", "t", "SE",
                      "PPslope", "pvalue", "pvalues", "pvaluet",
                      "ayst", "cst", "smax", "smsyst", "rmsyst", "msyst", "umsyst", "deviance")

chilcotin_jm <- jags.model(textConnection(chilcotin_model_string), data = chilcotin_jags_data,
                           inits = chilcotin_inits, n.chains = 2, n.adapt = 1000)
update(chilcotin_jm, n.iter = 2000)
chilcotin_dic <- dic.samples(chilcotin_jm, n.iter = 2000, type = "pD")
chilcotin_samples <- coda.samples(chilcotin_jm, variable.names = chilcotin_params, n.iter = 30000, thin = 1)
chilcotin_samples_kept <- window(chilcotin_samples, start = 2001)

chilcotin_scalar_params <- c("a", "b", "s", "t", "SE", "smax", "deviance")
cat("\n==== Chilcotin ====\n")
print(gelman.diag(chilcotin_samples_kept[, chilcotin_scalar_params], multivariate = FALSE))
print(effectiveSize(chilcotin_samples_kept[, chilcotin_scalar_params]))
print(chilcotin_dic)

# ============================================================
# 2. THOMPSON STEELHEAD -- Spawners, SL, SST (unchanged)
# ============================================================

thompson_data <- steelhead_data %>%
  filter(Stock == "Thompson", !is.na(Spawners), !is.na(lnRS), !is.na(SST), !is.na(SL)) %>%
  arrange(Year)

stopifnot(nrow(thompson_data) == 46)

thompson_model_string <- "
model{
  intercept ~ dnorm(0, 0.000001)
  LSE       ~ dunif(LminSE, LmaxSE)
  LminSE    <- log(0.01)
  LmaxSE    <- log(10)
  SE        <- exp(LSE)
  tau       <- 1 / (SE * SE)

  # Bounded Ricker b prior -- same methodology as Chilcotin/sockeye
  ms    <- log(1 / prSmax)
  msmin <- 1 / Smaxmax
  taus  <- 1 / (prCV * prCV)
  b     ~ dlnorm(ms, taus) T(msmin, )

  s ~ dnorm(0, 0.0001)
  t ~ dnorm(0, 0.0001)

  a <- intercept

  PPslope <- step(intercept)
  pvalues <- step(-s)
  pvaluet <- step(-t)

  smax <- 1 / b

  for (i in 1:ndata) {
    lnrs_pred[i] <- a - b * sp[i] + s * SL[i] + t * SST[i]
    lnrs[i]      ~ dnorm(lnrs_pred[i], tau)
    lnrs_rep[i]  ~ dnorm(lnrs_pred[i], tau)
    pvalue[i]    <- step(lnrs_rep[i] - lnrs[i])

    ayst[i]    <- a + s * SL[i] + t * SST[i]
    cstx[i]    <- ayst[i] / b
    cst[i]     <- max(cstx[i], 0.0001)
    smsystx[i] <- cst[i] * (0.5 - 0.07 * ayst[i])
    smsyst[i]  <- max(smsystx[i], 0.0001)
    rmsyst[i]  <- smsyst[i] * exp(ayst[i] - b * smsyst[i])
    msystx[i]  <- rmsyst[i] - smsyst[i]
    msyst[i]   <- max(msystx[i], 0.00000001)
    umsystx[i] <- msyst[i] / rmsyst[i]
    umsysty[i] <- max(umsystx[i], 0.0001)
    umsyst[i]  <- min(umsysty[i], 0.999)
  }
}
"

thompson_ricker_priors <- list(prSmax = 1.2567, Smaxmax = 6.2836, prCV = 2)

thompson_jags_data <- list(
  ndata = nrow(thompson_data), sp = thompson_data$Spawners, lnrs = thompson_data$lnRS,
  SL = thompson_data$SL, SST = thompson_data$SST,
  prSmax = thompson_ricker_priors$prSmax, Smaxmax = thompson_ricker_priors$Smaxmax,
  prCV = thompson_ricker_priors$prCV
)

thompson_b_init_1 <- 1 / thompson_ricker_priors$prSmax
thompson_b_init_2 <- 1.3 / thompson_ricker_priors$prSmax

thompson_inits <- list(
  list(LSE = 0.9, intercept = 1.1, s = -0.75, t = -0.127, b = thompson_b_init_1),
  list(LSE = 0.8, intercept = 1.2, s = -0.61, t = -0.137, b = thompson_b_init_2)
)

thompson_params <- c("intercept", "a", "b", "s", "t", "SE",
                     "PPslope", "pvalue", "pvalues", "pvaluet",
                     "ayst", "cst", "smax", "smsyst", "rmsyst", "msyst", "umsyst", "deviance")

thompson_jm <- jags.model(textConnection(thompson_model_string), data = thompson_jags_data,
                          inits = thompson_inits, n.chains = 2, n.adapt = 1000)
update(thompson_jm, n.iter = 2000)
thompson_dic <- dic.samples(thompson_jm, n.iter = 2000, type = "pD")
thompson_samples <- coda.samples(thompson_jm, variable.names = thompson_params, n.iter = 30000, thin = 1)
thompson_samples_kept <- window(thompson_samples, start = 2001)

thompson_scalar_params <- c("a", "b", "s", "t", "SE", "smax", "deviance")
cat("\n==== Thompson ====\n")
print(gelman.diag(thompson_samples_kept[, thompson_scalar_params], multivariate = FALSE))
print(effectiveSize(thompson_samples_kept[, thompson_scalar_params]))
print(thompson_dic)

# ============================================================
# 3. CHUM -- NEW model. Spawners, PDO, SL 
# ============================================================

chum_data_raw <- read_csv("data/chum_data.csv")

chum_data <- chum_data_raw %>%
  filter(!is.na(spawners), !is.na(lnrs), !is.na(PDO), !is.na(SL)) %>%
  arrange(Year)

stopifnot(nrow(chum_data) == 66)   # matches the original WinBUGS ndata=66

# SCALING FIX: the original WinBUGS data (data_m1_1b-v4.odc) has sp[]
# values like 1.7173, 1.7625, 1.3763 -- NOT raw spawner counts. Checked
# directly against chum_data.csv: 171725 / 100000 = 1.71725, an exact
# match. The model's slope prior (dunif(0.00001, 1000), "max smax of
# 100,000" per the original comment) was calibrated entirely around
# spawners scaled this way. Feeding it raw spawners (up to ~2,000,000)
# pins b against its lower bound and forces the intercept to compensate,
# which is exactly what was producing ayst drifting up to ~20 -- not a
# coincidental fitting issue, a real unit mismatch.
chum_spawners_scaled <- chum_data$spawners / 100000

chum_model_string <- "
model{
  intercept ~ dnorm(0, 0.000001)
  LSE       ~ dunif(LminSE, LmaxSE)
  LminSE    <- log(0.01)
  LmaxSE    <- log(10)
  SE        <- exp(LSE)
  tau       <- 1 / (SE * SE)

  # Bounded Ricker b prior -- same methodology as Chilcotin/Thompson/
  # sockeye. prSmax/Smaxmax here are on the SAME scaled-spawners unit as
  # `sp` below (spawners / 100000), not raw spawner counts.
  ms    <- log(1 / prSmax)
  msmin <- 1 / Smaxmax
  taus  <- 1 / (prCV * prCV)
  b     ~ dlnorm(ms, taus) T(msmin, )

  s ~ dnorm(0, 0.0001)   # SL coefficient
  t ~ dnorm(0, 0.0001)   # PDO coefficient

  a <- intercept

  PPslope <- step(intercept)
  pvalues <- step(-s)
  pvaluet <- step(-t)

  smax <- 1 / b

  for (i in 1:ndata) {
    lnrs_pred[i] <- a - b * sp[i] + s * SL[i] + t * PDO[i]
    lnrs[i]      ~ dnorm(lnrs_pred[i], tau)
    lnrs_rep[i]  ~ dnorm(lnrs_pred[i], tau)
    pvalue[i]    <- step(lnrs_rep[i] - lnrs[i])

    ayst[i]    <- a + s * SL[i] + t * PDO[i]
    cstx[i]    <- ayst[i] / b
    cst[i]     <- max(cstx[i], 0.0001)
    smsystx[i] <- cst[i] * (0.5 - 0.07 * ayst[i])
    smsyst[i]  <- max(smsystx[i], 0.0001)
    rmsyst[i]  <- smsyst[i] * exp(ayst[i] - b * smsyst[i])
    msystx[i]  <- rmsyst[i] - smsyst[i]
    msyst[i]   <- max(msystx[i], 0.00000001)
    umsystx[i] <- msyst[i] / rmsyst[i]
    umsysty[i] <- max(umsystx[i], 0.0001)
    umsyst[i]  <- min(umsysty[i], 0.999)
  }
}
"

# Computed from the dredge Spawners coefficient (RAW-spawner scale,
# -4.359408204003911e-7), converted to the /100000-scaled unit used for
# `sp` above: b_mle_scaled = -coef * 100000 = 0.0436. Smax_mle_scaled
# implies ~2.29 million RAW spawners -- checked sensible against chum's
# actual historical spawner range (up to ~2 million).
chum_ricker_priors <- list(prSmax = 22.9389, Smaxmax = 114.6945, prCV = 2)

chum_jags_data <- list(
  ndata = nrow(chum_data), sp = chum_spawners_scaled, lnrs = chum_data$lnrs,
  SL = chum_data$SL, PDO = chum_data$PDO,
  prSmax = chum_ricker_priors$prSmax, Smaxmax = chum_ricker_priors$Smaxmax,
  prCV = chum_ricker_priors$prCV
)

chum_b_init_1 <- 1 / chum_ricker_priors$prSmax
chum_b_init_2 <- 1.3 / chum_ricker_priors$prSmax

chum_inits <- list(
  list(LSE = 0.9, intercept = 1.1, s = -0.3, t = -0.1, b = chum_b_init_1),
  list(LSE = 0.8, intercept = 1.3, s = -0.2, t = -0.2, b = chum_b_init_2)
)

chum_params <- c("intercept", "a", "b", "s", "t", "SE",
                 "PPslope", "pvalue", "pvalues", "pvaluet",
                 "ayst", "cst", "smax", "smsyst", "rmsyst", "msyst", "umsyst", "deviance")

chum_jm <- jags.model(textConnection(chum_model_string), data = chum_jags_data,
                      inits = chum_inits, n.chains = 2, n.adapt = 1000)
update(chum_jm, n.iter = 2000)
chum_dic <- dic.samples(chum_jm, n.iter = 2000, type = "pD")
chum_samples <- coda.samples(chum_jm, variable.names = chum_params, n.iter = 30000, thin = 1)
chum_samples_kept <- window(chum_samples, start = 2001)

chum_scalar_params <- c("a", "b", "s", "t", "SE", "smax", "deviance")
cat("\n==== Chum ====\n")
print(gelman.diag(chum_samples_kept[, chum_scalar_params], multivariate = FALSE))
print(effectiveSize(chum_samples_kept[, chum_scalar_params]))
print(chum_dic)

# ============================================================
# 4. SAVE POSTERIOR SAMPLES -- one RDS per stock
# ============================================================

saveRDS(chilcotin_samples_kept, "chilcotin_posterior_samples.rds")
saveRDS(thompson_samples_kept,  "thompson_posterior_samples.rds")
saveRDS(chum_samples_kept,      "chum_posterior_samples.rds")

# ============================================================
# 5. POSTERIOR SUMMARY TABLE -- all three models together
# ============================================================

posterior_summary <- bind_rows(
  { s <- summary(chilcotin_samples_kept[, chilcotin_scalar_params])
  tibble(Stock = "Chilcotin", parameter = rownames(s$statistics), mean = s$statistics[, "Mean"],
         sd = s$statistics[, "SD"], q2.5 = s$quantiles[, "2.5%"], median = s$quantiles[, "50%"],
         q97.5 = s$quantiles[, "97.5%"]) },
  { s <- summary(thompson_samples_kept[, thompson_scalar_params])
  tibble(Stock = "Thompson", parameter = rownames(s$statistics), mean = s$statistics[, "Mean"],
         sd = s$statistics[, "SD"], q2.5 = s$quantiles[, "2.5%"], median = s$quantiles[, "50%"],
         q97.5 = s$quantiles[, "97.5%"]) },
  { s <- summary(chum_samples_kept[, chum_scalar_params])
  tibble(Stock = "Chum", parameter = rownames(s$statistics), mean = s$statistics[, "Mean"],
         sd = s$statistics[, "SD"], q2.5 = s$quantiles[, "2.5%"], median = s$quantiles[, "50%"],
         q97.5 = s$quantiles[, "97.5%"]) }
)
print(posterior_summary, n = Inf)

# Quick check: how close is each stock's posterior b to msmin (the
# truncation bound)? Same diagnostic as the sockeye bounded-b models --
# a posterior 2.5th percentile close to msmin means the bound is
# actively constraining the posterior, not just a formality.
bound_check <- bind_rows(
  { priors <- chilcotin_ricker_priors; msmin <- 1/priors$Smaxmax
  b_summary <- posterior_summary %>% filter(Stock == "Chilcotin", parameter == "b")
  tibble(Stock = "Chilcotin", msmin = msmin, b_q2.5 = b_summary$q2.5,
         pct_above_bound = 100 * (b_summary$q2.5 - msmin) / msmin) },
  { priors <- thompson_ricker_priors; msmin <- 1/priors$Smaxmax
  b_summary <- posterior_summary %>% filter(Stock == "Thompson", parameter == "b")
  tibble(Stock = "Thompson", msmin = msmin, b_q2.5 = b_summary$q2.5,
         pct_above_bound = 100 * (b_summary$q2.5 - msmin) / msmin) },
  { priors <- chum_ricker_priors; msmin <- 1/priors$Smaxmax
  b_summary <- posterior_summary %>% filter(Stock == "Chum", parameter == "b")
  tibble(Stock = "Chum", msmin = msmin, b_q2.5 = b_summary$q2.5,
         pct_above_bound = 100 * (b_summary$q2.5 - msmin) / msmin) }
)
print(bound_check, n = Inf)

# ============================================================
# 6. TIME-VARYING PRODUCTIVITY PLOTS, all three
# ============================================================

chilcotin_combined <- as.matrix(do.call(rbind, chilcotin_samples_kept))
thompson_combined  <- as.matrix(do.call(rbind, thompson_samples_kept))
chum_combined      <- as.matrix(do.call(rbind, chum_samples_kept))

ayst_all <- bind_rows(
  extract_year_series(chilcotin_combined, "ayst", chilcotin_data$Year) %>% mutate(Stock = "Chilcotin"),
  extract_year_series(thompson_combined,  "ayst", thompson_data$Year)  %>% mutate(Stock = "Thompson"),
  extract_year_series(chum_combined,      "ayst", chum_data$Year)      %>% mutate(Stock = "Chum")
)

ggplot(ayst_all, aes(Year, mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "#4682B4", alpha = 0.25) +
  geom_line(color = "#4682B4", linewidth = 1) +
  facet_wrap(~ Stock, scales = "free", ncol = 1) +
  labs(x = "Year", y = "Productivity (ayst)",
       title = "Time-varying productivity: Chilcotin & Thompson steelhead, Chum") +
  theme_minimal()

ggsave("figures/steelhead_chum_productivity.png", width = 8, height = 10, dpi = 600)


# SAVE trace plots -- intercept + spawners beta + covariate coefficients,
# steelhead & chum -- with panels relabeled to the real covariate names

steelhead_chum_fits <- list(
  Chilcotin = chilcotin_samples_kept,   # s = SL, t = SST
  Thompson  = thompson_samples_kept,    # s = SL, t = SST
  Chum      = chum_samples_kept         # s = SL, t = PDO
)

# What each raw parameter name (intercept, b, s, t) actually represents,
# per stock -- used to relabel the trace/density panel titles directly,
# not just the filename
param_labels <- list(
  Chilcotin = c(intercept = "Intercept", b = "Spawners beta (b)", s = "SL", t = "SST"),
  Thompson  = c(intercept = "Intercept", b = "Spawners beta (b)", s = "SL", t = "SST"),
  Chum      = c(intercept = "Intercept", b = "Spawners beta (b)", s = "SL", t = "PDO")
)

for (stock_name in names(steelhead_chum_fits)) {
  
  samples <- steelhead_chum_fits[[stock_name]]
  all_names <- varnames(samples)
  
  coef_names <- c("intercept", "b", "s", "t")
  coef_names <- coef_names[coef_names %in% all_names]
  
  sub_samples <- samples[, coef_names]
  
  # Relabel each chain's columns to the real covariate names, so the
  # actual plot panels (not just the filename) show what's being traced
  new_labels <- unname(param_labels[[stock_name]][coef_names])
  for (j in seq_along(sub_samples)) {
    colnames(sub_samples[[j]]) <- new_labels
  }
  
  covariate_suffix <- paste(param_labels[[stock_name]][c("s", "t")], collapse = "_")
  filename <- paste0(stock_name, "_intercept_spawnersB_", covariate_suffix, ".png")
  
  png(
    filename = file.path("Figures/diagnostics", filename),
    width = 2000,
    height = 3000,
    res = 300
  )
  
  plot(sub_samples)
  
  dev.off()
}
