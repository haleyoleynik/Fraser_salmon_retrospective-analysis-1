# ============================================================
# Sockeye stock-recruit JAGS models -- one per stock (15 stocks, Late
# Shuswap excluded), using each stock's top-performing model
# (deltaAIC = 0) from Fraser_sockeye_dredge-results_alph05.csv (stricter
# model-selection alpha, matching the coho v2 update).
#
# COVARIATE CHANGES vs. the original alpha~0.10 results:
#   Bowron: PDO, SeaLions, pink, NPGO  ->  SeaLions, pink
#   Scotch: SeaLions, pink, NPGO       ->  pink, NPGO
#   (all other 13 stocks unchanged)
#
# BOUNDED RICKER b PRIOR -- per your Skeena methodology:
#   1. b_mle = -spawners coefficient from the dredge table (frequentist
#      MLE, already computed via lm()/dredge outside this script)
#   2. Smax_mle = 1 / b_mle
#   3. Smaxmax = 5 * Smax_mle  (your stated factor, absent lake-capacity
#      information)
#   4. Lower prior bound for b = 1 / Smaxmax
#   5. b ~ dlnorm(log(1/prSmax), 1/prCV^2) T(1/Smaxmax, )  -- lognormal,
#      median at 1/prSmax (i.e. at b_mle), truncated below at 1/Smaxmax
#      so Smax can't wander above Smaxmax. prCV = 2 (vague), matching
#      the Skeena WinBUGS convention exactly.
#
# This REPLACES the previous b_spawners ~ dnorm(0, 0.0001) (unbounded)
# prior. The linear predictor also changes form to match: previously
# `+ b_spawners * spawners[i]` (raw coefficient, sign absorbed);
# now `- b * spawners[i]` (b is a positive slope by construction, same
# convention as the Skeena/coho bounded-b models).
#
# DOWNSTREAM IMPACT: any script reading "b_spawners" from the posterior
# (build_stock_draws(), compute_smsy_benchmark(), etc.) needs updating to
# read "b" instead, and WITHOUT negating it -- b is already the positive
# Ricker slope. This is also the fix for the sign bug found in the
# existing ratio/Smsy scripts (see chat) -- see
# sockeye_ratio_smsy_full_script_v2.R for the updated downstream script.
# ============================================================

library(rjags)
library(coda)
library(dplyr)
library(readr)

set.seed(1)

sockeye_data_raw <- read_csv("data/sockeye_data.csv") %>% rename(BroodYear = yr)

# Top-model (deltaAIC = 0) extra covariates per stock, from
# Fraser_sockeye_dredge-results_alph05.csv. spawners is in every stock's
# top model, so it's not listed here -- added automatically below.
STOCK_COVARIATES <- list(
  Birkenhead     = c("PDO", "SeaLions"),
  Bowron         = c("SeaLions", "pink"),
  Chilko         = c("PDO", "SeaLions"),
  Cultus         = c("SeaLions", "smolt.sst"),
  `Early Stuart` = c("seal", "adult.sst"),
  Gates          = c("smolt.sst"),
  `Late Stuart`  = c("pink"),
  Pitt           = c("SeaLions", "smolt.sst"),
  Portage        = c("pink"),
  Quesnel        = c("adult.sst"),
  Raft           = c("SeaLions", "smolt.sst"),
  Scotch         = c("pink", "NPGO"),
  Seymour        = c("pink", "smolt.sst"),
  Stellako       = c("SeaLions", "pink"),
  Weaver         = c("pink")
)

# Bounded Ricker b priors per stock, derived from the alpha<0.05 dredge
# table's own spawners MLE (b_mle = -spawners coefficient), per the
# Skeena methodology: prSmax = Smax_mle = 1/b_mle, Smaxmax = 5*prSmax,
# prCV = 2. These are NOT independent expert priors -- they're built
# directly from this same fit's own point estimate, same as your WinBUGS
# example. If you have actual prior beliefs about plausible Smax per
# stock (habitat capacity etc.), replace these.
RICKER_B_PRIORS <- list(
  Birkenhead     = list(prSmax = 174771.0,   Smaxmax = 873854.8,    prCV = 2),
  Bowron         = list(prSmax = 20404.9,    Smaxmax = 102024.7,    prCV = 2),
  Chilko         = list(prSmax = 1085592.7,  Smaxmax = 5427963.7,   prCV = 2),
  Cultus         = list(prSmax = 45043.9,    Smaxmax = 225219.4,    prCV = 2),
  `Early Stuart` = list(prSmax = 916408.2,   Smaxmax = 4582040.8,   prCV = 2),
  Gates          = list(prSmax = 87073.0,    Smaxmax = 435364.9,    prCV = 2),
  `Late Stuart`  = list(prSmax = 2025911.1,  Smaxmax = 10129555.4,  prCV = 2),
  Pitt           = list(prSmax = 51945.7,    Smaxmax = 259728.5,    prCV = 2),
  Portage        = list(prSmax = 40242.5,    Smaxmax = 201212.6,    prCV = 2),
  Quesnel        = list(prSmax = 6466792.4,  Smaxmax = 32333962.1,  prCV = 2),
  Raft           = list(prSmax = 46988.4,    Smaxmax = 234941.9,    prCV = 2),
  Scotch         = list(prSmax = 256122.0,   Smaxmax = 1280609.9,   prCV = 2),
  Seymour        = list(prSmax = 348909.2,   Smaxmax = 1744546.1,   prCV = 2),
  Stellako       = list(prSmax = 430842.1,   Smaxmax = 2154210.4,   prCV = 2),
  Weaver         = list(prSmax = 203076.4,   Smaxmax = 1015382.2,   prCV = 2)
)

jags_name <- function(cv) gsub("\\.", "_", cv)

build_sockeye_model_string <- function(extra_covariates) {
  
  extra_priors <- paste(sprintf("  b_%s ~ dnorm(0, 0.0001)   # %s coefficient",
                                jags_name(extra_covariates), extra_covariates),
                        collapse = "\n")
  
  extra_terms <- paste(sprintf(" + b_%s * %s[i]",
                               jags_name(extra_covariates), jags_name(extra_covariates)),
                       collapse = "")
  
  sprintf('
model{
  intercept ~ dnorm(0, 0.000001)
  LSE       ~ dunif(LminSE, LmaxSE)
  LminSE    <- log(0.01)
  LmaxSE    <- log(10)
  SE        <- exp(LSE)
  tau       <- 1 / (SE * SE)

  ms    <- log(1 / prSmax)
  msmin <- 1 / Smaxmax
  taus  <- 1 / (prCV * prCV)
  b     ~ dlnorm(ms, taus) T(msmin, )
  Smax  <- 1 / b

%s

  for (i in 1:ndata) {
    lnrs_pred[i] <- intercept - b * spawners[i]%s
    lnrs[i]      ~ dnorm(lnrs_pred[i], tau)
    lnrs_rep[i]  ~ dnorm(lnrs_pred[i], tau)
    pvalue[i]    <- step(lnrs_rep[i] - lnrs[i])
  }
}
', extra_priors, extra_terms)
}

fit_sockeye_stock <- function(stock_name, extra_covariates) {
  
  message("Fitting: ", stock_name, " -- extra covariates: ",
          if (length(extra_covariates) == 0) "(none)" else paste(extra_covariates, collapse = ", "))
  
  stock_data <- sockeye_data_raw %>%
    filter(Stock == stock_name) %>%
    arrange(BroodYear)
  
  needed_cols <- c("lnrs", "spawners", extra_covariates)
  stock_data <- stock_data %>% filter(if_all(all_of(needed_cols), ~ !is.na(.)))
  
  model_string <- build_sockeye_model_string(extra_covariates)
  
  priors <- RICKER_B_PRIORS[[stock_name]]
  
  jags_data <- c(
    list(ndata = nrow(stock_data), lnrs = stock_data$lnrs, spawners = stock_data$spawners,
         prSmax = priors$prSmax, Smaxmax = priors$Smaxmax, prCV = priors$prCV),
    setNames(as.list(stock_data[extra_covariates]), jags_name(extra_covariates))
  )
  
  # b starting values: the prior median (1/prSmax, i.e. the MLE itself)
  # for chain 1, a nearby but distinct value for chain 2 -- both
  # comfortably above msmin (= b_mle/5) since Smaxmax is 5x Smax_mle
  b_init_1 <- 1 / priors$prSmax
  b_init_2 <- 1.3 / priors$prSmax
  
  inits <- list(
    list(LSE = 0.9, intercept = 1.0, b = b_init_1),
    list(LSE = 0.8, intercept = 1.2, b = b_init_2)
  )
  
  params <- c("intercept", "b", "Smax", paste0("b_", jags_name(extra_covariates)),
              "SE", "pvalue", "deviance")
  
  jm <- jags.model(textConnection(model_string), data = jags_data, inits = inits,
                   n.chains = 2, n.adapt = 1000)
  update(jm, n.iter = 2000)
  
  dic <- dic.samples(jm, n.iter = 2000, type = "pD")
  
  samples <- coda.samples(jm, variable.names = params, n.iter = 30000, thin = 1)
  samples_kept <- window(samples, start = 2001)
  
  list(stock = stock_name, extra_covariates = extra_covariates,
       jags_model = jm, samples = samples_kept, dic = dic, data = stock_data)
}

# ------------------------------------------------------------
# FIT ALL 15 STOCKS
# ------------------------------------------------------------

sockeye_fits <- lapply(names(STOCK_COVARIATES), function(stock_name) {
  fit_sockeye_stock(stock_name, STOCK_COVARIATES[[stock_name]])
})
names(sockeye_fits) <- names(STOCK_COVARIATES)

# ------------------------------------------------------------
# DIAGNOSTICS -- same checks as before. Watch b and Smax specifically:
# if a posterior mass piles up right against msmin (the truncation
# point), that's the bound actually binding -- worth knowing, since it
# means the data alone couldn't rule out an even-larger Smax and the
# prior is doing real work, not just a formality.
# ------------------------------------------------------------

for (stock_name in names(sockeye_fits)) {
  fit <- sockeye_fits[[stock_name]]
  all_names <- varnames(fit$samples)
  scalar_params <- setdiff(all_names, c("deviance", grep("^pvalue\\[", all_names, value = TRUE)))
  
  cat("\n====", stock_name, "====\n")
  print(gelman.diag(fit$samples[, scalar_params], multivariate = FALSE))
  print(effectiveSize(fit$samples[, scalar_params]))
  print(fit$dic)
}

# ------------------------------------------------------------
# SAVE POSTERIOR SAMPLES -- one RDS per stock
# ------------------------------------------------------------

for (stock_name in names(sockeye_fits)) {
  saveRDS(sockeye_fits[[stock_name]]$samples,
          paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds"))
}

# Posterior summary table, all stocks -- names pulled from samples directly
sockeye_posterior_summary <- bind_rows(lapply(names(sockeye_fits), function(stock_name) {
  fit <- sockeye_fits[[stock_name]]
  all_names <- varnames(fit$samples)
  coef_names <- setdiff(all_names, c("deviance", grep("^pvalue\\[", all_names, value = TRUE)))
  s <- summary(fit$samples[, coef_names])
  tibble(
    Stock     = stock_name,
    parameter = rownames(s$statistics),
    mean      = s$statistics[, "Mean"],
    sd        = s$statistics[, "SD"],
    q2.5      = s$quantiles[, "2.5%"],
    median    = s$quantiles[, "50%"],
    q97.5     = s$quantiles[, "97.5%"]
  )
}))
print(sockeye_posterior_summary, n = Inf)

# Quick check: how close is each stock's posterior b to msmin (the
# truncation bound)? A posterior 2.5th percentile very close to msmin
# suggests the bound is actively constraining the posterior, not just a
# formality.
bound_check <- bind_rows(lapply(names(sockeye_fits), function(stock_name) {
  priors <- RICKER_B_PRIORS[[stock_name]]
  msmin <- 1 / priors$Smaxmax
  b_summary <- sockeye_posterior_summary %>% filter(Stock == stock_name, parameter == "b")
  tibble(Stock = stock_name, msmin = msmin, b_q2.5 = b_summary$q2.5,
         pct_above_bound = 100 * (b_summary$q2.5 - msmin) / msmin)
}))
print(bound_check, n = Inf)

# ============================================================
# Sockeye: time-varying productivity (intercept + covariate effects,
# EXCLUDING the spawners/density-dependence term), faceted by stock,
# with 95% credible intervals -- same approach as the steelhead/coho
# productivity plots.
#
# Uses each stock's own real historical covariate values (not part of
# the low-point/MC retrospective work) -- this is just the fitted
# productivity trajectory itself, for sanity-checking the JAGS models.
#
# Requires sockeye_fits and STOCK_COVARIATES from sockeye_jags_models.R
# (or readRDS the saved posterior samples per stock if sockeye_fits
# isn't in your session).
# ============================================================

library(tidyverse)
library(coda)

STOCK_COVARIATES <- list(
  Birkenhead     = c("PDO", "SeaLions"),
  Bowron         = c("SeaLions", "pink"),
  Chilko         = c("PDO", "SeaLions"),
  Cultus         = c("SeaLions", "smolt.sst"),
  `Early Stuart` = c("seal", "adult.sst"),
  Gates          = c("smolt.sst"),
  `Late Stuart`  = c("pink"),
  Pitt           = c("SeaLions", "smolt.sst"),
  Portage        = c("pink"),
  Quesnel        = c("adult.sst"),
  Raft           = c("SeaLions", "smolt.sst"),
  Scotch         = c("pink", "NPGO"),
  Seymour        = c("pink", "smolt.sst"),
  Stellako       = c("SeaLions", "pink"),
  Weaver         = c("pink")
)

jags_name <- function(cv) gsub("\\.", "_", cv)

compute_ayst_sockeye <- function(stock_name) {
  
  fit <- if (exists("sockeye_fits")) {
    sockeye_fits[[stock_name]]
  } else {
    list(samples = readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds")))
  }
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  
  # stock_data: needs BroodYear + the real covariate columns for this
  # stock. If sockeye_fits[[stock]]$data exists (from this session's
  # fitting run), use it directly; otherwise fall back to re-deriving
  # the same complete-case slice from sockeye_data.csv.
  stock_data <- if (!is.null(fit$data)) {
    fit$data
  } else {
    sockeye_data_raw <- read_csv("sockeye_data.csv") %>% rename(BroodYear = yr)
    needed_cols <- c("lnrs", "spawners", extra_covariates)
    sockeye_data_raw %>%
      filter(Stock == stock_name) %>%
      arrange(BroodYear) %>%
      filter(if_all(all_of(needed_cols), ~ !is.na(.)))
  }
  
  all_names <- varnames(fit$samples)
  b_names <- setdiff(grep("^b_", all_names, value = TRUE), "b_spawners")  # exclude density dependence
  coef_names <- c("intercept", b_names)
  
  post <- as.matrix(fit$samples[, coef_names])
  
  # map each b_<name> back to the original data column (handles the
  # dot -> underscore conversion, e.g. b_adult_sst -> adult.sst)
  covariate_lookup <- setNames(extra_covariates, jags_name(extra_covariates))
  data_cols <- unname(covariate_lookup[sub("^b_", "", b_names)])
  
  X <- cbind(1, as.matrix(stock_data[, data_cols, drop = FALSE]))
  
  ayst_draws <- X %*% t(post)   # n_years x n_kept_draws
  
  tibble(
    Stock = stock_name,
    Year  = stock_data$BroodYear,
    mean  = rowMeans(ayst_draws),
    q2.5  = apply(ayst_draws, 1, quantile, 0.025),
    q97.5 = apply(ayst_draws, 1, quantile, 0.975)
  )
}

sockeye_ayst_summary <- bind_rows(lapply(names(STOCK_COVARIATES), compute_ayst_sockeye))

ggplot(sockeye_ayst_summary, aes(Year, mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "#4682B4", alpha = 0.25) +
  geom_line(color = "#4682B4", linewidth = 1) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 4) +
  labs(x = "Year", y = "Productivity (alpha, excl. density dependence)") +
  theme_minimal(base_size = 9)

ggsave("figures/sockeye_productivity_by_stock_bayesian.png", width = 14, height = 10, dpi = 600)



# SAVE trace plots --------------------------------------------------
dir.create("Figures", showWarnings = FALSE)

for(stock_name in names(sockeye_fits)) {
  
  fit <- sockeye_fits[[stock_name]]
  
  all_names <- varnames(fit$samples)
  
  # Intercept + spawner coefficient + all covariate slopes
  coef_names <- c(
    "intercept",
    "b",
    grep("^b_", all_names, value = TRUE)
  )
  
  # Keep only parameters that actually exist
  coef_names <- coef_names[coef_names %in% all_names]
  
  png(
    filename = file.path("Figures", paste0(stock_name, "_coefficients.png")),
    width = 2000,
    height = 3000,
    res = 300
  )
  
  
  plot(fit$samples[, coef_names])
  
  dev.off()
}

