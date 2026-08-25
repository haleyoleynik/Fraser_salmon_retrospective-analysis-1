# ============================================================
# Sockeye stock-recruit JAGS models -- one per stock (15 stocks, Late
# Shuswap excluded -- no usable top model), using each stock's
# top-performing model (deltaAIC = 0) from Fraser_sockeye_dredge-results.csv
#
# Same generic-builder approach as coho_jags_models.R: every top model
# includes spawners (density dependence); stocks differ only in which of
# NPGO/PDO/SeaLions/seal/adult.sst/pink/smolt.sst their top model adds.
#
# lnrs_pred[i] <- intercept + b_spawners*spawners[i] + <extra covariates>
# lnrs[i] ~ dnorm(lnrs_pred[i], tau)
#
# Maps directly onto the sockeye retrospective script's "terms" object:
#   ra = intercept, rb = -b_spawners, sel_covs/cov_coefs = the extras
# (get_top_model_terms() computes rb the same way from the frequentist
# spawners coefficient -- see sockeye_mc_probabilistic_retrospective.R)
# ============================================================

library(rjags)
library(coda)
library(dplyr)
library(readr)

set.seed(1)

sockeye_data_raw <- read_csv("data/sockeye_data.csv") %>% rename(BroodYear = yr)

# Top-model (deltaAIC = 0) extra covariates per stock, from
# Fraser_sockeye_dredge-results.csv. spawners is in every stock's top
# model, so it's not listed here -- added automatically below.
STOCK_COVARIATES <- list(
  Birkenhead   = c("PDO", "SeaLions"),
  Bowron       = c("PDO", "SeaLions", "pink", "NPGO"),
  Chilko       = c("PDO", "SeaLions"),
  Cultus       = c("SeaLions", "smolt.sst"),
  `Early Stuart` = c("seal", "adult.sst"),
  Gates        = c("smolt.sst"),
  `Late Stuart`  = c("pink"),
  Pitt         = c("SeaLions", "smolt.sst"),
  Portage      = c("pink"),
  Quesnel      = c("adult.sst"),
  Raft         = c("SeaLions", "smolt.sst"),
  Scotch       = c("SeaLions", "pink", "NPGO"),
  Seymour      = c("pink", "smolt.sst"),
  Stellako     = c("SeaLions", "pink"),
  Weaver       = c("pink")
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

  b_spawners ~ dnorm(0, 0.0001)
%s

  for (i in 1:ndata) {
    lnrs_pred[i] <- intercept + b_spawners * spawners[i]%s
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
  
  jags_data <- c(
    list(ndata = nrow(stock_data), lnrs = stock_data$lnrs, spawners = stock_data$spawners),
    setNames(as.list(stock_data[extra_covariates]), jags_name(extra_covariates))
  )
  
  inits <- list(
    list(LSE = 0.9, intercept = 1.0, b_spawners = 0),
    list(LSE = 0.8, intercept = 1.2, b_spawners = 0)
  )
  
  params <- c("intercept", "b_spawners", paste0("b_", jags_name(extra_covariates)),
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
# DIAGNOSTICS -- same checks as coho/steelhead: trace, Gelman-Rubin,
# effective size, DIC. Names pulled directly from the samples (not
# reconstructed), same fix already needed twice for coho.
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


# plot trace plots ------------------------------------------
stock_name <- "Birkenhead"  # change as needed
fit <- sockeye_fits[[stock_name]]
coef_names <- setdiff(varnames(fit$samples), c("deviance", grep("^pvalue\\[", varnames(fit$samples), value = TRUE)))
plot(fit$samples[, coef_names])

# Quick sanity check against the frequentist dredge coefficients --
# compare sockeye_posterior_summary means to the "(Intercept)"/spawners/
# covariate values in Fraser_sockeye_dredge-results.csv the same way we
# validated the coho models