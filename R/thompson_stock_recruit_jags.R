# ============================================================
# Thompson steelhead stock-recruit model -- JAGS
#
# Built analogously to the Chilcotin model (chilcotin_stock_recruit_jags.R),
# which was translated directly from the supplied WinBUGS files. No
# Thompson-specific .odc files were provided, so this follows the same
# prior structure and script logic, using your requested covariate set:
# Spawners (density dependence), SL, and SST -- no max_flow.
#
# See chilcotin_stock_recruit_jags.R for the full WinBUGS -> JAGS mapping
# notes (update/dic.set/beg etc.) -- same translation logic applies here.
#
# Requires: JAGS itself installed (https://mcmc-jags.sourceforge.io/) plus
# the R packages below.
# ============================================================

library(rjags)
library(coda)
library(dplyr)
library(readr)
library(ggplot2)

set.seed(1)

# ============================================================
# DATA
# ============================================================

steelhead_data <- read_csv("data/steelhead_data.csv")

thompson_data <- steelhead_data %>%
  filter(Stock == "Thompson",
         !is.na(Spawners), !is.na(lnRS), !is.na(SST), !is.na(SL)) %>%
  arrange(Year)

# 41 complete rows in the current data
stopifnot(nrow(thompson_data) == 41)

# ============================================================
# MODEL
# ============================================================

thompson_model_string <- "
model{
  # ---- priors (same structure as the Chilcotin model) ----
  slope     ~ dunif(0.0001, 1000)
  intercept ~ dnorm(0, 0.000001)
  LSE       ~ dunif(LminSE, LmaxSE)
  LminSE    <- log(0.01)
  LmaxSE    <- log(10)
  SE        <- exp(LSE)
  tau       <- 1 / (SE * SE)

  # ---- covariate coefficients: Spawners, SL, SST only ----
  s ~ dnorm(0, 0.0001)   # SL (sea lion index) coefficient
  t ~ dnorm(0, 0.0001)   # SST coefficient

  a <- intercept
  b <- slope

  PPslope <- step(intercept)   # P(intercept > 0)
  pvalues <- step(-s)          # P(s < 0)
  pvaluet <- step(-t)          # P(t < 0)

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

# ============================================================
# DATA LIST, INITS, MONITORED PARAMETERS
# ============================================================

thompson_jags_data <- list(
  ndata = nrow(thompson_data),
  sp    = thompson_data$Spawners,
  lnrs  = thompson_data$lnRS,
  SL    = thompson_data$SL,
  SST   = thompson_data$SST
)

# s/t starting values carried over from the Chilcotin inits_m1_4 files as a
# reasonable starting point (same covariates, same rough prior belief about
# sign/magnitude); slope/intercept/LSE nudged slightly for chain separation
thompson_inits <- list(
  list(LSE = 0.9, slope = 1.05, intercept = 1.1, s = -0.75, t = -0.127),
  list(LSE = 0.8, slope = 1.09, intercept = 1.2, s = -0.61, t = -0.137)
)

thompson_params <- c("slope", "intercept", "a", "b", "s", "t", "SE",
                      "PPslope", "pvalue", "pvalues", "pvaluet",
                      "ayst", "cst", "smax", "smsyst", "rmsyst", "msyst", "umsyst",
                      "deviance")

# ============================================================
# RUN
# ============================================================

thompson_jm <- jags.model(textConnection(thompson_model_string),
                           data = thompson_jags_data,
                           inits = thompson_inits,
                           n.chains = 2, n.adapt = 1000)

update(thompson_jm, n.iter = 2000)

thompson_dic <- dic.samples(thompson_jm, n.iter = 2000, type = "pD")

thompson_samples <- coda.samples(thompson_jm,
                                  variable.names = thompson_params,
                                  n.iter = 30000, thin = 1)

thompson_samples_kept <- window(thompson_samples, start = 2001)

# ============================================================
# DIAGNOSTICS
# ============================================================

scalar_params <- c("a", "b", "s", "t", "SE", "smax", "deviance")

print(gelman.diag(thompson_samples_kept[, scalar_params], multivariate = FALSE))
print(effectiveSize(thompson_samples_kept[, scalar_params]))
print(thompson_dic)

plot(thompson_samples_kept[, scalar_params])

# ============================================================
# POSTERIOR SUMMARY
# ============================================================

thompson_summary <- summary(thompson_samples_kept[, scalar_params])
thompson_param_table <- data.frame(
  parameter = rownames(thompson_summary$statistics),
  mean      = thompson_summary$statistics[, "Mean"],
  sd        = thompson_summary$statistics[, "SD"],
  q2.5      = thompson_summary$quantiles[, "2.5%"],
  median    = thompson_summary$quantiles[, "50%"],
  q97.5     = thompson_summary$quantiles[, "97.5%"]
)
print(thompson_param_table)

combined <- as.matrix(do.call(rbind, thompson_samples_kept))

extract_year_series <- function(prefix, years) {
  cols <- grep(paste0("^", prefix, "\\["), colnames(combined), value = TRUE)
  cols <- cols[order(as.numeric(gsub(paste0(prefix, "\\[|\\]"), "", cols)))]
  data.frame(
    Year    = years,
    metric  = prefix,
    mean    = colMeans(combined[, cols, drop = FALSE]),
    q2.5    = apply(combined[, cols, drop = FALSE], 2, quantile, 0.025),
    q97.5   = apply(combined[, cols, drop = FALSE], 2, quantile, 0.975)
  )
}

thompson_year_summary <- bind_rows(
  extract_year_series("ayst",    thompson_data$Year),
  extract_year_series("smsyst",  thompson_data$Year),
  extract_year_series("rmsyst",  thompson_data$Year),
  extract_year_series("msyst",   thompson_data$Year),
  extract_year_series("umsyst",  thompson_data$Year)
)

ggplot(thompson_year_summary %>% filter(metric == "ayst"),
       aes(Year, mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "#F8766D", alpha = 0.25) +
  geom_line(color = "#F8766D", linewidth = 1) +
  labs(x = "Year", y = "Productivity (ayst)",
       title = "Thompson steelhead: time-varying productivity (Spawners, SL, SST)") +
  theme_minimal()
