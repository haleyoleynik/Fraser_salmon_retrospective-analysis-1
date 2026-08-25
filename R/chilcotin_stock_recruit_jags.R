# ============================================================
# Chilcotin steelhead stock-recruit model -- WinBUGS -> JAGS
#
# Source: model_cow_m1_2-v4.odc / data_m1_1a.odc / data_m1_1b-v4.odc /
#         inits_1_m1_4.odc / inits_2_m1_4.odc (Chilcotin folder)
#
# The original WinBUGS model regressed ln(R/S) on spawners (Ricker
# density-dependence) plus SSL, SST, PDOJ, and NPGO as covariates on
# productivity. Per your request, this version drops PDOJ and NPGO and
# adds max_flow, so the covariate set is: Spawners (density dependence),
# max_flow, SL, SST.
#
# WinBUGS script -> JAGS mapping (script_cow_m1_2-v4.odc):
#   check(...) / data(...) / compile(2)      -> jags.model()
#   inits(1,...) / inits(2,...) / gen.inits()-> inits list passed to jags.model()
#                                                (JAGS auto-fills any node
#                                                you don't give a value for)
#   update(2000)                             -> update(jm, 2000)          [burn-in]
#   dic.set() / dic.stats()                  -> dic.samples(jm, ...)      [post-hoc DIC]
#   update(30000)                            -> coda.samples(jm, ..., n.iter = 30000)
#   beg(2001)                                -> window(samples, start = 2001)
#                                                [extra burn-in discard within the
#                                                 monitored run, as in the original]
#   set(<param>)                             -> variable.names in coda.samples()
#
# Requires: JAGS itself installed (https://mcmc-jags.sourceforge.io/) plus
# the R packages below. JAGS is a standalone program, not just an R package.
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

chilcotin_data <- steelhead_data %>%
  filter(Stock == "Chilcotin",
         !is.na(Spawners), !is.na(lnRS), !is.na(SST), !is.na(SL), !is.na(max_flow)) %>%
  arrange(Year)

# Should match ndata = 45 from the original data_m1_1a.odc
stopifnot(nrow(chilcotin_data) == 45)

# ============================================================
# MODEL
# ============================================================

chilcotin_model_string <- "
model{
  # ---- priors (unchanged from the WinBUGS original) ----
  slope     ~ dunif(0.0001, 1000)      # max Smax of 100,000 (sp already in 1000s)
  intercept ~ dnorm(0, 0.000001)       # prior SD = 1000
  LSE       ~ dunif(LminSE, LmaxSE)    # quasi-noninformative prior on SE
  LminSE    <- log(0.01)               # uniform on the log of SE
  LmaxSE    <- log(10)                 # range spans orders of magnitude
  SE        <- exp(LSE)                # SE is a nuisance scale parameter
  tau       <- 1 / (SE * SE)

  # ---- covariate coefficients (replaces PDOJ/NPGO with max_flow) ----
  s ~ dnorm(0, 0.0001)   # SL (sea lion index) coefficient
  t ~ dnorm(0, 0.0001)   # SST coefficient
  f ~ dnorm(0, 0.0001)   # max_flow coefficient

  a <- intercept
  b <- slope

  PPslope <- step(intercept)   # P(intercept > 0)
  pvalues <- step(-s)          # P(s < 0)
  pvaluet <- step(-t)          # P(t < 0)
  pvaluef <- step(-f)          # P(f < 0) -- direction not pre-specified for flow;
                                # look at both pvaluef and 1-pvaluef when interpreting

  smax <- 1 / b   # spawners producing maximum recruitment

  for (i in 1:ndata) {
    lnrs_pred[i] <- a - b * sp[i] + s * SL[i] + t * SST[i] + f * max_flow[i]
    lnrs[i]      ~ dnorm(lnrs_pred[i], tau)          # likelihood
    lnrs_rep[i]  ~ dnorm(lnrs_pred[i], tau)          # posterior predictive replicate
    pvalue[i]    <- step(lnrs_rep[i] - lnrs[i])      # Bayesian p-value, per point

    ayst[i]    <- a + s * SL[i] + t * SST[i] + f * max_flow[i]   # time-varying productivity
    cstx[i]    <- ayst[i] / b                                    # carrying capacity
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

chilcotin_jags_data <- list(
  ndata    = nrow(chilcotin_data),
  sp       = chilcotin_data$Spawners,
  lnrs     = chilcotin_data$lnRS,
  SL       = chilcotin_data$SL,
  SST      = chilcotin_data$SST,
  max_flow = chilcotin_data$max_flow
)

# s/t starting values taken directly from inits_1_m1_4.odc / inits_2_m1_4.odc
# (the covariate-model inits). f (max_flow) has no historical starting value
# since it's a new covariate here -- I picked small, chain-distinguishing
# starting values (0 and 0.1); adjust if you have a better prior guess.
chilcotin_inits <- list(
  list(LSE = 0.9, slope = 1.05, intercept = 1.1, s = -0.75, t = -0.127, f = 0.0),
  list(LSE = 0.8, slope = 1.09, intercept = 1.2, s = -0.61, t = -0.137, f = 0.1)
)

chilcotin_params <- c("slope", "intercept", "a", "b", "s", "t", "f", "SE",
                       "PPslope", "pvalue", "pvalues", "pvaluet", "pvaluef",
                       "ayst", "cst", "smax", "smsyst", "rmsyst", "msyst", "umsyst",
                       "deviance")

# ============================================================
# RUN
# ============================================================

chilcotin_jm <- jags.model(textConnection(chilcotin_model_string),
                            data = chilcotin_jags_data,
                            inits = chilcotin_inits,
                            n.chains = 2, n.adapt = 1000)

update(chilcotin_jm, n.iter = 2000)   # burn-in, matches update(2000)

chilcotin_dic <- dic.samples(chilcotin_jm, n.iter = 2000, type = "pD")  # matches dic.set()/dic.stats()

chilcotin_samples <- coda.samples(chilcotin_jm,
                                   variable.names = chilcotin_params,
                                   n.iter = 30000, thin = 1)            # matches update(30000)

chilcotin_samples_kept <- window(chilcotin_samples, start = 2001)      # matches beg(2001)

# ============================================================
# DIAGNOSTICS
# ============================================================

scalar_params <- c("a", "b", "s", "t", "f", "SE", "smax", "deviance")

print(gelman.diag(chilcotin_samples_kept[, scalar_params], multivariate = FALSE))
print(effectiveSize(chilcotin_samples_kept[, scalar_params]))
print(chilcotin_dic)

# traceplots for the key parameters
plot(chilcotin_samples_kept[, scalar_params])

# ============================================================
# POSTERIOR SUMMARY
# ============================================================

chilcotin_summary <- summary(chilcotin_samples_kept[, scalar_params])
chilcotin_param_table <- data.frame(
  parameter = rownames(chilcotin_summary$statistics),
  mean      = chilcotin_summary$statistics[, "Mean"],
  sd        = chilcotin_summary$statistics[, "SD"],
  q2.5      = chilcotin_summary$quantiles[, "2.5%"],
  median    = chilcotin_summary$quantiles[, "50%"],
  q97.5     = chilcotin_summary$quantiles[, "97.5%"]
)
print(chilcotin_param_table)

# per-year reference points (ayst, cst, smsyst, rmsyst, msyst, umsyst), joined
# back to Year, for downstream plotting/comparison work
combined <- as.matrix(do.call(rbind, chilcotin_samples_kept))

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

chilcotin_year_summary <- bind_rows(
  extract_year_series("ayst",    chilcotin_data$Year),
  extract_year_series("smsyst",  chilcotin_data$Year),
  extract_year_series("rmsyst",  chilcotin_data$Year),
  extract_year_series("msyst",   chilcotin_data$Year),
  extract_year_series("umsyst",  chilcotin_data$Year)
)

# Quick look: time-varying productivity (ayst) with 95% CI
ggplot(chilcotin_year_summary %>% filter(metric == "ayst"),
       aes(Year, mean)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5), fill = "#4682B4", alpha = 0.25) +
  geom_line(color = "#4682B4", linewidth = 1) +
  labs(x = "Year", y = "Productivity (ayst)",
       title = "Chilcotin steelhead: time-varying productivity (Spawners, max_flow, SL, SST)") +
  theme_minimal()
