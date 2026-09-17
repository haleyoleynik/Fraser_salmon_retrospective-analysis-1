# Sockeye retrospective model — scenario comparisons
# Adapted from Carl Walters Excel model
# Haley Oleynik
#
# Three covariate-driven scenarios: pinniped, SST, pink salmon. Each runs
# through the same pipeline so catch-lost and abundance-recovery estimates
# are directly comparable across drivers.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# ============================================================
# CONFIG
# ============================================================

FIT_YEARS      <- 1952:2019
LAG_YEARS      <- 4
COVARIATE_COLS <- c("NPGO", "PDO", "SeaLions", "seal", "adult.sst", "pink", "smolt.sst")

FREEZE_YEAR        <- 1970        # pinniped covariates frozen at this year's level, forward
SST_BASELINE_YEARS <- 1950:1975   # years used for stock-specific SST long-term mean

retroU_default   <- 0.3   # retrospective harvest rate cap (*_retro model runs)
useretro_default <- TRUE
yrretro_default  <- 1990  # year the retrospective harvest rate cap kicks in

STOCKS <- c("Birkenhead", "Bowron", "Chilko", "Cultus", "Early Stuart", "Gates",
            "Late Shuswap", "Late Stuart", "Pitt", "Portage", "Quesnel", "Raft",
            "Scotch", "Seymour", "Stellako", "Weaver")

# ============================================================
# LOAD DATA
# ============================================================

obs_raw <- read_csv("R/Sockeye Retrospective Shiny App/Walters_model_all-stocks.csv") %>%
  rename(
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year), Stock %in% STOCKS) %>%
  arrange(Stock, Year)

covariates_main <- read_csv("Data/sockeye_standardized_covariates.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, any_of(COVARIATE_COLS))

covariates_pink_wild <- read_csv("Data/sockeye_standardized_covariates_pink-wild.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, pink_wild)

dredge_models <- read_csv("Data/sockeye_top_models_dredge_wo-aquaculture.csv")

top_models <- dredge_models %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )

# ============================================================
# BUILD SCENARIO COVARIATES
# Each scenario gets its own "<var>_scenario" column(s), joined onto a
# single `obs` table so any driver can be swapped in at model-run time.
# ============================================================

## Pinniped: SeaLions/seal frozen at FREEZE_YEAR levels, forward
no_freeze_year <- covariates_main %>%
  group_by(Stock) %>%
  summarise(has_freeze_year = any(Year == FREEZE_YEAR & if_any(c(SeaLions, seal), ~ !is.na(.))),
            .groups = "drop") %>%
  filter(!has_freeze_year)

if (nrow(no_freeze_year) > 0) {
  warning("No ", FREEZE_YEAR, " SeaLions/seal value for stock(s): ",
          paste(no_freeze_year$Stock, collapse = ", "),
          " -- pinniped scenario will equal actual for these stocks.")
}

pinniped_scenario_cov <- covariates_main %>%
  group_by(Stock) %>%
  mutate(across(c(SeaLions, seal), function(x) {
    freeze_val <- x[Year == FREEZE_YEAR]
    if (length(freeze_val) != 1 || is.na(freeze_val)) x
    else ifelse(Year > FREEZE_YEAR, freeze_val, x)
  })) %>%
  ungroup() %>%
  select(Stock, Year, SeaLions_scenario = SeaLions, seal_scenario = seal)

## Pink: substituted with the wild-only pink series
no_alt_series <- covariates_pink_wild %>%
  group_by(Stock) %>%
  summarise(has_alt = any(!is.na(pink_wild)), .groups = "drop") %>%
  filter(!has_alt)

if (nrow(no_alt_series) > 0) {
  warning("No pink_wild values for stock(s): ",
          paste(no_alt_series$Stock, collapse = ", "),
          " -- pink scenario will equal actual for these stocks.")
}

pink_scenario_cov <- covariates_pink_wild %>%
  select(Stock, Year, pink_scenario = pink_wild)

## SST: adult.sst / smolt.sst frozen at each stock's 1950-1975 mean
sst_means <- covariates_main %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  group_by(Stock) %>%
  summarise(
    adult.sst_scenario = mean(adult.sst, na.rm = TRUE),
    smolt.sst_scenario = mean(smolt.sst, na.rm = TRUE),
    .groups = "drop"
  )

no_sst_mean <- sst_means %>% filter(is.na(adult.sst_scenario) | is.na(smolt.sst_scenario))
if (nrow(no_sst_mean) > 0) {
  warning("No usable ", min(SST_BASELINE_YEARS), "-", max(SST_BASELINE_YEARS),
          " SST mean for stock(s): ", paste(no_sst_mean$Stock, collapse = ", "))
}

sst_scenario_cov <- covariates_main %>%
  select(Stock, Year) %>%
  left_join(sst_means, by = "Stock")

obs <- obs_raw %>%
  left_join(covariates_main,       by = c("Stock", "Year")) %>%
  left_join(pinniped_scenario_cov, by = c("Stock", "Year")) %>%
  left_join(pink_scenario_cov,     by = c("Stock", "Year")) %>%
  left_join(sst_scenario_cov,      by = c("Stock", "Year"))

# ============================================================
# MODEL FUNCTIONS
# ============================================================

# Late Shuswap: modeled by 4-year cycle line (cycle = Year %% 4, per your
# cycle_df construction), not the single stock-wide intercept/spawners-
# slope structure every other stock uses -- so it needs its own lookup
# rather than fitting into top_models. Read directly from the cycle-
# stratified dredge output rather than hardcoding values, so a re-run of
# your dredge (Shuswap_sockeye_dredge-results.csv) is picked up
# automatically instead of requiring a manual re-transcription here.
late_shuswap_dredge <- read_csv("Data/Shuswap_sockeye_dredge-results.csv")

# Currently exactly one top model (deltaAIC=0). If a future re-dredge
# produces more than one model within your deltaAIC<=2 cutoff, this stops
# loudly rather than silently picking one arbitrarily or averaging them
# incorrectly -- decide how you want to combine them if/when that happens.
late_shuswap_top <- late_shuswap_dredge %>% filter(deltaAIC == min(deltaAIC))
stopifnot(
  "Expected exactly one Late Shuswap top model (deltaAIC=0) -- got more than one; decide how to combine them before proceeding." =
    nrow(late_shuswap_top) == 1
)

late_shuswap_intercept_base <- late_shuswap_top[["(Intercept)"]]

# Intercept: cycle 0 is the reference level, absorbed into (Intercept) --
# cycles 1-3 are (Intercept) + factor(cycle)N (0 if a cycle's delta
# wasn't selected/significant in a future re-dredge).
# Spawners slope: NO reference-level folding here, since the model has no
# bare "spawners" main effect, only the interaction -- all four cycles
# get their own explicit factor(cycle)N:spawners column directly.
late_shuswap_cycle_terms <- tibble(
  cycle = 0:3,
  ra = late_shuswap_intercept_base + c(
    0,
    coalesce(late_shuswap_top[["factor(cycle)1"]], 0),
    coalesce(late_shuswap_top[["factor(cycle)2"]], 0),
    coalesce(late_shuswap_top[["factor(cycle)3"]], 0)
  ),
  rb = -c(
    late_shuswap_top[["factor(cycle)0:spawners"]],
    late_shuswap_top[["factor(cycle)1:spawners"]],
    late_shuswap_top[["factor(cycle)2:spawners"]],
    late_shuswap_top[["factor(cycle)3:spawners"]]
  )
)

# Covariates: whichever of COVARIATE_COLS are actually present as columns
# in this file AND non-NA in the top model -- generalized rather than
# hardcoded to smolt.sst, so a future re-dredge selecting adult.sst (or
# anything else in COVARIATE_COLS) is picked up correctly too.
late_shuswap_covariate_cols_present <- intersect(COVARIATE_COLS, colnames(late_shuswap_top))
late_shuswap_sel_covs <- late_shuswap_covariate_cols_present[
  !is.na(as.numeric(late_shuswap_top[late_shuswap_covariate_cols_present]))
]
late_shuswap_cov_coefs <- if (length(late_shuswap_sel_covs) > 0) {
  as.numeric(late_shuswap_top[late_shuswap_sel_covs])
} else {
  numeric(0)
}

# AIC-weighted top-model coefficients for a stock; falls back to a plain
# Ricker fit if the stock has no dredge top model.
get_top_model_terms <- function(stock_name, fit_df_for_fallback) {
  
  if (stock_name == "Late Shuswap") {
    # ra/rb are cycle-varying, not scalar -- NA here deliberately, since
    # run_retro_model() branches on stock_name to pull the right
    # cycle-specific value from cycle_terms instead of using these
    return(list(
      ra = NA_real_, rb = NA_real_,
      sel_covs = late_shuswap_sel_covs, cov_coefs = late_shuswap_cov_coefs,
      cycle_terms = late_shuswap_cycle_terms
    ))
  }
  
  top_model_row <- top_models %>% filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
    list(
      ra        = top_model_row[["(Intercept)"]],
      rb        = -top_model_row[["spawners"]],
      sel_covs  = sel_covs,
      cov_coefs = if (length(sel_covs) > 0) as.numeric(top_model_row[sel_covs]) else numeric(0)
    )
  } else {
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df_for_fallback)
    list(ra = unname(coef(fit)[1]), rb = -unname(coef(fit)[2]),
         sel_covs = character(0), cov_coefs = numeric(0))
  }
}

# Linear-predictor contribution from selected covariates. scenario_vars
# names which of sel_covs should be read from their "<var>_scenario"
# column instead of the actual observed column.
compute_cov_term <- function(dat, sel_covs, cov_coefs, scenario_vars = character(0)) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  cols <- ifelse(sel_covs %in% scenario_vars, paste0(sel_covs, "_scenario"), sel_covs)
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}

# Retrospective stock-recruit + harvest projection for one stock.
# scenario_vars swaps covariates to their scenario column in the forward
# projection only -- process-error residuals (wt) are always estimated
# against actual historical covariates, so a scenario comparison isolates
# the covariate effect rather than mixing it with a different noise draw.
run_retro_model <- function(dat, stock_name, retroU, useretro, yrretro,
                            scenario_vars = character(0)) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs = pmin(0.95, Catch / RunJacks),
      Ut_obs = if_else(
        stock_name == "Late Shuswap" & Year == 2012,  # manual correction to match historical U
        pmin(0.7, Catch / RunJacks),
        Ut_obs
      ),
      ENS = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS, is.finite(lnR_S), is.finite(AdultEscapement))
  
  terms <- get_top_model_terms(stock_name, fit_df)
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars)
  
  if (stock_name == "Late Shuswap") {
    # Cycle-specific ra/rb joined in as their own columns, matched to
    # each row's own cycle -- keeps everything below (wt, the recursive
    # loop) working off real per-row values instead of needing a lookup
    # inside the loop each iteration.
    obs3 <- obs2 %>%
      mutate(cycle = Year %% 4, cov_term = cov_term, cov_term_proj = cov_term_proj) %>%
      left_join(terms$cycle_terms, by = "cycle") %>%
      mutate(wt = lnR_S - (ra - rb * AdultEscapement + cov_term))
  } else {
    obs3 <- obs2 %>%
      mutate(
        cov_term      = cov_term,
        cov_term_proj = cov_term_proj,
        wt = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term)
      )
  }
  
  n <- nrow(obs3)
  retroR     <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroS     <- rep(NA_real_, n)
  retroC     <- rep(NA_real_, n)
  retroU_vec <- if (useretro) ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs) else obs3$Ut_obs
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * (1 - retroU_vec[1:LAG_YEARS]) * obs3$ENS[1:LAG_YEARS]
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    ra_use <- if (stock_name == "Late Shuswap") obs3$ra[j] else terms$ra
    rb_use <- if (stock_name == "Late Shuswap") obs3$rb[j] else terms$rb
    retro_lnRS[j] <- ra_use - rb_use * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(retroR = retroR, retro_lnRS = retro_lnRS, retroU = retroU_vec,
           retroS = retroS, retroC = retroC)
}

run_scenario_model <- function(dat, retroU, useretro, yrretro, scenario_vars = character(0)) {
  dat %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x, stock_name = .y$Stock, retroU = retroU,
                                   useretro = useretro, yrretro = yrretro,
                                   scenario_vars = scenario_vars)) %>%
    ungroup()
}

# ============================================================
# RUN MODELS
# *_retro : capped retrospective harvest rate (productivity/return/spawner comparisons)
# *_hist  : actual historical harvest rate ("catch lost to driver" comparisons)
# ============================================================

model_actual_retro <- run_scenario_model(obs, retroU_default, useretro_default, yrretro_default)
model_actual_hist  <- run_scenario_model(obs, retroU = 0, useretro = FALSE, yrretro = yrretro_default)

model_pinniped_retro <- run_scenario_model(obs, retroU_default, useretro_default, yrretro_default,
                                           scenario_vars = c("SeaLions", "seal"))
model_pinniped_hist  <- run_scenario_model(obs, retroU = 0, useretro = FALSE, yrretro = yrretro_default,
                                           scenario_vars = c("SeaLions", "seal"))

model_pink_retro <- run_scenario_model(obs, retroU_default, useretro_default, yrretro_default,
                                       scenario_vars = "pink")
model_pink_hist  <- run_scenario_model(obs, retroU = 0, useretro = FALSE, yrretro = yrretro_default,
                                       scenario_vars = "pink")

model_sst_retro <- run_scenario_model(obs, retroU_default, useretro_default, yrretro_default,
                                      scenario_vars = c("adult.sst", "smolt.sst"))
model_sst_hist  <- run_scenario_model(obs, retroU = 0, useretro = FALSE, yrretro = yrretro_default,
                                      scenario_vars = c("adult.sst", "smolt.sst"))

# ============================================================
# CATCH LOST TO EACH DRIVER
# catch_lost = scenario catch - actual catch, at the actual historical
# harvest rate (positive = catch available if the driver's covariate(s)
# had stayed at the scenario level)
# ============================================================

compute_catch_lost <- function(model_scenario_hist, scenario_label) {
  model_actual_hist %>%
    select(Stock, Year, catch_actual = retroC) %>%
    left_join(
      model_scenario_hist %>% select(Stock, Year, catch_scenario = retroC),
      by = c("Stock", "Year")
    ) %>%
    mutate(
      catch_lost = catch_scenario - catch_actual,
      scenario   = scenario_label
    ) %>%
    arrange(Stock, Year) %>%
    group_by(Stock) %>%
    mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
    ungroup()
}

catch_lost_pinniped <- compute_catch_lost(model_pinniped_hist, "Pinniped")
catch_lost_pink     <- compute_catch_lost(model_pink_hist,     "Pink salmon")
catch_lost_sst      <- compute_catch_lost(model_sst_hist,      "SST")

catch_lost_by_scenario <- bind_rows(catch_lost_pinniped, catch_lost_pink, catch_lost_sst)

catch_lost_totals <- catch_lost_by_scenario %>%
  group_by(scenario, Year) %>%
  summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
  arrange(scenario, Year) %>%
  group_by(scenario) %>%
  mutate(cum_catch_lost = cumsum(catch_lost)) %>%
  ungroup()

# ============================================================
# PLOTS
# ============================================================

plot_productivity_compare <- function(model_actual, model_scenario, scenario_label, stocks = NULL) {
  df <- model_actual %>%
    select(Stock, Year, Actual = retro_lnRS) %>%
    left_join(model_scenario %>% select(Stock, Year, Scenario = retro_lnRS),
              by = c("Stock", "Year")) %>%
    pivot_longer(c(Actual, Scenario), names_to = "series", values_to = "lnRS") %>%
    mutate(series = if_else(series == "Scenario", scenario_label, "Actual"))
  
  if (!is.null(stocks)) df <- df %>% filter(Stock %in% stocks)
  
  ggplot(df, aes(Year, lnRS, color = series, linetype = series)) +
    geom_line(linewidth = 1, alpha = 0.6) +
    facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
    labs(x = "Year", y = "ln(R/S)", color = "", linetype = "",
         title = paste("Productivity:", scenario_label, "vs. actual")) +
    scale_color_manual(values = c("#4682B4", "#FF4500")) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_catch_lost_cumulative <- function(catch_lost_df, scenario_label) {
  ggplot(catch_lost_df, aes(Year, cum_catch_lost)) +
    geom_area(fill = "#FF4500", alpha = 0.2) +
    geom_line(linewidth = 1, color = "#FF4500") +
    facet_wrap(~ Stock, scales = "free_y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Year", y = "Cumulative catch lost",
         title = paste("Catch lost to", scenario_label, "- by stock")) +
    theme_minimal()
}

plot_productivity_compare(model_actual_retro, model_pinniped_retro, "Pinniped controls")
plot_productivity_compare(model_actual_retro, model_pink_retro, "Wild pink")
plot_productivity_compare(model_actual_retro, model_sst_retro, "SST")

plot_catch_lost_cumulative(catch_lost_pinniped, "pinnipeds")
plot_catch_lost_cumulative(catch_lost_pink, "pink salmon")
plot_catch_lost_cumulative(catch_lost_sst, "SST")

ggplot(catch_lost_totals, aes(Year, cum_catch_lost, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost", color = "Scenario driver",
       title = "Cumulative catch lost across all stocks, by scenario driver") +
  theme_minimal()

ggsave("figures/catch_lost_by_scenario_driver.png", width = 10, height = 6, dpi = 600)

# ============================================================
# MSY UNDER HIGH VS. LOW PRODUCTIVITY -- by stock, then averaged across
# stocks, overlaid as reference lines on the catch-lost plots below.
#
# High productivity: 1950-1970 (same period convention used in the
# probabilistic retrospective's Smsy benchmark).
# Low productivity: each stock's own low-productivity window from
# low_periods (the SAME window Xmin is built from in the probabilistic
# retrospective) -- requires low_periods to already exist in this
# session (run the probabilistic script's low-period section first).
#
# MSY = Rmsy - Smsy, same formula used throughout this whole pipeline:
# Smsy = (a/b)*(0.5-0.07*a), Rmsy = Smsy*exp(a-b*Smsy). Computed PER
# YEAR using that year's real covariates (and, for Late Shuswap, that
# year's own cycle-matched a/b), then averaged across the years in the
# period -- same "real covariates per year, averaged output" principle
# used everywhere else in this pipeline.
# ============================================================

SMSY_PERIOD_START_YEAR <- 1950
SMSY_PERIOD_END_YEAR   <- 1970

compute_msy_by_stock <- function(stock_name, period_start, period_end) {
  
  terms <- get_top_model_terms(stock_name, obs %>% filter(Stock == stock_name))
  
  period_data <- obs %>%
    filter(Stock == stock_name, Year >= period_start, Year <= period_end) %>%
    arrange(Year)
  
  if (length(terms$sel_covs) > 0) {
    period_data <- period_data %>% filter(if_all(all_of(terms$sel_covs), ~ !is.na(.)))
  }
  if (nrow(period_data) == 0) return(NA_real_)
  
  if (stock_name == "Late Shuswap") {
    period_data <- period_data %>%
      mutate(cycle = Year %% 4) %>%
      left_join(terms$cycle_terms, by = "cycle")
    cov_term <- compute_cov_term(period_data, terms$sel_covs, terms$cov_coefs)
    a_vec <- period_data$ra + cov_term
    b_vec <- period_data$rb
  } else {
    cov_term <- compute_cov_term(period_data, terms$sel_covs, terms$cov_coefs)
    a_vec <- terms$ra + cov_term
    b_vec <- rep(terms$rb, nrow(period_data))
  }
  
  Smsy_vec <- (a_vec / b_vec) * (0.5 - 0.07 * a_vec)
  Rmsy_vec <- Smsy_vec * exp(a_vec - b_vec * Smsy_vec)
  MSY_vec  <- Rmsy_vec - Smsy_vec
  
  mean(MSY_vec, na.rm = TRUE)
}

msy_high_by_stock <- tibble(
  Stock = STOCKS,
  MSY = sapply(STOCKS, compute_msy_by_stock,
               period_start = SMSY_PERIOD_START_YEAR, period_end = SMSY_PERIOD_END_YEAR)
)

msy_low_by_stock <- tibble(
  Stock = STOCKS,
  MSY = sapply(STOCKS, function(s) {
    low_p <- low_periods %>% filter(Stock == s)
    if (nrow(low_p) == 0) return(NA_real_)
    compute_msy_by_stock(s, low_p$low_period_start, low_p$low_period_end)
  })
)

print(msy_high_by_stock, n = Inf)
print(msy_low_by_stock, n = Inf)

msy_high_avg <- mean(msy_high_by_stock$MSY, na.rm = TRUE)
msy_low_avg  <- mean(msy_low_by_stock$MSY, na.rm = TRUE)

cat("Average MSY across stocks -- High productivity (1950-1970):", msy_high_avg, "\n")
cat("Average MSY across stocks -- Low productivity (own low period):", msy_low_avg, "\n")

MSY_LINETYPES <- c("MSY (high productivity)" = "dashed", "MSY (low productivity)" = "dotted")

# ============================================================
# AVERAGE YEARLY CATCH LOST, 2000-present
# ============================================================

RECENT_YEARS <- 2000:max(catch_lost_totals$Year, na.rm = TRUE)
catch_lost_recent <- catch_lost_totals %>% filter(Year %in% RECENT_YEARS)

catch_lost_recent_summary <- catch_lost_recent %>%
  group_by(scenario) %>%
  summarise(
    mean_catch_lost = mean(catch_lost, na.rm = TRUE),
    se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p1 <- ggplot(catch_lost_recent_summary, aes(x = reorder(scenario, -mean_catch_lost),
                                            y = mean_catch_lost, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_catch_lost - se_catch_lost,
                    ymax = mean_catch_lost + se_catch_lost), width = 0.15) +
  geom_hline(aes(yintercept = msy_high_avg, linetype = "MSY (high productivity)"), color = "black", linewidth = 0.7) +
  geom_hline(aes(yintercept = msy_low_avg, linetype = "MSY (low productivity)"), color = "black", linewidth = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = NULL, values = MSY_LINETYPES) +
  labs(x = NULL, y = "Mean yearly catch lost", fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/mean_catch_lost_by_scenario_2000-presen_v2.png", width = 8, height = 5.5, dpi = 600)

p2 <- ggplot(catch_lost_recent, aes(x = reorder(scenario, catch_lost, FUN = median),
                                    y = catch_lost, fill = scenario)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  geom_hline(aes(yintercept = msy_high_avg, linetype = "MSY (high productivity)"), color = "black", linewidth = 0.7) +
  geom_hline(aes(yintercept = msy_low_avg, linetype = "MSY (low productivity)"), color = "black", linewidth = 0.7) +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = NULL, values = MSY_LINETYPES) +
  labs(x = NULL, y = "Yearly catch lost", fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/catch_lost_boxplot_by_scenario_2000-present_v2.png", width = 8, height = 5.5, dpi = 600)

p1 | p2

# p3 <- ggplot(catch_lost_recent, aes(x = catch_lost, fill = scenario)) +
#   geom_density(alpha = 0.5, color = NA) +
#   scale_x_continuous(labels = scales::comma) +
#   labs(x = "Yearly catch lost", y = "Density", fill = "Scenario driver") +
#   theme_minimal() +
#   theme(legend.position = "bottom")
#
# p1 | p2 | p3
#
# ggsave("figures/sockeye_catch_lost_alpha05.png", width = 11, height = 5.5, dpi = 600)

p2.2 <- ggplot(catch_lost_recent, aes(x = reorder(scenario, catch_lost, FUN = median),
                                      y = catch_lost, fill = scenario)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  geom_hline(aes(yintercept = msy_high_avg, linetype = "MSY (high productivity)"), color = "black", linewidth = 0.7) +
  geom_hline(aes(yintercept = msy_low_avg, linetype = "MSY (low productivity)"), color = "black", linewidth = 0.7) +
  scale_y_continuous(labels = scales::comma,
                     trans = scales::sqrt_trans(),
                     breaks = scales::breaks_pretty(n = 8)) +
  scale_linetype_manual(name = NULL, values = MSY_LINETYPES) +
  labs(x = NULL, y = "Yearly catch lost (sqrt transformed)", fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1)) 

#ggsave("figures/catch_lost_sqrt.png", width = 8, height = 5.5, dpi = 600)
# Note: 1 negative value (SST, 2010) excluded under sqrt transform

p2 | p2.2

ggsave("figures/catch_lost_sqrt2_v2.png", width = 11.5, height = 5.5, dpi = 600)

# ============================================================
# RETURN TRAJECTORIES: observed vs. all scenarios
# ============================================================

SCENARIO_COLORS <- c(
  "Observed"          = "black",
  "Pinniped scenario" = "#4682B4",
  "Pink scenario"     = "#FF4500",
  "SST scenario"      = "#2E8B57"
)

returns_total <- bind_rows(
  model_actual_retro %>% group_by(Year) %>%
    summarise(Return = sum(RunJacks, na.rm = TRUE), .groups = "drop") %>% mutate(series = "Observed"),
  model_pinniped_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>% mutate(series = "Pinniped scenario"),
  model_pink_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>% mutate(series = "Pink scenario"),
  model_sst_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>% mutate(series = "SST scenario")
) %>%
  mutate(series = factor(series, levels = names(SCENARIO_COLORS)))

ggplot() +
  geom_line(data = filter(returns_total, series == "Observed"),
            aes(Year, Return, color = series), linewidth = 1.3) +
  geom_line(data = filter(returns_total, series != "Observed"),
            aes(Year, Return, color = series), linewidth = 1, alpha = 0.6) +
  scale_color_manual(values = SCENARIO_COLORS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Total return (all stocks)", color = "") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/return_trajectories_all_scenarios_total.png", width = 10, height = 6, dpi = 600)

returns_by_stock <- bind_rows(
  model_actual_retro %>% select(Stock, Year, Return = RunJacks) %>% mutate(series = "Observed"),
  model_pinniped_retro %>% select(Stock, Year, Return = retroR) %>% mutate(series = "Pinniped scenario"),
  model_pink_retro %>% select(Stock, Year, Return = retroR) %>% mutate(series = "Pink scenario"),
  model_sst_retro %>% select(Stock, Year, Return = retroR) %>% mutate(series = "SST scenario")
) %>%
  mutate(series = factor(series, levels = names(SCENARIO_COLORS)))

ggplot() +
  geom_line(data = filter(returns_by_stock, series == "Observed"),
            aes(Year, Return, color = series), linewidth = 1) +
  geom_line(data = filter(returns_by_stock, series != "Observed"),
            aes(Year, Return, color = series), linewidth = 0.7, alpha = 0.7) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Return", color = "",
       title = "Observed vs. scenario-reconstructed returns, by stock") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/return_trajectories_all_scenarios_by_stock.png", width = 14, height = 10, dpi = 600)

# ============================================================
# COSEWIC-STYLE STATUS CLASSIFICATION
#
# Approximate re-application of COSEWIC's quantitative thresholds
# (Criteria A/C/D) to model output -- not the full COSEWIC process
# (which also uses Criteria B/E, qualitative evidence, expert judgement).
# Treat status labels as comparative indicators across scenarios, not a
# formal reassessment.
# ============================================================

# Generation time (years), per COSEWIC (2017) Technical Summaries.
# All stocks are 4 years except Pitt (5).
GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5

# Right-aligned trailing mean, requires a full window of non-NA values --
# approximates COSEWIC's generation-smoothed abundance, and is reused to
# find each stock's low-abundance period below.
trailing_mean <- function(x, width) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (i >= width) {
      window <- x[(i - width + 1):i]
      if (all(!is.na(window))) out[i] <- mean(window)
    }
  }
  out
}

# Most severe status implied by any criterion. Special Concern isn't
# distinguished from Not at Risk (it's qualitative, not threshold-based,
# in the source report).
classify_status <- function(decline_pct, current_abundance) {
  case_when(
    is.na(decline_pct) & is.na(current_abundance) ~ NA_character_,
    (!is.na(decline_pct) & decline_pct >= 0.50) |
      (!is.na(current_abundance) & current_abundance < 2500)  ~ "Endangered",
    (!is.na(decline_pct) & decline_pct >= 0.30) |
      (!is.na(current_abundance) & current_abundance < 10000) ~ "Threatened",
    TRUE ~ "Not at Risk / Special Concern"
  )
}

compute_status_from_returns <- function(returns_df, scenario_label) {
  returns_df %>%
    left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
    arrange(Stock, Year) %>%
    group_by(Stock) %>%
    mutate(
      gen_mean          = trailing_mean(Return, width = first(GT)),
      gen_mean_3gen_ago = lag(gen_mean, n = 3 * first(GT)),
      decline_pct       = 1 - gen_mean / gen_mean_3gen_ago
    ) %>%
    ungroup() %>%
    mutate(status = classify_status(decline_pct, gen_mean), scenario = scenario_label)
}

status_actual <- compute_status_from_returns(
  model_actual_retro %>% select(Stock, Year, Return = RunJacks), "Observed"
)
status_pinniped <- compute_status_from_returns(
  model_pinniped_retro %>% select(Stock, Year, Return = retroR), "Pinniped scenario"
)
status_sst <- compute_status_from_returns(
  model_sst_retro %>% select(Stock, Year, Return = retroR), "SST scenario"
)
status_pink <- compute_status_from_returns(
  model_pink_retro %>% select(Stock, Year, Return = retroR), "Pink scenario"
)

status_all <- bind_rows(status_actual, status_pinniped, status_sst, status_pink)

# All series capped at MAX_ASSESSMENT_YEAR for a fair, contemporaneous
# comparison -- otherwise Observed would be frozen at its last data year
# while scenarios run further, which isn't like-for-like. This answers
# "given what was known as of the 2017 report, would this driver have
# changed the designation." status_summary_latest_uncapped below instead
# answers "what would status be today under this scenario."
MAX_ASSESSMENT_YEAR <- 2016

status_summary_latest <- status_all %>%
  filter(Year <= MAX_ASSESSMENT_YEAR, !is.na(status)) %>%
  group_by(Stock, scenario) %>%
  slice_max(Year, n = 1) %>%
  ungroup() %>%
  select(Stock, scenario, Year, gen_mean, decline_pct, status)

status_summary_latest_uncapped <- status_all %>%
  filter(!is.na(status)) %>%
  group_by(Stock, scenario) %>%
  slice_max(Year, n = 1) %>%
  ungroup() %>%
  select(Stock, scenario, Year, gen_mean, decline_pct, status)

# Stocks COSEWIC (2017) actually designated Endangered -- used instead of
# our own classify_status() output, since our thresholds don't
# distinguish Special Concern and would otherwise pull in extra stocks.
COSEWIC_ENDANGERED_STOCKS <- c("Bowron", "Weaver", "Quesnel", "Early Stuart",
                               "Late Stuart", "Portage", "Cultus")

status_comparison_table <- status_summary_latest %>%
  select(Stock, scenario, status) %>%
  pivot_wider(names_from = scenario, values_from = status) %>%
  filter(Stock %in% COSEWIC_ENDANGERED_STOCKS) %>%
  arrange(Stock)

print(status_comparison_table)

# ============================================================
# % ABUNDANCE INCREASE OVER LOW-ABUNDANCE PERIODS, BY SCENARIO
# ============================================================

returns_by_stock_scenario <- bind_rows(
  model_actual_retro   %>% select(Stock, Year, Return = RunJacks) %>% mutate(scenario = "Observed"),
  model_pinniped_retro %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pinniped scenario"),
  model_sst_retro      %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "SST scenario"),
  model_pink_retro     %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pink scenario")
)

# Each stock's low-abundance PERIOD: the GT-year trailing window (in the
# observed series) with the lowest mean return, matching the same
# generation-based smoothing used in status classification above.
low_periods <- returns_by_stock_scenario %>%
  filter(scenario == "Observed", is.finite(Return)) %>%
  left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(obs_gen_mean = trailing_mean(Return, width = first(GT))) %>%
  filter(!is.na(obs_gen_mean)) %>%
  slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, GT, low_period_end = Year, low_period_start = Year - GT + 1,
            low_period_mean = obs_gen_mean)

pct_increase_by_year <- returns_by_stock_scenario %>%
  filter(scenario != "Observed") %>%
  select(Stock, Year, scenario, Return_scenario = Return) %>%
  left_join(
    returns_by_stock_scenario %>% filter(scenario == "Observed") %>%
      select(Stock, Year, Return_observed = Return),
    by = c("Stock", "Year")
  ) %>%
  mutate(pct_increase = 100 * (Return_scenario - Return_observed) / Return_observed)

# % increase in mean abundance over each stock's low-abundance period --
# how much better off the stock would have been under each counterfactual
pct_increase_over_low_period <- returns_by_stock_scenario %>%
  filter(scenario != "Observed") %>%
  inner_join(low_periods, by = "Stock") %>%
  filter(Year >= low_period_start, Year <= low_period_end) %>%
  group_by(Stock, scenario, low_period_start, low_period_end, low_period_mean) %>%
  summarise(scenario_period_mean = mean(Return, na.rm = TRUE), .groups = "drop") %>%
  mutate(pct_increase = 100 * (scenario_period_mean - low_period_mean) / low_period_mean) %>%
  arrange(Stock, scenario)

print(pct_increase_over_low_period)

recovery_summary <- status_comparison_table %>%
  left_join(
    pct_increase_over_low_period %>% select(Stock, scenario, pct_increase) %>%
      pivot_wider(names_from = scenario, values_from = pct_increase, names_prefix = "pct_increase_"),
    by = "Stock"
  )

print(recovery_summary)

# ============================================================
# PLOTS
# ============================================================

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

ggplot(pct_increase_over_low_period,
       aes(x = reorder(Stock, -pct_increase), y = pct_increase, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  labs(x = NULL, y = "% increase in mean abundance over low-abundance period",
       fill = "Scenario driver",
       title = "Potential abundance recovery over each stock's historical low-abundance period") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/pct_increase_over_low_period_by_scenario.png", width = 10, height = 6, dpi = 600)

ggplot(pct_increase_by_year, aes(Year, pct_increase, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.25) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS2) +
  labs(x = "Year", y = "% increase in abundance vs. observed", color = "Scenario driver",
       title = "Counterfactual abundance gain over time, by stock",
       subtitle = "Shaded band marks each stock's historical low-abundance period (worst generation-length window)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/pct_increase_over_time_by_scenario.png", width = 14, height = 10, dpi = 600)

status_plot_df <- status_summary_latest %>%
  filter(Stock %in% status_comparison_table$Stock) %>%
  mutate(
    scenario = factor(scenario, levels = c("Observed", "Pinniped scenario", "SST scenario", "Pink scenario")),
    status = factor(status, levels = c("Endangered", "Threatened", "Not at Risk / Special Concern"))
  )

ggplot(status_plot_df, aes(scenario, Stock, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "Endangered" = "#B22222",
    "Threatened" = "#E8A33D",
    "Not at Risk / Special Concern" = "#4C9A5B"
  )) +
  labs(x = NULL, y = NULL, fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("figures/status_by_scenario_heatmap.png", width = 9, height = 6, dpi = 600)

# ============================================================
# DIAGNOSTICS
#
# diagnose_stock: full spawners/return/productivity trajectory for one
# stock across scenarios, plus how far outside its historical covariate
# range the pinniped freeze-year value sits (extrapolation risk). Useful
# when a stock's classification worsens under a scenario despite a
# favorable covariate coefficient -- the retrospective model is
# recursive (spawners -> recruits -> spawners...), so a productivity
# boost raises escapement, which feeds the density-dependent term next
# generation; a stock with large |rb| relative to ra can overshoot and
# crash on the next cycle (Ricker overcompensation) rather than
# improving monotonically. Also worth checking whether
# status_summary_latest's slice_max(Year) happens to land on the
# down-swing of such a cycle.
#
# inspect_status_trajectory: year-by-year gen_mean/decline_pct/status for
# one stock/scenario, to see exactly which year the classification is
# evaluated at.
# ============================================================

diagnose_stock <- function(stock_name) {
  
  traj <- bind_rows(
    model_actual_retro   %>% filter(Stock == stock_name) %>% mutate(scenario = "Observed"),
    model_pinniped_retro %>% filter(Stock == stock_name) %>% mutate(scenario = "Pinniped scenario"),
    model_sst_retro      %>% filter(Stock == stock_name) %>% mutate(scenario = "SST scenario"),
    model_pink_retro     %>% filter(Stock == stock_name) %>% mutate(scenario = "Pink scenario")
  )
  
  p1 <- ggplot(traj, aes(Year, retroS, color = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Observed" = "black", SCENARIO_COLORS2)) +
    labs(y = "Spawners (retroS)", title = paste(stock_name, "- spawners by scenario")) +
    theme_minimal()
  
  p2 <- ggplot(traj, aes(Year, retroR, color = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Observed" = "black", SCENARIO_COLORS2)) +
    labs(y = "Return (retroR)", title = paste(stock_name, "- returns by scenario")) +
    theme_minimal()
  
  p3 <- ggplot(traj, aes(Year, retro_lnRS, color = scenario)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Observed" = "black", SCENARIO_COLORS2)) +
    labs(y = "ln(R/S)", title = paste(stock_name, "- productivity by scenario")) +
    theme_minimal()
  
  print((p1 / p2 / p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom"))
  
  ggsave(paste0("figures/diagnostic_", tolower(gsub(" ", "_", stock_name)), "_trajectories.png"),
         width = 9, height = 10, dpi = 600)
  
  hist_range <- covariates_main %>%
    filter(Stock == stock_name) %>%
    summarise(
      SeaLions_hist_min  = min(SeaLions, na.rm = TRUE),
      SeaLions_hist_max  = max(SeaLions, na.rm = TRUE),
      SeaLions_at_freeze = SeaLions[Year == FREEZE_YEAR],
      seal_hist_min      = min(seal, na.rm = TRUE),
      seal_hist_max      = max(seal, na.rm = TRUE),
      seal_at_freeze     = seal[Year == FREEZE_YEAR]
    )
  print(hist_range)
  
  invisible(traj)
}

diagnose_stock("Raft")

inspect_status_trajectory <- function(stock_name, scenario_label = "Observed") {
  status_all %>%
    filter(Stock == stock_name, scenario == scenario_label) %>%
    select(Year, gen_mean, decline_pct, status) %>%
    print(n = Inf)
}

inspect_status_trajectory("Late Stuart")