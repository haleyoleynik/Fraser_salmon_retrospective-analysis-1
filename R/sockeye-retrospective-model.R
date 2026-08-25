# Sockeye retrospective model — scenario comparisons
# Adapted from Carl Walters Excel model
# Haley Oleynik
#
# Cleaned version: keeps only the four covariate-driven retrospective
# scenarios (pinniped, pink vs. pink_wild, SST, combined) and runs each
# through the same pipeline so catch-lost estimates can be compared
# directly across drivers at the end.

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# ============================================================
# CONFIG
# ============================================================

FIT_YEARS          <- 1952:2019
LAG_YEARS          <- 4
COVARIATE_COLS      <- c("NPGO", "PDO", "SeaLions", "seal", "adult.sst", "pink", "smolt.sst")

FREEZE_YEAR         <- 1970        # pinniped covariates frozen at this year's level, forward
SST_BASELINE_YEARS  <- 1950:1975   # years used for stock-specific SST long-term mean

retroU_default   <- 0.3   # retrospective harvest rate cap (used for *_retro model runs)
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

# Main standardized covariates (NPGO, PDO, SeaLions, seal, adult.sst, pink, smolt.sst)
covariates_main <- read_csv("Data/sockeye_standardized_covariates.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, any_of(COVARIATE_COLS))

# Wild-pink alternative series (pink abundance had hatchery production not occurred)
covariates_pink_wild <- read_csv("Data/sockeye_standardized_covariates_pink-wild.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, pink_wild)

# Per-stock AIC-weighted top-model coefficients (from dredge/model-selection output)
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
# Each scenario gets its own "<var>_scenario" column(s), all joined onto
# a single `obs` table so any driver can be swapped in at model-run time.
# ============================================================

## Pinniped scenario: SeaLions/seal frozen at FREEZE_YEAR levels, forward
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

## Pink scenario: pink substituted with the wild-only pink series
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

## SST scenario: adult.sst / smolt.sst frozen at each stock's 1950-1975 mean
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

## Assemble obs with actual covariates + every scenario column attached.
## "Combined" reuses the pinniped + pink scenario columns already here --
## it doesn't need a table of its own.
obs <- obs_raw %>%
  left_join(covariates_main,        by = c("Stock", "Year")) %>%
  left_join(pinniped_scenario_cov,  by = c("Stock", "Year")) %>%
  left_join(pink_scenario_cov,      by = c("Stock", "Year")) %>%
  left_join(sst_scenario_cov,       by = c("Stock", "Year"))

# ============================================================
# MODEL FUNCTIONS
# ============================================================

# Look up a stock's AIC-weighted top-model coefficients; fall back to a
# plain Ricker fit if the stock has no dredge top model.
get_top_model_terms <- function(stock_name, fit_df_for_fallback) {
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
# scenario_vars controls which covariates get swapped to their scenario
# column *in the forward projection only* -- the process-error residuals
# (wt) are always estimated against the actual historical covariates, so
# a scenario comparison isolates the covariate effect rather than mixing
# it with a different noise realization.
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
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term)
    )
  
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
    retro_lnRS[j] <- terms$ra - terms$rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(retroR = retroR, retro_lnRS = retro_lnRS, retroU = retroU_vec,
           retroS = retroS, retroC = retroC)
}

# Run run_retro_model() across all stocks for a given harvest-rate
# assumption + scenario driver.
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
# Two harvest-rate assumptions per scenario:
#   *_retro : capped retrospective harvest rate (for productivity / return / spawner comparisons)
#   *_hist  : actual historical harvest rate (for "catch lost to driver" comparisons)
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

model_combined_retro <- run_scenario_model(obs, retroU_default, useretro_default, yrretro_default,
                                           scenario_vars = c("SeaLions", "seal", "pink"))
model_combined_hist  <- run_scenario_model(obs, retroU = 0, useretro = FALSE, yrretro = yrretro_default,
                                           scenario_vars = c("SeaLions", "seal", "pink"))

# ============================================================
# CATCH LOST TO EACH DRIVER
# catch_lost = scenario catch - actual catch, at the actual historical
# harvest rate (positive = catch that would have been available if the
# driver's covariate(s) had stayed at the scenario level)
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
catch_lost_combined <- compute_catch_lost(model_combined_hist, "Combined (pinniped + pink)")

# All scenarios together, for direct comparison
catch_lost_by_scenario <- bind_rows(
  catch_lost_pinniped, catch_lost_pink, catch_lost_sst, catch_lost_combined
)

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
plot_productivity_compare(model_actual_retro, model_combined_retro, "Combined")

plot_catch_lost_cumulative(catch_lost_pinniped, "pinnipeds")
plot_catch_lost_cumulative(catch_lost_pink, "pink salmon")
plot_catch_lost_cumulative(catch_lost_sst, "SST")
plot_catch_lost_cumulative(catch_lost_combined, "combined drivers")

# Cross-scenario comparison: total cumulative catch lost, all drivers overlaid
ggplot(catch_lost_totals, aes(Year, cum_catch_lost, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost", color = "Scenario driver",
       title = "Cumulative catch lost across all stocks, by scenario driver") +
  theme_minimal()

ggsave("figures/catch_lost_by_scenario_driver.png", width = 10, height = 6, dpi = 600)


# ============================================================
# COMPARE SCENARIOS: average yearly catch lost, 2000-present
# Uses catch_lost_totals (total catch lost across all stocks, per year,
# per scenario) from sockeye_retrospective_scenarios.R
# ============================================================

RECENT_YEARS <- 2000:max(catch_lost_totals$Year, na.rm = TRUE)

catch_lost_recent <- catch_lost_totals %>%
  filter(Year %in% RECENT_YEARS)

## Bar chart: mean yearly catch lost per scenario, 2000-present -----------

catch_lost_recent_summary <- catch_lost_recent %>%
  group_by(scenario) %>%
  summarise(
    mean_catch_lost = mean(catch_lost, na.rm = TRUE),
    se_catch_lost    = sd(catch_lost, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

catch_lost_recent_summary %>%
  filter(scenario != "Combined (pinniped + pink)") %>%
ggplot(aes(x = reorder(scenario, -mean_catch_lost),
                                      y = mean_catch_lost, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_catch_lost - se_catch_lost,
                    ymax = mean_catch_lost + se_catch_lost),
                width = 0.15) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Mean yearly catch lost",
       fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/mean_catch_lost_by_scenario_2000-present.png", width = 8, height = 5.5, dpi = 600)

## Boxplot: yearly distribution of catch lost per scenario, 2000-present --

catch_lost_recent %>%
  filter(scenario != "Combined (pinniped + pink)") %>%
ggplot(aes(x = reorder(scenario, catch_lost, FUN = median),
                              y = catch_lost, fill = scenario)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Yearly catch lost",
       fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/catch_lost_boxplot_by_scenario_2000-present.png", width = 8, height = 5.5, dpi = 600)


# ============================================================
# COMPARE SCENARIOS: return trajectories, observed vs. all scenarios
# Uses model_actual_retro / model_pinniped_retro / model_pink_retro /
# model_sst_retro / model_combined_retro from
# sockeye_retrospective_scenarios.R (all run under the same retrospective
# harvest-rate assumption, so lines differ only by covariate driver)
# ============================================================

SCENARIO_COLORS <- c(
  "Observed"            = "black",
  "Pinniped scenario"   = "#4682B4",
  "Pink scenario" = "#FF4500",
  "SST scenario"         = "#2E8B57")

## --- Total return trajectory (summed across all stocks) -----------------

returns_total <- bind_rows(
  model_actual_retro %>% group_by(Year) %>%
    summarise(Return = sum(RunJacks, na.rm = TRUE), .groups = "drop") %>%
    mutate(series = "Observed"),
  model_pinniped_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>%
    mutate(series = "Pinniped scenario"),
  model_pink_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>%
    mutate(series = "Pink scenario"),
  model_sst_retro %>% group_by(Year) %>%
    summarise(Return = sum(retroR, na.rm = TRUE), .groups = "drop") %>%
    mutate(series = "SST scenario")
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

## --- Same comparison, faceted by stock -----------------------------------

returns_by_stock <- bind_rows(
  model_actual_retro %>% select(Stock, Year, Return = RunJacks) %>%
    mutate(series = "Observed"),
  model_pinniped_retro %>% select(Stock, Year, Return = retroR) %>%
    mutate(series = "Pinniped scenario"),
  model_pink_retro %>% select(Stock, Year, Return = retroR) %>%
    mutate(series = "Pink (wild) scenario"),
  model_sst_retro %>% select(Stock, Year, Return = retroR) %>%
    mutate(series = "SST scenario"),
  model_combined_retro %>% select(Stock, Year, Return = retroR) %>%
    mutate(series = "Combined scenario")
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
# Sockeye retrospective model — recovery status & low-period analysis
# Appends to: sockeye retrospective scenario script (pinniped / SST /
# pink models already run: model_actual_retro, model_pinniped_retro,
# model_sst_retro, model_pink_retro, plus catch_lost_totals etc.)
#
# Adds:
#   1. Approximate COSEWIC status classification (Criteria A/C/D) for
#      each stock under the Observed / Pinniped / SST / Pink scenarios,
#      so "would this stock still be Endangered" can be compared across
#      drivers.
#   2. %% increase in abundance, over each stock's historical
#      low-abundance period (worst generation-length window, not a
#      single year), that each scenario implies relative to what
#      actually happened -- plus the same comparison over the full
#      time series with that low period shaded.
#
# CAVEAT: This is a simplified re-application of COSEWIC's quantitative
# thresholds (Criteria A, C, D) to model output. It does not replicate
# the full COSEWIC process, which also incorporates Criteria B/E,
# qualitative lines of evidence, and expert judgement. Treat status
# labels here as approximate, comparative indicators across scenarios
# -- not as a formal reassessment.
# ============================================================

# ------------------------------------------------------------
# 1. COSEWIC-STYLE STATUS CLASSIFICATION
# ------------------------------------------------------------

# Generation time (years) per stock, per COSEWIC (2017) Technical
# Summaries. All stocks are 4 years except Pitt (5).
GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5

# Trailing mean over `width` years (right-aligned, requires a full
# window of non-NA values) -- used to approximate COSEWIC's
# generation-smoothed "number of mature individuals", and to find
# low-abundance periods below.
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

# Classify a stock-year as Endangered / Threatened / "Not at Risk or
# Special Concern" using the most severe status implied by any
# criterion (matching how COSEWIC applies e.g. "meets Endangered under
# A" even when C/D alone would only imply Threatened). Special Concern
# is not distinguished from Not at Risk here, since it is a
# qualitative (not threshold-based) designation in the source report.
classify_status <- function(decline_pct, current_abundance) {
  dplyr::case_when(
    is.na(decline_pct) & is.na(current_abundance) ~ NA_character_,
    (!is.na(decline_pct) & decline_pct >= 0.50) |
      (!is.na(current_abundance) & current_abundance < 2500)  ~ "Endangered",
    (!is.na(decline_pct) & decline_pct >= 0.30) |
      (!is.na(current_abundance) & current_abundance < 10000) ~ "Threatened",
    TRUE ~ "Not at Risk / Special Concern"
  )
}

# Compute generation-smoothed abundance, 3-generation decline, and
# resulting status for a Stock/Year/Return series under one scenario.
compute_status_from_returns <- function(returns_df, scenario_label) {
  returns_df %>%
    dplyr::left_join(
      tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"),
      by = "Stock"
    ) %>%
    dplyr::arrange(Stock, Year) %>%
    dplyr::group_by(Stock) %>%
    dplyr::mutate(
      gen_mean          = trailing_mean(Return, width = dplyr::first(GT)),
      gen_mean_3gen_ago = dplyr::lag(gen_mean, n = 3 * dplyr::first(GT)),
      decline_pct       = 1 - gen_mean / gen_mean_3gen_ago
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      status   = classify_status(decline_pct, gen_mean),
      scenario = scenario_label
    )
}

# "Observed" uses raw RunJacks (as the rest of the script does for the
# actual series); scenario models use their reconstructed retroR.
status_actual <- compute_status_from_returns(
  model_actual_retro %>% dplyr::select(Stock, Year, Return = RunJacks), "Observed"
)
status_pinniped <- compute_status_from_returns(
  model_pinniped_retro %>% dplyr::select(Stock, Year, Return = retroR), "Pinniped scenario"
)
status_sst <- compute_status_from_returns(
  model_sst_retro %>% dplyr::select(Stock, Year, Return = retroR), "SST scenario"
)
status_pink <- compute_status_from_returns(
  model_pink_retro %>% dplyr::select(Stock, Year, Return = retroR), "Pink scenario"
)

status_all <- dplyr::bind_rows(status_actual, status_pinniped, status_sst, status_pink)

# Most recent classifiable year, per stock and scenario
status_summary_latest <- status_all %>%
  dplyr::filter(!is.na(status)) %>%
  dplyr::group_by(Stock, scenario) %>%
  dplyr::slice_max(Year, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Stock, scenario, Year, gen_mean, decline_pct, status)

# Wide comparison table: observed status vs. each scenario's status,
# restricted to stocks currently Endangered or Threatened (i.e. the
# ones the recovery question is actually about)
status_comparison_table <- status_summary_latest %>%
  dplyr::select(Stock, scenario, status) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = status) %>%
  dplyr::filter(Observed %in% c("Endangered", "Threatened")) %>%
  dplyr::arrange(Stock)

print(status_comparison_table)

# ------------------------------------------------------------
# 2. %% ABUNDANCE INCREASE OVER LOW-ABUNDANCE PERIODS, BY SCENARIO
# ------------------------------------------------------------

returns_by_stock_scenario <- dplyr::bind_rows(
  model_actual_retro   %>% dplyr::select(Stock, Year, Return = RunJacks) %>% dplyr::mutate(scenario = "Observed"),
  model_pinniped_retro %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "Pinniped scenario"),
  model_sst_retro      %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "SST scenario"),
  model_pink_retro      %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "Pink scenario")
)

# Each stock's low-abundance PERIOD: the GT-year trailing window (in
# the *observed* series) with the lowest mean return, rather than a
# single low year. GT-length window matches the same generation-based
# smoothing used in the status classification above.
low_periods <- returns_by_stock_scenario %>%
  dplyr::filter(scenario == "Observed", is.finite(Return)) %>%
  dplyr::left_join(
    tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"),
    by = "Stock"
  ) %>%
  dplyr::arrange(Stock, Year) %>%
  dplyr::group_by(Stock) %>%
  dplyr::mutate(obs_gen_mean = trailing_mean(Return, width = dplyr::first(GT))) %>%
  dplyr::filter(!is.na(obs_gen_mean)) %>%
  dplyr::slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    Stock,
    GT,
    low_period_end   = Year,
    low_period_start = Year - GT + 1,
    low_period_mean  = obs_gen_mean
  )

# %% increase, scenario vs. observed, at every year (so the low period
# can be read in context of the whole trajectory)
pct_increase_by_year <- returns_by_stock_scenario %>%
  dplyr::filter(scenario != "Observed") %>%
  dplyr::select(Stock, Year, scenario, Return_scenario = Return) %>%
  dplyr::left_join(
    returns_by_stock_scenario %>%
      dplyr::filter(scenario == "Observed") %>%
      dplyr::select(Stock, Year, Return_observed = Return),
    by = c("Stock", "Year")
  ) %>%
  dplyr::mutate(pct_increase = 100 * (Return_scenario - Return_observed) / Return_observed)

# %% increase in mean abundance over each stock's low-abundance PERIOD
# -- i.e. how much better off the stock would have been, across its
# worst generation-length stretch, under each counterfactual driver
pct_increase_over_low_period <- returns_by_stock_scenario %>%
  dplyr::filter(scenario != "Observed") %>%
  dplyr::inner_join(low_periods, by = "Stock") %>%
  dplyr::filter(Year >= low_period_start, Year <= low_period_end) %>%
  dplyr::group_by(Stock, scenario, low_period_start, low_period_end, low_period_mean) %>%
  dplyr::summarise(scenario_period_mean = mean(Return, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pct_increase = 100 * (scenario_period_mean - low_period_mean) / low_period_mean) %>%
  dplyr::arrange(Stock, scenario)

print(pct_increase_over_low_period)

# Combined status + low-period recovery summary, one row per at-risk stock
recovery_summary <- status_comparison_table %>%
  dplyr::left_join(
    pct_increase_over_low_period %>%
      dplyr::select(Stock, scenario, pct_increase) %>%
      tidyr::pivot_wider(names_from = scenario, values_from = pct_increase,
                         names_prefix = "pct_increase_"),
    by = "Stock"
  )

print(recovery_summary)

# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

## Bar chart: % increase in mean abundance over each stock's
## low-abundance period, by scenario
ggplot(pct_increase_over_low_period,
       aes(x = reorder(Stock, -pct_increase), y = pct_increase, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  labs(x = NULL, y = "% increase in mean abundance over low-abundance period",
       fill = "Scenario driver") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/pct_increase_over_low_period_by_scenario.png", width = 10, height = 6, dpi = 600)

## Time series: % increase over time, faceted by stock, low-abundance
## period shaded
ggplot(pct_increase_by_year, aes(Year, pct_increase, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.25) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS2) +
  labs(x = "Year", y = "% increase in abundance vs. observed", color = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/pct_increase_over_time_by_scenario.png", width = 14, height = 10, dpi = 600)

## Status-comparison heatmap: does the stock stay Endangered/Threatened
## under each scenario, or does it drop out of at-risk status?
status_plot_df <- status_summary_latest %>%
  dplyr::filter(Stock %in% status_comparison_table$Stock) %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = c("Observed", "Pinniped scenario", "SST scenario", "Pink scenario")),
    status = factor(status, levels = c("Endangered", "Threatened",
                                       "Not at Risk / Special Concern"))
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


# do again with only stocks that are actually identified as thretened: 


# ============================================================
# Sockeye retrospective model — recovery status & low-period analysis
# Appends to: sockeye retrospective scenario script (pinniped / SST /
# pink models already run: model_actual_retro, model_pinniped_retro,
# model_sst_retro, model_pink_retro, plus catch_lost_totals etc.)
#
# Adds:
#   1. Approximate COSEWIC status classification (Criteria A/C/D) for
#      each stock under the Observed / Pinniped / SST / Pink scenarios,
#      so "would this stock still be Endangered" can be compared across
#      drivers.
#   2. %% increase in abundance, over each stock's historical
#      low-abundance period (worst generation-length window, not a
#      single year), that each scenario implies relative to what
#      actually happened -- plus the same comparison over the full
#      time series with that low period shaded.
#
# CAVEAT: This is a simplified re-application of COSEWIC's quantitative
# thresholds (Criteria A, C, D) to model output. It does not replicate
# the full COSEWIC process, which also incorporates Criteria B/E,
# qualitative lines of evidence, and expert judgement. Treat status
# labels here as approximate, comparative indicators across scenarios
# -- not as a formal reassessment.
# ============================================================

# ------------------------------------------------------------
# 1. COSEWIC-STYLE STATUS CLASSIFICATION
# ------------------------------------------------------------

# Generation time (years) per stock, per COSEWIC (2017) Technical
# Summaries. All stocks are 4 years except Pitt (5).
GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5

# Trailing mean over `width` years (right-aligned, requires a full
# window of non-NA values) -- used to approximate COSEWIC's
# generation-smoothed "number of mature individuals", and to find
# low-abundance periods below.
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

# Classify a stock-year as Endangered / Threatened / "Not at Risk or
# Special Concern" using the most severe status implied by any
# criterion (matching how COSEWIC applies e.g. "meets Endangered under
# A" even when C/D alone would only imply Threatened). Special Concern
# is not distinguished from Not at Risk here, since it is a
# qualitative (not threshold-based) designation in the source report.
classify_status <- function(decline_pct, current_abundance) {
  dplyr::case_when(
    is.na(decline_pct) & is.na(current_abundance) ~ NA_character_,
    (!is.na(decline_pct) & decline_pct >= 0.50) |
      (!is.na(current_abundance) & current_abundance < 2500)  ~ "Endangered",
    (!is.na(decline_pct) & decline_pct >= 0.30) |
      (!is.na(current_abundance) & current_abundance < 10000) ~ "Threatened",
    TRUE ~ "Not at Risk / Special Concern"
  )
}

# Compute generation-smoothed abundance, 3-generation decline, and
# resulting status for a Stock/Year/Return series under one scenario.
compute_status_from_returns <- function(returns_df, scenario_label) {
  returns_df %>%
    dplyr::left_join(
      tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"),
      by = "Stock"
    ) %>%
    dplyr::arrange(Stock, Year) %>%
    dplyr::group_by(Stock) %>%
    dplyr::mutate(
      gen_mean          = trailing_mean(Return, width = dplyr::first(GT)),
      gen_mean_3gen_ago = dplyr::lag(gen_mean, n = 3 * dplyr::first(GT)),
      decline_pct       = 1 - gen_mean / gen_mean_3gen_ago
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      status   = classify_status(decline_pct, gen_mean),
      scenario = scenario_label
    )
}

# "Observed" uses raw RunJacks (as the rest of the script does for the
# actual series); scenario models use their reconstructed retroR.
status_actual <- compute_status_from_returns(
  model_actual_retro %>% dplyr::select(Stock, Year, Return = RunJacks), "Observed"
)
status_pinniped <- compute_status_from_returns(
  model_pinniped_retro %>% dplyr::select(Stock, Year, Return = retroR), "Pinniped scenario"
)
status_sst <- compute_status_from_returns(
  model_sst_retro %>% dplyr::select(Stock, Year, Return = retroR), "SST scenario"
)
status_pink <- compute_status_from_returns(
  model_pink_retro %>% dplyr::select(Stock, Year, Return = retroR), "Pink scenario"
)

status_all <- dplyr::bind_rows(status_actual, status_pinniped, status_sst, status_pink)

# Most recent classifiable year, per stock and scenario.
#
# For a fair, contemporaneous comparison, ALL series (Observed AND
# scenario) are capped at MAX_ASSESSMENT_YEAR here -- otherwise
# Observed would be frozen at 2016 while scenarios ran through
# whatever year your data happens to end at, which isn't a like-for-
# like comparison. This answers "given what was known as of the 2017
# report, would this driver have changed the designation." An
# uncapped version (status_summary_latest_uncapped, below) is also
# provided if you instead want "what would status be today under this
# scenario" -- a different, forward-looking question.
MAX_ASSESSMENT_YEAR <- 2016

status_all_capped <- status_all %>%
  dplyr::filter(Year <= MAX_ASSESSMENT_YEAR)

status_summary_latest <- status_all_capped %>%
  dplyr::filter(!is.na(status)) %>%
  dplyr::group_by(Stock, scenario) %>%
  dplyr::slice_max(Year, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Stock, scenario, Year, gen_mean, decline_pct, status)

# Uncapped alternative: each series evaluated at its own latest
# available year (scenarios can run past 2016; Observed also reflects
# any post-report years in your data). Not used in the tables/plots
# below by default -- swap status_summary_latest for this if you want
# "status today" instead of "status as of the 2017 report."
status_summary_latest_uncapped <- status_all %>%
  dplyr::filter(!is.na(status)) %>%
  dplyr::group_by(Stock, scenario) %>%
  dplyr::slice_max(Year, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(Stock, scenario, Year, gen_mean, decline_pct, status)

# Wide comparison table: observed status vs. each scenario's status,
# restricted to stocks currently Endangered or Threatened (i.e. the
# ones the recovery question is actually about)
# Stocks COSEWIC (2017) actually designated Endangered (not our own
# simplified classify_status() output -- see explanation in chat: our
# thresholds don't distinguish Special Concern, so filtering on our
# own "Observed" status pulls in a few Special Concern stocks too).
COSEWIC_ENDANGERED_STOCKS <- c("Bowron", "Weaver", "Quesnel", "Early Stuart",
                               "Late Stuart", "Portage", "Cultus")

status_comparison_table <- status_summary_latest %>%
  dplyr::select(Stock, scenario, status) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = status) %>%
  dplyr::filter(Stock %in% COSEWIC_ENDANGERED_STOCKS) %>%
  dplyr::arrange(Stock)

print(status_comparison_table)

# ------------------------------------------------------------
# 2. %% ABUNDANCE INCREASE OVER LOW-ABUNDANCE PERIODS, BY SCENARIO
# ------------------------------------------------------------

returns_by_stock_scenario <- dplyr::bind_rows(
  model_actual_retro   %>% dplyr::select(Stock, Year, Return = RunJacks) %>% dplyr::mutate(scenario = "Observed"),
  model_pinniped_retro %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "Pinniped scenario"),
  model_sst_retro      %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "SST scenario"),
  model_pink_retro      %>% dplyr::select(Stock, Year, Return = retroR)   %>% dplyr::mutate(scenario = "Pink scenario")
)

# Each stock's low-abundance PERIOD: the GT-year trailing window (in
# the *observed* series) with the lowest mean return, rather than a
# single low year. GT-length window matches the same generation-based
# smoothing used in the status classification above.
low_periods <- returns_by_stock_scenario %>%
  dplyr::filter(scenario == "Observed", is.finite(Return)) %>%
  dplyr::left_join(
    tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"),
    by = "Stock"
  ) %>%
  dplyr::arrange(Stock, Year) %>%
  dplyr::group_by(Stock) %>%
  dplyr::mutate(obs_gen_mean = trailing_mean(Return, width = dplyr::first(GT))) %>%
  dplyr::filter(!is.na(obs_gen_mean)) %>%
  dplyr::slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    Stock,
    GT,
    low_period_end   = Year,
    low_period_start = Year - GT + 1,
    low_period_mean  = obs_gen_mean
  )

# %% increase, scenario vs. observed, at every year (so the low period
# can be read in context of the whole trajectory)
pct_increase_by_year <- returns_by_stock_scenario %>%
  dplyr::filter(scenario != "Observed") %>%
  dplyr::select(Stock, Year, scenario, Return_scenario = Return) %>%
  dplyr::left_join(
    returns_by_stock_scenario %>%
      dplyr::filter(scenario == "Observed") %>%
      dplyr::select(Stock, Year, Return_observed = Return),
    by = c("Stock", "Year")
  ) %>%
  dplyr::mutate(pct_increase = 100 * (Return_scenario - Return_observed) / Return_observed)

# %% increase in mean abundance over each stock's low-abundance PERIOD
# -- i.e. how much better off the stock would have been, across its
# worst generation-length stretch, under each counterfactual driver
pct_increase_over_low_period <- returns_by_stock_scenario %>%
  dplyr::filter(scenario != "Observed") %>%
  dplyr::inner_join(low_periods, by = "Stock") %>%
  dplyr::filter(Year >= low_period_start, Year <= low_period_end) %>%
  dplyr::group_by(Stock, scenario, low_period_start, low_period_end, low_period_mean) %>%
  dplyr::summarise(scenario_period_mean = mean(Return, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(pct_increase = 100 * (scenario_period_mean - low_period_mean) / low_period_mean) %>%
  dplyr::arrange(Stock, scenario)

print(pct_increase_over_low_period)

# Combined status + low-period recovery summary, one row per at-risk stock
recovery_summary <- status_comparison_table %>%
  dplyr::left_join(
    pct_increase_over_low_period %>%
      dplyr::select(Stock, scenario, pct_increase) %>%
      tidyr::pivot_wider(names_from = scenario, values_from = pct_increase,
                         names_prefix = "pct_increase_"),
    by = "Stock"
  )

print(recovery_summary)

# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

## Bar chart: % increase in mean abundance over each stock's
## low-abundance period, by scenario
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

## Time series: % increase over time, faceted by stock, low-abundance
## period shaded
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

## Status-comparison heatmap: does the stock stay Endangered/Threatened
## under each scenario, or does it drop out of at-risk status?
status_plot_df <- status_summary_latest %>%
  dplyr::filter(Stock %in% status_comparison_table$Stock) %>%
  dplyr::mutate(
    scenario = factor(scenario, levels = c("Observed", "Pinniped scenario", "SST scenario", "Pink scenario")),
    status = factor(status, levels = c("Endangered", "Threatened",
                                       "Not at Risk / Special Concern"))
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

# ------------------------------------------------------------
# DIAGNOSTIC: why does a stock's classification worsen under a
# scenario despite a favorable (negative) covariate coefficient?
#
# Two likely causes:
#   (a) The retrospective model is recursive (spawners -> recruits ->
#       spawners...), so a productivity boost raises escapement, which
#       feeds into the density-dependent (spawners) term next
#       generation. A stock with a comparatively large |rb| relative
#       to its intercept (ra) can "overshoot" and then crash on the
#       next cycle (Ricker overcompensation) rather than improving
#       monotonically.
#   (b) The status classification above uses only the single most
#       recent classifiable year (slice_max(Year)); if that year lands
#       on the down-swing of an overcompensation cycle, the stock
#       looks worse than its average scenario trajectory would suggest.
#
# This plots the full retroS / retroR / retro_lnRS trajectory for one
# stock across scenarios, and reports how far outside the stock's
# historical covariate range the pinniped freeze-year value sits
# (extrapolation risk).
# ------------------------------------------------------------

diagnose_stock <- function(stock_name) {
  
  traj <- dplyr::bind_rows(
    model_actual_retro   %>% dplyr::filter(Stock == stock_name) %>% dplyr::mutate(scenario = "Observed"),
    model_pinniped_retro %>% dplyr::filter(Stock == stock_name) %>% dplyr::mutate(scenario = "Pinniped scenario"),
    model_sst_retro       %>% dplyr::filter(Stock == stock_name) %>% dplyr::mutate(scenario = "SST scenario"),
    model_pink_retro       %>% dplyr::filter(Stock == stock_name) %>% dplyr::mutate(scenario = "Pink scenario")
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
  
  # How far outside its historical range does the pinniped freeze-year
  # covariate sit, for this stock specifically?
  hist_range <- covariates_main %>%
    dplyr::filter(Stock == stock_name) %>%
    dplyr::summarise(
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

# ------------------------------------------------------------
# DIAGNOSTIC: inspect a stock's full status trajectory year-by-year,
# to see exactly which year/terminal-year our decline_pct and
# gen_mean calc are evaluating -- and whether that year post-dates
# the data COSEWIC (2017) actually had available (their report was
# necessarily built on data through roughly 2014-2015, given
# publication lag).
# ------------------------------------------------------------

inspect_status_trajectory <- function(stock_name, scenario_label = "Observed") {
  status_all %>%
    dplyr::filter(Stock == stock_name, scenario == scenario_label) %>%
    dplyr::select(Year, gen_mean, decline_pct, status) %>%
    print(n = Inf)
}

inspect_status_trajectory("Late Stuart")


