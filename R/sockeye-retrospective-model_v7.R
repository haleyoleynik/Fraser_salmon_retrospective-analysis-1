# Sockeye retrospective model — scenario comparisons
# Adapted from Carl Walters Excel model
# Haley Oleynik
#
# Three covariate-driven scenarios: pinniped, SST, pink salmon. Each runs
# through the same pipeline so catch-lost and abundance-recovery estimates
# are directly comparable across drivers.
#
# v5 changes (vs v4):
#   - FIX: covariate selection for AIC-averaged stocks. v4 coalesced
#     unselected coefficients to 0 *before* testing !is.na(), so every stock
#     "selected" all 7 covariates. NA * 0 = NA in the matrix product, so an
#     NA in ANY covariate (e.g. pink pre-1950, seal/NPGO after 2016-17)
#     wiped out the covariate term for every stock in that year.
#   - FIX: low_periods is now built before the MSY section that uses it
#     (v4 referenced it ~300 lines before defining it).
#   - FIX: all-stock totals only use years where every stock has a value
#     (v4 summed with na.rm = TRUE, so years with a stock missing were
#     silently under-counted, e.g. post-2016 when covariates run out).
#   - FIX: run_retro_model() starts each stock at its first year with usable
#     covariates/residuals, so a leading NA can't break a cycle line for
#     the whole series.
#   - All model runs stored long-format (one tibble with a `scenario`
#     column) instead of 8 separately named objects.
#   - ggsave() always gets an explicit plot = .
#   - NEW: harvest-rate sweep figure (observed / historic harvest rate /
#     fixed 0-80%), all stocks combined, paired bars by scenario.

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
# TRUE  = baseline-mean SST replaces actual SST in EVERY year (v4 behaviour,
#         so the 1950-1975 years also differ from actual)
# FALSE = actual SST through the end of the baseline, mean afterwards
#         (parallels the pinniped freeze)
SST_SCENARIO_ALL_YEARS <- TRUE

retroU_default   <- 0.3   # retrospective harvest rate cap (*_retro model runs)
useretro_default <- TRUE
yrretro_default  <- 1990  # year the retrospective harvest rate cap kicks in

# Harvest-rate sweep (fixed-rate bars)
HR_SWEEP_RATES      <- seq(0, 0.80, by = 0.05)
HR_SWEEP_START_YEAR <- yrretro_default   # fixed rate applies from this year on

# Window for mean-yearly-catch summaries. Covariates end in 2016-2017, so
# 2016 is the last year every stock is reliably complete; this also
# matches the chum figure (2000-2016) and the COSEWIC cap below.
CATCH_WINDOW        <- 2000:2016
MAX_ASSESSMENT_YEAR <- 2016

STOCKS <- c("Birkenhead", "Bowron", "Chilko", "Cultus", "Early Stuart", "Gates",
            "Late Shuswap", "Late Stuart", "Pitt", "Portage", "Quesnel", "Raft",
            "Scotch", "Seymour", "Stellako", "Weaver")

# Walters-file stock names that differ from the names used everywhere else
# (covariates, dredge output). Renamed on load.
STOCK_NAME_MAP <- c("Upper Pitt River" = "Pitt", "Scotch Creek" = "Scotch")

# Stocks COSEWIC (2017) designated Endangered (flagged with * in the
# weak-stock figures)
COSEWIC_ENDANGERED_STOCKS <- c("Bowron", "Weaver", "Quesnel", "Early Stuart",
                               "Late Stuart", "Portage", "Cultus")

# Weak-stock analysis: harvest rates at which per-stock spawner
# trajectories are drawn
WEAK_STOCK_TRAJ_RATES <- c(0.4, 0.6)

# Generation time (years), per COSEWIC (2017) Technical Summaries.
# All stocks are 4 years except Pitt (5).
GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5
GT_TABLE <- tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT")

# Right-aligned trailing mean; requires a full window of non-NA values.
trailing_mean <- function(x, width) {
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    if (i >= width) {
      window <- x[(i - width + 1):i]
      if (all(!is.na(window))) out[i] <- mean(window)
    }
  }
  out
}

# Scenario name -> covariate(s) swapped to their "<var>_scenario" column
SCENARIOS <- list(
  "Pinniped scenario" = c("SeaLions", "seal"),
  "SST scenario"      = c("adult.sst", "smolt.sst"),
  "Pink scenario"     = "pink"
)

SCENARIO_COLORS <- c(
  "Observed"          = "black",
  "Pinniped scenario" = "#4682B4",
  "SST scenario"      = "#2E8B57",
  "Pink scenario"     = "#FF4500"
)
SCENARIO_COLORS_NO_OBS <- SCENARIO_COLORS[names(SCENARIOS)]

dir.create("figures", showWarnings = FALSE)

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
  mutate(Stock = trimws(Stock),
         Stock = dplyr::recode(Stock, !!!STOCK_NAME_MAP)) %>%
  filter(!is.na(Year), Stock %in% STOCKS) %>%
  # Fill in any missing years (as NA rows) BEFORE joining covariates, so a
  # missing row still gets its covariates and can be gap-filled below.
  group_by(Stock) %>%
  tidyr::complete(Year = seq(min(Year), max(Year))) %>%
  ungroup() %>%
  arrange(Stock, Year)

# Stocks in STOCKS that have no rows in the Walters file (usually a
# spelling mismatch). All-stock totals are taken over MODELLED_STOCKS.
missing_stocks <- setdiff(STOCKS, unique(obs_raw$Stock))
if (length(missing_stocks) > 0) {
  warning("No rows in the Walters file for: ", paste(missing_stocks, collapse = ", "),
          " -- check spelling. Stock names in the file are: ",
          paste(sort(unique(read_csv("R/Sockeye Retrospective Shiny App/Walters_model_all-stocks.csv",
                                     show_col_types = FALSE)$Stock)), collapse = ", "))
}
MODELLED_STOCKS <- intersect(STOCKS, unique(obs_raw$Stock))

covariates_main <- read_csv("Data/sockeye_standardized_covariates.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, any_of(COVARIATE_COLS))

covariates_pink_wild <- read_csv("Data/sockeye_standardized_covariates_pink-wild.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, pink_wild)

dredge_models <- read_csv("Data/sockeye_top_models_dredge_wo-aquaculture.csv") %>%
  mutate(across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS), deltaAIC), as.numeric))

# AIC-weighted average coefficients (unselected terms count as 0)
top_models <- dredge_models %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )

# Which covariates appear in at least one of a stock's top models. Must be
# taken from the RAW dredge rows -- after averaging, nothing is NA.
top_model_selected <- dredge_models %>%
  group_by(Stock) %>%
  summarise(across(all_of(COVARIATE_COLS), ~ any(!is.na(.x))), .groups = "drop")

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
    else if_else(Year > FREEZE_YEAR, freeze_val, x)
  })) %>%
  ungroup() %>%
  select(Stock, Year, SeaLions_scenario = SeaLions, seal_scenario = seal)

## Pink: substituted with the wild-only pink series (same standardized scale)
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

## SST: adult.sst / smolt.sst set to each stock's baseline-period mean
sst_means <- covariates_main %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  group_by(Stock) %>%
  summarise(
    adult.sst_mean = mean(adult.sst, na.rm = TRUE),
    smolt.sst_mean = mean(smolt.sst, na.rm = TRUE),
    .groups = "drop"
  )

no_sst_mean <- sst_means %>% filter(is.na(adult.sst_mean) | is.na(smolt.sst_mean))
if (nrow(no_sst_mean) > 0) {
  warning("No usable ", min(SST_BASELINE_YEARS), "-", max(SST_BASELINE_YEARS),
          " SST mean for stock(s): ", paste(no_sst_mean$Stock, collapse = ", "))
}

sst_scenario_cov <- covariates_main %>%
  select(Stock, Year, adult.sst, smolt.sst) %>%
  left_join(sst_means, by = "Stock") %>%
  mutate(
    swap = SST_SCENARIO_ALL_YEARS | Year > max(SST_BASELINE_YEARS),
    adult.sst_scenario = if_else(swap, adult.sst_mean, adult.sst),
    smolt.sst_scenario = if_else(swap, smolt.sst_mean, smolt.sst)
  ) %>%
  select(Stock, Year, adult.sst_scenario, smolt.sst_scenario)

obs <- obs_raw %>%
  left_join(covariates_main,       by = c("Stock", "Year")) %>%
  left_join(pinniped_scenario_cov, by = c("Stock", "Year")) %>%
  left_join(pink_scenario_cov,     by = c("Stock", "Year")) %>%
  left_join(sst_scenario_cov,      by = c("Stock", "Year"))

# ============================================================
# LATE SHUSWAP (cycle-line model)
# Modeled by 4-year cycle line (cycle = Year %% 4), not the single
# stock-wide intercept/slope every other stock uses. Read directly from the
# cycle-stratified dredge output so a re-dredge is picked up automatically.
# ============================================================

late_shuswap_dredge <- read_csv("Data/Shuswap_sockeye_dredge-results.csv")

# Stops loudly if a re-dredge produces >1 top model -- decide how to
# combine them before proceeding.
late_shuswap_top <- late_shuswap_dredge %>% filter(deltaAIC == min(deltaAIC))
stopifnot(
  "Expected exactly one Late Shuswap top model (deltaAIC=0) -- got more than one; decide how to combine them before proceeding." =
    nrow(late_shuswap_top) == 1
)

# Intercept: cycle 0 is the reference level; cycles 1-3 add their
# factor(cycle)N offset (0 if not selected).
# Spawner slope: no bare "spawners" main effect, only the interaction, so
# each cycle has its own factor(cycle)N:spawners column.
late_shuswap_cycle_terms <- tibble(
  cycle = 0:3,
  ra = late_shuswap_top[["(Intercept)"]] + c(
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

late_shuswap_cov_present <- intersect(COVARIATE_COLS, colnames(late_shuswap_top))
late_shuswap_sel_covs <- late_shuswap_cov_present[
  !is.na(as.numeric(late_shuswap_top[late_shuswap_cov_present]))
]
late_shuswap_cov_coefs <- as.numeric(late_shuswap_top[late_shuswap_sel_covs])

# ============================================================
# MODEL FUNCTIONS
# ============================================================

# Stock-recruit terms for a stock: AIC-averaged top-model coefficients, or
# a plain Ricker fit if the stock has no dredge top model.
get_top_model_terms <- function(stock_name, fit_df_for_fallback = NULL) {
  
  if (stock_name == "Late Shuswap") {
    # ra/rb vary by cycle -- run_retro_model() pulls them from cycle_terms
    return(list(ra = NA_real_, rb = NA_real_,
                sel_covs = late_shuswap_sel_covs, cov_coefs = late_shuswap_cov_coefs,
                cycle_terms = late_shuswap_cycle_terms))
  }
  
  top_model_row <- top_models %>% filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_row  <- top_model_selected %>% filter(Stock == stock_name)
    sel_covs <- COVARIATE_COLS[unlist(sel_row[COVARIATE_COLS])]
    list(
      ra        = top_model_row[["(Intercept)"]],
      rb        = -top_model_row[["spawners"]],
      sel_covs  = sel_covs,
      cov_coefs = as.numeric(top_model_row[sel_covs])
    )
  } else {
    if (is.null(fit_df_for_fallback) || !"lnR_S" %in% names(fit_df_for_fallback)) {
      stop("No dredge top model for ", stock_name, " and no fit data for a fallback Ricker fit.")
    }
    warning("No dredge top model for ", stock_name, " -- using plain Ricker fit.")
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df_for_fallback)
    list(ra = unname(coef(fit)[1]), rb = -unname(coef(fit)[2]),
         sel_covs = character(0), cov_coefs = numeric(0))
  }
}

# Linear-predictor contribution from selected covariates. scenario_vars
# names which of sel_covs are read from their "<var>_scenario" column.
compute_cov_term <- function(dat, sel_covs, cov_coefs, scenario_vars = character(0)) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  cols <- ifelse(sel_covs %in% scenario_vars, paste0(sel_covs, "_scenario"), sel_covs)
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}

# Adds per-row ra/rb columns (cycle-specific for Late Shuswap, constant
# otherwise) so every downstream calculation can work row-wise.
add_ricker_terms <- function(dat, stock_name, terms) {
  if (stock_name == "Late Shuswap") {
    dat %>% mutate(cycle = Year %% 4) %>% left_join(terms$cycle_terms, by = "cycle")
  } else {
    dat %>% mutate(ra = terms$ra, rb = terms$rb)
  }
}

# Retrospective stock-recruit + harvest projection for one stock.
# scenario_vars swaps covariates to their scenario column in the forward
# projection only -- process-error residuals (wt) are always estimated
# against actual covariates, so a scenario comparison isolates the
# covariate effect.
run_retro_model <- function(dat, stock_name, retroU, useretro, yrretro,
                            scenario_vars = character(0)) {
  
  # (Year sequence is already complete -- see obs_raw -- so lead() by rows
  # lines recruits up with the right brood year.)
  obs2 <- dat %>%
    arrange(Year) %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch    = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs   = pmin(0.95, Catch / RunJacks),
      Ut_obs   = if_else(
        stock_name == "Late Shuswap" & Year == 2012,  # manual correction to match historical U
        pmin(0.7, Catch / RunJacks),
        Ut_obs
      ),
      ENS         = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort     = 1 - ENS,
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S       = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>% filter(Year %in% FIT_YEARS, is.finite(lnR_S), is.finite(AdultEscapement))
  terms  <- get_top_model_terms(stock_name, fit_df)
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs),
      cov_term_proj = compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars)
    ) %>%
    add_ricker_terms(stock_name, terms) %>%
    mutate(wt = lnR_S - (ra - rb * AdultEscapement + cov_term))
  
  # Gap-filling. One unusable brood year (zero/missing escapement or return,
  # or a missing row) used to give NA for that year's recruits, and because
  # the model is recursive the NA then propagated down that whole cycle line
  # (e.g. Portage lost every 4th year from ~2000 on, which dropped those
  # years from ALL all-stock totals). Instead:
  #   wt  (process error) -> 0, i.e. the deterministic model prediction
  #   Ut_obs / ENS        -> the stock's median
  # Only gaps after the first usable brood year and while covariates are
  # still available are filled; filled rows are flagged in `gap_filled`.
  first_ok <- which(is.finite(obs3$wt) & is.finite(obs3$cov_term_proj))[1]
  if (is.na(first_ok)) stop("No usable brood years for ", stock_name)
  in_range <- seq_len(nrow(obs3)) >= first_ok & is.finite(obs3$cov_term_proj) &
    seq_len(nrow(obs3)) <= nrow(obs3) - LAG_YEARS
  obs3 <- obs3 %>%
    mutate(
      gap_filled = in_range & (!is.finite(wt) | !is.finite(Ut_obs) | !is.finite(ENS)),
      wt     = if_else(in_range & !is.finite(wt), 0, wt),
      Ut_obs = if_else(is.finite(Ut_obs), Ut_obs, median(Ut_obs[is.finite(Ut_obs)])),
      ENS    = if_else(is.finite(ENS), ENS, median(ENS[is.finite(ENS)])),
      Catch  = if_else(is.finite(RunJacks), Catch, NA_real_)
    )
  
  # Start at the first brood year with a usable residual + projection term.
  obs3 <- obs3[first_ok:nrow(obs3), ]
  
  n <- nrow(obs3)
  retroR     <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroS     <- rep(NA_real_, n)
  retroC     <- rep(NA_real_, n)
  retroU_vec <- if (useretro) if_else(obs3$Year >= yrretro, retroU, obs3$Ut_obs) else obs3$Ut_obs
  
  # Each cycle line is seeded from its first observed run. Normally that is
  # just the first LAG_YEARS rows; if a cycle has no data there (missing
  # rows), it is seeded from its first observed run later on instead of
  # staying NA for the whole series.
  cycle_started <- rep(FALSE, LAG_YEARS)
  for (i in seq_len(n)) {
    cyc <- (i - 1) %% LAG_YEARS + 1
    if (i > LAG_YEARS && cycle_started[cyc]) {
      j <- i - LAG_YEARS
      retro_lnRS[j] <- obs3$ra[j] - obs3$rb[j] * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
      retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    } else if (is.finite(obs3$RunJacks[i])) {
      retroR[i] <- obs3$RunJacks[i]
      cycle_started[cyc] <- TRUE
    }
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

# Runs the actual-covariate model (scenario = "Actual") plus every
# scenario in SCENARIOS, stacked long.
run_all_scenarios <- function(retroU, useretro, yrretro, include_actual = TRUE) {
  scen <- if (include_actual) c(list(Actual = character(0)), SCENARIOS) else SCENARIOS
  imap_dfr(scen, ~ run_scenario_model(obs, retroU, useretro, yrretro, scenario_vars = .x) %>%
             mutate(scenario = .y))
}

# All-stock total of `value_col` by the grouping variables. A year only
# gets a total if every modelled stock has a finite value -- otherwise NA,
# so a missing stock can't masquerade as a drop in catch/abundance.
sum_complete_stocks <- function(df, value_col, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n_stocks = n_distinct(Stock[is.finite(.data[[value_col]])]),
      total    = sum(.data[[value_col]], na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(total = if_else(n_stocks == length(MODELLED_STOCKS), total, NA_real_))
}

# ============================================================
# RUN MODELS
# runs_retro : capped retrospective harvest rate (productivity/return/spawner comparisons)
# runs_hist  : actual historical harvest rate ("catch lost to driver" comparisons)
# ============================================================

runs_retro <- run_all_scenarios(retroU_default, useretro_default, yrretro_default)
runs_hist  <- run_all_scenarios(retroU = 0, useretro = FALSE, yrretro = yrretro_default)

# Brood years where an unusable residual/harvest rate/ENS was gap-filled
gap_filled_years <- runs_hist %>%
  filter(scenario == "Actual", gap_filled) %>%
  group_by(Stock) %>%
  summarise(filled_brood_years = paste(Year, collapse = ","), .groups = "drop")
if (nrow(gap_filled_years) > 0) {
  message("Gap-filled brood years (process error set to 0 / median U, ENS):")
  print(gap_filled_years, n = Inf)
}

# Coverage check: first/last projected year per stock (actual covariates)
projection_coverage <- runs_hist %>%
  filter(scenario == "Actual", is.finite(retroC)) %>%
  group_by(Stock) %>%
  summarise(first_year = min(Year), last_year = max(Year),
            n_gaps = (last_year - first_year + 1) - n(), .groups = "drop")
print(projection_coverage, n = Inf)

if (any(projection_coverage$n_gaps > 0)) {
  warning("Mid-series gaps in projected catch for: ",
          paste(projection_coverage$Stock[projection_coverage$n_gaps > 0], collapse = ", "),
          " -- an NA residual/covariate breaks that cycle line from that year on.")
}

# ============================================================
# CATCH LOST TO EACH DRIVER (historical harvest rate)
# catch_lost = scenario catch - actual catch (positive = catch that would
# have been available had the driver stayed at its scenario level)
# ============================================================

catch_lost_by_stock <- runs_hist %>%
  filter(scenario != "Actual") %>%
  select(Stock, Year, scenario, catch_scenario = retroC) %>%
  left_join(runs_hist %>% filter(scenario == "Actual") %>%
              select(Stock, Year, catch_actual = retroC),
            by = c("Stock", "Year")) %>%
  mutate(catch_lost = catch_scenario - catch_actual) %>%
  arrange(scenario, Stock, Year) %>%
  group_by(scenario, Stock) %>%
  mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
  ungroup()

catch_lost_totals <- catch_lost_by_stock %>%
  sum_complete_stocks("catch_lost", scenario, Year) %>%
  rename(catch_lost = total) %>%
  filter(!is.na(catch_lost)) %>%
  arrange(scenario, Year) %>%
  group_by(scenario) %>%
  mutate(cum_catch_lost = cumsum(catch_lost)) %>%
  ungroup()

# ============================================================
# HARVEST-RATE SWEEP: mean yearly catch & catch lost, all stocks combined
#   Observed               : observed catch (Below + Above Mission)
#   Historic harvest rate  : each scenario at the actual historical U
#   0%-80%                 : each scenario at a fixed U from HR_SWEEP_START_YEAR,
#                            plus "Historical conditions" (grey; actual
#                            covariates at that fixed U -- top panel only)
# Catch lost = scenario catch - catch under historical conditions at the
# SAME harvest rate (grey bar for fixed rates; actual-covariate model at the
# actual U for the historic group). This isolates the driver effect: at 0%
# harvest it is exactly 0. (The chum figure subtracted OBSERVED catch,
# which mixes the driver effect with the harvest-rate change and gives
# large negatives at low rates.)
# No grey bar in the historic-rate group: actual covariates at the actual
# U reproduce the observed catch, so it would duplicate the black bar.
# ============================================================

sweep_runs <- map_dfr(HR_SWEEP_RATES, function(u) {
  run_all_scenarios(retroU = u, useretro = TRUE, yrretro = HR_SWEEP_START_YEAR,
                    include_actual = TRUE) %>%
    mutate(harvest_rate = u,
           scenario = if_else(scenario == "Actual", "Historical conditions", scenario))
})

HR_COLORS <- c(
  "Observed"              = "black",
  "Historical conditions" = "grey55",
  SCENARIO_COLORS_NO_OBS
)

HR_LEVELS <- c("Observed", "Historic\nharvest rate", paste0(round(HR_SWEEP_RATES * 100), "%"))

observed_catch_total <- runs_hist %>%
  filter(scenario == "Actual") %>%
  sum_complete_stocks("Catch", Year) %>%
  select(Year, observed_catch = total)

historic_baseline_total <- runs_hist %>%
  filter(scenario == "Actual") %>%
  sum_complete_stocks("retroC", Year) %>%
  transmute(Year, hr_group = "Historic\nharvest rate", baseline_catch = total)

fixed_baseline_total <- sweep_runs %>%
  filter(scenario == "Historical conditions") %>%
  sum_complete_stocks("retroC", harvest_rate, Year) %>%
  transmute(Year, hr_group = paste0(round(harvest_rate * 100), "%"), baseline_catch = total)

catch_baseline <- bind_rows(historic_baseline_total, fixed_baseline_total)

annual_catch_by_hr <- bind_rows(
  observed_catch_total %>%
    transmute(Year, scenario = "Observed", hr_group = "Observed", catch = observed_catch),
  runs_hist %>%
    filter(scenario != "Actual") %>%
    sum_complete_stocks("retroC", scenario, Year) %>%
    transmute(Year, scenario, hr_group = "Historic\nharvest rate", catch = total),
  sweep_runs %>%
    sum_complete_stocks("retroC", scenario, harvest_rate, Year) %>%
    transmute(Year, scenario, hr_group = paste0(round(harvest_rate * 100), "%"), catch = total)
) %>%
  filter(Year %in% CATCH_WINDOW) %>%
  left_join(catch_baseline, by = c("Year", "hr_group")) %>%
  mutate(catch_lost = if_else(scenario %in% c("Observed", "Historical conditions"),
                              NA_real_, catch - baseline_catch),
         hr_group   = factor(hr_group, levels = HR_LEVELS),
         scenario   = factor(scenario, levels = names(HR_COLORS)))

# Which stock is missing in which window year (empty = all good). Any
# stock/year listed here blanks that year's all-stock total.
window_coverage <- bind_rows(
  runs_hist %>% mutate(run = "historic rate"),
  sweep_runs %>% mutate(run = paste0(round(harvest_rate * 100), "% fixed")) %>% select(-harvest_rate)
) %>%
  select(run, scenario, Stock, Year, retroC) %>%
  right_join(tidyr::expand_grid(distinct(bind_rows(
    runs_hist %>% mutate(run = "historic rate"),
    sweep_runs %>% mutate(run = paste0(round(harvest_rate * 100), "% fixed"))), run, scenario),
    Stock = MODELLED_STOCKS, Year = CATCH_WINDOW),
    by = c("run", "scenario", "Stock", "Year")) %>%
  filter(!is.finite(retroC)) %>%
  group_by(Stock, scenario, run) %>%
  summarise(missing_years = paste(sort(Year), collapse = ","), .groups = "drop")

if (nrow(window_coverage) > 0) {
  warning("Some stock/years in CATCH_WINDOW have no projected catch, so those years drop ",
          "out of the all-stock means -- see `window_coverage`.")
  print(distinct(window_coverage, Stock, scenario, missing_years), n = 50)
}

if (all(is.na(annual_catch_by_hr$catch))) {
  stop("No year in CATCH_WINDOW has a complete all-stock total. ",
       "Check `window_coverage`, `projection_coverage` and the stock-name warning above.")
}

hr_sweep_summary <- annual_catch_by_hr %>%
  group_by(hr_group, scenario) %>%
  summarise(
    n_years         = sum(!is.na(catch)),
    mean_catch      = mean(catch, na.rm = TRUE),
    se_catch        = sd(catch, na.rm = TRUE) / sqrt(n_years),
    mean_catch_lost = mean(catch_lost, na.rm = TRUE),
    se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n_years),
    .groups = "drop"
  )

print(hr_sweep_summary, n = Inf)

# (The sweep figure itself is drawn after the MSY section, which it needs.)

# ============================================================
# RETURNS (long format) + COSEWIC helpers + LOW-ABUNDANCE PERIODS
# Built here because the MSY section below needs low_periods.
# ============================================================

returns_by_stock <- bind_rows(
  runs_retro %>% filter(scenario == "Actual") %>%
    transmute(Stock, Year, scenario = "Observed", Return = RunJacks),
  runs_retro %>% filter(scenario != "Actual") %>%
    transmute(Stock, Year, scenario, Return = retroR)
)

# (GENERATION_TIME, GT_TABLE and trailing_mean() are defined in CONFIG.)

# Each stock's low-abundance PERIOD: the GT-year trailing window (observed
# series) with the lowest mean return.
low_periods <- returns_by_stock %>%
  filter(scenario == "Observed", is.finite(Return)) %>%
  left_join(GT_TABLE, by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(obs_gen_mean = trailing_mean(Return, width = first(GT))) %>%
  filter(!is.na(obs_gen_mean)) %>%
  slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, GT, low_period_end = Year, low_period_start = Year - GT + 1,
            low_period_mean = obs_gen_mean)

# ============================================================
# MSY UNDER HIGH VS. LOW PRODUCTIVITY
# High: 1950-1970 (brood years).
# Low: each stock's own LOW-PRODUCTIVITY window -- the generation-length run
# of brood years with the lowest mean realized Ricker productivity,
# a_t = ln(R/S) + b*S (density-independent, so a big escapement doesn't
# masquerade as low productivity). Only brood years with complete
# covariates are eligible, so MSY can always be computed for the window.
# (v5 previously used low_periods, the lowest-ABUNDANCE window by return
# year; for most stocks that falls in 2017-2020, after the covariates end,
# so every low MSY came out NA. low_periods is still used for the
# abundance-recovery analysis below, where abundance is the right measure.)
# MSY = Rmsy - Smsy, Smsy = (a/b)*(0.5-0.07*a), Rmsy = Smsy*exp(a-b*Smsy),
# computed per year with that year's real covariates (and cycle-matched
# a/b for Late Shuswap), then averaged over the period.
# ============================================================

SMSY_PERIOD_START_YEAR <- 1950
SMSY_PERIOD_END_YEAR   <- 1970

compute_msy_by_stock <- function(stock_name, period_start, period_end) {
  
  terms <- get_top_model_terms(stock_name)
  
  period_data <- obs %>%
    filter(Stock == stock_name, Year >= period_start, Year <= period_end) %>%
    filter(if_all(all_of(terms$sel_covs), ~ !is.na(.))) %>%
    arrange(Year)
  if (nrow(period_data) == 0) return(NA_real_)
  
  period_data <- add_ricker_terms(period_data, stock_name, terms)
  a_vec <- period_data$ra + compute_cov_term(period_data, terms$sel_covs, terms$cov_coefs)
  b_vec <- period_data$rb
  
  Smsy_vec <- (a_vec / b_vec) * (0.5 - 0.07 * a_vec)
  Rmsy_vec <- Smsy_vec * exp(a_vec - b_vec * Smsy_vec)
  mean(Rmsy_vec - Smsy_vec, na.rm = TRUE)
}

msy_high_by_stock <- tibble(
  Stock = STOCKS,
  MSY   = map_dbl(STOCKS, compute_msy_by_stock,
                  period_start = SMSY_PERIOD_START_YEAR, period_end = SMSY_PERIOD_END_YEAR)
)

low_prod_periods <- runs_hist %>%
  filter(scenario == "Actual") %>%
  left_join(GT_TABLE, by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(
    prod          = if_else(is.finite(cov_term), lnR_S + rb * AdultEscapement, NA_real_),
    prod_gen_mean = trailing_mean(prod, width = first(GT))
  ) %>%
  filter(!is.na(prod_gen_mean)) %>%
  slice_min(prod_gen_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, low_prod_start = Year - GT + 1, low_prod_end = Year,
            low_prod_mean_a = prod_gen_mean)

print(low_prod_periods, n = Inf)

msy_low_by_stock <- tibble(
  Stock = STOCKS,
  MSY   = map_dbl(STOCKS, function(s) {
    low_p <- low_prod_periods %>% filter(Stock == s)
    if (nrow(low_p) == 0) return(NA_real_)
    compute_msy_by_stock(s, low_p$low_prod_start, low_p$low_prod_end)
  })
)

print(msy_high_by_stock, n = Inf)
print(msy_low_by_stock, n = Inf)

msy_high_avg <- mean(msy_high_by_stock$MSY, na.rm = TRUE)
msy_low_avg  <- mean(msy_low_by_stock$MSY, na.rm = TRUE)

cat("Average MSY across stocks -- High productivity (1950-1970):", msy_high_avg, "\n")
cat("Average MSY across stocks -- Low productivity (own low period):", msy_low_avg, "\n")

MSY_LINES <- tibble(
  label = c("MSY (high productivity)", "MSY (low productivity)"),
  value = c(msy_high_avg, msy_low_avg)
)
MSY_LINETYPES <- c("MSY (high productivity)" = "dashed", "MSY (low productivity)" = "dotted")

# Recruit-weighted average MSY: each stock weighted by its mean recruits
# (AdultReturn by brood year, actual data) over the full record, so the
# high and low lines share the same weights and differ only in MSY.
recruit_weights <- runs_hist %>%
  filter(scenario == "Actual") %>%
  group_by(Stock) %>%
  summarise(mean_recruits = mean(AdultReturn, na.rm = TRUE), .groups = "drop")

weighted_msy <- function(msy_by_stock) {
  d <- msy_by_stock %>%
    inner_join(recruit_weights, by = "Stock") %>%
    filter(is.finite(MSY), is.finite(mean_recruits))
  weighted.mean(d$MSY, d$mean_recruits)
}

msy_high_wavg <- weighted_msy(msy_high_by_stock)
msy_low_wavg  <- weighted_msy(msy_low_by_stock)

cat("Recruit-weighted MSY -- High productivity (1950-1970):", msy_high_wavg, "\n")
cat("Recruit-weighted MSY -- Low productivity (own low period):", msy_low_wavg, "\n")

msy_summary <- msy_high_by_stock %>%
  rename(MSY_high = MSY) %>%
  left_join(rename(msy_low_by_stock, MSY_low = MSY), by = "Stock") %>%
  left_join(low_prod_periods, by = "Stock") %>%
  left_join(recruit_weights, by = "Stock") %>%
  mutate(weight = mean_recruits / sum(mean_recruits, na.rm = TRUE))
print(msy_summary, n = Inf, width = Inf)

if (any(!STOCKS %in% MODELLED_STOCKS)) {
  warning("MSY lines exclude stock(s) with no rows in the Walters file: ",
          paste(setdiff(STOCKS, MODELLED_STOCKS), collapse = ", "),
          " -- they are also missing from every all-stock total.")
}

# ============================================================
# HARVEST-RATE SWEEP FIGURE (two versions: MSY lines as a simple mean
# across stocks, and as a recruit-weighted mean)
# ============================================================

hr_dodge <- position_dodge(width = 0.85)
hr_fill  <- scale_fill_manual(values = HR_COLORS, limits = names(HR_COLORS), name = NULL)

# Top panel: fixed-rate groups have 4 bars, the historic group 3. dodge2
# with preserve = "single" keeps every bar the same width and centred.
hr_dodge2_col <- position_dodge2(preserve = "single", padding = 0.1)
hr_dodge2_err <- position_dodge2(preserve = "single", padding = 0.6)
hr_top_bars   <- filter(hr_sweep_summary, scenario != "Observed")

# MSY reference lines for the top panel: one dashed line per productivity
# period, labelled at the right-hand end.
msy_line_layers <- function(msy_high, msy_low) {
  geom_hline(yintercept = c(msy_high, msy_low),
             linetype = "dashed", colour = "firebrick", linewidth = 0.7)
}

# lost_summary: table supplying the bottom panel (mean_catch_lost /
# se_catch_lost by hr_group and scenario). Defaults to catch lost relative
# to historical conditions at the same harvest rate.
make_hr_sweep_plot <- function(msy_high, msy_low, msy_label,
                               lost_summary = hr_sweep_summary) {
  
  # (separator lines go after the first bar layer so the x scale is discrete)
  p_hr_catch <- ggplot(mapping = aes(hr_group, fill = scenario)) +
    geom_col(data = filter(hr_sweep_summary, scenario == "Observed"),
             aes(y = mean_catch), width = 0.5) +
    geom_vline(xintercept = 2.5, linetype = "dashed", colour = "grey65", linewidth = 0.4) +
    geom_errorbar(data = filter(hr_sweep_summary, scenario == "Observed"),
                  aes(ymin = mean_catch - se_catch, ymax = mean_catch + se_catch),
                  width = 0.2, linewidth = 0.4) +
    geom_col(data = hr_top_bars,
             aes(y = mean_catch), position = hr_dodge2_col, width = 0.85) +
    geom_errorbar(data = hr_top_bars,
                  aes(ymin = mean_catch - se_catch, ymax = mean_catch + se_catch),
                  position = hr_dodge2_err, width = 0.85, linewidth = 0.4) +
    msy_line_layers(msy_high, msy_low) +
    hr_fill +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, NA)) +
    labs(x = NULL, y = "Mean yearly catch") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  
  p_hr_lost <- ggplot(filter(lost_summary, !scenario %in% c("Observed", "Historical conditions")),
                      aes(hr_group, mean_catch_lost, fill = scenario)) +
    geom_col(position = hr_dodge, width = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_vline(xintercept = 1.5, linetype = "dashed", colour = "grey65", linewidth = 0.4) +
    geom_errorbar(aes(ymin = mean_catch_lost - se_catch_lost,
                      ymax = mean_catch_lost + se_catch_lost),
                  position = hr_dodge, width = 0.25, linewidth = 0.4) +
    hr_fill +
    scale_y_continuous(labels = scales::comma) +
    guides(fill = "none") +   # legend comes from the top panel
    labs(x = "Harvest rate scenario", y = "Mean yearly catch lost") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
  
  (p_hr_catch / p_hr_lost) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

p_hr_sweep_msy_mean <- make_hr_sweep_plot(
  msy_high_avg, msy_low_avg, "MSY averaged across stocks (simple mean)")
p_hr_sweep_msy_wmean <- make_hr_sweep_plot(
  msy_high_wavg, msy_low_wavg, "MSY averaged across stocks, weighted by mean recruits")

p_hr_sweep_msy_mean
p_hr_sweep_msy_wmean

sweep_file_stub <- paste0("figures/sockeye_mean_catch_observed_vs_scenarios_harvest_sweep_",
                          min(CATCH_WINDOW), "-", max(CATCH_WINDOW))
ggsave(paste0(sweep_file_stub, "_MSY-mean.png"), plot = p_hr_sweep_msy_mean,
       width = 15, height = 10, dpi = 600, bg = "white")
ggsave(paste0(sweep_file_stub, "_MSY-recruit-weighted.png"), plot = p_hr_sweep_msy_wmean,
       width = 15, height = 10, dpi = 600, bg = "white")

# ============================================================
# ALTERNATIVE SWEEP FIGURE: bottom panel = scenario catch - OBSERVED catch
# (the earlier definition). Same top panel as p_hr_sweep_msy_mean. Note
# this mixes the driver effect with the change in harvest rate, so at 0%
# every scenario shows -(observed catch).
# ============================================================

hr_sweep_summary_vs_obs <- annual_catch_by_hr %>%
  filter(!scenario %in% c("Observed", "Historical conditions")) %>%
  left_join(observed_catch_total, by = "Year") %>%
  mutate(catch_lost_obs = catch - observed_catch) %>%
  group_by(hr_group, scenario) %>%
  summarise(
    n_years         = sum(!is.na(catch_lost_obs)),
    mean_catch_lost = mean(catch_lost_obs, na.rm = TRUE),
    se_catch_lost   = sd(catch_lost_obs, na.rm = TRUE) / sqrt(n_years),
    .groups = "drop"
  )

p_hr_sweep_msy_mean_vs_obs <- make_hr_sweep_plot(
  msy_high_avg, msy_low_avg, "MSY averaged across stocks (simple mean)",
  lost_summary = hr_sweep_summary_vs_obs)

p_hr_sweep_msy_mean_vs_obs

ggsave(paste0(sweep_file_stub, "_MSY-mean_lost-vs-observed.png"),
       plot = p_hr_sweep_msy_mean_vs_obs, width = 15, height = 10, dpi = 600, bg = "white")

# ============================================================
# WEAK STOCKS: per-stock status across the harvest-rate grid, by scenario
#
# Question: could the weaker stocks have held up better -- tolerated
# higher harvest rates -- under the alternative ecological scenarios?
#
# Benchmarks (Wild Salmon Policy style), computed per stock for TWO
# productivity periods -- the same periods as the MSY lines:
#   High productivity : brood years SMSY_PERIOD_START_YEAR-SMSY_PERIOD_END_YEAR
#   Low productivity  : each stock's lowest-productivity generation
#                       (low_prod_periods, from the MSY section)
# For each period:
#   a      = ra + the stock's mean ACTUAL covariate term over that period
#            (cycle-specific ra/rb for Late Shuswap)
#   Smsy   = (a/b)(0.5 - 0.07a)
#   Sgen   = spawners that produce Smsy recruits in one generation,
#            i.e. Sgen * exp(a - b*Sgen) = Smsy
# Within a benchmark period, every scenario and harvest rate is judged
# against the SAME benchmarks -- only the spawner trajectories differ.
#
# Status = geometric mean of S_t / benchmark_t over the last generation of
# CATCH_WINDOW (2013-2016; 2012-2016 for Pitt):
#   Red   : spawners below Sgen            ("crashed" / weak)
#   Amber : between Sgen and 0.8 * Smsy
#   Green : at or above 0.8 * Smsy
# Every figure is produced once per benchmark period (file suffix
# _bm-high / _bm-low).
# ============================================================

STATUS_LEVELS <- c("Red (< Sgen)", "Amber", "Green (>= 0.8 Smsy)")
STATUS_COLORS <- setNames(c("#B22222", "#E8A33D", "#4C9A5B"), STATUS_LEVELS)
WEAK_SCENARIO_LEVELS <- c("Historical conditions", names(SCENARIOS))
WEAK_SCENARIO_COLORS <- HR_COLORS[WEAK_SCENARIO_LEVELS]

BM_PERIODS <- c(
  high = paste0("High productivity (", SMSY_PERIOD_START_YEAR, "–", SMSY_PERIOD_END_YEAR, ")"),
  low  = "Low productivity (each stock's lowest-productivity generation)"
)

solve_sgen <- function(a, b, smsy) {
  if (!is.finite(a) || !is.finite(b) || a <= 0 || b <= 0 || !is.finite(smsy) || smsy <= 0) {
    return(NA_real_)
  }
  uniroot(function(S) S * exp(a - b * S) - smsy,
          lower = smsy * 1e-8, upper = smsy, tol = 1e-6)$root
}

# Brood years defining each benchmark period, per stock
bm_period_years <- bind_rows(
  tibble(Stock = MODELLED_STOCKS, bm_period = BM_PERIODS[["high"]],
         period_start = SMSY_PERIOD_START_YEAR, period_end = SMSY_PERIOD_END_YEAR),
  low_prod_periods %>%
    transmute(Stock, bm_period = BM_PERIODS[["low"]],
              period_start = low_prod_start, period_end = low_prod_end)
)

actual_rows <- runs_hist %>%
  filter(scenario == "Actual") %>%
  mutate(cycle = Year %% 4) %>%
  select(Stock, Year, cycle, ra, rb, cov_term)

stock_benchmarks <- bm_period_years %>%
  left_join(actual_rows, by = "Stock", relationship = "many-to-many") %>%
  group_by(Stock, bm_period) %>%
  mutate(mean_cov = mean(cov_term[Year >= period_start & Year <= period_end], na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Stock, bm_period, period_start, period_end, cycle, ra, rb, mean_cov) %>%
  mutate(
    a_bm     = ra + mean_cov,
    Smsy     = (a_bm / rb) * (0.5 - 0.07 * a_bm),
    Sgen     = pmap_dbl(list(a_bm, rb, Smsy), solve_sgen),
    upper_bm = 0.8 * Smsy
  ) %>%
  select(Stock, bm_period, period_start, period_end, cycle, a_bm, Smsy, Sgen, upper_bm)

# One row per stock and period (Late Shuswap averaged over its 4 cycles)
stock_benchmarks_summary <- stock_benchmarks %>%
  group_by(bm_period, Stock, period_start, period_end) %>%
  summarise(a = mean(a_bm), Smsy = mean(Smsy), Sgen = mean(Sgen), .groups = "drop") %>%
  arrange(bm_period, Stock)
print(stock_benchmarks_summary, n = Inf, width = Inf)

# Per-stock runs: historic harvest rate + every fixed rate, all scenarios
WEAK_HR_LEVELS <- c("Historic", paste0(round(HR_SWEEP_RATES * 100), "%"))

stock_runs_base <- bind_rows(
  runs_hist %>% mutate(hr_group = "Historic", harvest_rate = NA_real_),
  sweep_runs %>% mutate(hr_group = paste0(round(harvest_rate * 100), "%"))
) %>%
  mutate(scenario = if_else(scenario == "Actual", "Historical conditions", scenario),
         cycle    = Year %% 4) %>%
  select(Stock, Year, cycle, scenario, hr_group, harvest_rate, retroS, retroR, retroC)

# Duplicated once per benchmark period
stock_runs <- stock_runs_base %>%
  left_join(stock_benchmarks %>% select(Stock, cycle, bm_period, Smsy, Sgen, upper_bm),
            by = c("Stock", "cycle"), relationship = "many-to-many")

stock_status <- stock_runs %>%
  left_join(GT_TABLE, by = "Stock") %>%
  filter(Year > max(CATCH_WINDOW) - GT, Year <= max(CATCH_WINDOW)) %>%
  group_by(bm_period, Stock, scenario, hr_group, harvest_rate) %>%
  summarise(
    n_years      = sum(is.finite(retroS)),
    gen_spawners = exp(mean(log(pmax(retroS, 1)), na.rm = TRUE)),
    ratio_sgen   = exp(mean(log(pmax(retroS, 1) / Sgen), na.rm = TRUE)),
    ratio_upper  = exp(mean(log(pmax(retroS, 1) / upper_bm), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    status = case_when(
      n_years == 0 | !is.finite(ratio_sgen) ~ NA_character_,
      ratio_sgen < 1  ~ STATUS_LEVELS[1],
      ratio_upper < 1 ~ STATUS_LEVELS[2],
      TRUE            ~ STATUS_LEVELS[3]
    ),
    status   = factor(status, levels = STATUS_LEVELS),
    scenario = factor(scenario, levels = WEAK_SCENARIO_LEVELS),
    hr_group = factor(hr_group, levels = WEAK_HR_LEVELS)
  )

stock_red_freq <- stock_runs %>%
  filter(Year %in% CATCH_WINDOW, is.finite(retroS), is.finite(Sgen)) %>%
  group_by(bm_period, Stock, scenario, hr_group, harvest_rate) %>%
  summarise(prop_below_sgen = mean(retroS < Sgen), .groups = "drop") %>%
  mutate(scenario = factor(scenario, levels = WEAK_SCENARIO_LEVELS),
         hr_group = factor(hr_group, levels = WEAK_HR_LEVELS))

n_red_by_hr <- stock_status %>%
  filter(!is.na(harvest_rate)) %>%
  group_by(bm_period, scenario, harvest_rate) %>%
  summarise(n_red    = sum(status == STATUS_LEVELS[1], na.rm = TRUE),
            n_amber  = sum(status == STATUS_LEVELS[2], na.rm = TRUE),
            n_stocks = sum(!is.na(status)), .groups = "drop")

# Highest fixed harvest rate each stock tolerates: above Sgen at that rate
# AND every lower rate
hr_step <- HR_SWEEP_RATES[2] - HR_SWEEP_RATES[1]
weak_thresholds <- stock_status %>%
  filter(!is.na(harvest_rate), !is.na(status)) %>%
  group_by(bm_period, Stock, scenario) %>%
  summarise(first_red = suppressWarnings(min(harvest_rate[status == STATUS_LEVELS[1]])),
            .groups = "drop") %>%
  mutate(
    outcome = case_when(
      !is.finite(first_red) ~ "Above Sgen at all rates tested",
      first_red == 0        ~ "Below Sgen even at 0%",
      TRUE                  ~ "Highest rate above Sgen"
    ),
    max_tolerated_hr = case_when(
      !is.finite(first_red) ~ max(HR_SWEEP_RATES),
      TRUE                  ~ pmax(first_red - hr_step, 0)
    )
  )

# Stock order (same in every figure): lowest productivity at the top, using
# productivity at average actual conditions over the full record.
# COSEWIC (2017) Endangered stocks are starred.
stock_order <- runs_hist %>%
  filter(scenario == "Actual", Year %in% FIT_YEARS) %>%
  group_by(Stock) %>%
  summarise(a_avg = mean(ra + cov_term, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(a_avg)) %>%
  pull(Stock)
stock_label <- function(x) ifelse(x %in% COSEWIC_ENDANGERED_STOCKS, paste0(x, "*"), x)

# Which stocks each scenario can affect: a scenario only changes a stock if
# at least one of its covariates is in that stock's top model(s). Rows for
# unaffected stocks are greyed out in the scenario panels of the heatmaps
# (their results are identical to historical conditions by construction).
scenario_applies <- expand_grid(Stock = MODELLED_STOCKS, scenario = names(SCENARIOS)) %>%
  mutate(applies = map2_lgl(Stock, scenario, function(st, sc) {
    any(SCENARIOS[[sc]] %in% get_top_model_terms(st)$sel_covs)
  }))
print(scenario_applies %>% pivot_wider(names_from = scenario, values_from = applies), n = Inf)

NOT_LINKED_LABEL <- "Stock not linked to scenario covariate(s)"

# Grey overlay for (stock, scenario) panels the scenario doesn't touch
not_linked_layer <- function(stock_levels, hr_levels) {
  df <- scenario_applies %>%
    filter(!applies) %>%
    select(Stock, scenario) %>%
    expand_grid(hr_group = factor(hr_levels, levels = hr_levels)) %>%
    mutate(Stock    = factor(Stock, levels = stock_levels),
           scenario = factor(scenario, levels = WEAK_SCENARIO_LEVELS),
           shade    = NOT_LINKED_LABEL)
  list(
    geom_tile(data = df, aes(hr_group, Stock, alpha = shade), inherit.aes = FALSE,
              fill = "grey65", colour = "white", linewidth = 0.4),
    scale_alpha_manual(values = setNames(1, NOT_LINKED_LABEL), name = NULL)
  )
}

# ------------------------------------------------------------
# Figures + tables for one benchmark period
# ------------------------------------------------------------
make_weak_stock_outputs <- function(bm_key) {
  
  bm      <- BM_PERIODS[[bm_key]]
  sfx     <- paste0("_bm-", bm_key)
  status  <- filter(stock_status, bm_period == bm)
  cat("\n==== Weak-stock results, benchmarks from:", bm, "====\n")
  
  # W1: status heatmap
  p_heat <- status %>%
    mutate(Stock = factor(Stock, levels = stock_order)) %>%
    ggplot(aes(hr_group, Stock, fill = status)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    not_linked_layer(stock_order, WEAK_HR_LEVELS) +
    facet_wrap(~ scenario, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = STATUS_COLORS, na.value = "grey90", drop = FALSE, name = NULL) +
    scale_y_discrete(labels = stock_label) +
    labs(x = "Harvest rate", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(), legend.position = "bottom")
  print(p_heat)
  ggsave(paste0("figures/weak_stocks_status_heatmap_by_scenario", sfx, ".png"),
         plot = p_heat, width = 16, height = 7, dpi = 600, bg = "white")
  
  # W1b: share of years in CATCH_WINDOW below Sgen (less sensitive to cycle
  # timing than the last-generation snapshot -- read alongside W1)
  p_freq <- stock_red_freq %>%
    filter(bm_period == bm) %>%
    mutate(Stock = factor(Stock, levels = stock_order)) %>%
    ggplot(aes(hr_group, Stock, fill = prop_below_sgen)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    not_linked_layer(stock_order, WEAK_HR_LEVELS) +
    facet_wrap(~ scenario, nrow = 1, drop = FALSE) +
    scale_fill_gradient(low = "white", high = "#B22222", limits = c(0, 1),
                        labels = scales::percent,
                        name = paste0("Years below Sgen, ", min(CATCH_WINDOW), "-", max(CATCH_WINDOW))) +
    scale_y_discrete(labels = stock_label) +
    labs(x = "Harvest rate", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank(), legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),
          panel.background = element_rect(fill = NA, colour = "grey85"))
  print(p_freq)
  ggsave(paste0("figures/weak_stocks_prop_years_below_sgen_by_scenario", sfx, ".png"),
         plot = p_freq, width = 16, height = 7, dpi = 600, bg = "white")
  
  # W2: number of stocks below Sgen vs fixed harvest rate
  nred <- filter(n_red_by_hr, bm_period == bm)
  print(nred %>% select(-bm_period, -n_amber, -n_stocks) %>%
          pivot_wider(names_from = scenario, values_from = n_red), n = Inf)
  p_nred <- ggplot(nred, aes(harvest_rate, n_red, colour = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = WEAK_SCENARIO_COLORS, name = NULL) +
    scale_x_continuous(labels = scales::percent, breaks = HR_SWEEP_RATES[c(TRUE, FALSE)]) +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    labs(x = "Fixed harvest rate", y = "Number of stocks below Sgen") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  print(p_nred)
  ggsave(paste0("figures/weak_stocks_n_below_sgen_by_harvest_rate", sfx, ".png"),
         plot = p_nred, width = 9, height = 6, dpi = 600, bg = "white")
  
  # W3: highest tolerated fixed harvest rate, by stock and scenario
  thr <- filter(weak_thresholds, bm_period == bm)
  thr %>%
    mutate(value = case_when(
      outcome == "Below Sgen even at 0%" ~ "red at 0%",
      outcome == "Above Sgen at all rates tested" ~ paste0(">=", round(max(HR_SWEEP_RATES) * 100), "%"),
      TRUE ~ paste0(round(max_tolerated_hr * 100), "%"))) %>%
    select(Stock, scenario, value) %>%
    pivot_wider(names_from = scenario, values_from = value) %>%
    mutate(Stock = factor(Stock, levels = rev(stock_order))) %>%   # weakest first
    arrange(Stock) %>%
    print(n = Inf, width = Inf)
  
  p_thresh <- thr %>%
    mutate(Stock = factor(Stock, levels = rev(stock_order))) %>%
    ggplot(aes(max_tolerated_hr, Stock, colour = scenario, shape = outcome)) +
    geom_point(size = 3, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = WEAK_SCENARIO_COLORS, name = NULL) +
    scale_shape_manual(values = c("Highest rate above Sgen" = 16,
                                  "Above Sgen at all rates tested" = 17,
                                  "Below Sgen even at 0%" = 4), name = NULL) +
    scale_x_continuous(labels = scales::percent, limits = c(0, max(HR_SWEEP_RATES)),
                       breaks = HR_SWEEP_RATES[c(TRUE, FALSE)]) +
    scale_y_discrete(labels = stock_label) +
    labs(x = "Highest fixed harvest rate keeping spawners above Sgen", y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom", legend.box = "vertical",
          panel.grid.minor = element_blank())
  print(p_thresh)
  ggsave(paste0("figures/weak_stocks_max_tolerated_harvest_rate", sfx, ".png"),
         plot = p_thresh, width = 9, height = 8, dpi = 600, bg = "white")
  
  # W4: spawner trajectories, COSEWIC Endangered stocks, with Sgen (dashed)
  walk(WEAK_STOCK_TRAJ_RATES, function(rate) {
    df <- stock_runs %>%
      filter(bm_period == bm, round(harvest_rate * 100) == round(rate * 100),
             Stock %in% COSEWIC_ENDANGERED_STOCKS,
             Year >= HR_SWEEP_START_YEAR, Year <= max(CATCH_WINDOW)) %>%
      mutate(scenario = factor(scenario, levels = WEAK_SCENARIO_LEVELS))
    sgen_df <- df %>% distinct(Stock, Year, Sgen)
    p_traj <- ggplot(df, aes(Year, retroS, colour = scenario)) +
      geom_line(data = sgen_df, aes(Year, Sgen), inherit.aes = FALSE,
                linetype = "dashed", colour = "firebrick", linewidth = 0.5) +
      # historical conditions drawn wider underneath, so it stays visible
      # where a scenario's covariates don't apply to that stock
      geom_line(data = filter(df, scenario == "Historical conditions"), linewidth = 2) +
      geom_line(data = filter(df, scenario != "Historical conditions"), linewidth = 0.7) +
      facet_wrap(~ Stock, scales = "free_y") +
      scale_colour_manual(values = WEAK_SCENARIO_COLORS, name = NULL, drop = FALSE) +
      scale_y_continuous(labels = scales::comma) +
      labs(x = "Year", y = paste0("Spawners at ", round(rate * 100), "% harvest rate")) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    print(p_traj)
    ggsave(paste0("figures/weak_stocks_spawner_trajectories_", round(rate * 100), "pct", sfx, ".png"),
           plot = p_traj, width = 12, height = 8, dpi = 600, bg = "white")
  })
  
  invisible(NULL)
}

walk(names(BM_PERIODS), make_weak_stock_outputs)

# ============================================================
# DENSITY-INDEPENDENT PRODUCTIVITY: observed vs. model vs. scenarios
# ln(R/S) with the density-dependent term removed (+ b*S), by brood year:
#   Observed           : ln(R/S) + b*S      (realized productivity)
#   Model predicted    : a + covariate term (actual covariates)
#   Scenario predicted : a + covariate term (scenario covariates)
# Only stocks linked to at least one scenario are shown; within a stock's
# panel, only the scenarios that affect that stock are drawn.
# ============================================================

PROD_SCENARIO_LABELS <- c(
  "Pinniped scenario" = "Low pinnipeds",
  "SST scenario"      = "Low SST",
  "Pink scenario"     = "No hatchery pinks"
)
PROD_COLORS <- c(
  "Observed"        = "black",
  "Model predicted" = "grey50",
  setNames(unname(SCENARIO_COLORS[names(PROD_SCENARIO_LABELS)]), PROD_SCENARIO_LABELS)
)

prod_actual <- runs_hist %>%
  filter(scenario == "Actual", Year %in% FIT_YEARS) %>%
  transmute(Stock, Year,
            Observed          = lnR_S + rb * AdultEscapement,
            `Model predicted` = ra + cov_term) %>%
  pivot_longer(c(Observed, `Model predicted`), names_to = "series", values_to = "value")

prod_scen <- runs_hist %>%
  filter(scenario %in% names(PROD_SCENARIO_LABELS), Year %in% FIT_YEARS) %>%
  semi_join(filter(scenario_applies, applies), by = c("Stock", "scenario")) %>%
  transmute(Stock, Year,
            series = unname(PROD_SCENARIO_LABELS[scenario]),
            value  = ra + cov_term_proj)

prod_df <- bind_rows(prod_actual, prod_scen) %>%
  filter(Stock %in% unique(prod_scen$Stock), is.finite(value)) %>%
  mutate(series = factor(series, levels = names(PROD_COLORS)),
         Stock  = factor(Stock, levels = rev(stock_order)))   # weakest first

p_prod <- ggplot(prod_df, aes(Year, value, colour = series)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = filter(prod_df, series == "Observed"), size = 0.8, alpha = 0.6) +
  geom_line(data = filter(prod_df, series != "Observed"), linewidth = 0.8) +
  facet_wrap(~ Stock, scales = "free_y", labeller = as_labeller(stock_label)) +
  scale_colour_manual(values = PROD_COLORS, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(
    linetype = c(0, 1, 1, 1, 1), shape = c(16, NA, NA, NA, NA)))) +
  labs(x = "Brood year", y = "Productivity, ln(R/S) + bS") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_prod
ggsave("figures/productivity_observed_model_scenarios_by_stock.png", plot = p_prod,
       width = 14, height = 10, dpi = 600, bg = "white")

PROD_OBS_LABEL <- "Observed"

PROD_COLORS <- c(
  setNames("grey60", PROD_OBS_LABEL),
  "Model predicted" = "black",
  setNames(unname(SCENARIO_COLORS[names(PROD_SCENARIO_LABELS)]), PROD_SCENARIO_LABELS)
)
PROD_LINETYPES <- c(
  setNames("dotted", PROD_OBS_LABEL),
  "Model predicted" = "solid",
  setNames(rep("solid", length(PROD_SCENARIO_LABELS)), PROD_SCENARIO_LABELS)
)

p_prod <- ggplot(prod_df, aes(Year, value, colour = series, linetype = series)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_line(data = filter(prod_df, series == PROD_OBS_LABEL), linewidth = 0.6) +
  geom_line(data = filter(prod_df, series != PROD_OBS_LABEL), linewidth = 0.8) +
  facet_wrap(~ Stock, scales = "free_y", labeller = as_labeller(stock_label)) +
  scale_colour_manual(values = PROD_COLORS, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = PROD_LINETYPES, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.9))) +
  labs(x = "Brood year", y = "Productivity, ln(R/S) + bS") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

p_prod
ggsave("figures/productivity_observed_model_scenarios_by_stock.png", plot = p_prod,
       width = 14, height = 10, dpi = 600, bg = "white")

# ============================================================
# PRODUCTIVITY: observed vs. model vs. scenarios (full Ricker model)
#   Observed           : ln(R/S)
#   Model predicted    : a - b*S + covariate term (actual covariates)
#   Scenario predicted : a - b*S + covariate term (scenario covariates)
# Scenario lines use the OBSERVED spawners, so the gap between a scenario
# line and the model line is the covariate effect alone.
# Only stocks linked to at least one scenario are shown; within a stock's
# panel, only the scenarios that affect that stock are drawn.
# Requires scenario_applies, stock_order, stock_label (WEAK STOCKS section).
# ============================================================

PROD_SCENARIO_LABELS <- c(
  "Pinniped scenario" = "Low pinnipeds",
  "SST scenario"      = "Low SST",
  "Pink scenario"     = "No hatchery pinks"
)
PROD_OBS_LABEL <- "Observed"

PROD_COLORS <- c(
  setNames("grey60", PROD_OBS_LABEL),
  "Model predicted" = "black",
  setNames(unname(SCENARIO_COLORS[names(PROD_SCENARIO_LABELS)]), PROD_SCENARIO_LABELS)
)
PROD_LINETYPES <- c(
  setNames("dashed", PROD_OBS_LABEL),
  "Model predicted" = "solid",
  setNames(rep("solid", length(PROD_SCENARIO_LABELS)), PROD_SCENARIO_LABELS)
)

prod_actual <- runs_hist %>%
  filter(scenario == "Actual", Year %in% FIT_YEARS) %>%
  transmute(Stock, Year,
            obs               = lnR_S,
            `Model predicted` = ra - rb * AdultEscapement + cov_term) %>%
  pivot_longer(c(obs, `Model predicted`), names_to = "series", values_to = "value") %>%
  mutate(series = if_else(series == "obs", PROD_OBS_LABEL, series))

prod_scen <- runs_hist %>%
  filter(scenario %in% names(PROD_SCENARIO_LABELS), Year %in% FIT_YEARS) %>%
  semi_join(filter(scenario_applies, applies), by = c("Stock", "scenario")) %>%
  transmute(Stock, Year,
            series = unname(PROD_SCENARIO_LABELS[scenario]),
            value  = ra - rb * AdultEscapement + cov_term_proj)

prod_df <- bind_rows(prod_actual, prod_scen) %>%
  filter(Stock %in% unique(prod_scen$Stock), is.finite(value)) %>%
  mutate(series = factor(series, levels = names(PROD_COLORS)),
         Stock  = factor(Stock, levels = rev(stock_order)))   # weakest first

p_prod <- ggplot(prod_df, aes(Year, value, colour = series, linetype = series)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_line(data = filter(prod_df, series == PROD_OBS_LABEL), linewidth = 0.6) +
  geom_line(data = filter(prod_df, series != PROD_OBS_LABEL), linewidth = 0.8) +
  facet_wrap(~ Stock, scales = "free_y", labeller = as_labeller(stock_label), ncol=3) +
  scale_colour_manual(values = PROD_COLORS, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = PROD_LINETYPES, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(linewidth = 0.9))) +
  labs(x = "Brood year", y = "ln(R/S)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"))

p_prod
ggsave("figures/productivity_observed_model_scenarios_by_stock.png", plot = p_prod,
       width = 10, height = 14, dpi = 600, bg = "white")

# ============================================================
# PLOTS: productivity and catch lost through time
# ============================================================

plot_productivity_compare <- function(scenario_label, stocks = NULL) {
  df <- runs_retro %>%
    filter(scenario %in% c("Actual", scenario_label)) %>%
    select(Stock, Year, scenario, lnRS = retro_lnRS)
  if (!is.null(stocks)) df <- df %>% filter(Stock %in% stocks)
  
  ggplot(df, aes(Year, lnRS, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1, alpha = 0.6) +
    facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(Actual = "grey30", SCENARIO_COLORS_NO_OBS)) +
    labs(x = "Year", y = "ln(R/S)", color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_catch_lost_cumulative <- function(scenario_label) {
  ggplot(filter(catch_lost_by_stock, scenario == scenario_label), aes(Year, cum_catch_lost)) +
    geom_area(fill = SCENARIO_COLORS[[scenario_label]], alpha = 0.2) +
    geom_line(linewidth = 1, color = SCENARIO_COLORS[[scenario_label]]) +
    facet_wrap(~ Stock, scales = "free_y") +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Year", y = "Cumulative catch lost") +
    theme_minimal()
}

walk(names(SCENARIOS), ~ print(plot_productivity_compare(.x)))
walk(names(SCENARIOS), ~ print(plot_catch_lost_cumulative(.x)))

p_cum_lost <- ggplot(catch_lost_totals, aes(Year, cum_catch_lost, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = SCENARIO_COLORS_NO_OBS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost", color = "Scenario driver") +
  theme_minimal()
p_cum_lost
ggsave("figures/catch_lost_by_scenario_driver.png", plot = p_cum_lost,
       width = 10, height = 6, dpi = 600, bg = "white")

# ============================================================
# AVERAGE YEARLY CATCH LOST over CATCH_WINDOW (historical harvest rate)
# ============================================================

catch_lost_recent <- catch_lost_totals %>% filter(Year %in% CATCH_WINDOW)

catch_lost_recent_summary <- catch_lost_recent %>%
  group_by(scenario) %>%
  summarise(
    mean_catch_lost = mean(catch_lost),
    se_catch_lost   = sd(catch_lost) / sqrt(n()),
    .groups = "drop"
  )

msy_hlines <- geom_hline(data = MSY_LINES, aes(yintercept = value, linetype = label),
                         color = "black", linewidth = 0.7)

p1 <- ggplot(catch_lost_recent_summary,
             aes(x = reorder(scenario, -mean_catch_lost), y = mean_catch_lost, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_catch_lost - se_catch_lost,
                    ymax = mean_catch_lost + se_catch_lost), width = 0.15) +
  msy_hlines +
  scale_fill_manual(values = SCENARIO_COLORS_NO_OBS) +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = NULL, values = MSY_LINETYPES) +
  labs(x = NULL, y = "Mean yearly catch lost", fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/mean_catch_lost_by_scenario_2000-present_v2.png", plot = p1,
       width = 8, height = 5.5, dpi = 600, bg = "white")

p2 <- ggplot(catch_lost_recent,
             aes(x = reorder(scenario, catch_lost, FUN = median), y = catch_lost, fill = scenario)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  msy_hlines +
  scale_fill_manual(values = SCENARIO_COLORS_NO_OBS) +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual(name = NULL, values = MSY_LINETYPES) +
  labs(x = NULL, y = "Yearly catch lost", fill = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/catch_lost_boxplot_by_scenario_2000-present_v2.png", plot = p2,
       width = 8, height = 5.5, dpi = 600, bg = "white")

# Note: negative values are dropped under the sqrt transform
p2.2 <- p2 +
  scale_y_continuous(labels = scales::comma, trans = scales::sqrt_trans(),
                     breaks = scales::breaks_pretty(n = 8)) +
  labs(y = "Yearly catch lost (sqrt transformed)") +
  theme(legend.position = "right")

p_box_pair <- p2 | p2.2
p_box_pair
ggsave("figures/catch_lost_sqrt2_v2.png", plot = p_box_pair,
       width = 11.5, height = 5.5, dpi = 600, bg = "white")

# ============================================================
# RETURN TRAJECTORIES: observed vs. all scenarios
# ============================================================

returns_total <- returns_by_stock %>%
  sum_complete_stocks("Return", scenario, Year) %>%
  rename(Return = total) %>%
  mutate(scenario = factor(scenario, levels = names(SCENARIO_COLORS)))

p_ret_total <- ggplot(mapping = aes(Year, Return, color = scenario)) +
  geom_line(data = filter(returns_total, scenario == "Observed"), linewidth = 1.3) +
  geom_line(data = filter(returns_total, scenario != "Observed"), linewidth = 1, alpha = 0.6) +
  scale_color_manual(values = SCENARIO_COLORS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Total return (all stocks)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
p_ret_total
ggsave("figures/return_trajectories_all_scenarios_total.png", plot = p_ret_total,
       width = 10, height = 6, dpi = 600, bg = "white")

returns_by_stock_plot <- returns_by_stock %>%
  mutate(scenario = factor(scenario, levels = names(SCENARIO_COLORS)))

p_ret_stock <- ggplot(mapping = aes(Year, Return, color = scenario)) +
  geom_line(data = filter(returns_by_stock_plot, scenario == "Observed"), linewidth = 1) +
  geom_line(data = filter(returns_by_stock_plot, scenario != "Observed"), linewidth = 0.7, alpha = 0.7) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Return", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")
p_ret_stock
ggsave("figures/return_trajectories_all_scenarios_by_stock.png", plot = p_ret_stock,
       width = 14, height = 10, dpi = 600, bg = "white")

# ============================================================
# COSEWIC-STYLE STATUS CLASSIFICATION
# Approximate re-application of COSEWIC's quantitative thresholds
# (Criteria A/C/D) -- comparative indicators across scenarios, not a
# formal reassessment.
# ============================================================

# Most severe status implied by any criterion. Special Concern isn't
# distinguished from Not at Risk (qualitative in the source report).
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

status_all <- returns_by_stock %>%
  left_join(GT_TABLE, by = "Stock") %>%
  arrange(scenario, Stock, Year) %>%
  group_by(scenario, Stock) %>%
  mutate(
    gen_mean          = trailing_mean(Return, width = first(GT)),
    gen_mean_3gen_ago = lag(gen_mean, n = 3 * first(GT)),
    decline_pct       = 1 - gen_mean / gen_mean_3gen_ago
  ) %>%
  ungroup() %>%
  mutate(status = classify_status(decline_pct, gen_mean))

# Capped at MAX_ASSESSMENT_YEAR for a contemporaneous comparison ("given
# what was known as of the 2017 report, would this driver have changed the
# designation"). The uncapped version answers "status today".
latest_status <- function(df) {
  df %>%
    filter(!is.na(status)) %>%
    group_by(Stock, scenario) %>%
    slice_max(Year, n = 1) %>%
    ungroup() %>%
    select(Stock, scenario, Year, gen_mean, decline_pct, status)
}
status_summary_latest          <- latest_status(filter(status_all, Year <= MAX_ASSESSMENT_YEAR))
status_summary_latest_uncapped <- latest_status(status_all)

# (COSEWIC_ENDANGERED_STOCKS is defined in CONFIG.)

status_comparison_table <- status_summary_latest %>%
  filter(Stock %in% COSEWIC_ENDANGERED_STOCKS) %>%
  select(Stock, scenario, status) %>%
  pivot_wider(names_from = scenario, values_from = status) %>%
  arrange(Stock)

print(status_comparison_table)

# ============================================================
# % ABUNDANCE INCREASE OVER LOW-ABUNDANCE PERIODS, BY SCENARIO
# ============================================================

pct_increase_by_year <- returns_by_stock %>%
  filter(scenario != "Observed") %>%
  rename(Return_scenario = Return) %>%
  left_join(returns_by_stock %>% filter(scenario == "Observed") %>%
              select(Stock, Year, Return_observed = Return),
            by = c("Stock", "Year")) %>%
  mutate(pct_increase = 100 * (Return_scenario - Return_observed) / Return_observed)

pct_increase_over_low_period <- returns_by_stock %>%
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
    pct_increase_over_low_period %>%
      select(Stock, scenario, pct_increase) %>%
      pivot_wider(names_from = scenario, values_from = pct_increase, names_prefix = "pct_increase_"),
    by = "Stock"
  )

print(recovery_summary)

p_pct_low <- ggplot(pct_increase_over_low_period,
                    aes(x = reorder(Stock, -pct_increase), y = pct_increase, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS_NO_OBS) +
  labs(x = NULL, y = "% increase in mean abundance over low-abundance period",
       fill = "Scenario driver") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
p_pct_low
ggsave("figures/pct_increase_over_low_period_by_scenario.png", plot = p_pct_low,
       width = 10, height = 6, dpi = 600, bg = "white")

p_pct_time <- ggplot(pct_increase_by_year, aes(Year, pct_increase, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.25) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS_NO_OBS) +
  labs(x = "Year", y = "% increase in abundance vs. observed", color = "Scenario driver") +
  theme_minimal() +
  theme(legend.position = "bottom")
p_pct_time
ggsave("figures/pct_increase_over_time_by_scenario.png", plot = p_pct_time,
       width = 14, height = 10, dpi = 600, bg = "white")

status_plot_df <- status_summary_latest %>%
  filter(Stock %in% COSEWIC_ENDANGERED_STOCKS) %>%
  mutate(
    scenario = factor(scenario, levels = names(SCENARIO_COLORS)),
    status   = factor(status, levels = c("Endangered", "Threatened", "Not at Risk / Special Concern"))
  )

p_status <- ggplot(status_plot_df, aes(scenario, Stock, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "Endangered" = "#B22222",
    "Threatened" = "#E8A33D",
    "Not at Risk / Special Concern" = "#4C9A5B"
  )) +
  labs(x = NULL, y = NULL, fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
p_status
ggsave("figures/status_by_scenario_heatmap.png", plot = p_status,
       width = 9, height = 6, dpi = 600, bg = "white")

# ============================================================
# DIAGNOSTICS
#
# diagnose_stock: spawners/return/productivity trajectory for one stock
# across scenarios, plus where the pinniped freeze-year value sits in the
# historical covariate range (extrapolation risk). Useful when a stock's
# status worsens under a scenario despite a favourable coefficient: the
# model is recursive, so a productivity boost raises escapement, which
# feeds density dependence next generation; a stock with large |rb|
# relative to ra can overshoot and crash (Ricker overcompensation). Also
# check whether slice_max(Year) lands on the down-swing of such a cycle.
#
# inspect_status_trajectory: year-by-year gen_mean/decline_pct/status.
# ============================================================

diagnose_stock <- function(stock_name) {
  
  traj <- runs_retro %>%
    filter(Stock == stock_name) %>%
    mutate(scenario = if_else(scenario == "Actual", "Observed", scenario))
  diag_colors <- scale_color_manual(values = SCENARIO_COLORS)
  
  p_s <- ggplot(traj, aes(Year, retroS, color = scenario)) +
    geom_line(linewidth = 0.8) + diag_colors +
    labs(y = "Spawners (retroS)") +
    theme_minimal()
  p_r <- ggplot(traj, aes(Year, retroR, color = scenario)) +
    geom_line(linewidth = 0.8) + diag_colors +
    labs(y = "Return (retroR)") +
    theme_minimal()
  p_p <- ggplot(traj, aes(Year, retro_lnRS, color = scenario)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    diag_colors +
    labs(y = "ln(R/S)") +
    theme_minimal()
  
  p_diag <- (p_s / p_r / p_p) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  print(p_diag)
  ggsave(paste0("figures/diagnostic_", tolower(gsub(" ", "_", stock_name)), "_trajectories.png"),
         plot = p_diag, width = 9, height = 10, dpi = 600, bg = "white")
  
  covariates_main %>%
    filter(Stock == stock_name) %>%
    summarise(
      SeaLions_hist_min  = min(SeaLions, na.rm = TRUE),
      SeaLions_hist_max  = max(SeaLions, na.rm = TRUE),
      SeaLions_at_freeze = SeaLions[Year == FREEZE_YEAR],
      seal_hist_min      = min(seal, na.rm = TRUE),
      seal_hist_max      = max(seal, na.rm = TRUE),
      seal_at_freeze     = seal[Year == FREEZE_YEAR]
    ) %>%
    print()
  
  invisible(traj)
}

inspect_status_trajectory <- function(stock_name, scenario_label = "Observed") {
  status_all %>%
    filter(Stock == stock_name, scenario == scenario_label) %>%
    select(Year, gen_mean, decline_pct, status) %>%
    print(n = Inf)
}

diagnose_stock("Raft")
inspect_status_trajectory("Late Stuart")