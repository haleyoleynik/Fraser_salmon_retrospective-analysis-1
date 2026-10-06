# Chum / steelhead retrospective model — scenario comparisons
# Haley Oleynik, Murdoch McAllister
#
# Merged version, combining two parallel copies of this script into one.
# Keeps the three core scenario analyses:
#   1. Historic harvest rate: SSL control vs. no control
#   2. Harvest-rate sweep (byrate held fixed): SSL control vs. no control
#   3. Bycatch-rate sweep (harvest held at historic): SSL control vs. no control
# ...plus chum catch-lost-to-pinnipeds averaging, 2000-present, and a
# steelhead low-abundance-period % increase summary (pinniped-control
# scenario), styled to match the sockeye retrospective bar chart.
#
# Dropped from the original: dead/commented-out code, duplicate function
# definitions, the SSL-predation-rate exploration section, the
# lost-value/CPI section, and the standalone COSEWIC status-replication
# section (validation against the 2020 assessment's ~82%/~80% decline
# figures) — the latter lives in a separate copy of this script if needed
# again; it was dropped here in favor of the low-abundance-period summary
# below, which isn't meant to coexist with it in the same run.
#
# FIXES applied (both now folded into this version permanently):
#   1. run_thompson_scenario() / run_chilcotin_scenario() / run_chum_scenario():
#      the observed-residual calibration term (exp(ln_obs_pred)) now falls
#      back to 1 (i.e. pure alpha-model prediction) whenever the observed
#      value for that year is NA, instead of propagating NA through
#      recruits_alt and every downstream age-class/lag that depends on it.
#      This was producing the disconnected gaps in the chum/Thompson/
#      Chilcotin productivity plots.
#   2. plot_scenario_comparison() / plot_scenario_sweep(): now actually
#      pass `title` into labs() — it was accepted as an argument but never
#      used, so plots were rendering without titles.
#   3. Added the missing chum_sweep_with_obs construction (model + observed,
#      pivoted long) ahead of the "recruits per spawner & returns across
#      harvest rates, with historical observed overlaid" section — it was
#      referenced but never built, which would error out and stop
#      everything after it (including the whole steelhead bycatch-grid
#      section) from running.
#   4. run_thompson_scenario() / run_chilcotin_scenario(): sh_*_spawners and
#      sh_*_recruits are already reported in thousands in the input data
#      (e.g. 1.67 = 1,670 fish). sh_*_model_recruits and sh_*_ln_obs_pred
#      were dividing them by 1000 a second time, corrupting the residual-
#      calibration term from the very first modeled brood year. That error
#      compounded through the recursive age-lag structure and, combined
#      with FN/sport catch being subtracted as fixed historical headcounts
#      from an increasingly understated modeled return, drove escapement to
#      exactly 0 in many years -- the source of the disconnected segments
#      in the Thompson/Chilcotin "No SSL control" productivity/returns
#      lines. Removed the extra /1000.
#   5. run_chum_scenario() / run_thompson_scenario() / run_chilcotin_scenario():
#      "No SSL control" with harvest reconstructed to match history
#      (U_historic == 1, and for steelhead byrate == default_byrate) has no
#      counterfactual left to simulate -- it's meant to reproduce the
#      observed record exactly. Rather than rely on the recursive
#      Ricker/lag simulation (which can still drift from minor residual/
#      interpolation effects even with fix #4 above), these cases now pull
#      recruits/spawners/returns/catch directly from the observed columns.
#      The SSL-control counterfactual and the harvest/bycatch sweeps
#      (U_historic == 0, or byrate swept away from default) are genuine
#      hypotheticals with no historical record to copy, so they still run
#      the full simulation.
#   6. Replaced the single shared MODEL_YEARS plotting window (1980-2016,
#      which truncated chum and Chilcotin to Thompson's shorter record)
#      with per-stock ranges: chum 1951-2016, Thompson 1978-2016, Chilcotin
#      1973-2016 (each stock's own observed-data start; 2016 kept as the
#      shared upper bound since covariates only run through then). Faceted
#      plots now use scales = "free" so each stock's panel shows its own
#      x-axis range instead of being padded to match the others.
#   7. Steelhead recruitment (Thompson + Chilcotin) no longer multiplies by
#      the observed residual exp(ln(Robs/Rpred)). The spreadsheet
#      ('Retrospective Run1' CQ / EU) uses the pure Ricker prediction for
#      steelhead; only CHUM recruits carry the residual (AW: *EXP(N)). With
#      the residual, year-specific noise was fed back through spawners and
#      the age lags, roughly doubling peak recruits and year-to-year CV
#      under SSL control. Toggle with APPLY_SH_RESIDUALS.
#   8. SSL control year is now one setting (SSL_CONTROL_YEAR = 1978) used by
#      all three stocks, the predation block and the catch-lost summaries.
#      Chilcotin previously switched to the 1978 control value from 1974,
#      while its actual SSLs were still in use through 1978 (the spreadsheet
#      applies control from 1979 on, EW42); it now matches the spreadsheet.
#   9. Added the predation-mortality accounting block (spreadsheet cols
#      FB:GL, v12): modeled scenario spawners, scenario SSL, rundays = 30,
#      baseline = zero SSLs, as in the spreadsheet (PRED_BASELINE).

library(tidyverse)
library(readr)
library(patchwork)
library(scales)

# ============================================================
# CONFIG
# ============================================================

MODEL_YEARS <- 1980:2016   # legacy shared window; kept as fallback default only
harvest_rates <- seq(0, 0.8, by = 0.1)
bycatch_rates <- seq(0, 1, by = 0.1)
default_byrate <- 0.69     # current bycatch/FN-mortality proxy rate for steelhead

# Sea lion control counterfactual: SSLs held at this brood year's covariate
# value for every later year (spreadsheet rows 41: $CZ$41, $CS$41, $EW$41).
SSL_CONTROL_YEAR <- 1978

# Standardized SSL covariate (z) at SSL_CONTROL_YEAR for each stock. NA =
# look it up from the input data (Year == SSL_CONTROL_YEAR). Spreadsheet
# values for 1978: chum -0.6884925755025276, Thompson -0.904618196,
# Chilcotin -0.877048729. Only needed if SSL_CONTROL_YEAR is moved to a year
# the input data don't cover.
SSL_CONTROL_Z <- c(chum = NA_real_, thompson = NA_real_, chilcotin = NA_real_)

# Spreadsheet uses pure Ricker predictions for steelhead (no obs residual).
APPLY_SH_RESIDUALS <- FALSE

get_control_ssl <- function(df, col, stock) {
  if (!is.na(SSL_CONTROL_Z[[stock]])) return(unname(SSL_CONTROL_Z[[stock]]))
  val <- df %>% filter(Year == SSL_CONTROL_YEAR) %>% pull(all_of(col)) %>% as.numeric()
  val <- val[!is.na(val)]
  if (length(val) == 0) {
    stop(sprintf(paste0("No %s value for %d in the input data. Set ",
                        "SSL_CONTROL_Z[\"%s\"] in CONFIG to the %d SSL z-score."),
                 col, SSL_CONTROL_YEAR, stock, SSL_CONTROL_YEAR), call. = FALSE)
  }
  val[1]
}

# Per-stock plotting windows -- each stock's own reliable data range, rather
# than forcing every panel to the shortest-history stock's window. Chum
# observed data starts 1951 (covariates run only through 2016, capping the
# upper end for all three); Thompson steelhead observed data starts 1978;
# Chilcotin steelhead observed data starts 1973.
CHUM_YEARS      <- 1951:2016
THOMPSON_YEARS  <- 1978:2016
CHILCOTIN_YEARS <- 1973:2016

stock_year_bounds <- tibble::tibble(
  Stock  = c("Chum", "Thompson steelhead", "Chilcotin steelhead"),
  yr_min = c(min(CHUM_YEARS), min(THOMPSON_YEARS), min(CHILCOTIN_YEARS)),
  yr_max = c(max(CHUM_YEARS), max(THOMPSON_YEARS), max(CHILCOTIN_YEARS))
)

## Filters a long-format df (must have Stock + Year columns) down to each
## stock's own reliable year range, instead of one shared range for all.
filter_stock_years <- function(df) {
  df %>%
    left_join(stock_year_bounds, by = "Stock") %>%
    filter(Year >= yr_min, Year <= yr_max) %>%
    select(-yr_min, -yr_max)
}

# ============================================================
# LOAD DATA
# ============================================================

data       <- read_csv("R/Chum Steelhead Retrospective Shiny App/s-r_data.csv") %>% select(-SSL)
sh_data    <- read_csv("R/Chum Steelhead Retrospective Shiny App/sh_s-r_data.csv")
covariates <- read_csv("R/Chum Steelhead Retrospective Shiny App/covariates.csv")

# ============================================================
# MODEL FUNCTIONS
# ============================================================

## --- Chum -----------------------------------------------------------------

run_chum_scenario <- function(data, covariates, U_apply, SSL_control,
                              U_historic = 0,
                              start_year = 1978,
                              intercept = 1.03737862843252,
                              pdo_adult_coef = 0.0929480015696915,
                              npgo_coef = 0.102617626088713,
                              pdo_smolt_coef = -0.105950753725792,
                              ssl_coef = -0.224246284912103,
                              spawners_coef = -4.95136622626478E-07) {
  
  new.data <- data %>%
    left_join(covariates, by = "Year") %>%
    arrange(Year) %>%
    mutate(
      catch = chum_total_stock - chum_spawners,
      U_chum = catch / chum_total_stock,
      chum_base_alpha = intercept + (PDO_adult * pdo_adult_coef + NPGO * npgo_coef +
                                       PDO_smolt * pdo_smolt_coef + SSL * ssl_coef),
      chum_model_recruits = chum_spawners * exp(chum_base_alpha + chum_spawners * spawners_coef),
      chum_ln_obs_pred = log(chum_recruits_obs / chum_model_recruits),
      Nage3_obs = chum_total_stock * prop3,
      Nage4_obs = chum_total_stock * prop4,
      Nage5_obs = chum_total_stock * prop5,
      Nage6_obs = chum_total_stock * prop6,
      Nage3_pred = case_when(
        Year <= 1953 ~ Nage3_obs,
        Year >= 1954 ~ lag(chum_model_recruits, 3) * prop3 * exp(lag(chum_ln_obs_pred, 3))),
      Nage4_pred = case_when(
        Year <= 1954 ~ Nage4_obs,
        Year >= 1955 ~ lag(chum_model_recruits, 4) * prop4 * exp(lag(chum_ln_obs_pred, 4))),
      Nage5_pred = case_when(
        Year <= 1955 ~ Nage5_obs,
        Year >= 1956 ~ lag(chum_model_recruits, 5) * prop5 * exp(lag(chum_ln_obs_pred, 5))),
      Nage6_pred = case_when(
        Year <= 1956 ~ Nage6_obs,
        Year >= 1957 ~ lag(chum_model_recruits, 6) * prop6 * exp(lag(chum_ln_obs_pred, 6))),
      recruits_pred = rowSums(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)), na.rm = TRUE),
      U_chum_pred = catch / recruits_pred,
      chum_commercial_harvest_uapply = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 & U_historic == 0 ~ U_apply,
        Year >= 1991 & U_historic == 1 ~ U_chum_pred)
    )
  
  SSL_ctrl <- get_control_ssl(new.data, "SSL", "chum")
  
  df <- new.data %>%
    arrange(Year) %>%
    mutate(
      chum_recruits_alt = NA_real_,
      Nage3_alt = Nage3_pred, Nage4_alt = Nage4_pred,
      Nage5_alt = Nage5_pred, Nage6_alt = Nage6_pred,
      sum_alt = NA_real_,
      catch_alt = NA_real_,
      chum_spawners_pred = chum_spawners,
      chum_SSL_alt = NA_real_,
      chum_SSL_alpha = NA_real_
    )
  
  for (i in seq_len(nrow(df))) {
    
    df$chum_SSL_alt[i] <- if (df$Year[i] <= SSL_CONTROL_YEAR) {
      df$SSL[i]
    } else {
      (1 - SSL_control) * df$SSL[i] + SSL_control * SSL_ctrl
    }
    
    df$chum_SSL_alpha[i] <-
      intercept + df$PDO_adult[i] * pdo_adult_coef + df$NPGO[i] * npgo_coef +
      df$PDO_smolt[i] * pdo_smolt_coef + df$chum_SSL_alt[i] * ssl_coef
    
    # age structure from PAST recruits_alt -- can run before this year's
    # recruits are computed
    if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
    if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
    if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
    if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
    
    df$sum_alt[i] <- sum(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i], na.rm = TRUE)
    
    # apply the scenario harvest rate -> catch & escapement, BEFORE this
    # year's recruits are calculated
    df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
    df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i])
    
    # FIX: fall back to the pure alpha-model prediction (residual adjustment
    # = 1) when the observed value for this year is missing, so a gap in
    # chum_recruits_obs can't inject an NA that then propagates through
    # future age-classes via the lag(recruits_alt) terms above.
    resid_adj <- if (is.na(df$chum_ln_obs_pred[i])) 1 else exp(df$chum_ln_obs_pred[i])
    
    df$chum_recruits_alt[i] <-
      if (SSL_control == 0) {
        df$chum_spawners_pred[i] *
          exp(df$chum_base_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
          resid_adj
      } else {
        df$chum_spawners_pred[i] *
          exp(df$chum_SSL_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
          resid_adj
      }
  }
  
  # FIX: "No SSL control" with harvest reconstructed to match history
  # (U_historic == 1) is meant to reproduce the actual observed record --
  # there's no counterfactual left to simulate once SSL_control == 0 and
  # harvest is historic. Rather than rely on the recursive Ricker/lag
  # simulation to land on the same numbers (it can drift from small
  # residual/interpolation effects even when correctly specified), pull
  # the key outputs directly from the observed columns so this scenario
  # is exact. The SSL-control counterfactual (SSL_control == 1) and the
  # harvest/bycatch sweeps (U_historic == 0) still need the simulation,
  # since those are genuine hypotheticals with no historical record to
  # copy.
  if (SSL_control == 0 && U_historic == 1) {
    df <- df %>%
      mutate(
        chum_recruits_alt   = chum_recruits_obs,
        chum_spawners_pred  = chum_spawners,
        catch_alt           = catch
      )
  }
  
  df %>%
    mutate(
      ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control")
    )
}

## --- Steelhead: Thompson ---------------------------------------------------

run_thompson_scenario <- function(all_data, SSL_control, U_historic, byrate,
                                  start_year = 1978,
                                  sh_thompson_intercept      = 1.572107637,
                                  sh_thompson_sst_coef       = -0.203463091,
                                  sh_thompson_ssl_coef       = -0.764277677,
                                  sh_thompson_npgo_coef      = -0.0402,
                                  sh_thompson_spawners_coef  = -0.804438793) {
  
  sh_thompson_SSL_ctrl <- get_control_ssl(all_data, "sh_thompson_SL", "thompson")
  FN_thompson_2018 <- all_data %>% filter(Year == 2018) %>% pull(sh_thompson_FN_mortalities) %>% as.numeric()
  
  df <- all_data %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_thompson_base_alpha = sh_thompson_intercept +
        sh_thompson_SST  * sh_thompson_sst_coef +
        sh_thompson_SL   * sh_thompson_ssl_coef +
        sh_thompson_NPGO * sh_thompson_npgo_coef,
      # FIX: sh_thompson_spawners/sh_thompson_recruits are already reported in
      # thousands in the input data (e.g. 1.67 = 1,670 fish) -- dividing by
      # 1000 again here silently corrupted the residual-calibration term
      # (sh_thompson_ln_obs_pred), understating recruits_alt from the very
      # first modeled brood year and compounding through the recursive
      # age-lag structure. This is the primary driver of the "No SSL
      # control" line collapsing to disconnected zero-runs in the
      # productivity/returns plots.
      sh_thompson_model_recruits =
        sh_thompson_spawners *
        exp(sh_thompson_base_alpha +
              sh_thompson_spawners * sh_thompson_spawners_coef),
      sh_thompson_ln_obs_pred =
        log(sh_thompson_recruits / sh_thompson_model_recruits),
      sh_thompson_pred_bycatch = sh_thompson_prefishery_N - sh_thompson_sport_mortalities -
        sh_thompson_FN_mortalities - 1000 * sh_thompson_spawners,
      sh_thompson_U = sh_thompson_pred_bycatch / sh_thompson_prefishery_N,
      sh_thompson_recruits_alt      = NA_real_,
      sh_thompson_spawners_pred     = sh_thompson_spawners,
      sh_thompson_Nage4_pred = NA_real_, sh_thompson_Nage5_pred = NA_real_,
      sh_thompson_Nage6_pred = NA_real_, sh_thompson_Nage7_pred = NA_real_,
      sh_thompson_Nage8_pred = NA_real_,
      sh_thompson_SSL_alt           = NA_real_,
      sh_thompson_SSL_alpha         = NA_real_,
      sh_thompson_alpha_CN          = NA_real_,
      sh_thompson_sum_pred          = NA_real_,
      sh_thompson_bycatch_pred      = NA_real_,
      sh_thompson_FN_catch_pred     = NA_real_,
      sh_thompson_total_catch_pred  = NA_real_,
      sh_thompson_U_comm            = NA_real_
    )
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_thompson_U_comm[i] <-
      if (yr <= 1990) {
        df$sh_thompson_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_thompson_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_thompson_SSL_alt[i] <- if (yr <= SSL_CONTROL_YEAR) {
      df$sh_thompson_SL[i]
    } else {
      (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_ctrl
    }
    
    df$sh_thompson_SSL_alpha[i] <-
      sh_thompson_intercept +
      df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
      df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
      df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    df$sh_thompson_alpha_CN[i] <-
      sh_thompson_intercept +
      df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
      df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
      df$sh_thompson_SL[i]   * sh_thompson_ssl_coef
    
    # within-year fixed-point iteration to resolve the spawner/catch circularity
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      # FIX: same NA-residual fallback as chum -- a missing observed value
      # in a given year (real gaps in the Thompson survey record) used to
      # set recruits_alt to NA for that year, which then propagated
      # through every downstream age-class (ages 4-8, i.e. up to 5 years
      # of lags) via lag_recruits() below. That's what produced the
      # disconnected multi-year gaps in the productivity plot.
      # FIX 7: spreadsheet CQ has no residual term for steelhead
      resid_adj <- if (!APPLY_SH_RESIDUALS || is.na(df$sh_thompson_ln_obs_pred[i])) 1 else exp(df$sh_thompson_ln_obs_pred[i])
      
      df$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_thompson_alpha_CN[i] +
                      sh_thompson_spawners_coef * spk) *
            resid_adj
        } else {
          spk * exp(df$sh_thompson_SSL_alpha[i] +
                      sh_thompson_spawners_coef * spk) *
            resid_adj
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df$sh_thompson_recruits_alt[j]
      }
      
      df$sh_thompson_Nage4_pred[i] <- if (yr < start_year + 4) df$sh_thompson_prefishery_N[i] * df$sh_thompson_p4[i] else lag_recruits(yr - 4) * df$sh_thompson_p4[i] * 1000
      df$sh_thompson_Nage5_pred[i] <- if (yr < start_year + 5) df$sh_thompson_prefishery_N[i] * df$sh_thompson_p5[i] else lag_recruits(yr - 5) * df$sh_thompson_p5[i] * 1000
      df$sh_thompson_Nage6_pred[i] <- if (yr < start_year + 6) df$sh_thompson_prefishery_N[i] * df$sh_thompson_p6[i] else lag_recruits(yr - 6) * df$sh_thompson_p6[i] * 1000
      df$sh_thompson_Nage7_pred[i] <- if (yr < start_year + 7) df$sh_thompson_prefishery_N[i] * df$sh_thompson_p7[i] else lag_recruits(yr - 7) * df$sh_thompson_p7[i] * 1000
      df$sh_thompson_Nage8_pred[i] <- if (yr < start_year + 8) df$sh_thompson_prefishery_N[i] * df$sh_thompson_p8[i] else lag_recruits(yr - 8) * df$sh_thompson_p8[i] * 1000
      
      df$sh_thompson_sum_pred[i] <- sum(df$sh_thompson_Nage4_pred[i], df$sh_thompson_Nage5_pred[i],
                                        df$sh_thompson_Nage6_pred[i], df$sh_thompson_Nage7_pred[i],
                                        df$sh_thompson_Nage8_pred[i], na.rm = TRUE)
      
      df$sh_thompson_bycatch_pred[i] <- df$sh_thompson_sum_pred[i] * df$sh_thompson_U_comm[i]
      
      if (yr <= 2018) {
        df$sh_thompson_FN_catch_pred[i] <- df$sh_thompson_FN_mortalities[i]
      } else {
        denom_2018 <- df$sh_thompson_sum_pred[df$Year == 2018] - df$sh_thompson_bycatch_pred[df$Year == 2018]
        df$sh_thompson_FN_catch_pred[i] <- FN_thompson_2018 / denom_2018 *
          (df$sh_thompson_sum_pred[i] - df$sh_thompson_bycatch_pred[i])
      }
      
      df$sh_thompson_total_catch_pred[i] <- df$sh_thompson_FN_catch_pred[i] +
        df$sh_thompson_sport_mortalities[i] + df$sh_thompson_bycatch_pred[i]
      
      S_new <- max(df$sh_thompson_sum_pred[i] - df$sh_thompson_total_catch_pred[i], 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= 1e-8 * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    
    df$sh_thompson_spawners_pred[i] <- S_old
  }
  
  # FIX: same rationale as chum above. "No SSL control" with harvest
  # reconstructed to match history (U_historic == 1) AND byrate at its
  # default (i.e. not being swept away from the historic proxy) has no
  # counterfactual left -- pull straight from the observed record instead
  # of trusting the recursive fixed-point simulation to reproduce it
  # exactly. sh_thompson_pred_bycatch/U are already the observed bycatch/
  # exploitation-rate reconstruction computed above, so this is an exact
  # accounting identity: sum_pred - total_catch_pred = spawners_pred.
  if (SSL_control == 0 && U_historic == 1 && byrate == default_byrate) {
    df <- df %>%
      mutate(
        sh_thompson_recruits_alt     = sh_thompson_recruits,
        sh_thompson_spawners_pred    = sh_thompson_spawners * 1000,
        sh_thompson_sum_pred         = sh_thompson_prefishery_N,
        sh_thompson_bycatch_pred     = sh_thompson_pred_bycatch,
        sh_thompson_FN_catch_pred    = sh_thompson_FN_mortalities,
        sh_thompson_total_catch_pred = sh_thompson_FN_mortalities +
          sh_thompson_sport_mortalities + sh_thompson_pred_bycatch,
        sh_thompson_U_comm           = sh_thompson_U
      )
  }
  
  df
}

## --- Steelhead: Chilcotin ---------------------------------------------------

run_chilcotin_scenario <- function(df, SSL_control, U_historic, byrate,
                                   start_year = 1973,
                                   sh_chilcotin_intercept = 1.053608979,
                                   sh_chilcotin_sst_coef = -0.127949278,
                                   sh_chilcotin_ssl_coef = -0.792741195,
                                   sh_chilcotin_npgo_coef = 0.152526045,
                                   sh_chilcotin_pdo_coef = 0.202708011,
                                   sh_chilcotin_spawners_coef = -1.022467631) {
  
  sh_chilcotin_SSL_ctrl <- get_control_ssl(df, "sh_chilcotin_SL", "chilcotin")
  FN_chilcotin_2018 <- df %>% filter(Year == 2018) %>% pull(sh_chilcotin_FN_mortalities) %>% as.numeric()
  
  df <- df %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_chilcotin_base_alpha = sh_chilcotin_intercept +
        sh_chilcotin_SST  * sh_chilcotin_sst_coef +
        sh_chilcotin_SL   * sh_chilcotin_ssl_coef +
        sh_chilcotin_NPGO * sh_chilcotin_npgo_coef +
        sh_chilcotin_PDO * sh_chilcotin_pdo_coef,
      # FIX: same double-unit-division bug as Thompson above -- sh_chilcotin_
      # spawners/recruits are already in thousands, so the extra /1000 here
      # corrupted sh_chilcotin_ln_obs_pred the same way.
      sh_chilcotin_model_recruits =
        sh_chilcotin_spawners *
        exp(sh_chilcotin_base_alpha +
              sh_chilcotin_spawners * sh_chilcotin_spawners_coef),
      sh_chilcotin_ln_obs_pred =
        log(sh_chilcotin_recruits / sh_chilcotin_model_recruits),
      sh_chilcotin_pred_bycatch = sh_chilcotin_prefishery_N - sh_chilcotin_sport_mortalities -
        sh_chilcotin_FN_mortalities - 1000 * sh_chilcotin_spawners,
      sh_chilcotin_U = sh_chilcotin_pred_bycatch / sh_chilcotin_prefishery_N,
      sh_chilcotin_recruits_alt      = NA_real_,
      sh_chilcotin_spawners_pred     = sh_chilcotin_spawners,
      sh_chilcotin_Nage4_pred = NA_real_, sh_chilcotin_Nage5_pred = NA_real_,
      sh_chilcotin_Nage6_pred = NA_real_, sh_chilcotin_Nage7_pred = NA_real_,
      sh_chilcotin_Nage8_pred = NA_real_,
      sh_chilcotin_SSL_alt           = NA_real_,
      sh_chilcotin_SSL_alpha         = NA_real_,
      sh_chilcotin_alpha_CN          = NA_real_,
      sh_chilcotin_sum_pred          = NA_real_,
      sh_chilcotin_bycatch_pred      = NA_real_,
      sh_chilcotin_FN_catch_pred     = NA_real_,
      sh_chilcotin_total_catch_pred  = NA_real_,
      sh_chilcotin_U_comm            = NA_real_
    )
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_chilcotin_U_comm[i] <-
      if (yr <= 1990) {
        df$sh_chilcotin_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_chilcotin_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= SSL_CONTROL_YEAR) {
      df$sh_chilcotin_SL[i]
    } else {
      (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_ctrl
    }
    
    df$sh_chilcotin_SSL_alpha[i] <-
      sh_chilcotin_intercept +
      df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
      df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef +
      df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    df$sh_chilcotin_alpha_CN[i] <-
      sh_chilcotin_intercept +
      df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
      df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df$sh_chilcotin_SL[i]   * sh_chilcotin_ssl_coef +
      df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef
    
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      # FIX: same NA-residual fallback as Thompson/chum above.
      # FIX 7: spreadsheet EU has no residual term for steelhead
      resid_adj <- if (!APPLY_SH_RESIDUALS || is.na(df$sh_chilcotin_ln_obs_pred[i])) 1 else exp(df$sh_chilcotin_ln_obs_pred[i])
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_chilcotin_alpha_CN[i] +
                      sh_chilcotin_spawners_coef * spk) *
            resid_adj
        } else {
          spk * exp(df$sh_chilcotin_SSL_alpha[i] +
                      sh_chilcotin_spawners_coef * spk) *
            resid_adj
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df$sh_chilcotin_recruits_alt[j]
      }
      
      df$sh_chilcotin_Nage4_pred[i] <- if (yr < start_year + 4) df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p4[i] else lag_recruits(yr - 4) * df$sh_chilcotin_p4[i] * 1000
      df$sh_chilcotin_Nage5_pred[i] <- if (yr < start_year + 5) df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p5[i] else lag_recruits(yr - 5) * df$sh_chilcotin_p5[i] * 1000
      df$sh_chilcotin_Nage6_pred[i] <- if (yr < start_year + 6) df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p6[i] else lag_recruits(yr - 6) * df$sh_chilcotin_p6[i] * 1000
      df$sh_chilcotin_Nage7_pred[i] <- if (yr < start_year + 7) df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p7[i] else lag_recruits(yr - 7) * df$sh_chilcotin_p7[i] * 1000
      df$sh_chilcotin_Nage8_pred[i] <- if (yr < start_year + 8) df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p8[i] else lag_recruits(yr - 8) * df$sh_chilcotin_p8[i] * 1000
      
      df$sh_chilcotin_sum_pred[i] <- sum(df$sh_chilcotin_Nage4_pred[i], df$sh_chilcotin_Nage5_pred[i],
                                         df$sh_chilcotin_Nage6_pred[i], df$sh_chilcotin_Nage7_pred[i],
                                         df$sh_chilcotin_Nage8_pred[i], na.rm = TRUE)
      
      df$sh_chilcotin_bycatch_pred[i] <- df$sh_chilcotin_sum_pred[i] * df$sh_chilcotin_U_comm[i]
      
      if (yr <= 2018) {
        df$sh_chilcotin_FN_catch_pred[i] <- df$sh_chilcotin_FN_mortalities[i]
      } else {
        denom_2018 <- df$sh_chilcotin_sum_pred[df$Year == 2018] - df$sh_chilcotin_bycatch_pred[df$Year == 2018]
        df$sh_chilcotin_FN_catch_pred[i] <- FN_chilcotin_2018 / denom_2018 *
          (df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_bycatch_pred[i])
      }
      
      df$sh_chilcotin_total_catch_pred[i] <- df$sh_chilcotin_FN_catch_pred[i] +
        df$sh_chilcotin_sport_mortalities[i] + df$sh_chilcotin_bycatch_pred[i]
      
      S_new <- max(df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_total_catch_pred[i], 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= 1e-8 * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    
    df$sh_chilcotin_spawners_pred[i] <- S_old
  }
  
  # FIX: same rationale as Thompson above.
  if (SSL_control == 0 && U_historic == 1 && byrate == default_byrate) {
    df <- df %>%
      mutate(
        sh_chilcotin_recruits_alt     = sh_chilcotin_recruits,
        sh_chilcotin_spawners_pred    = sh_chilcotin_spawners * 1000,
        sh_chilcotin_sum_pred         = sh_chilcotin_prefishery_N,
        sh_chilcotin_bycatch_pred     = sh_chilcotin_pred_bycatch,
        sh_chilcotin_FN_catch_pred    = sh_chilcotin_FN_mortalities,
        sh_chilcotin_total_catch_pred = sh_chilcotin_FN_mortalities +
          sh_chilcotin_sport_mortalities + sh_chilcotin_pred_bycatch,
        sh_chilcotin_U_comm           = sh_chilcotin_U
      )
  }
  
  df
}

## --- Combined pipeline: chum -> Thompson -> Chilcotin ----------------------

run_full_scenario <- function(chum_data, covariates, sh_data,
                              U_apply, SSL_control, U_historic = 0, byrate = default_byrate) {
  
  chum_df     <- run_chum_scenario(chum_data, covariates, U_apply, SSL_control, U_historic)
  all_data    <- chum_df %>% right_join(sh_data, by = "Year")
  thompson_df <- run_thompson_scenario(all_data, SSL_control, U_historic, byrate)
  full_df     <- run_chilcotin_scenario(thompson_df, SSL_control, U_historic, byrate)
  
  full_df %>%
    mutate(
      harvest_rate = U_apply,
      byrate       = byrate,
      U_historic   = U_historic,
      ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control")
    )
}

# ============================================================
# RUN SCENARIOS
# ============================================================

## 1. Historic harvest rate: SSL control vs. no control ----------------------
scenarios_historic_harvest <- purrr::map_dfr(c(0, 1), function(sc) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = 0, SSL_control = sc, U_historic = 1, byrate = default_byrate)
})

## 2. Harvest-rate sweep (byrate fixed): SSL control vs. no control ----------
scenarios_harvest_rate <- purrr::pmap_dfr(
  expand_grid(SSL_control = c(0, 1), U_apply = harvest_rates),
  function(SSL_control, U_apply) {
    run_full_scenario(data, covariates, sh_data,
                      U_apply = U_apply, SSL_control = SSL_control,
                      U_historic = 0, byrate = default_byrate)
  }
)

## 3. Bycatch-rate sweep (harvest held at historic): SSL control vs. no control
scenarios_bycatch_rate <- purrr::pmap_dfr(
  expand_grid(SSL_control = c(0, 1), byrate = bycatch_rates),
  function(SSL_control, byrate) {
    run_full_scenario(data, covariates, sh_data,
                      U_apply = 0, SSL_control = SSL_control,
                      U_historic = 1, byrate = byrate)
  }
)

# ============================================================
# PREDATION MORTALITY ACCOUNTING  (spreadsheet 'Retrospective Run1' FB:GL)
#
# "Run before predation"   = Ricker prediction with ZERO sea lions
#                            (spreadsheet FD / FU / GC, via SSLz / thompSSLz /
#                            chilcotSSLz). PRED_BASELINE <- "control" uses the
#                            SSL_CONTROL_YEAR (1978) level instead.
# "Recruits after predation" = same prediction at the scenario SSL.
# Killed = before - after, i.e. all predation (relative to zero sea lions).
# Spawners are the modeled scenario spawners (chum AP, Thompson CO,
# Chilcotin ES); predictions carry no obs residual, as in the spreadsheet.
#
# With PRED_BASELINE = "control", kills under "SSL control" are ~0 after
# SSL_CONTROL_YEAR by construction.
# ============================================================

# SSL standardization constants from the 'covariates' tab, used to
# back-transform z-scores to abundance (FK / FY / GI)
SSL_SCALE <- list(
  chum      = c(mean = 16088.708673270821, sd = 12075.157324801276),  # V3/V4
  thompson  = c(mean = 23610.490624783917, sd = 14522.985143567075),  # AC3/AC4
  chilcotin = c(mean = 21919.907643703482, sd = 14447.561481342822)   # AD3/AD4
)
ssl_backtransform <- function(z, stock) z * SSL_SCALE[[stock]][["sd"]] + SSL_SCALE[[stock]][["mean"]]

RUNDAYS <- 30   # 'rundays' (FQ8); was 60 in v9

# "control" = SSL_CONTROL_YEAR level (consistent with catch-lost);
# "zero"    = zero sea lions, as in spreadsheet FD / FU / GC.
PRED_BASELINE <- "zero"
baseline_z <- function(df, col, stock) {
  if (PRED_BASELINE == "zero") {
    unname(-SSL_SCALE[[stock]][["mean"]] / SSL_SCALE[[stock]][["sd"]])  # covariates!V83:V85
  } else {
    get_control_ssl(df, col, stock)
  }
}

# 6-yr trailing mean, expanding over the first 5 rows (FP column)
running_mean6 <- function(x) sapply(seq_along(x), function(i) mean(x[max(1, i - 5):i], na.rm = TRUE))

add_predation_block <- function(stock_df, stock_label, stock_key, run_before, run_after, ssl_z) {
  stock_df %>%
    mutate(
      Stock                  = stock_label,
      run_before_pred        = run_before,                                   # FD / FU / GC
      recruits_after         = run_after,                                    # FF / FV / GE
      killed_by_pred         = run_before_pred - recruits_after,             # FG / FW / GF
      frac_killed            = killed_by_pred / run_before_pred,             # FI / FX / GH
      SSL_abundance          = ssl_backtransform(ssl_z, stock_key),          # FK / FY / GI
      killed_per_SSL         = killed_by_pred / SSL_abundance,               # FM / FZ / GJ
      killed_per_SSL_day     = killed_per_SSL / RUNDAYS,                     # FO / GA / GK
      run_avg_killed_per_SSL = running_mean6(killed_per_SSL),                # FP
      delta_SSL              = (SSL_abundance - lag(SSL_abundance)) / lag(SSL_abundance),  # FQ
      SSL_rel_min            = SSL_abundance / min(SSL_abundance, na.rm = TRUE)            # FR / GL
    )
}

compute_predation_mortality <- function(scenario_df,
                                        chum_intercept = 1.03737862843252,
                                        pdo_adult_coef = 0.0929480015696915,
                                        npgo_coef      = 0.102617626088713,
                                        pdo_smolt_coef = -0.105950753725792,
                                        chum_ssl_coef  = -0.224246284912103,
                                        chum_S_coef    = -4.95136622626478E-07,
                                        th_intercept = 1.572107637, th_sst = -0.203463091,
                                        th_ssl = -0.764277677, th_npgo = -0.0402,
                                        th_S = -0.804438793,
                                        ch_intercept = 1.053608979, ch_sst = -0.127949278,
                                        ch_ssl = -0.792741195, ch_npgo = 0.152526045,
                                        ch_pdo = 0.202708011, ch_S = -1.022467631,
                                        th_include_npgo    = FALSE,
                                        th_ssl_denominator = c("actual", "scenario")) {
  th_ssl_denominator <- match.arg(th_ssl_denominator)
  
  ## --- Chum (FC:FR) ---------------------------------------------------------
  chum <- scenario_df %>% filter(Year %in% CHUM_YEARS) %>% distinct(Year, .keep_all = TRUE) %>% arrange(Year)
  z0   <- baseline_z(covariates, "SSL", "chum")   # scenario_df may not include SSL_CONTROL_YEAR for chum (right_join on sh_data)
  S    <- chum$chum_spawners_pred
  env  <- chum_intercept + pdo_adult_coef * chum$PDO_adult + npgo_coef * chum$NPGO +
    pdo_smolt_coef * chum$PDO_smolt
  chum <- add_predation_block(
    chum, "Chum", "chum",
    run_before = S * exp(env + chum_S_coef * S + chum_ssl_coef * z0),
    run_after  = S * exp(env + chum_S_coef * S + chum_ssl_coef * chum$chum_SSL_alt),
    ssl_z      = chum$chum_SSL_alt)
  
  ## --- Thompson steelhead (FT:GA) -------------------------------------------
  th  <- scenario_df %>% filter(Year %in% THOMPSON_YEARS) %>% distinct(Year, .keep_all = TRUE) %>% arrange(Year)
  z0  <- baseline_z(scenario_df, "sh_thompson_SL", "thompson")
  S   <- th$sh_thompson_spawners_pred                       # fish, not thousands
  # Spreadsheet FU/FV omit NPGO although the Thompson recruitment model
  # (BO/CV) includes it; default matches the spreadsheet.
  env <- th_intercept + th_sst * th$sh_thompson_SST +
    (if (th_include_npgo) th_npgo * th$sh_thompson_NPGO else 0)
  # Spreadsheet FY divides by ACTUAL SSL (BL) while FV uses scenario SSL (CS);
  # default matches the spreadsheet, "scenario" matches chum/Chilcotin.
  th_denom_z <- if (th_ssl_denominator == "actual") th$sh_thompson_SL else th$sh_thompson_SSL_alt
  th <- add_predation_block(
    th, "Thompson steelhead", "thompson",
    run_before = S * exp(env + th_S * S / 1000 + th_ssl * z0),
    run_after  = S * exp(env + th_S * S / 1000 + th_ssl * th$sh_thompson_SSL_alt),
    ssl_z      = th_denom_z)
  
  ## --- Chilcotin steelhead (GB:GK) ------------------------------------------
  ch  <- scenario_df %>% filter(Year %in% CHILCOTIN_YEARS) %>% distinct(Year, .keep_all = TRUE) %>% arrange(Year)
  z0  <- baseline_z(scenario_df, "sh_chilcotin_SL", "chilcotin")
  S   <- ch$sh_chilcotin_spawners_pred
  env <- ch_intercept + ch_sst * ch$sh_chilcotin_SST + ch_npgo * ch$sh_chilcotin_NPGO +
    ch_pdo * ch$sh_chilcotin_PDO
  ch <- add_predation_block(
    ch, "Chilcotin steelhead", "chilcotin",
    run_before = S * exp(env + ch_S * S / 1000 + ch_ssl * z0),
    run_after  = S * exp(env + ch_S * S / 1000 + ch_ssl * ch$sh_chilcotin_SSL_alt),
    ssl_z      = ch$sh_chilcotin_SSL_alt)
  
  keep <- c("Year", "ssl_scenario", "harvest_rate", "byrate", "Stock",
            "run_before_pred", "recruits_after", "killed_by_pred", "frac_killed",
            "SSL_abundance", "killed_per_SSL", "killed_per_SSL_day",
            "run_avg_killed_per_SSL", "delta_SSL", "SSL_rel_min")
  bind_rows(select(chum, all_of(keep)), select(th, all_of(keep)), select(ch, all_of(keep)))
}

## Run per SSL scenario so lags / running means don't mix scenarios
predation_historic <- scenarios_historic_harvest %>%
  group_split(ssl_scenario) %>%
  map_dfr(compute_predation_mortality)

## Thompson steelhead: total run before predation (FU) vs. total recruits
## after predation (FV) under historical pinniped conditions -- the
## "No SSL control" / historic-harvest run, i.e. observed spawners and
## actual SSL abundance. The shaded gap is fish killed by sea lions (FW).
th_run_before_after <- predation_historic %>%
  filter(Stock == "Thompson steelhead", ssl_scenario == "No SSL control") %>%
  filter_stock_years()

ggplot(th_run_before_after, aes(Year)) +
  geom_ribbon(aes(ymin = recruits_after, ymax = run_before_pred),
              fill = "#eb6834", alpha = 0.12) +
  geom_line(aes(y = run_before_pred, color = "Total run before predation"), linewidth = 0.9) +
  geom_line(aes(y = recruits_after,  color = "Total recruits after predation"), linewidth = 0.9) +
  scale_color_manual(values = c("Total run before predation"     = "#2a78d6",
                                "Total recruits after predation" = "#eb6834")) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  labs(x = "Brood year", y = "Thompson steelhead (fish)", color = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("figures/thompson_run_before_after_historic.png", width = 10, height = 5, dpi = 600)

## Same quantities at the spreadsheet's saved settings (SSLcontrol = 1,
## Uhistoric = 0, Uapply = 0.2, byrate = 0.69) -- reproduces the Excel chart
## of FT (spawners), FU (before), FV (after) and FW (killed). Here SSLs are
## held at the 1978 level, so the fraction killed is constant (~42%).
predation_u02 <- scenarios_harvest_rate %>%
  filter(abs(harvest_rate - 0.2) < 1e-9, ssl_scenario == "SSL control") %>%
  compute_predation_mortality()

th_spawners_u02 <- scenarios_harvest_rate %>%
  filter(abs(harvest_rate - 0.2) < 1e-9, ssl_scenario == "SSL control") %>%
  distinct(Year, .keep_all = TRUE) %>%
  select(Year, spawners = sh_thompson_spawners_pred)

th_sheet_equiv <- predation_u02 %>%
  filter(Stock == "Thompson steelhead") %>%
  filter_stock_years() %>%
  left_join(th_spawners_u02, by = "Year") %>%
  transmute(Year,
            `Total TH run before predation`     = run_before_pred,   # FU
            `Total TH recruits after predation` = recruits_after,    # FV
            `TH brood year spawners`            = spawners,          # FT
            `Total TH run killed by predators`  = killed_by_pred) %>% # FW
  pivot_longer(-Year, names_to = "series", values_to = "fish") %>%
  mutate(series = factor(series, levels = c("Total TH run before predation",
                                            "Total TH recruits after predation",
                                            "TH brood year spawners",
                                            "Total TH run killed by predators")))

ggplot(th_sheet_equiv, aes(Year, fish, color = series, linetype = series)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#6b6a63")) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed")) +
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
  labs(x = "Brood year", y = "Total Thompson steelhead", color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave("figures/thompson_run_before_after_SSLcontrol_U02.png", width = 10, height = 5, dpi = 600)

## Long-term summaries (spreadsheet AG11 / FP10 = avg chums killed per SSL)
predation_historic %>%
  filter(ssl_scenario == "No SSL control", Year > SSL_CONTROL_YEAR) %>%
  group_by(Stock) %>%
  summarise(mean_killed_per_SSL = mean(killed_per_SSL, na.rm = TRUE),
            mean_frac_killed    = mean(frac_killed, na.rm = TRUE),
            .groups = "drop") %>%
  print()

# ============================================================
# METRICS EXTRACTION: tidy productivity + returns for all 3 stocks
# ============================================================

extract_stock_metrics <- function(scenario_df) {
  bind_rows(
    scenario_df %>%
      transmute(Year, ssl_scenario,
                harvest_rate = harvest_rate, byrate = byrate,
                Stock = "Chum",
                productivity = log(chum_recruits_alt / chum_spawners_pred),
                returns = chum_recruits_alt),
    scenario_df %>%
      transmute(Year, ssl_scenario,
                harvest_rate = harvest_rate, byrate = byrate,
                Stock = "Thompson steelhead",
                productivity = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
                returns = sh_thompson_recruits_alt),
    scenario_df %>%
      transmute(Year, ssl_scenario,
                harvest_rate = harvest_rate, byrate = byrate,
                Stock = "Chilcotin steelhead",
                productivity = log(sh_chilcotin_recruits_alt / (sh_chilcotin_spawners_pred / 1000)),
                returns = sh_chilcotin_recruits_alt)
  )
}

metrics_historic_harvest <- extract_stock_metrics(scenarios_historic_harvest)
metrics_harvest_rate     <- extract_stock_metrics(scenarios_harvest_rate)
metrics_bycatch_rate     <- extract_stock_metrics(scenarios_bycatch_rate)

# ============================================================
# PLOTS
# ============================================================

plot_scenario_comparison <- function(long_df, metric_col, y_label, title) {
  df <- filter_stock_years(long_df)
  ggplot(df, aes(Year, .data[[metric_col]], color = ssl_scenario)) +
    geom_line(linewidth = 1, alpha = 0.85) +
    facet_wrap(~ Stock, scales = "free") +
    scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
    labs(x = "Year", y = y_label, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_scenario_sweep <- function(long_df, metric_col, group_var, group_label, y_label, title) {
  df <- filter_stock_years(long_df)
  ggplot(df, aes(Year, .data[[metric_col]], color = factor(.data[[group_var]]), linetype = ssl_scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    facet_wrap(~ Stock, scales = "free") +
    scale_color_viridis_d(name = group_label) +
    labs(x = "Year", y = y_label, linetype = "SSL scenario") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

## 1. Historic harvest: productivity & returns --------------------------------
(p1 <- plot_scenario_comparison(metrics_historic_harvest, "productivity", "Productivity (alpha)", NULL))


ggsave("figures/chum_sh_productivity_historic.png", width = 11, height = 5, dpi = 600)

p2 <- plot_scenario_comparison(metrics_historic_harvest, "returns", "Returns", NULL)

(p1 + theme(legend.position = "none") + xlab("")) / p2

ggsave("figures/chum_sh_returns_productivity_historic.png", width = 10, height = 5, dpi = 600)

## 2. Harvest-rate sweep: productivity & returns ------------------------------
#plot_scenario_sweep(metrics_harvest_rate, "productivity", "harvest_rate", "Harvest rate",
#                    "Productivity (alpha)", "Productivity across chum harvest-rate scenarios")
#ggsave("figures/chum_sh_productivity_harvest_sweep.png", width = 11, height = 5, dpi = 600)

plot_scenario_sweep(metrics_harvest_rate, "returns", "harvest_rate", "Harvest rate",
                    "Returns")
ggsave("figures/chum_sh_returns_harvest_sweep.png", width = 11, height = 5, dpi = 600)

## 3. Bycatch-rate sweep: productivity & returns ------------------------------
plot_scenario_sweep(metrics_bycatch_rate, "productivity", "byrate", "Bycatch rate",
                    "")
ggsave("figures/chum_sh_productivity_bycatch_sweep.png", width = 11, height = 5, dpi = 600)

plot_scenario_sweep(metrics_bycatch_rate, "returns", "byrate", "Bycatch rate",
                    "Returns", "Returns across steelhead bycatch-rate scenarios")
ggsave("figures/chum_sh_returns_bycatch_sweep.png", width = 11, height = 5, dpi = 600)

# ============================================================
# CHUM CATCH LOST TO PINNIPEDS (historic harvest scenario)
# "SSL control" holds SSLs at the SSL_CONTROL_YEAR (1978) level, so
# catch_lost is measured against 1978 sea lions.
# ============================================================

chum_catch_lost_by_year <- scenarios_historic_harvest %>%
  filter(Year >= 1978, Year < 2017) %>%   # reliable model window
  select(Year, ssl_scenario, catch_alt) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  mutate(catch_lost = `No SSL control` - `SSL control`)

## Average yearly catch lost, 2000-present ------------------------------------

RECENT_YEARS_CHUM <- 2000:max(chum_catch_lost_by_year$Year, na.rm = TRUE)
chum_catch_lost_recent <- chum_catch_lost_by_year %>% filter(Year %in% RECENT_YEARS_CHUM)

chum_mean_catch_lost <- chum_catch_lost_recent %>%
  summarise(
    mean_catch_lost = mean(catch_lost, na.rm = TRUE),
    se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n())
  )
print(chum_mean_catch_lost)

# Bar chart
chum_lost <- ggplot(chum_mean_catch_lost %>% mutate(label = "Chum"), aes(label, -mean_catch_lost)) +
  geom_col(width = 0.4, fill = "#FF4500") +
  geom_errorbar(aes(ymin = -mean_catch_lost - se_catch_lost,
                    ymax = -mean_catch_lost + se_catch_lost), width = 0.1) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Mean yearly catch lost") +
  theme_minimal()
ggsave("figures/chum_mean_catch_lost_2000-present.png", width = 5, height = 5.5, dpi = 600)

# Boxplot (yearly distribution)
ggplot(chum_catch_lost_recent, aes(x = "Chum", y = catch_lost)) +
  geom_boxplot(width = 0.3, fill = "#FF4500", alpha = 0.7, outlier.shape = 21) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Yearly catch lost",
       title = paste0("Chum: distribution of yearly catch lost to pinnipeds, ",
                      min(RECENT_YEARS_CHUM), "-", max(RECENT_YEARS_CHUM))) +
  theme_minimal()
ggsave("figures/chum_catch_lost_boxplot_2000-present.png", width = 5, height = 5.5, dpi = 600)

# ============================================================
# CHUM: mean yearly catch lost to pinnipeds, by harvest rate
# Uses scenarios_harvest_rate (the U_apply sweep, byrate fixed) instead of
# the single historic-harvest scenario -- same "No SSL control minus SSL
# control" logic, just repeated for each harvest rate in the sweep.
# ============================================================

chum_catch_lost_by_harvest_rate <- scenarios_harvest_rate %>%
  filter(Year >= 1978, Year < 2017) %>%   # reliable model window
  select(Year, harvest_rate, ssl_scenario, catch_alt) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  mutate(catch_lost = `No SSL control` - `SSL control`)

chum_catch_lost_harvest_rate_recent <- chum_catch_lost_by_harvest_rate %>%
  filter(Year %in% RECENT_YEARS_CHUM)   # 2000-present, defined earlier

chum_mean_catch_lost_by_harvest_rate <- chum_catch_lost_harvest_rate_recent %>%
  group_by(harvest_rate) %>%
  summarise(
    mean_catch_lost = mean(catch_lost, na.rm = TRUE),
    se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(chum_mean_catch_lost_by_harvest_rate,
       aes(factor(harvest_rate), -mean_catch_lost, fill = factor(harvest_rate))) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = -mean_catch_lost - se_catch_lost,
                    ymax = -mean_catch_lost + se_catch_lost), width = 0.2) +
  scale_fill_viridis_d(guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Harvest rate", y = "Mean yearly catch lost") +
  theme_minimal()

ggsave("figures/chum_mean_catch_lost_by_harvest_rate_2000-present.png", width = 9, height = 5.5, dpi = 600)

# ============================================================
# CHUM: mean yearly catch lost, historic harvest + harvest-rate sweep together
# Combines chum_mean_catch_lost (historic scenario) and
# chum_mean_catch_lost_by_harvest_rate (the U_apply sweep) into one bar
# chart, with "Historic" placed first.
# ============================================================

chum_catch_lost_all_bars <- bind_rows(
  chum_mean_catch_lost %>%
    mutate(scenario_label = "Historic"),
  chum_mean_catch_lost_by_harvest_rate %>%
    mutate(scenario_label = as.character(harvest_rate))
) %>%
  mutate(scenario_label = factor(scenario_label,
                                 levels = c("Historic", as.character(sort(unique(harvest_rates))))))

ggplot(chum_catch_lost_all_bars,
       aes(scenario_label, -mean_catch_lost, fill = scenario_label == "Historic")) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = -mean_catch_lost - se_catch_lost,
                    ymax = -mean_catch_lost + se_catch_lost), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#FF4500", "FALSE" = "grey50"), guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Harvest rate scenario", y = "Mean yearly catch lost") +
  theme_minimal()

ggsave("figures/chum_mean_catch_lost_historic_and_harvest_sweep_2000-present.png",
       width = 10, height = 5.5, dpi = 600)

# ============================================================
# Productivity & returns across harvest rates 0.1-0.8
# Same style as the historic-vs-pinniped comparison, but faceted into a
# 2-column (productivity / returns) x 8-row (harvest rate) grid, one plot
# per stock. Uses metrics_harvest_rate from the main scenario script.
# ============================================================

plot_data_sweep <- metrics_harvest_rate %>%
  filter(harvest_rate %in% seq(0.1, 0.8, by = 0.1)) %>%
  pivot_longer(cols = c(productivity, returns), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, productivity = "ln(R/S)", returns = "Returns"))

plot_harvest_rate_grid <- function(long_df, stock_name) {
  df <- long_df %>%
    filter(Stock == stock_name,
           Year >= stock_year_bounds$yr_min[stock_year_bounds$Stock == stock_name],
           Year <= stock_year_bounds$yr_max[stock_year_bounds$Stock == stock_name])
  
  ggplot(df, aes(Year, value, color = ssl_scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    facet_grid(harvest_rate ~ metric, scales = "free_y") +
    scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
    labs(x = "Year", y = NULL, color = NULL,
         title = paste0(stock_name, ": productivity & returns")) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
}

stocks <- c("Chum", "Thompson steelhead", "Chilcotin steelhead")

for (s in stocks) {
  p <- plot_harvest_rate_grid(plot_data_sweep, s)
  print(p)
  ggsave(
    filename = paste0("figures/", tolower(gsub(" ", "_", s)), "_productivity_returns_harvest_grid.png"),
    plot = p, width = 7, height = 14, dpi = 600
  )
}

# ============================================================
# CHUM: recruits per spawner & returns across harvest rates 0-0.8,
# with historical observed overlaid
#
# The previous facet_grid(harvest_rate ~ metric, scales = "free_y")
# only frees the y-axis PER ROW, and both columns in a row still share
# that one scale. Since Returns is in the millions and Recruits per
# spawner is ~1-5, the shared per-row scale is dominated by Returns and
# Recruits-per-spawner flatlines at ~0.
#
# Fix: build each metric as its own facet_wrap(~harvest_rate, ncol = 1,
# scales = "free_y") column -- facet_wrap frees scales fully per panel,
# not per row/column -- then place the two columns side by side with
# patchwork.
# ============================================================

## Model output: recruits-per-spawner + returns, by harvest rate & SSL scenario

chum_harvest_metrics <- scenarios_harvest_rate %>%
  filter(harvest_rate %in% seq(0.1, 0.8, by = 0.1), Year %in% CHUM_YEARS) %>%
  transmute(Year, harvest_rate, ssl_scenario,
            `Recruits per spawner` = chum_recruits_alt / chum_spawners_pred,
            Returns = chum_recruits_alt) %>%
  pivot_longer(cols = c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value")

## Historical observed: doesn't vary by scenario/harvest_rate -- pull once
## (from any single harvest-rate/SSL-scenario combo, since chum_recruits_obs
## and chum_spawners are the same input data regardless of scenario) and
## repeat across every harvest_rate row

chum_observed <- scenarios_harvest_rate %>%
  filter(ssl_scenario == "No SSL control", harvest_rate == harvest_rates[1], Year %in% CHUM_YEARS) %>%
  transmute(Year,
            `Recruits per spawner` = chum_recruits_obs / chum_spawners,
            Returns = chum_recruits_obs) %>%
  pivot_longer(cols = c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value") %>%
  mutate(ssl_scenario = "Observed")

chum_observed_repeated <- tidyr::crossing(
  harvest_rate = unique(chum_harvest_metrics$harvest_rate),
  chum_observed
)

chum_sweep_with_obs <- bind_rows(chum_harvest_metrics, chum_observed_repeated)

color_values    <- c("Observed" = "black", "No SSL control" = "#F8766D", "SSL control" = "#00BFC4")
linetype_values <- c("Observed" = "solid", "No SSL control" = "solid", "SSL control" = "dashed")

plot_chum_metric_stack <- function(df, metric_name, y_label) {
  df %>%
    filter(metric == metric_name) %>%
    ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    facet_wrap(~ harvest_rate, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = color_values) +
    scale_linetype_manual(values = linetype_values) +
    labs(x = "Year", y = y_label, color = NULL, linetype = NULL, title = metric_name) +
    theme_minimal(base_size = 9) +
    theme(strip.placement = "outside")
}

p_recruits_per_spawner <- plot_chum_metric_stack(chum_sweep_with_obs, "Recruits per spawner", "Recruits per spawner")
p_returns              <- plot_chum_metric_stack(chum_sweep_with_obs, "Returns", "Returns")

(p_recruits_per_spawner + theme(legend.position = "none")) | p_returns +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/chum_harvest_rate_pred_control_with_observed.png", width = 9, height = 15, dpi = 600)

# check numbers 
scenarios_historic_harvest %>%
  filter(ssl_scenario == "No SSL control", Year >= 1978, Year < 2017) %>%
  summarise(min_U = min(U_chum, na.rm = TRUE),
            max_U = max(U_chum, na.rm = TRUE),
            mean_U = mean(U_chum, na.rm = TRUE))

# ============================================================
# STEELHEAD (Thompson + Chilcotin): productivity & returns across
# bycatch rates 0-0.8, SSL control vs. none, with historical observed
# overlaid. Same patchwork-column approach as the chum grid, to avoid
# the facet_grid shared-scale issue.
#
#
# NOTE on "Returns": recruits_alt (used in the earlier harvest-rate grid)
# is on a spawners/1000-normalized scale, not directly comparable to the
# observed historical run size. Here "Returns" uses sh_*_sum_pred (the
# modeled total return, raw-count scale) instead, so it lines up with
# the observed sh_*_prefishery_N series.
#
# NOTE on "Productivity": sh_*_alpha_CN (the "No SSL control" alpha) is
# calculated directly from the actual historical SL/SST/NPGO covariates
# -- it IS the observed productivity trajectory. So "Observed" and
# "No SSL control" overlap exactly on the productivity panels; they're
# both plotted for a consistent 3-line legend across panels, but expect
# them to sit on top of each other.
# ============================================================

## Model output: productivity + returns, by bycatch rate & SSL scenario ------

sh_bycatch_metrics <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(byrate <= 0.8, Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead",
              productivity = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
              returns = sh_thompson_sum_pred),
  scenarios_bycatch_rate %>%
    filter(byrate <= 0.8, Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Chilcotin steelhead",
              productivity = log(sh_chilcotin_recruits_alt / (sh_chilcotin_spawners_pred / 1000)),
              returns = sh_chilcotin_sum_pred)
) %>%
  pivot_longer(cols = c(productivity, returns), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, productivity = "ln(R/S)", returns = "Returns"))

## Historical observed: doesn't vary by scenario/byrate -- pull once and
## repeat across every byrate row

sh_observed <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% THOMPSON_YEARS) %>%
    transmute(Year, Stock = "Thompson steelhead",
              productivity = sh_thompson_alpha_CN,
              returns = sh_thompson_prefishery_N),
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, Stock = "Chilcotin steelhead",
              productivity = sh_chilcotin_alpha_CN,
              returns = sh_chilcotin_prefishery_N)
) %>%
  pivot_longer(cols = c(productivity, returns), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, productivity = "Productivity (alpha)", returns = "Returns"),
         ssl_scenario = "Observed")

sh_observed_repeated <- tidyr::crossing(byrate = unique(sh_bycatch_metrics$byrate), sh_observed)

sh_bycatch_with_obs <- bind_rows(sh_bycatch_metrics, sh_observed_repeated)

## Plot: 2 metric columns x 9 byrate rows, per stock ---------------------------

color_values_sh    <- c("Observed" = "black", "No SSL control" = "#F8766D", "SSL control" = "#00BFC4")
linetype_values_sh <- c("Observed" = "solid", "No SSL control" = "solid", "SSL control" = "dashed")

plot_sh_metric_stack <- function(df, stock_name, metric_name, y_label) {
  df %>%
    filter(Stock == stock_name, metric == metric_name) %>%
    ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    facet_wrap(~ byrate, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = color_values_sh) +
    scale_linetype_manual(values = linetype_values_sh) +
    labs(x = "Year", y = y_label, color = NULL, linetype = NULL, title = metric_name) +
    theme_minimal(base_size = 9) +
    theme(strip.placement = "outside")
}

for (s in c("Thompson steelhead", "Chilcotin steelhead")) {
  
  p_prod <- plot_sh_metric_stack(sh_bycatch_with_obs, s, "ln(R/S)", "ln(R/S)")
  p_ret  <- plot_sh_metric_stack(sh_bycatch_with_obs, s, "Returns", "Returns")
  
  p <- (p_prod + theme(legend.position = "none")) | p_ret +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom") &
    plot_annotation(title = paste0(s, ": productivity & returns across bycatch rates (0-0.8)"))
  
  print(p)
  ggsave(
    filename = paste0("figures/", tolower(gsub(" ", "_", s)), "_productivity_returns_bycatch_grid_with_observed.png"),
    plot = p, width = 9, height = 12, dpi = 600
  )
}

### byrate exploration 
## Isolate the byrate effect: % difference in returns vs. the byrate = 0
## baseline, "No SSL control" only (drops the SSL-control spike that
## dominates the shared y-axis and swamps the byrate signal).

byrate_effect <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, Stock = "Thompson steelhead", returns = sh_thompson_sum_pred),
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, Stock = "Chilcotin steelhead", returns = sh_chilcotin_sum_pred)
)

baseline <- byrate_effect %>%
  filter(byrate == 0) %>%
  select(Year, Stock, baseline_returns = returns)

byrate_pct_diff <- byrate_effect %>%
  filter(byrate > 0) %>%
  left_join(baseline, by = c("Year", "Stock")) %>%
  mutate(pct_diff = 100 * (returns - baseline_returns) / baseline_returns)

ggplot(byrate_pct_diff, aes(Year, pct_diff, color = factor(byrate))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  facet_wrap(~ Stock, scales = "free", ncol = 1) +
  scale_color_viridis_d(name = "Bycatch rate") +
  labs(x = "Year", y = "% difference in returns vs. byrate = 0",
       title = "No SSL control: sensitivity of returns to bycatch rate") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/byrate_returns_pct_diff.png", width = 9, height = 8, dpi = 300)

## Compare bycatch-rate sensitivity between SSL scenarios: same % diff vs.
## byrate = 0 idea, computed separately within each ssl_scenario, then
## shown side by side (columns = ssl_scenario, rows = Stock).

byrate_effect_both <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead", returns = sh_thompson_sum_pred),
  scenarios_bycatch_rate %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Chilcotin steelhead", returns = sh_chilcotin_sum_pred)
)

baseline_both <- byrate_effect_both %>%
  filter(byrate == 0) %>%
  select(Year, Stock, ssl_scenario, baseline_returns = returns)

byrate_pct_diff_both <- byrate_effect_both %>%
  filter(byrate > 0) %>%
  left_join(baseline_both, by = c("Year", "Stock", "ssl_scenario")) %>%
  mutate(pct_diff = 100 * (returns - baseline_returns) / baseline_returns)

ggplot(byrate_pct_diff_both, aes(Year, pct_diff, color = factor(byrate))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  facet_grid(Stock ~ ssl_scenario, scales = "free_y") +
  scale_color_viridis_d(name = "Bycatch rate") +
  labs(x = "Year", y = "% difference in returns vs. byrate = 0",
       title = "Sensitivity of returns to bycatch rate: SSL control vs. none") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/byrate_returns_pct_diff_both_scenarios.png", width = 11, height = 8, dpi = 300)

## Baseline run at byrate ~ 0.69 (default_byrate), via the SAME fully-
## simulated pathway as every other point in the sweep. A tiny epsilon
## offset avoids the `byrate == default_byrate` exact-match check that
## would otherwise trigger the observed-data passthrough for "No SSL
## control" -- the offset is numerically negligible to the model itself
## (byrate only enters as a multiplier on U_comm), but keeps this baseline
## on the same simulated footing as the rest of the sweep, so the
## comparison isolates *only* the byrate effect, not model-vs-observed drift.
scenarios_byrate_baseline <- purrr::map_dfr(c(0, 1), function(sc) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = 0, SSL_control = sc, U_historic = 1,
                    byrate = default_byrate + 1e-9)
})

byrate_effect_both <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead", returns = sh_thompson_sum_pred),
  scenarios_bycatch_rate %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Chilcotin steelhead", returns = sh_chilcotin_sum_pred)
)

baseline_sim <- bind_rows(
  scenarios_byrate_baseline %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, ssl_scenario, Stock = "Thompson steelhead", baseline_returns = sh_thompson_sum_pred),
  scenarios_byrate_baseline %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, ssl_scenario, Stock = "Chilcotin steelhead", baseline_returns = sh_chilcotin_sum_pred)
)

byrate_pct_diff_sim <- byrate_effect_both %>%
  left_join(baseline_sim, by = c("Year", "Stock", "ssl_scenario")) %>%
  mutate(pct_diff = 100 * (returns - baseline_returns) / baseline_returns)

ggplot(byrate_pct_diff_sim, aes(Year, pct_diff, color = factor(byrate))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  facet_grid(Stock ~ ssl_scenario, scales = "free_y") +
  scale_color_viridis_d(name = "Bycatch rate") +
  labs(x = "Year", y = "% difference in returns vs. byrate = 0.69 (simulated baseline)") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/byrate_returns_pct_diff_vs_069_simbaseline.png", width = 11, height = 8, dpi = 300)

## Modeled recruits (brood-year total, both SSL scenarios) across bycatch
## rates, plus observed recruits (doesn't vary by scenario/byrate) overlaid
## on every byrate row, faceted by stock.

recruits_modeled <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead",
              recruits = sh_thompson_recruits_alt),
  scenarios_bycatch_rate %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Chilcotin steelhead",
              recruits = sh_chilcotin_recruits_alt)
)

recruits_observed <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% THOMPSON_YEARS) %>%
    transmute(Year, Stock = "Thompson steelhead", recruits = sh_thompson_recruits),
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, Stock = "Chilcotin steelhead", recruits = sh_chilcotin_recruits)
) %>%
  mutate(ssl_scenario = "Observed")

## Repeat observed across every byrate row so it appears in every facet
recruits_observed_repeated <- tidyr::crossing(
  byrate = unique(recruits_modeled$byrate),
  recruits_observed
)

recruits_with_obs <- bind_rows(recruits_modeled, recruits_observed_repeated) %>%
  mutate(ssl_scenario = factor(ssl_scenario, levels = c("No SSL control", "Observed", "SSL control")))

color_values_sh <- c("No SSL control" = "#E41A1C", "Observed" = "darkgrey", "SSL control" = "#377EB8")
linetype_values_sh <- c("No SSL control" = "solid", "Observed" = "solid", "SSL control" = "dashed")

ggplot(recruits_with_obs, aes(Year, recruits, color = ssl_scenario, linetype = ssl_scenario)) +
  geom_line(linewidth = 1, alpha = 0.7) +
  facet_grid(byrate ~ Stock) +
  scale_color_manual(values = color_values_sh) +
  scale_linetype_manual(values = linetype_values_sh) +
  labs(x = "Year", y = "Recruits (thousands, brood-year total)", color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))

ggsave("figures/steelhead_recruits_byrate_grid.png", width = 9, height = 16, dpi = 300)

recruits_modeled <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead",
              recruits = sh_thompson_recruits_alt * 1000),
  scenarios_bycatch_rate %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Chilcotin steelhead",
              recruits = sh_chilcotin_recruits_alt * 1000)
)

recruits_observed <- bind_rows(
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% THOMPSON_YEARS) %>%
    transmute(Year, Stock = "Thompson steelhead", recruits = sh_thompson_recruits * 1000),
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, Stock = "Chilcotin steelhead", recruits = sh_chilcotin_recruits * 1000)
) %>%
  mutate(ssl_scenario = "Observed")

recruits_observed_repeated <- tidyr::crossing(
  byrate = unique(recruits_modeled$byrate),
  recruits_observed
)

recruits_with_obs <- bind_rows(recruits_modeled, recruits_observed_repeated) %>%
  mutate(ssl_scenario = factor(ssl_scenario, levels = c("No SSL control", "Observed", "SSL control")))

color_values_sh <- c("No SSL control" = "#E41A1C", "Observed" = "black", "SSL control" = "#377EB8")
linetype_values_sh <- c("No SSL control" = "solid", "Observed" = "dashed", "SSL control" = "solid")

ggplot(recruits_with_obs, aes(Year, recruits, color = ssl_scenario, linetype = ssl_scenario)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_grid(byrate ~ Stock, scales = "free_y") +
  scale_color_manual(values = color_values_sh) +
  scale_linetype_manual(values = linetype_values_sh) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Recruits", color = NULL, linetype = NULL) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))

ggsave("figures/steelhead_recruits_byrate_grid_realscale.png", width = 9, height = 16, dpi = 300)

recruits_with_obs %>%
  filter(Year > 1989) %>%
  filter(ssl_scenario != "SSL control") %>%
  ggplot(aes(Year, recruits, color = ssl_scenario, linetype = ssl_scenario)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_grid(byrate ~ Stock, scales = "free_y") +
  scale_color_manual(values = color_values_sh) +
  scale_linetype_manual(values = linetype_values_sh) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Recruits", color = NULL, linetype = NULL) +
  ylim(0,3000) + 
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))

# Steelhead: % increase in mean abundance over each stock's low-abundance
# period, under the pinniped-control ("SSL control") scenario -- styled to
# match the sockeye retrospective bar chart.
#
# Assumptions (I don't have your exact sockeye "low-abundance period"
# definition, so adjust as needed):
#  - Baseline = "No SSL control" (historic reconstruction) within each
#    stock's own year range (THOMPSON_YEARS / CHILCOTIN_YEARS).
#  - "Low-abundance period" for each stock = the years where that stock's
#    baseline abundance is <= its own 25th percentile within its own range
#    (i.e. bottom quartile of years, not a single low year).
#  - Abundance = sh_thompson_sum_pred / sh_chilcotin_sum_pred, same metric
#    used in the earlier line-chart version.
#  - pct_increase = (mean abundance under SSL control - mean abundance
#    under No SSL control) / mean abundance under No SSL control * 100,
#    both means taken over that stock's own low-abundance years only.
#
# Requires scenarios_historic_harvest, THOMPSON_YEARS, and CHILCOTIN_YEARS
# to already exist (run the main scenario script first).

steelhead_abundance <- bind_rows(
  scenarios_historic_harvest %>%
    filter(Year %in% THOMPSON_YEARS) %>%
    transmute(Year, ssl_scenario, Stock = "Thompson steelhead", abundance = sh_thompson_sum_pred),
  scenarios_historic_harvest %>%
    filter(Year %in% CHILCOTIN_YEARS) %>%
    transmute(Year, ssl_scenario, Stock = "Chilcotin steelhead", abundance = sh_chilcotin_sum_pred)
)

# Identify each stock's low-abundance years from the baseline ("No SSL
# control") series only
low_period_years <- steelhead_abundance %>%
  filter(ssl_scenario == "No SSL control") %>%
  group_by(Stock) %>%
  filter(abundance <= quantile(abundance, 0.25, na.rm = TRUE)) %>%
  distinct(Stock, Year)

pct_increase_over_low_period <- steelhead_abundance %>%
  inner_join(low_period_years, by = c("Stock", "Year")) %>%
  group_by(Stock, ssl_scenario) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = ssl_scenario, values_from = mean_abundance) %>%
  mutate(
    pct_increase = (`SSL control` - `No SSL control`) / `No SSL control` * 100,
    scenario     = "Pinniped scenario"
  )

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

ggsave("figures/steelhead_pct_increase_low_period_pinniped_scenario.png", width = 7, height = 5.5, dpi = 600)