# Chum / steelhead deterministic retrospective model
# Haley Oleynik, Murdoch McAllister
#
# Reproduces the 'Retrospective Run1' tab of
# Fraser_Chum_data_v19_alpha_SSL_est_fin_yrs_v12.xlsx and produces:
#   1. Run before vs. after sea lion predation (spawners, run before,
#      recruits after, killed) under SSL control at chum U = 0.2 -- one
#      figure per stock (chum, Thompson, Chilcotin)
#   2. The same for both steelhead stocks across steelhead bycatch rates
#   3. chum_sh_returns_productivity_historic
#   4. chum_mean_catch_lost_historic_and_harvest_sweep_2000-present
#   5. chum_harvest_rate_pred_control_with_observed
#
# Key modeling choices (all match the spreadsheet):
#   * SSL control holds sea lions at the 1978 covariate value for every
#     later brood year (SSL_CONTROL_YEAR).
#   * Steelhead recruits are pure Ricker predictions (no observed residual);
#     chum recruits keep the observed residual (spreadsheet AW: *EXP(N)).
#   * "No SSL control" with historic harvest reproduces the observed record.
#   * Run before predation = Ricker prediction at ZERO sea lions
#     (thompSSLz etc.); after predation = same spawners at scenario SSLs.

library(tidyverse)
library(patchwork)
library(scales)

# ============================================================
# CONFIG
# ============================================================

harvest_rates  <- seq(0, 0.8, by = 0.1)   # chum harvest-rate sweep
harvest_rates  <- round(seq(0, 0.8, by = 0.05), 2)   # chum harvest-rate sweep
bycatch_rates  <- seq(0, 1, by = 0.1)     # steelhead bycatch-rate sweep
default_byrate <- 0.69                    # steelhead bycatch/FN-mortality proxy rate
sheet_U_apply  <- 0.2                     # spreadsheet's saved chum harvest rate (Uapply)

SSL_CONTROL_YEAR <- 1978
# Optional override of the SSL z-score at SSL_CONTROL_YEAR (NA = read from data).
# Spreadsheet 1978 values: chum -0.68849, Thompson -0.90462, Chilcotin -0.87705.
SSL_CONTROL_Z <- c(chum = NA_real_, thompson = NA_real_, chilcotin = NA_real_)

APPLY_SH_RESIDUALS <- FALSE   # spreadsheet uses pure Ricker for steelhead

# Per-stock plotting windows (covariates end 2016)
CHUM_YEARS      <- 1951:2016
THOMPSON_YEARS  <- 1978:2016
CHILCOTIN_YEARS <- 1973:2016
stock_year_bounds <- tibble(
  Stock  = c("Chum", "Thompson steelhead", "Chilcotin steelhead"),
  yr_min = c(min(CHUM_YEARS), min(THOMPSON_YEARS), min(CHILCOTIN_YEARS)),
  yr_max = c(max(CHUM_YEARS), max(THOMPSON_YEARS), max(CHILCOTIN_YEARS))
)

filter_stock_years <- function(df) {
  df %>%
    left_join(stock_year_bounds, by = "Stock") %>%
    filter(Year >= yr_min, Year <= yr_max) %>%
    select(-yr_min, -yr_max)
}

get_control_ssl <- function(df, col, stock) {
  if (!is.na(SSL_CONTROL_Z[[stock]])) return(unname(SSL_CONTROL_Z[[stock]]))
  val <- df %>% filter(Year == SSL_CONTROL_YEAR) %>% pull(all_of(col)) %>% as.numeric()
  val <- val[!is.na(val)]
  if (length(val) == 0) {
    stop(sprintf("No %s value for %d in the input data; set SSL_CONTROL_Z[\"%s\"].",
                 col, SSL_CONTROL_YEAR, stock), call. = FALSE)
  }
  val[1]
}

dir.create("figures", showWarnings = FALSE)

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

## Historic harvest: SSL control vs. no control
scenarios_historic_harvest <- map_dfr(c(0, 1), function(sc) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = 0, SSL_control = sc, U_historic = 1, byrate = default_byrate)
})

## Chum harvest-rate sweep (steelhead byrate fixed at default)
scenarios_harvest_rate <- pmap_dfr(
  expand_grid(SSL_control = c(0, 1), U_apply = harvest_rates),
  function(SSL_control, U_apply) {
    run_full_scenario(data, covariates, sh_data,
                      U_apply = U_apply, SSL_control = SSL_control,
                      U_historic = 0, byrate = default_byrate)
  }
)

## Steelhead bycatch-rate sweep under SSL control at the spreadsheet's
## chum harvest rate (bycatch rate = byrate x chum U after 1990)
scenarios_bycatch_sslcontrol <- map_dfr(bycatch_rates, function(br) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = sheet_U_apply, SSL_control = 1,
                    U_historic = 0, byrate = br)
})

# ============================================================
# PREDATION MORTALITY ACCOUNTING  (spreadsheet cols FB:GL)
#   run_before_pred = Ricker prediction at zero sea lions   (FD / FU / GC)
#   recruits_after  = same spawners at scenario SSLs        (FF / FV / GE)
#   killed_by_pred  = before - after                        (FG / FW / GF)
# ============================================================

## SSL standardization (covariates tab): zero sea lions = -mean / sd
SSL_SCALE <- list(
  chum      = c(mean = 16088.708673270821, sd = 12075.157324801276),  # V3/V4
  thompson  = c(mean = 23610.490624783917, sd = 14522.985143567075),  # AC3/AC4
  chilcotin = c(mean = 21919.907643703482, sd = 14447.561481342822)   # AD3/AD4
)
ssl_zero_z <- function(stock) unname(-SSL_SCALE[[stock]][["mean"]] / SSL_SCALE[[stock]][["sd"]])

compute_predation <- function(scenario_df) {
  keep_cols <- c("Year", "ssl_scenario", "harvest_rate", "byrate")
  
  ## Chum (FC, FD, FF)
  chum <- scenario_df %>% filter(Year %in% CHUM_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- chum$chum_spawners_pred
  env <- 1.03737862843252 + 0.0929480015696915 * chum$PDO_adult +
    0.102617626088713 * chum$NPGO - 0.105950753725792 * chum$PDO_smolt +
    -4.95136622626478E-07 * S
  chum_out <- chum %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Chum", spawners = S,
           run_before_pred = S * exp(env - 0.224246284912103 * ssl_zero_z("chum")),
           recruits_after  = S * exp(env - 0.224246284912103 * chum$chum_SSL_alt))
  
  ## Thompson (FT, FU, FV) -- spreadsheet omits NPGO here
  th  <- scenario_df %>% filter(Year %in% THOMPSON_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- th$sh_thompson_spawners_pred
  env <- 1.572107637 - 0.804438793 * S / 1000 - 0.203463091 * th$sh_thompson_SST
  th_out <- th %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Thompson steelhead", spawners = S,
           run_before_pred = S * exp(env - 0.764277677 * ssl_zero_z("thompson")),
           recruits_after  = S * exp(env - 0.764277677 * th$sh_thompson_SSL_alt))
  
  ## Chilcotin (GB, GC, GE)
  ch  <- scenario_df %>% filter(Year %in% CHILCOTIN_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- ch$sh_chilcotin_spawners_pred
  env <- 1.053608979 - 1.022467631 * S / 1000 - 0.127949278 * ch$sh_chilcotin_SST +
    0.152526045 * ch$sh_chilcotin_NPGO + 0.202708011 * ch$sh_chilcotin_PDO
  ch_out <- ch %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Chilcotin steelhead", spawners = S,
           run_before_pred = S * exp(env - 0.792741195 * ssl_zero_z("chilcotin")),
           recruits_after  = S * exp(env - 0.792741195 * ch$sh_chilcotin_SSL_alt))
  
  bind_rows(chum_out, th_out, ch_out) %>%
    mutate(killed_by_pred = run_before_pred - recruits_after)
}

## Long format for plotting: spawners, before, after, killed
predation_long <- function(pred_df) {
  pred_df %>%
    filter_stock_years() %>%
    transmute(Year, Stock, byrate,
              `Total run before predation`     = run_before_pred,
              `Total recruits after predation` = recruits_after,
              `Brood year spawners`            = spawners,
              `Total run killed by predators`  = killed_by_pred) %>%
    pivot_longer(-c(Year, Stock, byrate), names_to = "series", values_to = "fish") %>%
    mutate(series = factor(series, levels = c("Total run before predation",
                                              "Total recruits after predation",
                                              "Brood year spawners",
                                              "Total run killed by predators")))
}

predation_theme <- list(
  geom_line(linewidth = 0.8),
  scale_color_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#6b6a63")),
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed")),
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)),
  theme_minimal(),
  theme(legend.position = "bottom", panel.grid.minor = element_blank()),
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
)

# ============================================================
# FIGURE 1: run before vs. after predation, SSL control, chum U = 0.2
# (spreadsheet's saved settings), one figure per stock
# ============================================================

predation_sheet <- scenarios_harvest_rate %>%
  filter(abs(harvest_rate - sheet_U_apply) < 1e-9, ssl_scenario == "SSL control") %>%
  compute_predation() %>%
  predation_long()

for (stk in c("Chum", "Thompson steelhead", "Chilcotin steelhead")) {
  p <- ggplot(filter(predation_sheet, Stock == stk),
              aes(Year, fish, color = series, linetype = series)) +
    predation_theme +
    labs(x = "Brood year", y = paste("Total", stk), color = NULL, linetype = NULL,
         title = paste0(stk, " under SSL control (SSLs held at ", SSL_CONTROL_YEAR, " level)"),
         subtitle = paste0("Chum harvest rate ", sheet_U_apply,
                           ", steelhead bycatch rate ", default_byrate))
  print(p)
  ggsave(paste0("figures/", tolower(gsub(" ", "_", stk)), "_run_before_after_predation_SSLcontrol.png"),
         p, width = 10, height = 5, dpi = 600)
}

# ============================================================
# FIGURE 2: same, both steelhead stocks across bycatch rates
# ============================================================

predation_bycatch <- scenarios_bycatch_sslcontrol %>%
  group_split(byrate) %>%
  map_dfr(compute_predation) %>%
  filter(Stock != "Chum") %>%
  predation_long()

p <- ggplot(predation_bycatch, aes(Year, fish, color = series, linetype = series)) +
  predation_theme +
  facet_grid(byrate ~ Stock, scales = "free_y") +
  labs(x = "Brood year", y = "Total steelhead", color = NULL, linetype = NULL) +
  theme(strip.text.y = element_text(angle = 0))
print(p)


ggsave("figures/steelhead_run_before_after_predation_SSLcontrol_byrate_grid.png",
       p, width = 9, height = 16, dpi = 300)

# ============================================================
# FIGURE 3: chum_sh_returns_productivity_historic
# ============================================================

# extract_stock_metrics <- function(scenario_df) {
#   bind_rows(
#     scenario_df %>%
#       transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chum",
#                 productivity = log(chum_recruits_alt / chum_spawners_pred),
#                 returns = chum_recruits_alt),
#     scenario_df %>%
#       transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Thompson steelhead",
#                 productivity = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
#                 returns = sh_thompson_recruits_alt),
#     scenario_df %>%
#       transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chilcotin steelhead",
#                 productivity = log(sh_chilcotin_recruits_alt / (sh_chilcotin_spawners_pred / 1000)),
#                 returns = sh_chilcotin_recruits_alt)
#   )
# }
# 
# plot_scenario_comparison <- function(long_df, metric_col, y_label) {
#   ggplot(filter_stock_years(long_df), aes(Year, .data[[metric_col]], color = ssl_scenario)) +
#     geom_line(linewidth = 1, alpha = 0.85) +
#     facet_wrap(~ Stock, scales = "free") +
#     scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
#     labs(x = "Year", y = y_label, color = NULL) +
#     theme_minimal() +
#     theme(legend.position = "bottom")
# }
# 
# metrics_historic_harvest <- extract_stock_metrics(scenarios_historic_harvest)
# 
# p1 <- plot_scenario_comparison(metrics_historic_harvest, "productivity", "Productivity (alpha)")
# p2 <- plot_scenario_comparison(metrics_historic_harvest, "returns", "Returns")
# (p1 + theme(legend.position = "none") + xlab("")) / p2
#ggsave("figures/chum_sh_returns_productivity_historic.png", width = 10, height = 5, dpi = 600)

# ============================================================
# FIGURE 3: chum_sh_returns_productivity_historic
#
# Productivity is plotted as recruits per spawner with the effect of
# spawner density removed (Ricker: R = S * exp(a + covariates + b*S)):
#   observed  : R(y) / S(y) / exp(b * S(y))
#   predicted : exp(a + covariate effects)   -- scenario sea lion levels
# on a log axis labelled in R/S. Returns are in 000s of fish.
#
# Units: chum S and R are in fish (b is per fish). Steelhead S is stored
# in fish but the Ricker b is per 1000 spawners and R is already in 000s.
# ============================================================

# pull density coefficients from the model function defaults so they stay in sync
coef_of     <- function(fn, arg) eval(formals(fn)[[arg]])
b_chum      <- coef_of(run_chum_scenario,      "spawners_coef")
b_thompson  <- coef_of(run_thompson_scenario,  "sh_thompson_spawners_coef")
b_chilcotin <- coef_of(run_chilcotin_scenario, "sh_chilcotin_spawners_coef")

extract_stock_metrics <- function(scenario_df) {
  bind_rows(
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chum",
                rs_observed  = chum_recruits_alt / chum_spawners_pred /
                  exp(b_chum * chum_spawners_pred),
                rs_predicted = exp(chum_SSL_alpha),       # = base alpha when no SSL control
                returns_000s = chum_recruits_alt / 1000),
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Thompson steelhead",
                S_000s       = sh_thompson_spawners_pred / 1000,
                rs_observed  = sh_thompson_recruits_alt / S_000s / exp(b_thompson * S_000s),
                rs_predicted = exp(sh_thompson_SSL_alpha),
                returns_000s = sh_thompson_recruits_alt) %>%
      select(-S_000s),
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chilcotin steelhead",
                S_000s       = sh_chilcotin_spawners_pred / 1000,
                rs_observed  = sh_chilcotin_recruits_alt / S_000s / exp(b_chilcotin * S_000s),
                rs_predicted = exp(sh_chilcotin_SSL_alpha),
                returns_000s = sh_chilcotin_recruits_alt) %>%
      select(-S_000s)
  )
}

# panel order: chum first, then the two steelhead stocks
# (applied after filter_stock_years(), whose join would drop factor levels)
order_stocks <- function(df) {
  mutate(df, Stock = factor(Stock, levels = c("Chum", "Thompson steelhead", "Chilcotin steelhead")))
}

ssl_colours <- c("No SSL control" = "black", "SSL control" = "#4682B4")

metrics_historic_harvest <- extract_stock_metrics(scenarios_historic_harvest)

## --- Productivity: density-independent R/S, log axis in R/S units ---------
prod_long <- metrics_historic_harvest %>%
  pivot_longer(c(rs_observed, rs_predicted), names_to = "type", values_to = "rs") %>%
  mutate(type = factor(recode(type, rs_observed = "Observed",
                              rs_predicted = "Model predicted"),
                       levels = c("Observed", "Model predicted"))) %>%
  filter(is.finite(rs), rs > 0) %>%
  filter_stock_years() %>%
  order_stocks() %>%
  # draw "No SSL control" last so it isn't hidden where the scenarios coincide
  # (identical up to SSL_CONTROL_YEAR)
  mutate(ssl_scenario = factor(ssl_scenario, levels = c("SSL control", "No SSL control")))

p1 <- ggplot(prod_long, aes(Year, rs, color = ssl_scenario, linetype = type)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  facet_wrap(~ Stock, scales = "free") +
  scale_y_log10(breaks = breaks_log(n = 6),
                labels = label_number(drop0trailing = TRUE, big.mark = ",")) +
  scale_color_manual(values = ssl_colours, breaks = names(ssl_colours)) +
  scale_linetype_manual(values = c("Observed" = "solid", "Model predicted" = "dashed")) +
  labs(x = NULL,
       y = "Recruits per spawner",
       color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

## --- Returns in 000s of fish ----------------------------------------------
p2 <- ggplot(order_stocks(filter_stock_years(metrics_historic_harvest)),
             aes(Year, returns_000s, color = ssl_scenario)) +
  geom_line(linewidth = 1, alpha = 0.85) +
  facet_wrap(~ Stock, scales = "free") +
  scale_color_manual(values = ssl_colours) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Brood year", y = "Returns (000s of fish)", color = NULL) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

p1 / p2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/chum_sh_returns_productivity_historic.png", width = 11, height = 7, dpi = 600)

# ============================================================
# FIGURE 4: chum_mean_catch_lost_historic_and_harvest_sweep_2000-present
# catch lost = "No SSL control" minus "SSL control" (1978 SSL level)
# ============================================================

RECENT_YEARS_CHUM <- 2000:2016

chum_catch_lost <- function(scenario_df) {
  scenario_df %>%
    filter(Year %in% RECENT_YEARS_CHUM) %>%
    select(Year, harvest_rate, ssl_scenario, catch_alt) %>%
    pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
    mutate(catch_lost = `No SSL control` - `SSL control`) %>%
    group_by(harvest_rate) %>%
    summarise(mean_catch_lost = mean(catch_lost, na.rm = TRUE),
              se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")
}

chum_catch_lost_all_bars <- bind_rows(
  chum_catch_lost(scenarios_historic_harvest) %>% mutate(scenario_label = "Historic"),
  chum_catch_lost(scenarios_harvest_rate)     %>% mutate(scenario_label = as.character(harvest_rate))
) %>%
  mutate(scenario_label = factor(scenario_label,
                                 levels = c("Historic", as.character(sort(harvest_rates)))))
print(chum_catch_lost_all_bars)

ggplot(chum_catch_lost_all_bars,
       aes(scenario_label, -mean_catch_lost, fill = scenario_label == "Historic")) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = -mean_catch_lost - se_catch_lost,
                    ymax = -mean_catch_lost + se_catch_lost), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "#FF4500", "FALSE" = "grey50"), guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Harvest rate scenario", y = "Mean yearly catch lost") +
  theme_minimal()
#ggsave("figures/chum_mean_catch_lost_historic_and_harvest_sweep_2000-present.png",
#      width = 10, height = 5.5, dpi = 600)


# ============================================================
# FIGURE 4: chum mean yearly catch, observed vs. SSL control (2000-2016)
#   top    : mean yearly catch -- observed historical record, then SSL
#            control at the historic harvest rate and at each fixed rate
#   bottom : mean yearly catch lost to SSLs = SSL control catch - observed
#            historical catch, for the historic harvest rate (red) and
#            each fixed harvest rate
# Bars are means over RECENT_YEARS_CHUM; intervals are +/- 1 SE across
# years (bottom: SE of the paired yearly differences).
# Requires harvest_rates at 0.05 steps (see CONFIG).
# ============================================================

RECENT_YEARS_CHUM <- 2000:2016

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

fill_colours <- c("Observed"          = "black",
                  "SSL control, historic harvest rate" = "#FF4500",
                  "SSL control, fixed harvest rate"    = "grey50")

lvl_obs  <- "Observed"
lvl_hist <- "Historic\nharvest rate"
rate_lvls <- percent(sort(harvest_rates), accuracy = 1)

## observed historical catch (historic harvest, no SSL control = observed record)
obs_catch <- scenarios_historic_harvest %>%
  filter(ssl_scenario == "No SSL control", Year %in% RECENT_YEARS_CHUM) %>%
  distinct(Year, .keep_all = TRUE) %>%
  select(Year, catch_obs = catch_alt)

## SSL control catch: historic harvest rate + each fixed harvest rate
ctrl_catch <- bind_rows(
  scenarios_historic_harvest %>%
    filter(ssl_scenario == "SSL control") %>%
    mutate(scenario_label = lvl_hist),
  scenarios_harvest_rate %>%
    filter(ssl_scenario == "SSL control") %>%
    mutate(scenario_label = percent(harvest_rate, accuracy = 1))
) %>%
  filter(Year %in% RECENT_YEARS_CHUM) %>%
  distinct(Year, scenario_label, .keep_all = TRUE) %>%
  select(Year, scenario_label, catch = catch_alt) %>%
  left_join(obs_catch, by = "Year") %>%
  mutate(catch_lost = catch - catch_obs)

ctrl_summary <- ctrl_catch %>%
  group_by(scenario_label) %>%
  summarise(mean_catch = mean(catch, na.rm = TRUE),      se_catch = se(catch),
            mean_lost  = mean(catch_lost, na.rm = TRUE), se_lost  = se(catch_lost),
            n_years    = sum(!is.na(catch)),
            .groups = "drop") %>%
  mutate(fill_group = if_else(scenario_label == lvl_hist,
                              "SSL control, historic harvest rate",
                              "SSL control, fixed harvest rate"))

catch_top <- bind_rows(
  tibble(scenario_label = lvl_obs,
         mean_catch = mean(obs_catch$catch_obs, na.rm = TRUE),
         se_catch   = se(obs_catch$catch_obs),
         fill_group = "Observed"),
  ctrl_summary %>% select(scenario_label, mean_catch, se_catch, fill_group)
) %>%
  mutate(scenario_label = factor(scenario_label, levels = c(lvl_obs, lvl_hist, rate_lvls)))

catch_bottom <- ctrl_summary %>%
  mutate(scenario_label = factor(scenario_label, levels = c(lvl_hist, rate_lvls)))

print(catch_top, n = Inf)
print(catch_bottom %>% select(scenario_label, mean_lost, se_lost, n_years), n = Inf)

bar_theme <- list(
  geom_col(width = 0.75),
  # identical limits in both panels so patchwork merges them into one legend
  scale_fill_manual(values = fill_colours, limits = names(fill_colours)),
  scale_y_continuous(labels = label_comma()),
  theme_minimal(),
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
)

## --- Top: mean yearly catch -------------------------------------------------
p_catch <- ggplot(catch_top, aes(scenario_label, mean_catch, fill = fill_group)) +
  bar_theme +
  geom_errorbar(aes(ymin = mean_catch - se_catch, ymax = mean_catch + se_catch), width = 0.25) +
  labs(x = NULL, y = "Mean yearly catch", fill = NULL)

## --- Bottom: mean yearly catch lost (SSL control - observed) ----------------
p_lost <- ggplot(catch_bottom, aes(scenario_label, mean_lost, fill = fill_group)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  bar_theme +
  geom_errorbar(aes(ymin = mean_lost - se_lost, ymax = mean_lost + se_lost), width = 0.25) +
  labs(x = "Harvest rate scenario",
       y = "Mean yearly catch lost", fill = NULL) +
  # the top panel's legend already has all three keys; hide this one
  # (guides(), not theme(), so the patchwork '&' below can't override it)
  guides(fill = "none")

p_catch / p_lost +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(sprintf("figures/chum_mean_catch_observed_vs_SSLcontrol_harvest_sweep_%d-%d.png",
               min(RECENT_YEARS_CHUM), max(RECENT_YEARS_CHUM)),
       width = 11, height = 8, dpi = 600)


# ============================================================
# FIGURE 4: chum mean yearly catch, historical conditions vs. SSL control
# (layout matches the sockeye catch-lost figure)
#   top    : mean yearly catch -- Observed; SSL control at the historic
#            harvest rate; then paired bars at each fixed harvest rate:
#            historical conditions (no SSL control, grey) vs SSL control
#   bottom : mean yearly catch lost = SSL control catch (historic harvest
#            rate and each fixed rate) - OBSERVED historical catch.
#            Grey bars appear only in the top panel, for reference.
# Bars are means over RECENT_YEARS_CHUM; intervals are +/- 1 SE across
# years (bottom: SE of the paired yearly differences).
# Requires harvest_rates at 0.05 steps (see CONFIG).
# ============================================================

RECENT_YEARS_CHUM <- 2000:2016

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

fill_colours <- c("Observed"              = "black",
                  "Historical conditions" = "grey50",
                  "SSL control"           = "#4682B4")

lvl_obs   <- "Observed"
lvl_hist  <- "Historic\nharvest rate"
rate_lvls <- percent(sort(harvest_rates), accuracy = 1)
x_lvls    <- c(lvl_obs, lvl_hist, rate_lvls)

## yearly catch for every scenario x harvest-rate slot ------------------------
## historic harvest: "No SSL control" is the observed record
yearly_catch <- bind_rows(
  scenarios_historic_harvest %>%
    mutate(scenario_label = lvl_hist),
  scenarios_harvest_rate %>%
    mutate(scenario_label = percent(harvest_rate, accuracy = 1))
) %>%
  filter(Year %in% RECENT_YEARS_CHUM) %>%
  distinct(Year, scenario_label, ssl_scenario, .keep_all = TRUE) %>%
  transmute(Year, scenario_label,
            condition = if_else(ssl_scenario == "SSL control", "SSL control", "Historical conditions"),
            catch = catch_alt)

## top panel -----------------------------------------------------------------
catch_top <- yearly_catch %>%
  mutate(fill_group = case_when(
    scenario_label == lvl_hist & condition == "Historical conditions" ~ "Observed",
    TRUE ~ condition),
    scenario_label = if_else(fill_group == "Observed", lvl_obs, scenario_label)) %>%
  group_by(scenario_label, fill_group) %>%
  summarise(mean_catch = mean(catch, na.rm = TRUE), se_catch = se(catch),
            n_years = sum(!is.na(catch)), .groups = "drop") %>%
  mutate(scenario_label = factor(scenario_label, levels = x_lvls),
         fill_group     = factor(fill_group, levels = names(fill_colours))) %>%
  arrange(scenario_label, fill_group)

## bottom panel: SSL control - OBSERVED catch, paired by year ----------------
obs_catch <- yearly_catch %>%
  filter(scenario_label == lvl_hist, condition == "Historical conditions") %>%
  select(Year, catch_obs = catch)

catch_bottom <- yearly_catch %>%
  filter(condition == "SSL control") %>%
  left_join(obs_catch, by = "Year") %>%
  mutate(catch_lost = catch - catch_obs) %>%
  group_by(scenario_label) %>%
  summarise(mean_lost = mean(catch_lost, na.rm = TRUE), se_lost = se(catch_lost),
            n_years = sum(!is.na(catch_lost)), .groups = "drop") %>%
  mutate(fill_group     = factor("SSL control", levels = names(fill_colours)),
         scenario_label = factor(scenario_label, levels = c(lvl_hist, rate_lvls))) %>%
  arrange(scenario_label)

print(catch_top, n = Inf)
print(catch_bottom, n = Inf)

## shared bar styling --------------------------------------------------------
BAR_W <- 0.8
dodge_bar <- position_dodge2(preserve = "single", padding = 0.05)
dodge_err <- position_dodge2(preserve = "single", padding = 0.6)

bar_layers <- list(
  geom_col(width = BAR_W, position = dodge_bar),
  scale_fill_manual(values = fill_colours, limits = names(fill_colours)),
  scale_y_continuous(labels = label_comma()),
  theme_minimal(),
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
)

## --- Top: mean yearly catch -------------------------------------------------
p_catch <- ggplot(catch_top, aes(scenario_label, mean_catch, fill = fill_group)) +
  bar_layers +
  geom_errorbar(aes(ymin = mean_catch - se_catch, ymax = mean_catch + se_catch),
                width = BAR_W, position = dodge_err) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey60") +
  labs(x = NULL, y = "Mean yearly catch", fill = NULL)

## --- Bottom: mean yearly catch lost -----------------------------------------
p_lost <- ggplot(catch_bottom, aes(scenario_label, mean_lost, fill = fill_group)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  bar_layers +
  # one bar per slot here, so no dodging needed -- a fixed narrow cap
  geom_errorbar(aes(ymin = mean_lost - se_lost, ymax = mean_lost + se_lost),
                width = 0.25) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey60") +
  labs(x = "Chum harvest rate scenario",
       y = "Mean yearly catch lost\n(SSL control − observed)", fill = NULL) +
  guides(fill = "none")   # top panel's legend covers all three keys

p_catch / p_lost +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(sprintf("figures/chum_mean_catch_historical_vs_SSLcontrol_harvest_sweep_%d-%d.png",
               min(RECENT_YEARS_CHUM), max(RECENT_YEARS_CHUM)),
       width = 12, height = 8, dpi = 600)

# ============================================================
# FIGURE 4: chum mean yearly catch, historical conditions vs. SSL control -- WITH MSY 
# (layout matches the sockeye catch-lost figure)
#   top    : mean yearly catch -- Observed; SSL control at the historic
#            harvest rate; then paired bars at each fixed harvest rate:
#            historical conditions (no SSL control, grey) vs SSL control
#   bottom : mean yearly catch lost = SSL control catch (historic harvest
#            rate and each fixed rate) - OBSERVED historical catch.
#            Grey bars appear only in the top panel, for reference.
# Bars are means over RECENT_YEARS_CHUM; intervals are +/- 1 SE across
# years (bottom: SE of the paired yearly differences).
# Requires harvest_rates at 0.05 steps (see CONFIG).
# ============================================================

RECENT_YEARS_CHUM <- 2000:2016

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

fill_colours <- c("Observed"              = "black",
                  "Historical conditions" = "grey50",
                  "SSL control"           = "#4682B4")

lvl_obs   <- "Observed"
lvl_hist  <- "Historic\nharvest rate"
rate_lvls <- percent(sort(harvest_rates), accuracy = 1)
x_lvls    <- c(lvl_obs, lvl_hist, rate_lvls)

## yearly catch for every scenario x harvest-rate slot ------------------------
## historic harvest: "No SSL control" is the observed record
yearly_catch <- bind_rows(
  scenarios_historic_harvest %>%
    mutate(scenario_label = lvl_hist),
  scenarios_harvest_rate %>%
    mutate(scenario_label = percent(harvest_rate, accuracy = 1))
) %>%
  filter(Year %in% RECENT_YEARS_CHUM) %>%
  distinct(Year, scenario_label, ssl_scenario, .keep_all = TRUE) %>%
  transmute(Year, scenario_label,
            condition = if_else(ssl_scenario == "SSL control", "SSL control", "Historical conditions"),
            catch = catch_alt)

## top panel -----------------------------------------------------------------
catch_top <- yearly_catch %>%
  mutate(fill_group = case_when(
    scenario_label == lvl_hist & condition == "Historical conditions" ~ "Observed",
    TRUE ~ condition),
    scenario_label = if_else(fill_group == "Observed", lvl_obs, scenario_label)) %>%
  group_by(scenario_label, fill_group) %>%
  summarise(mean_catch = mean(catch, na.rm = TRUE), se_catch = se(catch),
            n_years = sum(!is.na(catch)), .groups = "drop") %>%
  mutate(scenario_label = factor(scenario_label, levels = x_lvls),
         fill_group     = factor(fill_group, levels = names(fill_colours))) %>%
  arrange(scenario_label, fill_group)

## bottom panel: SSL control - OBSERVED catch, paired by year ----------------
obs_catch <- yearly_catch %>%
  filter(scenario_label == lvl_hist, condition == "Historical conditions") %>%
  select(Year, catch_obs = catch)

catch_bottom <- yearly_catch %>%
  filter(condition == "SSL control") %>%
  left_join(obs_catch, by = "Year") %>%
  mutate(catch_lost = catch - catch_obs) %>%
  group_by(scenario_label) %>%
  summarise(mean_lost = mean(catch_lost, na.rm = TRUE), se_lost = se(catch_lost),
            n_years = sum(!is.na(catch_lost)), .groups = "drop") %>%
  mutate(fill_group     = factor("SSL control", levels = names(fill_colours)),
         scenario_label = factor(scenario_label, levels = c(lvl_hist, rate_lvls))) %>%
  arrange(scenario_label)

print(catch_top, n = Inf)
print(catch_bottom, n = Inf)

## MSY reference lines: lowest vs highest productivity periods ---------------
## Productivity = density-independent ln(R/S) = covariate alpha + observed
## residual (what the simulation itself applies each year). Smoothed with a
## rolling MSY_WINDOW-year mean; the lowest and highest windows give a_low
## and a_high. MSY = max over S of [S*exp(a - beta*S) - S] for the Ricker
## (equilibrium surplus production), solved numerically.
MSY_WINDOW <- 10
beta_chum  <- -eval(formals(run_chum_scenario)$spawners_coef)   # positive, per fish

chum_productivity <- scenarios_historic_harvest %>%
  filter(ssl_scenario == "No SSL control") %>%
  distinct(Year, .keep_all = TRUE) %>%
  transmute(Year, a_obs = chum_base_alpha + chum_ln_obs_pred) %>%
  filter(is.finite(a_obs)) %>%
  arrange(Year) %>%
  mutate(a_roll   = zoo::rollmeanr(a_obs, k = MSY_WINDOW, fill = NA),
         win_start = Year - MSY_WINDOW + 1)

ricker_msy <- function(a, beta) {
  if (a <= 0) return(c(Smsy = 0, MSY = 0, Umsy = 0))   # no surplus production
  surplus <- function(S) S * exp(a - beta * S) - S
  opt  <- optimize(surplus, interval = c(0, a / beta), maximum = TRUE)
  Smsy <- opt$maximum
  c(Smsy = Smsy, MSY = opt$objective, Umsy = opt$objective / (opt$objective + Smsy))
}

msy_ref <- bind_rows(
  chum_productivity %>% slice_min(a_roll, n = 1, with_ties = FALSE) %>% mutate(level = "Low productivity"),
  chum_productivity %>% slice_max(a_roll, n = 1, with_ties = FALSE) %>% mutate(level = "High productivity")
) %>%
  rowwise() %>%
  mutate(as_tibble_row(ricker_msy(a_roll, beta_chum))) %>%
  ungroup() %>%
  mutate(label = sprintf("MSY, %s (%d–%d): %s",
                         tolower(level), win_start, Year, comma(round(MSY))))
print(msy_ref %>% select(level, win_start, win_end = Year, a_roll, Smsy, Umsy, MSY))

## shared bar styling --------------------------------------------------------
BAR_W <- 0.8
dodge_bar <- position_dodge2(preserve = "single", padding = 0.05)
dodge_err <- position_dodge2(preserve = "single", padding = 0.6)

bar_layers <- list(
  geom_col(width = BAR_W, position = dodge_bar),
  scale_fill_manual(values = fill_colours, limits = names(fill_colours)),
  scale_y_continuous(labels = label_comma()),
  theme_minimal(),
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
)

## --- Top: mean yearly catch -------------------------------------------------
p_catch <- ggplot(catch_top, aes(scenario_label, mean_catch, fill = fill_group)) +
  bar_layers +
  geom_errorbar(aes(ymin = mean_catch - se_catch, ymax = mean_catch + se_catch),
                width = BAR_W, position = dodge_err) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey60") +
  geom_hline(data = msy_ref, aes(yintercept = MSY), inherit.aes = FALSE,
             linetype = "dashed", color = "firebrick", linewidth = 0.6) +
  labs(x = NULL, y = "Mean yearly catch", fill = NULL)

## --- Bottom: mean yearly catch lost -----------------------------------------
p_lost <- ggplot(catch_bottom, aes(scenario_label, mean_lost, fill = fill_group)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  bar_layers +
  # one bar per slot here, so no dodging needed -- a fixed narrow cap
  geom_errorbar(aes(ymin = mean_lost - se_lost, ymax = mean_lost + se_lost),
                width = 0.25) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey60") +
  labs(x = "Chum harvest rate scenario",
       y = "Mean yearly catch lost\n(SSL control − observed)", fill = NULL) +
  guides(fill = "none")   # top panel's legend covers all three keys

p_catch / p_lost +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(sprintf("figures/chum_mean_catch_historical_vs_SSLcontrol_harvest_sweep_%d-%d.png",
               min(RECENT_YEARS_CHUM), max(RECENT_YEARS_CHUM)),
       width = 12, height = 8, dpi = 600)


# ============================================================
# FIGURE 5: chum_harvest_rate_pred_control_with_observed
# recruits per spawner & returns across harvest rates, observed overlaid.
# Each metric is its own facet_wrap column (fully free y) joined by patchwork.
# ============================================================

chum_harvest_metrics <- scenarios_harvest_rate %>%
  filter(harvest_rate %in% seq(0.1, 0.8, by = 0.1), Year %in% CHUM_YEARS) %>%
  transmute(Year, harvest_rate, ssl_scenario,
            `Recruits per spawner` = chum_recruits_alt / chum_spawners_pred,
            Returns = chum_recruits_alt) %>%
  pivot_longer(c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value")

chum_observed <- scenarios_harvest_rate %>%
  filter(ssl_scenario == "No SSL control", harvest_rate == harvest_rates[1], Year %in% CHUM_YEARS) %>%
  transmute(Year,
            `Recruits per spawner` = chum_recruits_obs / chum_spawners,
            Returns = chum_recruits_obs) %>%
  pivot_longer(c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value") %>%
  mutate(ssl_scenario = "Observed")

chum_sweep_with_obs <- bind_rows(
  chum_harvest_metrics,
  crossing(harvest_rate = unique(chum_harvest_metrics$harvest_rate), chum_observed)
)

plot_chum_metric_stack <- function(df, metric_name) {
  df %>%
    filter(metric == metric_name) %>%
    ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    facet_wrap(~ harvest_rate, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = c("Observed" = "black", "No SSL control" = "#F8766D",
                                  "SSL control" = "#00BFC4")) +
    scale_linetype_manual(values = c("Observed" = "solid", "No SSL control" = "solid",
                                     "SSL control" = "dashed")) +
    labs(x = "Year", y = metric_name, color = NULL, linetype = NULL, title = metric_name) +
    theme_minimal(base_size = 9) +
    theme(strip.placement = "outside")
}

(plot_chum_metric_stack(chum_sweep_with_obs, "Recruits per spawner") + theme(legend.position = "none")) |
  plot_chum_metric_stack(chum_sweep_with_obs, "Returns") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/chum_harvest_rate_pred_control_with_observed.png", width = 9, height = 15, dpi = 600)


# ============================================================
# FIGURE 5: chum_harvest_rate_pred_control_with_observed
# Recruits per spawner (density effect removed) & returns across harvest
# rates, observed overlaid. Each metric is its own facet_wrap column
# (fully free y) joined by patchwork.
#
# Recruits per spawner = R(y) / S(y) / exp(b * S(y)), same index as
# Figure 3, on a log axis labelled in R/S. Because the Ricker density term
# is divided out, this index does not depend on spawner abundance -- so it
# is the same at every harvest rate, and "No SSL control" equals Observed
# exactly. Only the SSL control line differs (after SSL_CONTROL_YEAR).
# ============================================================

b_chum <- eval(formals(run_chum_scenario)$spawners_coef)

fig5_rates <- round(seq(0.1, 0.8, by = 0.1), 2)   # keep 0.1 steps here

rs_label <- "Recruits per spawner"

chum_harvest_metrics <- scenarios_harvest_rate %>%
  filter(round(harvest_rate, 2) %in% fig5_rates, Year %in% CHUM_YEARS) %>%
  transmute(Year, harvest_rate, ssl_scenario,
            !!rs_label := chum_recruits_alt / chum_spawners_pred /
              exp(b_chum * chum_spawners_pred),
            Returns = chum_recruits_alt) %>%
  pivot_longer(c(all_of(rs_label), Returns), names_to = "metric", values_to = "value")

chum_observed <- scenarios_harvest_rate %>%
  filter(ssl_scenario == "No SSL control", harvest_rate == harvest_rates[1], Year %in% CHUM_YEARS) %>%
  transmute(Year,
            !!rs_label := chum_recruits_obs / chum_spawners / exp(b_chum * chum_spawners),
            Returns = chum_recruits_obs) %>%
  pivot_longer(c(all_of(rs_label), Returns), names_to = "metric", values_to = "value") %>%
  mutate(ssl_scenario = "Observed")

chum_sweep_with_obs <- bind_rows(
  chum_harvest_metrics,
  crossing(harvest_rate = unique(chum_harvest_metrics$harvest_rate), chum_observed)
) %>%
  filter(is.finite(value), value > 0) %>%
  # draw order: Observed first, "No SSL control" last (dashed) so it stays
  # visible where it coincides with Observed
  mutate(ssl_scenario = factor(ssl_scenario, levels = c("Observed", "SSL control", "No SSL control")))

scenario_colours   <- c("Observed" = "grey", "No SSL control" = "#F8766D", "SSL control" = "#4682B4")
scenario_linetypes <- c("Observed" = "solid", "No SSL control" = "dashed", "SSL control" = "solid")

plot_chum_metric_stack <- function(df, metric_name, log_y = FALSE) {
  p <- df %>%
    filter(metric == metric_name) %>%
    ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    facet_wrap(~ harvest_rate, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = scenario_colours, breaks = names(scenario_colours)) +
    scale_linetype_manual(values = scenario_linetypes, breaks = names(scenario_colours)) +
    labs(x = "Year", y = metric_name, color = NULL, linetype = NULL, title = metric_name) +
    theme_minimal(base_size = 9) +
    theme(strip.placement = "outside")
  if (log_y) {
    p <- p + scale_y_log10(breaks = breaks_log(n = 4),
                           labels = label_number(drop0trailing = TRUE, big.mark = ","))
  }
  p
}

(plot_chum_metric_stack(chum_sweep_with_obs, rs_label, log_y = TRUE) + theme(legend.position = "none")) |
  plot_chum_metric_stack(chum_sweep_with_obs, "Returns") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/chum_harvest_rate_pred_control_with_observed.png", width = 9, height = 15, dpi = 600)




## NEW with one difference: 
# this one is difference -- FIGURE 2 under historic chum harvest (0.2 is )
# Chum / steelhead deterministic retrospective model
# Haley Oleynik, Murdoch McAllister
#
# Reproduces the 'Retrospective Run1' tab of
# Fraser_Chum_data_v19_alpha_SSL_est_fin_yrs_v12.xlsx and produces:
#   1. Run before vs. after sea lion predation (spawners, run before,
#      recruits after, killed) under SSL control at chum U = 0.2 -- one
#      figure per stock (chum, Thompson, Chilcotin)
#   2. The same for both steelhead stocks across steelhead bycatch rates,
#      with historic chum harvest
#   3. chum_sh_returns_productivity_historic
#   4. chum_mean_catch_lost_historic_and_harvest_sweep_2000-present
#   5. chum_harvest_rate_pred_control_with_observed
#
# Key modeling choices (all match the spreadsheet):
#   * SSL control holds sea lions at the 1978 covariate value for every
#     later brood year (SSL_CONTROL_YEAR).
#   * Steelhead recruits are pure Ricker predictions (no observed residual);
#     chum recruits keep the observed residual (spreadsheet AW: *EXP(N)).
#   * "No SSL control" with historic harvest reproduces the observed record.
#   * Run before predation = Ricker prediction at ZERO sea lions
#     (thompSSLz etc.); after predation = same spawners at scenario SSLs.

library(tidyverse)
library(patchwork)
library(scales)

# ============================================================
# CONFIG
# ============================================================

harvest_rates  <- seq(0, 0.8, by = 0.1)   # chum harvest-rate sweep
bycatch_rates  <- seq(0, 1, by = 0.1)     # steelhead bycatch-rate sweep
default_byrate <- 0.69                    # steelhead bycatch/FN-mortality proxy rate
sheet_U_apply  <- 0.2                     # spreadsheet's saved chum harvest rate (Uapply)

SSL_CONTROL_YEAR <- 1978
# Optional override of the SSL z-score at SSL_CONTROL_YEAR (NA = read from data).
# Spreadsheet 1978 values: chum -0.68849, Thompson -0.90462, Chilcotin -0.87705.
SSL_CONTROL_Z <- c(chum = NA_real_, thompson = NA_real_, chilcotin = NA_real_)

APPLY_SH_RESIDUALS <- FALSE   # spreadsheet uses pure Ricker for steelhead

# Per-stock plotting windows (covariates end 2016)
CHUM_YEARS      <- 1951:2016
THOMPSON_YEARS  <- 1978:2016
CHILCOTIN_YEARS <- 1973:2016
stock_year_bounds <- tibble(
  Stock  = c("Chum", "Thompson steelhead", "Chilcotin steelhead"),
  yr_min = c(min(CHUM_YEARS), min(THOMPSON_YEARS), min(CHILCOTIN_YEARS)),
  yr_max = c(max(CHUM_YEARS), max(THOMPSON_YEARS), max(CHILCOTIN_YEARS))
)

filter_stock_years <- function(df) {
  df %>%
    left_join(stock_year_bounds, by = "Stock") %>%
    filter(Year >= yr_min, Year <= yr_max) %>%
    select(-yr_min, -yr_max)
}

get_control_ssl <- function(df, col, stock) {
  if (!is.na(SSL_CONTROL_Z[[stock]])) return(unname(SSL_CONTROL_Z[[stock]]))
  val <- df %>% filter(Year == SSL_CONTROL_YEAR) %>% pull(all_of(col)) %>% as.numeric()
  val <- val[!is.na(val)]
  if (length(val) == 0) {
    stop(sprintf("No %s value for %d in the input data; set SSL_CONTROL_Z[\"%s\"].",
                 col, SSL_CONTROL_YEAR, stock), call. = FALSE)
  }
  val[1]
}

dir.create("figures", showWarnings = FALSE)

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

## Historic harvest: SSL control vs. no control
scenarios_historic_harvest <- map_dfr(c(0, 1), function(sc) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = 0, SSL_control = sc, U_historic = 1, byrate = default_byrate)
})

## Chum harvest-rate sweep (steelhead byrate fixed at default)
scenarios_harvest_rate <- pmap_dfr(
  expand_grid(SSL_control = c(0, 1), U_apply = harvest_rates),
  function(SSL_control, U_apply) {
    run_full_scenario(data, covariates, sh_data,
                      U_apply = U_apply, SSL_control = SSL_control,
                      U_historic = 0, byrate = default_byrate)
  }
)

## Steelhead bycatch-rate sweep under SSL control with historic chum harvest
## (steelhead bycatch rate = byrate x observed U after 1990; observed before)
scenarios_bycatch_sslcontrol <- map_dfr(bycatch_rates, function(br) {
  run_full_scenario(data, covariates, sh_data,
                    U_apply = 0, SSL_control = 1,
                    U_historic = 1, byrate = br)
})

# ============================================================
# PREDATION MORTALITY ACCOUNTING  (spreadsheet cols FB:GL)
#   run_before_pred = Ricker prediction at zero sea lions   (FD / FU / GC)
#   recruits_after  = same spawners at scenario SSLs        (FF / FV / GE)
#   killed_by_pred  = before - after                        (FG / FW / GF)
# ============================================================

## SSL standardization (covariates tab): zero sea lions = -mean / sd
SSL_SCALE <- list(
  chum      = c(mean = 16088.708673270821, sd = 12075.157324801276),  # V3/V4
  thompson  = c(mean = 23610.490624783917, sd = 14522.985143567075),  # AC3/AC4
  chilcotin = c(mean = 21919.907643703482, sd = 14447.561481342822)   # AD3/AD4
)
ssl_zero_z <- function(stock) unname(-SSL_SCALE[[stock]][["mean"]] / SSL_SCALE[[stock]][["sd"]])

compute_predation <- function(scenario_df) {
  keep_cols <- c("Year", "ssl_scenario", "harvest_rate", "byrate")
  
  ## Chum (FC, FD, FF)
  chum <- scenario_df %>% filter(Year %in% CHUM_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- chum$chum_spawners_pred
  env <- 1.03737862843252 + 0.0929480015696915 * chum$PDO_adult +
    0.102617626088713 * chum$NPGO - 0.105950753725792 * chum$PDO_smolt +
    -4.95136622626478E-07 * S
  chum_out <- chum %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Chum", spawners = S,
           run_before_pred = S * exp(env - 0.224246284912103 * ssl_zero_z("chum")),
           recruits_after  = S * exp(env - 0.224246284912103 * chum$chum_SSL_alt))
  
  ## Thompson (FT, FU, FV) -- spreadsheet omits NPGO here
  th  <- scenario_df %>% filter(Year %in% THOMPSON_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- th$sh_thompson_spawners_pred
  env <- 1.572107637 - 0.804438793 * S / 1000 - 0.203463091 * th$sh_thompson_SST
  th_out <- th %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Thompson steelhead", spawners = S,
           run_before_pred = S * exp(env - 0.764277677 * ssl_zero_z("thompson")),
           recruits_after  = S * exp(env - 0.764277677 * th$sh_thompson_SSL_alt))
  
  ## Chilcotin (GB, GC, GE)
  ch  <- scenario_df %>% filter(Year %in% CHILCOTIN_YEARS) %>% distinct(Year, .keep_all = TRUE)
  S   <- ch$sh_chilcotin_spawners_pred
  env <- 1.053608979 - 1.022467631 * S / 1000 - 0.127949278 * ch$sh_chilcotin_SST +
    0.152526045 * ch$sh_chilcotin_NPGO + 0.202708011 * ch$sh_chilcotin_PDO
  ch_out <- ch %>%
    select(all_of(keep_cols)) %>%
    mutate(Stock = "Chilcotin steelhead", spawners = S,
           run_before_pred = S * exp(env - 0.792741195 * ssl_zero_z("chilcotin")),
           recruits_after  = S * exp(env - 0.792741195 * ch$sh_chilcotin_SSL_alt))
  
  bind_rows(chum_out, th_out, ch_out) %>%
    mutate(killed_by_pred = run_before_pred - recruits_after)
}

## Long format for plotting: spawners, before, after, killed
predation_long <- function(pred_df) {
  pred_df %>%
    filter_stock_years() %>%
    transmute(Year, Stock, byrate,
              `Total run before predation`     = run_before_pred,
              `Total recruits after predation` = recruits_after,
              `Brood year spawners`            = spawners,
              `Total run killed by predators`  = killed_by_pred) %>%
    pivot_longer(-c(Year, Stock, byrate), names_to = "series", values_to = "fish") %>%
    mutate(series = factor(series, levels = c("Total run before predation",
                                              "Total recruits after predation",
                                              "Brood year spawners",
                                              "Total run killed by predators")))
}

predation_theme <- list(
  geom_line(linewidth = 0.8),
  scale_color_manual(values = c("#2a78d6", "#eb6834", "#1baf7a", "#6b6a63")),
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed")),
  scale_y_continuous(labels = scales::comma, limits = c(0, NA)),
  theme_minimal(),
  theme(legend.position = "bottom", panel.grid.minor = element_blank()),
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
)

# ============================================================
# FIGURE 1: run before vs. after predation, SSL control, chum U = 0.2
# (spreadsheet's saved settings), one figure per stock
# ============================================================

predation_sheet <- scenarios_harvest_rate %>%
  filter(abs(harvest_rate - sheet_U_apply) < 1e-9, ssl_scenario == "SSL control") %>%
  compute_predation() %>%
  predation_long()

for (stk in c("Chum", "Thompson steelhead", "Chilcotin steelhead")) {
  p <- ggplot(filter(predation_sheet, Stock == stk),
              aes(Year, fish, color = series, linetype = series)) +
    predation_theme +
    labs(x = "Brood year", y = paste("Total", stk), color = NULL, linetype = NULL,
         title = paste0(stk, " under SSL control (SSLs held at ", SSL_CONTROL_YEAR, " level)"),
         subtitle = paste0("Chum harvest rate ", sheet_U_apply,
                           ", steelhead bycatch rate ", default_byrate))
  print(p)
  ggsave(paste0("figures/", tolower(gsub(" ", "_", stk)), "_run_before_after_predation_SSLcontrol.png"),
         p, width = 10, height = 5, dpi = 600)
}

# ============================================================
# FIGURE 2: same, both steelhead stocks across bycatch rates,
# historic chum harvest
# ============================================================

predation_bycatch <- scenarios_bycatch_sslcontrol %>%
  group_split(byrate) %>%
  map_dfr(compute_predation) %>%
  filter(Stock != "Chum") %>%
  predation_long()

p <- ggplot(predation_bycatch, aes(Year, fish, color = series, linetype = series)) +
  predation_theme +
  facet_grid(byrate ~ Stock, scales = "free") +
  labs(x = "Brood year", y = "Total steelhead", color = NULL, linetype = NULL) +
  theme(strip.text.y = element_text(angle = 0))
print(p)
ggsave("figures/steelhead_run_before_after_predation_SSLcontrol_byrate_grid.png",
       p, width = 9, height = 16, dpi = 300)

# ============================================================
# FIGURE 3: chum_sh_returns_productivity_historic
# ============================================================

extract_stock_metrics <- function(scenario_df) {
  bind_rows(
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chum",
                productivity = log(chum_recruits_alt / chum_spawners_pred),
                returns = chum_recruits_alt),
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Thompson steelhead",
                productivity = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
                returns = sh_thompson_recruits_alt),
    scenario_df %>%
      transmute(Year, ssl_scenario, harvest_rate, byrate, Stock = "Chilcotin steelhead",
                productivity = log(sh_chilcotin_recruits_alt / (sh_chilcotin_spawners_pred / 1000)),
                returns = sh_chilcotin_recruits_alt)
  )
}

plot_scenario_comparison <- function(long_df, metric_col, y_label) {
  ggplot(filter_stock_years(long_df), aes(Year, .data[[metric_col]], color = ssl_scenario)) +
    geom_line(linewidth = 1, alpha = 0.85) +
    facet_wrap(~ Stock, scales = "free") +
    scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
    labs(x = "Year", y = y_label, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

metrics_historic_harvest <- extract_stock_metrics(scenarios_historic_harvest)

p1 <- plot_scenario_comparison(metrics_historic_harvest, "productivity", "Productivity (alpha)")
p2 <- plot_scenario_comparison(metrics_historic_harvest, "returns", "Returns")
(p1 + theme(legend.position = "none") + xlab("")) / p2
ggsave("figures/chum_sh_returns_productivity_historic.png", width = 10, height = 5, dpi = 600)

# ============================================================
# FIGURE 4: chum_mean_catch_lost_historic_and_harvest_sweep_2000-present
# catch lost = "No SSL control" minus "SSL control" (1978 SSL level)
# ============================================================

RECENT_YEARS_CHUM <- 2000:2016

chum_catch_lost <- function(scenario_df) {
  scenario_df %>%
    filter(Year %in% RECENT_YEARS_CHUM) %>%
    select(Year, harvest_rate, ssl_scenario, catch_alt) %>%
    pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
    mutate(catch_lost = `No SSL control` - `SSL control`) %>%
    group_by(harvest_rate) %>%
    summarise(mean_catch_lost = mean(catch_lost, na.rm = TRUE),
              se_catch_lost   = sd(catch_lost, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")
}

chum_catch_lost_all_bars <- bind_rows(
  chum_catch_lost(scenarios_historic_harvest) %>% mutate(scenario_label = "Historic"),
  chum_catch_lost(scenarios_harvest_rate)     %>% mutate(scenario_label = as.character(harvest_rate))
) %>%
  mutate(scenario_label = factor(scenario_label,
                                 levels = c("Historic", as.character(sort(harvest_rates)))))
print(chum_catch_lost_all_bars)

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
# FIGURE 5: chum_harvest_rate_pred_control_with_observed
# recruits per spawner & returns across harvest rates, observed overlaid.
# Each metric is its own facet_wrap column (fully free y) joined by patchwork.
# ============================================================

chum_harvest_metrics <- scenarios_harvest_rate %>%
  filter(harvest_rate %in% seq(0.1, 0.8, by = 0.1), Year %in% CHUM_YEARS) %>%
  transmute(Year, harvest_rate, ssl_scenario,
            `Recruits per spawner` = chum_recruits_alt / chum_spawners_pred,
            Returns = chum_recruits_alt) %>%
  pivot_longer(c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value")

chum_observed <- scenarios_harvest_rate %>%
  filter(ssl_scenario == "No SSL control", harvest_rate == harvest_rates[1], Year %in% CHUM_YEARS) %>%
  transmute(Year,
            `Recruits per spawner` = chum_recruits_obs / chum_spawners,
            Returns = chum_recruits_obs) %>%
  pivot_longer(c(`Recruits per spawner`, Returns), names_to = "metric", values_to = "value") %>%
  mutate(ssl_scenario = "Observed")

chum_sweep_with_obs <- bind_rows(
  chum_harvest_metrics,
  crossing(harvest_rate = unique(chum_harvest_metrics$harvest_rate), chum_observed)
)

plot_chum_metric_stack <- function(df, metric_name) {
  df %>%
    filter(metric == metric_name) %>%
    ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    facet_wrap(~ harvest_rate, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_color_manual(values = c("Observed" = "black", "No SSL control" = "#F8766D",
                                  "SSL control" = "#00BFC4")) +
    scale_linetype_manual(values = c("Observed" = "solid", "No SSL control" = "solid",
                                     "SSL control" = "dashed")) +
    labs(x = "Year", y = metric_name, color = NULL, linetype = NULL, title = metric_name) +
    theme_minimal(base_size = 9) +
    theme(strip.placement = "outside")
}

(plot_chum_metric_stack(chum_sweep_with_obs, "Recruits per spawner") + theme(legend.position = "none")) |
  plot_chum_metric_stack(chum_sweep_with_obs, "Returns") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/chum_harvest_rate_pred_control_with_observed.png", width = 9, height = 15, dpi = 600)