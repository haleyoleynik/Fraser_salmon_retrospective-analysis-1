# Chum / steelhead retrospective model — scenario comparisons
# Haley Oleynik, Murdoch McAllister

# USE steelhead_chum_retrospective_deterministic.R 














# Cleaned version. Keeps the three core scenario analyses:
#   1. Historic harvest rate: SSL control vs. no control
#   2. Harvest-rate sweep (byrate held fixed): SSL control vs. no control
#   3. Bycatch-rate sweep (harvest held at historic): SSL control vs. no control
# ...plus chum catch-lost-to-pinnipeds averaging, 2000-present.
#
# Dropped from the original: dead/commented-out code, duplicate function
# definitions, the SSL-predation-rate exploration section, and the
# lost-value/CPI section (not part of this request).

library(tidyverse)
library(readr)
library(patchwork)
library(scales)

# ============================================================
# CONFIG
# ============================================================

MODEL_YEARS <- 1980:2016   # reliable window used throughout for plotting
harvest_rates <- seq(0, 0.8, by = 0.1)
bycatch_rates <- seq(0, 1, by = 0.1)
default_byrate <- 0.69     # current bycatch/FN-mortality proxy rate for steelhead

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
  
  SSL_1978 <- new.data %>% filter(Year == 1978) %>% pull(SSL)
  
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
    
    df$chum_SSL_alt[i] <- if (df$Year[i] <= 1978) {
      df$SSL[i]
    } else {
      (1 - SSL_control) * df$SSL[i] + SSL_control * SSL_1978
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
    
    df$chum_recruits_alt[i] <-
      if (SSL_control == 0) {
        df$chum_spawners_pred[i] *
          exp(df$chum_base_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
          exp(df$chum_ln_obs_pred[i])
      } else {
        df$chum_spawners_pred[i] *
          exp(df$chum_SSL_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
          exp(df$chum_ln_obs_pred[i])
      }
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
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  FN_thompson_2018 <- all_data %>% filter(Year == 2018) %>% pull(sh_thompson_FN_mortalities) %>% as.numeric()
  
  df <- all_data %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_thompson_base_alpha = sh_thompson_intercept +
        sh_thompson_SST  * sh_thompson_sst_coef +
        sh_thompson_SL   * sh_thompson_ssl_coef +
        sh_thompson_NPGO * sh_thompson_npgo_coef,
      #sh_thompson_model_recruits = sh_thompson_spawners *
      #  exp(sh_thompson_base_alpha + sh_thompson_spawners * sh_thompson_spawners_coef),
      #sh_thompson_ln_obs_pred = log(sh_thompson_recruits / sh_thompson_model_recruits),
      sh_thompson_model_recruits =
        (sh_thompson_spawners / 1000) *
        exp(sh_thompson_base_alpha +
              (sh_thompson_spawners / 1000) * sh_thompson_spawners_coef),
      sh_thompson_ln_obs_pred =
        log((sh_thompson_recruits / 1000) / sh_thompson_model_recruits),
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
    
    df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) {
      df$sh_thompson_SL[i]
    } else {
      (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
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
      
      # df$sh_thompson_recruits_alt[i] <-
      #   if (SSL_control == 0) {
      #     spk * exp(df$sh_thompson_alpha_CN[i] + sh_thompson_spawners_coef * spk)
      #   } else {
      #     spk * exp(df$sh_thompson_SSL_alpha[i] + sh_thompson_spawners_coef * spk)
      #   }
      
      df$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_thompson_alpha_CN[i] +
                      sh_thompson_spawners_coef * spk) *
            exp(df$sh_thompson_ln_obs_pred[i])
        } else {
          spk * exp(df$sh_thompson_SSL_alpha[i] +
                      sh_thompson_spawners_coef * spk) *
            exp(df$sh_thompson_ln_obs_pred[i])
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
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
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
      # sh_chilcotin_model_recruits = sh_chilcotin_spawners *
      #   exp(sh_chilcotin_base_alpha + sh_chilcotin_spawners * sh_chilcotin_spawners_coef),
      # sh_chilcotin_ln_obs_pred = log(sh_chilcotin_recruits / sh_chilcotin_model_recruits),
      
      sh_chilcotin_model_recruits =
        (sh_chilcotin_spawners / 1000) *
        exp(sh_chilcotin_base_alpha +
              (sh_chilcotin_spawners / 1000) * sh_chilcotin_spawners_coef),
      
      sh_chilcotin_ln_obs_pred =
        log((sh_chilcotin_recruits / 1000) / sh_chilcotin_model_recruits),
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
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) {
      df$sh_chilcotin_SL[i]
    } else {
      (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
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
      
      # df$sh_chilcotin_recruits_alt[i] <-
      #   if (SSL_control == 0) {
      #     spk * exp(df$sh_chilcotin_alpha_CN[i] + sh_chilcotin_spawners_coef * spk)
      #   } else {
      #     spk * exp(df$sh_chilcotin_SSL_alpha[i] + sh_chilcotin_spawners_coef * spk)
      #   }
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_chilcotin_alpha_CN[i] +
                      sh_chilcotin_spawners_coef * spk) *
            exp(df$sh_chilcotin_ln_obs_pred[i])
        } else {
          spk * exp(df$sh_chilcotin_SSL_alpha[i] +
                      sh_chilcotin_spawners_coef * spk) *
            exp(df$sh_chilcotin_ln_obs_pred[i])
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

plot_scenario_comparison <- function(long_df, metric_col, y_label, title, year_range = MODEL_YEARS) {
  df <- long_df %>% filter(Year %in% year_range)
  ggplot(df, aes(Year, .data[[metric_col]], color = ssl_scenario)) +
    geom_line(linewidth = 1, alpha = 0.85) +
    facet_wrap(~ Stock, scales = "free_y") +
    scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
    labs(x = "Year", y = y_label, color = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_scenario_sweep <- function(long_df, metric_col, group_var, group_label, y_label, title, year_range = MODEL_YEARS) {
  df <- long_df %>% filter(Year %in% year_range)
  ggplot(df, aes(Year, .data[[metric_col]], color = factor(.data[[group_var]]), linetype = ssl_scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    facet_wrap(~ Stock, scales = "free_y") +
    scale_color_viridis_d(name = group_label) +
    labs(x = "Year", y = y_label, linetype = "SSL scenario") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

## 1. Historic harvest: productivity & returns --------------------------------
(p1 <- plot_scenario_comparison(metrics_historic_harvest, "productivity", "Productivity (alpha)",
                         "Productivity under historic harvest: SSL control vs. none"))


ggsave("figures/chum_sh_productivity_historic.png", width = 11, height = 5, dpi = 600)

p2 <- plot_scenario_comparison(metrics_historic_harvest, "returns", "Returns",
                         "Returns under historic harvest: SSL control vs. none")

(p1 + theme(legend.position = "none") + xlab("")) / p2

ggsave("figures/chum_sh_returns_productivity_historic.png", width = 10, height = 5, dpi = 600)

## 2. Harvest-rate sweep: productivity & returns ------------------------------
#plot_scenario_sweep(metrics_harvest_rate, "productivity", "harvest_rate", "Harvest rate",
#                    "Productivity (alpha)", "Productivity across chum harvest-rate scenarios")
#ggsave("figures/chum_sh_productivity_harvest_sweep.png", width = 11, height = 5, dpi = 600)

plot_scenario_sweep(metrics_harvest_rate, "returns", "harvest_rate", "Harvest rate",
                    "Returns", "Returns across chum harvest-rate scenarios")
ggsave("figures/chum_sh_returns_harvest_sweep.png", width = 11, height = 5, dpi = 600)

## 3. Bycatch-rate sweep: productivity & returns ------------------------------
plot_scenario_sweep(metrics_bycatch_rate, "productivity", "byrate", "Bycatch rate",
                    "Productivity (alpha)", "Productivity across steelhead bycatch-rate scenarios")
ggsave("figures/chum_sh_productivity_bycatch_sweep.png", width = 11, height = 5, dpi = 600)

plot_scenario_sweep(metrics_bycatch_rate, "returns", "byrate", "Bycatch rate",
                    "Returns", "Returns across steelhead bycatch-rate scenarios")
ggsave("figures/chum_sh_returns_bycatch_sweep.png", width = 11, height = 5, dpi = 600)

# ============================================================
# CHUM CATCH LOST TO PINNIPEDS (historic harvest scenario)
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

plot_harvest_rate_grid <- function(long_df, stock_name, year_range = MODEL_YEARS) {
  df <- long_df %>% filter(Stock == stock_name, Year %in% year_range)
  
  ggplot(df, aes(Year, value, color = ssl_scenario)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    facet_grid(harvest_rate ~ metric, scales = "free_y") +
    scale_color_manual(values = c("No SSL control" = "black", "SSL control" = "#4682B4")) +
    labs(x = "Year", y = NULL, color = NULL,
         title = paste0(stock_name, ": productivity & returns across harvest rates (0.1-0.8)")) +
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
# with historical observed overlaid -- FIXED VERSION
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
    filter(byrate <= 0.8, Year %in% MODEL_YEARS) %>%
    transmute(Year, byrate, ssl_scenario, Stock = "Thompson steelhead",
              productivity = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
              returns = sh_thompson_sum_pred),
  scenarios_bycatch_rate %>%
    filter(byrate <= 0.8, Year %in% MODEL_YEARS) %>%
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
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% MODEL_YEARS) %>%
    transmute(Year, Stock = "Thompson steelhead",
              productivity = sh_thompson_alpha_CN,
              returns = sh_thompson_prefishery_N),
  scenarios_bycatch_rate %>%
    filter(ssl_scenario == "No SSL control", byrate == 0, Year %in% MODEL_YEARS) %>%
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

# ============================================================
# STEELHEAD COSEWIC REPLICATION -- Thompson & Chilcotin, matched
# to the 2020 assessment. Chum components dropped entirely.
#
# Method notes (differs from the coho version):
#   - No smoothing. Figures 9/13 fit log-linear regression directly
#     on annual spawner counts.
#   - Generation time differs by DU: Thompson = 5 yr, Chilcotin = 6 yr.
#     "3 generations" is therefore a different window length for each.
#   - % change = 100*(exp(slope * (end_year - start_year)) - 1), using
#     elapsed years between window endpoints (not a fixed 10-yr
#     conversion like the coho script used).
#   - Metric = spawners (sh_thompson_spawners / sh_chilcotin_spawners),
#     matching "number of mature adults (spawners)" in Figs 9/13.
#
# Windows (from Figures 9 & 13):
#   Thompson:  entire 1978-2020, recent 2006-2020
#   Chilcotin: entire 1972-2020, recent 2003-2020
#
# Depends on:
#   - `sh_data` with sh_thompson_spawners / sh_chilcotin_spawners by
#     Year -- to VALIDATE against the report's 82% (Thompson) /
#     80% (Chilcotin) 3-generation decline numbers
#   - `run_thompson_scenario()` and `run_chilcotin_scenario()` from
#     your chum/steelhead retrospective script (unmodified)
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

STEELHEAD_WINDOWS <- list(
  Thompson = list(gen_time = 5, entire_start = 1978, entire_end = 2020,
                  recent_start = 2006, recent_end = 2020),
  Chilcotin = list(gen_time = 6, entire_start = 1972, entire_end = 2020,
                   recent_start = 2003, recent_end = 2020)
)

default_byrate <- 0.69   # unchanged from the original script

# ---- COSEWIC replication regression: fixed window, raw annual data --
cosewic_pct_change_steelhead <- function(df, year_col, value_col,
                                         start_year, end_year) {
  df <- df %>%
    rename(yr = {{ year_col }}, val = {{ value_col }}) %>%
    filter(yr >= start_year, yr <= end_year, !is.na(val), val > 0) %>%
    arrange(yr)
  
  if (nrow(df) < 4) {
    return(tibble(slope = NA_real_, p_value = NA_real_, r2 = NA_real_,
                  n = nrow(df), pct_change = NA_real_))
  }
  
  fit <- lm(log(val) ~ yr, data = df)
  b   <- coef(fit)[["yr"]]
  tibble(
    slope      = b,
    p_value    = summary(fit)$coefficients["yr", "Pr(>|t|)"],
    r2         = summary(fit)$r.squared,
    n          = nrow(df),
    pct_change = 100 * (exp(b * (end_year - start_year)) - 1)
  )
}

classify_cosewic_A <- function(pct_change) {
  case_when(
    is.na(pct_change) ~ NA_character_,
    pct_change <= -50 ~ "Endangered",
    pct_change <= -30 ~ "Threatened",
    TRUE              ~ "Not at Risk / Special Concern"
  )
}

# ============================================================
# STEP 1 -- VALIDATE against the report's own numbers
# Thompson recent-window pct_change should land near -82%;
# Chilcotin near -80%.
# ============================================================

validation <- bind_rows(
  cosewic_pct_change_steelhead(sh_data, Year, sh_thompson_spawners,
                               start_year = STEELHEAD_WINDOWS$Thompson$entire_start,
                               end_year   = STEELHEAD_WINDOWS$Thompson$entire_end) %>%
    mutate(Stock = "Thompson", period = "Entire series (1978-2020)"),
  cosewic_pct_change_steelhead(sh_data, Year, sh_thompson_spawners,
                               start_year = STEELHEAD_WINDOWS$Thompson$recent_start,
                               end_year   = STEELHEAD_WINDOWS$Thompson$recent_end) %>%
    mutate(Stock = "Thompson", period = "Recent 3 gen, 2006-2020 (expect ~ -82%)"),
  cosewic_pct_change_steelhead(sh_data, Year, sh_chilcotin_spawners,
                               start_year = STEELHEAD_WINDOWS$Chilcotin$entire_start,
                               end_year   = STEELHEAD_WINDOWS$Chilcotin$entire_end) %>%
    mutate(Stock = "Chilcotin", period = "Entire series (1972-2020)"),
  cosewic_pct_change_steelhead(sh_data, Year, sh_chilcotin_spawners,
                               start_year = STEELHEAD_WINDOWS$Chilcotin$recent_start,
                               end_year   = STEELHEAD_WINDOWS$Chilcotin$recent_end) %>%
    mutate(Stock = "Chilcotin", period = "Recent 3 gen, 2003-2020 (expect ~ -80%)")
)
print(validation)

# ============================================================
# STEP 2 -- Run the SSL-control scenario (steelhead only, no chum)
#
# U_historic = 1 for both stocks -> sh_*_U_comm uses each stock's own
# reconstructed historic exploitation rate, so chum_commercial_harvest_
# uapply is never evaluated. It's included as an unused placeholder
# column only because run_thompson_scenario()/run_chilcotin_scenario()
# reference it inside an untaken branch.
# ============================================================

run_full_steelhead_scenario <- function(sh_data, SSL_control, byrate = default_byrate) {
  all_data <- sh_data %>%
    mutate(chum_commercial_harvest_uapply = NA_real_)
  
  thompson_df <- run_thompson_scenario(all_data, SSL_control, U_historic = 1, byrate = byrate)
  full_df     <- run_chilcotin_scenario(thompson_df, SSL_control, U_historic = 1, byrate = byrate)
  
  full_df %>%
    mutate(ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control"))
}

scenarios_steelhead_historic <- purrr::map_dfr(c(0, 1), function(sc) {
  run_full_steelhead_scenario(sh_data, SSL_control = sc, byrate = default_byrate)
})

modeled_spawners <- bind_rows(
  scenarios_steelhead_historic %>%
    transmute(Year, ssl_scenario, Stock = "Thompson", spawners_pred = sh_thompson_spawners_pred),
  scenarios_steelhead_historic %>%
    transmute(Year, ssl_scenario, Stock = "Chilcotin", spawners_pred = sh_chilcotin_spawners_pred)
)

status_by_scenario <- modeled_spawners %>%
  group_by(Stock, ssl_scenario) %>%
  group_modify(~ {
    win <- STEELHEAD_WINDOWS[[unique(.y$Stock)]]
    cosewic_pct_change_steelhead(.x, Year, spawners_pred,
                                 start_year = win$recent_start, end_year = win$recent_end)
  }) %>%
  ungroup() %>%
  mutate(status = classify_cosewic_A(pct_change))

print(status_by_scenario)

status_comparison_table <- status_by_scenario %>%
  select(Stock, ssl_scenario, status) %>%
  pivot_wider(names_from = ssl_scenario, values_from = status)
print(status_comparison_table)

# ============================================================
# STEP 3 -- PLOTS
# ============================================================

SCENARIO_COLORS_SH <- c("Pinniped scenario" = "#4682B4")

scenario_recode <- c("No SSL control" = "Observed", "SSL control" = "Pinniped scenario")

modeled_spawners <- modeled_spawners %>%
  mutate(scenario = recode(ssl_scenario, !!!scenario_recode))

observed_spawners_by_year <- modeled_spawners %>%
  filter(scenario == "Observed") %>%
  select(Year, Stock, Observed = spawners_pred)

## % increase vs. observed, by year and stock -------------------------------
pct_increase_by_year <- modeled_spawners %>%
  filter(scenario != "Observed") %>%
  left_join(observed_spawners_by_year, by = c("Year", "Stock")) %>%
  mutate(pct_increase = 100 * (spawners_pred - Observed) / Observed) %>%
  select(Year, Stock, scenario, pct_increase)

## Assessment window per stock (for shading + the bar-chart summary) --------
assessment_windows <- bind_rows(
  tibble(Stock = "Thompson", win_start = STEELHEAD_WINDOWS$Thompson$recent_start,
         win_end = STEELHEAD_WINDOWS$Thompson$recent_end),
  tibble(Stock = "Chilcotin", win_start = STEELHEAD_WINDOWS$Chilcotin$recent_start,
         win_end = STEELHEAD_WINDOWS$Chilcotin$recent_end)
)

pct_increase_over_window <- pct_increase_by_year %>%
  inner_join(assessment_windows, by = "Stock") %>%
  filter(Year >= win_start, Year <= win_end) %>%
  group_by(Stock, scenario) %>%
  summarise(pct_increase = mean(pct_increase, na.rm = TRUE), .groups = "drop")

## Bar chart: % increase in mean spawners over the 3-generation
## COSEWIC assessment window, by stock and scenario
ggplot(pct_increase_over_window,
       aes(x = reorder(Stock, -pct_increase), y = pct_increase, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS_SH) +
  labs(x = NULL, y = "% increase in mean spawners over the 3-generation assessment window",
       fill = "Scenario driver",
       title = "Steelhead: potential spawner recovery under a pinniped-control scenario") +
  theme_minimal()
ggsave("figures/steelhead_pct_increase_over_window_by_scenario.png", width = 7, height = 6, dpi = 600)

## Time series: % increase over time, faceted by stock, assessment
## window shaded
ggplot(pct_increase_by_year, aes(Year, pct_increase, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_rect(data = assessment_windows,
            aes(xmin = win_start, xmax = win_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey70", alpha = 0.25) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_color_manual(values = SCENARIO_COLORS_SH) +
  labs(x = "Year", y = "% increase in spawners vs. observed", color = "Scenario driver",
       title = "Steelhead: counterfactual spawner gain over time, by stock",
       subtitle = "Shaded band marks each DU's COSEWIC 3-generation assessment window") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("figures/steelhead_pct_increase_over_time_by_scenario.png", width = 11, height = 6, dpi = 600)

## Status-comparison heatmap
status_plot_df <- status_by_scenario %>%
  mutate(
    scenario = recode(ssl_scenario, !!!scenario_recode),
    scenario = factor(scenario, levels = c("Observed", "Pinniped scenario")),
    status   = factor(status, levels = c("Endangered", "Threatened",
                                         "Not at Risk / Special Concern"))
  )

ggplot(status_plot_df, aes(scenario, Stock, fill = status)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "Endangered" = "#B22222",
    "Threatened" = "#E8A33D",
    "Not at Risk / Special Concern" = "#4C9A5B"
  )) +
  labs(x = NULL, y = NULL, fill = "Status",
       title = "Steelhead: modeled COSEWIC status band by DU and scenario\n(3-generation assessment window)") +
  theme_minimal()
ggsave("figures/steelhead_status_by_scenario_heatmap.png", width = 7, height = 4.5, dpi = 600)
