# ============================================================
# Probabilistic (Monte Carlo) retrospective analysis --
# Chum, Thompson steelhead, Chilcotin steelhead
#
# CHANGES FROM THE PREVIOUS VERSION:
#  1. Chum is now drawn PROBABILISTICALLY, same as Thompson/Chilcotin --
#     previously it ran once outside the draw loop with fixed point-
#     estimate coefficients (intercept/pdo_adult_coef/npgo_coef/
#     pdo_smolt_coef/ssl_coef/spawners_coef), identical across every
#     draw. Now (a, b, s, t) are sampled from chum_posterior_samples.rds
#     each draw, same as the other two stocks.
#  2. Chum's covariate set changed to match the new dredge results:
#     Spawners + SL + a SINGLE PDO covariate (dropping PDO_adult,
#     PDO_smolt, and NPGO entirely -- per your instruction, "just use
#     one PDO covariate ... forget about smolt and adult"). The single
#     PDO column comes from chum_data.csv (the same source the new JAGS
#     fit used), joined in by Year -- covariates.csv only ever had
#     PDO_adult/PDO_smolt, no plain PDO, so this is a new join, not a
#     column that already existed in the old covariate data.
#  3. run_chilcotin_scenario_v2() no longer includes max_flow -- the
#     updated dredge results dropped it entirely for Chilcotin, making
#     its covariate set (Spawners + SL + SST) identical in structure to
#     Thompson now. This was still present in the previous version of
#     this script; removed here.
#  4. Because chum feeds into Thompson/Chilcotin via
#     chum_commercial_harvest_uapply (used when U_historic == 0),
#     making chum probabilistic means Thompson/Chilcotin's results now
#     inherit chum's posterior uncertainty too, not just their own --
#     worth remembering when interpreting Chilcotin/Thompson results
#     specifically. This retrospective doesn't use U_historic == 0
#     scenarios (see the historic-harvest-only design below), so this
#     doesn't currently bite, but would if the harvest-rate/bycatch-rate
#     sweeps from the original scenario script were revived here.
#  5. Added the X/Xmin analysis (task 3), matching the coho/sockeye
#     design: Xmin = fixed, from real observed abundance's own 10-year
#     low-period mean per stock; X = mean SCENARIO-reconstructed
#     abundance over the final 10 years, per posterior draw; ratio =
#     X / Xmin, summarized via pointrange + natural log scale.
#
# DEPENDENCIES -- run the ORIGINAL chum_steelhead_retrospective_
# scenarios.R (the cleaned version, through its DATA LOADING section)
# first, so `data`, `covariates`, `sh_data` exist. Also run
# steelhead_chum_jags_models.R (or readRDS the saved posterior samples)
# so `chilcotin_samples_kept` / `thompson_samples_kept` /
# `chum_samples_kept` exist.
# ============================================================

# need to run steelhead_retrospective-model.R first 

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50   # start small, matching the pattern used throughout this pipeline
LOW_PERIOD_YEARS <- 10
FINAL_PERIOD_YEARS <- 10

# ------------------------------------------------------------
# 1. LOAD POSTERIOR SAMPLES, ALL THREE STOCKS
# ------------------------------------------------------------

if (!exists("chilcotin_samples_kept")) chilcotin_samples_kept <- readRDS("chilcotin_posterior_samples.rds")
if (!exists("thompson_samples_kept"))  thompson_samples_kept  <- readRDS("thompson_posterior_samples.rds")
if (!exists("chum_samples_kept"))      chum_samples_kept      <- readRDS("chum_posterior_samples.rds")

chilcotin_post <- as.data.frame(as.matrix(chilcotin_samples_kept[, c("a", "b", "s", "t")]))
thompson_post  <- as.data.frame(as.matrix(thompson_samples_kept[,  c("a", "b", "s", "t")]))
chum_post      <- as.data.frame(as.matrix(chum_samples_kept[,      c("a", "b", "s", "t")]))

# ------------------------------------------------------------
# 2. SUBSAMPLE N_DRAWS FROM EACH POSTERIOR
#    (fit independently -- draw i of each is paired purely to keep the
#    loop to N_DRAWS total runs rather than N_DRAWS^3)
# ------------------------------------------------------------

draw_idx_chilcotin <- sample(seq_len(nrow(chilcotin_post)), N_DRAWS, replace = nrow(chilcotin_post) < N_DRAWS)
draw_idx_thompson  <- sample(seq_len(nrow(thompson_post)),  N_DRAWS, replace = nrow(thompson_post)  < N_DRAWS)
draw_idx_chum      <- sample(seq_len(nrow(chum_post)),      N_DRAWS, replace = nrow(chum_post)      < N_DRAWS)

chilcotin_draws <- chilcotin_post[draw_idx_chilcotin, ]
thompson_draws  <- thompson_post[draw_idx_thompson, ]
chum_draws      <- chum_post[draw_idx_chum, ]

# ------------------------------------------------------------
# 3. CHUM'S OWN COVARIATES -- from chum_data.csv directly (it already
#    has both SL and PDO), NOT the separate covariates.csv the original
#    run_chum_scenario() needed. That file was only necessary for the
#    old covariate set (PDO_adult, PDO_smolt, SSL/SL) -- now that chum's
#    model is just SL + a single PDO, and chum_data.csv already has
#    exactly those two columns, there's no reason to bring in a second
#    external file at all.
# ------------------------------------------------------------

chum_data_raw <- read_csv("data/chum_data.csv")
chum_covariates <- chum_data_raw %>% select(Year, SL, PDO)

# ------------------------------------------------------------
# 4. run_chum_scenario_v2() -- Spawners + SL + single PDO, coefficients
#    as arguments (posterior draws), not hardcoded. Age-structured
#    recursion, NA-residual fallback, and SSL-freeze mechanism are
#    otherwise UNCHANGED from your cleaned run_chum_scenario() --
#    none of that logic depends on which/how-many covariates feed
#    chum_base_alpha.
# ------------------------------------------------------------

run_chum_scenario_v2 <- function(data, chum_covariates, U_apply, SSL_control,
                                 U_historic = 0, start_year = 1978,
                                 chum_intercept, chum_pdo_coef, chum_ssl_coef,
                                 chum_spawners_coef) {
  
  # SCALING: chum_spawners_coef (from the bounded-b JAGS fit) was fit on
  # spawners / 100000, matching the original WinBUGS convention -- same
  # reasoning as Thompson/Chilcotin's / 1000 scaling. The exponential
  # regression (chum_model_recruits, chum_recruits_alt) works entirely in
  # this scaled unit; raw-count quantities (Nage*_pred/alt, sum_alt,
  # catch_alt, chum_spawners_pred) convert back via *100000 at the point
  # each age-class is built from lagged recruits -- mirrors exactly where
  # Thompson's *1000 conversion happens (in Nage4_pred etc.), not baked
  # into recruits_alt itself.
  CHUM_SCALE <- 100000
  
  new.data <- data %>%
    left_join(chum_covariates, by = "Year") %>%
    arrange(Year) %>%
    mutate(
      catch = chum_total_stock - chum_spawners,
      U_chum = catch / chum_total_stock,
      chum_base_alpha = chum_intercept + (PDO * chum_pdo_coef + SL * chum_ssl_coef),
      chum_model_recruits =
        (chum_spawners / CHUM_SCALE) *
        exp(chum_base_alpha + (chum_spawners / CHUM_SCALE) * chum_spawners_coef),
      chum_ln_obs_pred = log((chum_recruits_obs / CHUM_SCALE) / chum_model_recruits),
      Nage3_obs = chum_total_stock * prop3,
      Nage4_obs = chum_total_stock * prop4,
      Nage5_obs = chum_total_stock * prop5,
      Nage6_obs = chum_total_stock * prop6,
      Nage3_pred = case_when(
        Year <= 1953 ~ Nage3_obs,
        Year >= 1954 ~ lag(chum_model_recruits, 3) * prop3 * exp(lag(chum_ln_obs_pred, 3)) * CHUM_SCALE),
      Nage4_pred = case_when(
        Year <= 1954 ~ Nage4_obs,
        Year >= 1955 ~ lag(chum_model_recruits, 4) * prop4 * exp(lag(chum_ln_obs_pred, 4)) * CHUM_SCALE),
      Nage5_pred = case_when(
        Year <= 1955 ~ Nage5_obs,
        Year >= 1956 ~ lag(chum_model_recruits, 5) * prop5 * exp(lag(chum_ln_obs_pred, 5)) * CHUM_SCALE),
      Nage6_pred = case_when(
        Year <= 1956 ~ Nage6_obs,
        Year >= 1957 ~ lag(chum_model_recruits, 6) * prop6 * exp(lag(chum_ln_obs_pred, 6)) * CHUM_SCALE),
      recruits_pred = rowSums(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)), na.rm = TRUE),
      U_chum_pred = catch / recruits_pred,
      chum_commercial_harvest_uapply = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 & U_historic == 0 ~ U_apply,
        Year >= 1991 & U_historic == 1 ~ U_chum_pred)
    )
  
  SSL_1978 <- new.data %>% filter(Year == 1978) %>% pull(SL)
  
  df <- new.data %>%
    arrange(Year) %>%
    mutate(
      chum_recruits_alt = NA_real_,   # stays in SCALED (/100000) units throughout
      Nage3_alt = Nage3_pred, Nage4_alt = Nage4_pred,
      Nage5_alt = Nage5_pred, Nage6_alt = Nage6_pred,
      sum_alt = NA_real_,
      catch_alt = NA_real_,
      chum_spawners_pred = chum_spawners,   # RAW units
      chum_SSL_alt = NA_real_,
      chum_SSL_alpha = NA_real_
    )
  
  for (i in seq_len(nrow(df))) {
    
    df$chum_SSL_alt[i] <- if (df$Year[i] <= 1978) {
      df$SL[i]
    } else {
      (1 - SSL_control) * df$SL[i] + SSL_control * SSL_1978
    }
    
    df$chum_SSL_alpha[i] <- chum_intercept + df$PDO[i] * chum_pdo_coef + df$chum_SSL_alt[i] * chum_ssl_coef
    
    # Nage*_alt built from LAGGED chum_recruits_alt (scaled units) --
    # *CHUM_SCALE converts back to raw here, same point Thompson does it
    if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i] * CHUM_SCALE
    if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i] * CHUM_SCALE
    if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i] * CHUM_SCALE
    if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i] * CHUM_SCALE
    
    df$sum_alt[i] <- sum(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i], na.rm = TRUE)
    
    # RAW-unit bookkeeping: catch, spawners_pred -- unaffected by the
    # scaling fix, same as before
    df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
    df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i])
    
    resid_adj <- if (is.na(df$chum_ln_obs_pred[i])) 1 else exp(df$chum_ln_obs_pred[i])
    
    # spk: spawners_pred converted to the SAME scaled unit the JAGS b was
    # fit on -- this is the actual fix. recruits_alt stays in scaled
    # units (converted back to raw only when building the NEXT age-class
    # above, via *CHUM_SCALE)
    spk <- df$chum_spawners_pred[i] / CHUM_SCALE
    
    df$chum_recruits_alt[i] <-
      if (SSL_control == 0) {
        spk * exp(df$chum_base_alpha[i] + chum_spawners_coef * spk) * resid_adj
      } else {
        spk * exp(df$chum_SSL_alpha[i] + chum_spawners_coef * spk) * resid_adj
      }
  }
  
  df %>% mutate(ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control"))
}

# ------------------------------------------------------------
# 5. run_thompson_scenario_v2() -- unchanged from your existing version
#    (Spawners + SL + SST, already matches the updated dredge results)
# ------------------------------------------------------------

run_thompson_scenario_v2 <- function(all_data, SSL_control, U_historic, byrate,
                                     start_year = 1978,
                                     sh_thompson_intercept,
                                     sh_thompson_sst_coef,
                                     sh_thompson_ssl_coef,
                                     sh_thompson_spawners_coef) {
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  FN_thompson_2018 <- all_data %>% filter(Year == 2018) %>% pull(sh_thompson_FN_mortalities) %>% as.numeric()
  
  df <- all_data %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_thompson_base_alpha = sh_thompson_intercept +
        sh_thompson_SST * sh_thompson_sst_coef +
        sh_thompson_SL  * sh_thompson_ssl_coef,
      sh_thompson_model_recruits =
        (sh_thompson_spawners / 1000) *
        exp(sh_thompson_base_alpha + (sh_thompson_spawners / 1000) * sh_thompson_spawners_coef),
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
      if (yr <= 1990) df$sh_thompson_U[i]
    else if (U_historic == 1) byrate * df$sh_thompson_U[i]
    else byrate * df$chum_commercial_harvest_uapply[i]
    
    df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) df$sh_thompson_SL[i]
    else (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
    
    df$sh_thompson_SSL_alpha[i] <-
      sh_thompson_intercept + df$sh_thompson_SST[i] * sh_thompson_sst_coef + df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    df$sh_thompson_alpha_CN[i] <-
      sh_thompson_intercept + df$sh_thompson_SST[i] * sh_thompson_sst_coef + df$sh_thompson_SL[i] * sh_thompson_ssl_coef
    
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_thompson_ln_obs_pred[i])) 1 else exp(df$sh_thompson_ln_obs_pred[i])
      
      df$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) spk * exp(df$sh_thompson_alpha_CN[i] + sh_thompson_spawners_coef * spk) * resid_adj
      else spk * exp(df$sh_thompson_SSL_alpha[i] + sh_thompson_spawners_coef * spk) * resid_adj
      
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
        df$sh_thompson_FN_catch_pred[i] <- FN_thompson_2018 / denom_2018 * (df$sh_thompson_sum_pred[i] - df$sh_thompson_bycatch_pred[i])
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

# ------------------------------------------------------------
# 6. run_chilcotin_scenario_v2() -- max_flow REMOVED. Now Spawners +
#    SL + SST only, structurally identical to Thompson, matching the
#    updated dredge results.
# ------------------------------------------------------------

run_chilcotin_scenario_v2 <- function(df, SSL_control, U_historic, byrate,
                                      start_year = 1973,
                                      sh_chilcotin_intercept,
                                      sh_chilcotin_sst_coef,
                                      sh_chilcotin_ssl_coef,
                                      sh_chilcotin_spawners_coef) {
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  FN_chilcotin_2018 <- df %>% filter(Year == 2018) %>% pull(sh_chilcotin_FN_mortalities) %>% as.numeric()
  
  df <- df %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_chilcotin_base_alpha = sh_chilcotin_intercept +
        sh_chilcotin_SST * sh_chilcotin_sst_coef +
        sh_chilcotin_SL  * sh_chilcotin_ssl_coef,
      sh_chilcotin_model_recruits =
        (sh_chilcotin_spawners / 1000) *
        exp(sh_chilcotin_base_alpha + (sh_chilcotin_spawners / 1000) * sh_chilcotin_spawners_coef),
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
      if (yr <= 1990) df$sh_chilcotin_U[i]
    else if (U_historic == 1) byrate * df$sh_chilcotin_U[i]
    else byrate * df$chum_commercial_harvest_uapply[i]
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) df$sh_chilcotin_SL[i]
    else (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
    
    df$sh_chilcotin_SSL_alpha[i] <-
      sh_chilcotin_intercept + df$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef + df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    df$sh_chilcotin_alpha_CN[i] <-
      sh_chilcotin_intercept + df$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef + df$sh_chilcotin_SL[i] * sh_chilcotin_ssl_coef
    
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_chilcotin_ln_obs_pred[i])) 1 else exp(df$sh_chilcotin_ln_obs_pred[i])
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) spk * exp(df$sh_chilcotin_alpha_CN[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
      else spk * exp(df$sh_chilcotin_SSL_alpha[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
      
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
        df$sh_chilcotin_FN_catch_pred[i] <- FN_chilcotin_2018 / denom_2018 * (df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_bycatch_pred[i])
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

# ------------------------------------------------------------
# 7. EACH STOCK'S OWN 10-YEAR LOW PERIOD -- real observed abundance,
#    lowest mean. Fixed, computed once (doesn't vary by draw).
# ------------------------------------------------------------

chum_observed_abundance <- data %>%
  filter(is.finite(chum_recruits_obs)) %>%
  arrange(Year) %>%
  mutate(roll_mean = zoo::rollapply(chum_recruits_obs, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA))

thompson_observed_abundance <- sh_data %>%
  filter(is.finite(sh_thompson_prefishery_N)) %>%
  arrange(Year) %>%
  mutate(roll_mean = zoo::rollapply(sh_thompson_prefishery_N, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA))

chilcotin_observed_abundance <- sh_data %>%
  filter(is.finite(sh_chilcotin_prefishery_N)) %>%
  arrange(Year) %>%
  mutate(roll_mean = zoo::rollapply(sh_chilcotin_prefishery_N, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA))

get_low_period <- function(df) {
  df %>%
    filter(!is.na(roll_mean)) %>%
    slice_min(roll_mean, n = 1, with_ties = FALSE) %>%
    transmute(low_period_start = Year - LOW_PERIOD_YEARS + 1, low_period_end = Year, low_period_mean = roll_mean)
}

low_periods <- bind_rows(
  get_low_period(chum_observed_abundance)      %>% mutate(Stock = "Chum"),
  get_low_period(thompson_observed_abundance)  %>% mutate(Stock = "Thompson steelhead"),
  get_low_period(chilcotin_observed_abundance) %>% mutate(Stock = "Chilcotin steelhead")
)
print(low_periods)

# Diagnostic plot: rolling 10-year mean abundance per stock, with the
# selected low-period window shaded and its mean marked -- same style
# as the sockeye/coho low-period diagnostics
observed_all <- bind_rows(
  chum_observed_abundance      %>% transmute(Year, roll_mean, Stock = "Chum"),
  thompson_observed_abundance  %>% transmute(Year, roll_mean, Stock = "Thompson steelhead"),
  chilcotin_observed_abundance %>% transmute(Year, roll_mean, Stock = "Chilcotin steelhead")
)

ggplot(observed_all, aes(Year, roll_mean)) +
  geom_line(color = "#4682B4") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.15) +
  geom_point(data = low_periods, aes(x = low_period_end, y = low_period_mean), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "10-year rolling mean abundance") +
  theme_minimal()

ggsave("figures/steelhead_chum_low_period_diagnostic.png", width = 10, height = 4, dpi = 600)

natural_data_end <- list(
  Chum = max(chum_observed_abundance$Year),
  `Thompson steelhead` = max(thompson_observed_abundance$Year),
  `Chilcotin steelhead` = max(chilcotin_observed_abundance$Year)
)

# ------------------------------------------------------------
# 8. MONTE CARLO LOOP -- historic harvest scenario, SSL control vs. no
#    control, ALL THREE STOCKS drawn per iteration
# ------------------------------------------------------------

lowpoint_results <- vector("list", N_DRAWS * 2)
counter <- 0

for (d in seq_len(N_DRAWS)) {
  
  ch <- chilcotin_draws[d, ]
  th <- thompson_draws[d, ]
  cm <- chum_draws[d, ]
  
  for (ssl in c(0, 1)) {
    
    chum_df <- run_chum_scenario_v2(
      data, chum_covariates, U_apply = 0, SSL_control = ssl, U_historic = 1,
      chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
    )
    
    all_data <- chum_df %>% right_join(sh_data, by = "Year")
    
    thompson_df <- run_thompson_scenario_v2(
      all_data, SSL_control = ssl, U_historic = 1, byrate = default_byrate,
      sh_thompson_intercept = th$a, sh_thompson_sst_coef = th$t,
      sh_thompson_ssl_coef = th$s, sh_thompson_spawners_coef = -th$b
    )
    
    chilcotin_df <- run_chilcotin_scenario_v2(
      thompson_df, SSL_control = ssl, U_historic = 1, byrate = default_byrate,
      sh_chilcotin_intercept = ch$a, sh_chilcotin_sst_coef = ch$t,
      sh_chilcotin_ssl_coef = ch$s, sh_chilcotin_spawners_coef = -ch$b
    )
    
    ssl_label <- if (ssl == 1) "SSL control" else "No SSL control"
    
    for (stock_name in c("Chum", "Thompson steelhead", "Chilcotin steelhead")) {
      
      abundance_col <- switch(stock_name,
                              "Chum" = "sum_alt",
                              "Thompson steelhead" = "sh_thompson_sum_pred",
                              "Chilcotin steelhead" = "sh_chilcotin_sum_pred"
      )
      
      final_end <- natural_data_end[[stock_name]]
      final_start <- final_end - FINAL_PERIOD_YEARS + 1
      low_p <- low_periods %>% filter(Stock == stock_name)
      
      X <- chilcotin_df %>%
        filter(Year >= final_start, Year <= final_end) %>%
        summarise(m = mean(.data[[abundance_col]], na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = stock_name, ssl_scenario = ssl_label, draw = d,
        Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean
      )
    }
  }
  
  if (d %% 10 == 0) message("draw ", d, " / ", N_DRAWS)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "steelhead_chum_ratio_summary.rds")

# Diagnostic plot: distribution of X/Xmin draws per stock, colored by
# SSL scenario -- same style as the coho/sockeye ratio histograms
SSL_COLORS <- c("No SSL control" = "black", "SSL control" = "#4682B4")

ggplot(lowpoint_summary, aes(ratio, fill = ssl_scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SSL_COLORS, name = NULL) +
  facet_wrap(~ Stock, scales = "free") +
  labs(x = "X / Xmin", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_histogram.png", width = 10, height = 4, dpi = 600)

# ------------------------------------------------------------
# 10. 40% Smsy BENCHMARK -- per-year actual covariates, averaged output,
#     same approach as coho/sockeye. Reference period defaults to each
#     stock's FULL historical data range (no established high-
#     productivity reference year exists for steelhead/chum the way
#     1970 did for pinniped scenarios) -- say if you want a specific
#     window instead. Reuses the SAME draws already sampled for the
#     X/Xmin loop, for consistency.
#
# SCALING: Smsy/recruits come out in the SAME scaled unit each stock's
# b was fit on (Thompson/Chilcotin /1000, Chum /100000) -- converted
# back to raw via *SCALE before dividing by Xmin (which is raw).
# ------------------------------------------------------------

compute_smsy_benchmark <- function(draws_df, covariate_data, covariate_cols, coef_cols, Xmin, scale) {
  n_years <- nrow(covariate_data)
  a_mat <- outer(draws_df$a, rep(1, n_years))
  for (k in seq_along(covariate_cols)) {
    a_mat <- a_mat + outer(draws_df[[coef_cols[k]]], covariate_data[[covariate_cols[k]]])
  }
  beta_draw <- draws_df$b   # already the positive Ricker slope, bounded model -- no sign flip
  
  Smsy_mat <- (a_mat / beta_draw) * (0.5 - 0.07 * a_mat)
  R_mat    <- (0.4 * Smsy_mat) * exp(a_mat - beta_draw * (0.4 * Smsy_mat)) * scale
  
  R_draw_avg <- rowMeans(R_mat, na.rm = TRUE)
  ratio_draw <- R_draw_avg / Xmin
  tibble(benchmark_ratio = median(ratio_draw, na.rm = TRUE), n_years_used = n_years)
}

chum_covariates_complete <- chum_covariates %>% filter(!is.na(SL), !is.na(PDO))

smsy_benchmarks <- bind_rows(
  compute_smsy_benchmark(
    chilcotin_draws, chilcotin_data, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Chilcotin steelhead"], scale = 1000
  ) %>% mutate(Stock = "Chilcotin steelhead"),
  compute_smsy_benchmark(
    thompson_draws, thompson_data, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Thompson steelhead"], scale = 1000
  ) %>% mutate(Stock = "Thompson steelhead"),
  compute_smsy_benchmark(
    chum_draws, chum_covariates_complete, c("SL", "PDO"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Chum"], scale = 100000
  ) %>% mutate(Stock = "Chum")
)
print(smsy_benchmarks, n = Inf)

# ------------------------------------------------------------
# 11. SUMMARIZE + PLOT -- X/Xmin, pointrange, natural log scale, plus
#     the 40% Smsy benchmark, faceted by Chum vs. Steelhead (abundance
#     scales differ by orders of magnitude between the two)
# ------------------------------------------------------------

mc_summary_stats_ratio <- lowpoint_summary %>%
  group_by(Stock, ssl_scenario) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.05, na.rm = TRUE),
    q95    = quantile(ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species_type = if_else(Stock == "Chum", "Chum", "Steelhead"))
print(mc_summary_stats_ratio, n = Inf)

smsy_benchmarks <- smsy_benchmarks %>%
  mutate(species_type = if_else(Stock == "Chum", "Chum", "Steelhead"))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = ssl_scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.4), size = 0.6, linewidth = 1) +
  geom_point(data = smsy_benchmarks, aes(x = Stock, y = benchmark_ratio, shape = "40% Smsy"),
             inherit.aes = FALSE, color = "firebrick", size = 6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_wrap(~ species_type, scales = "free", ncol = 1) +
  scale_color_manual(values = SSL_COLORS, name = NULL) +
  scale_shape_manual(name = NULL, values = c("40% Smsy" = 95)) +   # pch 95 = horizontal dash
  scale_y_continuous(trans = "log", labels = scales::comma) +
  labs(x = NULL, y = "Recruits / min Recruits") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio.png", width = 9, height = 8, dpi = 600)
