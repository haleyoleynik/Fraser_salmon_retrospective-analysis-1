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

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50   # start small, matching the pattern used throughout this pipeline
LOW_PERIOD_YEARS <- 10
MIN_YEARS_BEFORE_LOW_PERIOD <- 10   # same buffer rule as sockeye -- a candidate low
# window is only considered if it's preceded by a
# full 10 years of real data, preventing selection
# right at the very start of a stock's record
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
# 3b. SHARED GAP-FILL HELPER -- same convention used throughout the
#     sockeye pipeline: any TRAILING gap in a covariate (including SL/
#     SST for years beyond real data -- confirmed Chilcotin's sh_data
#     has NO pre-filled future years at all, unlike Thompson's, which
#     happens to; relying on that inconsistency would be fragile) gets
#     filled with the mean of the 10 years immediately before the gap
#     starts. Applied to SL/SST inside both steelhead functions below,
#     BEFORE alpha_scenario is computed -- this is what was silently
#     going NA for Chilcotin and then getting masked by na.rm=TRUE in
#     the sum_pred calculation, eventually producing a hard 0 once every
#     age class's lag reference fell into the gap simultaneously.
# ------------------------------------------------------------

fill_trailing_gap <- function(vals) {
  n <- length(vals)
  if (all(is.na(vals))) return(vals)
  last_valid <- max(which(!is.na(vals)))
  if (last_valid == n) return(vals)
  window_start <- max(1, last_valid - 10 + 1)
  fill_val <- mean(vals[window_start:last_valid], na.rm = TRUE)
  vals[(last_valid + 1):n] <- fill_val
  vals
}

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

run_thompson_scenario_v2 <- function(all_data, SSL_control, SST_control, U_historic, byrate,
                                     start_year = 1978,
                                     sh_thompson_intercept,
                                     sh_thompson_sst_coef,
                                     sh_thompson_ssl_coef,
                                     sh_thompson_spawners_coef,
                                     sst_baseline_thompson,
                                     zero_harvest = FALSE) {
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  
  # FIX: fill any trailing gap in SL/SST (confirmed Chilcotin's sh_data
  # has none for years beyond real data; Thompson's happens to already
  # have repeated placeholder values, which is a fragile coincidence,
  # not something to rely on -- apply the same principled fill to both).
  all_data <- all_data %>%
    arrange(Year) %>%
    mutate(
      sh_thompson_SL  = fill_trailing_gap(sh_thompson_SL),
      sh_thompson_SST = fill_trailing_gap(sh_thompson_SST)
    )
  
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
      # FIX: FN_mortalities and sport_mortalities were previously used as
      # FIXED ABSOLUTE historical counts, applied directly regardless of
      # what the model itself reconstructed for sum_pred. Any posterior
      # draw with even slightly weaker productivity than the original
      # point estimate could then be asked to "pay" the same absolute
      # catch as the real, larger historical population -- overdrawing
      # the model's own smaller reconstruction and forcing spawners to
      # exactly 0, which then propagates forward permanently (0 in,
      # 0 out, forever). Converting to RATES (share of real historical
      # abundance) and applying them PROPORTIONALLY to each draw's own
      # sum_pred -- same principle already used correctly for bycatch_pred
      # via U_comm -- makes total_catch_pred structurally unable to
      # exceed sum_pred.
      sh_thompson_FN_rate    = sh_thompson_FN_mortalities / sh_thompson_prefishery_N,
      sh_thompson_sport_rate = sh_thompson_sport_mortalities / sh_thompson_prefishery_N,
      sh_thompson_recruits_alt      = NA_real_,
      sh_thompson_spawners_pred     = sh_thompson_spawners,
      sh_thompson_Nage4_pred = NA_real_, sh_thompson_Nage5_pred = NA_real_,
      sh_thompson_Nage6_pred = NA_real_, sh_thompson_Nage7_pred = NA_real_,
      sh_thompson_Nage8_pred = NA_real_,
      sh_thompson_SSL_alt           = NA_real_,
      sh_thompson_SST_alt           = NA_real_,
      sh_thompson_alpha_scenario    = NA_real_,
      sh_thompson_sum_pred          = NA_real_,
      sh_thompson_bycatch_pred      = NA_real_,
      sh_thompson_FN_catch_pred     = NA_real_,
      sh_thompson_sport_catch_pred  = NA_real_,
      sh_thompson_total_catch_pred  = NA_real_,
      sh_thompson_U_comm            = NA_real_
    )
  
  # Fallback rates -- mean of this stock's own last 10 REAL years' own
  # FN/sport harvest rate, same convention as the U_comm fallback
  real_FN_rate    <- df$sh_thompson_FN_rate[!is.na(df$sh_thompson_FN_rate)]
  real_sport_rate <- df$sh_thompson_sport_rate[!is.na(df$sh_thompson_sport_rate)]
  FN_rate_fallback_thompson    <- mean(tail(real_FN_rate, 10), na.rm = TRUE)
  sport_rate_fallback_thompson <- mean(tail(real_sport_rate, 10), na.rm = TRUE)
  
  # FIX: sh_thompson_U (real harvest rate) is only defined where real
  # prefishery_N/sport_mortalities/FN_mortalities exist -- for synthetic
  # (projected) years, these are all NA, so U_comm was silently going NA
  # too, cascading into permanent NA for every subsequent year. Fallback:
  # mean of this stock's own last 10 REAL years' harvest rate, same
  # "final-N-real-years mean" convention used everywhere else in this
  # pipeline (mirrors sockeye's retroU_default, just data-derived instead
  # of a fixed external constant).
  real_U_thompson <- df$sh_thompson_U[!is.na(df$sh_thompson_U)]
  retroU_fallback_thompson <- mean(tail(real_U_thompson, 10), na.rm = TRUE)
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_thompson_U_comm[i] <-
      if (zero_harvest) {
        0
      } else if (is.na(df$sh_thompson_U[i])) {
        byrate * retroU_fallback_thompson   # synthetic year -- no real harvest data
      } else if (yr <= 1990) {
        df$sh_thompson_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_thompson_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) df$sh_thompson_SL[i]
    else (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
    
    # SST freeze: baseline is the MEAN of this stock's own first 10 real
    # years (mirroring sockeye's period-mean baseline for SST, rather
    # than a single frozen year like SL uses) -- passed in as
    # sst_baseline_thompson, computed once outside this function.
    df$sh_thompson_SST_alt[i] <- (1 - SST_control) * df$sh_thompson_SST[i] + SST_control * sst_baseline_thompson
    
    # Unified alpha: independently reflects whichever of SL/SST is
    # frozen, generalizing the old alpha_CN/SSL_alpha binary branch to
    # any combination of the two controls (only "one on, one off" is
    # actually used below, but this handles all four combinations
    # correctly if ever needed).
    df$sh_thompson_alpha_scenario[i] <-
      sh_thompson_intercept + df$sh_thompson_SST_alt[i] * sh_thompson_sst_coef + df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_thompson_ln_obs_pred[i])) 1 else exp(df$sh_thompson_ln_obs_pred[i])
      
      df$sh_thompson_recruits_alt[i] <- spk * exp(df$sh_thompson_alpha_scenario[i] + sh_thompson_spawners_coef * spk) * resid_adj
      
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
      
      FN_rate_used_thompson    <- if (zero_harvest) 0 else if (is.na(df$sh_thompson_FN_rate[i]))    FN_rate_fallback_thompson    else df$sh_thompson_FN_rate[i]
      sport_rate_used_thompson <- if (zero_harvest) 0 else if (is.na(df$sh_thompson_sport_rate[i])) sport_rate_fallback_thompson else df$sh_thompson_sport_rate[i]
      
      df$sh_thompson_FN_catch_pred[i]    <- FN_rate_used_thompson    * df$sh_thompson_sum_pred[i]
      df$sh_thompson_sport_catch_pred[i] <- sport_rate_used_thompson * df$sh_thompson_sum_pred[i]
      
      df$sh_thompson_total_catch_pred[i] <- df$sh_thompson_FN_catch_pred[i] +
        df$sh_thompson_sport_catch_pred[i] + df$sh_thompson_bycatch_pred[i]
      
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

run_chilcotin_scenario_v2 <- function(df, SSL_control, SST_control, U_historic, byrate,
                                      start_year = 1973,
                                      sh_chilcotin_intercept,
                                      sh_chilcotin_sst_coef,
                                      sh_chilcotin_ssl_coef,
                                      sh_chilcotin_spawners_coef,
                                      sst_baseline_chilcotin,
                                      zero_harvest = FALSE) {
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  
  # FIX: fill any trailing gap in SL/SST -- this is the actual root
  # cause found for Chilcotin: sh_chilcotin_SL/SST were genuinely NA for
  # every year beyond 2017 (unlike Thompson's, which happened to already
  # have repeated placeholder values), making alpha_scenario NA there.
  # That NA didn't surface as NA in sum_pred -- na.rm=TRUE in the
  # rowSums-equivalent silently DROPPED each NA age-class contribution
  # as its lag reference entered the gap one at a time, understating the
  # population more each year, until all five age classes referenced the
  # gap simultaneously and sum(NA,NA,NA,NA,NA, na.rm=TRUE) returned an
  # honest-looking but wrong 0 (summing an empty set is defined as
  # zero) -- not a population collapse, a silently-masked missing-data
  # problem.
  df <- df %>%
    arrange(Year) %>%
    mutate(
      sh_chilcotin_SL  = fill_trailing_gap(sh_chilcotin_SL),
      sh_chilcotin_SST = fill_trailing_gap(sh_chilcotin_SST)
    )
  
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
      # Same FN/sport catch fix as Thompson -- rates instead of fixed
      # absolute counts, applied proportionally to this draw's own
      # sum_pred so total_catch_pred can never exceed it.
      sh_chilcotin_FN_rate    = sh_chilcotin_FN_mortalities / sh_chilcotin_prefishery_N,
      sh_chilcotin_sport_rate = sh_chilcotin_sport_mortalities / sh_chilcotin_prefishery_N,
      sh_chilcotin_recruits_alt      = NA_real_,
      sh_chilcotin_spawners_pred     = sh_chilcotin_spawners,
      sh_chilcotin_Nage4_pred = NA_real_, sh_chilcotin_Nage5_pred = NA_real_,
      sh_chilcotin_Nage6_pred = NA_real_, sh_chilcotin_Nage7_pred = NA_real_,
      sh_chilcotin_Nage8_pred = NA_real_,
      sh_chilcotin_SSL_alt           = NA_real_,
      sh_chilcotin_SST_alt           = NA_real_,
      sh_chilcotin_alpha_scenario    = NA_real_,
      sh_chilcotin_sum_pred          = NA_real_,
      sh_chilcotin_bycatch_pred      = NA_real_,
      sh_chilcotin_FN_catch_pred     = NA_real_,
      sh_chilcotin_sport_catch_pred  = NA_real_,
      sh_chilcotin_total_catch_pred  = NA_real_,
      sh_chilcotin_U_comm            = NA_real_
    )
  
  real_FN_rate_chilcotin    <- df$sh_chilcotin_FN_rate[!is.na(df$sh_chilcotin_FN_rate)]
  real_sport_rate_chilcotin <- df$sh_chilcotin_sport_rate[!is.na(df$sh_chilcotin_sport_rate)]
  FN_rate_fallback_chilcotin    <- mean(tail(real_FN_rate_chilcotin, 10), na.rm = TRUE)
  sport_rate_fallback_chilcotin <- mean(tail(real_sport_rate_chilcotin, 10), na.rm = TRUE)
  
  real_U_chilcotin <- df$sh_chilcotin_U[!is.na(df$sh_chilcotin_U)]
  retroU_fallback_chilcotin <- mean(tail(real_U_chilcotin, 10), na.rm = TRUE)
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_chilcotin_U_comm[i] <-
      if (zero_harvest) {
        0
      } else if (is.na(df$sh_chilcotin_U[i])) {
        byrate * retroU_fallback_chilcotin   # synthetic year -- no real harvest data
      } else if (yr <= 1990) {
        df$sh_chilcotin_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_chilcotin_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) df$sh_chilcotin_SL[i]
    else (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
    
    df$sh_chilcotin_SST_alt[i] <- (1 - SST_control) * df$sh_chilcotin_SST[i] + SST_control * sst_baseline_chilcotin
    
    df$sh_chilcotin_alpha_scenario[i] <-
      sh_chilcotin_intercept + df$sh_chilcotin_SST_alt[i] * sh_chilcotin_sst_coef + df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_chilcotin_ln_obs_pred[i])) 1 else exp(df$sh_chilcotin_ln_obs_pred[i])
      
      df$sh_chilcotin_recruits_alt[i] <- spk * exp(df$sh_chilcotin_alpha_scenario[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
      
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
      
      FN_rate_used_chilcotin    <- if (zero_harvest) 0 else if (is.na(df$sh_chilcotin_FN_rate[i]))    FN_rate_fallback_chilcotin    else df$sh_chilcotin_FN_rate[i]
      sport_rate_used_chilcotin <- if (zero_harvest) 0 else if (is.na(df$sh_chilcotin_sport_rate[i])) sport_rate_fallback_chilcotin else df$sh_chilcotin_sport_rate[i]
      
      df$sh_chilcotin_FN_catch_pred[i]    <- FN_rate_used_chilcotin    * df$sh_chilcotin_sum_pred[i]
      df$sh_chilcotin_sport_catch_pred[i] <- sport_rate_used_chilcotin * df$sh_chilcotin_sum_pred[i]
      
      df$sh_chilcotin_total_catch_pred[i] <- df$sh_chilcotin_FN_catch_pred[i] +
        df$sh_chilcotin_sport_catch_pred[i] + df$sh_chilcotin_bycatch_pred[i]
      
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
# 7. EACH STOCK'S OWN 10-YEAR LOW PERIOD -- now identified from real
#    observed ln(R/S) (productivity) instead of raw abundance, with the
#    same MIN_YEARS_BEFORE_LOW_PERIOD buffer rule used for sockeye
#    (prevents picking a window right at the very start of a stock's
#    record, before a rolling mean has real prior context). Xmin ITSELF
#    still comes from real observed ABUNDANCE over that window -- the
#    window is chosen via productivity, the value plugged into the
#    ratio stays abundance, same convention as the sockeye update.
#
#    ln(R/S) here is real recruits/real spawners directly (a ratio, so
#    it's scale-invariant regardless of any /1000 or /100000 scaling
#    used elsewhere in this script for the model fitting).
# ------------------------------------------------------------

chum_observed_lnrs <- data %>%
  filter(is.finite(chum_recruits_obs), is.finite(chum_spawners), chum_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(chum_recruits_obs / chum_spawners))

thompson_observed_lnrs <- sh_data %>%
  filter(is.finite(sh_thompson_recruits), is.finite(sh_thompson_spawners), sh_thompson_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(sh_thompson_recruits / sh_thompson_spawners))

chilcotin_observed_lnrs <- sh_data %>%
  filter(is.finite(sh_chilcotin_recruits), is.finite(sh_chilcotin_spawners), sh_chilcotin_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(sh_chilcotin_recruits / sh_chilcotin_spawners))

chum_observed_abundance <- data %>%
  filter(is.finite(chum_recruits_obs)) %>%
  select(Year, abundance = chum_recruits_obs)

thompson_observed_abundance <- sh_data %>%
  filter(is.finite(sh_thompson_prefishery_N)) %>%
  select(Year, abundance = sh_thompson_prefishery_N)

chilcotin_observed_abundance <- sh_data %>%
  filter(is.finite(sh_chilcotin_prefishery_N)) %>%
  select(Year, abundance = sh_chilcotin_prefishery_N)

get_low_period <- function(lnrs_df, abundance_df) {
  
  window <- lnrs_df %>%
    arrange(Year) %>%
    mutate(
      roll_mean = zoo::rollapply(lnRS, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA),
      roll_start = Year - LOW_PERIOD_YEARS + 1,
      min_year = min(Year)
    ) %>%
    filter(!is.na(roll_mean), roll_start >= min_year + MIN_YEARS_BEFORE_LOW_PERIOD) %>%
    slice_min(roll_mean, n = 1, with_ties = FALSE)
  
  low_period_mean <- abundance_df %>%
    filter(Year >= window$roll_start, Year <= window$Year) %>%
    summarise(m = mean(abundance, na.rm = TRUE)) %>% pull(m)
  
  tibble(low_period_start = window$roll_start, low_period_end = window$Year, low_period_mean = low_period_mean)
}

low_periods <- bind_rows(
  get_low_period(chum_observed_lnrs, chum_observed_abundance)           %>% mutate(Stock = "Chum"),
  get_low_period(thompson_observed_lnrs, thompson_observed_abundance)   %>% mutate(Stock = "Thompson steelhead"),
  get_low_period(chilcotin_observed_lnrs, chilcotin_observed_abundance) %>% mutate(Stock = "Chilcotin steelhead")
)
print(low_periods)

# Diagnostic plot: rolling 10-year mean ln(R/S) per stock, with the
# selected low-period window shaded -- matches the sockeye ln(R/S)
# diagnostic. Red dot marks the stock's own rolling lnRS at the window's
# end (not the abundance-based low_period_mean, since the axis here is
# productivity, not abundance).
observed_lnrs_rolling <- bind_rows(
  chum_observed_lnrs      %>% mutate(Stock = "Chum"),
  thompson_observed_lnrs  %>% mutate(Stock = "Thompson steelhead"),
  chilcotin_observed_lnrs %>% mutate(Stock = "Chilcotin steelhead")
) %>%
  group_by(Stock) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(roll_mean_lnRS = zoo::rollapply(lnRS, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA)) %>%
  ungroup()

lnrs_marker <- observed_lnrs_rolling %>%
  inner_join(low_periods %>% select(Stock, low_period_end), by = c("Stock", "Year" = "low_period_end"))

ggplot(observed_lnrs_rolling, aes(Year, roll_mean_lnRS)) +
  geom_line(color = "#4682B4") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.15) +
  geom_point(data = lnrs_marker, aes(x = Year, y = roll_mean_lnRS), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "10-year rolling mean ln(R/S)") +
  theme_minimal()

ggsave("figures/steelhead_chum_lnrs_low_period_v3.png", width = 10, height = 4, dpi = 600)

natural_data_end <- list(
  Chum = max(chum_observed_abundance$Year),
  `Thompson steelhead` = max(thompson_observed_abundance$Year),
  `Chilcotin steelhead` = max(chilcotin_observed_abundance$Year)
)

# ------------------------------------------------------------
# 7b. SST BASELINE, per steelhead stock -- 1950-1975 mean, matching
#     sockeye's convention. steelhead_data.csv's SST column is already
#     STANDARDIZED (z-scored), and we don't know the exact mean/SD used
#     to build it -- so rather than guess, recover the transformation
#     directly via linear regression against the raw SST anomaly series
#     (np_temp), over whichever years overlap with each stock's data.
#     Standardization is an affine transform, so a linear fit recovers
#     it exactly: standardized = raw*slope + intercept, where
#     slope = 1/SD and intercept = -mean/SD. R^2 should come out very
#     close to 1 if this is genuinely the same transform -- if it
#     isn't, don't trust the resulting baseline.
# ------------------------------------------------------------

np_temp <- read_csv("Data/ersst_MJJ_gulf_of_alaska_anomaly.csv") %>%
  rename(Year = year, sst = sst_anom_mjj)

recover_sst_standardization <- function(stock_data, np_temp) {
  merged <- stock_data %>%
    select(Year, SST_standardized = SST) %>%
    inner_join(np_temp, by = "Year")
  
  fit <- lm(SST_standardized ~ sst, data = merged)
  r_squared <- summary(fit)$r.squared
  
  cat("Recovered SST standardization -- n overlapping years:", nrow(merged),
      " | R^2:", round(r_squared, 6), "\n")
  if (r_squared < 0.999) {
    warning("R^2 is below 0.999 -- the recovered transform may not exactly match the ",
            "original standardization. Check before trusting the SST baseline below.")
  }
  
  list(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]), r_squared = r_squared)
}

thompson_sst_transform  <- recover_sst_standardization(thompson_data, np_temp)
chilcotin_sst_transform <- recover_sst_standardization(chilcotin_data, np_temp)

apply_sst_transform <- function(raw_sst, transform) transform$intercept + transform$slope * raw_sst

SST_BASELINE_YEARS <- 1950:1975   # matches the sockeye convention

sst_baseline_thompson <- np_temp %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  mutate(sst_standardized = apply_sst_transform(sst, thompson_sst_transform)) %>%
  summarise(m = mean(sst_standardized, na.rm = TRUE)) %>% pull(m)

sst_baseline_chilcotin <- np_temp %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  mutate(sst_standardized = apply_sst_transform(sst, chilcotin_sst_transform)) %>%
  summarise(m = mean(sst_standardized, na.rm = TRUE)) %>% pull(m)

cat("SST baseline (1950-1975, standardized units) -- Thompson:", sst_baseline_thompson,
    " | Chilcotin:", sst_baseline_chilcotin, "\n")
cat("For comparison, each stock's own historical SST range:\n")
print(summary(thompson_data$SST))
print(summary(chilcotin_data$SST))

# ------------------------------------------------------------
# 8. MONTE CARLO LOOP
#
#    CHUM: unchanged, own 2-way SSL control / no control comparison.
#
#    STEELHEAD (Thompson, Chilcotin): now TWO independent scenarios,
#    matching the sockeye design --
#      "Pinniped scenario": SSL frozen at 1978 value, SST real
#      "SST scenario":       SST frozen at each stock's own first-10-
#                            year baseline, SL real
#    Each compared against the SAME Xmin (fixed real data), same as
#    sockeye's Pinniped/SST/Pink comparison. Chum is still run once per
#    draw to build all_data (needed for the join), but its SSL_control
#    choice has NO effect on steelhead's results here -- this
#    retrospective always uses U_historic = 1, so the
#    chum_commercial_harvest_uapply interaction term never actually
#    fires (that branch only triggers when U_historic == 0). Chum is
#    run with SSL_control = 0 for this purpose, arbitrarily, since it's
#    inconsequential to the steelhead scenarios.
# ------------------------------------------------------------

lowpoint_results <- vector("list", N_DRAWS * 12)
counter <- 0

for (d in seq_len(N_DRAWS)) {
  
  ch <- chilcotin_draws[d, ]
  th <- thompson_draws[d, ]
  cm <- chum_draws[d, ]
  
  # ---- Chum: own SSL control / no control, as before ----
  X_actual_harvest <- NA_real_   # captured explicitly below, not relied on as loop leftover
  
  for (ssl in c(0, 1)) {
    
    chum_df <- run_chum_scenario_v2(
      data, chum_covariates, U_apply = 0, SSL_control = ssl, U_historic = 1,
      chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
    )
    
    ssl_label <- if (ssl == 1) "SSL control" else "No SSL control"
    final_end <- natural_data_end[["Chum"]]
    final_start <- final_end - FINAL_PERIOD_YEARS + 1
    low_p <- low_periods %>% filter(Stock == "Chum")
    
    X <- chum_df %>%
      filter(Year >= final_start, Year <= final_end) %>%
      summarise(m = mean(sum_alt, na.rm = TRUE)) %>% pull(m)
    
    if (ssl == 0) X_actual_harvest <- X   # this IS "actual harvest" (real, unfrozen SL, real harvest)
    
    counter <- counter + 1
    lowpoint_results[[counter]] <- tibble(
      Stock = "Chum", scenario = ssl_label, draw = d,
      Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean
    )
  }
  
  # ---- Chum: zero commercial harvest (1991+) vs actual harvest ----
  # SSL held at 0 (real, no pinniped counterfactual) for BOTH runs, to
  # isolate the harvest effect specifically -- same "vary one thing,
  # hold the other at its real value" principle used for steelhead's
  # Pinniped-vs-SST comparison. "Actual harvest" is IDENTICAL to
  # "No SSL control" above (U_historic=1, SSL_control=0) -- captured
  # directly from that run, not recomputed, just relabeled for this
  # comparison's clarity.
  chum_df_zero_harvest <- run_chum_scenario_v2(
    data, chum_covariates, U_apply = 0, SSL_control = 0, U_historic = 0,
    chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
  )
  
  X_zero_harvest <- chum_df_zero_harvest %>%
    filter(Year >= final_start, Year <= final_end) %>%
    summarise(m = mean(sum_alt, na.rm = TRUE)) %>% pull(m)
  
  counter <- counter + 1
  lowpoint_results[[counter]] <- tibble(
    Stock = "Chum", scenario = "Zero harvest", draw = d,
    Xmin = low_p$low_period_mean, X = X_zero_harvest, ratio = X_zero_harvest / low_p$low_period_mean
  )
  
  counter <- counter + 1
  lowpoint_results[[counter]] <- tibble(
    Stock = "Chum", scenario = "Actual harvest", draw = d,
    Xmin = low_p$low_period_mean, X = X_actual_harvest, ratio = X_actual_harvest / low_p$low_period_mean
  )
  
  # ---- Steelhead: Pinniped scenario and SST scenario, independently ----
  chum_df_for_steelhead <- run_chum_scenario_v2(
    data, chum_covariates, U_apply = 0, SSL_control = 0, U_historic = 1,
    chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
  )
  all_data <- chum_df_for_steelhead %>% right_join(sh_data, by = "Year")
  
  # ---- Steelhead: zero harvest (bycatch + FN catch + sport catch, ALL
  # zeroed via the new zero_harvest flag) vs actual harvest, SL/SST held
  # real (unfrozen) in both, to isolate the harvest effect specifically
  # -- same "vary one thing, hold the rest real" principle as everywhere
  # else in this comparison. Chilcotin's function only reads its OWN
  # sh_chilcotin_* columns from its input df, not anything Thompson
  # computed, so a single baseline thompson_df is a valid chaining input
  # regardless of which harvest scenario Chilcotin itself is testing.
  thompson_df_baseline <- run_thompson_scenario_v2(
    all_data, SSL_control = 0, SST_control = 0, U_historic = 1, byrate = default_byrate,
    sh_thompson_intercept = th$a, sh_thompson_sst_coef = th$t,
    sh_thompson_ssl_coef = th$s, sh_thompson_spawners_coef = -th$b,
    sst_baseline_thompson = sst_baseline_thompson, zero_harvest = FALSE
  )
  thompson_zero <- run_thompson_scenario_v2(
    all_data, SSL_control = 0, SST_control = 0, U_historic = 1, byrate = default_byrate,
    sh_thompson_intercept = th$a, sh_thompson_sst_coef = th$t,
    sh_thompson_ssl_coef = th$s, sh_thompson_spawners_coef = -th$b,
    sst_baseline_thompson = sst_baseline_thompson, zero_harvest = TRUE
  )
  chilcotin_actual <- run_chilcotin_scenario_v2(
    thompson_df_baseline, SSL_control = 0, SST_control = 0, U_historic = 1, byrate = default_byrate,
    sh_chilcotin_intercept = ch$a, sh_chilcotin_sst_coef = ch$t,
    sh_chilcotin_ssl_coef = ch$s, sh_chilcotin_spawners_coef = -ch$b,
    sst_baseline_chilcotin = sst_baseline_chilcotin, zero_harvest = FALSE
  )
  chilcotin_zero <- run_chilcotin_scenario_v2(
    thompson_df_baseline, SSL_control = 0, SST_control = 0, U_historic = 1, byrate = default_byrate,
    sh_chilcotin_intercept = ch$a, sh_chilcotin_sst_coef = ch$t,
    sh_chilcotin_ssl_coef = ch$s, sh_chilcotin_spawners_coef = -ch$b,
    sst_baseline_chilcotin = sst_baseline_chilcotin, zero_harvest = TRUE
  )
  
  steelhead_harvest_runs <- list(
    "Thompson steelhead"  = list(actual = thompson_df_baseline, zero = thompson_zero, abundance_col = "sh_thompson_sum_pred"),
    "Chilcotin steelhead" = list(actual = chilcotin_actual,      zero = chilcotin_zero, abundance_col = "sh_chilcotin_sum_pred")
  )
  
  for (stock_name in names(steelhead_harvest_runs)) {
    
    run_set <- steelhead_harvest_runs[[stock_name]]
    final_end_sh <- natural_data_end[[stock_name]]
    final_start_sh <- final_end_sh - FINAL_PERIOD_YEARS + 1
    low_p_sh <- low_periods %>% filter(Stock == stock_name)
    
    X_actual_sh <- run_set$actual %>%
      filter(Year >= final_start_sh, Year <= final_end_sh) %>%
      summarise(m = mean(.data[[run_set$abundance_col]], na.rm = TRUE)) %>% pull(m)
    X_zero_sh <- run_set$zero %>%
      filter(Year >= final_start_sh, Year <= final_end_sh) %>%
      summarise(m = mean(.data[[run_set$abundance_col]], na.rm = TRUE)) %>% pull(m)
    
    counter <- counter + 1
    lowpoint_results[[counter]] <- tibble(
      Stock = stock_name, scenario = "Actual harvest", draw = d,
      Xmin = low_p_sh$low_period_mean, X = X_actual_sh, ratio = X_actual_sh / low_p_sh$low_period_mean
    )
    counter <- counter + 1
    lowpoint_results[[counter]] <- tibble(
      Stock = stock_name, scenario = "Zero harvest", draw = d,
      Xmin = low_p_sh$low_period_mean, X = X_zero_sh, ratio = X_zero_sh / low_p_sh$low_period_mean
    )
  }
  
  steelhead_scenarios <- list(
    "SSL control" = list(SSL_control = 1, SST_control = 0),
    "SST scenario" = list(SSL_control = 0, SST_control = 1)
  )
  
  for (scenario_label in names(steelhead_scenarios)) {
    
    sc <- steelhead_scenarios[[scenario_label]]
    
    thompson_df <- run_thompson_scenario_v2(
      all_data, SSL_control = sc$SSL_control, SST_control = sc$SST_control,
      U_historic = 1, byrate = default_byrate,
      sh_thompson_intercept = th$a, sh_thompson_sst_coef = th$t,
      sh_thompson_ssl_coef = th$s, sh_thompson_spawners_coef = -th$b,
      sst_baseline_thompson = sst_baseline_thompson
    )
    
    chilcotin_df <- run_chilcotin_scenario_v2(
      thompson_df, SSL_control = sc$SSL_control, SST_control = sc$SST_control,
      U_historic = 1, byrate = default_byrate,
      sh_chilcotin_intercept = ch$a, sh_chilcotin_sst_coef = ch$t,
      sh_chilcotin_ssl_coef = ch$s, sh_chilcotin_spawners_coef = -ch$b,
      sst_baseline_chilcotin = sst_baseline_chilcotin
    )
    
    for (stock_name in c("Thompson steelhead", "Chilcotin steelhead")) {
      
      abundance_col <- switch(stock_name,
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
        Stock = stock_name, scenario = scenario_label, draw = d,
        Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean
      )
    }
  }
  
  if (d %% 10 == 0) message("draw ", d, " / ", N_DRAWS)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "steelhead_chum_ratio_summary.rds")

# Diagnostic plot: distribution of X/Xmin draws per stock, colored by
# scenario -- same style as the coho/sockeye ratio histograms. Chum uses
# "SSL control"/"No SSL control"; steelhead uses "SSL control"/
# "SST scenario" -- now the SAME "SSL control" label/color for both,
# naturally separated by facet.
SCENARIO_COLORS <- c(
  "No SSL control" = "black", "SSL control" = "#4682B4", "SST scenario" = "#2E8B57"
)

ggplot(lowpoint_summary, aes(ratio, fill = scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS, name = NULL) +
  facet_wrap(~ Stock, scales = "free") +
  labs(x = "X / Xmin", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_histogram.png", width = 10, height = 4, dpi = 600)

ggplot(lowpoint_summary, aes(ratio, fill = scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "white") +
  #geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS, name = NULL) +
  facet_wrap(scenario~ Stock, scales = "free") +
  labs(x = "X / Xmin", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_histogram_harvest_1.png", width = 10, height = 4, dpi = 600)

# ------------------------------------------------------------
# 10. 40% Smsy BENCHMARK -- per-year actual covariates, averaged output,
#     same approach as coho/sockeye. Reference periods, per your
#     instruction:
#       Chum: 1960-1970
#       Thompson/Chilcotin: beginning of each stock's own time series
#                           through 1980
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

CHUM_SMSY_YEARS <- 1960:1970
STEELHEAD_SMSY_END_YEAR <- 1980

chilcotin_data_smsy_period <- chilcotin_data %>% filter(Year <= STEELHEAD_SMSY_END_YEAR)
thompson_data_smsy_period  <- thompson_data %>% filter(Year <= STEELHEAD_SMSY_END_YEAR)
chum_covariates_smsy_period <- chum_covariates %>%
  filter(!is.na(SL), !is.na(PDO), Year %in% CHUM_SMSY_YEARS)

cat("Smsy reference period year counts -- Chilcotin:", nrow(chilcotin_data_smsy_period),
    " | Thompson:", nrow(thompson_data_smsy_period),
    " | Chum:", nrow(chum_covariates_smsy_period), "\n")

smsy_benchmarks <- bind_rows(
  compute_smsy_benchmark(
    chilcotin_draws, chilcotin_data_smsy_period, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Chilcotin steelhead"], scale = 1000
  ) %>% mutate(Stock = "Chilcotin steelhead"),
  compute_smsy_benchmark(
    thompson_draws, thompson_data_smsy_period, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Thompson steelhead"], scale = 1000
  ) %>% mutate(Stock = "Thompson steelhead"),
  compute_smsy_benchmark(
    chum_draws, chum_covariates_smsy_period, c("SL", "PDO"), c("s", "t"),
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
  filter(!scenario %in% c("No SSL control", "Actual harvest", "Zero harvest")) %>%   # main plot: SSL/SST scenarios only
  group_by(Stock, scenario) %>%
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

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.4), size = 0.6, linewidth = 1) +
  geom_point(data = smsy_benchmarks, aes(x = Stock, y = benchmark_ratio, shape = "40% Smsy"),
             inherit.aes = FALSE, color = "firebrick", size = 6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_wrap(~ species_type, scales = "free", ncol = 1) +
  scale_color_manual(values = SCENARIO_COLORS, name = NULL) +
  scale_shape_manual(name = NULL, values = c("40% Smsy" = 95)) +   # pch 95 = horizontal dash
  scale_y_continuous(trans = "log", labels = scales::comma, breaks = scales::breaks_log(n = 6)) +
  labs(x = NULL, y = "Recruits / min Recruits") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_5000_draws.png", width = 9, height = 8, dpi = 600)

# ------------------------------------------------------------
# 11. HARVEST COMPARISON, ALL THREE STOCKS -- actual (real, reconstructed)
#     harvest vs. zero harvest, SL/SST held real (unfrozen) in both, to
#     isolate the harvest effect specifically. For chum, "zero harvest"
#     means zero commercial harvest from 1991 onward (real before that).
#     For Thompson/Chilcotin, it means bycatch + First Nations catch +
#     sport catch ALL zeroed throughout, via the zero_harvest flag --
#     not just bycatch, since FN/sport catch are independent of byrate.
#     Both scenarios use the SAME Xmin and the SAME (real-years-only)
#     final window as everything else in this script.
# ------------------------------------------------------------

harvest_comparison <- lowpoint_summary %>%
  filter(scenario %in% c("Actual harvest", "Zero harvest"))

harvest_comparison_summary <- harvest_comparison %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.10, na.rm = TRUE),
    q95    = quantile(ratio, 0.90, na.rm = TRUE),
    .groups = "drop"
  )
print(harvest_comparison_summary, n = Inf)

HARVEST_COLORS <- c("Actual harvest" = "black", "Zero harvest" = "#2E8B57")

ggplot(harvest_comparison_summary, aes(x = scenario, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), size = 0.8, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_wrap(~ Stock, scales = "free") +
  scale_color_manual(values = HARVEST_COLORS, name = NULL, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Recruits / min Recruits ") +
  theme_minimal()

ggsave("figures/harvest_comparison_all_stocks.png", width = 12, height = 5, dpi = 600)



######################################
#######################################
# with introduction of harvest x scenario 
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

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50   # start small, matching the pattern used throughout this pipeline
LOW_PERIOD_YEARS <- 10
MIN_YEARS_BEFORE_LOW_PERIOD <- 10   # same buffer rule as sockeye -- a candidate low
# window is only considered if it's preceded by a
# full 10 years of real data, preventing selection
# right at the very start of a stock's record
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
# 3b. SHARED GAP-FILL HELPER -- same convention used throughout the
#     sockeye pipeline: any TRAILING gap in a covariate (including SL/
#     SST for years beyond real data -- confirmed Chilcotin's sh_data
#     has NO pre-filled future years at all, unlike Thompson's, which
#     happens to; relying on that inconsistency would be fragile) gets
#     filled with the mean of the 10 years immediately before the gap
#     starts. Applied to SL/SST inside both steelhead functions below,
#     BEFORE alpha_scenario is computed -- this is what was silently
#     going NA for Chilcotin and then getting masked by na.rm=TRUE in
#     the sum_pred calculation, eventually producing a hard 0 once every
#     age class's lag reference fell into the gap simultaneously.
# ------------------------------------------------------------

fill_trailing_gap <- function(vals) {
  n <- length(vals)
  if (all(is.na(vals))) return(vals)
  last_valid <- max(which(!is.na(vals)))
  if (last_valid == n) return(vals)
  window_start <- max(1, last_valid - 10 + 1)
  fill_val <- mean(vals[window_start:last_valid], na.rm = TRUE)
  vals[(last_valid + 1):n] <- fill_val
  vals
}

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

run_thompson_scenario_v2 <- function(all_data, SSL_control, SST_control, U_historic, byrate,
                                     start_year = 1978,
                                     sh_thompson_intercept,
                                     sh_thompson_sst_coef,
                                     sh_thompson_ssl_coef,
                                     sh_thompson_spawners_coef,
                                     sst_baseline_thompson,
                                     zero_harvest = FALSE) {
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  
  # FIX: fill any trailing gap in SL/SST (confirmed Chilcotin's sh_data
  # has none for years beyond real data; Thompson's happens to already
  # have repeated placeholder values, which is a fragile coincidence,
  # not something to rely on -- apply the same principled fill to both).
  all_data <- all_data %>%
    arrange(Year) %>%
    mutate(
      sh_thompson_SL  = fill_trailing_gap(sh_thompson_SL),
      sh_thompson_SST = fill_trailing_gap(sh_thompson_SST)
    )
  
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
      # FIX: FN_mortalities and sport_mortalities were previously used as
      # FIXED ABSOLUTE historical counts, applied directly regardless of
      # what the model itself reconstructed for sum_pred. Any posterior
      # draw with even slightly weaker productivity than the original
      # point estimate could then be asked to "pay" the same absolute
      # catch as the real, larger historical population -- overdrawing
      # the model's own smaller reconstruction and forcing spawners to
      # exactly 0, which then propagates forward permanently (0 in,
      # 0 out, forever). Converting to RATES (share of real historical
      # abundance) and applying them PROPORTIONALLY to each draw's own
      # sum_pred -- same principle already used correctly for bycatch_pred
      # via U_comm -- makes total_catch_pred structurally unable to
      # exceed sum_pred.
      sh_thompson_FN_rate    = sh_thompson_FN_mortalities / sh_thompson_prefishery_N,
      sh_thompson_sport_rate = sh_thompson_sport_mortalities / sh_thompson_prefishery_N,
      sh_thompson_recruits_alt      = NA_real_,
      sh_thompson_spawners_pred     = sh_thompson_spawners,
      sh_thompson_Nage4_pred = NA_real_, sh_thompson_Nage5_pred = NA_real_,
      sh_thompson_Nage6_pred = NA_real_, sh_thompson_Nage7_pred = NA_real_,
      sh_thompson_Nage8_pred = NA_real_,
      sh_thompson_SSL_alt           = NA_real_,
      sh_thompson_SST_alt           = NA_real_,
      sh_thompson_alpha_scenario    = NA_real_,
      sh_thompson_sum_pred          = NA_real_,
      sh_thompson_bycatch_pred      = NA_real_,
      sh_thompson_FN_catch_pred     = NA_real_,
      sh_thompson_sport_catch_pred  = NA_real_,
      sh_thompson_total_catch_pred  = NA_real_,
      sh_thompson_U_comm            = NA_real_
    )
  
  # Fallback rates -- mean of this stock's own last 10 REAL years' own
  # FN/sport harvest rate, same convention as the U_comm fallback
  real_FN_rate    <- df$sh_thompson_FN_rate[!is.na(df$sh_thompson_FN_rate)]
  real_sport_rate <- df$sh_thompson_sport_rate[!is.na(df$sh_thompson_sport_rate)]
  FN_rate_fallback_thompson    <- mean(tail(real_FN_rate, 10), na.rm = TRUE)
  sport_rate_fallback_thompson <- mean(tail(real_sport_rate, 10), na.rm = TRUE)
  
  # FIX: sh_thompson_U (real harvest rate) is only defined where real
  # prefishery_N/sport_mortalities/FN_mortalities exist -- for synthetic
  # (projected) years, these are all NA, so U_comm was silently going NA
  # too, cascading into permanent NA for every subsequent year. Fallback:
  # mean of this stock's own last 10 REAL years' harvest rate, same
  # "final-N-real-years mean" convention used everywhere else in this
  # pipeline (mirrors sockeye's retroU_default, just data-derived instead
  # of a fixed external constant).
  real_U_thompson <- df$sh_thompson_U[!is.na(df$sh_thompson_U)]
  retroU_fallback_thompson <- mean(tail(real_U_thompson, 10), na.rm = TRUE)
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_thompson_U_comm[i] <-
      if (zero_harvest) {
        0
      } else if (is.na(df$sh_thompson_U[i])) {
        byrate * retroU_fallback_thompson   # synthetic year -- no real harvest data
      } else if (yr <= 1990) {
        df$sh_thompson_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_thompson_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) df$sh_thompson_SL[i]
    else (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
    
    # SST freeze: baseline is the MEAN of this stock's own first 10 real
    # years (mirroring sockeye's period-mean baseline for SST, rather
    # than a single frozen year like SL uses) -- passed in as
    # sst_baseline_thompson, computed once outside this function.
    df$sh_thompson_SST_alt[i] <- (1 - SST_control) * df$sh_thompson_SST[i] + SST_control * sst_baseline_thompson
    
    # Unified alpha: independently reflects whichever of SL/SST is
    # frozen, generalizing the old alpha_CN/SSL_alpha binary branch to
    # any combination of the two controls (only "one on, one off" is
    # actually used below, but this handles all four combinations
    # correctly if ever needed).
    df$sh_thompson_alpha_scenario[i] <-
      sh_thompson_intercept + df$sh_thompson_SST_alt[i] * sh_thompson_sst_coef + df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_thompson_ln_obs_pred[i])) 1 else exp(df$sh_thompson_ln_obs_pred[i])
      
      df$sh_thompson_recruits_alt[i] <- spk * exp(df$sh_thompson_alpha_scenario[i] + sh_thompson_spawners_coef * spk) * resid_adj
      
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
      
      FN_rate_used_thompson    <- if (zero_harvest) 0 else if (is.na(df$sh_thompson_FN_rate[i]))    FN_rate_fallback_thompson    else df$sh_thompson_FN_rate[i]
      sport_rate_used_thompson <- if (zero_harvest) 0 else if (is.na(df$sh_thompson_sport_rate[i])) sport_rate_fallback_thompson else df$sh_thompson_sport_rate[i]
      
      df$sh_thompson_FN_catch_pred[i]    <- FN_rate_used_thompson    * df$sh_thompson_sum_pred[i]
      df$sh_thompson_sport_catch_pred[i] <- sport_rate_used_thompson * df$sh_thompson_sum_pred[i]
      
      df$sh_thompson_total_catch_pred[i] <- df$sh_thompson_FN_catch_pred[i] +
        df$sh_thompson_sport_catch_pred[i] + df$sh_thompson_bycatch_pred[i]
      
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

run_chilcotin_scenario_v2 <- function(df, SSL_control, SST_control, U_historic, byrate,
                                      start_year = 1973,
                                      sh_chilcotin_intercept,
                                      sh_chilcotin_sst_coef,
                                      sh_chilcotin_ssl_coef,
                                      sh_chilcotin_spawners_coef,
                                      sst_baseline_chilcotin,
                                      zero_harvest = FALSE) {
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  
  # FIX: fill any trailing gap in SL/SST -- this is the actual root
  # cause found for Chilcotin: sh_chilcotin_SL/SST were genuinely NA for
  # every year beyond 2017 (unlike Thompson's, which happened to already
  # have repeated placeholder values), making alpha_scenario NA there.
  # That NA didn't surface as NA in sum_pred -- na.rm=TRUE in the
  # rowSums-equivalent silently DROPPED each NA age-class contribution
  # as its lag reference entered the gap one at a time, understating the
  # population more each year, until all five age classes referenced the
  # gap simultaneously and sum(NA,NA,NA,NA,NA, na.rm=TRUE) returned an
  # honest-looking but wrong 0 (summing an empty set is defined as
  # zero) -- not a population collapse, a silently-masked missing-data
  # problem.
  df <- df %>%
    arrange(Year) %>%
    mutate(
      sh_chilcotin_SL  = fill_trailing_gap(sh_chilcotin_SL),
      sh_chilcotin_SST = fill_trailing_gap(sh_chilcotin_SST)
    )
  
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
      # Same FN/sport catch fix as Thompson -- rates instead of fixed
      # absolute counts, applied proportionally to this draw's own
      # sum_pred so total_catch_pred can never exceed it.
      sh_chilcotin_FN_rate    = sh_chilcotin_FN_mortalities / sh_chilcotin_prefishery_N,
      sh_chilcotin_sport_rate = sh_chilcotin_sport_mortalities / sh_chilcotin_prefishery_N,
      sh_chilcotin_recruits_alt      = NA_real_,
      sh_chilcotin_spawners_pred     = sh_chilcotin_spawners,
      sh_chilcotin_Nage4_pred = NA_real_, sh_chilcotin_Nage5_pred = NA_real_,
      sh_chilcotin_Nage6_pred = NA_real_, sh_chilcotin_Nage7_pred = NA_real_,
      sh_chilcotin_Nage8_pred = NA_real_,
      sh_chilcotin_SSL_alt           = NA_real_,
      sh_chilcotin_SST_alt           = NA_real_,
      sh_chilcotin_alpha_scenario    = NA_real_,
      sh_chilcotin_sum_pred          = NA_real_,
      sh_chilcotin_bycatch_pred      = NA_real_,
      sh_chilcotin_FN_catch_pred     = NA_real_,
      sh_chilcotin_sport_catch_pred  = NA_real_,
      sh_chilcotin_total_catch_pred  = NA_real_,
      sh_chilcotin_U_comm            = NA_real_
    )
  
  real_FN_rate_chilcotin    <- df$sh_chilcotin_FN_rate[!is.na(df$sh_chilcotin_FN_rate)]
  real_sport_rate_chilcotin <- df$sh_chilcotin_sport_rate[!is.na(df$sh_chilcotin_sport_rate)]
  FN_rate_fallback_chilcotin    <- mean(tail(real_FN_rate_chilcotin, 10), na.rm = TRUE)
  sport_rate_fallback_chilcotin <- mean(tail(real_sport_rate_chilcotin, 10), na.rm = TRUE)
  
  real_U_chilcotin <- df$sh_chilcotin_U[!is.na(df$sh_chilcotin_U)]
  retroU_fallback_chilcotin <- mean(tail(real_U_chilcotin, 10), na.rm = TRUE)
  
  start_i <- which(df$Year >= start_year)[1]
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_chilcotin_U_comm[i] <-
      if (zero_harvest) {
        0
      } else if (is.na(df$sh_chilcotin_U[i])) {
        byrate * retroU_fallback_chilcotin   # synthetic year -- no real harvest data
      } else if (yr <= 1990) {
        df$sh_chilcotin_U[i]
      } else if (U_historic == 1) {
        byrate * df$sh_chilcotin_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) df$sh_chilcotin_SL[i]
    else (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
    
    df$sh_chilcotin_SST_alt[i] <- (1 - SST_control) * df$sh_chilcotin_SST[i] + SST_control * sst_baseline_chilcotin
    
    df$sh_chilcotin_alpha_scenario[i] <-
      sh_chilcotin_intercept + df$sh_chilcotin_SST_alt[i] * sh_chilcotin_sst_coef + df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_chilcotin_ln_obs_pred[i])) 1 else exp(df$sh_chilcotin_ln_obs_pred[i])
      
      df$sh_chilcotin_recruits_alt[i] <- spk * exp(df$sh_chilcotin_alpha_scenario[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
      
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
      
      FN_rate_used_chilcotin    <- if (zero_harvest) 0 else if (is.na(df$sh_chilcotin_FN_rate[i]))    FN_rate_fallback_chilcotin    else df$sh_chilcotin_FN_rate[i]
      sport_rate_used_chilcotin <- if (zero_harvest) 0 else if (is.na(df$sh_chilcotin_sport_rate[i])) sport_rate_fallback_chilcotin else df$sh_chilcotin_sport_rate[i]
      
      df$sh_chilcotin_FN_catch_pred[i]    <- FN_rate_used_chilcotin    * df$sh_chilcotin_sum_pred[i]
      df$sh_chilcotin_sport_catch_pred[i] <- sport_rate_used_chilcotin * df$sh_chilcotin_sum_pred[i]
      
      df$sh_chilcotin_total_catch_pred[i] <- df$sh_chilcotin_FN_catch_pred[i] +
        df$sh_chilcotin_sport_catch_pred[i] + df$sh_chilcotin_bycatch_pred[i]
      
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
# 7. EACH STOCK'S OWN 10-YEAR LOW PERIOD -- now identified from real
#    observed ln(R/S) (productivity) instead of raw abundance, with the
#    same MIN_YEARS_BEFORE_LOW_PERIOD buffer rule used for sockeye
#    (prevents picking a window right at the very start of a stock's
#    record, before a rolling mean has real prior context). Xmin ITSELF
#    still comes from real observed ABUNDANCE over that window -- the
#    window is chosen via productivity, the value plugged into the
#    ratio stays abundance, same convention as the sockeye update.
#
#    ln(R/S) here is real recruits/real spawners directly (a ratio, so
#    it's scale-invariant regardless of any /1000 or /100000 scaling
#    used elsewhere in this script for the model fitting).
# ------------------------------------------------------------

chum_observed_lnrs <- data %>%
  filter(is.finite(chum_recruits_obs), is.finite(chum_spawners), chum_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(chum_recruits_obs / chum_spawners))

thompson_observed_lnrs <- sh_data %>%
  filter(is.finite(sh_thompson_recruits), is.finite(sh_thompson_spawners), sh_thompson_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(sh_thompson_recruits / sh_thompson_spawners))

chilcotin_observed_lnrs <- sh_data %>%
  filter(is.finite(sh_chilcotin_recruits), is.finite(sh_chilcotin_spawners), sh_chilcotin_spawners > 0) %>%
  arrange(Year) %>%
  mutate(lnRS = log(sh_chilcotin_recruits / sh_chilcotin_spawners))

chum_observed_abundance <- data %>%
  filter(is.finite(chum_recruits_obs)) %>%
  select(Year, abundance = chum_recruits_obs)

thompson_observed_abundance <- sh_data %>%
  filter(is.finite(sh_thompson_prefishery_N)) %>%
  select(Year, abundance = sh_thompson_prefishery_N)

chilcotin_observed_abundance <- sh_data %>%
  filter(is.finite(sh_chilcotin_prefishery_N)) %>%
  select(Year, abundance = sh_chilcotin_prefishery_N)

get_low_period <- function(lnrs_df, abundance_df) {
  
  window <- lnrs_df %>%
    arrange(Year) %>%
    mutate(
      roll_mean = zoo::rollapply(lnRS, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA),
      roll_start = Year - LOW_PERIOD_YEARS + 1,
      min_year = min(Year)
    ) %>%
    filter(!is.na(roll_mean), roll_start >= min_year + MIN_YEARS_BEFORE_LOW_PERIOD) %>%
    slice_min(roll_mean, n = 1, with_ties = FALSE)
  
  low_period_mean <- abundance_df %>%
    filter(Year >= window$roll_start, Year <= window$Year) %>%
    summarise(m = mean(abundance, na.rm = TRUE)) %>% pull(m)
  
  tibble(low_period_start = window$roll_start, low_period_end = window$Year, low_period_mean = low_period_mean)
}

low_periods <- bind_rows(
  get_low_period(chum_observed_lnrs, chum_observed_abundance)           %>% mutate(Stock = "Chum"),
  get_low_period(thompson_observed_lnrs, thompson_observed_abundance)   %>% mutate(Stock = "Thompson steelhead"),
  get_low_period(chilcotin_observed_lnrs, chilcotin_observed_abundance) %>% mutate(Stock = "Chilcotin steelhead")
)
print(low_periods)

# Diagnostic plot: rolling 10-year mean ln(R/S) per stock, with the
# selected low-period window shaded -- matches the sockeye ln(R/S)
# diagnostic. Red dot marks the stock's own rolling lnRS at the window's
# end (not the abundance-based low_period_mean, since the axis here is
# productivity, not abundance).
observed_lnrs_rolling <- bind_rows(
  chum_observed_lnrs      %>% mutate(Stock = "Chum"),
  thompson_observed_lnrs  %>% mutate(Stock = "Thompson steelhead"),
  chilcotin_observed_lnrs %>% mutate(Stock = "Chilcotin steelhead")
) %>%
  group_by(Stock) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(roll_mean_lnRS = zoo::rollapply(lnRS, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA)) %>%
  ungroup()

lnrs_marker <- observed_lnrs_rolling %>%
  inner_join(low_periods %>% select(Stock, low_period_end), by = c("Stock", "Year" = "low_period_end"))

ggplot(observed_lnrs_rolling, aes(Year, roll_mean_lnRS)) +
  geom_line(color = "#4682B4") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.15) +
  geom_point(data = lnrs_marker, aes(x = Year, y = roll_mean_lnRS), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "10-year rolling mean ln(R/S)") +
  theme_minimal()

ggsave("figures/steelhead_chum_lnrs_low_period_diagnostic.png", width = 10, height = 4, dpi = 600)

natural_data_end <- list(
  Chum = max(chum_observed_abundance$Year),
  `Thompson steelhead` = max(thompson_observed_abundance$Year),
  `Chilcotin steelhead` = max(chilcotin_observed_abundance$Year)
)

# ------------------------------------------------------------
# 7b. SST BASELINE, per steelhead stock -- 1950-1975 mean, matching
#     sockeye's convention. steelhead_data.csv's SST column is already
#     STANDARDIZED (z-scored), and we don't know the exact mean/SD used
#     to build it -- so rather than guess, recover the transformation
#     directly via linear regression against the raw SST anomaly series
#     (np_temp), over whichever years overlap with each stock's data.
#     Standardization is an affine transform, so a linear fit recovers
#     it exactly: standardized = raw*slope + intercept, where
#     slope = 1/SD and intercept = -mean/SD. R^2 should come out very
#     close to 1 if this is genuinely the same transform -- if it
#     isn't, don't trust the resulting baseline.
# ------------------------------------------------------------

np_temp <- read_csv("Data/ersst_MJJ_gulf_of_alaska_anomaly.csv") %>%
  rename(Year = year, sst = sst_anom_mjj)

recover_sst_standardization <- function(stock_data, np_temp) {
  merged <- stock_data %>%
    select(Year, SST_standardized = SST) %>%
    inner_join(np_temp, by = "Year")
  
  fit <- lm(SST_standardized ~ sst, data = merged)
  r_squared <- summary(fit)$r.squared
  
  cat("Recovered SST standardization -- n overlapping years:", nrow(merged),
      " | R^2:", round(r_squared, 6), "\n")
  if (r_squared < 0.999) {
    warning("R^2 is below 0.999 -- the recovered transform may not exactly match the ",
            "original standardization. Check before trusting the SST baseline below.")
  }
  
  list(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]), r_squared = r_squared)
}

thompson_sst_transform  <- recover_sst_standardization(thompson_data, np_temp)
chilcotin_sst_transform <- recover_sst_standardization(chilcotin_data, np_temp)

apply_sst_transform <- function(raw_sst, transform) transform$intercept + transform$slope * raw_sst

SST_BASELINE_YEARS <- 1950:1975   # matches the sockeye convention

sst_baseline_thompson <- np_temp %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  mutate(sst_standardized = apply_sst_transform(sst, thompson_sst_transform)) %>%
  summarise(m = mean(sst_standardized, na.rm = TRUE)) %>% pull(m)

sst_baseline_chilcotin <- np_temp %>%
  filter(Year %in% SST_BASELINE_YEARS) %>%
  mutate(sst_standardized = apply_sst_transform(sst, chilcotin_sst_transform)) %>%
  summarise(m = mean(sst_standardized, na.rm = TRUE)) %>% pull(m)

cat("SST baseline (1950-1975, standardized units) -- Thompson:", sst_baseline_thompson,
    " | Chilcotin:", sst_baseline_chilcotin, "\n")
cat("For comparison, each stock's own historical SST range:\n")
print(summary(thompson_data$SST))
print(summary(chilcotin_data$SST))

# ------------------------------------------------------------
# 8. MONTE CARLO LOOP
#
#    CHUM: unchanged, own 2-way SSL control / no control comparison.
#
#    STEELHEAD (Thompson, Chilcotin): now TWO independent scenarios,
#    matching the sockeye design --
#      "Pinniped scenario": SSL frozen at 1978 value, SST real
#      "SST scenario":       SST frozen at each stock's own first-10-
#                            year baseline, SL real
#    Each compared against the SAME Xmin (fixed real data), same as
#    sockeye's Pinniped/SST/Pink comparison. Chum is still run once per
#    draw to build all_data (needed for the join), but its SSL_control
#    choice has NO effect on steelhead's results here -- this
#    retrospective always uses U_historic = 1, so the
#    chum_commercial_harvest_uapply interaction term never actually
#    fires (that branch only triggers when U_historic == 0). Chum is
#    run with SSL_control = 0 for this purpose, arbitrarily, since it's
#    inconsequential to the steelhead scenarios.
# ------------------------------------------------------------

lowpoint_results <- vector("list", N_DRAWS * 12)
counter <- 0

for (d in seq_len(N_DRAWS)) {
  
  ch <- chilcotin_draws[d, ]
  th <- thompson_draws[d, ]
  cm <- chum_draws[d, ]
  
  # ---- CHUM: fully crossed SSL scenario x harvest level ----
  chum_scenario_defs <- list("No SSL control" = 0, "SSL control" = 1)
  chum_harvest_defs  <- list(
    "Actual harvest" = list(U_apply = 0, U_historic = 1),
    "Zero harvest"   = list(U_apply = 0, U_historic = 0)
  )
  
  final_end <- natural_data_end[["Chum"]]
  final_start <- final_end - FINAL_PERIOD_YEARS + 1
  low_p <- low_periods %>% filter(Stock == "Chum")
  
  for (scenario_name in names(chum_scenario_defs)) {
    for (harvest_name in names(chum_harvest_defs)) {
      
      ssl_val <- chum_scenario_defs[[scenario_name]]
      hv <- chum_harvest_defs[[harvest_name]]
      
      chum_df <- run_chum_scenario_v2(
        data, chum_covariates, U_apply = hv$U_apply, SSL_control = ssl_val, U_historic = hv$U_historic,
        chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
      )
      
      X <- chum_df %>%
        filter(Year >= final_start, Year <= final_end) %>%
        summarise(m = mean(sum_alt, na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = "Chum", scenario = scenario_name, harvest = harvest_name, draw = d,
        Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean
      )
    }
  }
  
  # ---- STEELHEAD: fully crossed scenario x harvest level ----
  chum_df_for_steelhead <- run_chum_scenario_v2(
    data, chum_covariates, U_apply = 0, SSL_control = 0, U_historic = 1,
    chum_intercept = cm$a, chum_pdo_coef = cm$t, chum_ssl_coef = cm$s, chum_spawners_coef = -cm$b
  )
  all_data <- chum_df_for_steelhead %>% right_join(sh_data, by = "Year")
  
  steelhead_scenario_defs <- list(
    "SSL control" = list(SSL_control = 1, SST_control = 0),
    "SST scenario" = list(SSL_control = 0, SST_control = 1)
  )
  harvest_levels <- list("Actual harvest" = FALSE, "Zero harvest" = TRUE)
  
  for (scenario_name in names(steelhead_scenario_defs)) {
    
    sc <- steelhead_scenario_defs[[scenario_name]]
    
    for (harvest_name in names(harvest_levels)) {
      
      zh <- harvest_levels[[harvest_name]]
      
      thompson_df <- run_thompson_scenario_v2(
        all_data, SSL_control = sc$SSL_control, SST_control = sc$SST_control,
        U_historic = 1, byrate = default_byrate,
        sh_thompson_intercept = th$a, sh_thompson_sst_coef = th$t,
        sh_thompson_ssl_coef = th$s, sh_thompson_spawners_coef = -th$b,
        sst_baseline_thompson = sst_baseline_thompson, zero_harvest = zh
      )
      
      chilcotin_df <- run_chilcotin_scenario_v2(
        thompson_df, SSL_control = sc$SSL_control, SST_control = sc$SST_control,
        U_historic = 1, byrate = default_byrate,
        sh_chilcotin_intercept = ch$a, sh_chilcotin_sst_coef = ch$t,
        sh_chilcotin_ssl_coef = ch$s, sh_chilcotin_spawners_coef = -ch$b,
        sst_baseline_chilcotin = sst_baseline_chilcotin, zero_harvest = zh
      )
      
      for (stock_name in c("Thompson steelhead", "Chilcotin steelhead")) {
        
        abundance_col <- switch(stock_name,
                                "Thompson steelhead" = "sh_thompson_sum_pred",
                                "Chilcotin steelhead" = "sh_chilcotin_sum_pred"
        )
        
        final_end_sh <- natural_data_end[[stock_name]]
        final_start_sh <- final_end_sh - FINAL_PERIOD_YEARS + 1
        low_p_sh <- low_periods %>% filter(Stock == stock_name)
        
        X <- chilcotin_df %>%
          filter(Year >= final_start_sh, Year <= final_end_sh) %>%
          summarise(m = mean(.data[[abundance_col]], na.rm = TRUE)) %>% pull(m)
        
        counter <- counter + 1
        lowpoint_results[[counter]] <- tibble(
          Stock = stock_name, scenario = scenario_name, harvest = harvest_name, draw = d,
          Xmin = low_p_sh$low_period_mean, X = X, ratio = X / low_p_sh$low_period_mean
        )
      }
    }
  }
  
  if (d %% 10 == 0) message("draw ", d, " / ", N_DRAWS)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "steelhead_chum_ratio_summary.rds")

# Diagnostic plot: distribution of X/Xmin draws, faceted by stock,
# scenario, AND harvest level now that harvest is its own crossed
# dimension -- same style as the coho/sockeye ratio histograms.
SCENARIO_COLORS <- c(
  "No SSL control" = "black", "SSL control" = "#4682B4", "SST scenario" = "#2E8B57"
)

ggplot(lowpoint_summary, aes(ratio, fill = scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS, name = NULL) +
  facet_wrap(~ Stock + scenario + harvest, scales = "free", ncol = 4) +
  labs(x = "Recruits / min Recruits", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_histogram.png", width = 16, height = 10, dpi = 600)

# ------------------------------------------------------------
# 10. 40% Smsy BENCHMARK -- per-year actual covariates, averaged output,
#     same approach as coho/sockeye. Reference periods, per your
#     instruction:
#       Chum: 1960-1970
#       Thompson/Chilcotin: beginning of each stock's own time series
#                           through 1980
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

CHUM_SMSY_YEARS <- 1960:1970
STEELHEAD_SMSY_END_YEAR <- 1980

chilcotin_data_smsy_period <- chilcotin_data %>% filter(Year <= STEELHEAD_SMSY_END_YEAR)
thompson_data_smsy_period  <- thompson_data %>% filter(Year <= STEELHEAD_SMSY_END_YEAR)
chum_covariates_smsy_period <- chum_covariates %>%
  filter(!is.na(SL), !is.na(PDO), Year %in% CHUM_SMSY_YEARS)

cat("Smsy reference period year counts -- Chilcotin:", nrow(chilcotin_data_smsy_period),
    " | Thompson:", nrow(thompson_data_smsy_period),
    " | Chum:", nrow(chum_covariates_smsy_period), "\n")

smsy_benchmarks <- bind_rows(
  compute_smsy_benchmark(
    chilcotin_draws, chilcotin_data_smsy_period, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Chilcotin steelhead"], scale = 1000
  ) %>% mutate(Stock = "Chilcotin steelhead"),
  compute_smsy_benchmark(
    thompson_draws, thompson_data_smsy_period, c("SL", "SST"), c("s", "t"),
    Xmin = low_periods$low_period_mean[low_periods$Stock == "Thompson steelhead"], scale = 1000
  ) %>% mutate(Stock = "Thompson steelhead"),
  compute_smsy_benchmark(
    chum_draws, chum_covariates_smsy_period, c("SL", "PDO"), c("s", "t"),
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
  filter(scenario != "No SSL control") %>%   # main plot: alternative scenarios only
  group_by(Stock, scenario, harvest) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.10, na.rm = TRUE),
    q95    = quantile(ratio, 0.90, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species_type = if_else(Stock == "Chum", "Chum", "Steelhead"))
print(mc_summary_stats_ratio, n = Inf)

# Smsy is a fixed biological reference, independent of the harvest axis
# -- duplicated across both harvest levels so it renders in both facet
# rows for comparison, rather than only appearing in one.
smsy_benchmarks <- bind_rows(
  smsy_benchmarks %>% mutate(harvest = "Actual harvest"),
  smsy_benchmarks %>% mutate(harvest = "Zero harvest")
) %>%
  mutate(species_type = if_else(Stock == "Chum", "Chum", "Steelhead"))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.4), size = 0.6, linewidth = 1) +
  geom_point(data = smsy_benchmarks, aes(x = Stock, y = benchmark_ratio, shape = "40% Smsy"),
             inherit.aes = FALSE, color = "firebrick", size = 6) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  facet_grid(harvest ~ species_type, scales = "free") +
  scale_color_manual(values = SCENARIO_COLORS, name = NULL) +
  scale_shape_manual(name = NULL, values = c("40% Smsy" = 95)) +   # pch 95 = horizontal dash
  scale_y_continuous(trans = "log", labels = scales::comma, breaks = scales::breaks_log(n = 6)) +
  labs(x = NULL, y = "Recruits / min Recruits") +
  theme_light() +
  theme(legend.position = "bottom")

ggsave("figures/steelhead_chum_ratio_harvest.png", width = 10, height = 10, dpi = 600)
