# ============================================================
# Probabilistic (Monte Carlo) retrospective analysis --
# Thompson & Chilcotin steelhead
#
# For each of N_DRAWS posterior draws of (a, b, s, t[, f]) from the JAGS
# stock-recruit models, re-runs the age-structured retrospective scenario
# script with those draw values substituted for the fixed point-estimate
# coefficients, for the historic-harvest scenario (SSL control vs. none,
# byrate = default_byrate, harvest held at the historic reconstructed
# rate i.e. U_historic = 1) -- scenario 1 from the original retrospective
# script. No bycatch-rate or harvest-rate sweep.
#
# Chum is NOT drawn probabilistically here -- it runs with U_historic = 1
# (harvest pinned to the historic reconstructed rate) and its own fixed
# point-estimate coefficients, identical across every draw.
#
# DEPENDENCIES -- run chum_steelhead_retrospective_scenarios.R (or at
# least its data-loading + run_chum_scenario() section) in this session
# first, so `data`, `covariates`, `sh_data`, `MODEL_YEARS`, and
# `default_byrate` already exist. Also run the two JAGS scripts (or
# readRDS their saved posterior samples) so
# `chilcotin_samples_kept` / `thompson_samples_kept` exist.
#
# CAVEATS -- please check these against your actual setup:
#  1. The max_flow join below assumes sh_data is keyed by Year with
#     columns like sh_chilcotin_SL/SST etc., as in the original
#     retrospective script. Adjust if your actual sh_data differs.
#  2. Runtime: 5000 draws x 2 SSL scenarios x ~50-iteration fixed-point
#     solve per year, per stock, in a row-by-row data.frame loop -- slow
#     in base R. Test with a small N_DRAWS (e.g. 50) first to gauge
#     runtime before committing to the full 5000. Draws are independent,
#     so this also parallelizes well -- see the note at the bottom.
#  3. Some posterior draws imply weak enough productivity that the
#     no-control reconstruction hits the model's spawner floor of 0
#     during the stock's low-abundance years -- a real feature of that
#     draw, not a bug, but it makes % increase (division by ~0) explode.
#     Primary reported metric is now the ABSOLUTE increase in mean
#     abundance (never blows up); the fraction of draws that collapsed
#     this way is reported separately (collapse_fraction) as its own
#     diagnostic, and a %-increase summary excluding collapsed draws is
#     kept for reference only. Adjust `collapse_threshold` in
#     summarize_pct_increase() (default 1 fish) if you want a different
#     cutoff for "collapsed."
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 5000
N_DRAWS <- 50

# ------------------------------------------------------------
# 1. LOAD POSTERIOR SAMPLES
# ------------------------------------------------------------

if (!exists("chilcotin_samples_kept")) {
  chilcotin_samples_kept <- readRDS("chilcotin_posterior_samples.rds")
}
if (!exists("thompson_samples_kept")) {
  thompson_samples_kept <- readRDS("thompson_posterior_samples.rds")
}

chilcotin_post <- as.data.frame(as.matrix(chilcotin_samples_kept[, c("a", "b", "s", "t", "f")]))
thompson_post  <- as.data.frame(as.matrix(thompson_samples_kept[, c("a", "b", "s", "t")]))

# ------------------------------------------------------------
# 2. SUBSAMPLE N_DRAWS FROM EACH POSTERIOR
# ------------------------------------------------------------
# Thompson and Chilcotin posteriors were fit independently -- draw i of
# one is paired with draw i of the other purely to keep the loop to
# N_DRAWS total runs rather than N_DRAWS^2.

draw_idx_chilcotin <- sample(seq_len(nrow(chilcotin_post)), N_DRAWS,
                             replace = nrow(chilcotin_post) < N_DRAWS)
draw_idx_thompson  <- sample(seq_len(nrow(thompson_post)), N_DRAWS,
                             replace = nrow(thompson_post) < N_DRAWS)

chilcotin_draws <- chilcotin_post[draw_idx_chilcotin, ]
thompson_draws  <- thompson_post[draw_idx_thompson, ]

# ------------------------------------------------------------
# 3. UPDATED run_thompson_scenario() / run_chilcotin_scenario(),
#    matching the covariate sets fit in JAGS: Thompson = Spawners, SL,
#    SST (drops NPGO); Chilcotin = Spawners, max_flow, SL, SST (drops PDO)
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
      df$sh_thompson_SST[i] * sh_thompson_sst_coef +
      df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    df$sh_thompson_alpha_CN[i] <-
      sh_thompson_intercept +
      df$sh_thompson_SST[i] * sh_thompson_sst_coef +
      df$sh_thompson_SL[i]  * sh_thompson_ssl_coef
    
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_thompson_ln_obs_pred[i])) 1 else exp(df$sh_thompson_ln_obs_pred[i])
      
      df$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_thompson_alpha_CN[i] + sh_thompson_spawners_coef * spk) * resid_adj
        } else {
          spk * exp(df$sh_thompson_SSL_alpha[i] + sh_thompson_spawners_coef * spk) * resid_adj
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

run_chilcotin_scenario_v2 <- function(df, SSL_control, U_historic, byrate,
                                      start_year = 1973,
                                      sh_chilcotin_intercept,
                                      sh_chilcotin_sst_coef,
                                      sh_chilcotin_ssl_coef,
                                      sh_chilcotin_maxflow_coef,
                                      sh_chilcotin_spawners_coef) {
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  FN_chilcotin_2018 <- df %>% filter(Year == 2018) %>% pull(sh_chilcotin_FN_mortalities) %>% as.numeric()
  
  df <- df %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_chilcotin_base_alpha = sh_chilcotin_intercept +
        sh_chilcotin_SST * sh_chilcotin_sst_coef +
        sh_chilcotin_SL  * sh_chilcotin_ssl_coef +
        sh_chilcotin_max_flow * sh_chilcotin_maxflow_coef,
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
      df$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef +
      df$sh_chilcotin_max_flow[i] * sh_chilcotin_maxflow_coef +
      df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    df$sh_chilcotin_alpha_CN[i] <-
      sh_chilcotin_intercept +
      df$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef +
      df$sh_chilcotin_SL[i]  * sh_chilcotin_ssl_coef +
      df$sh_chilcotin_max_flow[i] * sh_chilcotin_maxflow_coef
    
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    for (iter in seq_len(50)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      resid_adj <- if (is.na(df$sh_chilcotin_ln_obs_pred[i])) 1 else exp(df$sh_chilcotin_ln_obs_pred[i])
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_chilcotin_alpha_CN[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
        } else {
          spk * exp(df$sh_chilcotin_SSL_alpha[i] + sh_chilcotin_spawners_coef * spk) * resid_adj
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

# ------------------------------------------------------------
# 4. JOIN max_flow INTO sh_data -- VERIFY THIS AGAINST YOUR ACTUAL sh_data
# ------------------------------------------------------------

steelhead_data <- read_csv("data/steelhead_data.csv")

chilcotin_maxflow_by_year <- steelhead_data %>%
  filter(Stock == "Chilcotin") %>%
  transmute(Year, sh_chilcotin_max_flow = max_flow)

sh_data <- sh_data %>%
  left_join(chilcotin_maxflow_by_year, by = "Year")

# ------------------------------------------------------------
# 5. PRECOMPUTE CHUM ONCE PER SSL_control -- valid for BOTH scenario
#    types here, since both use U_historic = 1 (chum harvest pinned to
#    the historic rate regardless of byrate)
# ------------------------------------------------------------

chum_df_no_control <- run_chum_scenario(data, covariates, U_apply = 0, SSL_control = 0, U_historic = 1)
chum_df_control     <- run_chum_scenario(data, covariates, U_apply = 0, SSL_control = 1, U_historic = 1)

all_data_no_control <- chum_df_no_control %>% right_join(sh_data, by = "Year")
all_data_control     <- chum_df_control     %>% right_join(sh_data, by = "Year")

# ------------------------------------------------------------
# 6. HELPER: run Thompson + Chilcotin for one draw / SSL_control / byrate
# ------------------------------------------------------------

run_draw_scenario <- function(th, ch, ssl_control, byrate) {
  
  base_all_data <- if (ssl_control == 0) all_data_no_control else all_data_control
  
  thompson_df <- run_thompson_scenario_v2(
    base_all_data, SSL_control = ssl_control, U_historic = 1, byrate = byrate,
    sh_thompson_intercept     = th$a,
    sh_thompson_sst_coef      = th$t,
    sh_thompson_ssl_coef      = th$s,
    sh_thompson_spawners_coef = -th$b
  )
  
  chilcotin_df <- run_chilcotin_scenario_v2(
    thompson_df, SSL_control = ssl_control, U_historic = 1, byrate = byrate,
    sh_chilcotin_intercept      = ch$a,
    sh_chilcotin_sst_coef       = ch$t,
    sh_chilcotin_ssl_coef       = ch$s,
    sh_chilcotin_maxflow_coef   = ch$f,
    sh_chilcotin_spawners_coef  = -ch$b
  )
  
  ssl_label <- if (ssl_control == 1) "SSL control" else "No SSL control"
  
  chilcotin_df %>%
    filter(Year %in% MODEL_YEARS) %>%
    transmute(
      Year,
      ssl_scenario = ssl_label,
      thompson_productivity  = log(sh_thompson_recruits_alt / (sh_thompson_spawners_pred / 1000)),
      thompson_returns       = sh_thompson_sum_pred,
      chilcotin_productivity = log(sh_chilcotin_recruits_alt / (sh_chilcotin_spawners_pred / 1000)),
      chilcotin_returns      = sh_chilcotin_sum_pred
    )
}

summarize_pct_increase <- function(draw_df, draw_num, scenario_type, byrate_val,
                                   collapse_threshold = 1) {
  
  long_abund <- bind_rows(
    draw_df %>% transmute(Year, ssl_scenario, Stock = "Thompson steelhead", abundance = thompson_returns),
    draw_df %>% transmute(Year, ssl_scenario, Stock = "Chilcotin steelhead", abundance = chilcotin_returns)
  )
  
  low_years <- long_abund %>%
    filter(ssl_scenario == "No SSL control") %>%
    group_by(Stock) %>%
    filter(abundance <= quantile(abundance, 0.25, na.rm = TRUE)) %>%
    distinct(Stock, Year)
  
  long_abund %>%
    inner_join(low_years, by = c("Stock", "Year")) %>%
    group_by(Stock, ssl_scenario) %>%
    summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = ssl_scenario, values_from = mean_abundance) %>%
    mutate(
      abs_increase = `SSL control` - `No SSL control`,
      # % increase kept for reference, but it's unstable whenever the
      # no-control baseline is at/near zero (see collapse flag below) --
      # don't average this across draws without excluding those.
      pct_increase = (`SSL control` - `No SSL control`) / `No SSL control` * 100,
      # TRUE when this draw's no-control reconstruction has collapsed to
      # (near) zero abundance during the stock's own low-abundance years --
      # i.e. the population hit the spawner floor of 0 in the model. This
      # is a real feature of that posterior draw, not a numerical error,
      # but it makes pct_increase meaningless for that draw.
      collapse = `No SSL control` <= collapse_threshold,
      draw = draw_num, scenario_type = scenario_type, byrate = byrate_val
    )
}

# ------------------------------------------------------------
# 7. MONTE CARLO LOOP -- historic harvest scenario only
# ------------------------------------------------------------

mc_results_full    <- vector("list", N_DRAWS)
mc_results_summary <- vector("list", N_DRAWS)

for (d in seq_len(N_DRAWS)) {
  
  th <- thompson_draws[d, ]
  ch <- chilcotin_draws[d, ]
  
  historic_pair <- bind_rows(lapply(c(0, 1), function(ssl) {
    run_draw_scenario(th, ch, ssl_control = ssl, byrate = default_byrate)
  })) %>%
    mutate(draw = d, scenario_type = "Historic harvest", byrate = default_byrate)
  
  mc_results_full[[d]]    <- historic_pair
  mc_results_summary[[d]] <- summarize_pct_increase(historic_pair, d, "Historic harvest", default_byrate)
  
  if (d %% 500 == 0) message("draw ", d, " / ", N_DRAWS)
}

mc_full    <- bind_rows(mc_results_full)
mc_summary <- bind_rows(mc_results_summary)

saveRDS(mc_full,    "mc_retrospective_full_timeseries.rds")
saveRDS(mc_summary, "mc_retrospective_pct_increase_summary.rds")

# ------------------------------------------------------------
# 8. POSTERIOR SUMMARIES
# ------------------------------------------------------------

# Per-year posterior of productivity/returns, by stock and SSL scenario
mc_year_summary <- mc_full %>%
  pivot_longer(cols = c(thompson_productivity, thompson_returns,
                        chilcotin_productivity, chilcotin_returns),
               names_to = "metric", values_to = "value") %>%
  group_by(Year, ssl_scenario, metric) %>%
  summarise(
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q2.5   = quantile(value, 0.025, na.rm = TRUE),
    q97.5  = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Fraction of draws where the no-control reconstruction collapsed to
# (near) zero abundance during the stock's own low-abundance years --
# reported separately since it makes pct_increase meaningless for those
# draws, but is itself a meaningful finding about the posterior
collapse_fraction <- mc_summary %>%
  group_by(Stock) %>%
  summarise(
    n_draws          = n(),
    n_collapsed      = sum(collapse, na.rm = TRUE),
    pct_collapsed    = 100 * n_collapsed / n_draws,
    .groups = "drop"
  )
print(collapse_fraction)

# Posterior of the ABSOLUTE increase in mean abundance over the
# low-abundance period -- primary metric, since it doesn't blow up when
# the no-control baseline is near zero
mc_summary_stats <- mc_summary %>%
  group_by(Stock) %>%
  summarise(
    mean   = mean(abs_increase, na.rm = TRUE),
    median = median(abs_increase, na.rm = TRUE),
    q2.5   = quantile(abs_increase, 0.025, na.rm = TRUE),
    q97.5  = quantile(abs_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats)

ggplot(mc_summary_stats, aes(x = reorder(Stock, -mean), y = mean, fill = "Pinniped scenario")) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Pinniped scenario" = "#4682B4"), guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL,
       y = "Increase in mean abundance over low-abundance period (fish)\n(posterior mean, 95% credible interval)",
       title = "Steelhead abundance recovery under the pinniped-control scenario -- full posterior") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/steelhead_abs_increase_low_period_posterior.png", width = 7, height = 5.5, dpi = 600)

# Reference only: % increase, excluding collapsed draws, so it isn't
# dominated by division-by-near-zero. Still treat cautiously -- excluding
# collapsed draws changes what population of draws this represents.
mc_summary_stats_pct_noncollapsed <- mc_summary %>%
  filter(!collapse) %>%
  group_by(Stock) %>%
  summarise(
    n_used = n(),
    mean   = mean(pct_increase, na.rm = TRUE),
    median = median(pct_increase, na.rm = TRUE),
    q2.5   = quantile(pct_increase, 0.025, na.rm = TRUE),
    q97.5  = quantile(pct_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_pct_noncollapsed)

# ------------------------------------------------------------
# OPTIONAL: parallelize across draws
# ------------------------------------------------------------
# Each draw's computation is independent of every other draw, so if this
# is too slow as a plain for loop, wrap the body of section 7 in
# future.apply::future_lapply() or furrr::future_map() over
# seq_len(N_DRAWS) with a multisession/multicore plan. Worth doing once
# you've confirmed correctness and timed a small N_DRAWS run.