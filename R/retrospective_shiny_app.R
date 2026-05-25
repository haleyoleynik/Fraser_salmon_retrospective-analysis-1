# app.R
library(shiny)
library(tidyverse)
library(slider)

#-------------------------------
# Model runner (wraps  script)
#-------------------------------
run_model <- function(
    U_apply = 0.2,
    byrate = 0.69,
    SSL_control = 0,
    U_historic = 1,
    data_fp = "data/s-r_data.csv",
    sh_fp   = "data/sh_s-r_data.csv",
    cov_fp  = "data/covariates.csv"
) {
  
  # read data
  data <- readr::read_csv(data_fp, show_col_types = FALSE)
  sh_data <- readr::read_csv(sh_fp, show_col_types = FALSE)
  covariates <- readr::read_csv(cov_fp, show_col_types = FALSE)
  
  #-------------------------------
  # CHUM
  #-------------------------------
  intercept       <- 1.03737862843252
  pdo_adult_coef  <- 0.093
  npgo_coef       <- 0.103
  pdo_smolt_coef  <- -0.106
  ssl_coef        <- -0.224
  spawners_coef   <- -4.95E-07
  
  # your "new.data"
  new.data <- data %>%
    left_join(covariates, by = "Year") %>%
    arrange(Year) %>%
    mutate(
      catch = chum_total_stock - chum_spawners,
      U_chum = catch / chum_total_stock,
      chum_base_alpha = intercept + PDO_adult * pdo_adult_coef + NPGO * npgo_coef +
        PDO_smolt * pdo_smolt_coef + SSL * ssl_coef,
      alpha_running_avg = slide_dbl(chum_base_alpha, mean, .before = 9, .complete = TRUE),
      chum_model_recruits = chum_spawners * exp(chum_base_alpha + chum_spawners * spawners_coef),
      chum_ln_obs_pred = log(chum_recruits_obs / chum_model_recruits),
      Nage3_obs = chum_total_stock * prop3,
      Nage4_obs = chum_total_stock * prop4,
      Nage5_obs = chum_total_stock * prop5,
      Nage6_obs = chum_total_stock * prop6,
      Nage3_pred = case_when(
        Year <= 1953 ~ Nage3_obs,
        Year >= 1954 ~ lag(chum_model_recruits, 3) * prop3 * exp(lag(chum_ln_obs_pred, 3))
      ),
      Nage4_pred = case_when(
        Year <= 1954 ~ Nage4_obs,
        Year >= 1955 ~ lag(chum_model_recruits, 4) * prop4 * exp(lag(chum_ln_obs_pred, 4))
      ),
      Nage5_pred = case_when(
        Year <= 1955 ~ Nage5_obs,
        Year >= 1956 ~ lag(chum_model_recruits, 5) * prop5 * exp(lag(chum_ln_obs_pred, 5))
      ),
      Nage6_pred = case_when(
        Year <= 1956 ~ Nage6_obs,
        Year >= 1957 ~ lag(chum_model_recruits, 6) * prop6 * exp(lag(chum_ln_obs_pred, 6))
      ),
      recruits_pred = rowSums(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)), na.rm = TRUE),
      recruits_dif = chum_total_stock - recruits_pred,
      U_chum_pred = catch / recruits_pred,
      U_chum_dif = U_chum_pred - U_chum,
      chum_commercial_harvest = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 ~ U_chum
      ),
      chum_commercial_harvest_uapply = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 ~ U_apply
      )
    )
  
  SSL_1978 <- new.data %>% filter(Year == 1978) %>% pull(SSL)
  if (length(SSL_1978) == 0) SSL_1978 <- NA_real_
  
  df_chum <- new.data %>%
    arrange(Year) %>%
    mutate(
      chum_recruits_alt = NA_real_,
      Nage3_alt = Nage3_pred,
      Nage4_alt = Nage4_pred,
      Nage5_alt = Nage5_pred,
      Nage6_alt = Nage6_pred,
      sum_alt = NA_real_,
      catch_alt = NA_real_,
      chum_spawners_pred = chum_spawners
    )
  
  for (i in seq_len(nrow(df_chum))) {
    
    SSL_alt <- if (df_chum$Year[i] <= 1978) {
      df_chum$SSL[i]
    } else {
      (1 - SSL_control) * df_chum$SSL[i] + SSL_control * SSL_1978
    }
    
    SSL_control_alpha <-
      intercept +
      df_chum$PDO_adult[i] * pdo_adult_coef +
      df_chum$NPGO[i] * npgo_coef +
      df_chum$PDO_smolt[i] * pdo_smolt_coef +
      SSL_alt * ssl_coef
    
    df_chum$chum_recruits_alt[i] <-
      if (SSL_control == 0) {
        df_chum$chum_spawners_pred[i] *
          exp(df_chum$chum_base_alpha[i] + spawners_coef * df_chum$chum_spawners_pred[i]) *
          exp(df_chum$chum_ln_obs_pred[i])
      } else {
        df_chum$chum_spawners_pred[i] *
          exp(SSL_control_alpha + spawners_coef * df_chum$chum_spawners_pred[i]) *
          exp(df_chum$chum_ln_obs_pred[i])
      }
    
    if (i > 3) df_chum$Nage3_alt[i] <- df_chum$chum_recruits_alt[i - 3] * df_chum$prop3[i]
    if (i > 4) df_chum$Nage4_alt[i] <- df_chum$chum_recruits_alt[i - 4] * df_chum$prop4[i]
    if (i > 5) df_chum$Nage5_alt[i] <- df_chum$chum_recruits_alt[i - 5] * df_chum$prop5[i]
    if (i > 6) df_chum$Nage6_alt[i] <- df_chum$chum_recruits_alt[i - 6] * df_chum$prop6[i]
    
    df_chum$sum_alt[i] <- sum(df_chum$Nage3_alt[i], df_chum$Nage4_alt[i], df_chum$Nage5_alt[i], df_chum$Nage6_alt[i], na.rm = TRUE)
    df_chum$catch_alt[i] <- df_chum$sum_alt[i] * df_chum$U_chum_pred[i]
    df_chum$chum_spawners_pred[i] <- df_chum$sum_alt[i] * (1 - df_chum$U_chum_pred[i])
  }
  
  #-------------------------------
  # Merge chum -> steelhead input
  #-------------------------------
  all_data <- df_chum %>% right_join(sh_data, by = "Year")
  
  #-------------------------------
  # STEELHEAD: Thompson
  #-------------------------------
  sh_thompson_intercept     <- 1.572107637
  sh_thompson_sst_coef      <- -0.203463091
  sh_thompson_ssl_coef      <- -0.764277677
  sh_thompson_npgo_coef     <- -0.0402
  sh_thompson_spawners_coef <- -0.804438793
  start_year_thompson       <- 1978
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  if (length(sh_thompson_SSL_1978) == 0) sh_thompson_SSL_1978 <- NA_real_
  
  FN_thompson_2018 <- all_data %>% filter(Year == 2018) %>% pull(sh_thompson_FN_mortalities) %>% as.numeric()
  if (length(FN_thompson_2018) == 0) FN_thompson_2018 <- NA_real_
  
  df <- all_data %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_thompson_base_alpha =
        sh_thompson_intercept +
        sh_thompson_SST  * sh_thompson_sst_coef +
        sh_thompson_SL   * sh_thompson_ssl_coef +
        sh_thompson_NPGO * sh_thompson_npgo_coef,
      
      sh_thompson_model_recruits =
        sh_thompson_spawners *
        exp(sh_thompson_base_alpha + sh_thompson_spawners * sh_thompson_spawners_coef),
      
      sh_thompson_ln_obs_pred =
        log(sh_thompson_recruits / sh_thompson_model_recruits),
      
      sh_thompson_pred_bycatch =
        sh_thompson_prefishery_N -
        sh_thompson_sport_mortalities -
        sh_thompson_FN_mortalities -
        1000 * sh_thompson_spawners,
      
      sh_thompson_U =
        sh_thompson_pred_bycatch / sh_thompson_prefishery_N,
      
      sh_thompson_recruits_alt      = NA_real_,
      sh_thompson_spawners_pred     = sh_thompson_spawners,
      
      sh_thompson_Nage4_pred = NA_real_,
      sh_thompson_Nage5_pred = NA_real_,
      sh_thompson_Nage6_pred = NA_real_,
      sh_thompson_Nage7_pred = NA_real_,
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
  
  start_i <- which(df$Year >= start_year_thompson)[1]
  if (length(start_i) == 0 || is.na(start_i)) start_i <- 1
  
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_thompson_U_comm[i] <-
      if (yr <= 1990) {
        df$sh_thompson_U[i]
      } else if (U_historic == 1) {
        df$sh_thompson_U[i]
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
    
    # fixed-point iteration
    S_old <- df$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    max_iter <- 50
    tol <- 1e-8
    
    for (iter in seq_len(max_iter)) {
      
      df$sh_thompson_spawners_pred[i] <- S_old
      spk <- df$sh_thompson_spawners_pred[i] / 1000
      
      df$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_thompson_alpha_CN[i] + sh_thompson_spawners_coef * spk)
        } else {
          spk * exp(df$sh_thompson_SSL_alpha[i] + sh_thompson_spawners_coef * spk)
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df$sh_thompson_recruits_alt[j]
      }
      
      if (yr < start_year_thompson + 4) {
        df$sh_thompson_Nage4_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p4[i]
      } else {
        df$sh_thompson_Nage4_pred[i] <- lag_recruits(yr - 4) * df$sh_thompson_p4[i] * 1000
      }
      
      if (yr < start_year_thompson + 5) {
        df$sh_thompson_Nage5_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p5[i]
      } else {
        df$sh_thompson_Nage5_pred[i] <- lag_recruits(yr - 5) * df$sh_thompson_p5[i] * 1000
      }
      
      if (yr < start_year_thompson + 6) {
        df$sh_thompson_Nage6_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p6[i]
      } else {
        df$sh_thompson_Nage6_pred[i] <- lag_recruits(yr - 6) * df$sh_thompson_p6[i] * 1000
      }
      
      if (yr < start_year_thompson + 7) {
        df$sh_thompson_Nage7_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p7[i]
      } else {
        df$sh_thompson_Nage7_pred[i] <- lag_recruits(yr - 7) * df$sh_thompson_p7[i] * 1000
      }
      
      if (yr < start_year_thompson + 8) {
        df$sh_thompson_Nage8_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p8[i]
      } else {
        df$sh_thompson_Nage8_pred[i] <- lag_recruits(yr - 8) * df$sh_thompson_p8[i] * 1000
      }
      
      df$sh_thompson_sum_pred[i] <-
        sum(df$sh_thompson_Nage4_pred[i],
            df$sh_thompson_Nage5_pred[i],
            df$sh_thompson_Nage6_pred[i],
            df$sh_thompson_Nage7_pred[i],
            df$sh_thompson_Nage8_pred[i],
            na.rm = TRUE)
      
      df$sh_thompson_bycatch_pred[i] <- df$sh_thompson_sum_pred[i] * df$sh_thompson_U_comm[i]
      
      if (yr <= 2018) {
        df$sh_thompson_FN_catch_pred[i] <- df$sh_thompson_FN_mortalities[i]
      } else {
        denom_2018 <- (df$sh_thompson_sum_pred[df$Year == 2018] -
                         df$sh_thompson_bycatch_pred[df$Year == 2018])
        df$sh_thompson_FN_catch_pred[i] <-
          FN_thompson_2018 / denom_2018 *
          (df$sh_thompson_sum_pred[i] - df$sh_thompson_bycatch_pred[i])
      }
      
      df$sh_thompson_total_catch_pred[i] <-
        df$sh_thompson_FN_catch_pred[i] +
        df$sh_thompson_sport_mortalities[i] +
        df$sh_thompson_bycatch_pred[i]
      
      S_new <- df$sh_thompson_sum_pred[i] - df$sh_thompson_total_catch_pred[i]
      S_new <- max(S_new, 0)
      
      if (is.finite(S_old) && is.finite(S_new) &&
          abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    
    df$sh_thompson_spawners_pred[i] <- S_old
  }
  
  #-------------------------------
  # STEELHEAD: Chilcotin
  #-------------------------------
  sh_chilcotin_intercept     <- 1.053608979
  sh_chilcotin_sst_coef      <- -0.127949278
  sh_chilcotin_ssl_coef      <- -0.792741195
  sh_chilcotin_npgo_coef     <- 0.152526045
  sh_chilcotin_pdo_coef      <- 0.202708011
  sh_chilcotin_spawners_coef <- -1.022467631
  start_year_chilcotin       <- 1973
  
  sh_chilcotin_SSL_1978 <- df %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  if (length(sh_chilcotin_SSL_1978) == 0) sh_chilcotin_SSL_1978 <- NA_real_
  
  FN_chilcotin_2018 <- df %>% filter(Year == 2018) %>% pull(sh_chilcotin_FN_mortalities) %>% as.numeric()
  if (length(FN_chilcotin_2018) == 0) FN_chilcotin_2018 <- NA_real_
  
  df <- df %>%
    arrange(Year) %>%
    mutate(
      sh_chilcotin_base_alpha =
        sh_chilcotin_intercept +
        sh_chilcotin_SST  * sh_chilcotin_sst_coef +
        sh_chilcotin_SL   * sh_chilcotin_ssl_coef +
        sh_chilcotin_NPGO * sh_chilcotin_npgo_coef +
        sh_chilcotin_PDO  * sh_chilcotin_pdo_coef,
      
      sh_chilcotin_model_recruits =
        sh_chilcotin_spawners *
        exp(sh_chilcotin_base_alpha + sh_chilcotin_spawners * sh_chilcotin_spawners_coef),
      
      sh_chilcotin_ln_obs_pred =
        log(sh_chilcotin_recruits / sh_chilcotin_model_recruits),
      
      sh_chilcotin_pred_bycatch =
        sh_chilcotin_prefishery_N -
        sh_chilcotin_sport_mortalities -
        sh_chilcotin_FN_mortalities -
        1000 * sh_chilcotin_spawners,
      
      sh_chilcotin_U =
        sh_chilcotin_pred_bycatch / sh_chilcotin_prefishery_N,
      
      sh_chilcotin_recruits_alt      = NA_real_,
      sh_chilcotin_spawners_pred     = sh_chilcotin_spawners,
      
      sh_chilcotin_Nage4_pred = NA_real_,
      sh_chilcotin_Nage5_pred = NA_real_,
      sh_chilcotin_Nage6_pred = NA_real_,
      sh_chilcotin_Nage7_pred = NA_real_,
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
  
  start_i <- which(df$Year >= start_year_chilcotin)[1]
  if (length(start_i) == 0 || is.na(start_i)) start_i <- 1
  
  year_to_i <- setNames(seq_len(nrow(df)), df$Year)
  
  for (i in seq(from = start_i, to = nrow(df))) {
    
    yr <- df$Year[i]
    
    df$sh_chilcotin_U_comm[i] <-
      if (yr <= 1990) {
        df$sh_chilcotin_U[i]
      } else if (U_historic == 1) {
        df$sh_chilcotin_U[i]
      } else {
        byrate * df$chum_commercial_harvest_uapply[i]
      }
    
    df$sh_chilcotin_SSL_alt[i] <- if (yr <= start_year_chilcotin) {
      df$sh_chilcotin_SL[i]
    } else {
      (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
    }
    
    df$sh_chilcotin_SSL_alpha[i] <-
      sh_chilcotin_intercept +
      df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
      df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df$sh_chilcotin_PDO[i]  * sh_chilcotin_pdo_coef +
      df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    df$sh_chilcotin_alpha_CN[i] <-
      sh_chilcotin_intercept +
      df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
      df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df$sh_chilcotin_PDO[i]  * sh_chilcotin_pdo_coef +
      df$sh_chilcotin_SL[i]   * sh_chilcotin_ssl_coef
    
    # fixed-point iteration
    S_old <- df$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    max_iter <- 50
    tol <- 1e-8
    
    for (iter in seq_len(max_iter)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_chilcotin_alpha_CN[i] + sh_chilcotin_spawners_coef * spk)
        } else {
          spk * exp(df$sh_chilcotin_SSL_alpha[i] + sh_chilcotin_spawners_coef * spk)
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df$sh_chilcotin_recruits_alt[j]
      }
      
      if (yr < start_year_chilcotin + 4) {
        df$sh_chilcotin_Nage4_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p4[i]
      } else {
        df$sh_chilcotin_Nage4_pred[i] <- lag_recruits(yr - 4) * df$sh_chilcotin_p4[i] * 1000
      }
      
      if (yr < start_year_chilcotin + 5) {
        df$sh_chilcotin_Nage5_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p5[i]
      } else {
        df$sh_chilcotin_Nage5_pred[i] <- lag_recruits(yr - 5) * df$sh_chilcotin_p5[i] * 1000
      }
      
      if (yr < start_year_chilcotin + 6) {
        df$sh_chilcotin_Nage6_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p6[i]
      } else {
        df$sh_chilcotin_Nage6_pred[i] <- lag_recruits(yr - 6) * df$sh_chilcotin_p6[i] * 1000
      }
      
      if (yr < start_year_chilcotin + 7) {
        df$sh_chilcotin_Nage7_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p7[i]
      } else {
        df$sh_chilcotin_Nage7_pred[i] <- lag_recruits(yr - 7) * df$sh_chilcotin_p7[i] * 1000
      }
      
      if (yr < start_year_chilcotin + 8) {
        df$sh_chilcotin_Nage8_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p8[i]
      } else {
        df$sh_chilcotin_Nage8_pred[i] <- lag_recruits(yr - 8) * df$sh_chilcotin_p8[i] * 1000
      }
      
      df$sh_chilcotin_sum_pred[i] <-
        sum(df$sh_chilcotin_Nage4_pred[i],
            df$sh_chilcotin_Nage5_pred[i],
            df$sh_chilcotin_Nage6_pred[i],
            df$sh_chilcotin_Nage7_pred[i],
            df$sh_chilcotin_Nage8_pred[i],
            na.rm = TRUE)
      
      df$sh_chilcotin_bycatch_pred[i] <- df$sh_chilcotin_sum_pred[i] * df$sh_chilcotin_U_comm[i]
      
      if (yr <= 2018) {
        df$sh_chilcotin_FN_catch_pred[i] <- df$sh_chilcotin_FN_mortalities[i]
      } else {
        denom_2018 <- (df$sh_chilcotin_sum_pred[df$Year == 2018] -
                         df$sh_chilcotin_bycatch_pred[df$Year == 2018])
        df$sh_chilcotin_FN_catch_pred[i] <-
          FN_chilcotin_2018 / denom_2018 *
          (df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_bycatch_pred[i])
      }
      
      df$sh_chilcotin_total_catch_pred[i] <-
        df$sh_chilcotin_FN_catch_pred[i] +
        df$sh_chilcotin_sport_mortalities[i] +
        df$sh_chilcotin_bycatch_pred[i]
      
      S_new <- df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_total_catch_pred[i]
      S_new <- max(S_new, 0)
      
      if (is.finite(S_old) && is.finite(S_new) &&
          abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    
    df$sh_chilcotin_spawners_pred[i] <- S_old
  }
  
  # Attach scenario settings as columns (handy for exports)
  df %>%
    mutate(
      scenario_U_apply = U_apply,
      scenario_byrate = byrate,
      scenario_SSL_control = SSL_control,
      scenario_U_historic = U_historic
    )
}

# Map a "rate" (0-1) onto the alpha axis, and provide inverse for sec_axis
rate_to_alpha_transform <- function(alpha_min, alpha_max) {
  scale <- (alpha_max - alpha_min)
  list(
    to_alpha = function(rate) alpha_min + rate * scale,
    to_rate  = function(alpha_y) (alpha_y - alpha_min) / scale
  )
}


# add plot 
#left axis: alpha series (Chum base alpha; SH Thompson base/SSLC alpha; SH Chilcotin base/SSLC alpha)
#right axis: rates (Chum commercial harvest rate; SH commercial mortality rate)

plot_alpha_and_rates <- function(d, title = "Retrospective Chum and SH bycatch harvest rates and alpha values") {
  
  # pick alpha series you want present (only keep columns that exist)
  alpha_cols <- c(
    "chum_base_alpha",
    "sh_thompson_alpha_CN",
    "sh_thompson_SSL_alpha",
    "sh_chilcotin_alpha_CN",
    "sh_chilcotin_SSL_alpha"
  )
  alpha_cols <- alpha_cols[alpha_cols %in% names(d)]
  
  # pick rate series you want present (0-1)
  # chum commercial fishery harvest rate: use chum_commercial_harvest_uapply if you want the scenario rate,
  # or U_chum_pred if you want the model-implied rate.
  rate_cols <- c(
    "chum_commercial_harvest_uapply",
    "sh_thompson_U_comm",
    "sh_chilcotin_U_comm"
  )
  rate_cols <- rate_cols[rate_cols %in% names(d)]
  
  # alpha range for scaling (use observed alpha range, with a little padding)
  alpha_vals <- d %>% select(any_of(alpha_cols)) %>% unlist(use.names = FALSE)
  alpha_vals <- alpha_vals[is.finite(alpha_vals)]
  alpha_min <- if (length(alpha_vals)) min(alpha_vals) else -2
  alpha_max <- if (length(alpha_vals)) max(alpha_vals) else  3
  pad <- 0.05 * (alpha_max - alpha_min)
  alpha_min <- alpha_min - pad
  alpha_max <- alpha_max + pad
  
  tr <- rate_to_alpha_transform(alpha_min, alpha_max)
  
  alpha_long <- d %>%
    select(Year, any_of(alpha_cols)) %>%
    pivot_longer(-Year, names_to = "series", values_to = "value") %>%
    mutate(kind = "alpha")
  
  rate_long <- d %>%
    select(Year, any_of(rate_cols)) %>%
    pivot_longer(-Year, names_to = "series", values_to = "value") %>%
    mutate(
      kind = "rate",
      value = tr$to_alpha(pmin(pmax(value, 0), 1))  # clamp 0-1 then map to alpha axis
    )
  
  plot_df <- bind_rows(alpha_long, rate_long) %>%
    mutate(
      series = recode(series,
                      chum_base_alpha                = "Chum base alpha",
                      sh_thompson_alpha_CN           = "SH TH base alpha",
                      sh_thompson_SSL_alpha          = "SH TH SSLC alpha",
                      sh_chilcotin_alpha_CN          = "SH Chilc base alpha",
                      sh_chilcotin_SSL_alpha         = "SH Chilc SSLC alpha",
                      chum_commercial_harvest_uapply = "Chum commercial fishery harvest rate",
                      sh_thompson_U_comm             = "Thompson SH mortality rate in commercial fishery",
                      sh_chilcotin_U_comm            = "Chilcotin SH mortality rate in commercial fishery",
                      .default = series
      ),
      lty = case_when(
        kind == "alpha" ~ "alpha",
        TRUE ~ "rate"
      )
    )
  
  ggplot(plot_df, aes(x = Year, y = value, group = series)) +
    geom_line(aes(linetype = series), linewidth = 0.9, na.rm = TRUE) +
    scale_y_continuous(
      name = "Alpha",
      limits = c(alpha_min, alpha_max),
      sec.axis = sec_axis(~ tr$to_rate(.), name = "Harvest or fishing mortality rate")
    ) +
    labs(title = title, x = "Year", linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

# plot 
#left axis: chum spawners (observed points + predicted line)
#right axis: steelhead spawners (Thompson and Chilcotin; observed points + predicted lines)
plot_spawners_dual_axis <- function(
    d,
    title = "Predicted and Observed Chum and Thompson and Chilcotin River Steelhead Spawner Abundance",
    sh_pred_divisor = 1000  # <- key: predicted SH spawners are in fish; observed are in thousands
) {
  
  needed <- c("Year", "chum_spawners", "chum_spawners_pred",
              "sh_thompson_spawners", "sh_thompson_spawners_pred",
              "sh_chilcotin_spawners", "sh_chilcotin_spawners_pred")
  
  d <- d %>% select(any_of(needed)) %>% filter(!is.na(Year))
  
  # --- put SH predicted into same units as observed (thousands) ---
  if ("sh_thompson_spawners_pred" %in% names(d)) {
    d <- d %>% mutate(sh_thompson_spawners_pred_k = sh_thompson_spawners_pred / sh_pred_divisor)
  }
  if ("sh_chilcotin_spawners_pred" %in% names(d)) {
    d <- d %>% mutate(sh_chilcotin_spawners_pred_k = sh_chilcotin_spawners_pred / sh_pred_divisor)
  }
  
  # ranges for axis mapping
  chum_vals <- c(d$chum_spawners, d$chum_spawners_pred)
  chum_vals <- chum_vals[is.finite(chum_vals)]
  if (!length(chum_vals)) chum_vals <- c(0, 1)
  chum_lim <- c(0, max(chum_vals, na.rm = TRUE))
  
  # steelhead values: observed + (corrected) predicted
  sh_vals <- c(
    d$sh_thompson_spawners,
    d$sh_thompson_spawners_pred_k,
    d$sh_chilcotin_spawners,
    d$sh_chilcotin_spawners_pred_k
  )
  sh_vals <- sh_vals[is.finite(sh_vals)]
  if (!length(sh_vals)) sh_vals <- c(0, 1)
  sh_lim <- c(0, max(sh_vals, na.rm = TRUE))
  
  # map steelhead (right) -> chum scale (left)
  a <- diff(chum_lim) / diff(sh_lim)
  b <- chum_lim[1] - a * sh_lim[1]
  to_left  <- function(x) a * x + b
  to_right <- function(y) (y - b) / a
  
  p <- ggplot(d, aes(x = Year)) +
    labs(title = title, x = "Year") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # ---- CHUM (left axis) ----
  if ("chum_spawners" %in% names(d)) {
    p <- p + geom_point(aes(y = chum_spawners, shape = "Observed Chum spawners"),
                        size = 2.6, stroke = 1.1, na.rm = TRUE, color = "firebrick")
  }
  if ("chum_spawners_pred" %in% names(d)) {
    p <- p + geom_line(aes(y = chum_spawners_pred, linetype = "Predicted Chum spawners"),
                       linewidth = 1.0, na.rm = TRUE, color = "firebrick")
  }
  
  # ---- STEELHEAD (right axis, scaled onto left) ----
  if ("sh_thompson_spawners" %in% names(d)) {
    p <- p + geom_point(aes(y = to_left(sh_thompson_spawners), shape = "Observed Thompson spawners"),
                        size = 2.6, stroke = 1.1, na.rm = TRUE, color = "purple")
  }
  if ("sh_thompson_spawners_pred_k" %in% names(d)) {
    p <- p + geom_line(aes(y = to_left(sh_thompson_spawners_pred_k), linetype = "Predicted Thompson spawners"),
                       linewidth = 1.0, na.rm = TRUE, color = "purple")
  }
  
  if ("sh_chilcotin_spawners" %in% names(d)) {
    p <- p + geom_point(aes(y = to_left(sh_chilcotin_spawners), shape = "Observed Chilcotin spawners"),
                        size = 2.6, stroke = 1.1, na.rm = TRUE, color = "darkgreen")
  }
  if ("sh_chilcotin_spawners_pred_k" %in% names(d)) {
    p <- p + geom_line(aes(y = to_left(sh_chilcotin_spawners_pred_k), linetype = "Predicted Chilcotin spawners"),
                       linewidth = 1.0, na.rm = TRUE, color = "darkgreen")
  }
  
  p +
    scale_y_continuous(
      name = "Chum Salmon Spawners",
      limits = chum_lim,
      sec.axis = sec_axis(~ to_right(.), name = "Steelhead spawner abundance")
    ) +
    labs(shape = NULL, linetype = NULL)
}




#-------------------------------
# Shiny UI
#-------------------------------
ui <- fluidPage(
  titlePanel("Retrospective Model Scenario Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("U_apply", "U_apply", value = 0.2, min = 0, max = 1, step = 0.01),
      numericInput("byrate", "byrate", value = 0.69, min = 0, max = 2, step = 0.01),
      selectInput("SSL_control", "SSL_control", choices = c("0" = 0, "1" = 1), selected = 0),
      selectInput("U_historic", "U_historic", choices = c("0" = 0, "1" = 1), selected = 1),
      
      tags$hr(),
      
      actionButton("run", "Run scenario", class = "btn-primary"),
      downloadButton("download", "Download results CSV")
    ),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Chum",
          plotOutput("plot_chum_spawners", height = 280),
          plotOutput("plot_chum_recruits", height = 280)
        ),
        tabPanel(
          "Steelhead (Thompson)",
          plotOutput("plot_th_spawners", height = 280),
          plotOutput("plot_th_ucomm", height = 280)
        ),
        tabPanel(
          "Steelhead (Chilcotin)",
          plotOutput("plot_ch_spawners", height = 280),
          plotOutput("plot_ch_ucomm", height = 280)
        ),
        tabPanel("Spawner Abundance", plotOutput("plot_spawners_dual", height = 520)),
        tabPanel("Alpha + Rates", plotOutput("plot_alpha_rates", height = 520)),
        tabPanel(
          "Data",
          helpText("Showing a subset of columns. Download CSV for everything."),
          tableOutput("table_head")
        )
      )
    )
  )
)



#-------------------------------
# Shiny server
#-------------------------------
server <- function(input, output, session) {
  
  # run on button click (prevents re-running on every tiny input change)
  results <- eventReactive(input$run, {
    run_model(
      U_apply = input$U_apply,
      byrate = input$byrate,
      SSL_control = as.integer(input$SSL_control),
      U_historic  = as.integer(input$U_historic)
    )
  }, ignoreInit = TRUE)
  
  # convenient reactive for plots
  dfp <- reactive({
    req(results())
    results()
  })
  
  output$plot_chum_spawners <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year)) +
      geom_line(aes(y = chum_spawners, linetype = "Observed")) +
      geom_line(aes(y = chum_spawners_pred, linetype = "Predicted")) +
      labs(y = "Chum spawners", linetype = NULL) +
      theme_minimal()
  })
  
  output$plot_chum_recruits <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year)) +
      geom_line(aes(y = chum_recruits_obs, linetype = "Observed")) +
      geom_line(aes(y = chum_recruits_alt, linetype = "Scenario recruits")) +
      labs(y = "Chum recruits", linetype = NULL) +
      theme_minimal()
  })
  
  output$plot_th_spawners <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year)) +
      geom_line(aes(y = sh_thompson_spawners, linetype = "Observed")) +
      geom_line(aes(y = sh_thompson_spawners_pred, linetype = "Predicted")) +
      labs(y = "Thompson steelhead spawners", linetype = NULL) +
      theme_minimal()
  })
  
  output$plot_th_ucomm <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year, y = sh_thompson_U_comm)) +
      geom_line() +
      labs(y = "Thompson U_comm") +
      theme_minimal()
  })
  
  output$plot_ch_spawners <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year)) +
      geom_line(aes(y = sh_chilcotin_spawners, linetype = "Observed")) +
      geom_line(aes(y = sh_chilcotin_spawners_pred, linetype = "Predicted")) +
      labs(y = "Chilcotin steelhead spawners", linetype = NULL) +
      theme_minimal()
  })
  
  output$plot_ch_ucomm <- renderPlot({
    d <- dfp()
    ggplot(d, aes(x = Year, y = sh_chilcotin_U_comm)) +
      geom_line() +
      labs(y = "Chilcotin U_comm") +
      theme_minimal()
  })
  
  output$plot_alpha_rates <- renderPlot({
    req(dfp())
    plot_alpha_and_rates(dfp())
  }, res = 120)
  
  output$plot_spawners_dual <- renderPlot({
    req(dfp())
    plot_spawners_dual_axis(dfp())
  }, res = 120)
  
  
  output$table_head <- renderTable({
    d <- dfp()
    d %>%
      select(
        Year,
        chum_spawners, chum_spawners_pred, chum_recruits_obs, chum_recruits_alt,
        sh_thompson_spawners, sh_thompson_spawners_pred, sh_thompson_U_comm,
        sh_chilcotin_spawners, sh_chilcotin_spawners_pred, sh_chilcotin_U_comm,
        scenario_U_apply, scenario_byrate, scenario_SSL_control, scenario_U_historic
      ) %>%
      tail(30)
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0(
        "scenario_",
        "Uapply", input$U_apply, "_",
        "byrate", input$byrate, "_",
        "SSL", input$SSL_control, "_",
        "Uhist", input$U_historic,
        ".csv"
      )
    },
    content = function(file) {
      readr::write_csv(dfp(), file)
    }
  )
}

shinyApp(ui, server)
