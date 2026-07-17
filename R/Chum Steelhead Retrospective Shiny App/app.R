# Retrospective Model - Shiny App
# Based on Haley Oleynik Murdoch McAllister's chum/steelhead retrospective mode
# Toggle U_apply, bycatch_rate, and SSL_control reactively

library(shiny)
library(tidyverse)
library(slider)

# ---- Fixed coefficients (not exposed as controls) ----------------------
intercept          <- 1.03737862843252
pdo_adult_coef     <- 0.093
npgo_coef          <- 0.103
pdo_smolt_coef     <- -0.106
ssl_coef           <- -0.224
spawners_coef      <- -4.95E-07

sh_thompson_intercept     <- 1.572107637
sh_thompson_sst_coef      <- -0.203463091
sh_thompson_ssl_coef      <- -0.764277677
sh_thompson_npgo_coef     <- -0.0402
sh_thompson_spawners_coef <- -0.804438793

sh_chilcotin_intercept     <- 1.053608979
sh_chilcotin_sst_coef      <- -0.127949278
sh_chilcotin_ssl_coef      <- -0.792741195
sh_chilcotin_npgo_coef     <- 0.152526045
sh_chilcotin_pdo_coef      <- 0.202708011
sh_chilcotin_spawners_coef <- -1.022467631


# ---- Data load -----------------------------------------------------------
data_raw       <- read_csv("s-r_data.csv", show_col_types = FALSE) %>%
  select(-SSL)  
sh_data        <- read_csv("sh_s-r_data.csv", show_col_types = FALSE)
covariates_raw <- read_csv("covariates.csv", show_col_types = FALSE) 

# ---- Core model function --------------------------------------------------
# Wraps the full chum -> steelhead (Thompson + Chilcotin) pipeline.
# U_apply, bycatch_rate, SSL_control are the reactive controls.
run_model <- function(U_apply, bycatch_rate, SSL_control, U_historic) {
  
  byrate <- bycatch_rate
  
  ## ---- CHUM -----------------------------------------------------------
  new.data <- data_raw %>%
    left_join(covariates_raw, by = "Year") %>%
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
      recruits_pred = if_else(
        rowSums(!is.na(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)))) == 0,
        NA_real_,
        rowSums(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)), na.rm = TRUE)
      ),
      recruits_dif = chum_total_stock - recruits_pred,
      U_chum_pred = catch / recruits_pred,
      U_chum_dif = U_chum_pred - U_chum,
      chum_commercial_harvest = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 ~ U_chum),
      chum_commercial_harvest_uapply = case_when(
        Year <= 1990 ~ U_chum_pred,
        Year >= 1991 ~ U_apply)
    )
  
  SSL_1978 <- new.data %>% filter(Year == 1978) %>% pull(SSL)
  
  df <- new.data %>%
    arrange(Year) %>%
    mutate(
      chum_recruits_alt  = NA_real_,
      Nage3_alt = Nage3_pred,
      Nage4_alt = Nage4_pred,
      Nage5_alt = Nage5_pred,
      Nage6_alt = Nage6_pred,
      sum_alt   = NA_real_,
      catch_alt = NA_real_,
      chum_spawners_pred = chum_spawners
    )
  
  for (i in seq_len(nrow(df))) {
    
    SSL_alt <- if (df$Year[i] <= 1978) {
      df$SSL[i]
    } else {
      (1 - SSL_control) * df$SSL[i] + SSL_control * SSL_1978
    }
    
    SSL_control_alpha <-
      intercept +
      df$PDO_adult[i] * pdo_adult_coef +
      df$NPGO[i] * npgo_coef +
      df$PDO_smolt[i] * pdo_smolt_coef +
      SSL_alt * ssl_coef
    
    df$chum_recruits_alt[i] <-
      if (SSL_control == 0) {
        df$chum_spawners_pred[i] *
          exp(df$chum_base_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
          exp(df$chum_ln_obs_pred[i])
      } else {
        df$chum_spawners_pred[i] *
          exp(SSL_control_alpha + spawners_coef * df$chum_spawners_pred[i]) *
          exp(df$chum_ln_obs_pred[i])
      }
    
    if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
    if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
    if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
    if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
    
    df$sum_alt[i] <- if (all(is.na(c(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i]))))
      NA_real_ else
        sum(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i], na.rm = TRUE)
    df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
    df$chum_spawners_pred[i] <- max(df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i]), 0)
  }
  
  ## ---- Bind chum to steelhead ------------------------------------------
  all_data <- df %>% right_join(sh_data, by = "Year")
  
  ## ---- STEELHEAD: Thompson ---------------------------------------------
  start_year_t <- 1978
  
  sh_thompson_SSL_1978 <- all_data %>% filter(Year == 1978) %>% pull(sh_thompson_SL) %>% as.numeric()
  FN_thompson_2018 <- all_data %>% filter(Year == 2018) %>% pull(sh_thompson_FN_mortalities) %>% as.numeric()
  
  df_t <- all_data %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_thompson_base_alpha = sh_thompson_intercept +
        sh_thompson_SST * sh_thompson_sst_coef +
        sh_thompson_SL * sh_thompson_ssl_coef +
        sh_thompson_NPGO * sh_thompson_npgo_coef,
      sh_thompson_model_recruits = (sh_thompson_spawners / 1000) *
        exp(sh_thompson_base_alpha + sh_thompson_spawners_coef * (sh_thompson_spawners / 1000)),
      sh_thompson_ln_obs_pred = log(sh_thompson_recruits / sh_thompson_model_recruits),
      sh_thompson_pred_bycatch = sh_thompson_prefishery_N - sh_thompson_sport_mortalities -
        sh_thompson_FN_mortalities - 1000 * sh_thompson_spawners,
      sh_thompson_U = sh_thompson_pred_bycatch / sh_thompson_prefishery_N,
      sh_thompson_recruits_alt  = NA_real_,
      sh_thompson_spawners_pred = sh_thompson_spawners,
      sh_thompson_Nage4_pred = NA_real_,
      sh_thompson_Nage5_pred = NA_real_,
      sh_thompson_Nage6_pred = NA_real_,
      sh_thompson_Nage7_pred = NA_real_,
      sh_thompson_Nage8_pred = NA_real_,
      sh_thompson_SSL_alt          = NA_real_,
      sh_thompson_SSL_alpha        = NA_real_,
      sh_thompson_alpha_CN         = NA_real_,
      sh_thompson_sum_pred         = NA_real_,
      sh_thompson_bycatch_pred     = NA_real_,
      sh_thompson_FN_catch_pred    = NA_real_,
      sh_thompson_total_catch_pred = NA_real_,
      sh_thompson_U_comm           = NA_real_
    )
  
  start_i <- which(df_t$Year >= start_year_t)[1]
  year_to_i <- setNames(seq_len(nrow(df_t)), df_t$Year)
  
  for (i in seq(from = start_i, to = nrow(df_t))) {
    yr <- df_t$Year[i]
    
    df_t$sh_thompson_U_comm[i] <-
      if (yr <= 1990) {
        df_t$sh_thompson_U[i]
      } else if (U_historic == 1) {
        df_t$sh_thompson_U[i]
      } else {
        byrate * df_t$chum_commercial_harvest_uapply[i]
      }
    
    df_t$sh_thompson_SSL_alt[i] <- if (yr <= 1978) {
      df_t$sh_thompson_SL[i]
    } else {
      (1 - SSL_control) * df_t$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
    }
    
    df_t$sh_thompson_SSL_alpha[i] <-
      sh_thompson_intercept +
      df_t$sh_thompson_SST[i] * sh_thompson_sst_coef +
      df_t$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
      df_t$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
    
    df_t$sh_thompson_alpha_CN[i] <-
      sh_thompson_intercept +
      df_t$sh_thompson_SST[i] * sh_thompson_sst_coef +
      df_t$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
      df_t$sh_thompson_SL[i] * sh_thompson_ssl_coef
    
    S_old <- df_t$sh_thompson_spawners_pred[i]
    if (is.na(S_old)) S_old <- df_t$sh_thompson_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    max_iter <- 50; tol <- 1e-8
    
    for (iter in seq_len(max_iter)) {
      df_t$sh_thompson_spawners_pred[i] <- S_old
      spk <- df_t$sh_thompson_spawners_pred[i]
      
      df_t$sh_thompson_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df_t$sh_thompson_alpha_CN[i] + sh_thompson_spawners_coef * spk)
        } else {
          spk * exp(df_t$sh_thompson_SSL_alpha[i] + sh_thompson_spawners_coef * spk)
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df_t$sh_thompson_recruits_alt[j]
      }
      
      df_t$sh_thompson_Nage4_pred[i] <- if (yr < start_year_t + 4) df_t$sh_thompson_prefishery_N[i] * df_t$sh_thompson_p4[i] else lag_recruits(yr - 4) * df_t$sh_thompson_p4[i] * 1000
      df_t$sh_thompson_Nage5_pred[i] <- if (yr < start_year_t + 5) df_t$sh_thompson_prefishery_N[i] * df_t$sh_thompson_p5[i] else lag_recruits(yr - 5) * df_t$sh_thompson_p5[i] * 1000
      df_t$sh_thompson_Nage6_pred[i] <- if (yr < start_year_t + 6) df_t$sh_thompson_prefishery_N[i] * df_t$sh_thompson_p6[i] else lag_recruits(yr - 6) * df_t$sh_thompson_p6[i] * 1000
      df_t$sh_thompson_Nage7_pred[i] <- if (yr < start_year_t + 7) df_t$sh_thompson_prefishery_N[i] * df_t$sh_thompson_p7[i] else lag_recruits(yr - 7) * df_t$sh_thompson_p7[i] * 1000
      df_t$sh_thompson_Nage8_pred[i] <- if (yr < start_year_t + 8) df_t$sh_thompson_prefishery_N[i] * df_t$sh_thompson_p8[i] else lag_recruits(yr - 8) * df_t$sh_thompson_p8[i] * 1000
      
      df_t$sh_thompson_sum_pred[i] <- sum(
        df_t$sh_thompson_Nage4_pred[i], df_t$sh_thompson_Nage5_pred[i],
        df_t$sh_thompson_Nage6_pred[i], df_t$sh_thompson_Nage7_pred[i],
        df_t$sh_thompson_Nage8_pred[i], na.rm = TRUE)
      
      df_t$sh_thompson_bycatch_pred[i] <- df_t$sh_thompson_sum_pred[i] * df_t$sh_thompson_U_comm[i]
      
      if (yr <= 2018) {
        df_t$sh_thompson_FN_catch_pred[i] <- df_t$sh_thompson_FN_mortalities[i]
      } else {
        denom_2018 <- (df_t$sh_thompson_sum_pred[df_t$Year == 2018] - df_t$sh_thompson_bycatch_pred[df_t$Year == 2018])
        df_t$sh_thompson_FN_catch_pred[i] <- FN_thompson_2018 / denom_2018 *
          (df_t$sh_thompson_sum_pred[i] - df_t$sh_thompson_bycatch_pred[i])
      }
      
      df_t$sh_thompson_total_catch_pred[i] <-
        df_t$sh_thompson_FN_catch_pred[i] + df_t$sh_thompson_sport_mortalities[i] + df_t$sh_thompson_bycatch_pred[i]
      
      S_new <- max((df_t$sh_thompson_sum_pred[i] - df_t$sh_thompson_total_catch_pred[i]) / 1000, 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    df_t$sh_thompson_spawners_pred[i] <- S_old
  }
  
  ## ---- STEELHEAD: Chilcotin --------------------------------------------
  start_year_c <- 1973
  
  sh_chilcotin_SSL_1978 <- df_t %>% filter(Year == 1978) %>% pull(sh_chilcotin_SL) %>% as.numeric()
  FN_chilcotin_2018 <- df_t %>% filter(Year == 2018) %>% pull(sh_chilcotin_FN_mortalities) %>% as.numeric()
  
  df_c <- df_t %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_chilcotin_base_alpha = sh_chilcotin_intercept +
        sh_chilcotin_SST * sh_chilcotin_sst_coef +
        sh_chilcotin_SL * sh_chilcotin_ssl_coef +
        sh_chilcotin_NPGO * sh_chilcotin_npgo_coef +
        sh_chilcotin_PDO * sh_chilcotin_pdo_coef,
      sh_chilcotin_model_recruits = (sh_chilcotin_spawners / 1000) *
        exp(sh_chilcotin_base_alpha + sh_chilcotin_spawners_coef * (sh_chilcotin_spawners / 1000)),
      sh_chilcotin_ln_obs_pred = log(sh_chilcotin_recruits / sh_chilcotin_model_recruits),
      sh_chilcotin_pred_bycatch = sh_chilcotin_prefishery_N - sh_chilcotin_sport_mortalities -
        sh_chilcotin_FN_mortalities - 1000 * sh_chilcotin_spawners,
      sh_chilcotin_U = sh_chilcotin_pred_bycatch / sh_chilcotin_prefishery_N,
      sh_chilcotin_recruits_alt  = NA_real_,
      sh_chilcotin_spawners_pred = sh_chilcotin_spawners,
      sh_chilcotin_Nage4_pred = NA_real_,
      sh_chilcotin_Nage5_pred = NA_real_,
      sh_chilcotin_Nage6_pred = NA_real_,
      sh_chilcotin_Nage7_pred = NA_real_,
      sh_chilcotin_Nage8_pred = NA_real_,
      sh_chilcotin_SSL_alt          = NA_real_,
      sh_chilcotin_SSL_alpha        = NA_real_,
      sh_chilcotin_alpha_CN         = NA_real_,
      sh_chilcotin_sum_pred         = NA_real_,
      sh_chilcotin_bycatch_pred     = NA_real_,
      sh_chilcotin_FN_catch_pred    = NA_real_,
      sh_chilcotin_total_catch_pred = NA_real_,
      sh_chilcotin_U_comm           = NA_real_
    )
  
  start_i_c <- which(df_c$Year >= start_year_c)[1]
  year_to_i_c <- setNames(seq_len(nrow(df_c)), df_c$Year)
  
  for (i in seq(from = start_i_c, to = nrow(df_c))) {
    yr <- df_c$Year[i]
    
    df_c$sh_chilcotin_U_comm[i] <-
      if (yr <= 1990) {
        df_c$sh_chilcotin_U[i]
      } else if (U_historic == 1) {
        df_c$sh_chilcotin_U[i]
      } else {
        byrate * df_c$chum_commercial_harvest_uapply[i]
      }
    
    df_c$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) {
      df_c$sh_chilcotin_SL[i]
    } else {
      (1 - SSL_control) * df_c$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
    }
    
    df_c$sh_chilcotin_SSL_alpha[i] <-
      sh_chilcotin_intercept +
      df_c$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef +
      df_c$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df_c$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef +
      df_c$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
    
    df_c$sh_chilcotin_alpha_CN[i] <-
      sh_chilcotin_intercept +
      df_c$sh_chilcotin_SST[i] * sh_chilcotin_sst_coef +
      df_c$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
      df_c$sh_chilcotin_SL[i] * sh_chilcotin_ssl_coef +
      df_c$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef
    
    S_old <- df_c$sh_chilcotin_spawners_pred[i]
    if (is.na(S_old)) S_old <- df_c$sh_chilcotin_spawners[i]
    if (is.na(S_old)) S_old <- 0
    
    max_iter <- 50; tol <- 1e-8
    
    for (iter in seq_len(max_iter)) {
      df_c$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df_c$sh_chilcotin_spawners_pred[i]
      
      df_c$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df_c$sh_chilcotin_alpha_CN[i] + sh_chilcotin_spawners_coef * spk)
        } else {
          spk * exp(df_c$sh_chilcotin_SSL_alpha[i] + sh_chilcotin_spawners_coef * spk)
        }
      
      lag_recruits <- function(lag_year) {
        j <- year_to_i_c[as.character(lag_year)]
        if (is.na(j)) NA_real_ else df_c$sh_chilcotin_recruits_alt[j]
      }
      
      df_c$sh_chilcotin_Nage4_pred[i] <- if (yr < start_year_c + 4) df_c$sh_chilcotin_prefishery_N[i] * df_c$sh_chilcotin_p4[i] else lag_recruits(yr - 4) * df_c$sh_chilcotin_p4[i] * 1000
      df_c$sh_chilcotin_Nage5_pred[i] <- if (yr < start_year_c + 5) df_c$sh_chilcotin_prefishery_N[i] * df_c$sh_chilcotin_p5[i] else lag_recruits(yr - 5) * df_c$sh_chilcotin_p5[i] * 1000
      df_c$sh_chilcotin_Nage6_pred[i] <- if (yr < start_year_c + 6) df_c$sh_chilcotin_prefishery_N[i] * df_c$sh_chilcotin_p6[i] else lag_recruits(yr - 6) * df_c$sh_chilcotin_p6[i] * 1000
      df_c$sh_chilcotin_Nage7_pred[i] <- if (yr < start_year_c + 7) df_c$sh_chilcotin_prefishery_N[i] * df_c$sh_chilcotin_p7[i] else lag_recruits(yr - 7) * df_c$sh_chilcotin_p7[i] * 1000
      df_c$sh_chilcotin_Nage8_pred[i] <- if (yr < start_year_c + 8) df_c$sh_chilcotin_prefishery_N[i] * df_c$sh_chilcotin_p8[i] else lag_recruits(yr - 8) * df_c$sh_chilcotin_p8[i] * 1000
      
      df_c$sh_chilcotin_sum_pred[i] <- sum(
        df_c$sh_chilcotin_Nage4_pred[i], df_c$sh_chilcotin_Nage5_pred[i],
        df_c$sh_chilcotin_Nage6_pred[i], df_c$sh_chilcotin_Nage7_pred[i],
        df_c$sh_chilcotin_Nage8_pred[i], na.rm = TRUE)
      
      df_c$sh_chilcotin_bycatch_pred[i] <- df_c$sh_chilcotin_sum_pred[i] * df_c$sh_chilcotin_U_comm[i]
      
      if (yr <= 2018) {
        df_c$sh_chilcotin_FN_catch_pred[i] <- df_c$sh_chilcotin_FN_mortalities[i]
      } else {
        denom_2018 <- (df_c$sh_chilcotin_sum_pred[df_c$Year == 2018] - df_c$sh_chilcotin_bycatch_pred[df_c$Year == 2018])
        df_c$sh_chilcotin_FN_catch_pred[i] <- FN_chilcotin_2018 / denom_2018 *
          (df_c$sh_chilcotin_sum_pred[i] - df_c$sh_chilcotin_bycatch_pred[i])
      }
      
      df_c$sh_chilcotin_total_catch_pred[i] <-
        df_c$sh_chilcotin_FN_catch_pred[i] + df_c$sh_chilcotin_sport_mortalities[i] + df_c$sh_chilcotin_bycatch_pred[i]
      
      S_new <- max((df_c$sh_chilcotin_sum_pred[i] - df_c$sh_chilcotin_total_catch_pred[i]) / 1000, 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      S_old <- S_new
    }
    df_c$sh_chilcotin_spawners_pred[i] <- S_old
  }
  
  list(chum = df, thompson = df_t, chilcotin = df_c)
}

# ---- UI --------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Fraser Chum & Steelhead Retrospective Model"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("U_apply", "U_apply (post-1990 chum exploitation rate)",
                  min = 0, max = 1, value = 0.2, step = 0.01),
      sliderInput("bycatch_rate", "bycatch_rate (steelhead bycatch scalar)",
                  min = 0, max = 2, value = 0.67, step = 0.01),
      checkboxInput("SSL_control", "Hold SSL at 1978 level (SSL_control)",
                    value = FALSE),
      checkboxInput("U_historic", "Use historic steelhead U (U_historic)",
                    value = TRUE),
      hr(),
      helpText("When 'Use historic steelhead U' is unchecked, steelhead U_comm switches to ",
               "bycatch_rate * chum_commercial_harvest_uapply for years after 1990.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Chum",
                 plotOutput("chumPlot"),
                 downloadButton("dlChum", "Download chum results (.csv)")#,
                 #tableOutput("chumTable")
        ),
        tabPanel("Steelhead - Thompson",
                 plotOutput("thompsonPlot"),
                 downloadButton("dlThompson", "Download Thompson results (.csv)")#,
                 #tableOutput("thompsonTable")
        ),
        tabPanel("Steelhead - Chilcotin",
                 plotOutput("chilcotinPlot"),
                 downloadButton("dlChilcotin", "Download Chilcotin results (.csv)")#,
                 #tableOutput("chilcotinTable")
        )
      )
    )
  )
)

# ---- Server ------------------------------------------------------------
server <- function(input, output, session) {
  
  model_out <- reactive({
    run_model(
      U_apply      = input$U_apply,
      bycatch_rate = input$bycatch_rate,
      SSL_control  = as.numeric(input$SSL_control),
      U_historic   = as.numeric(input$U_historic)
    )
  })
  
  output$chumPlot <- renderPlot({
    d <- model_out()$chum
    d %>%
      select(Year, chum_spawners, chum_spawners_pred, chum_recruits_obs, chum_recruits_alt) %>%
      pivot_longer(-Year, names_to = "series", values_to = "value") %>%
      mutate(
        series = factor(series, levels = c("chum_spawners", "chum_spawners_pred",
                                           "chum_recruits_obs", "chum_recruits_alt"))
      ) %>%
      ggplot(aes(Year, value, color = series, linetype = series)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = c(
        chum_spawners      = "#24492e",
        chum_spawners_pred = "#015b58",
        chum_recruits_obs  = "#e69b99",
        chum_recruits_alt  = "#ba7999"
      )) +
      scale_linetype_manual(values = c(
        chum_spawners      = "dashed",
        chum_spawners_pred = "solid",
        chum_recruits_obs  = "dashed",
        chum_recruits_alt  = "solid"
      )) +
      labs(title = "Chum: observed vs. modeled spawners & recruits",
           y = "Fish", color = NULL, linetype = NULL) +
      theme_minimal()
  })
  
  #output$chumTable <- renderTable({
  #  model_out()$chum %>%
  #    select(Year, chum_spawners, chum_spawners_pred, chum_recruits_obs,
  #           chum_recruits_alt, U_chum, U_chum_pred, chum_commercial_harvest_uapply) %>%
  #    tail(15)
  #})
  
  output$dlChum <- downloadHandler(
    filename = function() "chum_results.csv",
    content = function(file) write_csv(model_out()$chum, file)
  )
  
  output$thompsonPlot <- renderPlot({
    d <- model_out()$thompson
    d %>%
      select(Year, sh_thompson_spawners, sh_thompson_spawners_pred,
             sh_thompson_recruits, sh_thompson_recruits_alt) %>%
      pivot_longer(-Year, names_to = "series", values_to = "value") %>%
      mutate(
        series = factor(series, levels = c("sh_thompson_spawners", "sh_thompson_spawners_pred",
                                           "sh_thompson_recruits", "sh_thompson_recruits_alt"))
      ) %>%
      ggplot(aes(Year, value, color = series, linetype = series)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = c(
        sh_thompson_spawners      = "#24492e",
        sh_thompson_spawners_pred = "#015b58",
        sh_thompson_recruits      = "#e69b99",
        sh_thompson_recruits_alt  = "#ba7999"
      )) +
      scale_linetype_manual(values = c(
        sh_thompson_spawners      = "dashed",
        sh_thompson_spawners_pred = "solid",
        sh_thompson_recruits      = "dashed",
        sh_thompson_recruits_alt  = "solid"
      )) +
      labs(title = "Steelhead (Thompson): observed vs. modeled spawners & recruits",
           y = "Fish", color = NULL, linetype = NULL) +
      theme_minimal()
  })
  
  #output$thompsonTable <- renderTable({
  #  model_out()$thompson %>%
  #    select(Year, sh_thompson_spawners, sh_thompson_spawners_pred,
  #           sh_thompson_U_comm, sh_thompson_total_catch_pred) %>%
  #    tail(15)
  #})
  
  output$dlThompson <- downloadHandler(
    filename = function() "thompson_steelhead_results.csv",
    content = function(file) write_csv(model_out()$thompson, file)
  )
  
  output$chilcotinPlot <- renderPlot({
    d <- model_out()$chilcotin
    d %>%
      select(Year, sh_chilcotin_spawners, sh_chilcotin_spawners_pred,
             sh_chilcotin_recruits, sh_chilcotin_recruits_alt) %>%
      pivot_longer(-Year, names_to = "series", values_to = "value") %>%
      mutate(
        series = factor(series, levels = c("sh_chilcotin_spawners", "sh_chilcotin_spawners_pred",
                                           "sh_chilcotin_recruits", "sh_chilcotin_recruits_alt"))
      ) %>%
      ggplot(aes(Year, value, color = series, linetype = series)) +
      geom_line(linewidth = 0.8) +
      scale_color_manual(values = c(
        sh_chilcotin_spawners      = "#24492e",
        sh_chilcotin_spawners_pred = "#015b58",
        sh_chilcotin_recruits      = "#e69b99",
        sh_chilcotin_recruits_alt  = "#ba7999"
      )) +
      scale_linetype_manual(values = c(
        sh_chilcotin_spawners      = "dashed",
        sh_chilcotin_spawners_pred = "solid",
        sh_chilcotin_recruits      = "dashed",
        sh_chilcotin_recruits_alt  = "solid"
      )) +
      labs(title = "Steelhead (Chilcotin): observed vs. modeled spawners & recruits",
           y = "Fish", color = NULL, linetype = NULL) +
      theme_minimal()
  })
  
  #output$chilcotinTable <- renderTable({
  #  model_out()$chilcotin %>%
  #    select(Year, sh_chilcotin_spawners, sh_chilcotin_spawners_pred,
  #           sh_chilcotin_U_comm, sh_chilcotin_total_catch_pred) %>%
  #    tail(15)
  #})
  
  output$dlChilcotin <- downloadHandler(
    filename = function() "chilcotin_steelhead_results.csv",
    content = function(file) write_csv(model_out()$chilcotin, file)
  )
}

shinyApp(ui, server)