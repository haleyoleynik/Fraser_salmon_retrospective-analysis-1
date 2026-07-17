# Retrospective Model 
# Haley Oleynik Murdoch McAllister
# October 2025

# load libraries 
require(tidyverse)
require(ggplot2)
require(readr)
require(slider)
require(zoo)

# read data ------------------
data <- read_csv("R/Chum Steelhead Retrospective Shiny App/s-r_data.csv") %>%
  select(-SSL)
sh_data <- read_csv("R/Chum Steelhead Retrospective Shiny App/sh_s-r_data.csv")
covariates <- read_csv("R/Chum Steelhead Retrospective Shiny App/covariates.csv") 
# Scenario 1 --------------------------------------
# no predator controls 
## CHUM ----------------
# from Fraser Chum data v19_alpha_SSL_est_fin_yrs_v8 & v9 sheet 

# estimate the coefficients 
# calculate alpha with coefficients using this formula: 

# coefficients (from lnrs model)
# can estimate these directly in r and then input (to do)
intercept = 1.03737862843252
pdo_adult_coef = 0.0929480015696915
npgo_coef = 0.102617626088713
pdo_smolt_coef = -0.105950753725792
ssl_coef = -0.224246284912103
spawners_coef = -4.95136622626478E-07

# control settings to change 
U_apply = 0.2
SSL_control = 0 
U_historic  <- 1
byrate      <- 0.69
start_year  <- 1978

# mutate dataframe with calculations from exel 
new.data <- data %>%
  left_join(covariates, by = "Year") %>%
  arrange(Year) %>%  # ensure chronological order
  mutate(
    catch = chum_total_stock - chum_spawners,
    U_chum = catch / chum_total_stock,
    chum_base_alpha = intercept + (PDO_adult * pdo_adult_coef + NPGO * npgo_coef +
      PDO_smolt * pdo_smolt_coef + SSL * ssl_coef),
    alpha_running_avg = slide_dbl(chum_base_alpha, mean, .before = 9, .complete = TRUE),
    chum_model_recruits = chum_spawners * exp(chum_base_alpha + chum_spawners * spawners_coef),
    chum_ln_obs_pred = log(chum_recruits_obs / chum_model_recruits),
    Nage3_obs = chum_total_stock * prop3,
    Nage4_obs = chum_total_stock * prop4,
    Nage5_obs = chum_total_stock * prop5,
    Nage6_obs = chum_total_stock * prop6,
    Nage3_pred = case_when(          # column X 
      Year <= 1953 ~ Nage3_obs,
      Year >= 1954 ~ lag(chum_model_recruits, 3) * prop3 * exp(lag(chum_ln_obs_pred, 3))),
    Nage4_pred = case_when(
      Year <= 1954 ~ Nage4_obs,
      Year >= 1955 ~ lag(chum_model_recruits, 4) * prop4 * exp(lag(chum_ln_obs_pred, 4))),
    Nage5_pred = case_when(
      Year <= 1955 ~ Nage5_obs,
      Year >= 1956 ~ lag(chum_model_recruits, 5) * prop5 * exp(lag(chum_ln_obs_pred, 5))), # why chum_ln_obs_pred here? 
    Nage6_pred = case_when(
      Year <= 1956 ~ Nage6_obs,
      Year >= 1957 ~ lag(chum_model_recruits, 6) * prop6 * exp(lag(chum_ln_obs_pred, 6))),
    recruits_pred = rowSums(across(c(Nage3_pred, Nage4_pred, Nage5_pred, Nage6_pred)), na.rm = TRUE),
    recruits_dif = chum_total_stock - recruits_pred, # deviation from the obs total stock
    U_chum_pred = catch / recruits_pred,     # U calc - column AE 
    U_chum_dif = U_chum_pred - U_chum,       # U dif 
    chum_commercial_harvest = case_when(
      Year <= 1990 ~ U_chum_pred,
      Year >= 1991 ~ U_chum),
    chum_commercial_harvest_uapply = case_when( # column AH 
      Year <= 1990 ~ U_chum_pred,
      Year >= 1991 ~ U_apply))

# set 1978 value 
SSL_1978 <- new.data %>%
  filter(Year == 1978) %>%
  pull(SSL)

# arrange and initialize values: 
df <- new.data %>%
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

# forward simulation loop: 
for (i in seq_len(nrow(df))) {
  
  # --- SSL logic ---
  SSL_alt <- if (df$Year[i] <= 1978) {
    df$SSL[i]
  } else {
    (1 - SSL_control) * df$SSL[i] +
      SSL_control * SSL_1978
  }
  
  SSL_control_alpha <-
    intercept +
    df$PDO_adult[i] * pdo_adult_coef +
    df$NPGO[i] * npgo_coef +
    df$PDO_smolt[i] * pdo_smolt_coef +
    SSL_alt * ssl_coef
  
  # --- Recruits ---
  df$chum_recruits_alt[i] <-
    if (SSL_control == 0) {
      df$chum_spawners_pred[i] *
        exp(df$chum_base_alpha[i] +
              spawners_coef * df$chum_spawners_pred[i]) *
        exp(df$chum_ln_obs_pred[i])
    } else {
      df$chum_spawners_pred[i] *
        exp(SSL_control_alpha +
              spawners_coef * df$chum_spawners_pred[i]) *
        exp(df$chum_ln_obs_pred[i])
    }
  
  # --- Age structure (lagged recruits) ---
  if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
  if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
  if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
  if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
  
  # --- Totals ---
  df$sum_alt[i] <-
    sum(df$Nage3_alt[i],
        df$Nage4_alt[i],
        df$Nage5_alt[i],
        df$Nage6_alt[i],
        na.rm = TRUE)
  
  df$catch_alt[i] <- df$sum_alt[i] * df$U_chum_pred[i]
  
  # biologically consistent spawners
  df$chum_spawners_pred[i] <-
    df$sum_alt[i] * (1 - df$U_chum_pred[i])
}


### Bind chum to steelhead data --------------------
# df to sh_data

df$Year
sh_data$Year

all_data <- df %>%
  right_join(sh_data, by = "Year")


## STEELHEAD  --------------------
### Thompson ---------------
# coefficients (from lnrs model)
sh_thompson_intercept      <- 1.572107637
sh_thompson_sst_coef       <- -0.203463091
sh_thompson_ssl_coef       <- -0.764277677
sh_thompson_npgo_coef      <- -0.0402
sh_thompson_spawners_coef  <- -0.804438793

# constants pulled from data
sh_thompson_SSL_1978 <- all_data %>%
  filter(Year == 1978) %>%
  pull(sh_thompson_SL) %>%
  as.numeric()

FN_thompson_2018 <- all_data %>%
  filter(Year == 2018) %>%
  pull(sh_thompson_FN_mortalities) %>%
  as.numeric()

# initialize state
df <- all_data %>%
  arrange(Year) %>%
  filter(!is.na(Year)) %>%
  mutate(
    sh_thompson_base_alpha =                     # might not need this 
      sh_thompson_intercept +
      sh_thompson_SST  * sh_thompson_sst_coef +
      sh_thompson_SL   * sh_thompson_ssl_coef +
      sh_thompson_NPGO * sh_thompson_npgo_coef,
    
    sh_thompson_model_recruits =
      sh_thompson_spawners *
      exp(sh_thompson_base_alpha +
            sh_thompson_spawners * sh_thompson_spawners_coef),
    
    sh_thompson_ln_obs_pred =
      log(sh_thompson_recruits / sh_thompson_model_recruits),
    
    sh_thompson_pred_bycatch =
      sh_thompson_prefishery_N -
      sh_thompson_sport_mortalities -
      sh_thompson_FN_mortalities -
      1000 * sh_thompson_spawners,
    
    sh_thompson_U =
      sh_thompson_pred_bycatch / sh_thompson_prefishery_N,
    
    # state variables (empty)
    sh_thompson_recruits_alt      = NA_real_,
    sh_thompson_spawners_pred     = sh_thompson_spawners,
    
    sh_thompson_Nage4_pred = NA_real_,
    sh_thompson_Nage5_pred = NA_real_,
    sh_thompson_Nage6_pred = NA_real_,
    sh_thompson_Nage7_pred = NA_real_,
    sh_thompson_Nage8_pred = NA_real_,
    
    sh_thompson_SSL_alt           = NA_real_,
    sh_thompson_SSL_alpha         = NA_real_,
    sh_thompson_alpha_CN          = NA_real_,   # <-- correct way to add
    sh_thompson_sum_pred          = NA_real_,
    sh_thompson_bycatch_pred      = NA_real_,
    sh_thompson_FN_catch_pred     = NA_real_,
    sh_thompson_total_catch_pred  = NA_real_,
    sh_thompson_U_comm            = NA_real_
  )

# start index
start_i <- which(df$Year >= start_year)[1]

# build a fast year->row lookup (assumes one row per Year)
year_to_i <- setNames(seq_len(nrow(df)), df$Year)

# forward simulation loop
for (i in seq(from = start_i, to = nrow(df))) {
  
  yr <- df$Year[i]
  
  # --- commercial U - column BZ 
  df$sh_thompson_U_comm[i] <-
    if (yr <= 1990) {
      df$sh_thompson_U[i]
    } else if (U_historic == 1) {
      df$sh_thompson_U[i]
    } else {
      byrate * df$chum_commercial_harvest_uapply[i]
    }
  
  # --- SSL 
  df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) {
    df$sh_thompson_SL[i]
  } else {
    (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
  }
  
  # --- SSL alpha - column CT
  df$sh_thompson_SSL_alpha[i] <-
    sh_thompson_intercept +
    df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
    df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
    df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
  
  # --- base alpha (replaces above before loop)
  df$sh_thompson_alpha_CN[i] <-
    sh_thompson_intercept +
    df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
    df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
    df$sh_thompson_SL[i]   * sh_thompson_ssl_coef
  
  # --- within-year fixed-point iteration to resolve circularity ---
  S_old <- df$sh_thompson_spawners_pred[i]
  if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]     # fallback
  if (is.na(S_old)) S_old <- 0
  
  max_iter <- 50
  tol <- 1e-8
  
  for (iter in seq_len(max_iter)) {
    
    # 1) set current guess
    df$sh_thompson_spawners_pred[i] <- S_old
    
    # 2) recruits (Excel formula)
    spk <- df$sh_thompson_spawners_pred[i] / 1000
    
    df$sh_thompson_recruits_alt[i] <-
      if (SSL_control == 0) {
        spk * exp(df$sh_thompson_alpha_CN[i] +
                    sh_thompson_spawners_coef * spk)
      } else {
        spk * exp(df$sh_thompson_SSL_alpha[i] +
                    sh_thompson_spawners_coef * spk)
      }
    
    # 3) ages -> stock (your existing lag logic, but computed inside iter)
    lag_recruits <- function(lag_year) {
      j <- year_to_i[as.character(lag_year)]
      if (is.na(j)) NA_real_ else df$sh_thompson_recruits_alt[j]
    }
    
    if (yr < start_year + 4) {
      df$sh_thompson_Nage4_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p4[i]
    } else {
      df$sh_thompson_Nage4_pred[i] <- lag_recruits(yr - 4) * df$sh_thompson_p4[i] * 1000
    }
    
    if (yr < start_year + 5) {
      df$sh_thompson_Nage5_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p5[i]
    } else {
      df$sh_thompson_Nage5_pred[i] <- lag_recruits(yr - 5) * df$sh_thompson_p5[i] * 1000
    }
    
    if (yr < start_year + 6) {
      df$sh_thompson_Nage6_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p6[i]
    } else {
      df$sh_thompson_Nage6_pred[i] <- lag_recruits(yr - 6) * df$sh_thompson_p6[i] * 1000
    }
    
    if (yr < start_year + 7) {
      df$sh_thompson_Nage7_pred[i] <- df$sh_thompson_prefishery_N[i] * df$sh_thompson_p7[i]
    } else {
      df$sh_thompson_Nage7_pred[i] <- lag_recruits(yr - 7) * df$sh_thompson_p7[i] * 1000
    }
    
    if (yr < start_year + 8) {
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
    
    # 4) catch components
    df$sh_thompson_bycatch_pred[i] <-
      df$sh_thompson_sum_pred[i] * df$sh_thompson_U_comm[i]
    
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
    
    # 5) implied new spawners
    S_new <- df$sh_thompson_sum_pred[i] - df$sh_thompson_total_catch_pred[i]
    
    # optional: keep spawners nonnegative (Excel often implicitly does)
    S_new <- max(S_new, 0)
    
    # 6) convergence check
    if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
      S_old <- S_new
      break
    }
    
    S_old <- S_new
  }
  
  # after iteration finishes, store final spawners
  df$sh_thompson_spawners_pred[i] <- S_old
  
}

### Chilcotin ------------------------ 

# coefficients (from lnrs model)
sh_chilcotin_intercept = 1.053608979
sh_chilcotin_sst_coef = -0.127949278
sh_chilcotin_ssl_coef = -0.792741195
sh_chilcotin_npgo_coef = 0.152526045
sh_chilcotin_pdo_coef = 0.202708011
sh_chilcotin_spawners_coef = -1.022467631

start_year  <- 1973

# constants pulled from data
sh_chilcotin_SSL_1978 <- df %>%
  filter(Year == 1978) %>%
  pull(sh_chilcotin_SL) %>%
  as.numeric()

FN_chilcotin_2018 <- df %>%
  filter(Year == 2018) %>%
  pull(sh_chilcotin_FN_mortalities) %>%
  as.numeric()

# initialize state
df <- df %>%
  arrange(Year) %>%
  filter(!is.na(Year)) %>%
  mutate(
    sh_chilcotin_base_alpha =                     # might not need this 
      sh_chilcotin_intercept +
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
    
    sh_chilcotin_pred_bycatch =
      sh_chilcotin_prefishery_N -
      sh_chilcotin_sport_mortalities -
      sh_chilcotin_FN_mortalities -
      1000 * sh_chilcotin_spawners,
    
    sh_chilcotin_U =
      sh_chilcotin_pred_bycatch / sh_chilcotin_prefishery_N,
    
    # state variables (empty)
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

# start index
start_i <- which(df$Year >= start_year)[1]

# build a fast year->row lookup (assumes one row per Year)
year_to_i <- setNames(seq_len(nrow(df)), df$Year)

# forward simulation loop
for (i in seq(from = start_i, to = nrow(df))) {
  
  yr <- df$Year[i]
  
  # --- commercial U - column BZ 
  df$sh_chilcotin_U_comm[i] <-
    if (yr <= 1990) {
      df$sh_chilcotin_U[i]
    } else if (U_historic == 1) {
      df$sh_chilcotin_U[i]
    } else {
      byrate * df$chum_commercial_harvest_uapply[i]
    }
  
  # --- SSL 
  df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) {
    df$sh_chilcotin_SL[i]
  } else {
    (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
  }
  
  # --- SSL alpha - column CT
  df$sh_chilcotin_SSL_alpha[i] <-
    sh_chilcotin_intercept +
    df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
    df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
    df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef +
    df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
  
  # --- base alpha (replaces above before loop)
  df$sh_chilcotin_alpha_CN[i] <-
    sh_chilcotin_intercept +
    df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
    df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
    df$sh_chilcotin_SL[i]   * sh_chilcotin_ssl_coef +
    df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef
  
  # --- within-year fixed-point iteration to resolve circularity ---
  S_old <- df$sh_chilcotin_spawners_pred[i]
  if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]     # fallback
  if (is.na(S_old)) S_old <- 0
  
  max_iter <- 50
  tol <- 1e-8
  
  for (iter in seq_len(max_iter)) {
    
    # 1) set current guess
    df$sh_chilcotin_spawners_pred[i] <- S_old
    
    # 2) recruits (Excel formula)
    spk <- df$sh_chilcotin_spawners_pred[i] / 1000
    
    df$sh_chilcotin_recruits_alt[i] <-
      if (SSL_control == 0) {
        spk * exp(df$sh_chilcotin_alpha_CN[i] +
                    sh_chilcotin_spawners_coef * spk)
      } else {
        spk * exp(df$sh_chilcotin_SSL_alpha[i] +
                    sh_chilcotin_spawners_coef * spk)
      }
    
    # 3) ages -> stock (your existing lag logic, but computed inside iter)
    lag_recruits <- function(lag_year) {
      j <- year_to_i[as.character(lag_year)]
      if (is.na(j)) NA_real_ else df$sh_chilcotin_recruits_alt[j]
    }
    
    if (yr < start_year + 4) {
      df$sh_chilcotin_Nage4_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p4[i]
    } else {
      df$sh_chilcotin_Nage4_pred[i] <- lag_recruits(yr - 4) * df$sh_chilcotin_p4[i] * 1000
    }
    
    if (yr < start_year + 5) {
      df$sh_chilcotin_Nage5_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p5[i]
    } else {
      df$sh_chilcotin_Nage5_pred[i] <- lag_recruits(yr - 5) * df$sh_chilcotin_p5[i] * 1000
    }
    
    if (yr < start_year + 6) {
      df$sh_chilcotin_Nage6_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p6[i]
    } else {
      df$sh_chilcotin_Nage6_pred[i] <- lag_recruits(yr - 6) * df$sh_chilcotin_p6[i] * 1000
    }
    
    if (yr < start_year + 7) {
      df$sh_chilcotin_Nage7_pred[i] <- df$sh_chilcotin_prefishery_N[i] * df$sh_chilcotin_p7[i]
    } else {
      df$sh_chilcotin_Nage7_pred[i] <- lag_recruits(yr - 7) * df$sh_chilcotin_p7[i] * 1000
    }
    
    if (yr < start_year + 8) {
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
    
    # 4) catch components
    df$sh_chilcotin_bycatch_pred[i] <-
      df$sh_chilcotin_sum_pred[i] * df$sh_chilcotin_U_comm[i]
    
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
    
    # 5) implied new spawners
    S_new <- df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_total_catch_pred[i]
    
    # optional: keep spawners nonnegative (Excel often implicitly does)
    S_new <- max(S_new, 0)
    
    # 6) convergence check
    if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
      S_old <- S_new
      break
    }
    
    S_old <- S_new
  }
  
  # after iteration finishes, store final spawners
  df$sh_chilcotin_spawners_pred[i] <- S_old
  
}

# no control 
scenario_1 <- df 
# with controls 
#scenario_2 <- df 

## Plots & calculations ------------------------------
# check columns 
df %>% select(Year, chum_base_alpha) %>%
  print(n=100)

# plot harvest base alphas 
scenario_1 %>%
  select(Year, chum_base_alpha, sh_thompson_base_alpha, sh_chilcotin_base_alpha, 
         U_chum, sh_chilcotin_U_comm, sh_thompson_U_comm) %>%
  filter(Year < 2017) %>%
  filter(Year > 1979) %>%
  pivot_longer(cols = 2:6, names_to = "stock", values_to = "base_alpha") %>%
  mutate(
    species = case_when(
      str_detect(stock, "chum") ~ "Chum",
      str_starts(stock, "sh_thompson") ~ "Steelhead",
      str_starts(stock, "sh_chilcotin") ~ "Steelhead",
      TRUE ~ NA_character_
    ),
    stock = recode(stock, 
                   "chum_base_alpha" = "Chum (base alpha)", 
                   "sh_thompson_base_alpha" = "Thompson steelhead (base alpha)",
                   "sh_chilcotin_base_alpha" = "Chilcotin steelhead (base alpha)",
                   "U_chum" = "Chum (commercial U)",
                   "sh_chilcotin_U_comm" = "Chilcotin steelhead (commercial U)",
                   "sh_thompson_U_comm" = "Thompson steelhead (commercial U)")
  ) %>%
  ggplot(aes(Year, base_alpha, color = stock, linetype = species)) + 
  geom_line() + 
  labs(x = "Year", y = "alpha / harvest rate", color = "Metric") + 
  theme_minimal()

scenario_2 %>%
  select(Year, chum_base_alpha, sh_thompson_base_alpha, sh_chilcotin_base_alpha, 
         U_chum, sh_chilcotin_U_comm, sh_thompson_U_comm) %>%
  filter(Year < 2017) %>%
  filter(Year > 1979) %>%
  pivot_longer(cols = 2:6, names_to = "stock", values_to = "base_alpha") %>%
  mutate(
    species = case_when(
      str_detect(stock, "chum") ~ "Chum",
      str_starts(stock, "sh_thompson") ~ "Steelhead",
      str_starts(stock, "sh_chilcotin") ~ "Steelhead",
      TRUE ~ NA_character_
    ),
    stock = recode(stock, 
                   "chum_base_alpha" = "Chum (base alpha)", 
                   "sh_thompson_base_alpha" = "Thompson steelhead (base alpha)",
                   "sh_chilcotin_base_alpha" = "Chilcotin steelhead (base alpha)",
                   "U_chum" = "Chum (commercial U)",
                   "sh_chilcotin_U_comm" = "Chilcotin steelhead (commercial U)",
                   "sh_thompson_U_comm" = "Thompson steelhead (commercial U)")
  ) %>%
  ggplot(aes(Year, base_alpha, color = stock, linetype = species)) + 
  geom_line() + 
  labs(x = "Year", y = "alpha / harvest rate", color = "Metric", linetype = "Species") + 
  theme_minimal()

# stopped here! figure out which to use to show diff between controls and no controls 




# other 
scenario_2 %>%
  select(Year, sh_thompson_SSL_alpha, sh_chilcotin_SSL_alpha, U_chum, sh_chilcotin_U_comm, sh_thompson_U_comm) %>%
  filter(Year < 2017) %>%
  filter(Year > 1979) %>%
  pivot_longer(cols = 2:5, names_to = "stock", values_to = "base_alpha") %>%
  ggplot(aes(Year, base_alpha, color = stock)) + 
  geom_line() + 
  labs(x = "Year", y = "alpha / harvest rate", color = "Metric") + 
  theme_light()



# Projections ------------------------------------------
## Thompson average alphas -------------------------
# 5 year 
alpha_5 <- sh %>% 
  filter(Year > 2014) %>%
  summarize(mean_alpha = mean(sh_thompson_base_alpha)) %>%
  pull(mean_alpha)

# 10 year 
alpha_10 <- sh %>% 
  filter(Year > 2008) %>%
  summarize(mean_alpha = mean(sh_thompson_base_alpha)) %>%
  pull(mean_alpha)

# 20 year 
alpha_20 <- sh %>% 
  filter(Year > 1998) %>%
  summarize(mean_alpha = mean(sh_thompson_base_alpha)) %>%
  pull(mean_alpha)

## Thompson average alphas -------------------------



# Predation rates ------------------------------------------------
## chum -------------------------
# toggles 
run_days <- 80

SSL_number <- read_csv("data/SSL_numbers.csv") 

# apply age props 
chum_age_3 = 0.107415165345185
chum_age_4 = 0.759302338214896
chum_age_5 = 0.131836042641171
chum_age_6 = 0.00144645379874856

# stopped column N in covariates sheet ! 
SSL_number <- SSL_number  %>% 
  mutate(SL_3 = lead(SSL_number, 3),
           SL_4 = lead(SSL_number, 4),
           SL_5 = lead(SSL_number, 5),
           SL_6 = lead(SSL_number, 6)) %>% 
  mutate(SSL_number_weighted = (chum_age_3*SL_3)+(chum_age_4*SL_4)+(chum_age_5*SL_5)+(chum_age_6*SL_6)) %>%
  filter(Year >= 1951, Year <= 2016)

df <- df %>%
  full_join(SSL_number, by = "Year")

# calculate 
chum_SSLZ = (0-mean(SSL_number$SSL_number_weighted))/sd(SSL_number$SSL_number_weighted)
sh_thompson_SSLZ
sh_chilcotin_SSLZ

# predation calculations 
df <- df %>% 
  mutate(chum_run_before_predation = chum_spawners * exp(intercept+spawners_coef*chum_spawners+pdo_adult_coef*PDO_adult+npgo_coef*NPGO+pdo_smolt_coef*PDO_smolt+ssl_coef*chum_SSLZ)) %>% # FB - SSLZ at the end (calculated above) 
  mutate(chum_return_after_predation = chum_spawners * exp(intercept+spawners_coef*chum_spawners+pdo_adult_coef*PDO_adult+npgo_coef*NPGO+pdo_smolt_coef*PDO_smolt+ssl_coef*SSL)) %>% # FC 
  mutate(chum_total_killed = chum_run_before_predation - chum_return_after_predation) %>% #FD
  mutate(chum_fraction_killed = chum_total_killed / chum_run_before_predation, #FE
         chum_percent_killed = chum_fraction_killed* 100) %>% #FF
  mutate(chum_killed_per_predator = chum_total_killed / SSL_number_weighted,
         chum_killed_per_predator_per_day = chum_killed_per_predator / run_days) %>%
  arrange(Year) %>%
  mutate(chum_killed_per_predator_rolling_avg = rollmean(chum_killed_per_predator, k = 5, fill = NA, align = "right"),
         delta_SSL= ((SSL_number_weighted - lag(SSL_number_weighted)) / lag(SSL_number_weighted))*100,
         SSL_relative_abundance = SSL_number_weighted / min(SSL_number_weighted, na.rm=T))

# check values 
df %>%
  select(Year, chum_killed_per_predator,chum_killed_per_predator_per_day, delta_SSL,SSL_relative_abundance) %>%
  print(n=100)

###  plots -------------
# chums killed per predator (60 days)
ggplot(df, aes(Year, chum_killed_per_predator)) +
  geom_line() +
  theme_light()


# different day lengths 
chum_killed_60 <- df %>% select(Year, chum_killed_60 = chum_killed_per_predator)
chum_killed_70 <- df %>% select(Year, chum_killed_70 = chum_killed_per_predator)
chum_killed_80 <- df %>% select(Year, chum_killed_70 = chum_killed_per_predator)

chum_killed <- chum_killed_60 %>% 
  mutate(run_days = "60 Days") %>% 
  rename(chum_killed_per_predator = chum_killed_60) %>%
  bind_rows(chum_killed_70 %>% 
              mutate(run_days = "70 Days") %>% 
              rename(chum_killed_per_predator = chum_killed_70)) %>%
  bind_rows(chum_killed_80 %>% 
              mutate(run_days = "80 Days") %>% 
              rename(chum_killed_per_predator = chum_killed_70))

# plot 
ggplot(chum_killed, aes(Year, chum_killed_per_predator, color = run_days)) +
  geom_line() +
  theme_light()

## steelhead --------------------
# toggles 
run_days <- 60

SSL_number <- read_csv("data/SSL_numbers.csv") 

# apply age props 
SH_chilcotin_age_4 = 0.00761421319796954
SH_chilcotin_age_5 = 0.295121697253677
SH_chilcotin_age_6 = 0.606023688663283
SH_chilcotin_age_7 = 0.0887023298190811
SH_chilcotin_age_8 = 0.00253807106598985

SH_thompson_age_4 = 0.0357506950880445
SH_thompson_age_5 = 0.835875810936052
SH_thompson_age_6 = 0.119569045412419
SH_thompson_age_7 = 0.00834105653382762
SH_thompson_age_8 = 0.00046339202965709

SSL_number <- SSL_number  %>% 
  mutate(SL_4 = lead(SSL_number, 4),
         SL_5 = lead(SSL_number, 5),
         SL_6 = lead(SSL_number, 6),
         SL_7 = lead(SSL_number, 7),
         SL_8 = lead(SSL_number, 8)) %>% 
  mutate(SSL_number_weighted_thompson = (SH_thompson_age_4*SL_4)+(SH_thompson_age_5*SL_5)+(SH_thompson_age_6*SL_6)+(SH_thompson_age_7*SL_7)+(SH_thompson_age_8*SL_8),
         SSL_number_weighted_chilcotin = (SH_chilcotin_age_4*SL_4)+(SH_chilcotin_age_5*SL_5)+(SH_chilcotin_age_6*SL_6)+(SH_chilcotin_age_7*SL_7)+(SH_chilcotin_age_8*SL_8)) %>%
  filter(Year >= 1951, Year <= 2016)

df <- df %>%
  full_join(SSL_number, by = "Year")

# calculate SL at 0 
sh_thompson_SSLZ = (0-mean(SSL_number$SSL_number_weighted_thompson, na.rm = T))/sd(SSL_number$SSL_number_weighted_thompson, na.rm = T)
sh_chilcotin_SSLZ = (0-mean(SSL_number$SSL_number_weighted_chilcotin, na.rm = T))/sd(SSL_number$SSL_number_weighted_chilcotin, na.rm = T)

# predation calculations 
# thompson 
df <- df %>% 
  # thompson
  mutate(sh_thompson_run_before_predation = sh_thompson_spawners * exp(sh_thompson_intercept+sh_thompson_spawners_coef*sh_thompson_spawners+sh_thompson_sst_coef*sh_thompson_SST+sh_thompson_npgo_coef*sh_thompson_NPGO+sh_thompson_ssl_coef*sh_thompson_SSLZ)) %>% # FB 
  mutate(sh_thompson_return_after_predation = sh_thompson_spawners * exp(sh_thompson_intercept+sh_thompson_spawners_coef*sh_thompson_spawners+sh_thompson_sst_coef*sh_thompson_SST+sh_thompson_npgo_coef*sh_thompson_NPGO+sh_thompson_ssl_coef*sh_thompson_SL)) %>% # FC 
  mutate(sh_thompson_total_killed = sh_thompson_run_before_predation - sh_thompson_return_after_predation) %>% #FD
  mutate(sh_thompson_fraction_killed = sh_thompson_total_killed / sh_thompson_run_before_predation, #FE
         sh_thompson_percent_killed = sh_thompson_fraction_killed* 100) %>% #FF
  mutate(sh_thompson_killed_per_predator = sh_thompson_total_killed / SSL_number_weighted_thompson,
         sh_thompson_killed_per_predator_per_day = sh_thompson_killed_per_predator / run_days) %>%
  arrange(Year) %>%
  mutate(sh_thompson_killed_per_predator_rolling_avg = rollmean(sh_thompson_killed_per_predator, k = 5, fill = NA, align = "right"),
         delta_SSL= ((SSL_number_weighted_thompson - lag(SSL_number_weighted_thompson)) / lag(SSL_number_weighted_thompson))*100,
         SSL_relative_abundance_thompson = SSL_number_weighted_thompson / min(SSL_number_weighted_thompson, na.rm=T)) %>%
  # chilcotin
  mutate(sh_chilcotin_run_before_predation = sh_chilcotin_spawners * exp(sh_chilcotin_intercept+sh_chilcotin_spawners_coef*sh_chilcotin_spawners+sh_chilcotin_sst_coef*sh_chilcotin_SST+sh_chilcotin_npgo_coef*sh_chilcotin_NPGO+sh_chilcotin_PDO*sh_chilcotin_pdo_coef+sh_chilcotin_ssl_coef*sh_chilcotin_SSLZ)) %>% # FB 
  mutate(sh_chilcotin_return_after_predation = sh_chilcotin_spawners * exp(sh_chilcotin_intercept+sh_chilcotin_spawners_coef*sh_chilcotin_spawners+sh_chilcotin_sst_coef*sh_chilcotin_SST+sh_chilcotin_npgo_coef*sh_chilcotin_NPGO+sh_chilcotin_PDO*sh_chilcotin_pdo_coef+sh_chilcotin_ssl_coef*sh_chilcotin_SL)) %>% # FC 
  mutate(sh_chilcotin_total_killed = sh_chilcotin_run_before_predation - sh_chilcotin_return_after_predation) %>% #FD
  mutate(sh_chilcotin_fraction_killed = sh_chilcotin_total_killed / sh_chilcotin_run_before_predation, #FE
         sh_chilcotin_percent_killed = sh_chilcotin_fraction_killed* 100) %>% #FF
  mutate(sh_chilcotin_killed_per_predator = sh_chilcotin_total_killed / SSL_number_weighted_chilcotin,
         sh_chilcotin_killed_per_predator_per_day = sh_chilcotin_killed_per_predator / run_days) %>%
  arrange(Year) %>%
  mutate(sh_chilcotin_killed_per_predator_rolling_avg = rollmean(sh_chilcotin_killed_per_predator, k = 5, fill = NA, align = "right"),
         delta_SSL= ((SSL_number_weighted_chilcotin - lag(SSL_number_weighted_chilcotin)) / lag(SSL_number_weighted_chilcotin))*100,
         SSL_relative_abundance_chilcotin = SSL_number_weighted_chilcotin / min(SSL_number_weighted_chilcotin, na.rm=T)) 

# check values 
df %>%
  select(Year, sh_chilcotin_total_killed) %>%
  print(n=100)


### plots ---------------
#  killed per predator (60 days)
ggplot(df, aes(Year)) +
  geom_line(aes(y=sh_thompson_killed_per_predator)) +
  geom_line(aes(y=sh_chilcotin_killed_per_predator), col = "red") +
  theme_light()







