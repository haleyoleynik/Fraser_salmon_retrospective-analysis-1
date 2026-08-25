# Retrospective Model 
# Haley Oleynik Murdoch McAllister
# October 2025

## original chum retrospective model 
# see chum-retrospective-model.R for cleaned and updated code. 

# load libraries 
require(tidyverse)
require(ggplot2)
require(readr)
require(slider)
require(zoo)
require(patchwork)
require(scales)

# read data ------------------
data <- read_csv("R/Chum Steelhead Retrospective Shiny App/s-r_data.csv") %>%
  select(-SSL)
sh_data <- read_csv("R/Chum Steelhead Retrospective Shiny App/sh_s-r_data.csv")
covariates <- read_csv("R/Chum Steelhead Retrospective Shiny App/covariates.csv") 
# Scenario 1 --------------------------------------
# no predator controls 
## CHUM ----------------
# from Fraser Chum data v19_alpha_SSL_est_fin_yrs_v8 & v9 sheet 
# calculate covariates in situ !!!!!!!!!!!!!!!!!

# estimate the coefficients from regression models in situ !!!!!!!!!!!!!!!
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
    catch = chum_total_stock - chum_spawners, # catch 
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
    chum_spawners_pred = chum_spawners,
    chum_SSL_alt = NA_real_,
    chum_SSL_alpha = NA_real_ 
  )

# forward simulation loop: 
for (i in seq_len(nrow(df))) {
  
  # SSL logic 
  df$chum_SSL_alt[i] <- if (df$Year[i] <= 1978) {
    df$SSL[i]
  } else {
    (1 - SSL_control) * df$SSL[i] +
      SSL_control * SSL_1978
  }
  
  df$chum_SSL_alpha[i] <-
    intercept +
    df$PDO_adult[i] * pdo_adult_coef +
    df$NPGO[i] * npgo_coef +
    df$PDO_smolt[i] * pdo_smolt_coef +
    df$chum_SSL_alt[i] * ssl_coef
  
  # Recruits
  df$chum_recruits_alt[i] <-
    if (SSL_control == 0) {
      df$chum_spawners_pred[i] *
        exp(df$chum_base_alpha[i] +
              spawners_coef * df$chum_spawners_pred[i]) *
        exp(df$chum_ln_obs_pred[i])
    } else {
      df$chum_spawners_pred[i] *
        exp(df$chum_SSL_alpha[i] +
              spawners_coef * df$chum_spawners_pred[i]) *
        exp(df$chum_ln_obs_pred[i])
    }
  
# Age structure (lagged recruits
  if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
  if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
  if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
  if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
  
  #  Totals
  df$sum_alt[i] <-
    sum(df$Nage3_alt[i],
        df$Nage4_alt[i],
        df$Nage5_alt[i],
        df$Nage6_alt[i],
        na.rm = TRUE)
  
  #df$catch_alt[i] <- df$sum_alt[i] * df$U_chum_pred[i]
  df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
  
  # biologically consistent spawners
  #df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$U_chum_pred[i])
  df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i]) # add Uapply
  
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
  
  #  commercial U - column BZ 
  df$sh_thompson_U_comm[i] <-
    if (yr <= 1990) {
      df$sh_thompson_U[i]
    } else if (U_historic == 1) {
      df$sh_thompson_U[i]
    } else {
      byrate * df$chum_commercial_harvest_uapply[i]
    }
  
  #  SSL 
  df$sh_thompson_SSL_alt[i] <- if (yr <= 1978) {
    df$sh_thompson_SL[i]
  } else {
    (1 - SSL_control) * df$sh_thompson_SL[i] + SSL_control * sh_thompson_SSL_1978
  }
  
  #  SSL alpha - column CT
  df$sh_thompson_SSL_alpha[i] <-
    sh_thompson_intercept +
    df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
    df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
    df$sh_thompson_SSL_alt[i] * sh_thompson_ssl_coef
  
  #  base alpha (replaces above before loop)
  df$sh_thompson_alpha_CN[i] <-
    sh_thompson_intercept +
    df$sh_thompson_SST[i]  * sh_thompson_sst_coef +
    df$sh_thompson_NPGO[i] * sh_thompson_npgo_coef +
    df$sh_thompson_SL[i]   * sh_thompson_ssl_coef
  
  # within-year fixed-point iteration to resolve circularity ---
  S_old <- df$sh_thompson_spawners_pred[i]
  if (is.na(S_old)) S_old <- df$sh_thompson_spawners[i]     # fallback
  if (is.na(S_old)) S_old <- 0
  
  max_iter <- 50
  tol <- 1e-8
  
  for (iter in seq_len(max_iter)) {
    
    # set current guess
    df$sh_thompson_spawners_pred[i] <- S_old
    
    # recruits (Excel formula)
    spk <- df$sh_thompson_spawners_pred[i] / 1000
    
    df$sh_thompson_recruits_alt[i] <-
      if (SSL_control == 0) {
        spk * exp(df$sh_thompson_alpha_CN[i] +
                    sh_thompson_spawners_coef * spk)
      } else {
        spk * exp(df$sh_thompson_SSL_alpha[i] +
                    sh_thompson_spawners_coef * spk)
      }
    
    # ages -> stock (your existing lag logic, but computed inside iter)
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
    
    # catch components
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
    
    # implied new spawners
    S_new <- df$sh_thompson_sum_pred[i] - df$sh_thompson_total_catch_pred[i]
    
    # optional: keep spawners nonnegative (Excel often implicitly does)
    S_new <- max(S_new, 0)
    
    # convergence check
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
  
  # commercial U - column BZ 
  df$sh_chilcotin_U_comm[i] <-
    if (yr <= 1990) {
      df$sh_chilcotin_U[i]
    } else if (U_historic == 1) {
      df$sh_chilcotin_U[i]
    } else {
      byrate * df$chum_commercial_harvest_uapply[i]
    }
  
  # SSL 
  df$sh_chilcotin_SSL_alt[i] <- if (yr <= 1973) {
    df$sh_chilcotin_SL[i]
  } else {
    (1 - SSL_control) * df$sh_chilcotin_SL[i] + SSL_control * sh_chilcotin_SSL_1978
  }
  
  # SSL alpha - column CT
  df$sh_chilcotin_SSL_alpha[i] <-
    sh_chilcotin_intercept +
    df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
    df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
    df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef +
    df$sh_chilcotin_SSL_alt[i] * sh_chilcotin_ssl_coef
  
  # base alpha (replaces above before loop)
  df$sh_chilcotin_alpha_CN[i] <-
    sh_chilcotin_intercept +
    df$sh_chilcotin_SST[i]  * sh_chilcotin_sst_coef +
    df$sh_chilcotin_NPGO[i] * sh_chilcotin_npgo_coef +
    df$sh_chilcotin_SL[i]   * sh_chilcotin_ssl_coef +
    df$sh_chilcotin_PDO[i] * sh_chilcotin_pdo_coef
  
  #  within-year fixed-point iteration to resolve circularity
  S_old <- df$sh_chilcotin_spawners_pred[i]
  if (is.na(S_old)) S_old <- df$sh_chilcotin_spawners[i]     # fallback
  if (is.na(S_old)) S_old <- 0
  
  max_iter <- 50
  tol <- 1e-8
  
  for (iter in seq_len(max_iter)) {
    
    # set current guess
    df$sh_chilcotin_spawners_pred[i] <- S_old
    
    # recruits (Excel formula)
    spk <- df$sh_chilcotin_spawners_pred[i] / 1000
    
    df$sh_chilcotin_recruits_alt[i] <-
      if (SSL_control == 0) {
        spk * exp(df$sh_chilcotin_alpha_CN[i] +
                    sh_chilcotin_spawners_coef * spk)
      } else {
        spk * exp(df$sh_chilcotin_SSL_alpha[i] +
                    sh_chilcotin_spawners_coef * spk)
      }
    
    # ages -> stock (your existing lag logic, but computed inside iter)
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
    
    # catch components
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
    
    # implied new spawners
    S_new <- df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_total_catch_pred[i]
    
    # optional: keep spawners nonnegative (Excel often implicitly does)
    S_new <- max(S_new, 0)
    
    # convergence check
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

# Plots & calculations ------------------------------
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

# parameters of interest 
sh_chilcotin_SSL_alpha
sh_thompson_SSL_alpha
chum_SSL_alpha

# SSL controls on!  
p1 <- scenario_2 %>%
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
  geom_line(size = 1) + 
  labs(x = "Year", y = "alpha / harvest rate", color = "Metric", linetype = "Species") + 
  ylim(0,3) + 
  scale_color_manual(values = hcl.colors(5, palette = "Dark 2")) +
  theme_minimal()


p2 <- scenario_2 %>%
  select(Year, chum_SSL_alpha, sh_thompson_SSL_alpha, sh_chilcotin_SSL_alpha, 
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
                   "chum_SSL_alpha" = "Chum alpha", 
                   "sh_thompson_SSL_alpha" = "Thompson steelhead alpha",
                   "sh_chilcotin_SSL_alpha" = "Chilcotin steelhead alpha",
                   "U_chum" = "Chum (commercial U)",
                   "sh_chilcotin_U_comm" = "Chilcotin steelhead (commercial U)",
                   "sh_thompson_U_comm" = "Thompson steelhead (commercial U)")
  ) %>%
  ggplot(aes(Year, base_alpha, color = stock, linetype = species)) + 
  geom_line(size = 1) + 
  labs(x = "Year", y = "alpha / harvest rate", color = "Metric", linetype = "Species") + 
  ylim(0,3) + 
  scale_color_manual(values = hcl.colors(5, palette = "Dark 2")) +
  theme_minimal()

(p1 + theme(legend.position = "none")) | p2 

## alpha plots ---------------------------------------
# SSL controls on!  
(p1 <- scenario_2 %>%
  select(Year, chum_base_alpha, sh_thompson_base_alpha, sh_chilcotin_base_alpha) %>%
  filter(Year < 2017) %>%
  filter(Year > 1979) %>%
  pivot_longer(cols = 2:4, names_to = "stock", values_to = "base_alpha") %>%
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
                   "sh_chilcotin_base_alpha" = "Chilcotin steelhead (base alpha)")
  ) %>%
  ggplot(aes(Year, base_alpha, color = stock, linetype = species)) + 
  geom_line(size = 1) + 
  labs(x = "Year", y = "Alpha", color = "Metric", linetype = "Species") + 
  ylim(0,3) + 
  scale_color_manual(values = hcl.colors(3, palette = "Dark 2")) +
  theme_minimal())


(p2 <- scenario_2 %>%
  select(Year, chum_SSL_alpha, sh_thompson_SSL_alpha, sh_chilcotin_SSL_alpha) %>%
  filter(Year < 2017) %>%
  filter(Year > 1979) %>%
  pivot_longer(cols = 2:4, names_to = "stock", values_to = "base_alpha") %>%
  mutate(
    species = case_when(
      str_detect(stock, "chum") ~ "Chum",
      str_starts(stock, "sh_thompson") ~ "Steelhead",
      str_starts(stock, "sh_chilcotin") ~ "Steelhead",
      TRUE ~ NA_character_
    ),
    stock = recode(stock, 
                   "chum_SSL_alpha" = "Chum alpha", 
                   "sh_thompson_SSL_alpha" = "Thompson steelhead alpha",
                   "sh_chilcotin_SSL_alpha" = "Chilcotin steelhead alpha")
  ) %>%
  ggplot(aes(Year, base_alpha, color = stock, linetype = species)) + 
  geom_line(size = 1) + 
  labs(x = "Year", y = "", color = "Metric", linetype = "Species") + 
  ylim(0,3) + 
  scale_color_manual(values = hcl.colors(3, palette = "Dark 2")) +
  theme_minimal())

(p1 + theme(legend.position = "none")) | p2 + ylab("") + theme(legend.position = "right")

ggsave(file = "figures/chum_steelhead_alpha.png", height = 5, width = 11, dpi = 600)

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


# Run scenarios -----------------------------
# wrap original retrospective model in a function 
# chum
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
      #chum_commercial_harvest_uapply = case_when(
      #  Year <= 1990 ~ U_chum_pred,
      #  Year >= 1991 ~ U_apply)
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
    
    # SSL logic (unchanged)
    df$chum_SSL_alt[i] <- if (df$Year[i] <= 1978) {
      df$SSL[i]
    } else {
      (1 - SSL_control) * df$SSL[i] + SSL_control * SSL_1978
    }
    
    df$chum_SSL_alpha[i] <-
      intercept + df$PDO_adult[i] * pdo_adult_coef + df$NPGO[i] * npgo_coef +
      df$PDO_smolt[i] * pdo_smolt_coef + df$chum_SSL_alt[i] * ssl_coef
    
    # age structure — from PAST recruits_alt only, so this can run first
    if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
    if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
    if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
    if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
    
    # total returning stock this year
    df$sum_alt[i] <- sum(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i], na.rm = TRUE)
    
    # apply the scenario harvest rate -> catch & escapement, BEFORE recruits are calculated
    df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
    df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i])
    
    # compute this year's recruits, using the just-finalized post-harvest escapement
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
  
  # for (i in seq_len(nrow(df))) {
  #   
  #   df$chum_SSL_alt[i] <- if (df$Year[i] <= 1978) {
  #     df$SSL[i]
  #   } else {
  #     (1 - SSL_control) * df$SSL[i] + SSL_control * SSL_1978
  #   }
  #   
  #   df$chum_SSL_alpha[i] <-
  #     intercept + df$PDO_adult[i] * pdo_adult_coef + df$NPGO[i] * npgo_coef +
  #     df$PDO_smolt[i] * pdo_smolt_coef + df$chum_SSL_alt[i] * ssl_coef
  #   
  #   df$chum_recruits_alt[i] <-
  #     if (SSL_control == 0) {
  #       df$chum_spawners_pred[i] *
  #         exp(df$chum_base_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
  #         exp(df$chum_ln_obs_pred[i])
  #     } else {
  #       df$chum_spawners_pred[i] *
  #         exp(df$chum_SSL_alpha[i] + spawners_coef * df$chum_spawners_pred[i]) *
  #         exp(df$chum_ln_obs_pred[i])
  #     }
  #   
  #   if (i > 3) df$Nage3_alt[i] <- df$chum_recruits_alt[i - 3] * df$prop3[i]
  #   if (i > 4) df$Nage4_alt[i] <- df$chum_recruits_alt[i - 4] * df$prop4[i]
  #   if (i > 5) df$Nage5_alt[i] <- df$chum_recruits_alt[i - 5] * df$prop5[i]
  #   if (i > 6) df$Nage6_alt[i] <- df$chum_recruits_alt[i - 6] * df$prop6[i]
  #   
  #   df$sum_alt[i] <- sum(df$Nage3_alt[i], df$Nage4_alt[i], df$Nage5_alt[i], df$Nage6_alt[i], na.rm = TRUE)
  #   #df$catch_alt[i] <- df$sum_alt[i] * df$U_chum_pred[i]
  #   df$catch_alt[i] <- df$sum_alt[i] * df$chum_commercial_harvest_uapply[i]
  #   
  #   #df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$U_chum_pred[i])
  #   df$chum_spawners_pred[i] <- df$sum_alt[i] * (1 - df$chum_commercial_harvest_uapply[i])
  #   
  # }
  
  df %>%
    mutate(
      harvest_rate = chum_commercial_harvest_uapply,
      ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control"),
      productivity = if (SSL_control == 1) chum_SSL_alpha else chum_base_alpha,
      returns = chum_recruits_alt,
      recruits_per_spawner = chum_recruits_alt / chum_spawners_pred
    )
}

# steelhead 
run_thompson_scenario <- function(all_data, SSL_control, U_historic, byrate,
                                  start_year = 1978,
                                  sh_thompson_intercept      = 1.572107637,
                                  sh_thompson_sst_coef       = -0.203463091,
                                  sh_thompson_ssl_coef       = -0.764277677,
                                  sh_thompson_npgo_coef      = -0.0402,
                                  sh_thompson_spawners_coef  = -0.804438793) {
  
  sh_thompson_SSL_1978 <- all_data %>%
    filter(Year == 1978) %>%
    pull(sh_thompson_SL) %>%
    as.numeric()
  
  FN_thompson_2018 <- all_data %>%
    filter(Year == 2018) %>%
    pull(sh_thompson_FN_mortalities) %>%
    as.numeric()
  
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
          spk * exp(df$sh_thompson_alpha_CN[i] +
                      sh_thompson_spawners_coef * spk)
        } else {
          spk * exp(df$sh_thompson_SSL_alpha[i] +
                      sh_thompson_spawners_coef * spk)
        }
      
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
      
      S_new <- df$sh_thompson_sum_pred[i] - df$sh_thompson_total_catch_pred[i]
      S_new <- max(S_new, 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      
      S_old <- S_new
    }
    
    df$sh_thompson_spawners_pred[i] <- S_old
  }
  
  df
}

run_chilcotin_scenario <- function(df, SSL_control, U_historic, byrate,
                                   start_year = 1973,
                                   sh_chilcotin_intercept = 1.053608979,
                                   sh_chilcotin_sst_coef = -0.127949278,
                                   sh_chilcotin_ssl_coef = -0.792741195,
                                   sh_chilcotin_npgo_coef = 0.152526045,
                                   sh_chilcotin_pdo_coef = 0.202708011,
                                   sh_chilcotin_spawners_coef = -1.022467631) {
  
  sh_chilcotin_SSL_1978 <- df %>%
    filter(Year == 1978) %>%
    pull(sh_chilcotin_SL) %>%
    as.numeric()
  
  FN_chilcotin_2018 <- df %>%
    filter(Year == 2018) %>%
    pull(sh_chilcotin_FN_mortalities) %>%
    as.numeric()
  
  df <- df %>%
    arrange(Year) %>%
    filter(!is.na(Year)) %>%
    mutate(
      sh_chilcotin_base_alpha =
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
  
  start_i <- which(df$Year >= start_year)[1]
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
    
    max_iter <- 50
    tol <- 1e-8
    
    for (iter in seq_len(max_iter)) {
      
      df$sh_chilcotin_spawners_pred[i] <- S_old
      spk <- df$sh_chilcotin_spawners_pred[i] / 1000
      
      df$sh_chilcotin_recruits_alt[i] <-
        if (SSL_control == 0) {
          spk * exp(df$sh_chilcotin_alpha_CN[i] +
                      sh_chilcotin_spawners_coef * spk)
        } else {
          spk * exp(df$sh_chilcotin_SSL_alpha[i] +
                      sh_chilcotin_spawners_coef * spk)
        }
      
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
      
      S_new <- df$sh_chilcotin_sum_pred[i] - df$sh_chilcotin_total_catch_pred[i]
      S_new <- max(S_new, 0)
      
      if (is.finite(S_old) && is.finite(S_new) && abs(S_new - S_old) <= tol * max(1, abs(S_old))) {
        S_old <- S_new
        break
      }
      
      S_old <- S_new
    }
    
    df$sh_chilcotin_spawners_pred[i] <- S_old
  }
  
  df
}

## run scenarios (set harvest rates) --------------------------
harvest_rates <- seq(0, 0.8, by = 0.1)   
bycatch_rates <- seq(0, 0.8, by = 0.1)   

# run full sceanrio
run_full_scenario <- function(chum_data, covariates, sh_data,
                              U_apply, SSL_control, U_historic = 0, byrate = 0.69) {
  
  chum_df <- run_chum_scenario(chum_data, covariates, U_apply, SSL_control)
  all_data <- chum_df %>% right_join(sh_data, by = "Year")
  thompson_df <- run_thompson_scenario(all_data, SSL_control, U_historic, byrate)
  full_df <- run_chilcotin_scenario(thompson_df, SSL_control, U_historic, byrate)
  
  full_df %>%
    mutate(harvest_rate = U_apply,
           ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control"))
}

# make a scenario grid 
#scenario_grid <- expand_grid(SSL_control = c(0, 1), U_apply = harvest_rates)

# with alternative bycatch rates 
scenario_grid <- expand_grid(
  SSL_control = c(0, 1),
  U_apply = harvest_rates,
  byrate = bycatch_rates
)

#all_scenarios <- purrr::pmap_dfr(
#  scenario_grid,
#  function(SSL_control, U_apply) {
#    run_full_scenario(chum_data = data, covariates = covariates, sh_data = sh_data,
#                      U_apply = U_apply, SSL_control = SSL_control)
#  }
#)

# with bycatch rates 
all_scenarios <- purrr::pmap_dfr(
  scenario_grid,
  function(SSL_control, U_apply, byrate) {
    run_full_scenario(
      chum_data = data,
      covariates = covariates,
      sh_data = sh_data,
      U_apply = U_apply,
      SSL_control = SSL_control,
      U_historic = 0,
      byrate = byrate
    )
  }
)



# historic scenarios ---------------------------------------

# run full sceanrio
run_full_scenario <- function(chum_data, covariates, sh_data,
                              U_apply, SSL_control, U_historic = 0, byrate = 0.69) {
  
  #chum_df <- run_chum_scenario(chum_data, covariates, U_apply, SSL_control)
  chum_df <- run_chum_scenario(
    chum_data,
    covariates,
    U_apply,
    SSL_control,
    U_historic
  )
  all_data <- chum_df %>% right_join(sh_data, by = "Year")
  thompson_df <- run_thompson_scenario(all_data, SSL_control, U_historic, byrate)
  full_df <- run_chilcotin_scenario(thompson_df, SSL_control, U_historic, byrate)
  
  full_df %>%
    mutate(
      ssl_scenario = if_else(
        SSL_control == 1,
        "SSL control",
        "No SSL control"
      )
    )
}

historic_scenarios <- purrr::map_dfr(
  c(0, 1),
  function(SSL_control) {
    run_full_scenario(
      chum_data = data,
      covariates = covariates,
      sh_data = sh_data,
      U_apply = 0,
      SSL_control = SSL_control,
      U_historic = 1,
      byrate = 0.69
    )
  }
)

# plot historical conditions ---------------------------------
historic_data <- all_scenarios %>%
  filter(harvest_rate == harvest_rates[1], ssl_scenario == "No SSL control") %>%
  select(Year, catch, chum_spawners, chum_total_stock, U_chum) %>%
  distinct()

historic_long <- historic_data %>%
  pivot_longer(cols = c(catch, chum_spawners, chum_total_stock, U_chum),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         catch = "Historic catch",
                         chum_spawners = "Historic spawners",
                         chum_total_stock = "Historic total stock",
                         U_chum = "Historic commercial harvest rate"))

ggplot(historic_long, aes(Year, value)) +
  geom_line(color = "#5E4FA2", linewidth = 1) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  labs(x = "Year", y = "") +
  theme_minimal()

# find the max of the left-axis metrics to set a sensible scaling factor
left_max <- max(c(historic_data$catch, historic_data$chum_spawners, 
                  historic_data$chum_total_stock), na.rm = TRUE)

scale_factor <- left_max / max(historic_data$U_chum, na.rm = TRUE)

stock_long <- historic_data %>%
  pivot_longer(cols = c(catch, chum_spawners, chum_total_stock),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         catch = "Catch",
                         chum_spawners = "Spawners",
                         chum_total_stock = "Total stock"))

ggplot() +
  geom_line(data = stock_long, aes(Year, value, color = metric), linewidth = 1) +
  geom_line(data = historic_data, aes(Year, U_chum * scale_factor), 
            color = "black", linewidth = 1, linetype = "dashed") +
  scale_y_continuous(
    name = "Fish",
    labels = scales::label_comma(),
    sec.axis = sec_axis(~ . / scale_factor, name = "Commercial harvest rate", 
                        labels = scales::label_percent())
  ) +
  scale_color_manual(values = palette.colors(n = 3, palette = "Okabe-Ito")) +
  labs(x = "Year", color = NULL) +
  theme_minimal() +
  theme(axis.title.y.right = element_text(angle = 90))


left_max <- max(c(historic_data$catch, historic_data$chum_spawners, 
                  historic_data$chum_total_stock), na.rm = TRUE)

scale_factor <- left_max / max(historic_data$U_chum, na.rm = TRUE)

stock_long <- historic_data %>%
  pivot_longer(cols = c(catch, chum_spawners, chum_total_stock),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         catch = "Catch",
                         chum_spawners = "Spawners",
                         chum_total_stock = "Total stock"))

harvest_long <- historic_data %>%
  select(Year, U_chum) %>%
  mutate(value = U_chum * scale_factor,
         metric = "Commercial harvest rate")

combined <- bind_rows(stock_long, harvest_long %>% select(Year, metric, value))

ggplot(combined, aes(Year, value, color = metric, 
                     linetype = metric == "Commercial harvest rate")) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    name = "Fish",
    labels = scales::label_comma(),
    sec.axis = sec_axis(~ . / scale_factor, name = "Commercial harvest rate", 
                        labels = scales::label_percent())
  ) +
  scale_color_manual(values = palette.colors(n = 4, palette = "Okabe-Ito")) +
  scale_linetype_manual(values = c("solid", "dashed"), guide = "none") +
  labs(x = "Year", color = NULL) +
  geom_vline(xintercept = 1978, color = "darkgrey", linetype = "dashed") +
  theme_minimal() + 
  theme(legend.position = "top")

ggsave(file = "figures/chum_historical_conditions.png", height = 6, width = 10, dpi = 600)


# plot productivity under different scenarios ------------------------------
plot_data <- all_scenarios %>%
  select(Year, harvest_rate, ssl_scenario, productivity, returns) %>%
  pivot_longer(cols = c(productivity, returns), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         productivity = "Productivity (alpha)",
                         returns = "Recruits/returns"))

table(plot_data$ssl_scenario)

plot_data %>%
  filter(Year > 1979) %>%
  filter(Year < 2017) %>%
ggplot(aes(Year, value, color = ssl_scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(harvest_rate~metric,  ncol = 2, scales = "free_y") +
  labs(x = "Year", y = "") +
  theme_minimal(base_size = 9)

ggsave(file = "figures/chum_harvest_rate_pred_control.png", height = 9, width = 6, dpi = 600)

# plot with predicted lnrs
plot_data <- all_scenarios %>%
  select(Year, harvest_rate, ssl_scenario, productivity, returns, recruits_per_spawner) %>%
  pivot_longer(cols = c(productivity, returns, recruits_per_spawner), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         productivity = "Productivity (alpha)",
                         returns = "Recruits/returns",
                         recruits_per_spawner = "Recruits per spawner"))

# plot with predicted lnrs
plot_data <- all_scenarios %>%
  select(Year, harvest_rate, ssl_scenario, returns, recruits_per_spawner) %>%
  pivot_longer(cols = c(returns, recruits_per_spawner), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         returns = "Returns",
                         recruits_per_spawner = "Recruits per spawner"))

table(plot_data$ssl_scenario)

plot_data %>%
  filter(Year > 1979) %>%
  filter(Year < 2017) %>%
  ggplot(aes(Year, value, color = ssl_scenario, linetype = ssl_scenario)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  facet_wrap(harvest_rate~metric, ncol = 2, scales = "free_y") +
  labs(x = "Year", y = "") +
  theme_minimal(base_size = 9)

ggsave(file = "figures/chum_harvest_rate_pred_control_v2.png", height = 10, width = 6, dpi = 600)

# lost catch historic harvest rate scenario ------------------------------

catch_lost_to_pinnipeds <- historic_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Stock, Year, ssl_scenario, catch_alt) %>%
  pivot_wider(
    names_from = ssl_scenario,
    values_from = catch_alt
  ) %>%
  mutate(
    catch_lost = `No SSL control` - `SSL control`
  )


catch_lost_annual <- catch_lost_to_pinnipeds %>%
  group_by(Year) %>%
  summarise(
    catch_lost = sum(catch_lost, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(catch_lost_annual, aes(Year, catch_lost)) +
  geom_area(alpha = 0.2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Year",
    y = "Catch lost to pinnipeds"
  ) +
  theme_minimal()


catch_lost_cum <- catch_lost_annual %>%
  arrange(Year) %>%
  mutate(
    cum_catch_lost = cumsum(catch_lost)
  )

ggplot(catch_lost_cum, aes(Year, cum_catch_lost)) +
  geom_area(alpha = 0.2, fill = "#FF4500") +
  geom_line(linewidth = 1, color = "#FF4500") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Year",
    y = "Cumulative catch lost to sea lions"
  ) +
  theme_minimal()

ggsave(file = "figures/chum_catch_lost_Uhistoric.png", height = 6,width = 9, dpi = 600)

## lost catch Uapply ----------------------------
lost_catch_summary <- all_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, ssl_scenario, catch_alt) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  group_by(harvest_rate) %>%
  summarise(
    total_catch_no_control = sum(`No SSL control`, na.rm = TRUE),
    total_catch_ssl_control = sum(`SSL control`, na.rm = TRUE),
    total_lost_catch = total_catch_no_control - total_catch_ssl_control,
    .groups = "drop"
  )

lost_catch_summary

# with bycatch rates 
catch_lost_summary <- all_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, byrate, ssl_scenario, catch_alt) %>%
  pivot_wider(
    names_from = ssl_scenario,
    values_from = catch_alt
  ) %>%
  group_by(harvest_rate, byrate) %>%
  summarise(
    total_catch_no_control = sum(`No SSL control`, na.rm = TRUE),
    total_catch_ssl_control = sum(`SSL control`, na.rm = TRUE),
    total_lost_catch = total_catch_no_control - total_catch_ssl_control,
    .groups = "drop"
  )

# lost catch year to year 
lost_catch_by_year <- all_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, ssl_scenario, catch_alt) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  mutate(lost_catch = `No SSL control` - `SSL control`)

# cumulative lost catch 
lost_catch_cumulative <- lost_catch_by_year %>%
  group_by(harvest_rate) %>%
  arrange(Year) %>%
  mutate(cumulative_lost_catch = cumsum(lost_catch)) %>%
  ungroup()

(p1 <- ggplot(lost_catch_cumulative, aes(Year, cumulative_lost_catch, color = factor(harvest_rate))) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Cumulative Lost Catch (number of fish)",
       color = "Harvest rate") +
  scale_color_viridis_d() + 
  #scale_color_manual(values = palette.colors(n = 8, palette = "Okabe-Ito")[-1][1:7]) +
  theme_minimal()+
    theme(legend.position = "none") )


p2 <- ggplot(lost_catch_summary, aes(factor(harvest_rate), -total_lost_catch)) +
  geom_col(fill = "#440154") +
  labs(x = "Harvest rate", y = "Total lost catch since 1978 (No SSL control − SSL control)") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") 

p1 +  p2 + plot_layout(widths = c(2, 1))

p2 <- ggplot(lost_catch_summary, aes(factor(harvest_rate), -total_lost_catch, fill = factor(harvest_rate))) +
  geom_col() +
  labs(x = "Harvest rate", y = "Total lost catch since 1978") +
  scale_fill_viridis_d() + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none") 

p1 +  p2 + plot_layout(widths = c(2, 1))

ggsave(file = "figures/chum_lost_catch.png", height = 5, width = 10, dpi = 600)



## lost money historic U -------------------------------------
cpi_table <- read_csv("data/cpi_table_world-bank.csv")
value_per_fish <- read_csv("data/value-per-fish_2020_v2.csv")

cpi_2026 <- cpi_table %>%
  filter(Year == 2025) %>%
  pull(CPI)

# value added 
value_per_fish_table <- cpi_table %>%
  mutate(value_per_fish = 439 * (CPI / cpi_2026))

value_per_fish_table

# raw value 
value_per_fish_table_raw <- value_per_fish %>%
  mutate(value_per_fish = Chum * 5) %>%
  select(Year, value_per_fish)

value_per_fish_table_raw_2 <- cpi_table %>%
  mutate(value_per_fish = 3.50 * 11 * (CPI/cpi_2026)) %>%
  select(Year, value_per_fish)

#historic_scenarios <- historic_scenarios %>%
  #select(-value_per_fish_added) %>%
  #select(-value_per_fish_raw) %>%
#  select(-value_per_fish)


historic_scenarios <- historic_scenarios %>%
  left_join(value_per_fish_table_raw %>% 
              select(Year, value_per_fish) %>% 
              rename(value_per_fish_raw = value_per_fish),
            by = "Year") %>%
  left_join(value_per_fish_table %>% select(Year, value_per_fish) %>% 
              rename(value_per_fish_added = value_per_fish),
            by = "Year") %>%
  left_join(value_per_fish_table_raw_2 %>% select(Year, value_per_fish) %>% 
              rename(value_per_fish_raw_2 = value_per_fish),
            by = "Year")


# Lost catch and money year to year
lost_catch_by_year <- historic_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, ssl_scenario, catch_alt, value_per_fish_raw, value_per_fish_raw_2, value_per_fish_added) %>%
  pivot_wider(
    names_from = ssl_scenario,
    values_from = catch_alt
  ) %>%
  mutate(
    lost_catch = `No SSL control` - `SSL control`,
    lost_money_raw = lost_catch * value_per_fish_raw,
    lost_money_raw_2 = lost_catch * value_per_fish_raw_2,
    lost_money_added = lost_catch * value_per_fish_added
  )

# Cumulative lost catch and money
lost_catch_cumulative <- lost_catch_by_year %>%
  arrange(Year) %>%
  mutate(
    cumulative_lost_catch = cumsum(lost_catch),
    cumulative_lost_money_raw = cumsum(lost_money_raw),
    cumulative_lost_money_raw_2 = cumsum(lost_money_raw_2),
    cumulative_lost_money_added = cumsum(lost_money_added)
  )

# Total lost catch and money
lost_catch_summary <- lost_catch_by_year %>%
  summarise(
    total_catch_no_control = sum(`No SSL control`, na.rm = TRUE),
    total_catch_ssl_control = sum(`SSL control`, na.rm = TRUE),
    total_lost_catch = total_catch_no_control - total_catch_ssl_control,
    total_lost_money = sum(lost_money, na.rm = TRUE)
  )

# Cumulative lost money through time
ggplot(
  lost_catch_cumulative,
  aes(Year, -cumulative_lost_money)
) +
  geom_area(alpha = 0.2, fill = "#FF4500") +
  geom_line(linewidth = 1, color = "#FF4500") +
  labs(
    x = "Year",
    y = "Cumulative lost value (CAD)"
  ) +
  ylim(0,1500000000) +
  scale_y_continuous(labels = scales::label_dollar(scale = 1)) +
  theme_minimal()

# Cumulative lost money through time
lost_catch_cumulative %>%
  pivot_longer(
    cols = c(
      "cumulative_lost_money_raw",
      "cumulative_lost_money_added",
      "cumulative_lost_money_raw_2"
    ),
    names_to = "value",
    values_to = "lost_money"
  ) %>%
  mutate(
    value = recode(
      value,
      cumulative_lost_money_raw = "Raw Landed (historic, DFO)",
      cumulative_lost_money_added = "Added",
      cumulative_lost_money_raw_2 = "Raw Landed (current)"
    )
  ) %>%
  ggplot(aes(Year, lost_money, color = value, fill = value)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Year",
    y = "Cumulative lost value (CAD)",
    color = "Value ",
    fill = "Value"
  ) +
  scale_y_continuous(labels = scales::label_dollar(scale = 1)) +
  theme_minimal()

ggsave(file = "figures/chum_catch_money_Uhistoric_value_all.png", height = 5,width = 9, dpi = 600)

## lost money Uapply-------------------------------------
cpi_table <- read_csv("data/cpi_table_world-bank.csv")

cpi_2026 <- cpi_table %>% filter(Year == 2025) %>% pull(CPI)

value_per_fish_table <- cpi_table %>%
  mutate(value_per_fish = 439 * (CPI / cpi_2026))

value_per_fish_table

all_scenarios <- all_scenarios %>%
  left_join(value_per_fish_table %>% select(Year, value_per_fish), by = "Year")

# lost catch year to year 
lost_catch_by_year <- all_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, ssl_scenario, catch_alt, value_per_fish) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  mutate(lost_catch = `No SSL control` - `SSL control`,
         lost_money = lost_catch*value_per_fish)

# cumulative lost catch 
lost_catch_cumulative <- lost_catch_by_year %>%
  group_by(harvest_rate) %>%
  arrange(Year) %>%
  mutate(cumulative_lost_catch = cumsum(lost_catch),
         cumulative_lost_money = cumsum(lost_money)) %>%
  ungroup()

lost_catch_summary <- all_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  select(Year, harvest_rate, ssl_scenario, catch_alt, value_per_fish) %>%
  pivot_wider(names_from = ssl_scenario, values_from = catch_alt) %>%
  mutate(lost_catch = `No SSL control` - `SSL control`,
         lost_money = lost_catch * value_per_fish) %>%
  group_by(harvest_rate) %>%
  summarise(
    total_catch_no_control = sum(`No SSL control`, na.rm = TRUE),
    total_catch_ssl_control = sum(`SSL control`, na.rm = TRUE),
    total_lost_catch = total_catch_no_control - total_catch_ssl_control,
    total_lost_money = sum(lost_money, na.rm = TRUE),   # <-- add
    .groups = "drop"
  )


(p1 <- ggplot(lost_catch_cumulative, aes(Year, cumulative_lost_money, color = factor(harvest_rate))) +
    geom_line(linewidth = 1) +
    labs(x = "Year", y = "Cumulative Lost Money (CAD)",
         color = "Harvest rate") +
    scale_y_continuous(labels = label_dollar(scale = 1)) +
    scale_color_viridis_d() + 
    theme_minimal()+
    theme(legend.position = "none") )

(p2 <- ggplot(lost_catch_summary, aes(factor(harvest_rate), -total_lost_money, fill = factor(harvest_rate))) +
  geom_col() +
  labs(x = "Harvest rate", y = "Total lost money (CAD) since 1978") +
  scale_y_continuous(labels = label_dollar(scale = 1)) +
  scale_fill_viridis_d() + 
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "none"))

p1 +  p2 + plot_layout(widths = c(2, 1.5))

ggsave(file = "figures/chum_lost_money.png", height = 5, width = 10, dpi = 600)

# bycatch scenarios --------------------------------
bycatch_rates <- seq(0, 1, by = 0.1)

# run full sceanrio
run_full_scenario <- function(chum_data, covariates, sh_data,
                              U_apply, SSL_control, U_historic = 0, byrate = 0.69) {
    chum_df <- run_chum_scenario(
    chum_data,
    covariates,
    U_apply,
    SSL_control,
    U_historic
  )
  all_data <- chum_df %>% right_join(sh_data, by = "Year")
  thompson_df <- run_thompson_scenario(all_data, SSL_control, U_historic, byrate)
  full_df <- run_chilcotin_scenario(thompson_df, SSL_control, U_historic, byrate)
  
  full_df %>%
    mutate(
      byrate = byrate,
      ssl_scenario = if_else(
        SSL_control == 1,
        "SSL control",
        "No SSL control"
      )
    )
}

steelhead_scenarios <- purrr::pmap_dfr(
  expand_grid(
    SSL_control = c(0, 1),
    byrate = bycatch_rates
  ),
  function(SSL_control, byrate) {
    run_full_scenario(
      chum_data = data,
      covariates = covariates,
      sh_data = sh_data,
      U_apply = 0,
      SSL_control = SSL_control,
      U_historic = 1,
      byrate = byrate
    )
  }
)

# check 
steelhead_scenarios %>%
  distinct(byrate, ssl_scenario)

steelhead_productivity <- steelhead_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  mutate(
    lnRS_thompson = log(sh_thompson_recruits_alt / sh_thompson_spawners_pred)
  )

steelhead_productivity <- steelhead_scenarios %>%
  filter(Year >= 1978, Year < 2017) %>%
  mutate(
    lnRS_thompson = log(
      sh_thompson_recruits_alt / sh_thompson_spawners_pred
    ),
    lnRS_chilcotin = log(
      sh_chilcotin_recruits_alt / sh_chilcotin_spawners_pred
    )
  ) %>%
  select(
    Year,
    byrate,
    ssl_scenario,
    lnRS_thompson,
    lnRS_chilcotin
  ) %>%
  pivot_longer(
    cols = c(lnRS_thompson, lnRS_chilcotin),
    names_to = "Stock",
    values_to = "lnRS"
  ) %>%
  mutate(
    Stock = recode(
      Stock,
      lnRS_thompson = "Thompson",
      lnRS_chilcotin = "Chilcotin"
    )
  )

ggplot(
  steelhead_productivity,
  aes(
    x = Year,
    y = lnRS,
    color = factor(byrate),
    linetype = ssl_scenario
  )
) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(
    x = "Year",
    y = "ln(R/S)",
    color = "Bycatch rate",
    linetype = "SSL control"
  ) +
  theme_minimal()

### plots together -------------
make_scenario_plot <- function(df, hr, scenario) {
  d <- df %>% filter(harvest_rate == hr, ssl_scenario == scenario)
  
  ggplot(d, aes(Year, value)) +
    geom_line(color = "#5E4FA2", linewidth = 1) +
    facet_wrap(~metric, ncol = 1, scales = "free_y") +
    labs(title = paste0(scenario, " — U = ", hr), x = NULL, y = NULL) +
    theme_minimal(base_size = 9)
}

plots <- scenario_grid %>%
  mutate(ssl_scenario = if_else(SSL_control == 1, "SSL control", "No SSL control")) %>%
  mutate(plot = purrr::map2(U_apply, ssl_scenario, ~make_scenario_plot(plot_data, .x, .y)))

no_control_plots <- plots %>% filter(ssl_scenario == "No SSL control") %>% pull(plot)
control_plots     <- plots %>% filter(ssl_scenario == "SSL control") %>% pull(plot)

wrap_plots(c(no_control_plots, control_plots), ncol = length(harvest_rates))


make_prod_plot <- function(df, hr, scenario, yrange) {
  d <- df %>% filter(harvest_rate == hr, ssl_scenario == scenario, metric == "Productivity (alpha)")
  ggplot(d, aes(Year, value)) +
    geom_line(color = "#5E4FA2", linewidth = 1) +
    ylim(yrange) +
    labs(title = paste0(scenario, " — U = ", hr), x = NULL, y = "alpha") +
    theme_minimal(base_size = 9)
}

make_returns_plot <- function(df, hr, scenario, yrange) {
  d <- df %>% filter(harvest_rate == hr, ssl_scenario == scenario, metric == "Recruits/returns")
  ggplot(d, aes(Year, value)) +
    geom_line(color = "#D55E00", linewidth = 1) +
    ylim(yrange) +
    labs(title = paste0(scenario, " — U = ", hr), x = NULL, y = "recruits") +
    theme_minimal(base_size = 9)
}
wrap_plots(c(no_control_prod_plots, control_prod_plots), ncol = length(harvest_rates)) +
  plot_annotation(title = "Productivity across chum harvest rate scenarios")




















