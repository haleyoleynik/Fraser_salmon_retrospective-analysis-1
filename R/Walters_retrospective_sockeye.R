# Sockeye retrospective model from Carl Walters excel model 
# January 2026
# Haley Oleynik 

#   * Observed exploitation (Ut) and ENS
#   * Fitting ln(R/S) = ra - rb * S
#   * Retrospective predictionwith 4-year lag using historical residuals

# load packages 
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
require(ggplot2)
require(patchwork)

# Configure  -----------------------------

# Years used in the Excel regression ranges:
FIT_YEARS  <- 1952:2019

# Lag in rows for R/S in the spreadsheet (X4 = O8, i.e., +4 rows)
LAG_YEARS  <- 4

# Read observed data -----------------------------

# INPUT_XLSX <- "Walters_model_simplified.xlsx"
# SHEET      <- "Sheet1"
# 
# # The main data table starts at row 3 (headers) and runs through year 2023.
# obs <- read_excel(
#   path      = INPUT_XLSX,
#   sheet     = SHEET,
#   range     = "A3:I75",
#   col_names = TRUE
# ) %>%
#   rename(
#     Year            = `Year`,
#     AdultEscapement = `Sum of Adult Escapement`,
#     JackEscapement  = `Sum of Jack Escapement`,
#     TotalEscapement = `Sum of Total Escapement`,
#     DBE             = `Sum of DBE`,
#     BelowMissionC   = `Sum of Below Mission Catch (excluding Alaska catch)`,
#     AboveMissionC   = `Sum of Above Mission Catch`,
#     AlaskaCatch     = `Sum of Alaska Catch`,
#     RunSize         = `Sum of Run Size`
#   ) %>%
#   mutate(Year = as.integer(Year)) %>%
#   filter(!is.na(Year)) %>%
#   arrange(Year)

# new read observed data with all stocks 
obs <- read_csv("R/Sockeye Retrospective Shiny App/Walters_model_all-stocks.csv")

obs <- obs %>%
  rename(
    Year            = `Year`,
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    DBE             = `DBE`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  arrange(Stock,Year)

# Set parameters ---------------------------

retroU <- 0.3
useretro <- 1 # 1 for true, 0 for false 
yrretro <- 1990 


#  Reproduce the observed-derived columns ---------------
obs2 <- obs %>%
  mutate(
    RunJacks = RunSize - JackEscapement,  
    
    # sum catches but ignore NAs
    Catch = rowSums(
      cbind(BelowMissionC, AboveMissionC),
      na.rm = TRUE
    ),
    
    Ut_obs  = pmin(0.999, Catch / RunJacks),
    ENS     = pmin(1, pmax(0.0001,
                           AdultEscapement / RunJacks / (1 - Ut_obs))),
    migmort = 1 - ENS
  ) %>%
  mutate(
    AdultReturn = lead(RunJacks, n = LAG_YEARS),
    lnR_S       = log(AdultReturn / AdultEscapement)
  )

# Fit a & b parameters ----------------------------------------
# using 1952-2019
fit_df <- obs2 %>%
  filter(Year %in% FIT_YEARS) %>%
  filter(is.finite(lnR_S), is.finite(AdultEscapement))

fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)

ra <- unname(coef(fit)[[1]])
rb <- unname(-coef(fit)[[2]])  

obs3 <- obs2 %>%
  mutate(
    wt = lnR_S - (ra - rb * AdultEscapement)       # residuals
  )

# Retrospective predictions  -----------------------------

n <- nrow(obs3)

retroR <- rep(NA_real_, n)  # AA
retroU_vec <- if (useretro) ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs) else obs3$Ut_obs  # AB
retroS <- rep(NA_real_, n)  # AC
retroC <- rep(NA_real_, n)  # AD

# Seed with observed RunJacks for first LAG_YEARS years
retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]

# Compute spawners/catch for those seeded rows
retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * (1 - retroU_vec[1:LAG_YEARS]) * obs3$ENS[1:LAG_YEARS]
retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * retroU_vec[1:LAG_YEARS]

# Iterate forward
for (i in (LAG_YEARS + 1):n) {
  j <- i - LAG_YEARS
  # recruit/run from spawners 4 years prior
  retroR[i] <- retroS[j] * exp(ra - rb * retroS[j] + obs3$wt[j])
  retroS[i] <- retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i]
  retroC[i] <- retroR[i] * retroU_vec[i]
}

out <- obs3 %>%
  transmute(
    Year,
    AdultEscapement,
    JackEscapement,
    TotalEscapement,
    DBE,
    BelowMissionC,
    AboveMissionC,
    AlaskaCatch,
    RunSize,
    RunJacks,
    Catch,
    Ut_obs,
    ENS,
    migmort,
    AdultReturn,
    lnR_S,
    wt,
    retroR = retroR,
    retroU = retroU_vec,
    retroS = retroS,
    retroC = retroC
  )


hist_catch <- sum(out$Catch, na.rm=T)
retro_catch <- sum(out$retroC, na.rm=T)
lost_Catch <- retro_catch-hist_catch

# Save outputs -----------------------------

# message("Fitted parameters:")
# message(sprintf("  ra = %.8f", ra))
# message(sprintf("  rb = %.12f", rb))
# message("Retro toggles:")
# message(sprintf("  useretro = %s", useretro))
# message(sprintf("  retroU   = %.4f", retroU))
# message(sprintf("  yrretro  = %d", yrretro))
# 
# write.csv(out, "walters_model_outputs.csv", row.names = FALSE)
# message("Wrote: walters_model_outputs.csv")

# NEW retro with covariates -----------------------------------------------

# read in standardized covariates 
covariates <- read_csv("Data/sockeye_standardized_covariates.csv")
dredge_models <- read_csv("Data/sockeye_top_models_dredge_wo-aquaculture.csv")

dredge_models_wild_pink <- read_csv("Data/sockeye_top_models_dredge_wo-aquaculture_wild-pink.csv")

## configure ----------------------
FIT_YEARS <- 1952:2019
LAG_YEARS <- 4

#  Parameters that were previously Shiny inputs
retroU       <- 0.3          # retrospective harvest rate
useretro     <- TRUE         # use retrospective harvest rate?
yrretro      <- 1990         # retrospective start year
stock_choice <- "All Stocks" # "All Stocks" or a specific Stock name

## load data -----------------------------
obs <- read_csv("R/Sockeye Retrospective Shiny App/Walters_model_all-stocks.csv") %>%
  rename(
    Year            = `Year`,
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    DBE             = `DBE`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  arrange(Stock, Year) %>%
  filter(Stock %in% c( "Birkenhead" ,
                       "Bowron"   ,
                       "Chilko"    ,
                       "Cultus"     ,
                       "Early Stuart" ,
                       "Gates"  ,
                       "Late Shuswap" ,
                       "Late Stuart" ,
                       "Pitt"   ,
                       "Portage"    ,
                       "Quesnel"   ,
                       "Raft"      ,
                       "Scotch"    ,
                       "Seymour"    ,
                       "Stellako"   ,
                       "Weaver"  ))

# # Per-stock top models from dredge/model-selection output
# top_models <- dredge_models %>%
#   group_by(Stock) %>%
#   slice_min(deltaAIC, n = 1, with_ties = FALSE) %>%
#   ungroup()
# 
# # Candidate covariate columns that can appear in a stock's top model
# COVARIATE_COLS <- c("NPGO", "PDO", "SeaLions", "seal",
#                     "adult.sst", "pink", "smolt.sst")

# specify covariates 
COVARIATE_COLS <- c("NPGO", "PDO", "SeaLions", "seal",
                    "adult.sst", "pink", "smolt.sst")

# AIC weighted coefficients for each covariate 
top_models <- dredge_models %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )

top_models_pink <- dredge_models_wild_pink %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )


# Yearly covariate time series, keyed by Stock + Year
covariates <- read_csv("Data/sockeye_standardized_covariates.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, any_of(COVARIATE_COLS))

obs <- obs %>%
  left_join(covariates, by = c("Stock", "Year"))

# Pinniped control scenario -------------------------------------
FREEZE_YEAR <- 1970
FREEZE_VARS <- c("SeaLions", "seal")

no_freeze_year <- covariates %>%
  group_by(Stock) %>%
  summarise(has_freeze_year = any(Year == FREEZE_YEAR & if_any(all_of(FREEZE_VARS), ~ !is.na(.))),
            .groups = "drop") %>%
  filter(!has_freeze_year)

if (nrow(no_freeze_year) > 0) {
  warning("No ", FREEZE_YEAR, " SeaLions/seal value for stock(s): ",
          paste(no_freeze_year$Stock, collapse = ", "),
          " -- scenario will equal actual for these stocks.")
}

scenario_covariates <- covariates %>%
  group_by(Stock) %>%
  mutate(across(all_of(FREEZE_VARS), function(x) {
    freeze_val <- x[Year == FREEZE_YEAR]
    if (length(freeze_val) != 1 || is.na(freeze_val)) {
      x
    } else {
      ifelse(Year > FREEZE_YEAR, freeze_val, x)
    }
  })) %>%
  ungroup() %>%
  select(Stock, Year, all_of(FREEZE_VARS)) %>%
  rename_with(~ paste0(., "_scenario"), all_of(FREEZE_VARS))

obs_with_scenario <- obs %>%
  left_join(scenario_covariates, by = c("Stock", "Year"))

#  TOP-MODEL / COVARIATE HELPERS

# Look up a stock's top-model coefficients from the dredge table.
get_top_model_terms <- function(stock_name) {
  top_model_row <- top_models %>% filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
    list(
      ra            = top_model_row[["(Intercept)"]],
      rb            = -top_model_row[["spawners"]],
      sel_covs      = sel_covs,
      cov_coefs     = if (length(sel_covs) > 0) as.numeric(top_model_row[sel_covs]) else numeric(0),
      has_top_model = TRUE
    )
  } else {
    list(ra = NA_real_, rb = NA_real_, sel_covs = character(0),
         cov_coefs = numeric(0), has_top_model = FALSE)
  }
}

# Linear-predictor contribution from the selected covariates.
# use_scenario = TRUE pulls the "<var>_scenario" (pinniped-frozen) column
# for any covariate in FREEZE_VARS instead of its actual observed column.
compute_cov_term <- function(dat, sel_covs, cov_coefs, use_scenario = FALSE) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  cols <- if (use_scenario) {
    ifelse(sel_covs %in% FREEZE_VARS, paste0(sel_covs, "_scenario"), sel_covs)
  } else {
    sel_covs
  }
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}

#  RETRO FUNCTION
run_retro_model <- function(dat,
                            stock_name,
                            retroU,
                            useretro,
                            yrretro,
                            run_scenario = FALSE) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs = pmin(0.95, Catch / RunJacks),
      Ut_obs = if_else(
        stock_name == "Late Shuswap" & Year == 2012,  # manually adjust late shuswap 2012 to match historical U 
        pmin(0.7, Catch / RunJacks),
        Ut_obs
      ),
      ENS = pmin(1, pmax(0.0001,
                         AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  terms <- get_top_model_terms(stock_name)
  
  if (terms$has_top_model) {
    ra <- terms$ra
    rb <- terms$rb
    sel_covs <- terms$sel_covs
    cov_coefs <- terms$cov_coefs
  } else {
    # Fallback: no dredge top model for this stock -> plain Ricker fit
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
    ra <- coef(fit)[1]
    rb <- -coef(fit)[2]
    sel_covs <- character(0)
    cov_coefs <- numeric(0)
  }
  
  # cov_term (actual/historical covariates) always drives wt, the estimated
  # process-error residual -- this stays fixed across scenarios so that a
  # scenario comparison isolates the covariate effect rather than mixing it
  # with a different noise realization.
  cov_term <- compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = FALSE)
  
  # cov_term_proj drives the forward recursive projection below, and is the
  # only thing that differs between the actual run and the pinniped scenario.
  cov_term_proj <- if (run_scenario) {
    compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = TRUE)
  } else {
    cov_term
  }
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (ra - rb * AdultEscapement + cov_term)
    )
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroU_vec <- if (useretro)
    ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
  else
    obs3$Ut_obs
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- ra - rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR     = retroR,
      retro_lnRS = retro_lnRS,
      retroU     = retroU_vec,
      retroS     = retroS,
      retroC     = retroC
    )
}

# RUN MODEL (all stocks
model_all <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()


# Pinniped scenario run: SeaLions/seal frozen at FREEZE_YEAR levels
model_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

# FILTER BY STOCK
filtered_model <- if (stock_choice == "All Stocks") {
  model_all
} else {
  model_all %>% filter(Stock == stock_choice)
}

#  STOCK SUMMARY
stock_table <- if (stock_choice == "All Stocks") {
  filtered_model %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch,
      .groups = "drop"
    )
} else {
  filtered_model %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch
    )
}
print(stock_table)

# TOTAL SUMMARY
hist_total  <- sum(filtered_model$Catch,  na.rm = TRUE)
retro_total <- sum(filtered_model$retroC, na.rm = TRUE)

cat("Historical Catch: ",
    format(round(hist_total, 0), big.mark = ","), "\n")
cat("Retrospective Catch: ",
    format(round(retro_total, 0), big.mark = ","), "\n")
cat("Difference: ",
    format(round(retro_total - hist_total, 0), big.mark = ","), "\n")

# AGGREGATE TIME SERIES (summed across stocks
summed_ts <- filtered_model %>%
  group_by(Year) %>%
  summarise(
    Catch = sum(Catch, na.rm = TRUE),
    retroC = sum(retroC, na.rm = TRUE),
    RunJacks = sum(RunJacks, na.rm = TRUE),
    retroR = sum(retroR, na.rm = TRUE),
    AdultEscapement = sum(AdultEscapement, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year) %>%
  mutate(
    AdultReturn = lead(RunJacks, LAG_YEARS)
  )

#  RETRO TIME SERIES TABLE
retro_ts_table <- summed_ts %>%
  select(Year, retroR, retroC, Run = RunJacks) %>%
  arrange(Year) %>%
  mutate(
    retroR = round(retroR, 0),
    retroC = round(retroC, 0),
    Run = round(Run, 0)
  )
print(retro_ts_table)

## plot retro catch, etc. ------------------------

# Catch plot
catch_plot <- ggplot(
  summed_ts %>% select(Year, Catch, retroC) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Catch", y = "Catch", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Adult return plot
return_plot <- ggplot(
  summed_ts %>% select(Year, AdultReturns = RunJacks, retroR) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Return", y = "Returns", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Escapement plot
esc_plot <- ggplot(
  summed_ts %>% select(Year, AdultEscapement, retroS) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Spawners", y = "Spawners", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

print(catch_plot)
print(return_plot)
print(esc_plot)

# ggsave("catch_plot.png", catch_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("return_plot.png", return_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("esc_plot.png", esc_plot, width = 8, height = 4.5, dpi = 300)

## pinniped control scenario comparison ---------------------------

filtered_scenario <- if (stock_choice == "All Stocks") {
  model_scenario
} else {
  model_scenario %>% filter(Stock == stock_choice)
}

# Per-stock, per-year comparison
scenario_compare <- filtered_model %>%
  select(Stock, Year, lnRS_actual = retro_lnRS, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    filtered_scenario %>%
      select(Stock, Year, lnRS_scenario = retro_lnRS, R_scenario = retroR, S_scenario = retroS),
    by = c("Stock", "Year")
  ) %>%
  mutate(
    lnRS_diff = lnRS_scenario - lnRS_actual,
    R_diff    = R_scenario - R_actual,
    S_diff    = S_scenario - S_actual
  )
print(scenario_compare)

# Aggregate (summed across stocks) comparison, for plotting
summed_scenario_ts <- filtered_scenario %>%
  group_by(Year) %>%
  summarise(
    retroR = sum(retroR, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year)

compare_ts <- summed_ts %>%
  select(Year, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    summed_scenario_ts %>% select(Year, R_scenario = retroR, S_scenario = retroS),
    by = "Year"
  )

# Reconstructed returns: actual vs pinniped-frozen scenario
pinniped_scenario_return_plot <- ggplot(
  compare_ts %>% select(Year, R_actual, R_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "pinnipeds",
       y = "Returns", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Reconstructed spawners: actual vs pinniped-frozen scenario
scenario_esc_plot <- ggplot(
  compare_ts %>% select(Year, S_actual, S_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(
       y = "Spawners", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

print(pinniped_scenario_return_plot)
print(scenario_esc_plot)



## productivity compare -------------------------------------
productivity_compare <- model_all %>%
  select(Stock, Year, lnRS_actual = retro_lnRS) %>%
  left_join(
    model_scenario %>% select(Stock, Year, lnRS_scenario = retro_lnRS),
    by = c("Stock", "Year")
  ) %>%
  pivot_longer(cols = c(lnRS_actual, lnRS_scenario),
               names_to = "scenario", values_to = "lnRS") %>%
  mutate(scenario = recode(scenario,
                           lnRS_actual   = "Actual",
                           lnRS_scenario = "Pinniped controls"))

ggplot(productivity_compare, aes(Year, lnRS, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_productivity.png", dpi = 600, width = 7, height = 10)

productivity_compare %>%
filter(Stock %in% c("Birkenhead", "Bowron", "Chilko","Cultus", "Early Stuart", "Raft", "Stellako")) %>%
  ggplot(aes(Year, lnRS, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_productivity_select-stocks.png", dpi = 600, width = 10, height = 6)

## catch lost by harvest rate scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.05)

# lost_catch = retro_catch - hist_catch (positive = more catch under the
# retro harvest rate than actually occurred; flip the sign if you want
# "catch lost" to read positive when the retro scenario catches less)
catch_lost_by_hr <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch, na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(lost_catch = retro_catch - hist_catch,
           retroU     = hr)
})

# Per-stock catch lost vs harvest rate
ggplot(catch_lost_by_hr, aes(retroU, lost_catch)) +
  geom_line(linewidth = 1, color = "#4682B4") +
  geom_point(size = 1.5, color = "#4682B4") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Catch lost by retrospective harvest rate, per stock",
       x = "Retrospective harvest rate", y = "Catch lost (retro - historical)") +
  theme_minimal()

# Cumulative catch lost across all stocks
catch_lost_cumulative <- catch_lost_by_hr %>%
  group_by(retroU) %>%
  summarise(total_lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop")

ggplot(catch_lost_cumulative, aes(retroU, total_lost_catch)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  geom_point(size = 2, color = "#FF4500") +
  labs(title = "Cumulative catch lost across all stocks, by harvest rate",
       x = "Retrospective harvest rate", y = "Total catch lost") +
  theme_minimal()

## catch lost by year scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.1)  # tweak step size as you like

catch_lost_ts <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    mutate(lost_catch = retroC - Catch,
           retroU     = hr) %>%
    select(Stock, Year, retroU, lost_catch)
})

# Cumulative catch lost over time, by stock
catch_lost_cum_by_stock <- catch_lost_ts %>%
  arrange(Stock, retroU, Year) %>%
  group_by(Stock, retroU) %>%
  mutate(cum_lost_catch = cumsum(replace_na(lost_catch, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Cumulative catch lost over time, by stock",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

# Cumulative catch lost over time, total across all stocks
catch_lost_cum_total <- catch_lost_ts %>%
  group_by(retroU, Year) %>%
  summarise(lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop") %>%
  arrange(retroU, Year) %>%
  group_by(retroU) %>%
  mutate(cum_lost_catch = cumsum(lost_catch)) %>%
  ungroup()

ggplot(catch_lost_cum_total, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cumulative catch lost over time, total across stocks",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

## what's the question 

#  Actual harvest rate, actual pinniped abundance -----------------
model_actual_harvest <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,       # unused when useretro = FALSE
                                 useretro     = FALSE,   # use actual historical harvest rate
                                 yrretro      = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()

# Actual harvest rate, SL/seal frozen at FREEZE_YEAR
model_actual_harvest_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,
                                 useretro     = FALSE,
                                 yrretro      = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

#  Catch lost to pinnipeds: scenario catch minus actual catch
catch_lost_to_pinnipeds <- model_actual_harvest %>%
  select(Stock, Year, catch_actual = retroC) %>%
  left_join(
    model_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
    by = c("Stock", "Year")
  ) %>%
  mutate(catch_lost = catch_scenario - catch_actual)  # positive = catch lost to pinnipeds

#  Cumulative catch lost to pinnipeds over time, by stock
catch_lost_cum_by_stock <- catch_lost_to_pinnipeds %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# filter out stocks that weren't impacted 
p1 <- catch_lost_cum_by_stock %>%
  filter(Stock %in% c("Birkenhead", "Bowron", "Chilko","Cultus", "Early Stuart", "Raft", "Stellako")) %>%
ggplot(aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, 
             #scales = "free_y"
             ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# Cumulative catch lost to pinnipeds over time, total across stocks
catch_lost_cum_total <- catch_lost_to_pinnipeds %>%
  group_by(Year) %>%
  summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cum_catch_lost = cumsum(catch_lost))

p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p1 | p2

ggsave("figures/sockeye_lost_catch.png", dpi = 600, width = 10, height = 6)


# NEW PINK retro ---------------------------------------------------------
# if hatchery production had not happened (so using wild only)

covariates <- read_csv("Data/sockeye_standardized_covariates_pink-wild.csv") %>%
  rename(Year = yr)

covariates %>%
  select(Stock, Year, pink, pink_wild) %>%
  pivot_longer(c(pink, pink_wild)) %>%
  ggplot(aes(Year, value, color = name)) +
  geom_line() +
  facet_wrap(~ Stock) +
  theme_minimal()

covariates %>%
  filter(Stock == "Late Stuart") %>%
  select(Stock, Year, pink, pink_wild) %>%
  pivot_longer(c(pink, pink_wild)) %>%
  ggplot(aes(Year, value, color = name)) +
  geom_line(size = 1) +
  labs(y = "Pink Abundance (standardized)") +
  theme_minimal()

ggsave("figures/pink-vs-pink_wild.png", dpi = 600, width = 10, height = 6)

dredge_models <- read_csv("Data/sockeye_top_models_dredge_wo-aquaculture.csv")

# only top (0 deltaAIC)
#top_models <- dredge_models %>%
#  group_by(Stock) %>%
#  slice_min(deltaAIC, n = 1, with_ties = FALSE) %>%
#  ungroup()

# AIC weighted coefficients for each covariate
top_models <- dredge_models %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )


FREEZE_VARS <- c("pink")

# alternative pink salmon timeseries 
no_alt_series <- covariates %>%
  group_by(Stock) %>%
  summarise(has_alt = any(!is.na(pink_wild)), .groups = "drop") %>%
  filter(!has_alt)

if (nrow(no_alt_series) > 0) {
  warning("No pink_wild values for stock(s): ",
          paste(no_alt_series$Stock, collapse = ", "),
          " -- scenario will equal actual for these stocks.")
}


scenario_covariates <- covariates %>%
  select(Stock, Year, pink_scenario = pink_wild)

obs_with_scenario <- obs %>%
  left_join(scenario_covariates, by = c("Stock", "Year"))

#  TOP-MODEL / COVARIATE HELPERS
# Look up a stock's top-model coefficients from the dredge table.
get_top_model_terms <- function(stock_name) {
  top_model_row <- top_models %>% filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
    list(
      ra            = top_model_row[["(Intercept)"]],
      rb            = -top_model_row[["spawners"]],
      sel_covs      = sel_covs,
      cov_coefs     = if (length(sel_covs) > 0) as.numeric(top_model_row[sel_covs]) else numeric(0),
      has_top_model = TRUE
    )
  } else {
    list(ra = NA_real_, rb = NA_real_, sel_covs = character(0),
         cov_coefs = numeric(0), has_top_model = FALSE)
  }
}

# Linear-predictor contribution from the selected covariates.
# use_scenario = TRUE pulls the "<var>_scenario" column
# for any covariate in FREEZE_VARS instead of its actual observed column.
compute_cov_term <- function(dat, sel_covs, cov_coefs, use_scenario = FALSE) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  cols <- if (use_scenario) {
    ifelse(sel_covs %in% FREEZE_VARS, paste0(sel_covs, "_scenario"), sel_covs)
  } else {
    sel_covs
  }
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}


# RUN MODEL (all stocks
model_all_pink <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()


# pink scenario run: pink salmon frozen at FREEZE_YEAR levels
model_scenario_pink <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

# filter by stock
filtered_model <- if (stock_choice == "All Stocks") {
  model_all
} else {
  model_all %>% filter(Stock == stock_choice)
}

#  stock summary
stock_table <- if (stock_choice == "All Stocks") {
  filtered_model %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch,
      .groups = "drop"
    )
} else {
  filtered_model %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch
    )
}
print(stock_table)

# TOTAL SUMMARY
hist_total  <- sum(filtered_model$Catch,  na.rm = TRUE)
retro_total <- sum(filtered_model$retroC, na.rm = TRUE)

cat("Historical Catch: ",
    format(round(hist_total, 0), big.mark = ","), "\n")
cat("Retrospective Catch: ",
    format(round(retro_total, 0), big.mark = ","), "\n")
cat("Difference: ",
    format(round(retro_total - hist_total, 0), big.mark = ","), "\n")

# AGGREGATE TIME SERIES (summed across stocks
summed_ts <- filtered_model %>%
  group_by(Year) %>%
  summarise(
    Catch = sum(Catch, na.rm = TRUE),
    retroC = sum(retroC, na.rm = TRUE),
    RunJacks = sum(RunJacks, na.rm = TRUE),
    retroR = sum(retroR, na.rm = TRUE),
    AdultEscapement = sum(AdultEscapement, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year) %>%
  mutate(
    AdultReturn = lead(RunJacks, LAG_YEARS)
  )

#  RETRO TIME SERIES TABLE
retro_ts_table <- summed_ts %>%
  select(Year, retroR, retroC, Run = RunJacks) %>%
  arrange(Year) %>%
  mutate(
    retroR = round(retroR, 0),
    retroC = round(retroC, 0),
    Run = round(Run, 0)
  )
print(retro_ts_table)

## plot retro catch, etc. ------------------------

# Catch plot
catch_plot <- ggplot(
  summed_ts %>% select(Year, Catch, retroC) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Catch", y = "Catch", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Adult return plot
return_plot <- ggplot(
  summed_ts %>% select(Year, AdultReturns = RunJacks, retroR) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Return", y = "Returns", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Escapement plot
esc_plot <- ggplot(
  summed_ts %>% select(Year, AdultEscapement, retroS) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Spawners", y = "Spawners", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

print(catch_plot)
print(return_plot)
print(esc_plot)

# ggsave("catch_plot.png", catch_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("return_plot.png", return_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("esc_plot.png", esc_plot, width = 8, height = 4.5, dpi = 300)

## pink control scenario comparison ---------------------------

filtered_scenario <- if (stock_choice == "All Stocks") {
  model_scenario
} else {
  model_scenario %>% filter(Stock == stock_choice)
}

# Per-stock, per-year comparison
scenario_compare <- filtered_model %>%
  select(Stock, Year, lnRS_actual = retro_lnRS, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    filtered_scenario %>%
      select(Stock, Year, lnRS_scenario = retro_lnRS, R_scenario = retroR, S_scenario = retroS),
    by = c("Stock", "Year")
  ) %>%
  mutate(
    lnRS_diff = lnRS_scenario - lnRS_actual,
    R_diff    = R_scenario - R_actual,
    S_diff    = S_scenario - S_actual
  )
print(scenario_compare)

# Aggregate (summed across stocks) comparison, for plotting
summed_scenario_ts <- filtered_scenario %>%
  group_by(Year) %>%
  summarise(
    retroR = sum(retroR, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year)

compare_ts <- summed_ts %>%
  select(Year, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    summed_scenario_ts %>% select(Year, R_scenario = retroR, S_scenario = retroS),
    by = "Year"
  )

# Reconstructed returns: actual vs pinniped-frozen scenario
pink_scenario_return_plot <- ggplot(
  compare_ts %>% select(Year, R_actual, R_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(y = "Returns", color = "", title = "pink salmon") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Reconstructed spawners: actual vs pinniped-frozen scenario
scenario_esc_plot <- ggplot(
  compare_ts %>% select(Year, S_actual, S_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(y = "Spawners", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

print(pink_scenario_return_plot)
print(scenario_esc_plot)

pinniped_scenario_return_plot + theme(legend.position = "none") | pink_scenario_return_plot + ylab("")

ggsave("figures/scenario_return_plot.png", width = 12, height = 6, dpi = 600)

# ggsave("scenario_return_plot.png", scenario_return_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("scenario_esc_plot.png", scenario_esc_plot, width = 8, height = 4.5, dpi = 300)

## productivity compare -------------------------------------
productivity_compare <- model_all %>%
  select(Stock, Year, lnRS_actual = retro_lnRS) %>%
  left_join(
    model_scenario %>% select(Stock, Year, lnRS_scenario = retro_lnRS),
    by = c("Stock", "Year")
  ) %>%
  pivot_longer(cols = c(lnRS_actual, lnRS_scenario),
               names_to = "scenario", values_to = "lnRS") %>%
  mutate(scenario = recode(scenario,
                           lnRS_actual   = "Actual",
                           lnRS_scenario = "pink salmon controlled"))

ggplot(productivity_compare, aes(Year, lnRS, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#ggsave("figures/sockeye_productivity_pink.png", dpi = 600, width = 7, height = 10)

productivity_compare %>%
  filter(Stock %in% c("Stellako", "Bowron", "Seymour", "Late Stuart", "Scotch", "Weaver", "Portage")) %>%
  ggplot(aes(Year, lnRS, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_productivity_select-stocks_pink.png", dpi = 600, width = 10, height = 6)

## catch lost by harvest rate scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.05)

# lost_catch = retro_catch - hist_catch (positive = more catch under the
# retro harvest rate than actually occurred; flip the sign if you want
# "catch lost" to read positive when the retro scenario catches less)
catch_lost_by_hr <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch, na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(lost_catch = retro_catch - hist_catch,
           retroU     = hr)
})

# Per-stock catch lost vs harvest rate
ggplot(catch_lost_by_hr, aes(retroU, lost_catch)) +
  geom_line(linewidth = 1, color = "#4682B4") +
  geom_point(size = 1.5, color = "#4682B4") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Catch lost by retrospective harvest rate, per stock",
       x = "Retrospective harvest rate", y = "Catch lost (retro - historical)") +
  theme_minimal()

# Cumulative catch lost across all stocks
catch_lost_cumulative <- catch_lost_by_hr %>%
  group_by(retroU) %>%
  summarise(total_lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop")

ggplot(catch_lost_cumulative, aes(retroU, total_lost_catch)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  geom_point(size = 2, color = "#FF4500") +
  labs(title = "Cumulative catch lost across all stocks, by harvest rate",
       x = "Retrospective harvest rate", y = "Total catch lost") +
  theme_minimal()

## catch lost by year scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.1)  # tweak step size as you like

catch_lost_ts <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    mutate(lost_catch = retroC - Catch,
           retroU     = hr) %>%
    select(Stock, Year, retroU, lost_catch)
})

# Cumulative catch lost over time, by stock
catch_lost_cum_by_stock <- catch_lost_ts %>%
  arrange(Stock, retroU, Year) %>%
  group_by(Stock, retroU) %>%
  mutate(cum_lost_catch = cumsum(replace_na(lost_catch, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Cumulative catch lost over time, by stock",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

# Cumulative catch lost over time, total across all stocks
catch_lost_cum_total <- catch_lost_ts %>%
  group_by(retroU, Year) %>%
  summarise(lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop") %>%
  arrange(retroU, Year) %>%
  group_by(retroU) %>%
  mutate(cum_lost_catch = cumsum(lost_catch)) %>%
  ungroup()

ggplot(catch_lost_cum_total, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cumulative catch lost over time, total across stocks",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

## what's the question 

##  Actual harvest rate, actual  -----------------
pink_actual_harvest <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,       # unused when useretro = FALSE
                                 useretro     = FALSE,   # use actual historical harvest rate
                                 yrretro      = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()

# Actual harvest rate, pink salmon frozen at FREEZE_YEAR
pink_actual_harvest_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,
                                 useretro     = FALSE,
                                 yrretro      = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

#  Catch lost to pinks: scenario catch minus actual catch
catch_lost_to_pink <- pink_actual_harvest %>%
  select(Stock, Year, catch_actual = retroC) %>%
  left_join(
    pink_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
    by = c("Stock", "Year")
  ) %>%
  mutate(catch_lost = catch_scenario - catch_actual)  # positive = catch lost to pinnipeds

#  Cumulative catch lost to pinnipeds over time, by stock
catch_lost_cum_by_stock_pink <- catch_lost_to_pink %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock_pink, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# filter out stocks that weren't impacted 
p1 <- catch_lost_cum_by_stock %>%
  filter(Stock %in% c("Seymour", "Bowron", "Late Stuart","Portage", "Weaver", "Stellako")) %>%
  ggplot(aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, 
             #scales = "free_y"
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# Cumulative catch lost to pinks over time, total across stocks
catch_lost_cum_total_pinks <- catch_lost_to_pink %>%
  group_by(Year) %>%
  summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cum_catch_lost = cumsum(catch_lost))

p2 <- ggplot(catch_lost_cum_total_pinks, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p2 <- ggplot(catch_lost_cum_total_pinks, aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p1 | p2

ggsave("figures/sockeye_lost_catch_pink.png", dpi = 600, width = 10, height = 6)


# exploration -----------------------
# distribution of run sizes 

all.sockeye.st %>%
  ggplot(aes(x= log(recruits), color = Stock, fill = Stock)) + 
  geom_density(alpha = 0.5) +
  theme_minimal()

covariates %>%
  ggplot(aes(x = pink, y = lnrs, color = Stock, fill = Stock)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_minimal()

covariates %>%
  filter(Stock == "Late Stuart") %>%
  ggplot(aes(x = pink, y = lnrs, color = Stock, fill = Stock)) + 
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_minimal()

# assess pink covariate -----------------------------------

top_models %>%
  select(Stock, pink) %>%
  mutate(model = "total") %>%
  bind_rows(top_models_pink %>%
              select(Stock, pink) %>%
              mutate(model = "wild")) %>%
  group_by(Stock) %>%
  pivot_wider(names_from = "model", values_from = "pink") %>%
  mutate(prop = wild / total)


top_models %>%
  select(Stock, pink) %>%
  mutate(model = "total") %>%
  bind_rows(top_models_pink %>%
              select(Stock, pink) %>%
              mutate(model = "wild")) %>%
  group_by(Stock) %>%
  pivot_wider(names_from = "model", values_from = "pink") %>%
  mutate(prop_wild = abs(wild) / abs(total))


# COMBINED scenario: pink_wild substituted + SeaLions/seal frozen at 1978 ----

FREEZE_YEAR  <- 1970
FREEZE_VARS  <- c("SeaLions", "seal")   # these get frozen forward from FREEZE_YEAR
DIRECT_VARS  <- c("pink")               # this gets a direct substitute series

# Warn if any stock is missing a 1978 value for the freeze vars
no_freeze_year <- covariates %>%
  group_by(Stock) %>%
  summarise(has_freeze_year = any(Year == FREEZE_YEAR & if_any(all_of(FREEZE_VARS), ~ !is.na(.))),
            .groups = "drop") %>%
  filter(!has_freeze_year)

if (nrow(no_freeze_year) > 0) {
  warning("No ", FREEZE_YEAR, " SeaLions/seal value for stock(s): ",
          paste(no_freeze_year$Stock, collapse = ", "),
          " -- scenario will equal actual for these stocks.")
}

# Frozen SeaLions/seal columns
frozen_covariates <- covariates %>%
  group_by(Stock) %>%
  mutate(across(all_of(FREEZE_VARS), function(x) {
    freeze_val <- x[Year == FREEZE_YEAR]
    if (length(freeze_val) != 1 || is.na(freeze_val)) x
    else ifelse(Year > FREEZE_YEAR, freeze_val, x)
  })) %>%
  ungroup() %>%
  select(Stock, Year, all_of(FREEZE_VARS)) %>%
  rename_with(~ paste0(., "_scenario"), all_of(FREEZE_VARS))

# Direct substitute pink series (already standardized against pink's own mean/SD)
direct_covariates <- covariates %>%
  select(Stock, Year, pink_scenario = pink_wild)

scenario_covariates <- frozen_covariates %>%
  left_join(direct_covariates, by = c("Stock", "Year"))

obs_with_scenario <- obs %>%
  left_join(scenario_covariates, by = c("Stock", "Year"))

# All base covariate names that have a "_scenario" counterpart -- used by
# compute_cov_term() to decide which columns to swap in scenario mode.
SCENARIO_VARS <- c(FREEZE_VARS, DIRECT_VARS)

compute_cov_term <- function(dat, sel_covs, cov_coefs, use_scenario = FALSE) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  cols <- if (use_scenario) {
    ifelse(sel_covs %in% SCENARIO_VARS, paste0(sel_covs, "_scenario"), sel_covs)
  } else {
    sel_covs
  }
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}

# Look up a stock's top-model coefficients from the dredge table.
get_top_model_terms <- function(stock_name) {
  top_model_row <- top_models %>% filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
    list(
      ra            = top_model_row[["(Intercept)"]],
      rb            = -top_model_row[["spawners"]],
      sel_covs      = sel_covs,
      cov_coefs     = if (length(sel_covs) > 0) as.numeric(top_model_row[sel_covs]) else numeric(0),
      has_top_model = TRUE
    )
  } else {
    list(ra = NA_real_, rb = NA_real_, sel_covs = character(0),
         cov_coefs = numeric(0), has_top_model = FALSE)
  }
}


# RUN MODEL (all stocks
model_all <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()


# pink scenario run: pink salmon frozen at FREEZE_YEAR levels
model_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name = .y$Stock,
                                 retroU     = retroU,
                                 useretro   = useretro,
                                 yrretro    = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

# filter by stock
filtered_model <- if (stock_choice == "All Stocks") {
  model_all
} else {
  model_all %>% filter(Stock == stock_choice)
}

#  stock summary
stock_table <- if (stock_choice == "All Stocks") {
  filtered_model %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch,
      .groups = "drop"
    )
} else {
  filtered_model %>%
    summarise(
      hist_catch  = sum(Catch,  na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      lost_catch  = retro_catch - hist_catch
    )
}
print(stock_table)

# TOTAL SUMMARY
hist_total  <- sum(filtered_model$Catch,  na.rm = TRUE)
retro_total <- sum(filtered_model$retroC, na.rm = TRUE)

cat("Historical Catch: ",
    format(round(hist_total, 0), big.mark = ","), "\n")
cat("Retrospective Catch: ",
    format(round(retro_total, 0), big.mark = ","), "\n")
cat("Difference: ",
    format(round(retro_total - hist_total, 0), big.mark = ","), "\n")

# AGGREGATE TIME SERIES (summed across stocks
summed_ts <- filtered_model %>%
  group_by(Year) %>%
  summarise(
    Catch = sum(Catch, na.rm = TRUE),
    retroC = sum(retroC, na.rm = TRUE),
    RunJacks = sum(RunJacks, na.rm = TRUE),
    retroR = sum(retroR, na.rm = TRUE),
    AdultEscapement = sum(AdultEscapement, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year) %>%
  mutate(
    AdultReturn = lead(RunJacks, LAG_YEARS)
  )

#  RETRO TIME SERIES TABLE
retro_ts_table <- summed_ts %>%
  select(Year, retroR, retroC, Run = RunJacks) %>%
  arrange(Year) %>%
  mutate(
    retroR = round(retroR, 0),
    retroC = round(retroC, 0),
    Run = round(Run, 0)
  )
print(retro_ts_table)

## plot retro catch, etc. ------------------------

# Catch plot
catch_plot <- ggplot(
  summed_ts %>% select(Year, Catch, retroC) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Catch", y = "Catch", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Adult return plot
return_plot <- ggplot(
  summed_ts %>% select(Year, AdultReturns = RunJacks, retroR) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Return", y = "Returns", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Escapement plot
esc_plot <- ggplot(
  summed_ts %>% select(Year, AdultEscapement, retroS) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(title = "Spawners", y = "Spawners", color = "") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

#print(catch_plot)
print(return_plot)
print(esc_plot)

return_plot | esc_plot


# ggsave("catch_plot.png", catch_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("return_plot.png", return_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("esc_plot.png", esc_plot, width = 8, height = 4.5, dpi = 300)

## pink control scenario comparison ---------------------------

filtered_scenario <- if (stock_choice == "All Stocks") {
  model_scenario
} else {
  model_scenario %>% filter(Stock == stock_choice)
}

# Per-stock, per-year comparison
scenario_compare <- filtered_model %>%
  select(Stock, Year, lnRS_actual = retro_lnRS, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    filtered_scenario %>%
      select(Stock, Year, lnRS_scenario = retro_lnRS, R_scenario = retroR, S_scenario = retroS),
    by = c("Stock", "Year")
  ) %>%
  mutate(
    lnRS_diff = lnRS_scenario - lnRS_actual,
    R_diff    = R_scenario - R_actual,
    S_diff    = S_scenario - S_actual
  )
print(scenario_compare)

# Aggregate (summed across stocks) comparison, for plotting
summed_scenario_ts <- filtered_scenario %>%
  group_by(Year) %>%
  summarise(
    retroR = sum(retroR, na.rm = TRUE),
    retroS = sum(retroS, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year)

compare_ts <- summed_ts %>%
  select(Year, R_actual = retroR, S_actual = retroS) %>%
  left_join(
    summed_scenario_ts %>% select(Year, R_scenario = retroR, S_scenario = retroS),
    by = "Year"
  )

# Reconstructed returns: actual vs pinniped-frozen scenario
scenario_return_plot <- ggplot(
  compare_ts %>% select(Year, R_actual, R_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(y = "Returns", color = "", title = "Returns") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

# Reconstructed spawners: actual vs pinniped-frozen scenario
scenario_esc_plot <- ggplot(
  compare_ts %>% select(Year, S_actual, S_scenario) %>% pivot_longer(-Year),
  aes(Year, value, color = name, linetype = name)
) +
  geom_line(linewidth = 1.2) +
  labs(y = "Spawners", color = "", title = "Spawners") +
  guides(linetype = "none") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal()

print(scenario_return_plot)
print(scenario_esc_plot)

scenario_return_plot + theme(legend.position = "none") | scenario_esc_plot 

ggsave("figures/scenario_return_plot_both.png", width = 12, height = 6, dpi = 600)

# ggsave("scenario_return_plot.png", scenario_return_plot, width = 8, height = 4.5, dpi = 300)
# ggsave("scenario_esc_plot.png", scenario_esc_plot, width = 8, height = 4.5, dpi = 300)

## productivity compare -------------------------------------
productivity_compare <- model_all %>%
  select(Stock, Year, lnRS_actual = retro_lnRS) %>%
  left_join(
    model_scenario %>% select(Stock, Year, lnRS_scenario = retro_lnRS),
    by = c("Stock", "Year")
  ) %>%
  pivot_longer(cols = c(lnRS_actual, lnRS_scenario),
               names_to = "scenario", values_to = "lnRS") %>%
  mutate(scenario = recode(scenario,
                           lnRS_actual   = "Actual",
                           lnRS_scenario = "pink salmon controlled"))

ggplot(productivity_compare, aes(Year, lnRS, color = scenario, linetype = scenario)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
  labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
  scale_color_manual(values = c("#4682B4", "#FF4500")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_productivity_full_retro.png", dpi = 600, width = 7, height = 10)

## catch lost by harvest rate scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.05)

# lost_catch = retro_catch - hist_catch (positive = more catch under the
# retro harvest rate than actually occurred; flip the sign if you want
# "catch lost" to read positive when the retro scenario catches less)
catch_lost_by_hr <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    group_by(Stock) %>%
    summarise(
      hist_catch  = sum(Catch, na.rm = TRUE),
      retro_catch = sum(retroC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(lost_catch = retro_catch - hist_catch,
           retroU     = hr)
})

# Per-stock catch lost vs harvest rate
ggplot(catch_lost_by_hr, aes(retroU, lost_catch)) +
  geom_line(linewidth = 1, color = "#4682B4") +
  geom_point(size = 1.5, color = "#4682B4") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Catch lost by retrospective harvest rate, per stock",
       x = "Retrospective harvest rate", y = "Catch lost (retro - historical)") +
  theme_minimal()

# Cumulative catch lost across all stocks
catch_lost_cumulative <- catch_lost_by_hr %>%
  group_by(retroU) %>%
  summarise(total_lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop")

ggplot(catch_lost_cumulative, aes(retroU, total_lost_catch)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  geom_point(size = 2, color = "#FF4500") +
  labs(x = "Retrospective harvest rate", y = "Total catch lost") +
  theme_minimal()

## catch lost by year scenarios  -------------------------------------

harvest_rates <- seq(0, 0.7, by = 0.1)  # tweak step size as you like

catch_lost_ts <- map_dfr(harvest_rates, function(hr) {
  obs %>%
    group_by(Stock) %>%
    group_modify(~ run_retro_model(.x,
                                   stock_name   = .y$Stock,
                                   retroU       = hr,
                                   useretro     = TRUE,
                                   yrretro      = yrretro,
                                   run_scenario = FALSE)) %>%
    ungroup() %>%
    mutate(lost_catch = retroC - Catch,
           retroU     = hr) %>%
    select(Stock, Year, retroU, lost_catch)
})

# Cumulative catch lost over time, by stock
catch_lost_cum_by_stock <- catch_lost_ts %>%
  arrange(Stock, retroU, Year) %>%
  group_by(Stock, retroU) %>%
  mutate(cum_lost_catch = cumsum(replace_na(lost_catch, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(title = "Cumulative catch lost over time, by stock",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

# Cumulative catch lost over time, total across all stocks
catch_lost_cum_total <- catch_lost_ts %>%
  group_by(retroU, Year) %>%
  summarise(lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop") %>%
  arrange(retroU, Year) %>%
  group_by(retroU) %>%
  mutate(cum_lost_catch = cumsum(lost_catch)) %>%
  ungroup()

ggplot(catch_lost_cum_total, aes(Year, cum_lost_catch, color = factor(retroU))) +
  geom_line(linewidth = 1.2) +
  labs(title = "Cumulative catch lost over time, total across stocks",
       x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
  theme_minimal()

## what's the question 

#  Actual harvest rate, actual  -----------------
model_actual_harvest <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,       # unused when useretro = FALSE
                                 useretro     = FALSE,   # use actual historical harvest rate
                                 yrretro      = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()

# Actual harvest rate, pink salmon frozen at FREEZE_YEAR
model_actual_harvest_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,
                                 useretro     = FALSE,
                                 yrretro      = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

#  Catch lost to pinnipeds: scenario catch minus actual catch
catch_lost_to_pinnipeds <- model_actual_harvest %>%
  select(Stock, Year, catch_actual = retroC) %>%
  left_join(
    model_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
    by = c("Stock", "Year")
  ) %>%
  mutate(catch_lost = catch_scenario - catch_actual)  # positive = catch lost to pinnipeds

#  Cumulative catch lost to pinnipeds over time, by stock
catch_lost_cum_by_stock <- catch_lost_to_pinnipeds %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
  ungroup()

ggplot(catch_lost_cum_by_stock, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# filter out stocks that weren't impacted 
p1 <- catch_lost_cum_by_stock %>%
  #filter(Stock %in% c("Seymour", "Bowron", "Late Stuart","Portage", "Weaver", "Stellako")) %>%
  ggplot(aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock, 
             #scales = "free_y"
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()

# Cumulative catch lost to pinks over time, total across stocks
catch_lost_cum_total <- catch_lost_to_pinnipeds %>%
  group_by(Year) %>%
  summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cum_catch_lost = cumsum(catch_lost))

p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

p1 | p2

ggsave("figures/sockeye_lost_catch_full_retro.png", dpi = 600, width = 12, height = 6)


# SST scenarios ---------------------------------
## load data -----------------------------
obs <- read_csv("R/Sockeye Retrospective Shiny App/Walters_model_all-stocks.csv") %>%
  rename(
    Year            = `Year`,
    AdultEscapement = `Adult Escapement`,
    JackEscapement  = `Jack Escapement`,
    TotalEscapement = `Total Escapement`,
    DBE             = `DBE`,
    BelowMissionC   = `Below Mission Catch`,
    AboveMissionC   = `Above Mission Catch`,
    AlaskaCatch     = `Alaska Catch`,
    RunSize         = `Run Size`
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year)) %>%
  arrange(Stock, Year) %>%
  filter(Stock %in% c("Birkenhead",
                      "Bowron",
                      "Chilko",
                      "Cultus",
                      "Early Stuart",
                      "Gates",
                      "Late Shuswap",
                      "Late Stuart",
                      "Pitt",
                      "Portage",
                      "Quesnel",
                      "Raft",
                      "Scotch",
                      "Seymour",
                      "Stellako",
                      "Weaver"))

# specify covariates
COVARIATE_COLS <- c("NPGO", "PDO", "SeaLions", "seal",
                    "adult.sst", "pink", "smolt.sst")

# AIC weighted coefficients for each covariate
top_models <- dredge_models %>%
  group_by(Stock) %>%
  mutate(aic_weight = exp(-0.5 * deltaAIC) / sum(exp(-0.5 * deltaAIC))) %>%
  summarise(
    across(c(`(Intercept)`, spawners, all_of(COVARIATE_COLS)),
           ~ sum(aic_weight * coalesce(.x, 0))),
    .groups = "drop"
  )



# Yearly covariate time series, keyed by Stock + Year
covariates <- read_csv("Data/sockeye_standardized_covariates.csv") %>%
  rename(Year = yr) %>%
  select(Stock, Year, any_of(COVARIATE_COLS))

obs <- obs %>%
  left_join(covariates, by = c("Stock", "Year"))

## SST run -------------------------------------
# Freeze adult.sst and smolt.sst at their 1950-1975
# long-term mean values for each stock.

FREEZE_VARS <- c("adult.sst", "smolt.sst")

sst_means <- covariates %>%
  filter(Year >= 1950, Year <= 1975) %>%
  group_by(Stock) %>%
  summarise(
    adult.sst_mean = mean(adult.sst, na.rm = TRUE),
    smolt.sst_mean = mean(smolt.sst, na.rm = TRUE),
    .groups = "drop"
  )

# Check the stock-specific SST means
print(sst_means)

# Check for stocks without usable SST values during 1950-1975
no_sst_mean <- sst_means %>%
  filter(is.na(adult.sst_mean) | is.na(smolt.sst_mean))

if (nrow(no_sst_mean) > 0) {
  warning(
    "No usable 1950-1975 SST mean for stock(s): ",
    paste(no_sst_mean$Stock, collapse = ", ")
  )
}

# Create scenario columns while retaining all historical covariates
scenario_covariates <- covariates %>%
  left_join(sst_means, by = "Stock") %>%
  mutate(
    adult.sst_scenario = adult.sst_mean,
    smolt.sst_scenario = smolt.sst_mean
  ) %>%
  select(-adult.sst_mean, -smolt.sst_mean)

# Add scenario covariates to obs
obs_with_scenario <- obs %>%
  left_join(
    scenario_covariates %>%
      select(Stock, Year, adult.sst_scenario, smolt.sst_scenario),
    by = c("Stock", "Year")
  )


# Look up a stock's AIC-weighted model coefficients
get_top_model_terms <- function(stock_name) {
  top_model_row <- top_models %>%
    filter(Stock == stock_name)
  
  if (nrow(top_model_row) == 1) {
    sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
    
    list(
      ra            = top_model_row[["(Intercept)"]],
      rb            = -top_model_row[["spawners"]],
      sel_covs      = sel_covs,
      cov_coefs     = if (length(sel_covs) > 0) {
        as.numeric(top_model_row[sel_covs])
      } else {
        numeric(0)
      },
      has_top_model = TRUE
    )
  } else {
    list(
      ra = NA_real_,
      rb = NA_real_,
      sel_covs = character(0),
      cov_coefs = numeric(0),
      has_top_model = FALSE
    )
  }
}

# Linear-predictor contribution from selected covariates.
# In the SST scenario, only adult.sst and smolt.sst
# are replaced by their scenario values.
compute_cov_term <- function(dat, sel_covs, cov_coefs, use_scenario = FALSE) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  
  cols <- if (use_scenario) {
    ifelse(
      sel_covs %in% FREEZE_VARS,
      paste0(sel_covs, "_scenario"),
      sel_covs
    )
  } else {
    sel_covs
  }
  
  as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
}

## retro function ---------------------------------------------

run_retro_model <- function(dat,
                            stock_name,
                            retroU,
                            useretro,
                            yrretro,
                            run_scenario = FALSE) {
  
  dat <- dat %>%
    arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs = pmin(0.95, Catch / RunJacks),
      Ut_obs = if_else(
        stock_name == "Late Shuswap" & Year == 2012,
        pmin(0.7, Catch / RunJacks),
        Ut_obs
      ),
      ENS = pmin(
        1,
        pmax(
          0.0001,
          AdultEscapement / RunJacks / (1 - Ut_obs)
        )
      ),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  terms <- get_top_model_terms(stock_name)
  
  if (terms$has_top_model) {
    ra <- terms$ra
    rb <- terms$rb
    sel_covs <- terms$sel_covs
    cov_coefs <- terms$cov_coefs
  } else {
    # Fallback: no dredge top model for this stock -> plain Ricker fit
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
    ra <- coef(fit)[1]
    rb <- -coef(fit)[2]
    sel_covs <- character(0)
    cov_coefs <- numeric(0)
  }
  
  # Historical covariate contribution.
  # This is always used to calculate wt so that the
  # process-error realization stays identical between scenarios.
  cov_term <- compute_cov_term(
    obs2,
    sel_covs,
    cov_coefs,
    use_scenario = FALSE
  )
  
  # Forward projection:
  # historical run = historical SST
  # SST scenario = 1950-1975 mean SST
  cov_term_proj <- if (run_scenario) {
    compute_cov_term(
      obs2,
      sel_covs,
      cov_coefs,
      use_scenario = TRUE
    )
  } else {
    cov_term
  }
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (ra - rb * AdultEscapement + cov_term)
    )
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  
  retroU_vec <- if (useretro) {
    ifelse(
      obs3$Year >= yrretro,
      retroU,
      obs3$Ut_obs
    )
  } else {
    obs3$Ut_obs
  }
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    
    retro_lnRS[j] <- ra -
      rb * retroS[j] +
      obs3$cov_term_proj[j] +
      obs3$wt[j]
    
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR     = retroR,
      retro_lnRS = retro_lnRS,
      retroU     = retroU_vec,
      retroS     = retroS,
      retroC     = retroC
    )
}

## historical model ---------------------------------------

model_all <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(
    .x,
    stock_name = .y$Stock,
    retroU = retroU,
    useretro = useretro,
    yrretro = yrretro,
    run_scenario = FALSE
  )) %>%
  ungroup()

## SST scenario -------------------------------------------
# adult.sst and smolt.sst are fixed at their
# stock-specific 1950-1975 means.
# All other covariates remain historical.

model_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(
    .x,
    stock_name = .y$Stock,
    retroU = retroU,
    useretro = useretro,
    yrretro = yrretro,
    run_scenario = TRUE
  )) %>%
  ungroup()

##  Actual harvest rate, actual  -----------------
sst_actual_harvest <- obs %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,       # unused when useretro = FALSE
                                 useretro     = FALSE,   # use actual historical harvest rate
                                 yrretro      = yrretro,
                                 run_scenario = FALSE)) %>%
  ungroup()

# Actual harvest rate, pink salmon frozen at FREEZE_YEAR
sst_actual_harvest_scenario <- obs_with_scenario %>%
  group_by(Stock) %>%
  group_modify(~ run_retro_model(.x,
                                 stock_name   = .y$Stock,
                                 retroU       = 0,
                                 useretro     = FALSE,
                                 yrretro      = yrretro,
                                 run_scenario = TRUE)) %>%
  ungroup()

#  Catch lost to sst: scenario catch minus actual catch
catch_lost_to_sst <- sst_actual_harvest %>%
  select(Stock, Year, catch_actual = retroC) %>%
  left_join(
    sst_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
    by = c("Stock", "Year")
  ) %>%
  mutate(catch_lost = catch_scenario - catch_actual)  # positive = catch lost to pinnipeds

#  Cumulative catch lost to pinnipeds over time, by stock
catch_lost_to_sst_cum <- catch_lost_to_sst %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
  ungroup()

ggplot(catch_lost_to_sst_cum, aes(Year, cum_catch_lost)) +
  geom_line(linewidth = 1, color = "#FF4500") +
  facet_wrap(~ Stock) +
  labs(x = "Year", y = "Cumulative catch lost") +
  theme_minimal()


# Cumulative catch lost to pinks over time, total across stocks
catch_lost_cum_total_sst <- catch_lost_to_sst %>%
  group_by(Year) %>%
  summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cum_catch_lost = cumsum(catch_lost))

ggplot(catch_lost_cum_total_sst, aes(Year, cum_catch_lost)) +
  geom_area(fill = "#FF4500", alpha = 0.2) +
  geom_line(linewidth = 1.2, color = "#FF4500") +
  labs(x = "Year", y = "") +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal()

# average yearly since 2000
catch_lost_to_sst <- sst_actual_harvest %>%
  select(Stock, Year, catch_actual = retroC) %>%
  left_join(
    sst_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
    by = c("Stock", "Year")
  ) %>%
  mutate(catch_lost = catch_scenario - catch_actual)  

# stats 
catch_lost_to_sst %>%
  filter(Year > 1999) %>%
  group_by(Stock) %>%
  summarize(average_catch_lost = mean(catch_lost,na.rm=T))

# boxplot 
catch_lost_to_sst %>%
  filter(Year > 1999) %>%
  ggplot(aes(x = Stock, y = log(catch_lost))) +
  geom_boxplot() +
  labs(
    x = "Stock",
    y = "Catch lost to pinnipeds"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

catch_lost_to_sst %>%
  filter(Year > 1999) %>%
  ggplot(aes(x = Stock, y = log(catch_lost))) +
  geom_boxplot() +
  labs(
    x = "Stock",
    y = "Catch lost to pinnipeds"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

catch_lost_to_sst %>%
  filter(Year > 1999) %>%
  ggplot(aes(x=Year,y=catch_lost, color = Stock)) + 
  geom_line() + 
  theme_minimal()

# compare average catch lost -----------------------
catch_lost_to_sst



####
####
# ## OLD PINK SALMON competition -------------------------------------------------
# ggplot(all.sockeye.st, aes(yr,pink)) + 
#   geom_line()
# 
# # how much catch lost if pink salmon abundance in north pacific was kept at 1980s levels 
# 
# FREEZE_YEAR <- 1980
# FREEZE_VARS <- c("pink")
# 
# no_freeze_year <- covariates %>%
#   group_by(Stock) %>%
#   summarise(has_freeze_year = any(Year == FREEZE_YEAR & if_any(all_of(FREEZE_VARS), ~ !is.na(.))),
#             .groups = "drop") %>%
#   filter(!has_freeze_year)
# 
# if (nrow(no_freeze_year) > 0) {
#   warning("No ", FREEZE_YEAR, " pink value for stock(s): ",
#           paste(no_freeze_year$Stock, collapse = ", "),
#           " -- scenario will equal actual for these stocks.")
# }
# 
# scenario_covariates <- covariates %>%
#   group_by(Stock) %>%
#   mutate(across(all_of(FREEZE_VARS), function(x) {
#     #freeze_val <- x[Year == FREEZE_YEAR]
#     freeze_val <- mean(x[Year >= 1960 & Year <= 1970], na.rm = TRUE)
#     if (length(freeze_val) != 1 || is.na(freeze_val)) {
#       x
#     } else {
#       ifelse(Year > FREEZE_YEAR, freeze_val, x)
#     }
#   })) %>%
#   ungroup() %>%
#   select(Stock, Year, all_of(FREEZE_VARS)) %>%
#   rename_with(~ paste0(., "_scenario"), all_of(FREEZE_VARS))
# 
# obs_with_scenario <- obs %>%
#   left_join(scenario_covariates, by = c("Stock", "Year"))
# 
# #  TOP-MODEL / COVARIATE HELPERS
# 
# # Look up a stock's top-model coefficients from the dredge table.
# get_top_model_terms <- function(stock_name) {
#   top_model_row <- top_models %>% filter(Stock == stock_name)
#   
#   if (nrow(top_model_row) == 1) {
#     sel_covs <- COVARIATE_COLS[!is.na(as.numeric(top_model_row[COVARIATE_COLS]))]
#     list(
#       ra            = top_model_row[["(Intercept)"]],
#       rb            = -top_model_row[["spawners"]],
#       sel_covs      = sel_covs,
#       cov_coefs     = if (length(sel_covs) > 0) as.numeric(top_model_row[sel_covs]) else numeric(0),
#       has_top_model = TRUE
#     )
#   } else {
#     list(ra = NA_real_, rb = NA_real_, sel_covs = character(0),
#          cov_coefs = numeric(0), has_top_model = FALSE)
#   }
# }
# 
# # Linear-predictor contribution from the selected covariates.
# # use_scenario = TRUE pulls the "<var>_scenario" (pinniped-frozen) column
# # for any covariate in FREEZE_VARS instead of its actual observed column.
# compute_cov_term <- function(dat, sel_covs, cov_coefs, use_scenario = FALSE) {
#   if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
#   cols <- if (use_scenario) {
#     ifelse(sel_covs %in% FREEZE_VARS, paste0(sel_covs, "_scenario"), sel_covs)
#   } else {
#     sel_covs
#   }
#   as.numeric(as.matrix(dat[cols]) %*% cov_coefs)
# }
# 
# #  RETRO FUNCTION
# run_retro_model <- function(dat,
#                             stock_name,
#                             retroU,
#                             useretro,
#                             yrretro,
#                             run_scenario = FALSE) {
#   
#   dat <- dat %>% arrange(Year)
#   
#   obs2 <- dat %>%
#     mutate(
#       RunJacks = RunSize - JackEscapement,
#       Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
#       Ut_obs  = pmin(0.999, Catch / RunJacks),
#       ENS     = pmin(1, pmax(0.0001,
#                              AdultEscapement / RunJacks / (1 - Ut_obs))),
#       migmort = 1 - ENS
#     ) %>%
#     mutate(
#       AdultReturn = lead(RunJacks, n = LAG_YEARS),
#       lnR_S       = log(AdultReturn / AdultEscapement)
#     )
#   
#   fit_df <- obs2 %>%
#     filter(Year %in% FIT_YEARS) %>%
#     filter(is.finite(lnR_S), is.finite(AdultEscapement))
#   
#   terms <- get_top_model_terms(stock_name)
#   
#   if (terms$has_top_model) {
#     ra <- terms$ra
#     rb <- terms$rb
#     sel_covs <- terms$sel_covs
#     cov_coefs <- terms$cov_coefs
#   } else {
#     # Fallback: no dredge top model for this stock -> plain Ricker fit
#     fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
#     ra <- coef(fit)[1]
#     rb <- -coef(fit)[2]
#     sel_covs <- character(0)
#     cov_coefs <- numeric(0)
#   }
#   
#   # cov_term (actual/historical covariates) always drives wt, the estimated
#   # process-error residual -- this stays fixed across scenarios so that a
#   # scenario comparison isolates the covariate effect rather than mixing it
#   # with a different noise realization.
#   cov_term <- compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = FALSE)
#   
#   # cov_term_proj drives the forward recursive projection below, and is the
#   # only thing that differs between the actual run and the pinniped scenario.
#   cov_term_proj <- if (run_scenario) {
#     compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = TRUE)
#   } else {
#     cov_term
#   }
#   
#   obs3 <- obs2 %>%
#     mutate(
#       cov_term      = cov_term,
#       cov_term_proj = cov_term_proj,
#       wt = lnR_S - (ra - rb * AdultEscapement + cov_term)
#     )
#   
#   n <- nrow(obs3)
#   
#   retroR <- rep(NA_real_, n)
#   retro_lnRS <- rep(NA_real_, n)
#   retroU_vec <- if (useretro)
#     ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
#   else
#     obs3$Ut_obs
#   
#   retroS <- rep(NA_real_, n)
#   retroC <- rep(NA_real_, n)
#   
#   retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
#   
#   retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
#     (1 - retroU_vec[1:LAG_YEARS]) *
#     obs3$ENS[1:LAG_YEARS]
#   
#   retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
#     retroU_vec[1:LAG_YEARS]
#   
#   for (i in (LAG_YEARS + 1):n) {
#     j <- i - LAG_YEARS
#     retro_lnRS[j] <- ra - rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
#     retroR[i] <- retroS[j] * exp(retro_lnRS[j])
#     retroS[i] <- retroR[i] *
#       (1 - retroU_vec[i]) *
#       obs3$ENS[i]
#     retroC[i] <- retroR[i] * retroU_vec[i]
#   }
#   
#   obs3 %>%
#     mutate(
#       retroR     = retroR,
#       retro_lnRS = retro_lnRS,
#       retroU     = retroU_vec,
#       retroS     = retroS,
#       retroC     = retroC
#     )
# }
# 
# # RUN MODEL (all stocks
# model_all <- obs %>%
#   group_by(Stock) %>%
#   group_modify(~ run_retro_model(.x,
#                                  stock_name = .y$Stock,
#                                  retroU     = retroU,
#                                  useretro   = useretro,
#                                  yrretro    = yrretro,
#                                  run_scenario = FALSE)) %>%
#   ungroup()
# 
# 
# # pink scenario run: pink salmon frozen at FREEZE_YEAR levels
# model_scenario <- obs_with_scenario %>%
#   group_by(Stock) %>%
#   group_modify(~ run_retro_model(.x,
#                                  stock_name = .y$Stock,
#                                  retroU     = retroU,
#                                  useretro   = useretro,
#                                  yrretro    = yrretro,
#                                  run_scenario = TRUE)) %>%
#   ungroup()
# 
# # filter by stock
# filtered_model <- if (stock_choice == "All Stocks") {
#   model_all
# } else {
#   model_all %>% filter(Stock == stock_choice)
# }
# 
# #  stock summary
# stock_table <- if (stock_choice == "All Stocks") {
#   filtered_model %>%
#     group_by(Stock) %>%
#     summarise(
#       hist_catch  = sum(Catch,  na.rm = TRUE),
#       retro_catch = sum(retroC, na.rm = TRUE),
#       lost_catch  = retro_catch - hist_catch,
#       .groups = "drop"
#     )
# } else {
#   filtered_model %>%
#     summarise(
#       hist_catch  = sum(Catch,  na.rm = TRUE),
#       retro_catch = sum(retroC, na.rm = TRUE),
#       lost_catch  = retro_catch - hist_catch
#     )
# }
# print(stock_table)
# 
# # TOTAL SUMMARY
# hist_total  <- sum(filtered_model$Catch,  na.rm = TRUE)
# retro_total <- sum(filtered_model$retroC, na.rm = TRUE)
# 
# cat("Historical Catch: ",
#     format(round(hist_total, 0), big.mark = ","), "\n")
# cat("Retrospective Catch: ",
#     format(round(retro_total, 0), big.mark = ","), "\n")
# cat("Difference: ",
#     format(round(retro_total - hist_total, 0), big.mark = ","), "\n")
# 
# # AGGREGATE TIME SERIES (summed across stocks
# summed_ts <- filtered_model %>%
#   group_by(Year) %>%
#   summarise(
#     Catch = sum(Catch, na.rm = TRUE),
#     retroC = sum(retroC, na.rm = TRUE),
#     RunJacks = sum(RunJacks, na.rm = TRUE),
#     retroR = sum(retroR, na.rm = TRUE),
#     AdultEscapement = sum(AdultEscapement, na.rm = TRUE),
#     retroS = sum(retroS, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(Year) %>%
#   mutate(
#     AdultReturn = lead(RunJacks, LAG_YEARS)
#   )
# 
# #  RETRO TIME SERIES TABLE
# retro_ts_table <- summed_ts %>%
#   select(Year, retroR, retroC, Run = RunJacks) %>%
#   arrange(Year) %>%
#   mutate(
#     retroR = round(retroR, 0),
#     retroC = round(retroC, 0),
#     Run = round(Run, 0)
#   )
# print(retro_ts_table)
# 
# ## plot retro catch, etc. ------------------------
# 
# # Catch plot
# catch_plot <- ggplot(
#   summed_ts %>% select(Year, Catch, retroC) %>% pivot_longer(-Year),
#   aes(Year, value, color = name, linetype = name)
# ) +
#   geom_line(linewidth = 1.2) +
#   labs(title = "Catch", y = "Catch", color = "") +
#   guides(linetype = "none") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal()
# 
# # Adult return plot
# return_plot <- ggplot(
#   summed_ts %>% select(Year, AdultReturns = RunJacks, retroR) %>% pivot_longer(-Year),
#   aes(Year, value, color = name, linetype = name)
# ) +
#   geom_line(linewidth = 1.2) +
#   labs(title = "Return", y = "Returns", color = "") +
#   guides(linetype = "none") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal()
# 
# # Escapement plot
# esc_plot <- ggplot(
#   summed_ts %>% select(Year, AdultEscapement, retroS) %>% pivot_longer(-Year),
#   aes(Year, value, color = name, linetype = name)
# ) +
#   geom_line(linewidth = 1.2) +
#   labs(title = "Spawners", y = "Spawners", color = "") +
#   guides(linetype = "none") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal()
# 
# print(catch_plot)
# print(return_plot)
# print(esc_plot)
# 
# # ggsave("catch_plot.png", catch_plot, width = 8, height = 4.5, dpi = 300)
# # ggsave("return_plot.png", return_plot, width = 8, height = 4.5, dpi = 300)
# # ggsave("esc_plot.png", esc_plot, width = 8, height = 4.5, dpi = 300)
# 
# ## pink control scenario comparison ---------------------------
# 
# filtered_scenario <- if (stock_choice == "All Stocks") {
#   model_scenario
# } else {
#   model_scenario %>% filter(Stock == stock_choice)
# }
# 
# # Per-stock, per-year comparison
# scenario_compare <- filtered_model %>%
#   select(Stock, Year, lnRS_actual = retro_lnRS, R_actual = retroR, S_actual = retroS) %>%
#   left_join(
#     filtered_scenario %>%
#       select(Stock, Year, lnRS_scenario = retro_lnRS, R_scenario = retroR, S_scenario = retroS),
#     by = c("Stock", "Year")
#   ) %>%
#   mutate(
#     lnRS_diff = lnRS_scenario - lnRS_actual,
#     R_diff    = R_scenario - R_actual,
#     S_diff    = S_scenario - S_actual
#   )
# print(scenario_compare)
# 
# # Aggregate (summed across stocks) comparison, for plotting
# summed_scenario_ts <- filtered_scenario %>%
#   group_by(Year) %>%
#   summarise(
#     retroR = sum(retroR, na.rm = TRUE),
#     retroS = sum(retroS, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(Year)
# 
# compare_ts <- summed_ts %>%
#   select(Year, R_actual = retroR, S_actual = retroS) %>%
#   left_join(
#     summed_scenario_ts %>% select(Year, R_scenario = retroR, S_scenario = retroS),
#     by = "Year"
#   )
# 
# # Reconstructed returns: actual vs pinniped-frozen scenario
# pink_scenario_return_plot <- ggplot(
#   compare_ts %>% select(Year, R_actual, R_scenario) %>% pivot_longer(-Year),
#   aes(Year, value, color = name, linetype = name)
# ) +
#   geom_line(linewidth = 1.2) +
#   labs(y = "Returns", color = "", title = "pink salmon") +
#   guides(linetype = "none") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal()
# 
# # Reconstructed spawners: actual vs pinniped-frozen scenario
# scenario_esc_plot <- ggplot(
#   compare_ts %>% select(Year, S_actual, S_scenario) %>% pivot_longer(-Year),
#   aes(Year, value, color = name, linetype = name)
# ) +
#   geom_line(linewidth = 1.2) +
#   labs(y = "Spawners", color = "") +
#   guides(linetype = "none") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal()
# 
# print(scenario_return_plot)
# print(scenario_esc_plot)
# 
# # ggsave("scenario_return_plot.png", scenario_return_plot, width = 8, height = 4.5, dpi = 300)
# # ggsave("scenario_esc_plot.png", scenario_esc_plot, width = 8, height = 4.5, dpi = 300)
# 
# ## productivity compare -------------------------------------
# productivity_compare <- model_all %>%
#   select(Stock, Year, lnRS_actual = retro_lnRS) %>%
#   left_join(
#     model_scenario %>% select(Stock, Year, lnRS_scenario = retro_lnRS),
#     by = c("Stock", "Year")
#   ) %>%
#   pivot_longer(cols = c(lnRS_actual, lnRS_scenario),
#                names_to = "scenario", values_to = "lnRS") %>%
#   mutate(scenario = recode(scenario,
#                            lnRS_actual   = "Actual",
#                            lnRS_scenario = "pink salmon controlled"))
# 
# ggplot(productivity_compare, aes(Year, lnRS, color = scenario, linetype = scenario)) +
#   geom_line(linewidth = 1, alpha = 0.6) +
#   facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
#   labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 
# ggsave("figures/sockeye_productivity_pink.png", dpi = 600, width = 7, height = 10)
# 
# productivity_compare %>%
#   filter(Stock %in% c("Birkenhead", "Bowron", "Chilko","Cultus", "Early Stuart", "Raft", "Stellako")) %>%
#   ggplot(aes(Year, lnRS, color = scenario, linetype = scenario)) +
#   geom_line(linewidth = 1, alpha = 0.6) +
#   facet_wrap(~ Stock, scales = "free_y", ncol = 2) +
#   labs(x = "Year", y = "ln(R/S)", color = "", linetype = "") +
#   scale_color_manual(values = c("#4682B4", "#FF4500")) +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 
# ggsave("figures/sockeye_productivity_select-stocks_pink.png", dpi = 600, width = 10, height = 6)
# 
# ## catch lost by harvest rate scenarios  -------------------------------------
# 
# harvest_rates <- seq(0, 0.7, by = 0.05)
# 
# # lost_catch = retro_catch - hist_catch (positive = more catch under the
# # retro harvest rate than actually occurred; flip the sign if you want
# # "catch lost" to read positive when the retro scenario catches less)
# catch_lost_by_hr <- map_dfr(harvest_rates, function(hr) {
#   obs %>%
#     group_by(Stock) %>%
#     group_modify(~ run_retro_model(.x,
#                                    stock_name   = .y$Stock,
#                                    retroU       = hr,
#                                    useretro     = TRUE,
#                                    yrretro      = yrretro,
#                                    run_scenario = FALSE)) %>%
#     ungroup() %>%
#     group_by(Stock) %>%
#     summarise(
#       hist_catch  = sum(Catch, na.rm = TRUE),
#       retro_catch = sum(retroC, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     mutate(lost_catch = retro_catch - hist_catch,
#            retroU     = hr)
# })
# 
# # Per-stock catch lost vs harvest rate
# ggplot(catch_lost_by_hr, aes(retroU, lost_catch)) +
#   geom_line(linewidth = 1, color = "#4682B4") +
#   geom_point(size = 1.5, color = "#4682B4") +
#   facet_wrap(~ Stock, scales = "free_y") +
#   labs(title = "Catch lost by retrospective harvest rate, per stock",
#        x = "Retrospective harvest rate", y = "Catch lost (retro - historical)") +
#   theme_minimal()
# 
# # Cumulative catch lost across all stocks
# catch_lost_cumulative <- catch_lost_by_hr %>%
#   group_by(retroU) %>%
#   summarise(total_lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop")
# 
# ggplot(catch_lost_cumulative, aes(retroU, total_lost_catch)) +
#   geom_line(linewidth = 1.2, color = "#FF4500") +
#   geom_point(size = 2, color = "#FF4500") +
#   labs(title = "Cumulative catch lost across all stocks, by harvest rate",
#        x = "Retrospective harvest rate", y = "Total catch lost") +
#   theme_minimal()
# 
# ## catch lost by year scenarios  -------------------------------------
# 
# harvest_rates <- seq(0, 0.7, by = 0.1)  # tweak step size as you like
# 
# catch_lost_ts <- map_dfr(harvest_rates, function(hr) {
#   obs %>%
#     group_by(Stock) %>%
#     group_modify(~ run_retro_model(.x,
#                                    stock_name   = .y$Stock,
#                                    retroU       = hr,
#                                    useretro     = TRUE,
#                                    yrretro      = yrretro,
#                                    run_scenario = FALSE)) %>%
#     ungroup() %>%
#     mutate(lost_catch = retroC - Catch,
#            retroU     = hr) %>%
#     select(Stock, Year, retroU, lost_catch)
# })
# 
# # Cumulative catch lost over time, by stock
# catch_lost_cum_by_stock <- catch_lost_ts %>%
#   arrange(Stock, retroU, Year) %>%
#   group_by(Stock, retroU) %>%
#   mutate(cum_lost_catch = cumsum(replace_na(lost_catch, 0))) %>%
#   ungroup()
# 
# ggplot(catch_lost_cum_by_stock, aes(Year, cum_lost_catch, color = factor(retroU))) +
#   geom_line(linewidth = 1) +
#   facet_wrap(~ Stock, scales = "free_y") +
#   labs(title = "Cumulative catch lost over time, by stock",
#        x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
#   theme_minimal()
# 
# # Cumulative catch lost over time, total across all stocks
# catch_lost_cum_total <- catch_lost_ts %>%
#   group_by(retroU, Year) %>%
#   summarise(lost_catch = sum(lost_catch, na.rm = TRUE), .groups = "drop") %>%
#   arrange(retroU, Year) %>%
#   group_by(retroU) %>%
#   mutate(cum_lost_catch = cumsum(lost_catch)) %>%
#   ungroup()
# 
# ggplot(catch_lost_cum_total, aes(Year, cum_lost_catch, color = factor(retroU))) +
#   geom_line(linewidth = 1.2) +
#   labs(title = "Cumulative catch lost over time, total across stocks",
#        x = "Year", y = "Cumulative catch lost", color = "Harvest rate") +
#   theme_minimal()
# 
# ## what's the question 
# 
# #  Actual harvest rate, actual pinniped abundance -----------------
# model_actual_harvest <- obs %>%
#   group_by(Stock) %>%
#   group_modify(~ run_retro_model(.x,
#                                  stock_name   = .y$Stock,
#                                  retroU       = 0,       # unused when useretro = FALSE
#                                  useretro     = FALSE,   # use actual historical harvest rate
#                                  yrretro      = yrretro,
#                                  run_scenario = FALSE)) %>%
#   ungroup()
# 
# # Actual harvest rate, pink salmon frozen at FREEZE_YEAR
# model_actual_harvest_scenario <- obs_with_scenario %>%
#   group_by(Stock) %>%
#   group_modify(~ run_retro_model(.x,
#                                  stock_name   = .y$Stock,
#                                  retroU       = 0,
#                                  useretro     = FALSE,
#                                  yrretro      = yrretro,
#                                  run_scenario = TRUE)) %>%
#   ungroup()
# 
# #  Catch lost to pinnipeds: scenario catch minus actual catch
# catch_lost_to_pinnipeds <- model_actual_harvest %>%
#   select(Stock, Year, catch_actual = retroC) %>%
#   left_join(
#     model_actual_harvest_scenario %>% select(Stock, Year, catch_scenario = retroC),
#     by = c("Stock", "Year")
#   ) %>%
#   mutate(catch_lost = catch_scenario - catch_actual)  # positive = catch lost to pinnipeds
# 
# #  Cumulative catch lost to pinnipeds over time, by stock
# catch_lost_cum_by_stock <- catch_lost_to_pinnipeds %>%
#   arrange(Stock, Year) %>%
#   group_by(Stock) %>%
#   mutate(cum_catch_lost = cumsum(replace_na(catch_lost, 0))) %>%
#   ungroup()
# 
# ggplot(catch_lost_cum_by_stock, aes(Year, cum_catch_lost)) +
#   geom_line(linewidth = 1, color = "#FF4500") +
#   facet_wrap(~ Stock, scales = "free_y") +
#   labs(x = "Year", y = "Cumulative catch lost") +
#   theme_minimal()
# 
# # filter out stocks that weren't impacted 
# p1 <- catch_lost_cum_by_stock %>%
#   filter(Stock %in% c("Seymour", "Bowron", "Late Stuart","Portage", "Weaver", "Stellako")) %>%
#   ggplot(aes(Year, cum_catch_lost)) +
#   geom_area(fill = "#FF4500", alpha = 0.2) +
#   geom_line(linewidth = 1, color = "#FF4500") +
#   facet_wrap(~ Stock, 
#              #scales = "free_y"
#   ) +
#   scale_y_continuous(labels = scales::comma) +
#   labs(x = "Year", y = "Cumulative catch lost") +
#   theme_minimal()
# 
# # Cumulative catch lost to pinks over time, total across stocks
# catch_lost_cum_total <- catch_lost_to_pinnipeds %>%
#   group_by(Year) %>%
#   summarise(catch_lost = sum(catch_lost, na.rm = TRUE), .groups = "drop") %>%
#   arrange(Year) %>%
#   mutate(cum_catch_lost = cumsum(catch_lost))
# 
# p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
#   geom_line(linewidth = 1.2, color = "#FF4500") +
#   labs(x = "Year", y = "") +
#   scale_y_continuous(labels = scales::comma) +
#   theme_minimal()
# 
# p2 <- ggplot(catch_lost_cum_total, aes(Year, cum_catch_lost)) +
#   geom_area(fill = "#FF4500", alpha = 0.2) +
#   geom_line(linewidth = 1.2, color = "#FF4500") +
#   labs(x = "Year", y = "") +
#   scale_y_continuous(labels = scales::comma) +
#   theme_minimal()
# 
# p1 | p2
# 
# ggsave("figures/sockeye_lost_catch_pink.png", dpi = 600, width = 10, height = 6)
# 
# 
# # exploration -----------------------
# # distribution of run sizes 
# 
# all.sockeye.st %>%
#   ggplot(aes(x= log(recruits), color = Stock, fill = Stock)) + 
#   geom_density(alpha = 0.5) +
#   theme_minimal()
# 
# # assess pink covariate -----------------------------------
# 
# top_models %>%
#   select(Stock, pink) %>%
#   mutate(model = "total") %>%
#   bind_rows(top_models_pink %>%
#               select(Stock, pink) %>%
#               mutate(model = "wild")) %>%
#   group_by(Stock) %>%
#   pivot_wider(names_from = "model", values_from = "pink") %>%
#   mutate(prop = wild / total)
# 
# 
# top_models %>%
#   select(Stock, pink) %>%
#   mutate(model = "total") %>%
#   bind_rows(top_models_pink %>%
#               select(Stock, pink) %>%
#               mutate(model = "wild")) %>%
#   group_by(Stock) %>%
#   pivot_wider(names_from = "model", values_from = "pink") %>%
#   mutate(prop_wild = abs(wild) / abs(total))
# 
# 
# 


# old retro functions -------------------------
#  RETRO FUNCTION
run_retro_model <- function(dat,
                            stock_name,
                            retroU,
                            useretro,
                            yrretro,
                            run_scenario = FALSE) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs  = pmin(0.999, Catch / RunJacks),
      ENS     = pmin(1, pmax(0.0001,
                             AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S       = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  terms <- get_top_model_terms(stock_name)
  
  if (terms$has_top_model) {
    ra <- terms$ra
    rb <- terms$rb
    sel_covs <- terms$sel_covs
    cov_coefs <- terms$cov_coefs
  } else {
    # Fallback: no dredge top model for this stock -> plain Ricker fit
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
    ra <- coef(fit)[1]
    rb <- -coef(fit)[2]
    sel_covs <- character(0)
    cov_coefs <- numeric(0)
  }
  
  # cov_term (actual/historical covariates) always drives wt, the estimated
  # process-error residual -- this stays fixed across scenarios so that a
  # scenario comparison isolates the covariate effect rather than mixing it
  # with a different noise realization.
  cov_term <- compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = FALSE)
  
  # cov_term_proj drives the forward recursive projection below, and is the
  # only thing that differs between the actual run and the pinniped scenario.
  cov_term_proj <- if (run_scenario) {
    compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = TRUE)
  } else {
    cov_term
  }
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (ra - rb * AdultEscapement + cov_term)
    )
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroU_vec <- if (useretro)
    ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
  else
    obs3$Ut_obs
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- ra - rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR     = retroR,
      retro_lnRS = retro_lnRS,
      retroU     = retroU_vec,
      retroS     = retroS,
      retroC     = retroC
    )
}

#  RETRO FUNCTION
run_retro_model <- function(dat,
                            stock_name,
                            retroU,
                            useretro,
                            yrretro,
                            run_scenario = FALSE) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs  = pmin(0.999, Catch / RunJacks),
      ENS     = pmin(1, pmax(0.0001,
                             AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S       = log(AdultReturn / AdultEscapement)
    )
  
  fit_df <- obs2 %>%
    filter(Year %in% FIT_YEARS) %>%
    filter(is.finite(lnR_S), is.finite(AdultEscapement))
  
  terms <- get_top_model_terms(stock_name)
  
  if (terms$has_top_model) {
    ra <- terms$ra
    rb <- terms$rb
    sel_covs <- terms$sel_covs
    cov_coefs <- terms$cov_coefs
  } else {
    # Fallback: no dredge top model for this stock -> plain Ricker fit
    fit <- lm(lnR_S ~ AdultEscapement, data = fit_df)
    ra <- coef(fit)[1]
    rb <- -coef(fit)[2]
    sel_covs <- character(0)
    cov_coefs <- numeric(0)
  }
  
  # cov_term (actual/historical covariates) always drives wt, the estimated
  # process-error residual -- this stays fixed across scenarios so that a
  # scenario comparison isolates the covariate effect rather than mixing it
  # with a different noise realization.
  cov_term <- compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = FALSE)
  
  # cov_term_proj drives the forward recursive projection below, and is the
  # only thing that differs between the actual run and the pinniped scenario.
  cov_term_proj <- if (run_scenario) {
    compute_cov_term(obs2, sel_covs, cov_coefs, use_scenario = TRUE)
  } else {
    cov_term
  }
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (ra - rb * AdultEscapement + cov_term)
    )
  
  n <- nrow(obs3)
  
  retroR <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroU_vec <- if (useretro)
    ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs)
  else
    obs3$Ut_obs
  
  retroS <- rep(NA_real_, n)
  retroC <- rep(NA_real_, n)
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    (1 - retroU_vec[1:LAG_YEARS]) *
    obs3$ENS[1:LAG_YEARS]
  
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] *
    retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- ra - rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] *
      (1 - retroU_vec[i]) *
      obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>%
    mutate(
      retroR     = retroR,
      retro_lnRS = retro_lnRS,
      retroU     = retroU_vec,
      retroS     = retroS,
      retroC     = retroC
    )
}