# ============================================================
# Sockeye exact-spec low-point retrospective ratio (X / Xmin) --
# ALL THREE SCENARIOS (Pinniped, SST, Pink), consolidating and fixing
# sockeye_pinniped_lowpoint_ratio_exact.R:
#
# FIX 1 -- real bug found in the pinniped-only version: the synthetic
# years' frozen scenario value was pulled from stock_data[[..._scenario]][1]
# (that stock's FIRST real row, typically 1940s-50s -- BEFORE the 1970
# freeze year even applies, so it's actually that stock's real early
# covariate value, not the intended frozen-1970 constant). If that early
# value happened to be NA for a stock, it injected NA into every
# synthetic year, cascading forward -- the likely cause of Birkenhead/
# Bowron/Cultus/Early Stuart/Late Stuart returning NA. Fixed here by
# pulling the scenario value from the LAST real year instead (guaranteed
# well past 1970, so guaranteed to already be the frozen constant).
#
# FIX 2 -- reinstates the runaway-draw guard (uncapped exponential
# recursion, see chat history) that existed in the earlier all-scenarios
# script but was never carried into the exact-spec version. The
# bootstrap-sampled synthetic residual makes this MORE likely to trigger
# now, not less -- an unlucky large historical residual can get sampled
# into a synthetic year and compound forward with nothing to check it.
#
# Otherwise identical design to the pinniped-only exact-spec script:
# fixed 10-year low window, 12-year minimum freeze with synthetic-year
# extrapolation (non-scenario covariates at their final-10-real-year
# mean, wt bootstrap-resampled from that draw's own real-year residuals),
# generalized here to loop over all three scenario drivers.
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50

LOW_PERIOD_YEARS      <- 10
POST_LOW_PERIOD_YEARS <- 12
X_YEARS                <- 10
FINAL_MEAN_YEARS        <- 10
RUNAWAY_MULTIPLE        <- 10   # a real population can't plausibly return at >10x its historical max

STOCKS_MC <- setdiff(STOCKS, "Late Shuswap")

STOCK_COVARIATES <- list(
  Birkenhead     = c("PDO", "SeaLions"),
  Bowron         = c("PDO", "SeaLions", "pink", "NPGO"),
  Chilko         = c("PDO", "SeaLions"),
  Cultus         = c("SeaLions", "smolt.sst"),
  `Early Stuart` = c("seal", "adult.sst"),
  Gates          = c("smolt.sst"),
  `Late Stuart`  = c("pink"),
  Pitt           = c("SeaLions", "smolt.sst"),
  Portage        = c("pink"),
  Quesnel        = c("adult.sst"),
  Raft           = c("SeaLions", "smolt.sst"),
  Scotch         = c("SeaLions", "pink", "NPGO"),
  Seymour        = c("pink", "smolt.sst"),
  Stellako       = c("SeaLions", "pink"),
  Weaver         = c("pink")
)

SCENARIOS_MC <- list(
  "Pinniped scenario" = c("SeaLions", "seal"),
  "SST scenario"       = c("adult.sst", "smolt.sst"),
  "Pink scenario"      = "pink"
)

# ------------------------------------------------------------
# 1. FIXED 10-YEAR LOW-ABUNDANCE WINDOW
# ------------------------------------------------------------

observed_returns <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  select(Stock, Year, RunJacks)

available_stocks <- intersect(STOCKS_MC, unique(observed_returns$Stock))
missing_stocks <- setdiff(STOCKS_MC, available_stocks)
if (length(missing_stocks) > 0) message("Not found in `obs`, excluding: ", paste(missing_stocks, collapse = ", "))
STOCKS_MC <- available_stocks

low_periods <- observed_returns %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  filter(n() >= LOW_PERIOD_YEARS) %>%
  mutate(roll_mean = zoo::rollapply(RunJacks, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA)) %>%
  filter(!is.na(roll_mean)) %>%
  slice_min(roll_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, low_period_end = Year, low_period_start = Year - LOW_PERIOD_YEARS + 1)

print(low_periods)

historical_max <- observed_returns %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 2. FREEZE WINDOW + SYNTHETIC-YEAR EXTENSION -- scenario-general,
#    built once per stock PER SCENARIO (the set of "non-scenario"
#    covariates being projected differs by scenario)
# ------------------------------------------------------------

build_stock_window <- function(stock_name, scenario_vars) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  non_scenario_covariates <- setdiff(extra_covariates, scenario_vars)
  
  stock_data <- obs %>% filter(Stock == stock_name) %>% arrange(Year)
  low_p <- low_periods %>% filter(Stock == stock_name)
  if (nrow(low_p) == 0) return(NULL)
  
  freeze_start <- low_p$low_period_end + 1
  natural_data_end <- max(stock_data$Year, na.rm = TRUE)
  freeze_end <- max(freeze_start + POST_LOW_PERIOD_YEARS - 1, natural_data_end)
  
  extended_data <- stock_data
  if (freeze_end > natural_data_end) {
    final_mean_start <- natural_data_end - FINAL_MEAN_YEARS + 1
    final_means <- stock_data %>%
      filter(Year >= final_mean_start, Year <= natural_data_end) %>%
      summarise(across(all_of(non_scenario_covariates), ~ mean(.x, na.rm = TRUE)))
    
    synthetic_years <- tibble(Year = (natural_data_end + 1):freeze_end, Stock = stock_name)
    for (cv in non_scenario_covariates) synthetic_years[[cv]] <- final_means[[cv]]
    
    # FIX: pull the frozen scenario value from the LAST real row (Year ==
    # natural_data_end), not the first -- guaranteed to already be the
    # frozen constant (well past FREEZE_YEAR/SST_BASELINE_YEARS), unlike
    # row 1 which could be an early, possibly-NA, unfrozen real value
    last_row <- stock_data %>% filter(Year == natural_data_end)
    for (cv in intersect(extra_covariates, scenario_vars)) {
      synthetic_years[[cv]] <- NA_real_
      synthetic_years[[paste0(cv, "_scenario")]] <- last_row[[paste0(cv, "_scenario")]][1]
    }
    extended_data <- bind_rows(stock_data, synthetic_years)
  }
  
  list(data = extended_data, low_p = low_p, freeze_start = freeze_start, freeze_end = freeze_end,
       natural_data_end = natural_data_end)
}

# ------------------------------------------------------------
# 3. compute_cov_term_from_year() + run_retro_model_lowpoint_draw()
#    (with the wt bootstrap-resampling fallback for synthetic years)
# ------------------------------------------------------------

compute_cov_term_from_year <- function(dat, sel_covs, cov_coefs, scenario_vars, scenario_start_year) {
  if (length(sel_covs) == 0) return(rep(0, nrow(dat)))
  mat <- matrix(NA_real_, nrow = nrow(dat), ncol = length(sel_covs))
  for (k in seq_along(sel_covs)) {
    cv <- sel_covs[k]
    if (cv %in% scenario_vars) {
      mat[, k] <- if_else(dat$Year >= scenario_start_year, dat[[paste0(cv, "_scenario")]], dat[[cv]])
    } else {
      mat[, k] <- dat[[cv]]
    }
  }
  as.numeric(mat %*% cov_coefs)
}

run_retro_model_lowpoint_draw <- function(dat, stock_name, terms, retroU, useretro, yrretro,
                                          scenario_vars, scenario_start_year) {
  
  dat <- dat %>% arrange(Year)
  
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs = pmin(0.95, Catch / RunJacks),
      ENS = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(AdultReturn = lead(RunJacks, n = LAG_YEARS), lnR_S = log(AdultReturn / AdultEscapement))
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term_from_year(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars, scenario_start_year)
  
  obs3 <- obs2 %>%
    mutate(cov_term = cov_term, cov_term_proj = cov_term_proj,
           wt_raw = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term))
  
  real_wt <- obs3$wt_raw[!is.na(obs3$wt_raw)]
  na_idx  <- which(is.na(obs3$wt_raw))
  wt_filled <- obs3$wt_raw
  if (length(na_idx) > 0) {
    wt_filled[na_idx] <- if (length(real_wt) == 0) 0 else sample(real_wt, length(na_idx), replace = TRUE)
  }
  obs3$wt <- wt_filled
  
  n <- nrow(obs3)
  retroR <- rep(NA_real_, n); retro_lnRS <- rep(NA_real_, n)
  retroS <- rep(NA_real_, n); retroC <- rep(NA_real_, n)
  retroU_vec <- if (useretro) ifelse(obs3$Year >= yrretro, retroU, obs3$Ut_obs) else obs3$Ut_obs
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * (1 - retroU_vec[1:LAG_YEARS]) * obs3$ENS[1:LAG_YEARS]
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- terms$ra - terms$rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- retroS[j] * exp(retro_lnRS[j])
    retroS[i] <- retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i]
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>% mutate(retroR = retroR, retro_lnRS = retro_lnRS, retroU = retroU_vec, retroS = retroS, retroC = retroC)
}

# ------------------------------------------------------------
# 4. LOAD POSTERIOR DRAWS
# ------------------------------------------------------------

build_stock_draws <- function(stock_name, n_draws) {
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  samples <- if (exists("sockeye_fits")) sockeye_fits[[stock_name]]$samples else
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds"))
  
  bnames <- grep("^b_", varnames(samples), value = TRUE)
  post <- as.data.frame(as.matrix(samples[, c("intercept", bnames)]))
  idx <- sample(seq_len(nrow(post)), n_draws, replace = nrow(post) < n_draws)
  post <- post[idx, ]
  
  normalize_name <- function(x) tolower(gsub("[._]", "", x))
  bnames_stripped <- sub("^b_", "", bnames)
  
  draws <- tibble(draw = seq_len(n_draws), ra = post$intercept, rb_raw = post$b_spawners)
  for (cv in extra_covariates) {
    match_idx <- which(normalize_name(bnames_stripped) == normalize_name(cv))
    if (length(match_idx) != 1) stop("Could not match covariate '", cv, "' for stock ", stock_name)
    draws[[cv]] <- post[[bnames[match_idx]]]
  }
  draws
}

stock_draws <- setNames(lapply(STOCKS_MC, build_stock_draws, n_draws = N_DRAWS), STOCKS_MC)

# ------------------------------------------------------------
# 5. MONTE CARLO LOOP -- 3 scenarios x N_DRAWS x stocks
# ------------------------------------------------------------

lowpoint_results <- vector("list", length(STOCKS_MC) * N_DRAWS * length(SCENARIOS_MC))
counter <- 0

for (stock_name in STOCKS_MC) {
  
  draws_table      <- stock_draws[[stock_name]]
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
  
  for (scenario_label in names(SCENARIOS_MC)) {
    
    scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
    w <- build_stock_window(stock_name, scenario_vars)
    if (is.null(w)) next
    
    X_start <- w$freeze_end - X_YEARS + 1
    X_end   <- w$freeze_end
    
    for (d in seq_len(N_DRAWS)) {
      
      draw_row <- draws_table[d, ]
      terms <- list(
        ra = draw_row$ra, rb = -draw_row$rb_raw, sel_covs = extra_covariates,
        cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
      )
      
      result <- run_retro_model_lowpoint_draw(
        w$data, stock_name, terms, retroU_default, useretro_default, yrretro_default,
        scenario_vars = scenario_vars, scenario_start_year = w$freeze_start
      )
      
      Xmin <- result %>%
        filter(Year >= w$low_p$low_period_start, Year <= w$low_p$low_period_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
      
      X <- result %>%
        filter(Year >= X_start, Year <= X_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = stock_name, scenario = scenario_label, draw = d, Xmin = Xmin, X = X, ratio = X / Xmin,
        runaway = !is.finite(X) | (X / hist_max_this_stock) > RUNAWAY_MULTIPLE
      )
    }
  }
  
  message("done: ", stock_name)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "sockeye_lowpoint_ratio_exact_all_scenarios.rds")

# ------------------------------------------------------------
# 6. DIAGNOSTICS
# ------------------------------------------------------------

runaway_fraction <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_runaway = sum(runaway, na.rm = TRUE),
            pct_runaway = 100 * n_runaway / n_draws, .groups = "drop") %>%
  arrange(desc(pct_runaway))
print(runaway_fraction, n = Inf)

na_check <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_na = sum(is.na(ratio)), .groups = "drop") %>%
  filter(n_na > 0)
if (nrow(na_check) > 0) {
  cat("\nStill seeing NA ratios (should be much rarer after the frozen-value fix) for:\n")
  print(na_check, n = Inf)
}

zero_variance_check <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(sd_ratio = sd(ratio, na.rm = TRUE), .groups = "drop") %>%
  arrange(sd_ratio)
print(zero_variance_check, n = Inf)

# ------------------------------------------------------------
# 7. SUMMARIZE (excluding runaway draws) + PAIRED, LOG-SCALE PLOT
# ------------------------------------------------------------

lowpoint_stats <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.05, na.rm = TRUE),
    q95    = quantile(ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
print(lowpoint_stats, n = Inf)

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

ggplot(lowpoint_stats, aes(x = Stock, y = median, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  scale_y_log10(labels = scales::comma) +
  labs(x = NULL, y = "X / Xmin (median, 90% CI, log scale)", fill = "Scenario driver") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")


ggsave("figures/sockeye_lowpoint_ratio_exact_paired_logscale.png", width = 13, height = 6.5, dpi = 600)


# ============================================================
# Diagnostic pass for sockeye_lowpoint_ratio_exact_all_scenarios.R
#
# Run after lowpoint_summary, lowpoint_stats, runaway_fraction,
# stock_draws, historical_max all exist in your session.
# ============================================================
# ------------------------------------------------------------
# 1. RUNAWAY PROPORTION -- full detail, not just the summary table.
#    Worth knowing not just THAT some stocks have runaway draws, but
#    whether it's a handful of outliers or a large chunk of the posterior
#    (which would mean the "clean" draws left over aren't a representative
#    sample anymore).
# ------------------------------------------------------------

cat("=== Runaway proportion, sorted worst first ===\n")
print(runaway_fraction %>% arrange(desc(pct_runaway)), n = Inf)

cat("\n=== Stock/scenario combos losing >20% of draws to the runaway filter ===\n")
print(runaway_fraction %>% filter(pct_runaway > 20), n = Inf)

# ------------------------------------------------------------
# 2. HOW MANY YEARS WERE ACTUALLY SYNTHETIC, PER STOCK/SCENARIO --
#    reconstructs what build_stock_window() produced, since the main
#    script doesn't save this directly. Runaway rate and synthetic-year
#    count should correlate if the bootstrap residual is the driver
#    (more synthetic years = more chances for a large resampled
#    residual to compound).
# ------------------------------------------------------------

years_used_diagnostic <- bind_rows(lapply(STOCKS_MC, function(stock_name) {
  bind_rows(lapply(names(SCENARIOS_MC), function(scenario_label) {
    scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], STOCK_COVARIATES[[stock_name]])
    w <- build_stock_window(stock_name, scenario_vars)
    if (is.null(w)) return(NULL)
    
    X_start <- w$freeze_end - X_YEARS + 1
    X_end   <- w$freeze_end
    X_total <- X_end - X_start + 1
    X_real  <- max(0, min(X_end, w$natural_data_end) - X_start + 1)
    
    tibble(
      Stock = stock_name, scenario = scenario_label,
      freeze_total_years = w$freeze_end - w$freeze_start + 1,
      freeze_synthetic_years = max(0, w$freeze_end - w$natural_data_end),
      X_window_years = X_total,
      X_window_synthetic_years = X_total - X_real
    )
  }))
}))
print(years_used_diagnostic, n = Inf)

cat("\n=== Correlation check: does synthetic-year count track with runaway rate? ===\n")
correlation_check <- years_used_diagnostic %>%
  left_join(runaway_fraction, by = c("Stock", "scenario")) %>%
  arrange(desc(pct_runaway))
print(correlation_check %>% select(Stock, scenario, X_window_synthetic_years, pct_runaway), n = Inf)

# ------------------------------------------------------------
# 3. HOW MANY USABLE DRAWS REMAIN AFTER FILTERING -- a stock/scenario
#    losing most of its draws to the runaway filter means the reported
#    median/CI are based on a much smaller (and non-random) sample than
#    N_DRAWS suggests
# ------------------------------------------------------------

cat("\n=== Draws remaining after runaway exclusion ===\n")
usable_draws <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_total = n(), n_usable = sum(!runaway), pct_usable = 100 * n_usable / n_total,
            .groups = "drop") %>%
  arrange(pct_usable)
print(usable_draws, n = Inf)

# ------------------------------------------------------------
# 4. SANITY RANGE CHECK ON Xmin/X -- Xmin should always be positive and
#    "reasonable" (order of magnitude similar to real historical
#    RunJacks for that stock); X (post-runaway-filter) should be too.
#    Flags anything that looks structurally off, not just numerically
#    extreme.
# ------------------------------------------------------------

cat("\n=== Xmin range per stock (should be positive, roughly in line with real RunJacks) ===\n")
xmin_range <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(Xmin_min = min(Xmin, na.rm = TRUE), Xmin_max = max(Xmin, na.rm = TRUE),
            Xmin_sd = sd(Xmin, na.rm = TRUE), .groups = "drop")
print(xmin_range, n = Inf)

if (any(xmin_range$Xmin_min <= 0, na.rm = TRUE)) {
  cat("\nWARNING: some Xmin values are <= 0 -- shouldn't happen for a real recruitment/return quantity:\n")
  print(xmin_range %>% filter(Xmin_min <= 0), n = Inf)
}

cat("\n=== X range (non-runaway draws only) vs. historical max, per stock/scenario ===\n")
x_range_check <- lowpoint_summary %>%
  filter(!runaway) %>%
  left_join(historical_max, by = "Stock") %>%
  group_by(Stock, scenario) %>%
  summarise(
    X_min = min(X, na.rm = TRUE), X_max = max(X, na.rm = TRUE),
    historical_max_return = first(historical_max_return),
    max_ratio_to_historical = X_max / first(historical_max_return),
    .groups = "drop"
  ) %>%
  arrange(desc(max_ratio_to_historical))
print(x_range_check, n = Inf)

# ------------------------------------------------------------
# 5. DID EVERY STOCK/SCENARIO ACTUALLY PRODUCE N_DRAWS ROWS? -- catches
#    silent early termination or a stock/scenario combo that got skipped
#    entirely (e.g. via the `if (is.null(w)) next` guard)
# ------------------------------------------------------------

cat("\n=== Row count check: every stock x scenario should have exactly N_DRAWS rows ===\n")
row_count_check <- lowpoint_summary %>%
  count(Stock, scenario, name = "n_rows")
print(row_count_check, n = Inf)

expected_combos <- expand.grid(Stock = STOCKS_MC, scenario = names(SCENARIOS_MC), stringsAsFactors = FALSE)
missing_combos <- anti_join(expected_combos, row_count_check, by = c("Stock", "scenario"))
if (nrow(missing_combos) > 0) {
  cat("\nMISSING ENTIRELY (build_stock_window returned NULL, or some other skip) --\n")
  print(missing_combos, n = Inf)
}

short_combos <- row_count_check %>% filter(n_rows != N_DRAWS)
if (nrow(short_combos) > 0) {
  cat("\nWrong row count (expected", N_DRAWS, "each):\n")
  print(short_combos, n = Inf)
}


