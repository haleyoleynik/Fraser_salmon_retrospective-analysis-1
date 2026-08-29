# ============================================================
# Sockeye X/Xmin ratio, all three scenarios, WITH SYNTHETIC-YEAR
# EXTENSION -- clean rewrite with ONE unified rule for missing
# covariate data:
#
#   Whenever a covariate (or its frozen "_scenario" counterpart) runs
#   out of real values -- whether that's because we've entered
#   synthetic projection years, OR because the real data itself has a
#   trailing gap before the freeze even starts (e.g. a stock's SeaLions
#   monitoring stopped a few years before present) -- everything from
#   that point forward is filled with the MEAN OF THE 10 YEARS
#   IMMEDIATELY BEFORE THE GAP STARTS. One rule, applied uniformly.
#
# THE THREE-STEP PROCESS, exactly as specified:
#   1. Find each stock's own 10-year low period: real observed Return,
#      lowest 10-year rolling mean.
#   2. Project forward: run the model through the low period using real
#      covariates. From the year after, freeze the scenario's covariates
#      at their historical reference value (1970 for pinniped, the
#      1950-1975 mean for SST, wild-pink substitution for pink -- all
#      already built as "_scenario" columns in your cleaned script) for
#      at least 3 generations (12 years for GT=4 stocks, 15 for Pitt),
#      longer if real data extends further; use synthetic years if real
#      data runs out first.
#   3. Calculate the increase: X (final 10 years of the frozen period) /
#      Xmin (the low period itself).
#
# DEPENDENCIES: obs, model_actual_retro/model_pinniped_retro/
# model_sst_retro/model_pink_retro, compute_cov_term, STOCKS, LAG_YEARS,
# retroU_default, useretro_default, yrretro_default (from your cleaned
# sockeye script), sockeye_fits (or saved RDS) from
# sockeye_jags_models_bounded_b.R.
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 5000  
FINAL_PERIOD_YEARS <- 10
LOW_PERIOD_YEARS <- 10
RUNAWAY_MULTIPLE <- 10
GAP_FILL_WINDOW <- 10   # width of the "mean of the previous N years" used to fill any gap
MIN_YEARS_BEFORE_LOW_PERIOD <- 10   # a candidate low window is only considered if it's
# preceded by a full 10 years of real data -- prevents
# picking a "low period" right at the very start of a
# stock's record, before a rolling mean even has context

STOCKS_MC <- setdiff(STOCKS, "Late Shuswap")

STOCK_COVARIATES <- list(
  Birkenhead     = c("PDO", "SeaLions"),
  Bowron         = c("SeaLions", "pink"),
  Chilko         = c("PDO", "SeaLions"),
  Cultus         = c("SeaLions", "smolt.sst"),
  `Early Stuart` = c("seal", "adult.sst"),
  Gates          = c("smolt.sst"),
  `Late Stuart`  = c("pink"),
  Pitt           = c("SeaLions", "smolt.sst"),
  Portage        = c("pink"),
  Quesnel        = c("adult.sst"),
  Raft           = c("SeaLions", "smolt.sst"),
  Scotch         = c("pink", "NPGO"),
  Seymour        = c("pink", "smolt.sst"),
  Stellako       = c("SeaLions", "pink"),
  Weaver         = c("pink")
)

SCENARIOS_MC <- list(
  "Pinniped scenario" = c("SeaLions", "seal"),
  "SST scenario"       = c("adult.sst", "smolt.sst"),
  "Pink scenario"      = "pink"
)

GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
#GENERATION_TIME["Pitt"] <- 5

SMSY_PERIOD_START_YEAR <- 1950
SMSY_PERIOD_END_YEAR <- 1970

# ------------------------------------------------------------
# 1. STOCK AVAILABILITY CHECK
# ------------------------------------------------------------

available_stocks <- intersect(STOCKS_MC, unique(obs$Stock))
missing_stocks <- setdiff(STOCKS_MC, available_stocks)
if (length(missing_stocks) > 0) message("Not found in `obs`, excluding: ", paste(missing_stocks, collapse = ", "))
STOCKS_MC <- available_stocks

# ------------------------------------------------------------
# 2. STEP 1 -- EACH STOCK'S OWN 10-YEAR LOW PERIOD
#    (real observed Return, lowest 10-year rolling mean)
# ------------------------------------------------------------

observed_returns <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  select(Stock, Year, RunJacks)

low_periods <- observed_returns %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(
    roll_mean = zoo::rollapply(RunJacks, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA),
    roll_start = Year - LOW_PERIOD_YEARS + 1,
    min_year = min(Year)
  ) %>%
  filter(!is.na(roll_mean), roll_start >= min_year + MIN_YEARS_BEFORE_LOW_PERIOD) %>%
  slice_min(roll_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, low_period_start = roll_start, low_period_end = Year, low_period_mean = roll_mean)

print(low_periods %>% arrange(Stock))

historical_max <- observed_returns %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 3. UNIFIED GAP-FILL: mean of the GAP_FILL_WINDOW years immediately
#    before a TRAILING gap starts (whether that gap begins within real
#    calendar years or at the synthetic boundary). Applied uniformly to
#    every covariate column actually used, after real+synthetic rows
#    are already combined.
# ------------------------------------------------------------

fill_trailing_gap <- function(vals) {
  n <- length(vals)
  if (all(is.na(vals))) return(vals)   # nothing to anchor a fill on
  last_valid <- max(which(!is.na(vals)))
  if (last_valid == n) return(vals)    # no trailing gap
  window_start <- max(1, last_valid - GAP_FILL_WINDOW + 1)
  fill_val <- mean(vals[window_start:last_valid], na.rm = TRUE)
  vals[(last_valid + 1):n] <- fill_val
  vals
}

# ------------------------------------------------------------
# 4. STEP 2 -- FREEZE WINDOW + SYNTHETIC-YEAR EXTENSION, per stock per
#    scenario. Synthetic rows start as Year/Stock only (everything else
#    NA); the unified gap-fill (step 3) handles populating every
#    covariate afterward, for BOTH the synthetic rows and any trailing
#    real-year gaps.
# ------------------------------------------------------------

build_stock_window <- function(stock_name, scenario_vars) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  scenario_covariates <- intersect(extra_covariates, scenario_vars)
  non_scenario_covariates <- setdiff(extra_covariates, scenario_vars)
  
  stock_data <- obs %>% filter(Stock == stock_name) %>% arrange(Year)
  low_p <- low_periods %>% filter(Stock == stock_name)
  if (nrow(low_p) == 0) return(NULL)
  
  gt <- GENERATION_TIME[[stock_name]]
  post_low_period_years <- 3 * gt   # 12 for GT=4, 15 for Pitt (GT=5)
  
  freeze_start <- low_p$low_period_end + 1
  natural_data_end <- max(stock_data$Year, na.rm = TRUE)
  freeze_end <- max(freeze_start + post_low_period_years - 1, natural_data_end)
  
  extended_data <- stock_data
  if (freeze_end > natural_data_end) {
    synthetic_years <- tibble(Year = (natural_data_end + 1):freeze_end, Stock = stock_name)
    for (cv in extra_covariates) synthetic_years[[cv]] <- NA_real_
    for (cv in scenario_covariates) synthetic_years[[paste0(cv, "_scenario")]] <- NA_real_
    extended_data <- bind_rows(stock_data, synthetic_years)
  }
  
  extended_data <- extended_data %>% arrange(Year)
  
  # Non-scenario covariates in synthetic years: held at the mean of the
  # most recent 10 REAL, valid years (whatever they are chronologically)
  # -- NOT specifically the low period's own mean. Reverted back to this
  # from the low-period-specific version: some stocks' low periods
  # themselves overlap a real data gap (e.g. Cultus's SeaLions/smolt.sst
  # have no valid data at all within its own low period), which made the
  # low-period mean come out NaN. This "most recent valid years" version
  # reaches back further when needed and always finds SOME usable value.
  cols_to_fill <- unique(c(extra_covariates, paste0(scenario_covariates, "_scenario")))
  cols_to_fill <- intersect(cols_to_fill, names(extended_data))
  for (col in cols_to_fill) {
    extended_data[[col]] <- fill_trailing_gap(extended_data[[col]])
  }
  
  list(data = extended_data, low_p = low_p, freeze_start = freeze_start, freeze_end = freeze_end,
       natural_data_end = natural_data_end, post_low_period_years = post_low_period_years)
}

# ------------------------------------------------------------
# 5. compute_cov_term_from_year() + run_retro_model_lowpoint_draw()
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
      Ut_obs = if_else(stock_name == "Late Shuswap" & Year == 2012, pmin(0.7, Catch / RunJacks), Ut_obs),
      ENS = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(AdultReturn = lead(RunJacks, n = LAG_YEARS), lnR_S = log(AdultReturn / AdultEscapement))
  
  # ENS has no real-data equivalent for synthetic years (RunSize/
  # AdultEscapement are never populated there) -- same unified gap-fill
  # rule applied here too, for consistency with the covariates.
  obs2$ENS <- fill_trailing_gap(obs2$ENS)
  
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
  
  NUMERIC_CEILING <- 1e15   # prevents floating-point overflow to Inf during long synthetic runs
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * (1 - retroU_vec[1:LAG_YEARS]) * obs3$ENS[1:LAG_YEARS]
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * retroU_vec[1:LAG_YEARS]
  
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- terms$ra - terms$rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- pmin(retroS[j] * exp(retro_lnRS[j]), NUMERIC_CEILING)
    retroS[i] <- pmin(retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i], NUMERIC_CEILING)
    retroC[i] <- retroR[i] * retroU_vec[i]
  }
  
  obs3 %>% mutate(retroR = retroR, retro_lnRS = retro_lnRS, retroU = retroU_vec, retroS = retroS, retroC = retroC)
}

# ------------------------------------------------------------
# 6. LOAD POSTERIOR DRAWS -- "b" directly (bounded model, no negation)
# ------------------------------------------------------------

build_stock_draws <- function(stock_name, n_draws) {
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  samples <- if (exists("sockeye_fits")) sockeye_fits[[stock_name]]$samples else
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds"))
  
  bnames <- grep("^b_", varnames(samples), value = TRUE)
  post <- as.data.frame(as.matrix(samples[, c("intercept", "b", bnames)]))
  idx <- sample(seq_len(nrow(post)), n_draws, replace = nrow(post) < n_draws)
  post <- post[idx, ]
  
  normalize_name <- function(x) tolower(gsub("[._]", "", x))
  bnames_stripped <- sub("^b_", "", bnames)
  
  draws <- tibble(draw = seq_len(n_draws), ra = post$intercept, rb = post$b)
  for (cv in extra_covariates) {
    match_idx <- which(normalize_name(bnames_stripped) == normalize_name(cv))
    if (length(match_idx) != 1) stop("Could not match covariate '", cv, "' for stock ", stock_name)
    draws[[cv]] <- post[[bnames[match_idx]]]
  }
  draws
}

stock_draws <- setNames(lapply(STOCKS_MC, build_stock_draws, n_draws = N_DRAWS), STOCKS_MC)

# ------------------------------------------------------------
# 7. MONTE CARLO LOOP -- STEP 3: X / Xmin, 3 scenarios x N_DRAWS x stocks
# ------------------------------------------------------------

lowpoint_results <- vector("list", length(STOCKS_MC) * N_DRAWS * length(SCENARIOS_MC))
counter <- 0
window_log <- list()

for (stock_name in STOCKS_MC) {
  
  draws_table      <- stock_draws[[stock_name]]
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
  
  for (scenario_label in names(SCENARIOS_MC)) {
    
    scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
    
    if (length(scenario_vars) == 0) {
      message("Skipping ", stock_name, " / ", scenario_label,
              " -- this stock's model has no covariates this scenario would freeze; ",
              "running it would be mechanically identical to no scenario at all.")
      next
    }
    
    w <- build_stock_window(stock_name, scenario_vars)
    if (is.null(w)) next
    
    window_log[[length(window_log) + 1]] <- tibble(
      Stock = stock_name, scenario = scenario_label,
      freeze_start = w$freeze_start, freeze_end = w$freeze_end,
      post_low_period_years = w$post_low_period_years,
      synthetic_years = max(0, w$freeze_end - w$natural_data_end)
    )
    
    X_start <- w$freeze_end - FINAL_PERIOD_YEARS + 1
    X_end   <- w$freeze_end
    
    for (d in seq_len(N_DRAWS)) {
      
      draw_row <- draws_table[d, ]
      terms <- list(
        ra = draw_row$ra, rb = draw_row$rb, sel_covs = extra_covariates,
        cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
      )
      
      result <- run_retro_model_lowpoint_draw(
        w$data, stock_name, terms, retroU_default, useretro_default, yrretro_default,
        scenario_vars = scenario_vars, scenario_start_year = w$freeze_start
      )
      
      X <- result %>%
        filter(Year >= X_start, Year <= X_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = stock_name, scenario = scenario_label, draw = d,
        Xmin = w$low_p$low_period_mean, X = X, ratio = X / w$low_p$low_period_mean,
        runaway = !is.finite(X) | (X / hist_max_this_stock) > RUNAWAY_MULTIPLE
      )
    }
  }
  
  message("done: ", stock_name)
}

window_log <- bind_rows(window_log)
print(window_log %>% arrange(desc(synthetic_years)), n = Inf)

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "sockeye_ratio_summary_synthetic_years.rds")

# ------------------------------------------------------------
# 8. DIAGNOSTICS
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
  cat("\nStill seeing NA ratios (should not happen with the unified gap-fill) for:\n")
  print(na_check, n = Inf)
}

zero_variance_check <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(sd_ratio = sd(ratio, na.rm = TRUE), .groups = "drop") %>%
  arrange(sd_ratio)
print(zero_variance_check, n = Inf)

# ------------------------------------------------------------
# 9. 40% Smsy BENCHMARK -- per-year actual covariates 1950-1970,
#    averaged output. Uses "b" directly, no sign flip.
# ------------------------------------------------------------

compute_smsy_benchmark <- function(stock_name) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  Xmin <- low_periods$low_period_mean[low_periods$Stock == stock_name]
  
  hp_data <- obs %>%
    filter(Stock == stock_name, Year >= SMSY_PERIOD_START_YEAR, Year <= SMSY_PERIOD_END_YEAR) %>%
    arrange(Year)
  if (length(extra_covariates) > 0) {
    hp_data <- hp_data %>% filter(if_all(all_of(extra_covariates), ~ !is.na(.)))
  }
  
  if (nrow(hp_data) == 0) {
    warning(stock_name, ": no usable years in ", SMSY_PERIOD_START_YEAR, "-", SMSY_PERIOD_END_YEAR,
            " -- NA benchmark.")
    return(tibble(Stock = stock_name, benchmark_ratio = NA_real_, n_years_used = 0))
  }
  
  samples <- if (exists("sockeye_fits")) sockeye_fits[[stock_name]]$samples else
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds"))
  
  bnames <- grep("^b_", varnames(samples), value = TRUE)
  post <- as.data.frame(as.matrix(samples[, c("intercept", "b", bnames)]))
  
  normalize_name <- function(x) tolower(gsub("[._]", "", x))
  bnames_stripped <- sub("^b_", "", bnames)
  
  n_years <- nrow(hp_data)
  a_mat <- outer(post$intercept, rep(1, n_years))
  
  for (cv in extra_covariates) {
    match_idx <- which(normalize_name(bnames_stripped) == normalize_name(cv))
    a_mat <- a_mat + outer(post[[bnames[match_idx]]], hp_data[[cv]])
  }
  
  beta_draw <- post$b
  
  Smsy_mat <- (a_mat / beta_draw) * (0.5 - 0.07 * a_mat)
  R_mat    <- (0.4 * Smsy_mat) * exp(a_mat - beta_draw * (0.4 * Smsy_mat))
  
  R_draw_avg <- rowMeans(R_mat, na.rm = TRUE)
  ratio_draw <- R_draw_avg / Xmin
  
  tibble(Stock = stock_name, benchmark_ratio = median(ratio_draw, na.rm = TRUE), n_years_used = n_years)
}

smsy_benchmarks <- bind_rows(lapply(STOCKS_MC, compute_smsy_benchmark))
print(smsy_benchmarks, n = Inf)

# ------------------------------------------------------------
# 10. SUMMARIZE + PLOT
# ------------------------------------------------------------

mc_summary_stats_ratio <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(n_used = n(), median = median(ratio, na.rm = TRUE),
            q05 = quantile(ratio, 0.10, na.rm = TRUE), 
            q95 = quantile(ratio, 0.90, na.rm = TRUE),
            .groups = "drop")
print(mc_summary_stats_ratio, n = Inf)

SCENARIO_COLORS2 <- c("Pinniped scenario" = "#4682B4", "SST scenario" = "#2E8B57", "Pink scenario" = "#FF4500")

stock_order <- mc_summary_stats_ratio %>%
  group_by(Stock) %>% summarise(m = max(median, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Stock)
mc_summary_stats_ratio <- mc_summary_stats_ratio %>% mutate(Stock = factor(Stock, levels = stock_order))
smsy_benchmarks <- smsy_benchmarks %>% mutate(Stock = factor(Stock, levels = stock_order))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.6), size = 0.5, linewidth = 1) +
  geom_segment(data = smsy_benchmarks,
               aes(x = as.numeric(Stock) - 0.4, xend = as.numeric(Stock) + 0.4,
                   y = benchmark_ratio, yend = benchmark_ratio, linetype = "40% Smsy"),
               inherit.aes = FALSE, color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  scale_linetype_manual(name = NULL, values = c("40% Smsy" = "dashed")) +
  scale_y_continuous(trans = "log", labels = scales::comma,
                     breaks = c(1, 10, 20, 30, 50, 100, 200)) +
  labs(x = NULL, y = "Recruits / min Recruits") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")

ggsave("figures/sockeye_ratio_v2.png", width = 13, height = 6.5, dpi = 600)

# histogram ---------------
lowpoint_summary %>%
  filter(!runaway) %>%
  ggplot(aes(ratio, fill = scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  facet_wrap(~ Stock, scales = "free") +
  labs(x = "Recruits / min Recruits", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_draws_hist.png", width = 12, height = 6.5, dpi = 600)

lowpoint_summary %>%
  filter(!runaway) %>%
  ggplot(aes(ratio, fill = scenario)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "white") +
  #geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  facet_wrap(scenario ~ Stock, scales = "free") +
  labs(x = "Recruits / min Recruits", y = "Number of draws") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("figures/sockeye_draws_hist_v2.png", width = 12, height = 6.5, dpi = 600)

# ------------------------------------------------------------
# 10. SUMMARIZE + PLOT
# ------------------------------------------------------------

mc_summary_stats_ratio <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(n_used = n(), median = median(ratio, na.rm = TRUE),
            q05 = quantile(ratio, 0.05, na.rm = TRUE), q95 = quantile(ratio, 0.95, na.rm = TRUE),
            .groups = "drop")
print(mc_summary_stats_ratio, n = Inf)

SCENARIO_COLORS2 <- c("Pinniped scenario" = "#4682B4", "SST scenario" = "#2E8B57", "Pink scenario" = "#FF4500")

stock_order <- mc_summary_stats_ratio %>%
  group_by(Stock) %>% summarise(m = max(median, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Stock)
mc_summary_stats_ratio <- mc_summary_stats_ratio %>% mutate(Stock = factor(Stock, levels = stock_order))
smsy_benchmarks <- smsy_benchmarks %>% mutate(Stock = factor(Stock, levels = stock_order))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.6), size = 0.5, linewidth = 1) +
  geom_segment(data = smsy_benchmarks,
               aes(x = as.numeric(Stock) - 0.4, xend = as.numeric(Stock) + 0.4,
                   y = benchmark_ratio, yend = benchmark_ratio, linetype = "40% Smsy"),
               inherit.aes = FALSE, color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  scale_linetype_manual(name = NULL, values = c("40% Smsy" = "dashed")) +
  scale_y_continuous(trans = "log", labels = scales::comma) +
  labs(x = NULL, y = "Recruits / min Recruits") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")

ggsave("figures/sockeye_ratio.png", width = 13, height = 6.5, dpi = 600)


MIN_YEARS_BEFORE_LOW_PERIOD <- 10

low_periods_v2 <- observed_returns %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(
    roll_mean = zoo::rollapply(RunJacks, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA),
    roll_start = Year - LOW_PERIOD_YEARS + 1,
    min_year = min(Year)
  ) %>%
  filter(!is.na(roll_mean), roll_start >= min_year + MIN_YEARS_BEFORE_LOW_PERIOD) %>%
  slice_min(roll_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, low_period_start = roll_start, low_period_end = Year, low_period_mean = roll_mean)

observed_returns %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(roll_mean = zoo::rollapply(RunJacks, width = LOW_PERIOD_YEARS, FUN = mean,
                                    align = "right", fill = NA)) %>%
  ggplot(aes(Year, roll_mean)) +
  geom_line(color = "#4682B4") +
  geom_rect(data = low_periods_v2,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.15) +
  geom_point(data = low_periods_v2, aes(x = low_period_end, y = low_period_mean), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Year", y = "10-year rolling mean return (recruits)") +
  theme_minimal()

ggsave("figures/sockeye_low_years.png", width = 13, height = 6.5, dpi = 600)


# ============================================================
# Diagnostic: visualize individual posterior draw trajectories through
# the freeze window (including any synthetic years), for one stock/
# scenario at a time -- colored by whether that draw got flagged
# runaway, with the freeze point and end-of-real-data marked.
#
# Run after sockeye_ratio_synthetic_years.R (needs STOCKS_MC,
# stock_draws, STOCK_COVARIATES, SCENARIOS_MC, build_stock_window(),
# run_retro_model_lowpoint_draw(), historical_max, RUNAWAY_MULTIPLE,
# retroU_default/useretro_default/yrretro_default).
# ============================================================

library(tidyverse)

plot_stock_trajectories <- function(stock_name, scenario_label, n_draws_to_plot = 40) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
  w <- build_stock_window(stock_name, scenario_vars)
  if (is.null(w)) stop("No window found for ", stock_name)
  
  draws_table <- stock_draws[[stock_name]]
  hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
  
  n_available <- nrow(draws_table)
  draw_ids <- sample(seq_len(n_available), min(n_draws_to_plot, n_available))
  
  traj_list <- lapply(draw_ids, function(d) {
    draw_row <- draws_table[d, ]
    terms <- list(
      ra = draw_row$ra, rb = draw_row$rb, sel_covs = extra_covariates,
      cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
    )
    result <- run_retro_model_lowpoint_draw(
      w$data, stock_name, terms, retroU_default, useretro_default, yrretro_default,
      scenario_vars = scenario_vars, scenario_start_year = w$freeze_start
    )
    final_X <- result %>%
      filter(Year >= w$freeze_end - 9, Year <= w$freeze_end) %>%
      summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
    is_runaway <- !is.finite(final_X) | (final_X / hist_max_this_stock) > RUNAWAY_MULTIPLE
    
    result %>% select(Year, retroR) %>% mutate(draw = d, runaway = is_runaway)
  })
  
  traj <- bind_rows(traj_list)
  
  n_runaway <- traj %>% distinct(draw, runaway) %>% pull(runaway) %>% sum()
  cat(stock_name, "/", scenario_label, ": ", n_runaway, "of", length(draw_ids),
      "plotted draws flagged runaway\n")
  cat("Freeze starts:", w$freeze_start, " | Real data ends:", w$natural_data_end,
      " | Freeze/window ends:", w$freeze_end,
      " | Synthetic years:", max(0, w$freeze_end - w$natural_data_end), "\n")
  
  ggplot(traj, aes(Year, retroR, group = draw, color = runaway)) +
    geom_line(alpha = 0.5, linewidth = 0.4) +
    geom_vline(xintercept = w$freeze_start, linetype = "dashed", color = "black") +
    geom_vline(xintercept = w$natural_data_end, linetype = "dotted", color = "grey30") +
    annotate("text", x = w$freeze_start, y = Inf, label = "freeze starts", vjust = 1.5, hjust = -0.05, size = 3) +
    annotate("text", x = w$natural_data_end, y = Inf, label = "real data ends", vjust = 3, hjust = 1.05, size = 3) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c(`FALSE` = "#4682B4", `TRUE` = "firebrick"),
                       labels = c(`FALSE` = "Not runaway", `TRUE` = "Runaway"), name = NULL) +
    labs(x = "Year", y = "Reconstructed return (retroR, log scale)",
         title = paste0(stock_name, " -- ", scenario_label),
         subtitle = "Dashed line = freeze starts (year after low window) | Dotted line = end of real data (synthetic years begin here if freeze extends further)") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Example usage -- swap in whichever stock/scenario had a high
# runaway_fraction in your last run:
plot_stock_trajectories("Birkenhead", "Pinniped scenario")
plot_stock_trajectories("Cultus", "Pinniped scenario")
# ggsave("figures/diagnostic_gates_pinniped_trajectories.png", width = 9, height = 6, dpi = 600)


# ============================================================
# Diagnostic: trace every intermediate quantity through the recursion
# for ONE draw, around the real/synthetic boundary -- to find exactly
# where (and in which quantity) things first go wrong, instead of
# inferring it from the plot.
#
# Run after sockeye_ratio_synthetic_years.R (needs the same objects as
# the trajectory diagnostic).
# ============================================================

library(tidyverse)

trace_boundary <- function(stock_name, scenario_label, draw_id = 1) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
  w <- build_stock_window(stock_name, scenario_vars)
  if (is.null(w)) stop("No window found for ", stock_name)
  
  cat("=== build_stock_window() output ===\n")
  cat("natural_data_end:", w$natural_data_end, " | freeze_start:", w$freeze_start,
      " | freeze_end:", w$freeze_end, "\n")
  cat("Rows in w$data:", nrow(w$data), " | Max Year in w$data:", max(w$data$Year), "\n")
  cat("Does w$data actually contain rows past natural_data_end? ",
      any(w$data$Year > w$natural_data_end), "\n\n")
  
  cat("=== w$data around the boundary (raw covariates) ===\n")
  boundary_cols <- c("Year", extra_covariates, paste0(intersect(extra_covariates, scenario_vars), "_scenario"))
  boundary_cols <- intersect(boundary_cols, names(w$data))
  print(w$data %>% filter(Year >= w$natural_data_end - 5, Year <= w$natural_data_end + 5) %>%
          select(all_of(boundary_cols)), n = Inf)
  
  draws_table <- stock_draws[[stock_name]]
  draw_row <- draws_table[draw_id, ]
  terms <- list(
    ra = draw_row$ra, rb = draw_row$rb, sel_covs = extra_covariates,
    cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
  )
  cat("\n=== Draw", draw_id, "terms ===\n")
  print(terms)
  
  # Re-run the SAME steps run_retro_model_lowpoint_draw() does, but keep
  # every intermediate column instead of discarding them, so we can see
  # exactly where things break.
  dat <- w$data %>% arrange(Year)
  obs2 <- dat %>%
    mutate(
      RunJacks = RunSize - JackEscapement,
      Catch = rowSums(cbind(BelowMissionC, AboveMissionC), na.rm = TRUE),
      Ut_obs = pmin(0.95, Catch / RunJacks),
      Ut_obs = if_else(stock_name == "Late Shuswap" & Year == 2012, pmin(0.7, Catch / RunJacks), Ut_obs),
      ENS = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(AdultReturn = lead(RunJacks, n = LAG_YEARS), lnR_S = log(AdultReturn / AdultEscapement))
  
  real_ens <- obs2$ENS[!is.na(obs2$ENS)]
  ens_fill <- mean(tail(real_ens, 10), na.rm = TRUE)
  cat("\nens_fill value:", ens_fill, " (NaN here means NO real ENS values were found at all)\n")
  obs2 <- obs2 %>% mutate(ENS = if_else(is.na(ENS), ens_fill, ENS))
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term_from_year(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars, w$freeze_start)
  
  obs3 <- obs2 %>%
    mutate(cov_term = cov_term, cov_term_proj = cov_term_proj,
           wt_raw = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term))
  
  real_wt <- obs3$wt_raw[!is.na(obs3$wt_raw)]
  na_idx  <- which(is.na(obs3$wt_raw))
  wt_filled <- obs3$wt_raw
  set.seed(1)  # deterministic for this diagnostic run
  if (length(na_idx) > 0) {
    wt_filled[na_idx] <- if (length(real_wt) == 0) 0 else sample(real_wt, length(na_idx), replace = TRUE)
  }
  obs3$wt <- wt_filled
  
  cat("\n=== ENS / cov_term_proj / wt around the boundary (should all be non-NA) ===\n")
  print(obs3 %>% filter(Year >= w$natural_data_end - 5, Year <= w$natural_data_end + 5) %>%
          select(Year, ENS, cov_term_proj, wt), n = Inf)
  
  n <- nrow(obs3)
  retroR <- rep(NA_real_, n); retro_lnRS <- rep(NA_real_, n)
  retroS <- rep(NA_real_, n); retroC <- rep(NA_real_, n)
  retroU_vec <- if (useretro_default) ifelse(obs3$Year >= yrretro_default, retroU_default, obs3$Ut_obs) else obs3$Ut_obs
  
  NUMERIC_CEILING <- 1e15
  
  retroR[1:LAG_YEARS] <- obs3$RunJacks[1:LAG_YEARS]
  retroS[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * (1 - retroU_vec[1:LAG_YEARS]) * obs3$ENS[1:LAG_YEARS]
  retroC[1:LAG_YEARS] <- retroR[1:LAG_YEARS] * retroU_vec[1:LAG_YEARS]
  
  first_bad <- NA_integer_
  for (i in (LAG_YEARS + 1):n) {
    j <- i - LAG_YEARS
    retro_lnRS[j] <- terms$ra - terms$rb * retroS[j] + obs3$cov_term_proj[j] + obs3$wt[j]
    retroR[i] <- pmin(retroS[j] * exp(retro_lnRS[j]), NUMERIC_CEILING)
    retroS[i] <- pmin(retroR[i] * (1 - retroU_vec[i]) * obs3$ENS[i], NUMERIC_CEILING)
    retroC[i] <- retroR[i] * retroU_vec[i]
    if (is.na(first_bad) && (!is.finite(retroR[i]) || !is.finite(retroS[i]))) first_bad <- i
  }
  
  cat("\n=== First row where retroR or retroS becomes non-finite ===\n")
  if (is.na(first_bad)) {
    cat("None -- the full recursion stayed finite through Year", obs3$Year[n], "\n")
  } else {
    cat("Row", first_bad, "-> Year", obs3$Year[first_bad], "\n")
    cat("j =", first_bad - LAG_YEARS, "-> that row's Year =", obs3$Year[first_bad - LAG_YEARS], "\n")
    cat("retroS[j]:", retroS[first_bad - LAG_YEARS], " | retro_lnRS[j]:", retro_lnRS[first_bad - LAG_YEARS], "\n")
    cat("cov_term_proj[j]:", obs3$cov_term_proj[first_bad - LAG_YEARS],
        " | wt[j]:", obs3$wt[first_bad - LAG_YEARS], "\n")
  }
  
  cat("\n=== Full retroR/retroS trajectory around the boundary ===\n")
  print(tibble(Year = obs3$Year, retroR = retroR, retroS = retroS) %>%
          filter(Year >= w$natural_data_end - 5, Year <= w$natural_data_end + 10), n = Inf)
  
  invisible(list(obs3 = obs3, retroR = retroR, retroS = retroS, retro_lnRS = retro_lnRS))
}

trace_boundary("Cultus", "Pinniped scenario", draw_id = 1)

# ============================================================
# Diagnostic: the 10-year low-abundance period per stock (real observed
# Return, lowest mean), and -- critically -- how many REAL years remain
# between the end of that window and the end of the actual data. This
# is what determines whether the 12-year/3-generation freeze requirement
# can be met with real data alone, or needs synthetic projection.
# ============================================================

library(tidyverse)

LOW_PERIOD_YEARS <- 10

observed_returns <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  select(Stock, Year, RunJacks)

low_periods_check <- observed_returns %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  filter(n() >= LOW_PERIOD_YEARS) %>%
  mutate(roll_mean = zoo::rollapply(RunJacks, width = LOW_PERIOD_YEARS, FUN = mean, align = "right", fill = NA)) %>%
  filter(!is.na(roll_mean)) %>%
  slice_min(roll_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, low_period_start = Year - LOW_PERIOD_YEARS + 1, low_period_end = Year,
            low_period_mean = roll_mean)

natural_data_end_by_stock <- observed_returns %>%
  group_by(Stock) %>%
  summarise(natural_data_end = max(Year, na.rm = TRUE), .groups = "drop")

low_period_diagnostic <- low_periods_check %>%
  left_join(natural_data_end_by_stock, by = "Stock") %>%
  mutate(
    real_years_after_low_period = natural_data_end - low_period_end,
    freeze_start = low_period_end + 1,
    is_current_low_period = (low_period_end == natural_data_end)
  ) %>%
  arrange(real_years_after_low_period)

print(low_period_diagnostic, n = Inf)

cat("\n=== Stocks whose low period is the MOST RECENT 10 years (real_years_after_low_period = 0) ===\n")
cat("These have zero real years before projection has to start -- worth confirming this is a\n")
cat("genuine ongoing decline in the data, not a data artifact.\n")
print(low_period_diagnostic %>% filter(real_years_after_low_period == 0), n = Inf)

