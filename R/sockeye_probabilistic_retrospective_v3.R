# ============================================================
# FULL SCRIPT -- sockeye X/Xmin recovery ratio + 40% Smsy benchmark,
# all three scenarios (Pinniped, SST, Pink)
#
# v2 CHANGES:
#  1. FIXED A REAL SIGN BUG from the previous version: rb/beta_draw were
#     read as post$b_spawners directly, with a comment claiming they
#     were "already positive" -- they weren't (b_spawners was a raw,
#     typically-negative coefficient, same convention as coho). This
#     flipped the sign of density dependence and broke the Smsy
#     denominator. Now reads sockeye_jags_models_bounded_b.R's "b"
#     parameter directly (no negation needed -- it's already the
#     positive Ricker slope by construction, thanks to the bounded
#     lognormal prior).
#  2. Uses the new bounded-b JAGS models (sockeye_jags_models_bounded_b.R)
#     and their updated alpha<0.05 covariate sets (Bowron and Scotch
#     changed -- see that script's header).
#
# DEPENDENCIES (must already exist in your session, from your cleaned
# sockeye retrospective script): obs, model_actual_retro,
# model_pinniped_retro, model_sst_retro, model_pink_retro,
# compute_cov_term, STOCKS, LAG_YEARS, retroU_default, useretro_default,
# yrretro_default, FIT_YEARS. Also sockeye_fits (or the saved
# sockeye_<stock>_posterior_samples.rds files) from
# sockeye_jags_models_bounded_b.R -- NOT the old sockeye_jags_models.R.
#
# SECTIONS:
#   1. low_periods setup (this was missing before -- GT-year trailing
#      window of lowest observed Return, per stock)
#   2. Diagnostic table + plot of the low periods found
#   3. Posterior draws + run_retro_model_draw()
#   4. Monte Carlo loop: X/Xmin ratio, 3 scenarios x N_DRAWS x stocks
#   5. 40% Smsy benchmark (per-year actual covariates through 1970,
#      averaged output)
#   6. Plot: X/Xmin pointrange (dodged by scenario) + Smsy benchmark
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50   # start small -- bump to 5000 once this runs cleanly
FINAL_PERIOD_YEARS <- 10
SMSY_PERIOD_END_YEAR <- 1970
RUNAWAY_MULTIPLE <- 10

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

# ------------------------------------------------------------
# 1. LOW_PERIODS SETUP -- was missing before. GT-year (generation-time)
#    trailing window of lowest mean OBSERVED Return, per stock.
# ------------------------------------------------------------

GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5

trailing_mean <- function(x, width) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (i >= width) {
      window <- x[(i - width + 1):i]
      if (all(!is.na(window))) out[i] <- mean(window)
    }
  }
  out
}

returns_by_stock_scenario <- bind_rows(
  model_actual_retro   %>% select(Stock, Year, Return = RunJacks) %>% mutate(scenario = "Observed"),
  model_pinniped_retro %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pinniped scenario"),
  model_sst_retro      %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "SST scenario"),
  model_pink_retro     %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pink scenario")
)

low_periods <- returns_by_stock_scenario %>%
  filter(scenario == "Observed", is.finite(Return)) %>%
  left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(obs_gen_mean = trailing_mean(Return, width = first(GT))) %>%
  filter(!is.na(obs_gen_mean)) %>%
  slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, GT, low_period_end = Year, low_period_start = Year - GT + 1,
            low_period_mean = obs_gen_mean)

# ------------------------------------------------------------
# 2. DIAGNOSTIC: table + plot of the low periods found
# ------------------------------------------------------------

print(low_periods %>% select(Stock, GT, low_period_start, low_period_end, low_period_mean) %>% arrange(Stock))

returns_by_stock_scenario %>%
  filter(scenario == "Observed", is.finite(Return)) %>%
  left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(roll_mean = trailing_mean(Return, width = first(GT))) %>%
  ggplot(aes(Year, roll_mean)) +
  geom_line(color = "#4682B4") +
  geom_point(data = low_periods, aes(x = low_period_end, y = low_period_mean), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "Generation-length (GT-year) rolling mean return",
       title = "Low-abundance window per stock (red = selected low period)") +
  theme_minimal()

# ------------------------------------------------------------
# Check for stocks missing from `obs` before proceeding
# ------------------------------------------------------------

available_stocks <- intersect(STOCKS_MC, unique(obs$Stock))
missing_stocks <- setdiff(STOCKS_MC, available_stocks)
if (length(missing_stocks) > 0) {
  message("Not found in `obs`, excluding: ", paste(missing_stocks, collapse = ", "),
          " -- worth checking whether this is a genuine data gap or a name mismatch.")
}
STOCKS_MC <- available_stocks

historical_max <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 3. run_retro_model_draw() + posterior draws
# ------------------------------------------------------------

run_retro_model_draw <- function(dat, stock_name, terms, retroU, useretro, yrretro,
                                 scenario_vars = character(0)) {
  
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
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars)
  
  obs3 <- obs2 %>%
    mutate(cov_term = cov_term, cov_term_proj = cov_term_proj,
           wt = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term))
  
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
  
  # FIX: "b" is the bounded Ricker slope (sockeye_jags_models_bounded_b.R)
  # -- already positive by construction (lognormal, truncated below at
  # zero), unlike the old unbounded "b_spawners" which was a raw,
  # typically-negative coefficient. NO negation here.
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
# 4. MONTE CARLO LOOP -- X/Xmin, 3 scenarios x N_DRAWS x stocks
# ------------------------------------------------------------

lowpoint_results <- vector("list", length(STOCKS_MC) * N_DRAWS * length(SCENARIOS_MC))
counter <- 0

for (stock_name in STOCKS_MC) {
  
  stock_data       <- obs %>% filter(Stock == stock_name)
  draws_table      <- stock_draws[[stock_name]]
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  low_p            <- low_periods %>% filter(Stock == stock_name)
  hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
  
  final_period_end   <- max(stock_data$Year, na.rm = TRUE)
  final_period_start <- final_period_end - FINAL_PERIOD_YEARS + 1
  
  for (scenario_label in names(SCENARIOS_MC)) {
    
    scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
    
    for (d in seq_len(N_DRAWS)) {
      
      draw_row <- draws_table[d, ]
      terms <- list(
        ra = draw_row$ra, rb = draw_row$rb, sel_covs = extra_covariates,
        cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
      )
      
      result <- run_retro_model_draw(stock_data, stock_name, terms, retroU_default, useretro_default,
                                     yrretro_default, scenario_vars = scenario_vars)
      
      X <- result %>%
        filter(Year >= final_period_start, Year <= final_period_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = stock_name, scenario = scenario_label, draw = d,
        Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean,
        runaway = !is.finite(X) | (X / hist_max_this_stock) > RUNAWAY_MULTIPLE
      )
    }
  }
  
  message("done: ", stock_name)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "sockeye_ratio_summary_all_scenarios.rds")

runaway_fraction <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_runaway = sum(runaway, na.rm = TRUE),
            pct_runaway = 100 * n_runaway / n_draws, .groups = "drop") %>%
  arrange(desc(pct_runaway))
print(runaway_fraction, n = Inf)

mc_summary_stats_ratio <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.05, na.rm = TRUE),
    q95    = quantile(ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_ratio, n = Inf)

# ------------------------------------------------------------
# 5. 40% Smsy BENCHMARK -- per-year actual covariates through 1970,
#    averaged output, one value per stock
# ------------------------------------------------------------

compute_smsy_benchmark <- function(stock_name) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  Xmin <- low_periods$low_period_mean[low_periods$Stock == stock_name]
  
  hp_data <- obs %>%
    filter(Stock == stock_name, Year <= SMSY_PERIOD_END_YEAR) %>%
    arrange(Year)
  if (length(extra_covariates) > 0) {
    hp_data <- hp_data %>% filter(if_all(all_of(extra_covariates), ~ !is.na(.)))
  }
  
  if (nrow(hp_data) == 0) {
    warning(stock_name, ": no usable years <= ", SMSY_PERIOD_END_YEAR, " -- NA benchmark.")
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
  
  # FIX: "b" (sockeye_jags_models_bounded_b.R) is already the positive
  # Ricker slope by construction -- no sign flip needed, unlike the old
  # unbounded "b_spawners" this previously (incorrectly) read as-is
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
# 6. PLOT: X/Xmin pointrange (dodged by scenario) + 40% Smsy benchmark
# ------------------------------------------------------------

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

stock_order <- mc_summary_stats_ratio %>%
  group_by(Stock) %>% summarise(m = max(median, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Stock)

mc_summary_stats_ratio <- mc_summary_stats_ratio %>% mutate(Stock = factor(Stock, levels = stock_order))
smsy_benchmarks <- smsy_benchmarks %>% mutate(Stock = factor(Stock, levels = stock_order))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 1) +
  geom_segment(data = smsy_benchmarks,
               aes(x = as.numeric(Stock) - 0.4, xend = as.numeric(Stock) + 0.4,
                   y = benchmark_ratio, yend = benchmark_ratio, linetype = "40% Smsy"),
               inherit.aes = FALSE, color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  scale_linetype_manual(name = NULL, values = c("40% Smsy" = "dashed")) +
  scale_y_continuous(trans = "log", labels = scales::comma) +
  labs(x = NULL, y = "X / Xmin (posterior median, 90% credible interval, natural log scale)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")

ggsave("figures/sockeye_ratio_pointrange_with_smsy_benchmark.png", width = 13, height = 6.5, dpi = 600)

############# NEW ##############################################



























# ============================================================
# FULL SCRIPT -- sockeye X/Xmin recovery ratio + 40% Smsy benchmark,
# all three scenarios (Pinniped, SST, Pink)
#
# v2 CHANGES:
#  1. FIXED A REAL SIGN BUG from the previous version: rb/beta_draw were
#     read as post$b_spawners directly, with a comment claiming they
#     were "already positive" -- they weren't (b_spawners was a raw,
#     typically-negative coefficient, same convention as coho). This
#     flipped the sign of density dependence and broke the Smsy
#     denominator. Now reads sockeye_jags_models_bounded_b.R's "b"
#     parameter directly (no negation needed -- it's already the
#     positive Ricker slope by construction, thanks to the bounded
#     lognormal prior).
#  2. Uses the new bounded-b JAGS models (sockeye_jags_models_bounded_b.R)
#     and their updated alpha<0.05 covariate sets (Bowron and Scotch
#     changed -- see that script's header).
#
# DEPENDENCIES (must already exist in your session, from your cleaned
# sockeye retrospective script): obs, model_actual_retro,
# model_pinniped_retro, model_sst_retro, model_pink_retro,
# compute_cov_term, STOCKS, LAG_YEARS, retroU_default, useretro_default,
# yrretro_default, FIT_YEARS. Also sockeye_fits (or the saved
# sockeye_<stock>_posterior_samples.rds files) from
# sockeye_jags_models_bounded_b.R -- NOT the old sockeye_jags_models.R.
#
# SECTIONS:
#   1. low_periods setup: COMMON 1994-2003 window across all stocks
#      (matching the coho approach), not each stock's own GT-based
#      window of lowest observed Return, per stock)
#   2. Diagnostic table + plot of the low periods found
#   3. Posterior draws + run_retro_model_draw()
#   4. Monte Carlo loop: X/Xmin ratio, 3 scenarios x N_DRAWS x stocks
#   5. 40% Smsy benchmark (per-year actual covariates through 1970,
#      averaged output)
#   6. Plot: X/Xmin pointrange (dodged by scenario) + Smsy benchmark
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 50   # start small -- bump to 5000 once this runs cleanly
FINAL_PERIOD_YEARS <- 10
SMSY_PERIOD_START_YEAR <- 1950   # high-productivity window: 1950-1970
SMSY_PERIOD_END_YEAR <- 1970
RUNAWAY_MULTIPLE <- 10

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

# ------------------------------------------------------------
# 1. LOW_PERIODS SETUP -- COMMON WINDOW ACROSS ALL STOCKS, matching the
#    coho approach: instead of each stock's own GT-year trailing-window
#    minimum, use a single shared 1994-2003 window (freeze/reference
#    point 2004) for every stock. This is a deliberate experiment ("try
#    it out") to see how results compare to the per-stock GT-based
#    version -- not a replacement decision yet. Xmin still comes from
#    real observed Return, just over the SAME calendar years for every
#    stock rather than each stock's own worst generation-length window.
# ------------------------------------------------------------

GENERATION_TIME <- setNames(rep(4, length(STOCKS)), STOCKS)
GENERATION_TIME["Pitt"] <- 5

trailing_mean <- function(x, width) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (i >= width) {
      window <- x[(i - width + 1):i]
      if (all(!is.na(window))) out[i] <- mean(window)
    }
  }
  out
}

returns_by_stock_scenario <- bind_rows(
  model_actual_retro   %>% select(Stock, Year, Return = RunJacks) %>% mutate(scenario = "Observed"),
  model_pinniped_retro %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pinniped scenario"),
  model_sst_retro      %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "SST scenario"),
  model_pink_retro     %>% select(Stock, Year, Return = retroR)   %>% mutate(scenario = "Pink scenario")
)

COMMON_LOW_PERIOD_START <- 1994
COMMON_LOW_PERIOD_END   <- 2003

low_periods <- returns_by_stock_scenario %>%
  filter(scenario == "Observed", is.finite(Return),
         Year >= COMMON_LOW_PERIOD_START, Year <= COMMON_LOW_PERIOD_END) %>%
  group_by(Stock) %>%
  summarise(low_period_mean = mean(Return, na.rm = TRUE), n_years_used = n(), .groups = "drop") %>%
  mutate(low_period_start = COMMON_LOW_PERIOD_START, low_period_end = COMMON_LOW_PERIOD_END)

# ------------------------------------------------------------
# 2. DIAGNOSTIC: table + plot of the common low period vs. each stock's
#    own GT-based rolling trajectory, for comparison
# ------------------------------------------------------------

print(low_periods %>% arrange(Stock))

if (any(low_periods$n_years_used < (COMMON_LOW_PERIOD_END - COMMON_LOW_PERIOD_START + 1))) {
  warning("Some stocks have fewer than 10 real years in the common window -- check n_years_used:")
  print(low_periods %>% filter(n_years_used < (COMMON_LOW_PERIOD_END - COMMON_LOW_PERIOD_START + 1)))
}

returns_by_stock_scenario %>%
  filter(scenario == "Observed", is.finite(Return)) %>%
  left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(roll_mean = trailing_mean(Return, width = first(GT))) %>%
  ggplot(aes(Year, roll_mean)) +
  geom_line(color = "#4682B4") +
  geom_rect(data = low_periods,
            aes(xmin = low_period_start, xmax = low_period_end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "red", alpha = 0.15) +
  geom_point(data = low_periods, aes(x = low_period_end, y = low_period_mean), color = "red", size = 2) +
  facet_wrap(~ Stock, scales = "free_y") +
  labs(x = "Year", y = "Generation-length (GT-year) rolling mean return",
       title = "Common low-period window (1994-2003, red) vs. each stock's own GT-based trajectory") +
  theme_minimal()

# ------------------------------------------------------------
# Check for stocks missing from `obs` before proceeding
# ------------------------------------------------------------

available_stocks <- intersect(STOCKS_MC, unique(obs$Stock))
missing_stocks <- setdiff(STOCKS_MC, available_stocks)
if (length(missing_stocks) > 0) {
  message("Not found in `obs`, excluding: ", paste(missing_stocks, collapse = ", "),
          " -- worth checking whether this is a genuine data gap or a name mismatch.")
}
STOCKS_MC <- available_stocks

historical_max <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 3. run_retro_model_draw() + posterior draws
# ------------------------------------------------------------

run_retro_model_draw <- function(dat, stock_name, terms, retroU, useretro, yrretro,
                                 scenario_vars = character(0)) {
  
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
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars)
  
  obs3 <- obs2 %>%
    mutate(cov_term = cov_term, cov_term_proj = cov_term_proj,
           wt = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term))
  
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
  
  # FIX: "b" is the bounded Ricker slope (sockeye_jags_models_bounded_b.R)
  # -- already positive by construction (lognormal, truncated below at
  # zero), unlike the old unbounded "b_spawners" which was a raw,
  # typically-negative coefficient. NO negation here.
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
# 4. MONTE CARLO LOOP -- X/Xmin, 3 scenarios x N_DRAWS x stocks
# ------------------------------------------------------------

lowpoint_results <- vector("list", length(STOCKS_MC) * N_DRAWS * length(SCENARIOS_MC))
counter <- 0

for (stock_name in STOCKS_MC) {
  
  stock_data       <- obs %>% filter(Stock == stock_name)
  draws_table      <- stock_draws[[stock_name]]
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  low_p            <- low_periods %>% filter(Stock == stock_name)
  hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
  
  final_period_end   <- max(stock_data$Year, na.rm = TRUE)
  final_period_start <- final_period_end - FINAL_PERIOD_YEARS + 1
  
  for (scenario_label in names(SCENARIOS_MC)) {
    
    scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
    
    for (d in seq_len(N_DRAWS)) {
      
      draw_row <- draws_table[d, ]
      terms <- list(
        ra = draw_row$ra, rb = draw_row$rb, sel_covs = extra_covariates,
        cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
      )
      
      result <- run_retro_model_draw(stock_data, stock_name, terms, retroU_default, useretro_default,
                                     yrretro_default, scenario_vars = scenario_vars)
      
      X <- result %>%
        filter(Year >= final_period_start, Year <= final_period_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>% pull(m)
      
      counter <- counter + 1
      lowpoint_results[[counter]] <- tibble(
        Stock = stock_name, scenario = scenario_label, draw = d,
        Xmin = low_p$low_period_mean, X = X, ratio = X / low_p$low_period_mean,
        runaway = !is.finite(X) | (X / hist_max_this_stock) > RUNAWAY_MULTIPLE
      )
    }
  }
  
  message("done: ", stock_name)
}

lowpoint_summary <- bind_rows(lowpoint_results)
saveRDS(lowpoint_summary, "sockeye_ratio_summary_all_scenarios.rds")

runaway_fraction <- lowpoint_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_runaway = sum(runaway, na.rm = TRUE),
            pct_runaway = 100 * n_runaway / n_draws, .groups = "drop") %>%
  arrange(desc(pct_runaway))
print(runaway_fraction, n = Inf)

mc_summary_stats_ratio <- lowpoint_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    median = median(ratio, na.rm = TRUE),
    q05    = quantile(ratio, 0.05, na.rm = TRUE),
    q95    = quantile(ratio, 0.95, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_ratio, n = Inf)

# ------------------------------------------------------------
# 5. 40% Smsy BENCHMARK -- per-year actual covariates through 1970,
#    averaged output, one value per stock
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
  
  # FIX: "b" (sockeye_jags_models_bounded_b.R) is already the positive
  # Ricker slope by construction -- no sign flip needed, unlike the old
  # unbounded "b_spawners" this previously (incorrectly) read as-is
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
# 6. PLOT: X/Xmin pointrange (dodged by scenario) + 40% Smsy benchmark
# ------------------------------------------------------------

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

stock_order <- mc_summary_stats_ratio %>%
  group_by(Stock) %>% summarise(m = max(median, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(m)) %>% pull(Stock)

mc_summary_stats_ratio <- mc_summary_stats_ratio %>% mutate(Stock = factor(Stock, levels = stock_order))
smsy_benchmarks <- smsy_benchmarks %>% mutate(Stock = factor(Stock, levels = stock_order))

ggplot(mc_summary_stats_ratio, aes(x = Stock, y = median, color = scenario)) +
  geom_pointrange(aes(ymin = q05, ymax = q95), position = position_dodge(width = 0.6),
                  size = 0.5, linewidth = 1) +
  geom_segment(data = smsy_benchmarks,
               aes(x = as.numeric(Stock) - 0.4, xend = as.numeric(Stock) + 0.4,
                   y = benchmark_ratio, yend = benchmark_ratio, linetype = "40% Smsy"),
               inherit.aes = FALSE, color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_color_manual(values = SCENARIO_COLORS2, name = "Scenario driver") +
  scale_linetype_manual(name = NULL, values = c("40% Smsy" = "dashed")) +
  scale_y_continuous(trans = "log", labels = scales::comma) +
  labs(x = NULL, y = "X / Xmin (posterior median, 90% credible interval, natural log scale)",
       title = paste0("Sockeye recovery ratio (X / Xmin) by scenario, vs. 40% Smsy benchmark\n",
                      "(per-year actual covariates, ", SMSY_PERIOD_START_YEAR, "-", SMSY_PERIOD_END_YEAR, ")")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "bottom")

ggsave("figures/sockeye_ratio_pointrange_with_smsy_benchmark.png", width = 13, height = 6.5, dpi = 600)


























