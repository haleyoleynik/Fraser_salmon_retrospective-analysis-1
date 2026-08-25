# ============================================================
# Probabilistic (Monte Carlo) retrospective analysis -- Fraser sockeye,
# 15 stocks (Late Shuswap excluded -- no usable top model)
#
# For each of N_DRAWS posterior draws per stock, reconstructs the
# Pinniped / SST / Pink scenarios (run_retro_model_draw(), a version of
# your run_retro_model() that takes `terms` directly instead of looking
# them up from the fixed top_models table) and computes % increase in
# mean abundance over that stock's low-abundance period vs. observed --
# the same metric as pct_increase_over_low_period in your cleaned
# script's COSEWIC section.
#
# SIMPLER than steelhead/coho in one respect: that metric compares each
# scenario directly against raw historical RunJacks (real data), not
# against a separately-reconstructed "actual, no-scenario" model run --
# so low_periods and the Observed baseline are fixed real data, computed
# ONCE, not per draw. No need to run an "actual" model per draw at all.
#
# DEPENDENCIES -- run your cleaned sockeye retrospective script first, so
# `obs`, `COVARIATE_COLS`, `retroU_default`, `useretro_default`,
# `yrretro_default`, `LAG_YEARS`, `FIT_YEARS` already exist. Also run
# sockeye_jags_models.R (or readRDS its saved posterior samples) so
# `sockeye_fits` / the `sockeye_<stock>_posterior_samples.rds` files
# exist.
#
# CAVEATS:
#  1. Late Shuswap is excluded throughout (STOCKS_MC below), since it has
#     no dredge top model to fit.
#  2. Uses retroU_default/useretro_default/yrretro_default (the *_retro
#     harvest assumption), matching how pct_increase_over_low_period was
#     originally built from model_pinniped_retro/model_sst_retro/
#     model_pink_retro -- not the *_hist catch-lost version.
#  3. Same safety net as steelhead/coho: reports BOTH % increase and an
#     absolute increase in mean abundance, plus a flag for draws where
#     the observed low-period baseline is near zero, rather than trusting
#     % alone. Given RunJacks is real historical run size (not a
#     model-floored reconstruction), near-zero baselines should be rare,
#     but worth checking rather than assuming.
#  4. Runtime: 15 stocks x N_DRAWS x 3 scenarios x ~65-year single-lag
#     loop per run. No age-structure/multi-lag chain (unlike steelhead),
#     so each run is cheap -- similar cost profile to coho.
# ============================================================

library(tidyverse)
library(coda)

set.seed(2026)
N_DRAWS <- 5000

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

jags_name <- function(cv) gsub("\\.", "_", cv)

# ------------------------------------------------------------
# 1. run_retro_model_draw() -- your run_retro_model(), modified to take
#    `terms` (a posterior draw's ra/rb/sel_covs/cov_coefs) directly
#    instead of calling get_top_model_terms(stock_name, fit_df). The
#    process-error residual (wt) is still recomputed fresh against THIS
#    draw's terms -- everything else is identical to your original.
# ------------------------------------------------------------

run_retro_model_draw <- function(dat, stock_name, terms, retroU, useretro, yrretro,
                                 scenario_vars = character(0)) {
  
  dat <- dat %>% arrange(Year)
  
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
      ENS = pmin(1, pmax(0.0001, AdultEscapement / RunJacks / (1 - Ut_obs))),
      migmort = 1 - ENS
    ) %>%
    mutate(
      AdultReturn = lead(RunJacks, n = LAG_YEARS),
      lnR_S = log(AdultReturn / AdultEscapement)
    )
  
  cov_term      <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs)
  cov_term_proj <- compute_cov_term(obs2, terms$sel_covs, terms$cov_coefs, scenario_vars)
  
  obs3 <- obs2 %>%
    mutate(
      cov_term      = cov_term,
      cov_term_proj = cov_term_proj,
      wt = lnR_S - (terms$ra - terms$rb * AdultEscapement + cov_term)
    )
  
  n <- nrow(obs3)
  retroR     <- rep(NA_real_, n)
  retro_lnRS <- rep(NA_real_, n)
  retroS     <- rep(NA_real_, n)
  retroC     <- rep(NA_real_, n)
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
  
  obs3 %>%
    mutate(retroR = retroR, retro_lnRS = retro_lnRS, retroU = retroU_vec,
           retroS = retroS, retroC = retroC)
}

# ------------------------------------------------------------
# 2. LOAD POSTERIOR SAMPLES AND SUBSAMPLE N_DRAWS PER STOCK
# ------------------------------------------------------------

build_stock_draws <- function(stock_name, n_draws) {
  
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  
  samples <- if (exists("sockeye_fits")) {
    sockeye_fits[[stock_name]]$samples
  } else {
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock_name)), "_posterior_samples.rds"))
  }
  
  bnames <- grep("^b_", varnames(samples), value = TRUE)
  post <- as.data.frame(as.matrix(samples[, c("intercept", bnames)]))
  
  idx <- sample(seq_len(nrow(post)), n_draws, replace = nrow(post) < n_draws)
  post <- post[idx, ]
  
  # match each b_<name> to its covariate, case/dot/underscore-insensitive
  # (same fix needed for coho -- see coho_mc_probabilistic_retrospective.R)
  normalize_name <- function(x) tolower(gsub("[._]", "", x))
  bnames_stripped <- sub("^b_", "", bnames)
  
  draws <- tibble(draw = seq_len(n_draws), ra = post$intercept, rb_raw = post$b_spawners)
  
  for (cv in extra_covariates) {
    match_idx <- which(normalize_name(bnames_stripped) == normalize_name(cv))
    if (length(match_idx) != 1) {
      stop("Could not uniquely match covariate '", cv, "' to a monitored coefficient for stock ",
           stock_name, ". Monitored: ", paste(bnames, collapse = ", "))
    }
    draws[[cv]] <- post[[bnames[match_idx]]]
  }
  
  draws
}

stock_draws <- setNames(
  lapply(STOCKS_MC, build_stock_draws, n_draws = N_DRAWS),
  STOCKS_MC
)

# ------------------------------------------------------------
# 3. LOW-ABUNDANCE PERIOD PER STOCK -- fixed, from real observed RunJacks
#    (doesn't depend on draws, so computed once). Reuses GENERATION_TIME/
#    trailing_mean from your cleaned script.
# ------------------------------------------------------------

observed_returns <- obs %>%
  mutate(RunJacks = RunSize - JackEscapement) %>%
  filter(Stock %in% STOCKS_MC, is.finite(RunJacks)) %>%
  select(Stock, Year, RunJacks)

low_periods <- observed_returns %>%
  left_join(tibble::enframe(GENERATION_TIME, name = "Stock", value = "GT"), by = "Stock") %>%
  arrange(Stock, Year) %>%
  group_by(Stock) %>%
  mutate(obs_gen_mean = trailing_mean(RunJacks, width = first(GT))) %>%
  filter(!is.na(obs_gen_mean)) %>%
  slice_min(obs_gen_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Stock, GT, low_period_end = Year, low_period_start = Year - GT + 1,
            low_period_mean = obs_gen_mean)

# Historical max RunJacks per stock -- used to flag "runaway" draws below.
# Sockeye's recursion (retroR[i] <- retroS[j] * exp(retro_lnRS[j])) has no
# floor or ceiling anywhere, unlike steelhead (spawner floor at 0) or coho
# (never goes negative). A posterior draw with weak density dependence
# (small rb) combined with a large residual in a low-escapement year (where
# log(R/S) is inherently volatile) can compound exponentially over 60+
# years with nothing to check it. mean()/quantile() treat Inf/extreme
# values as contaminating (na.rm doesn't remove them), so these need to be
# flagged and excluded explicitly, not just averaged in -- see chat.
RUNAWAY_MULTIPLE <- 10   # a real population can't plausibly return at
# >10x its all-time historical high

historical_max <- observed_returns %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 4. MONTE CARLO LOOP -- pinniped / SST / pink, per stock, per draw
# ------------------------------------------------------------

SCENARIOS_MC <- list(
  "Pinniped scenario" = c("SeaLions", "seal"),
  "SST scenario"       = c("adult.sst", "smolt.sst"),
  "Pink scenario"      = "pink"
)

collapse_threshold <- 1

mc_results_summary <- vector("list", length(STOCKS_MC) * N_DRAWS)
counter <- 0

for (stock_name in STOCKS_MC) {
  
  stock_data <- obs %>% filter(Stock == stock_name)
  draws_table <- stock_draws[[stock_name]]
  extra_covariates <- STOCK_COVARIATES[[stock_name]]
  low_p <- low_periods %>% filter(Stock == stock_name)
  
  if (nrow(low_p) == 0) {
    warning("No low_period found for ", stock_name, " -- skipping")
    next
  }
  
  for (d in seq_len(N_DRAWS)) {
    
    draw_row <- draws_table[d, ]
    terms <- list(
      ra = draw_row$ra,
      rb = -draw_row$rb_raw,
      sel_covs = extra_covariates,
      cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
    )
    
    for (scenario_label in names(SCENARIOS_MC)) {
      
      scenario_vars <- intersect(SCENARIOS_MC[[scenario_label]], extra_covariates)
      # if this stock's top model doesn't use any of this scenario's
      # covariates, the scenario is identical to "actual" for this stock
      # -- still run it (matches your original script's behavior with
      # the no_freeze_year/no_alt_series/no_sst_mean warnings) rather
      # than silently skipping
      
      result <- run_retro_model_draw(
        stock_data, stock_name, terms, retroU_default, useretro_default, yrretro_default,
        scenario_vars = scenario_vars
      )
      
      scenario_period_mean <- result %>%
        filter(Year >= low_p$low_period_start, Year <= low_p$low_period_end) %>%
        summarise(m = mean(retroR, na.rm = TRUE)) %>%
        pull(m)
      
      counter <- counter + 1
      hist_max_this_stock <- historical_max$historical_max_return[historical_max$Stock == stock_name]
      mc_results_summary[[counter]] <- tibble(
        Stock                 = stock_name,
        scenario               = scenario_label,
        draw                    = d,
        low_period_mean_observed = low_p$low_period_mean,
        scenario_period_mean    = scenario_period_mean,
        abs_increase             = scenario_period_mean - low_p$low_period_mean,
        pct_increase             = 100 * (scenario_period_mean - low_p$low_period_mean) / low_p$low_period_mean,
        collapse                 = low_p$low_period_mean <= collapse_threshold,
        runaway                  = !is.finite(scenario_period_mean) |
          (scenario_period_mean / hist_max_this_stock) > RUNAWAY_MULTIPLE
      )
    }
  }
  
  message("done: ", stock_name)
}

mc_summary <- bind_rows(mc_results_summary)
saveRDS(mc_summary, "sockeye_mc_pct_increase_summary.rds")

# ------------------------------------------------------------
# 5. DIAGNOSTICS
# ------------------------------------------------------------

collapse_fraction <- mc_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_collapsed = sum(collapse, na.rm = TRUE),
            pct_collapsed = 100 * n_collapsed / n_draws, .groups = "drop")
print(collapse_fraction, n = Inf)

runaway_fraction <- mc_summary %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_runaway = sum(runaway, na.rm = TRUE),
            pct_runaway = 100 * n_runaway / n_draws, .groups = "drop") %>%
  arrange(desc(pct_runaway))
print(runaway_fraction, n = Inf)

# ------------------------------------------------------------
# 6. POSTERIOR SUMMARIES -- excludes runaway draws (see RUNAWAY_MULTIPLE
#    above); mean()/quantile() would otherwise be dominated or made Inf
#    by a small number of exploded draws
# ------------------------------------------------------------

mc_summary_stats_abs <- mc_summary %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    mean = mean(abs_increase, na.rm = TRUE), median = median(abs_increase, na.rm = TRUE),
    q2.5 = quantile(abs_increase, 0.025, na.rm = TRUE), q97.5 = quantile(abs_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_abs, n = Inf)

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

ggplot(mc_summary_stats_abs, aes(x = reorder(Stock, -mean), y = mean, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  labs(x = NULL, y = "Increase in mean abundance over low-abundance period (fish)",
       fill = "Scenario driver",
       title = "Sockeye abundance recovery over each stock's low-abundance period -- full posterior") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/sockeye_abs_increase_low_period_posterior.png", width = 12, height = 6, dpi = 600)

# % increase, reference only, excluding collapsed AND runaway draws
mc_summary_stats_pct <- mc_summary %>%
  filter(!collapse, !runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    mean = mean(pct_increase, na.rm = TRUE), median = median(pct_increase, na.rm = TRUE),
    q2.5 = quantile(pct_increase, 0.025, na.rm = TRUE), q97.5 = quantile(pct_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_pct, n = Inf)

ggplot(mc_summary_stats_pct, aes(x = reorder(Stock, -mean), y = mean, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  labs(x = NULL, y = "% increase in mean abundance over low-abundance period",
       fill = "Scenario driver",
       title = "Sockeye abundance recovery over each stock's low-abundance period -- full posterior") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/sockeye_pct_increase_low_period_posterior.png", width = 12, height = 6, dpi = 600)

# ============================================================
# Diagnostic: identify runaway-growth draws in the sockeye MC results,
# and trace one to confirm the mechanism (uncapped exponential Ricker
# recursion: retroR[i] <- retroS[j] * exp(retro_lnRS[j]), no floor/ceiling
# unlike steelhead's spawner floor or coho's non-negative structure).
#
# Run after mc_summary exists (from
# sockeye_mc_probabilistic_retrospective.R).
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. How many draws are "runaway" per stock/scenario -- defined as
#    scenario_period_mean far exceeding that stock's historical maximum
#    observed RunJacks (a real population can't plausibly return at
#    100x its all-time historical high)
# ------------------------------------------------------------

historical_max <- observed_returns %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

runaway_check <- mc_summary %>%
  left_join(historical_max, by = "Stock") %>%
  mutate(ratio_to_historical_max = scenario_period_mean / historical_max_return) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_draws = n(),
    n_runaway_10x  = sum(ratio_to_historical_max > 10, na.rm = TRUE),
    n_runaway_100x = sum(ratio_to_historical_max > 100, na.rm = TRUE),
    max_ratio      = max(ratio_to_historical_max, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_ratio))
print(runaway_check, n = Inf)

# ------------------------------------------------------------
# 2. Trace the single worst draw for the worst stock/scenario found above
# ------------------------------------------------------------

worst <- runaway_check %>% slice_max(max_ratio, n = 1)
worst_stock <- worst$Stock[1]
worst_scenario <- worst$scenario[1]

worst_draw_row <- mc_summary %>%
  left_join(historical_max, by = "Stock") %>%
  filter(Stock == worst_stock, scenario == worst_scenario) %>%
  mutate(ratio = scenario_period_mean / historical_max_return) %>%
  slice_max(ratio, n = 1)

cat("Worst case:", worst_stock, "/", worst_scenario,
    "/ draw", worst_draw_row$draw,
    "-- scenario_period_mean =", worst_draw_row$scenario_period_mean,
    "vs. historical max =", worst_draw_row$historical_max_return, "\n\n")

# Re-run that exact draw and print the full year-by-year trajectory
draw_num <- worst_draw_row$draw
draw_row <- stock_draws[[worst_stock]][draw_num, ]
extra_covariates <- STOCK_COVARIATES[[worst_stock]]

terms <- list(
  ra = draw_row$ra,
  rb = -draw_row$rb_raw,
  sel_covs = extra_covariates,
  cov_coefs = if (length(extra_covariates) > 0) as.numeric(draw_row[extra_covariates]) else numeric(0)
)
cat("Terms used:\n")
print(terms)

scenario_vars <- intersect(SCENARIOS_MC[[worst_scenario]], extra_covariates)

traj <- run_retro_model_draw(
  obs %>% filter(Stock == worst_stock), worst_stock, terms,
  retroU_default, useretro_default, yrretro_default,
  scenario_vars = scenario_vars
)

traj %>%
  select(Year, AdultEscapement, RunJacks, retroS, retroR, retro_lnRS, wt) %>%
  print(n = Inf)

# Year-over-year growth ratio in retroR -- look for where it takes off
traj %>%
  arrange(Year) %>%
  mutate(growth_ratio = retroR / lag(retroR)) %>%
  select(Year, retroR, growth_ratio) %>%
  filter(is.finite(growth_ratio)) %>%
  slice_max(growth_ratio, n = 10) %>%
  print()



# ============================================================
# Post-hoc fix for the sockeye MC results: mean()/quantile() treat Inf as
# contaminating, not something na.rm removes -- a handful of "runaway"
# draws (uncapped exponential Ricker recursion, see chat) were making
# every stock/scenario's posterior mean/CI effectively Inf or dominated
# by a few extreme values. This re-aggregates your EXISTING mc_summary
# (no need to rerun the MC loop) excluding those draws, and reports what
# fraction they were as its own diagnostic -- same pattern as `collapse`.
#
# Run after mc_summary and observed_returns already exist in your session
# (from sockeye_mc_probabilistic_retrospective.R).
#
# "Runaway" defined the same way as in sockeye_runaway_diagnostic.R:
# scenario_period_mean > 10x that stock's actual historical maximum
# RunJacks -- a real population cannot plausibly return at >10x its
# all-time historical high, so a draw exceeding that is a numerical
# artifact of the recursion, not a real biological possibility.
# ============================================================

library(tidyverse)

RUNAWAY_MULTIPLE <- 10   # adjust if you want a stricter/looser cutoff

historical_max <- observed_returns %>%
  group_by(Stock) %>%
  summarise(historical_max_return = max(RunJacks, na.rm = TRUE), .groups = "drop")

mc_summary_flagged <- mc_summary %>%
  left_join(historical_max, by = "Stock") %>%
  mutate(
    runaway = !is.finite(scenario_period_mean) |
      (scenario_period_mean / historical_max_return) > RUNAWAY_MULTIPLE
  )

# ------------------------------------------------------------
# DIAGNOSTIC: what fraction of draws were runaway, per stock/scenario --
# report this alongside any results, the same way collapse_fraction was
# reported for steelhead/coho
# ------------------------------------------------------------

runaway_fraction <- mc_summary_flagged %>%
  group_by(Stock, scenario) %>%
  summarise(n_draws = n(), n_runaway = sum(runaway, na.rm = TRUE),
            pct_runaway = 100 * n_runaway / n_draws, .groups = "drop") %>%
  arrange(desc(pct_runaway))
print(runaway_fraction, n = Inf)

# ------------------------------------------------------------
# RE-AGGREGATED POSTERIOR SUMMARIES, excluding runaway draws
# ------------------------------------------------------------

mc_summary_stats_abs_fixed <- mc_summary_flagged %>%
  filter(!runaway) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    mean   = mean(abs_increase, na.rm = TRUE),
    median = median(abs_increase, na.rm = TRUE),
    q2.5   = quantile(abs_increase, 0.025, na.rm = TRUE),
    q97.5  = quantile(abs_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_abs_fixed, n = Inf)

SCENARIO_COLORS2 <- c(
  "Pinniped scenario" = "#4682B4",
  "SST scenario"       = "#2E8B57",
  "Pink scenario"      = "#FF4500"
)

ggplot(mc_summary_stats_abs_fixed, aes(x = reorder(Stock, -mean), y = mean, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), position = position_dodge(width = 0.7), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = SCENARIO_COLORS2) +
  labs(x = NULL, y = "Increase in mean abundance over low-abundance period (fish)",
       fill = "Scenario driver",
       title = "Sockeye abundance recovery over each stock's low-abundance period",
       subtitle = paste0("Excludes runaway draws (>", RUNAWAY_MULTIPLE, "x historical max) -- see runaway_fraction")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ggsave("figures/sockeye_abs_increase_low_period_posterior_fixed.png", width = 12, height = 6, dpi = 600)

mc_summary_stats_pct_fixed <- mc_summary_flagged %>%
  filter(!runaway, !collapse) %>%
  group_by(Stock, scenario) %>%
  summarise(
    n_used = n(),
    mean   = mean(pct_increase, na.rm = TRUE),
    median = median(pct_increase, na.rm = TRUE),
    q2.5   = quantile(pct_increase, 0.025, na.rm = TRUE),
    q97.5  = quantile(pct_increase, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
print(mc_summary_stats_pct_fixed, n = Inf)


















