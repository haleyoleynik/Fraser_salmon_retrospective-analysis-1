# ============================================================
# Diagnostic: single deterministic run with POSTERIOR MEAN coefficients
# (no Monte Carlo, no loop) -- to see exactly where and why the "No SSL
# control" reconstruction collapses to (near) zero abundance.
#
# Run this after everything through section 5 of
# steelhead_mc_probabilistic_retrospective.R has executed (needs
# all_data_no_control, run_thompson_scenario_v2, run_chilcotin_scenario_v2,
# thompson_post, chilcotin_post, default_byrate, MODEL_YEARS).
# ============================================================

th_mean <- colMeans(thompson_post)
ch_mean <- colMeans(chilcotin_post)

print(th_mean)   # sanity check: do a, b, s, t look like reasonable-magnitude numbers?
print(ch_mean)   # same for a, b, s, t, f

diag_thompson <- run_thompson_scenario_v2(
  all_data_no_control, SSL_control = 0, U_historic = 1, byrate = default_byrate,
  sh_thompson_intercept     = th_mean["a"],
  sh_thompson_sst_coef      = th_mean["t"],
  sh_thompson_ssl_coef      = th_mean["s"],
  sh_thompson_spawners_coef = -th_mean["b"]
)

diag_chilcotin <- run_chilcotin_scenario_v2(
  diag_thompson, SSL_control = 0, U_historic = 1, byrate = default_byrate,
  sh_chilcotin_intercept      = ch_mean["a"],
  sh_chilcotin_sst_coef       = ch_mean["t"],
  sh_chilcotin_ssl_coef       = ch_mean["s"],
  sh_chilcotin_maxflow_coef   = ch_mean["f"],
  sh_chilcotin_spawners_coef  = -ch_mean["b"]
)

# Thompson trace: actual spawners vs. reconstructed spawners_pred, total
# return, recruits_alt, the productivity term, and whether/where the
# observed-residual calibration was available (NA = no calibration that
# year, i.e. resid_adj fell back to 1)
diag_thompson %>%
  filter(Year %in% MODEL_YEARS) %>%
  transmute(
    Year,
    spawners_obs   = sh_thompson_spawners,
    spawners_pred  = sh_thompson_spawners_pred,
    sum_pred       = sh_thompson_sum_pred,
    recruits_alt   = sh_thompson_recruits_alt,
    alpha_CN       = sh_thompson_alpha_CN,
    ln_obs_pred    = sh_thompson_ln_obs_pred,     # NA => resid_adj fell back to 1 that year
    U_comm         = sh_thompson_U_comm,
    total_catch    = sh_thompson_total_catch_pred
  ) %>%
  print(n = 100)

cat("\n--- Chilcotin ---\n\n")

diag_chilcotin %>%
  filter(Year %in% MODEL_YEARS) %>%
  transmute(
    Year,
    spawners_obs   = sh_chilcotin_spawners,
    spawners_pred  = sh_chilcotin_spawners_pred,
    sum_pred       = sh_chilcotin_sum_pred,
    recruits_alt   = sh_chilcotin_recruits_alt,
    alpha_CN       = sh_chilcotin_alpha_CN,
    ln_obs_pred    = sh_chilcotin_ln_obs_pred,
    U_comm         = sh_chilcotin_U_comm,
    total_catch    = sh_chilcotin_total_catch_pred
  ) %>%
  print(n = 100)

# Quick flag: first year spawners_pred hits (near) zero, for each stock
cat("\nFirst near-zero spawners_pred year, Thompson:\n")
diag_thompson %>% filter(Year %in% MODEL_YEARS, sh_thompson_spawners_pred <= 1) %>% slice(1) %>% print()
diag_thompson %>%
  filter(Year %in% MODEL_YEARS, sh_thompson_spawners_pred <= 1) %>%
  slice(1) %>%
  select(Year, sh_thompson_spawners, sh_thompson_spawners_pred, sh_thompson_sum_pred,
         sh_thompson_recruits_alt, sh_thompson_alpha_CN, sh_thompson_ln_obs_pred,
         sh_thompson_U_comm, sh_thompson_total_catch_pred) %>%
  print()

cat("\nFirst near-zero spawners_pred year, Chilcotin:\n")
diag_chilcotin %>% filter(Year %in% MODEL_YEARS, sh_chilcotin_spawners_pred <= 1) %>% slice(1) %>% print()
diag_chilcotin %>%
  filter(Year >= 1980, Year <= 1995) %>%
  select(Year, sh_chilcotin_spawners, sh_chilcotin_spawners_pred, sh_chilcotin_sum_pred,
         sh_chilcotin_recruits_alt, sh_chilcotin_alpha_CN, sh_chilcotin_ln_obs_pred,
         sh_chilcotin_U_comm, sh_chilcotin_U, sh_chilcotin_total_catch_pred) %>%
  print(n = 20)

print(th_mean)
print(ch_mean)

diag_chilcotin %>% 
  filter(Year > 1984) %>%
  select(Year, sh_chilcotin_U_comm, sh_chilcotin_U, sh_chilcotin_total_catch_pred)%>%
  print(n = 20)


# ============================================================
# Diagnostic: shape of the abs_increase distribution, per stock
#
# Run this after mc_summary exists (from
# steelhead_mc_probabilistic_retrospective.R). Shows whether the wide
# credible interval is one smooth continuous spread (genuine compounding
# uncertainty) or being dragged out by a small number of extreme outlier
# draws (which would point more toward small-sample noise or a specific
# problematic draw worth inspecting individually).
# ============================================================

library(tidyverse)

# Detailed quantiles -- finer-grained than just the 95% CI endpoints
abs_increase_quantiles <- mc_summary %>%
  group_by(Stock) %>%
  summarise(
    n      = n(),
    min    = min(abs_increase, na.rm = TRUE),
    p10    = quantile(abs_increase, 0.10, na.rm = TRUE),
    p25    = quantile(abs_increase, 0.25, na.rm = TRUE),
    median = median(abs_increase, na.rm = TRUE),
    p75    = quantile(abs_increase, 0.75, na.rm = TRUE),
    p90    = quantile(abs_increase, 0.90, na.rm = TRUE),
    max    = max(abs_increase, na.rm = TRUE),
    .groups = "drop"
  )
print(abs_increase_quantiles)

# Which specific draws are the extremes -- worth a look if just 1-2 draws
# are far outside the rest, since that points to a specific coefficient
# combination rather than a smooth spread
mc_summary %>%
  group_by(Stock) %>%
  slice_max(abs_increase, n = 3) %>%
  select(Stock, draw, abs_increase, `No SSL control`, `SSL control`) %>%
  print()

mc_summary %>%
  group_by(Stock) %>%
  slice_min(abs_increase, n = 3) %>%
  select(Stock, draw, abs_increase, `No SSL control`, `SSL control`) %>%
  print()

# Histogram: is it one smooth spread, or a cluster + a few outliers?
ggplot(mc_summary, aes(abs_increase)) +
  geom_histogram(bins = 20, fill = "#4682B4", color = "white") +
  facet_wrap(~ Stock, scales = "free_x") +
  labs(x = "Increase in mean abundance over low-abundance period (fish)",
       y = "Number of draws",
       title = "Distribution of abs_increase across posterior draws") +
  theme_minimal()

# Same thing but on a log scale for the x-axis, which often makes the
# shape clearer when a Ricker-type model produces a right-skewed spread
# (a few very large values stretching out an otherwise tight cluster)
ggplot(mc_summary %>% filter(abs_increase > 0), aes(abs_increase)) +
  geom_histogram(bins = 20, fill = "#4682B4", color = "white") +
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~ Stock, scales = "free_x") +
  labs(x = "Increase in mean abundance over low-abundance period (fish, log scale)",
       y = "Number of draws",
       title = "Distribution of abs_increase across posterior draws (log x-axis)") +
  theme_minimal()

















