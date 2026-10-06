# MSY heatmaps: % change in MSY across 2-D covariate grids (top model per stock)
# Recreates the "seal density x hatchery releases" style figure for Fraser sockeye,
# averaged across stocks.
#
# POSTERIOR-MEDIAN VERSION: parameters (alpha, beta, covariate coefficients) are the
# posterior medians from each stock's JAGS top model, plugged into the same
# deterministic calculation as the frequentist version.
#
# requires: all.sockeye.st loaded (standardized covariates, by Stock and yr),
#           STOCK_COVARIATES and jags_name() from sockeye_jags_models_b-bound_alpha05_v2.R, and either
#           sockeye_fits in the session or the saved *_posterior_samples.rds files
library(tidyverse)
library(coda)
library(patchwork)

# settings -----------------------------------------------------------------
n_grid <- 60   # grid resolution per axis
include_late_shuswap <- TRUE   # cycle-structured model; each cycle enters as its own unit

# how to combine stock MSYs into one surface
# recruit-weighted mean (as in the partial-effect code); swap in the sum for Fraser-wide total MSY
agg_msy <- function(msy, w) weighted.mean(msy, w, na.rm = TRUE)
# agg_msy <- function(msy, w) sum(msy, na.rm = TRUE)

# composite axes: covariates in the same group move together across their observed range
axis_groups <- list(
  pinnipeds = c("seal", "SeaLions"),
  pink      = "pink",
  adult_sst = "adult.sst",
  juv_sst   = "smolt.sst"
)
axis_labels <- c(
  pinnipeds = "Pinniped abundance (seals + sea lions)",
  pink      = "Pink salmon abundance",
  adult_sst = "Adult SST",
  juv_sst   = "Juvenile (smolt) SST"
)
plot_pairs <- list(c("pinnipeds", "adult_sst"),
                   c("pinnipeds", "juv_sst"),
                   c("pinnipeds", "pink"),
                   c("pink",      "juv_sst"),
                   c("pink",      "adult_sst"))

all_covs <- c("seal", "SeaLions", "pink", "smolt.sst", "adult.sst", "PDO", "NPGO")

# top model per stock: posterior medians from the JAGS fits ----------------
get_post_median <- function(stock) {
  samp <- if (exists("sockeye_fits")) {
    sockeye_fits[[stock]]$samples
  } else {
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock)), "_posterior_samples.rds"))
  }
  m    <- as.matrix(samp)                                   # chains combined
  meds <- apply(m, 2, median)
  
  # b_<jags name> -> original covariate column name (e.g. b_adult_sst -> adult.sst)
  lookup <- setNames(all_covs, paste0("b_", jags_name(all_covs)))
  covs   <- setNames(rep(0, length(all_covs)), all_covs)   # not in model -> 0
  b_cov  <- intersect(names(meds), names(lookup))
  covs[lookup[b_cov]] <- meds[b_cov]
  
  if (stock == "Late Shuswap") {
    # cycle-structured model: one alpha and one b per cycle (cycle = BroodYear %% 4),
    # shared covariate coefficient(s). Each cycle becomes its own "unit".
    return(tibble(Stock = stock,
                  unit  = paste0(stock, " c", 0:3),
                  cycle = 0:3,
                  alpha = meds[paste0("ra_cycle[", 1:4, "]")],
                  beta  = meds[paste0("b_cycle[", 1:4, "]")],
                  !!!as.list(covs)))
  }
  
  # bounded-b model: lnrs = intercept - b * spawners, b > 0 -> beta = b (no sign flip)
  tibble(Stock = stock, unit = stock, cycle = NA_integer_,
         alpha = meds[["intercept"]],
         beta  = meds[["b"]],
         !!!as.list(covs))
}

stocks_to_plot <- names(STOCK_COVARIATES)   # the 15 single-alpha stocks
if (include_late_shuswap) stocks_to_plot <- c(stocks_to_plot, "Late Shuswap")

top_models <- map_dfr(stocks_to_plot, get_post_median)

if (any(top_models$beta <= 0)) warning("Some stocks have a posterior-median beta <= 0; Smsy is undefined for those.")

# optional: compare with the frequentist top-model coefficients
# read_csv("Results/Fraser_sockeye_dredge-results_alph05.csv") %>%
#   group_by(Stock) %>% slice_min(AIC, n = 1, with_ties = FALSE) %>% ungroup() %>%
#   select(Stock, `(Intercept)`, spawners, any_of(all_covs)) %>%
#   left_join(top_models, by = "Stock", suffix = c("_freq", "_bayes"))

# stock data, covariate ranges, recruit weights ----------------------------
all.sockeye.st <- read_csv("data/sockeye_data.csv")

stock_dat <- all.sockeye.st %>%
  filter(Stock %in% top_models$Stock)

missing_stocks <- setdiff(top_models$Stock, stock_dat$Stock)
if (length(missing_stocks) > 0) {
  warning("Not in all.sockeye.st, dropped: ", paste(missing_stocks, collapse = ", "))
  top_models <- filter(top_models, Stock %in% stock_dat$Stock)
}

# each stock's observed range of each standardized covariate (0-1 axis maps onto this)
cov_ranges <- stock_dat %>%
  select(Stock, all_of(all_covs)) %>%
  pivot_longer(-Stock, names_to = "covariate", values_to = "z") %>%
  group_by(Stock, covariate) %>%
  summarize(lo = min(z, na.rm = TRUE), hi = max(z, na.rm = TRUE), .groups = "drop")

# weights per unit: mean recruits for ordinary stocks; for Late Shuswap each cycle gets
# (mean recruits in that cycle's brood years) / 4, so the dominant cycle counts most and
# the four cycles together weigh about as much as the stock's overall mean recruits
stock_weights <- top_models %>%
  distinct(Stock, unit, cycle) %>%
  mutate(w = pmap_dbl(list(Stock, cycle), function(s, cy) {
    d <- filter(stock_dat, Stock == s)
    if (is.na(cy)) mean(d$recruits, na.rm = TRUE)
    else mean(d$recruits[d$yr %% 4 == cy], na.rm = TRUE) / 4
  })) %>%
  select(unit, w)

# MSY (Hilborn 1985 approximation); no surplus production when alpha_t <= 0
calc_msy <- function(alpha_t, beta) {
  a    <- pmax(alpha_t, 0)
  Umsy <- a * (0.5 - 0.07 * a)
  Smsy <- Umsy / beta
  Rmsy <- Smsy * exp(a - beta * Smsy)
  Rmsy - Smsy
}

# build surface + first/last-year points for one pair of axes --------------
# covariates not on either axis are held at their mean (0)
make_surface <- function(x_axis, y_axis, n = n_grid) {
  xv   <- axis_groups[[x_axis]]
  yv   <- axis_groups[[y_axis]]
  vars <- c(xv, yv)
  
  coefs <- top_models %>%
    select(Stock, unit, alpha, beta, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "covariate", values_to = "coef") %>%
    left_join(cov_ranges, by = c("Stock", "covariate")) %>%
    mutate(axis = if_else(covariate %in% xv, "x", "y"))
  
  # grid: position p in [0,1] -> each stock's own min..max of that covariate
  grid <- expand_grid(px = seq(0, 1, length.out = n),
                      py = seq(0, 1, length.out = n))
  
  surf_stock <- coefs %>%
    crossing(grid) %>%
    mutate(p   = if_else(axis == "x", px, py),
           eff = if_else(coef == 0, 0, coef * (lo + p * (hi - lo)))) %>%
    group_by(Stock, unit, px, py) %>%
    summarize(alpha_t = first(alpha) + sum(eff),
              beta    = first(beta), .groups = "drop") %>%
    mutate(MSY = calc_msy(alpha_t, beta))
  
  # observed years: first and last year with data for every stock
  obs <- stock_dat %>%
    select(Stock, yr, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "covariate", values_to = "z") %>%
    group_by(Stock, yr) %>%
    filter(!any(is.na(z))) %>%
    ungroup()
  
  common_yrs <- obs %>%
    distinct(Stock, yr) %>%
    count(yr) %>%
    filter(n == n_distinct(top_models$Stock)) %>%
    pull(yr)
  if (length(common_yrs) == 0) stop("No year has ", paste(vars, collapse = "/"), " data for every stock.")
  ref_yrs <- range(common_yrs)
  
  pts_stock <- obs %>%
    filter(yr %in% ref_yrs) %>%
    left_join(coefs %>% select(Stock, unit, covariate, alpha, beta, coef, lo, hi, axis),
              by = c("Stock", "covariate"), relationship = "many-to-many") %>%
    mutate(p   = (z - lo) / (hi - lo),
           eff = coef * z) %>%
    group_by(Stock, unit, yr) %>%
    summarize(px      = mean(p[axis == "x"]),
              py      = mean(p[axis == "y"]),
              alpha_t = first(alpha) + sum(eff),
              beta    = first(beta), .groups = "drop") %>%
    mutate(MSY = calc_msy(alpha_t, beta))
  
  # aggregate across stocks (Late Shuswap: its four cycles enter as separate units)
  surf <- surf_stock %>%
    left_join(stock_weights, by = "unit") %>%
    group_by(px, py) %>%
    summarize(MSY = agg_msy(MSY, w), .groups = "drop")
  
  pts <- pts_stock %>%
    left_join(stock_weights, by = "unit") %>%
    group_by(yr) %>%
    summarize(px  = weighted.mean(px, w, na.rm = TRUE),
              py  = weighted.mean(py, w, na.rm = TRUE),
              MSY = agg_msy(MSY, w), .groups = "drop")
  
  # % change relative to MSY in the first common year
  base     <- pts$MSY[pts$yr == ref_yrs[1]]
  surf$pct <- 100 * (surf$MSY / base - 1)
  pts$pct  <- 100 * (pts$MSY  / base - 1)
  
  list(surface = surf, points = pts, per_stock = surf_stock,
       ref_years = ref_yrs, x_axis = x_axis, y_axis = y_axis)
}

surfaces <- map(plot_pairs, ~ make_surface(.x[1], .x[2]))

# shared colour scale ------------------------------------------------------
all_pct <- unlist(map(surfaces, ~ c(.x$surface$pct, .x$points$pct)))
lims <- c(min(floor(min(all_pct, na.rm = TRUE) / 10) * 10, -10),
          max(ceiling(max(all_pct, na.rm = TRUE) / 10) * 10,  10))

msy_fill <- scale_fill_stepsn(
  colours = c("#D7191C", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#FFFFFF"),
  values  = scales::rescale(c(lims[1], 0.6 * lims[1], 0.3 * lims[1], 0.1 * lims[1], 0, lims[2])),
  limits  = lims,
  breaks  = seq(lims[1], lims[2], by = 10),
  name    = "% change\nin MSY\n(posterior\nmedians)"
)

# plot ---------------------------------------------------------------------
plot_msy_surface <- function(s) {
  ggplot(s$surface, aes(px, py)) +
    geom_raster(aes(fill = pct)) +
    geom_point(data = s$points, size = 2) +
    geom_text(data = s$points,
              aes(label = sprintf("%d\n(%+.0f%%)", as.integer(yr), pct)),
              vjust = -0.3, size = 3, lineheight = 0.9) +
    msy_fill +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    labs(x = axis_labels[[s$x_axis]], y = axis_labels[[s$y_axis]]) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())
}

msy_plots <- map(surfaces, plot_msy_surface)

(msy_fig <- wrap_plots(msy_plots, nrow = 2) + plot_layout(guides = "collect"))

ggsave("figures/sockeye_msy_heatmaps_posterior-median_5panel.png", msy_fig,
       width = 15, height = 10, dpi = 300)

##################################
# Three panel 

# MSY heatmaps: % change in MSY across 2-D covariate grids (top model per stock)
# Recreates the "seal density x hatchery releases" style figure for Fraser sockeye,
# averaged across stocks.
#
# POSTERIOR-MEDIAN VERSION: parameters (alpha, beta, covariate coefficients) are the
# posterior medians from each stock's JAGS top model, plugged into the same
# deterministic calculation as the frequentist version.
#
# requires: all.sockeye.st loaded (standardized covariates, by Stock and yr),
#           STOCK_COVARIATES and jags_name() from sockeye_jags_models_b-bound_alpha05_v2.R, and either
#           sockeye_fits in the session or the saved *_posterior_samples.rds files
library(tidyverse)
library(coda)
library(patchwork)

# settings -----------------------------------------------------------------
n_grid <- 60   # grid resolution per axis
include_late_shuswap <- FALSE   # cycle-structured model; each cycle enters as its own unit

# how to combine stock MSYs into one surface
# recruit-weighted mean (as in the partial-effect code); swap in the sum for Fraser-wide total MSY
agg_msy <- function(msy, w) weighted.mean(msy, w, na.rm = TRUE)
# agg_msy <- function(msy, w) sum(msy, na.rm = TRUE)

# composite axes: covariates in the same group move together across their observed range
axis_groups <- list(
  pinnipeds = c("seal", "SeaLions"),
  pink      = "pink",
  sst       = c("smolt.sst", "adult.sst")
)
axis_labels <- c(
  pinnipeds = "Pinniped abundance (seals + sea lions)",
  pink      = "Pink salmon abundance",
  sst       = "SST (juvenile + adult)"
)
plot_pairs <- list(c("pinnipeds", "sst"),
                   c("pinnipeds", "pink"),
                   c("pink",      "sst"))

all_covs <- c("seal", "SeaLions", "pink", "smolt.sst", "adult.sst", "PDO", "NPGO")

# top model per stock: posterior medians from the JAGS fits ----------------
get_post_median <- function(stock) {
  samp <- if (exists("sockeye_fits")) {
    sockeye_fits[[stock]]$samples
  } else {
    readRDS(paste0("sockeye_", tolower(gsub(" ", "_", stock)), "_posterior_samples.rds"))
  }
  m    <- as.matrix(samp)                                   # chains combined
  meds <- apply(m, 2, median)
  
  # b_<jags name> -> original covariate column name (e.g. b_adult_sst -> adult.sst)
  lookup <- setNames(all_covs, paste0("b_", jags_name(all_covs)))
  covs   <- setNames(rep(0, length(all_covs)), all_covs)   # not in model -> 0
  b_cov  <- intersect(names(meds), names(lookup))
  covs[lookup[b_cov]] <- meds[b_cov]
  
  if (stock == "Late Shuswap") {
    # cycle-structured model: one alpha and one b per cycle (cycle = BroodYear %% 4),
    # shared covariate coefficient(s). Each cycle becomes its own "unit".
    return(tibble(Stock = stock,
                  unit  = paste0(stock, " c", 0:3),
                  cycle = 0:3,
                  alpha = meds[paste0("ra_cycle[", 1:4, "]")],
                  beta  = meds[paste0("b_cycle[", 1:4, "]")],
                  !!!as.list(covs)))
  }
  
  # bounded-b model: lnrs = intercept - b * spawners, b > 0 -> beta = b (no sign flip)
  tibble(Stock = stock, unit = stock, cycle = NA_integer_,
         alpha = meds[["intercept"]],
         beta  = meds[["b"]],
         !!!as.list(covs))
}

stocks_to_plot <- names(STOCK_COVARIATES)   # the 15 single-alpha stocks
if (include_late_shuswap) stocks_to_plot <- c(stocks_to_plot, "Late Shuswap")

top_models <- map_dfr(stocks_to_plot, get_post_median)

if (any(top_models$beta <= 0)) warning("Some stocks have a posterior-median beta <= 0; Smsy is undefined for those.")

# optional: compare with the frequentist top-model coefficients
# read_csv("Results/Fraser_sockeye_dredge-results_alph05.csv") %>%
#   group_by(Stock) %>% slice_min(AIC, n = 1, with_ties = FALSE) %>% ungroup() %>%
#   select(Stock, `(Intercept)`, spawners, any_of(all_covs)) %>%
#   left_join(top_models, by = "Stock", suffix = c("_freq", "_bayes"))

# stock data, covariate ranges, recruit weights ----------------------------
stock_dat <- all.sockeye.st %>%
  filter(Stock %in% top_models$Stock)

missing_stocks <- setdiff(top_models$Stock, stock_dat$Stock)
if (length(missing_stocks) > 0) {
  warning("Not in all.sockeye.st, dropped: ", paste(missing_stocks, collapse = ", "))
  top_models <- filter(top_models, Stock %in% stock_dat$Stock)
}

# each stock's observed range of each standardized covariate (0-1 axis maps onto this)
cov_ranges <- stock_dat %>%
  select(Stock, all_of(all_covs)) %>%
  pivot_longer(-Stock, names_to = "covariate", values_to = "z") %>%
  group_by(Stock, covariate) %>%
  summarize(lo = min(z, na.rm = TRUE), hi = max(z, na.rm = TRUE), .groups = "drop")

# weights per unit: mean recruits for ordinary stocks; for Late Shuswap each cycle gets
# (mean recruits in that cycle's brood years) / 4, so the dominant cycle counts most and
# the four cycles together weigh about as much as the stock's overall mean recruits
stock_weights <- top_models %>%
  distinct(Stock, unit, cycle) %>%
  mutate(w = pmap_dbl(list(Stock, cycle), function(s, cy) {
    d <- filter(stock_dat, Stock == s)
    if (is.na(cy)) mean(d$recruits, na.rm = TRUE)
    else mean(d$recruits[d$yr %% 4 == cy], na.rm = TRUE) / 4
  })) %>%
  select(unit, w)

# MSY (Hilborn 1985 approximation); no surplus production when alpha_t <= 0
calc_msy <- function(alpha_t, beta) {
  a    <- pmax(alpha_t, 0)
  Umsy <- a * (0.5 - 0.07 * a)
  Smsy <- Umsy / beta
  Rmsy <- Smsy * exp(a - beta * Smsy)
  Rmsy - Smsy
}

# build surface + first/last-year points for one pair of axes --------------
# covariates not on either axis are held at their mean (0)
make_surface <- function(x_axis, y_axis, n = n_grid) {
  xv   <- axis_groups[[x_axis]]
  yv   <- axis_groups[[y_axis]]
  vars <- c(xv, yv)
  
  coefs <- top_models %>%
    select(Stock, unit, alpha, beta, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "covariate", values_to = "coef") %>%
    left_join(cov_ranges, by = c("Stock", "covariate")) %>%
    mutate(axis = if_else(covariate %in% xv, "x", "y"))
  
  # grid: position p in [0,1] -> each stock's own min..max of that covariate
  grid <- expand_grid(px = seq(0, 1, length.out = n),
                      py = seq(0, 1, length.out = n))
  
  surf_stock <- coefs %>%
    crossing(grid) %>%
    mutate(p   = if_else(axis == "x", px, py),
           eff = if_else(coef == 0, 0, coef * (lo + p * (hi - lo)))) %>%
    group_by(Stock, unit, px, py) %>%
    summarize(alpha_t = first(alpha) + sum(eff),
              beta    = first(beta), .groups = "drop") %>%
    mutate(MSY = calc_msy(alpha_t, beta))
  
  # observed years: first and last year with data for every stock
  obs <- stock_dat %>%
    select(Stock, yr, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "covariate", values_to = "z") %>%
    group_by(Stock, yr) %>%
    filter(!any(is.na(z))) %>%
    ungroup()
  
  common_yrs <- obs %>%
    distinct(Stock, yr) %>%
    count(yr) %>%
    filter(n == n_distinct(top_models$Stock)) %>%
    pull(yr)
  if (length(common_yrs) == 0) stop("No year has ", paste(vars, collapse = "/"), " data for every stock.")
  ref_yrs <- range(common_yrs)
  
  pts_stock <- obs %>%
    filter(yr %in% ref_yrs) %>%
    left_join(coefs %>% select(Stock, unit, covariate, alpha, beta, coef, lo, hi, axis),
              by = c("Stock", "covariate"), relationship = "many-to-many") %>%
    mutate(p   = (z - lo) / (hi - lo),
           eff = coef * z) %>%
    group_by(Stock, unit, yr) %>%
    summarize(px      = mean(p[axis == "x"]),
              py      = mean(p[axis == "y"]),
              alpha_t = first(alpha) + sum(eff),
              beta    = first(beta), .groups = "drop") %>%
    mutate(MSY = calc_msy(alpha_t, beta))
  
  # aggregate across stocks (Late Shuswap: its four cycles enter as separate units)
  surf <- surf_stock %>%
    left_join(stock_weights, by = "unit") %>%
    group_by(px, py) %>%
    summarize(MSY = agg_msy(MSY, w), .groups = "drop")
  
  pts <- pts_stock %>%
    left_join(stock_weights, by = "unit") %>%
    group_by(yr) %>%
    summarize(px  = weighted.mean(px, w, na.rm = TRUE),
              py  = weighted.mean(py, w, na.rm = TRUE),
              MSY = agg_msy(MSY, w), .groups = "drop")
  
  # % change relative to MSY in the first common year
  base     <- pts$MSY[pts$yr == ref_yrs[1]]
  surf$pct <- 100 * (surf$MSY / base - 1)
  pts$pct  <- 100 * (pts$MSY  / base - 1)
  
  list(surface = surf, points = pts, per_stock = surf_stock,
       ref_years = ref_yrs, x_axis = x_axis, y_axis = y_axis)
}

surfaces <- map(plot_pairs, ~ make_surface(.x[1], .x[2]))

# shared colour scale ------------------------------------------------------
all_pct <- unlist(map(surfaces, ~ c(.x$surface$pct, .x$points$pct)))
lims <- c(min(floor(min(all_pct, na.rm = TRUE) / 10) * 10, -10),
          max(ceiling(max(all_pct, na.rm = TRUE) / 10) * 10,  10))

msy_fill <- scale_fill_stepsn(
  colours = c("#D7191C", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#FFFFFF"),
  values  = scales::rescale(c(lims[1], 0.6 * lims[1], 0.3 * lims[1], 0.1 * lims[1], 0, lims[2])),
  limits  = lims,
  breaks  = seq(lims[1], lims[2], by = 10),
  name    = "% change\nin MSY\n(posterior\nmedians)"
)

# plot ---------------------------------------------------------------------
plot_msy_surface <- function(s) {
  ggplot(s$surface, aes(px, py)) +
    geom_raster(aes(fill = pct)) +
    geom_point(data = s$points, size = 2) +
    geom_text(data = s$points,
              aes(label = sprintf("%d\n(%+.0f%%)", as.integer(yr), pct)),
              vjust = -0.3, size = 3, lineheight = 0.9) +
    msy_fill +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    labs(x = axis_labels[[s$x_axis]], y = axis_labels[[s$y_axis]]) +
    theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank())
}

msy_plots <- map(surfaces, plot_msy_surface)

(msy_fig <- wrap_plots(msy_plots, nrow = 1) + plot_layout(guides = "collect"))

ggsave("Figures/sockeye_msy_heatmaps_posterior-median_3panel_wo_shuswap.png", msy_fig,
       width = 15, height = 5.5, dpi = 300)

#
# share of the first-year cross-stock MSY coming from each stock/cycle
s <- surfaces[[1]]   # any panel; baseline is the first common year
yr0 <- s$ref_years[1]

top_models %>%
  left_join(stock_weights, by = "unit") %>%
  rowwise() %>%
  mutate(MSY_yr0 = {
    d <- filter(stock_dat, Stock == .data$Stock, yr == yr0)
    a <- alpha + sum(c_across(all_of(all_covs)) * unlist(d[1, all_covs]), na.rm = TRUE)
    calc_msy(a, beta)
  }) %>%
  ungroup() %>%
  mutate(share = w * MSY_yr0 / sum(w * MSY_yr0, na.rm = TRUE)) %>%
  select(unit, w, alpha, beta, MSY_yr0, share) %>%
  arrange(desc(share))

