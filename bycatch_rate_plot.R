# harvest rates for bycatch calc 
# figure 1 in chapter 3


by <- read_csv("data/harvest_rates_for_bycatch_calc.csv")

by %>% 
  pivot_longer(cols = 2:3, names_to = "stock", values_to = "harvest_rate") %>%
ggplot(aes(x = U_chum, y = harvest_rate, color = stock)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Observed ln(R/S)", y = "Predicted ln(R/S)") +
  scale_color_manual(values = palette.colors(n = 2, palette = "Okabe-Ito")) +
  theme_minimal() +
  theme(legend.position = "none")

# 1997 to 2017
by %>% 
  filter(Year > 1996) %>%
  filter(Year < 2017) %>%
  pivot_longer(cols = 2:3, names_to = "stock", values_to = "harvest_rate") %>%
  ggplot(aes(x = U_chum, y = harvest_rate, color = stock)) +
  geom_point(alpha = 0.6, size = 3, position = "jitter") +
  geom_smooth(method = "lm", se = F) +
  labs(x = "Chum salmon harvest rate", y = "Steelhead harvest rate") +
  scale_color_manual(values = palette.colors(n = 2, palette = "Okabe-Ito")) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave("figures/bycatch_rate_set.png", width = 6, height = 5, dpi = 600)


