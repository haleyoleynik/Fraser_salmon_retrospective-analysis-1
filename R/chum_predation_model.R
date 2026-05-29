# Fraser Chum consumption model 

# predation models 

chum <- read_csv("data/chum_s-r_data_pred.csv")


chum <- chum %>%
  mutate(chum_lnrs = log(chum_recruits/chum_spawners))


mod <- lm(chum_lnrs~chum_spawners + SSL, data = chum)
summary(mod)

coef(mod)
coef(mod)["chum_spawners"]

chum <- chum %>%
  mutate(Mdensity = -chum_spawners*coef(mod)["chum_spawners"],
         Mssl = -SSL*coef(mod)["SSL"],
         M = Mdensity + Mssl,
         Wt = chum_lnrs-(coef(mod)["(Intercept)"]-Mssl-Mdensity))

max_pred_n = -coef(mod)["SSL"] * max(chum$SSL,na.rm=T)
unfished_S = -coef(mod)["(Intercept)"]/coef(mod)["chum_spawners"]
density_effect = -unfished_S*coef(mod)["chum_spawners"]

Mp = chum %>%
  filter(Year >= 2012 & Year <= 2017) %>%
  summarise(mean = mean(Mssl)) %>%
  pull(mean)

# why is z set to 2?
Z = 2

avg_rec = chum %>%
  filter(Year >= 2012 & Year <= 2017) %>%
  summarise(mean = mean(chum_recruits)) %>%
  pull(mean)

Q = avg_rec*Mp/Z*(1-exp(-Z))
Q_SSL = Q/max(chum$SSL,na.rm=T) # over 3 years? 
NPAFC_catches = 0.06 # fraser averages 6% of coastwide totals (NPAFC)

chum_SSL_year = Q_SSL / NPAFC_catches /3 # 3 years
chum_weight = 5 # average weight of adult chum = 5 kg 
SSL_Q_year = 365*17

# percent of the diet expected if distributed over 3 years
chum_totQ = chum_SSL_year*chum_weight/SSL_Q_year



# plot 
chum %>%
  select(Year, M, Mdensity, Mssl) %>%
  pivot_longer(cols = M:Mssl, values_to = "Mortality", names_to = "Source") %>%
  mutate(Source = recode(Source, 
                         "M" = "Total",
                         "Mdensity" = "Density dependent",
                         "Mssl" = "Sea lion predation")) %>%
ggplot(aes(x=Year)) + 
  geom_line(aes(y=Mortality, color = Source), size = 0.8) +
  labs(y = "Mortality") + 
  theme_light()






               