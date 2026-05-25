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
Q_SSL = Q/max(chum$SSL,na.rm=T)









               