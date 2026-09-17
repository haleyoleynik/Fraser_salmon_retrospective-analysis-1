require(tidyverse)
require(ggplot2)
require(readr)

df <- read_csv("data/price_per_fish_2020-2023.csv")

# price per fish 
df %>% group_by(species) %>%
  summarise(landed_value = sum(landed_value, na.rm=T),
            number = sum(number, na.rm=T),
            price = landed_value/number)


