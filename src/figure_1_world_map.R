library(tidyverse)
library(rworldmap)

here::here("data","seq_information_usa_complete.csv") %>% 
  read_csv() -> us_coord

here::here("data","seq_information_africa_easia_complete.csv") %>% 
  read_csv() -> africa_easia_coord

here::here("data","india_locations.csv") %>% 
  read_csv() %>% 
  dplyr::select(region, y, x) %>% 
  dplyr::rename(lat = y,
                lon = x,
                country = region) -> india_coord



us_coord %>% 
  bind_rows(africa_easia_coord) %>%
  #write_csv("us_africa_easia_sequence_info.csv")
  dplyr::select(country, lat, lon) %>% 
  bind_rows(india_coord) %>% 
  drop_na() %>% 
  group_by(country, lat, lon) %>% 
  summarise(count = n()) -> count_locations


mapWorld <- borders("world", colour="black", fill="#f7e5c2")

# complete_codes %>% 
#   drop_na() %>% 
#   group_by(query, lat, lon) %>% 
#   summarise(count = n()) -> count_locations

count_locations %>% 
  ggplot() +
  mapWorld +
  labs(x = "",
       y = "",
       size = "Count") +
  coord_sf(ylim = c(-50, 90), datum = NA) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.key=element_blank(),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 13, face = "bold")
  ) +
  geom_point(aes(lon, lat, size = count), colour = "blue") +
  ggsave(here::here("output","sequence_map.png"), scale = 2.5, dpi = 300)
