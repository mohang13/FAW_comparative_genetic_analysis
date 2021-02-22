library(tidyverse)
library(rgdal)
library(ggrepel)
library(patchwork)
library(Biostrings)
library(ggmsa)
library(janitor)


readOGR(here::here("data","India State wise shape files/state.shp")) -> india_shp

broom::tidy(india_shp) -> tidy_india_shp


here::here("data","india_locations.csv") %>% 
  read_csv() %>% 
  dplyr::select(region, haplotype_number, y, x) %>% 
  dplyr::rename(lat = y,
                lon = x) %>% 
  mutate(haplotype_number = haplotype_number %>% factor()) %>% 
  ggplot() +
  geom_path(data = tidy_india_shp, 
            aes(x = long, y = lat, group = group),
            color = 'black', fill = 'white', size = .2) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.key=element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.spacing.y = unit(0.5, 'cm')
  ) +
  geom_point(aes(lon, lat,
                 colour = haplotype_number
  )) +
  scale_colour_manual(values = c("#b22222", "#ff00ff", "#696969", "#ffa07a", "#808000",  "#008000", "#00008b", "#ff0000", "#ff8c00", "#ffd700", "#00ff7f",
                                 "#4169e1", "#00bfff",  "#db7093", "#f0e68c", "#ff1493", "#ffa07a", "#ee82ee", "#adff2f", "#0000ff", "#7fffd4")) +
  labs(x = "",
       y = "",
       colour = "Haplotypes") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size=3))) +
  geom_text_repel(aes(lon, lat, label = haplotype_number), 
                  segment.color = 'grey80', fontface = "bold", segment.alpha = 0.75,
                  point.padding = 0.05, size = 3, max.overlaps = 75) +
  theme(text=element_text(family="sans-serif")) -> india_plot


# haplotype count India segment plot----
here::here("data","india_locations.csv") %>% 
  read_csv() %>% 
  select(haplotype_number) %>% 
  count(haplotype_number) %>% 
  mutate(haplotype_number = haplotype_number %>% factor() %>% fct_rev()) %>% 
  ggplot(aes(haplotype_number,n)) +
  geom_segment(aes(xend = haplotype_number, yend = 0)) +
  geom_point(size=2) +
  geom_text(aes(label = n), nudge_y = 7.5, size = 3, fontface = "bold") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0), limits = c(0,160), breaks = seq(0,160,15)) +
  theme_bw() +
  labs(x = "Haplotypes", 
       y = "Frequency") +
  theme(text=element_text(family="sans-serif"),
        panel.grid = element_blank()) -> count_segment_plot

#coi rs plot

tibble::tribble(
  ~V1,     ~V2,     ~V3,
  0L, 0.75028, 0.77018,
  1L, 0.20569, 0.17700,
  2L, 0.03959, 0.04068,
  3L, 0.00423, 0.00935,
  4L, 0.00021, 0.00215,
  5L, 0.00000, 0.00049,
  6L, 0.00000, 0.00011,
  7L,       0, 0.00003,
  8L,       0, 0.00001,
  9L,       0,       0,
  10L,       0,       0,
  11L,       0,       0,
  12L,       0,       0,
  13L,       0,       0,
  14L,       0,       0,
  15L,       0,       0,
  16L,       0,       0,
  17L,       0,       0,
  18L,       0,       0,
  19L,       0,       0,
  20L,       0,       0
) -> strain_data

strain_data %>% 
  pivot_longer(-1) %>% 
  mutate(name = ifelse(name == "V2", "Freq. Obs.", "Freq. Exp.")) %>% 
  ggplot(aes(V1, value)) +
  geom_line(aes(color = name, linetype = name), size = 0.75) +
  geom_point(aes(color = name, shape = name), size = 2) + 
  scale_x_continuous( limits = c(0,22), breaks = seq(0,22,2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.8), breaks = seq(0,0.8,0.1)) +
  theme_bw() +
  scale_color_manual(values = c("red", "darkgreen")) +
  scale_shape_manual(values = c(16,8)) +
  labs(x = "Pairwise Differences",
       y = "",
       subtitle = "India COI-RS") +
  theme(
    legend.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = NA),
    text=element_text(family="sans-serif")
  ) -> coi_rs_plot
  

# seq color plot====



readDNAStringSet(here::here("data","Gene conversion.fas")) -> s

tibble(name = names(s),
       seq = paste(s)) %>% 
  separate(col = seq, into = as.character(1:491), sep = "") %>%
  mutate(`1` = `1` %>% na_if("")) -> df

t(apply(df, 1, function(x) c(x[!is.na(x)], x[is.na(x)]))) -> df[] 

all_na <- function(x){
  any(!is.na(x))
}

df %>%
  select_if(all_na) -> df1


df1 %>% 
  tidyr::pivot_longer(-1, names_to = "seq_no") %>% 
  mutate(name = name %>% str_trim() %>% factor(levels = c("India_haplotype2_CS", "India_haplotype21", "India_haplotype20","India_haplotype1_RS")),
         seq_no = seq_no %>% as.numeric()) -> df2

splitter <- function(seq_no){
  if(seq_no <= 100){
    "1-100"
  }else if(seq_no > 100 & seq_no <=200){
    "101-200"
  }else if(seq_no > 200 & seq_no <=300){
    "201-300"
  }else if(seq_no > 300 & seq_no <=400){
    "301-400"
  }else if(seq_no > 400 & seq_no <=500){
    "401-500"
  }
}

df2 %>% 
  mutate(cat = seq_no %>% map_chr(splitter)) %>% 
  ggplot() +
  geom_tile(aes(x = seq_no, y = name, fill = value), colour = "black") +
  facet_wrap(~cat, scales = "free_x", nrow = 5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,500,2))




rules <- function(s1, s2, s3){
  if(s1 == s2 & s1 == s3 ){
    "gray"
  } else if(s1 == s3 & s2 != s3){
    "red"
  } else if(s1 != s3 & s2 == s3){
    "green"
  } else if(s1 != s3 & s2 != s3){
    "black"
  }
}


rules1 <- function(s1, s2, s3, s4, colour){
  if(s1 == s2 & s1 == s3 & s1 == s4){
    "gray"
  } else{
    colour
  }
}


df1 %>% 
  mutate(name = name %>% str_trim()) %>% 
  t() %>% 
  as_tibble() %>%
  row_to_names(1) %>% 
  mutate(India_haplotype15_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype20), rules),
         India_haplotype16_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21), rules),
         India_haplotype1_RS_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21, India_haplotype20, "red"), rules1),
         India_haplotype2_CS_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21, India_haplotype20, "green"), rules1)) %>% 
  dplyr::select(-(1:4)) %>%
  mutate(r_no = row_number()) %>%
  pivot_longer(-5) %>%
  mutate(alph = ifelse(value == "gray70", 0.5, 1)) %>% 
  mutate(cat = r_no %>% map_chr(splitter)) %>% 
  mutate(name = name %>% 
           str_trim() %>% 
           str_remove("_colour") %>% 
           factor(levels = c("India_haplotype2_CS", "India_haplotype16", "India_haplotype15","India_haplotype1_RS"))) %>% 
  ggplot() +
  geom_tile(aes(x = r_no, y = name, fill = value, alpha = alph), colour = "black") +
  facet_wrap(~cat, scales = "free_x", nrow = 5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,500,3)) +
  scale_fill_identity() + 
  scale_alpha_identity() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8, colour = "black", face = "bold"),
    axis.text.x = element_text(size = 6, face = "bold"),
    text=element_text(family="sans-serif"),
    strip.background = element_rect(fill="white"),
    strip.text.x = element_text(margin = margin(0.07,0,0.07,0, "cm"))
  )  -> seq_color
  


(india_plot/(count_segment_plot+coi_rs_plot)/seq_color) +
  plot_layout(heights = c(3,1.3,2.7)) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"),
        text=element_text(family="sans-serif")) -> pp
ggsave(here::here("output","figure3.png"), pp, dpi = 300, height = 15,width = 12, units = "in")
  