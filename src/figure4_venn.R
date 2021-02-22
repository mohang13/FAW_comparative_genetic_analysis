library(tidyverse)

source(here::here("functions","ggvenn.R"))
source(here::here("functions","geom_venn.R"))


read_csv(here::here("data","Shared mutation Venn.csv")) -> df

gen_label_pos_4 <- function() {
  tribble(~name, ~x,   ~y,   ~hjust, ~vjust,
          "A",   -1.5, 0.95, 1,      1,
          "B",   -0.8,  1.2, 0.5,    0,
          "C",    0.8,  1.2, 0.5,    0,
          "D",    1.5, 0.95, 0,      1)
}

ggvenn <- function(data, columns = NULL,
                   show_elements = FALSE,
                   show_percentage = TRUE,
                   fill_color = c("blue", "yellow", "green", "red"),
                   fill_alpha = .5,
                   stroke_color = "black",
                   stroke_alpha = 1,
                   stroke_size = 1,
                   stroke_linetype = "solid",
                   set_name_color = "black",
                   set_name_size = 6,
                   text_color = "black",
                   text_size = 4,
                   label_sep = ",") {
  venn <- prepare_venn_data(data, columns, show_elements, show_percentage, label_sep)
  venn$shapes %>%
    mutate(group = LETTERS[group]) %>%
    ggplot() +
    geom_polygon(aes(x = x, y = y, group = group, fill = group),
                 alpha = fill_alpha) +
    geom_polygon(aes(x = x, y = y, group = group, fill = NA),
                 color = stroke_color,
                 size = stroke_size,
                 alpha = stroke_alpha,
                 linetype = stroke_linetype) +
    geom_text(data = venn$labels,
              aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
              color = set_name_color,
              size = set_name_size,
              fontface = "bold") +
    geom_text(data = venn$texts,
              aes(x = x, y = y, label = text, hjust = hjust, vjust = vjust),
              color = text_color,
              size = text_size) +
    scale_x_continuous(limits = c(-2, 2)) +
    scale_y_continuous(limits = c(-2, 2)) +
    scale_fill_manual(values = fill_color) +
    guides(fill = FALSE) +
    coord_fixed() +
    theme_void() +
    theme(text=element_text(family="sans-serif")) 
}

df %>% 
  gather() %>% 
  dplyr::select(value) %>% 
  mutate(value = value %>% str_remove("M") %>% as.numeric()) %>% 
  arrange(desc(value))

tibble(mutation = c(paste0("M", 1:55)),
       Africa = ifelse(mutation %in% (df %>% pull(Africa)), T, F),
       America = ifelse(mutation %in% (df %>% pull(America)), T, F),
       `East Asia` = ifelse(mutation %in% (df %>% pull(`East Asia`)), T, F),
       India = ifelse(mutation %in% (df %>% pull(India)), T, F)
) %>% 
  rename(`Asia-II` = `East Asia`) %>% 
  ggvenn(text_size = 5) +
  scale_fill_brewer(palette = "Dark2") +
  ggsave(here::here("output","Figure_4_Venn_diagram.png"), dpi = 400, scale = 2, height = 5, width = 5, units = "in")

