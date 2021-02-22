library(tidyverse)
library(patchwork)

df <-tribble(
  ~group, ~name, ~count,
  "COIA", "COI-RS", 174, 
  "COIA", "COI-CS", 17,
  "COIB", "COI-RS", 69,
  "COIB", "COI-CSh4", 4,
  "COIB", "COI-CSh1", 0,
  "COIB", "COI-CSh2", 0,
  "COIB", "COI-CSh3", 0,
  "Tpi", "Tpi-Ca1a", 39,
  "Tpi", "Tpi-Ca2a", 6,
  "Tpi", "Tpi-Ca2b", 4,
  "Tpi", "Tpi-Ca1a/2a", 7,
  "Tpi", "Tpi-Ca2a/2b", 2,
  "Tpi", "Tpi-Ca1a/2b", 4,
  "Tpi", "Tpi-R", 4
)


df %>% 
  group_by(group) %>% 
  mutate(percent = count/sum(count)*100) -> df1

df1 %>% 
  filter(group == "COIA") %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_col(width = 0.25) +
  geom_text(aes(label = paste0("n = ",count)), nudge_y = 5, fontface = "bold") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  labs(x = "",
       y = "Percentage") +
  annotate("label", x = 0.65, y = 95, label = "Total: 191", fontface = "bold") -> p1

df1 %>% 
  filter(group == "COIB") %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_col(width = 0.6) + 
  geom_text(aes(label = paste0("n = ",count)), nudge_y = 3, fontface = "bold") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  labs(x = "",
       y = "Percentage") +
  annotate("label", x = 0.95, y = 95, label = "Total: 73", fontface = "bold") -> p2

df1 %>% 
  filter(group == "Tpi") %>% 
  mutate(name = name %>% factor(
    levels = c("Tpi-Ca1a", "Tpi-Ca2a", "Tpi-Ca2b", "Tpi-Ca1a/2a", "Tpi-Ca1a/2b", "Tpi-Ca2a/2b", "Tpi-R"))) %>% 
  ggplot(aes(x = name, y = percent)) +
  geom_col(width = 0.4) +
  geom_text(aes(label = paste0("n = ",count)), nudge_y = 5, fontface = "bold") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  labs(x = "",
       y = "Percentage") +
  annotate("label", x = 0.75, y = 95, label = "Total: 66", fontface = "bold") -> p3


(p1 + p2)/ p3 +
  plot_annotation(tag_levels = 'a') &
  theme_classic() +
  theme(axis.text = element_text(face = "bold"), 
        axis.title = element_text(face = "bold"),
        plot.tag = element_text(face = "bold")) -> pp
ggsave(here::here("output","figure2.png"), pp, width = 11, height = 7.5, units = "in", dpi = 300)

