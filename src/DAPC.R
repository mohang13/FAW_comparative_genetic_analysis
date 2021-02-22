library(tidyverse)
library(seqinr)
library(ape)
library(adegenet)
library(patchwork)
library(huxtable)

# read fasta to DNAbin
fasta2DNAbin(here::here("data","Revised manuscript.fas")) -> seq

# convert DNAbin to genind
seq %>% 
  DNAbin2genind() -> seq_genind

# DAPC calculation ----

# Assigning population groups manually
c(rep("America", 163), rep("Africa", 150), rep("India", 183), rep("Asia-II", 76)) %>% 
  factor(levels = c("America", "Africa", "India", "Asia-II"))-> groups

pop(seq_genind) <- groups

dapc(seq_genind) -> dapc1 # Initial DAPC calculation

summary(dapc1)

optim.a.score(dapc1) -> oa_score # To find optimal PCs to retain

dapc(seq_genind, n.pca = 8, n.da = 3) -> dapc2 # DAPC calculation provided the parameters


# DAPC plot 1 -----
my_pal <- RColorBrewer::brewer.pal(n=4, name = "Dark2")

dapc2$ind.coord %>% 
  as_tibble() %>% 
  mutate(Group = dapc2$grp) %>%
  ggplot(aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  labs(fill = "Groups", colour = "Groups") +
  scale_color_manual(values=c(my_pal)) + 
  scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
  scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) +
  scale_y_continuous(breaks = seq(-5, 7, 1),  limits = c(-5, 7)) +
  stat_ellipse(type = "norm", level = 0.95) -> p1

# DAPC individual facet plot -----
dapc2$ind.coord %>% 
  as_tibble() %>% 
  mutate(Group = dapc2$grp) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  scale_color_manual(values=c(my_pal)) + 
  scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
  scale_x_continuous(breaks = seq(-5, 5, 1), limits = c(-5, 5)) + 
  scale_y_continuous(breaks = seq(-2, 7, 1), limits = c(-2, 7)) +
  facet_wrap(~Group) -> facet_plot


# DAPC after randomizing groups -----

seq %>% 
  DNAbin2genind() -> rand_seq

sample(c("America", "Africa", "India", "Asia-II"), 572, replace = T) %>% 
  factor(levels = c("America", "Africa", "India", "Asia-II")) -> rand_pop

pop(rand_seq) <- rand_pop


dapc(rand_seq, n.pca = 5, n.da = 3) -> dapc_rand


summary(dapc_rand)

optim.a.score(dapc_rand) 


# Randomized DAPC plot ----
my_pal <- RColorBrewer::brewer.pal(n=4, name = "Dark2")

dapc_rand$ind.coord %>% 
  as_tibble() %>% 
  mutate(Group = dapc_rand$grp %>% factor(levels = c("America", "Africa", "India", "Asia-II"))) %>% 
  ggplot(aes(x = LD1, y = LD2, color = Group, fill = Group)) +
  geom_point(size = 4, shape = 21) +
  theme_bw() +
  labs(fill = "Groups", colour = "Groups") +
  scale_color_manual(values=c(my_pal)) + 
  scale_fill_manual(values=c(paste(my_pal, "66", sep = ""))) +
  stat_ellipse(type = "norm", level = 0.95) -> p2


# Posterior probability plot -----

dapc2$posterior %>% 
  as_tibble() %>% 
  mutate(no = row_number(),
         Group = groups) %>% 
  pivot_longer(-c(5,6)) %>% 
  mutate(name = name %>% factor(levels = c("America", "Africa", "India", "Asia-II"))) %>% 
  ggplot(aes(x = no, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "fill", width = 1) +
  facet_grid(.~Group, scales = "free_x", space = "free_x") +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,536, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = my_pal) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(x = "Samples",
       y = "Posterior Membership Propability",
       fill = "Groups") -> p3



# Combining plots using patchwork ----

(p1 + p2) + 
  plot_layout(guides = "collect") -> pp

(pp / p3) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"),
        text=element_text(family="sans-serif")) -> dapc_final_plots

dapc_final_plots +
  ggsave(here::here("output","figure_6_dapc_plot.png"), width = 11, height = 7.5, units = "in", dpi = 400)


# Admixture table ----

dapc2$posterior %>% 
  as_tibble(rownames = "Sample") %>% 
  mutate(no = row_number(),
         prior_group = as.character(groups),
         assigned_group = as.character(dapc2$assign),
         Sample = Sample %>% str_trim()) %>% 
  mutate(posterior = ifelse(America >= 0.5| Africa >= 0.5| India >= 0.5| `Asia-II` >= 0.5, assigned_group, "Admixture"), # 50 % cut off
         `Same/other group` = ifelse(posterior == "Admixture", "Admixture", "Same")) %>% 
  dplyr::select(Sample, prior_group, posterior, `Same/other group`) %>% 
  rename(Prior = prior_group,
         Posterior= posterior) -> admixture_table

admixture_table %>% 
  write_csv(here::here("output","supplementary_table1_admixture.csv"))

# membership assignment table------
here::here("functions","functions.R") %>% 
  source()

admixture_table %>% 
  crosstab(row.vars = "Prior", col.vars = "Posterior", type = "f") -> output_crosstab

output_crosstab$table %>% 
  as.matrix() %>% 
  as_hux() %>% 
  as_tibble() %>% 
  rename(`Prior/Posterior` = rownames) %>% 
  slice(-1) %>% 
  mutate(`Prior/Posterior` = `Prior/Posterior` %>% factor(levels = c("America", "Africa", "India", "Asia-II", "Sum"))) %>% 
  arrange(`Prior/Posterior`) %>% 
  mutate(across(-1, as.numeric)) %>% 
  mutate(`% admixture detected` = (Admixture/Sum)*100) -> membership_table


membership_table %>% 
  write_rds(here::here("rds","membership_table.rds"))

membership_table %>% 
  as_hux() %>% 
  set_all_borders() %>% 
  set_bold(1,everywhere) %>% 
  set_bold(everywhere, 1) %>% 
  set_align(everywhere, -1, "centre") %>% 
  quick_docx(file = here::here("output","membership_table_revised.docx"))


