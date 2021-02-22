library(tidyverse)
library(adegenet)
library(poppr) 
library(ape)
library(huxtable)


fasta2DNAbin(here::here("data","Revised manuscript.fas")) -> seq



#insect::subset.DNAbin(seq,subset = c(rep(T, 163), rep(F,150, rep(F, 183), rep(F, 76))))

# convert DNAbin to genind
seq %>% 
  DNAbin2genind() -> seq_genind


tibble(ind = 1:572,
       pop = c(rep("America", 163), rep("Africa", 150), rep("India", 183), rep("East Asia", 76))
) -> hierarchy_info

tibble(ind = 1:572,
       pop = sample(c("America", "Africa", "India", "East Asia"), 572, replace = T)
) -> hierarchy_random

strata(seq_genind) <- hierarchy_info


nameStrata(seq_genind) <- ~individual/group

table(strata(seq_genind, ~group, combine = F)) 

poppr.amova(seq_genind, ~group) -> seq_amova

randtest(seq_amova, nrepet = 1000) -> seq_amova_test

plot(seq_amova_test)


# US vs Ind -----------------------------------------------------------------

# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(T, 163), rep(F,150), rep(T, 183), rep(F, 76))) -> seq_ind_us

# convert DNAbin to genind
seq_ind_us %>% 
  DNAbin2genind() -> genind_ind_us

tibble(ind = 1:346,
       pop = c(rep("America", 163), rep("India", 183))
) -> hierarchy_info_ind_us

strata(genind_ind_us) <- hierarchy_info_ind_us

nameStrata(genind_ind_us) <- ~individual/group

poppr.amova(genind_ind_us, ~group) -> seq_amova_us_india

randtest(seq_amova_us_india, nrepet = 1000) -> seq_amova_test_us_india

plot(seq_amova_test_us_india)


# us vs africa -----------

# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(T, 163), rep(T,150), rep(F, 183), rep(F, 76))) -> seq_us_africa


#tibble(ID=labels(seq_us_africa)) %>% view()


# convert DNAbin to genind
seq_us_africa %>% 
  DNAbin2genind() -> genind_us_africa

tibble(ind = 1:313,
       pop = c(rep("America", 163), rep("Africa", 150))
) -> hierarchy_info_us_africa

strata(genind_us_africa) <- hierarchy_info_us_africa

nameStrata(genind_us_africa) <- ~individual/group

poppr.amova(genind_us_africa, ~group) -> seq_amova_us_africa

randtest(seq_amova_us_africa, nrepet = 1000) -> seq_amova_test_us_africa

plot(seq_amova_test_us_africa)


# africa vs india -----------

# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(F, 163), rep(T,150), rep(T, 183), rep(F, 76))) -> seq_africa_india

# convert DNAbin to genind
seq_africa_india %>% 
  DNAbin2genind() -> genind_africa_india

tibble(ind = 1:333,
       pop = c(rep("Africa", 150), rep("India", 183))
) -> hierarchy_info_africa_india

strata(genind_africa_india) <- hierarchy_info_africa_india

nameStrata(genind_africa_india) <- ~individual/group

poppr.amova(genind_africa_india, ~group) -> seq_amova_africa_india

randtest(seq_amova_africa_india, nrepet = 1000) -> seq_amova_test_africa_india

plot(seq_amova_test_africa_india)


# US vs E Asia ------
# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(T, 163), rep(F,150), rep(F, 183), rep(T, 76))) -> seq_us_easia

# convert DNAbin to genind
seq_us_easia %>% 
  DNAbin2genind() -> genind_us_easia

tibble(ind = 1:239,
       pop = c(rep("America", 163), rep("East Asia", 76))
) -> hierarchy_info_us_easia

strata(genind_us_easia) <- hierarchy_info_us_easia

nameStrata(genind_us_easia) <- ~individual/group

poppr.amova(genind_us_easia, ~group) -> seq_amova_us_easia

randtest(seq_amova_us_easia, nrepet = 1000) -> seq_amova_test_us_easia

plot(seq_amova_test_us_easia)

# Africa vs E Asia ------
# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(F, 163), rep(T,150), rep(F, 183), rep(T, 76))) -> seq_africa_easia

# convert DNAbin to genind
seq_africa_easia %>% 
  DNAbin2genind() -> genind_africa_easia

tibble(ind = 1:226,
       pop = c(rep("Africa", 150), rep("East Asia", 76))
) -> hierarchy_info_africa_easia

strata(genind_africa_easia) <- hierarchy_info_africa_easia

nameStrata(genind_africa_easia) <- ~individual/group

poppr.amova(genind_africa_easia, ~group) -> seq_amova_africa_easia

randtest(seq_amova_africa_easia, nrepet = 1000) -> seq_amova_test_africa_easia

plot(seq_amova_test_africa_easia)


# India vs E Asia ------
# read fasta to DNAbin
insect::subset.DNAbin(seq, subset = c(rep(F, 163), rep(F,150), rep(T, 183), rep(T, 76))) -> seq_india_easia
# convert DNAbin to genind
seq_india_easia %>% 
  DNAbin2genind() -> genind_india_easia

tibble(ind = 1:259,
       pop = c(rep("India", 183), rep("East Asia", 76))
) -> hierarchy_info_india_easia

strata(genind_india_easia) <- hierarchy_info_india_easia

nameStrata(genind_india_easia) <- ~individual/group

poppr.amova(genind_india_easia, ~group) -> seq_amova_india_easia

randtest(seq_amova_india_easia, nrepet = 1000) -> seq_amova_test_india_easia

plot(seq_amova_test_india_easia)


# table -----------

amova_table <- function(amova_result, significance_result, test_name){
  
  amova_result$results -> results
  amova_result$componentsofcovariance -> variance
  
  tibble(Groups = test_name,
         Source = row.names(results),
         Df = results$Df,
         SS = results$`Sum Sq`,
         VarianceComponents = variance$Sigma,
         `TotalVariance(%)` = variance$`%`,
         `P Value` = c(significance_result$pvalue, NA, NA))
}


amova_table(seq_amova, seq_amova_test, "All") %>% 
  bind_rows(amova_table(seq_amova_us_india, seq_amova_test_us_india, "America and India")) %>% 
  bind_rows(amova_table(seq_amova_us_africa, seq_amova_test_us_africa, "America and Africa")) %>% 
  bind_rows(amova_table(seq_amova_africa_india, seq_amova_test_africa_india, "Africa and India")) %>%
  bind_rows(amova_table(seq_amova_us_easia, seq_amova_test_us_easia, "America and E Asia")) %>% 
  bind_rows(amova_table(seq_amova_africa_easia, seq_amova_test_africa_easia, "Africa and E Asia")) %>% 
  bind_rows(amova_table(seq_amova_india_easia, seq_amova_test_india_easia, "India and E Asia")) %>% 
  mutate(Groups = replace(Groups, duplicated(Groups), NA)) %>% 
  mutate(across(where(is.numeric), round, 3),
         across(where(is.numeric), as.character)) %>% 
  mutate(Source = Source %>% str_replace_all("samples", "groups"))-> amova_final_table


amova_final_table %>% 
  write_rds(here::here("rds","amova_final_table.rds"))

amova_final_table %>% 
  as_hux() %>% 
  set_all_borders() %>% 
  set_bold(1, everywhere) %>% 
  set_bold(everywhere, 1) %>% 
  set_align(everywhere, -(1:2), "centre") %>% 
  merge_cells(2:4,1) %>% 
  merge_cells(5:7,1) %>% 
  merge_cells(8:10,1) %>% 
  merge_cells(11:13,1) %>% 
  merge_cells(14:16,1) %>%
  merge_cells(17:19,1) %>%
  merge_cells(20:22,1) %>%
  set_font_size(everywhere,everywhere, 9) %>% 
  quick_docx(file = here::here("output","Amova_table.docx"))
