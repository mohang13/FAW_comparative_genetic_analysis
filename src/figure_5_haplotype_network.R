library(tidyverse)
library(ape)
library(pegas)
library(RColorBrewer)


read.table(text = "Hap1 67,23,12,15
Hap2 2,0,0,0
Hap3 61,125,146,58
Hap4 4,0,0,0
Hap5 1,0,0,0
Hap6 1,0,0,0
Hap7 1,0,0,0
Hap8 2,0,0,0
Hap9 1,0,0,0
Hap10 3,0,0,0
Hap11 1,0,0,0
Hap12 1,0,0,0
Hap13 1,0,0,0
Hap14 1,0,0,0
Hap15 1,0,0,0
Hap16 1,0,0,0
Hap17 2,0,0,0
Hap18 3,0,0,0
Hap19 1,0,0,0
Hap20 2,0,0,0
Hap21 1,0,0,0
Hap22 1,0,0,0
Hap23 1,0,0,0
Hap24 1,0,0,0
Hap25 1,0,0,0
Hap26 1,0,0,0
Hap27 0,1,0,0
Hap28 0,1,0,0
Hap29 0,0,1,0
Hap30 0,0,12,0
Hap31 0,0,1,0
Hap32 0,0,1,0
Hap33 0,0,1,0
Hap34 0,0,1,0
Hap35 0,0,1,0
Hap36 0,0,1,0
Hap37 0,0,1,0
Hap38 0,0,1,0
Hap39 0,0,1,0
Hap40 0,0,1,0
Hap41 0,0,1,0
Hap42 0,0,1,0
Hap43 0,0,0,1
Hap44 0,0,0,1
Hap45 0,0,0,1") %>% 
  as_tibble() %>% 
  rename(TraitLabels = V1) %>% 
  separate(V2, c("America", "Africa", "India", "AsiaII"), convert = T) -> freq


read.dna(here::here("data","haplotype_seq.fasta"), format="fasta") -> seq

haplotype(seq, labels = 1:45) -> seq_haps

seq_haps <- sort(seq_haps, what = "label")
haploNet(seq_haps) -> net_haps

freq %>% 
  mutate(total = America+Africa+India+AsiaII) %>% 
  pull(total) -> freq_num

plot(net_haps, size = log2(real_freq_num)+1)

attr(net_haps, "labels") %>% 
  as.numeric() %>% 
  enframe(name = NULL, value = "pegas") %>% 
  left_join(enframe(freq_num, name = "pegas")) %>% 
  pull(value) -> real_freq_num

attr(net_haps, "labels") %>% 
  enframe(name = NULL, value = "TraitLabels") %>% 
  mutate(TraitLabels = paste0("Hap", TraitLabels)) %>% 
  left_join(freq) %>%
  dplyr::select(-1) %>% 
  as.matrix() -> fm



ind.hap<-with(
  stack(setNames(attr(seq_haps, "index"), rownames(seq_haps))),
  table(hap=ind, pop=rownames(seq)[values])
) 

ind.hap<-with(
  stack(setNames(attr(seq_haps, "index"), rownames(seq_haps))),
  table(hap=ind, pop=rownames(seq)[values])
)

plot(net_haps, size=log2(real_freq_num)+1, pie=fm)
legend("topright", colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=19, ncol=2)

brewer.pal(n = 4, name = "Dark2")




plot(net_haps, size =  log2(real_freq_num)+1, 
     pie=fm, legend=F, cex = 0.8, 
     show.mutation = 0, threshold = 0, 
     bg = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), 
     lwd = net_haps %>% as.matrix.data.frame() %>% as_tibble() %>% pull(V3)) 


myplot <- recordPlot()
saveRDS(myplot, here::here("rds","haplotype_plot.rds"))

read_rds(here::here("rds","haplotype_plot.rds"))
replot()
legend(3, 13, c("America", "Africa", "India", "Asia-II"), col= c("#1B9E77", "#D95F02", "#7570B3", "#E7298A"), 
       pch=19, ncol=1, cex = 1.1, box.lty=0, pt.cex = 2, 
       bg="transparent", y.intersp=1, x.intersp = 1)

myplot <- recordPlot()
#Cairo::Cairo("test.png", units="in", width=9, height=7, type = "png", dpi=300, bg = "white", fallback_resolution=300)
tiff(here::here("output","haplotype_network.tiff"), units="in", width=9, height=7, res = 300)
replayPlot(myplot)
dev.off()

# Checking sequences -----

library(Biostrings)

here::here("data","haplotype_seq.fasta") %>% 
  readDNAStringSet() -> s

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
  filter(name %in% c("Hap3", "Hap42")) %>% 
  tidyr::pivot_longer(-1, names_to = "seq_no") %>% 
  mutate(name = name %>% str_trim(),
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
  filter(name %in% c("Hap3", "Hap42")) %>% 
  t() %>% 
  as_tibble() %>% 
  janitor::row_to_names(1) %>% 
  mutate(Hap3_C = "gray",
         Hap42_C = ifelse(Hap3 == Hap42, "gray", "red")) %>% 
  # mutate(India_haplotype15_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype20), rules),
  #        India_haplotype16_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21), rules),
  #        India_haplotype1_RS_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21, India_haplotype20, "red"), rules1),
  #        India_haplotype2_CS_colour = pmap_chr(list(India_haplotype1_RS, India_haplotype2_CS, India_haplotype21, India_haplotype20, "green"), rules1)) %>% 
  dplyr::select(-(1:2)) %>%
  mutate(r_no = row_number()) %>% 
  pivot_longer(-3) %>% 
  #mutate(alph = ifelse(value == "gray70", 0.5, 1)) %>% 
  mutate(cat = r_no %>% map_chr(splitter)) %>% 
  mutate(name = name %>% 
           str_remove("_C")) -> ddf

ddf %>% 
  count(value) %>% 
  filter(value == "red") %>% 
  pull(n)

ddf %>% 
  ggplot() +
  geom_tile(aes(x = r_no, y = name, fill = value), colour = "black") +
  facet_wrap(~cat, scales = "free_x", nrow = 5) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), breaks = seq(1,500,3)) +
  scale_fill_identity() + 
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
  ) 
  