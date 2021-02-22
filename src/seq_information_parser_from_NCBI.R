library(tidyverse)
library(rentrez)
library(Biostrings)
library(janitor)
library(tmaptools)


all_na <- function(x){
  any(!is.na(x))
}

readDNAStringSet(here::here("data","Revised manuscript.fas")) -> seq

names(seq) %>% 
  enframe() %>% 
  mutate(acc_no = value %>% str_remove(".*_") %>% str_trim()) %>% 
dplyr::select(acc_no) -> df

df %>% 
  pull() -> acc_nos

accession_no_to_info <- function(acc_no){
  entrez_search(db = "Nucleotide", acc_no) -> ss
  Sys.sleep(2)
  if(ss$count == 1){
    entrez_summary(db = "Nucleotide", ss$ids) -> dd
    
    unlist(dd) %>% 
      enframe() -> dpt
    
    dpt %>% 
      filter(name %in% c("accessionversion", "title", "createdate", "updatedate", "organism"))  %>% 
      t() %>% 
      as_tibble() %>% 
      row_to_names(1)-> df1
    
    dpt %>% 
      filter(name %in% c("subtype", "subname")) %>% 
      dplyr::select(value) %>% 
      separate("value", sep = "\\|" , into = as.character(1:50)) %>% 
      select_if(all_na) %>% 
      row_to_names(1) -> df2
    
    df1 %>% 
      bind_cols(df2) -> df3
    
    tibble(accession_number = acc_no,
           accession_version = tryCatch({df3$accessionversion}, warning=function(cond){return(NA)}),
           title = tryCatch({df3$title}, warning=function(cond){return(NA)}),
           create_date = tryCatch({df3$createdate}, warning=function(cond){return(NA)}),
           update_date = tryCatch({df3$updatedate}, warning=function(cond){return(NA)}),
           organism = tryCatch({df3$organism}, warning=function(cond){return(NA)}),
           country = tryCatch({df3$country}, warning=function(cond){return(NA)}),
           lat_lon = tryCatch({df3$lat_lon}, warning=function(cond){return(NA)}),
           collection_date = tryCatch({df3$collection_date}, warning=function(cond){return(NA)})) -> df4
    
  } else {
    tibble(accession_number = acc_no,
           accession_version = as.character(NA),
           title = as.character(NA),
           create_date = as.character(NA),
           update_date = as.character(NA),
           organism = as.character(NA),
           country = as.character(NA),
           lat_lon = as.character(NA),
           collection_date = as.character(NA)) -> df4
  }
  
  Sys.sleep(2)
  
  return(df4)
  
}

acc_nos %>% 
  map_dfr(accession_no_to_info) -> ffg

ffg %>% 
  mutate(lat_lon = lat_lon %>% str_remove_all(" N| W| S| E") %>% str_trim()) %>% 
  separate(col = lat_lon, into = c("lat", "lon"), sep = " ") %>% 
  mutate(location_for_extraction = ifelse(is.na(lat), country, NA),
         location_for_extraction = ifelse((is.na(country) & is.na(location_for_extraction)), "unknown", location_for_extraction)) -> meta_info 

write_csv(meta_info,"meta_info.csv")

place_to_lat_lon <- function(place_name){
  place_name %>% 
    str_replace_all(":", ",") %>% 
    str_split(",") %>% 
    unlist() %>% 
    rev() %>% 
    sapply(str_trim) %>% 
    unname() -> query
  
  geocode_OSM(query) -> qq
  
  if(length(query) < 2){
    tibble(lat = qq$coords[2], lon = qq$coords[1]) -> df
    
  } else {
    qq %>% 
      as_tibble() -> df
  }
  
  print(place_name)
  print(df)
  
  response <- ""
  while (!is.integer(response)) {
    print("Please enter a numeric input")
    tryCatch({readline(prompt = "Select place: ") %>% as.integer()}, warning=function(cond){return("")}) -> response
  }
  
  df %>% 
    dplyr::slice(response) %>% 
    dplyr::select(lat,lon) %>% 
    dplyr::rename(lat_n = lat,
                  lon_n = lon) %>% 
    mutate(place = place_name)
}

meta_info %>% 
  dplyr::select(location_for_extraction) %>% 
  filter(location_for_extraction != "unknown") %>% 
  pull() %>% 
  unique() %>% 
  map_dfr(place_to_lat_lon) -> missing_lat_lon

write_csv(missing_lat_lon,"missing_lat_lon.csv")


meta_info %>% 
  left_join(missing_lat_lon %>% dplyr::rename(location_for_extraction = place)) %>% 
  write_csv("seq_information_africa_easia.csv")

"seq_information_africa_easia.csv" %>% 
  read_csv() %>% 
  mutate(lat = ifelse(is.na(lat), lat_n, lat)) %>% 
  mutate(lon = ifelse(is.na(lon), lon_n, lon)) %>% 
  drop_na(accession_version) %>% 
  dplyr::select(-c(location_for_extraction, lat_n, lon_n)) %>%
  distinct(accession_number, .keep_all = T) %>% 
  write_csv("seq_information_africa_easia_complete.csv")

readDNAStringSet("us_easia.fas") -> seq

names(seq) %>% 
  enframe() %>% 
  mutate(acc_no = value %>% str_remove(".*_") %>% str_trim()) %>% 
  dplyr::select(acc_no) -> df

df %>% 
  dplyr::slice(1:163) %>% 
  pull() -> acc_nos

acc_nos %>% 
  map_dfr(accession_no_to_info) -> ffg

lat_sighner <- function(x){
  if(!is.na(x)){
    if(is.numeric(x)){
      x
    } else {
      if(str_detect(x, "-")){
        x %>% 
          str_remove("-") -> x
        -as.numeric(x) 
      } else {
        as.numeric(x)
      }
    }
  } else {
    NA
  } 
}




ffg %>% 
  mutate(lat_lon = lat_lon %>% 
           str_replace_all(" N| E", "_") %>% 
           str_replace_all(" W| S", "-_") %>% 
           str_trim()) %>%
  #mutate(lat_lon = lat_lon %>% str_remove_all(" N| W| S| E") %>% str_trim()) %>% 
  separate(col = lat_lon, into = c("lat", "lon"), sep = "_") %>% 
  mutate(across(c(lat,lon), map_dbl, lat_sighner)) %>% 
  mutate(location_for_extraction = ifelse(is.na(lat), country, NA),
         location_for_extraction = ifelse((is.na(country) & is.na(location_for_extraction)), "unknown", location_for_extraction)) -> meta_info 

meta_info %>% 
  dplyr::select(location_for_extraction, everything()) %>% 
  write_csv("america_raw_info1.csv")
read_csv("america_raw_info1.csv") %>% 
  dplyr::select(location_for_extraction) %>% 
  filter(location_for_extraction != "unknown") %>% 
  pull() %>% 
  unique() %>% 
  map_dfr(place_to_lat_lon) -> missing_lat_lon



read_csv("america_raw_info1.csv") %>%
  dplyr::select(everything(), location_for_extraction) %>% 
  left_join(missing_lat_lon %>% dplyr::rename(location_for_extraction = place)) %>% 
  mutate(lat = ifelse(is.na(lat), lat_n, lat)) %>% 
  mutate(lon = ifelse(is.na(lon), lon_n, lon)) %>% 
  drop_na(accession_version) %>% 
  mutate(country = ifelse(is.na(country), location_for_extraction, country)) %>%
  dplyr::select(-c(location_for_extraction, lat_n, lon_n)) %>%
  distinct(accession_number, .keep_all = T) %>% 
  write_csv("seq_information_usa_complete.csv")

"seq_information_usa_complete.csv" %>% 
  read_csv() %>% 
  mutate(lon = ifelse(lon >= 0, -lon, lon))