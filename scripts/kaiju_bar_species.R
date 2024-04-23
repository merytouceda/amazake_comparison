# plot by species

library(tidyverse)
library(see)
library(ggplot2)

setwd("/Volumes/BunnyBike/amazake/")
amazake_counts <- read.csv("counts_separate_tax.csv")


# save viruses
virus_vec <- amazake_counts %>%
  filter(Domain == "Viruses") %>%
  dplyr::select(-c("Kingdom", "Phylum","Class", "Family", "Order", "Genus", "Species"))

colnames(virus_vec) <- c("Species", "AAO-0H-1",   "AAO-0H-2" ,  "AAO-0H-3" ,  "AAO-16H-1" , 
                         "AAO-16H-2" , "AAO-16H-3" , "BA-0H-1" , "BA-0H-2" ,   "BA-0H-3" ,
                         "BA-16H-1" , "BA-16H-2" ,  "BA-16H-3")


filtered_tax_count <- amazake_counts %>%
  dplyr::select(-c("Domain", "Kingdom", "Phylum","Class", "Order", "Family", "Genus")) %>%
  na.omit() #%>%
  rbind(., virus_vec)

write.csv(filtered_tax_count, "amazake_count_table_species.csv")
  
  
# plot

# inport metadata
metadata <- read_csv("/Volumes/Elements20/amazake/metadata.csv")



amazake_tax_plot <- filtered_tax_count  %>%
  #rownames_to_column(var = "phylum")%>%
  pivot_longer(!Family, names_to= "Sample", values_to = "count") %>%
  mutate(
    fermentation = case_when(
      str_detect(Sample, "^AAO") ~ "AAO",
      TRUE ~ "BA")) %>%
  mutate(
    timepoint = case_when(
      str_detect(Sample, "0H") ~ "0", 
      TRUE ~ "16")) %>%
  dplyr::select(-c("Sample")) %>%
  group_by(timepoint, fermentation, Family) %>%
  summarize(count = sum(count))

