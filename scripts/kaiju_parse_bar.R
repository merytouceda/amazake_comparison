# parse and plot the amazake species

# 0. set up
# install.packages("tidyverse")
library(tidyverse) # data wrangling data frames
library(see)
library(vegan)
library(ggalluvial)
library(car)
library(performance)
library(lme4)

# set working directory
setwd("/Volumes/Elements20/amazake/kaiju/unfiltered/tsv") # change to directory where you have your data
# inport metadata
metadata <- read.csv("/Volumes/Elements20/amazake/metadata.csv")



# take a look at one file
d <- read.table("AAO-0H-1.nr_euk.summary.tsv", sep = "\t", header = T)


# 1.Generate big table by uniting all sample tables

#initialize empty table
d0 <- data.frame(matrix(NA, nrow = 0, ncol = 3))
colnames <- c("sample_id", "taxon_name", "reads")
colnames(d0) <- colnames

# iterate over each file formating and joining to d0
# create list of files to iterate over
files <- list.files(pattern="*.tsv", full.names=F, recursive=FALSE)
for (file in files){
  filename <- basename(file) # make object with file name to retrieve sample name later
  d <- tryCatch(read_tsv(file), error=function(e) NULL) # open file
  # extract sample name from name of file
  samplename <- str_extract(filename, "[^.]*") # extract the string before the first "."
  
  d <- d %>%
    dplyr::select(c("taxon_name", "reads")) %>% #select columns of interest
    mutate(sample_id=samplename) # create sample name column
  
  # append to form a big dataset
  d0 <-rbind(d0, d)
  
}


# create abundance table (rows = microorganisms, columns = samples, values = reads of microorganim in that sample)
tax_count <- d0 %>%
  na.omit() %>%
  pivot_wider(id_cols = taxon_name, names_from = sample_id, values_from = reads) %>%
  replace(., is.na(.), 0) %>%
  #dplyr::select(-c("NA")) %>%
  separate(taxon_name, into = c("A1", "A2", "A3","A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14", "A15", "A16"), sep = ";")

write.csv(tax_count, "/Volumes/BunnyBike/amazake/amazake_count_table_absolute_unfiltered.csv")



tax_count <- read.csv("/Volumes/BunnyBike/amazake/amazake_count_table_absolute_unfiltered.csv")

tax_count <- tax_count %>%
  select(-X)
# 2. Calculate number of different species (also known as richness)
# first we calculate the mininum number of reads per sample and we use it to rarify the richness calculation
summary(colSums(tax_count[17:28]))
hist(colSums(tax_count[17:28]))
metadata$rich <- specnumber(rrarefy(t(tax_count[17:28]), 8027727))


# make variables into factors (works better with plotting)
metadata$filename <- as.factor(metadata$filename)
metadata$timepoint <- as.factor(metadata$timepoint)
metadata$organism <- as.factor(metadata$organism)
metadata$pair <- as.factor(metadata$pair)

# make a boxplot with  
ggplot(metadata, aes(x = timepoint, y = rich)) +
  geom_boxplot(aes(fill = organism, color = organism), alpha = 0.4) +
  geom_point(aes(color = organism), size = 2.5) + 
  geom_line(aes(group = pair), linewidth = 0.1) +
  xlab("Sample")+
  ylab("Number of different species")+
  scale_color_see()+
  scale_fill_see()+
  theme_bw() +
  facet_wrap(~organism)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("/Volumes/BunnyBike/amazake/amazake_rich.pdf", device = "pdf", width = 4, height = 5 , units = "in")

# linear model for organism
anova(lm(rich ~ organism, data = metadata))

# fixed effect for timepoint
metadata_AAO <- metadata %>%
  filter(organism == "AAO")
anova(lmer(rich ~ timepoint + (1|pair), data = metadata_AAO)) 
Anova(lmer(rich ~ timepoint + (1|pair), data = metadata_AAO))

metadata_ABA<- metadata %>%
  filter(!organism == "AAO")
anova(lmer(rich ~ timepoint + (1|pair), data = metadata_ABA)) 
Anova(lmer(rich ~ timepoint + (1|pair), data = metadata_ABA))





# 3. Calculate community composition similarity
counts.bray <- vegdist(t(tax_count[17:28]), method="bray") # using species level
#counts.bray <- vegdist(t(tax_count_genus), method="bray") # using genus level

#ordination (non-multidimensional scaling)
counts.nmds <- metaMDS(counts.bray, k=2, try = 100)
metadata$Axis01 = counts.nmds$points[,1]
metadata$Axis02 = counts.nmds$points[,2]
counts.nmds$stress #0.16


ggplot(metadata, aes(Axis01, Axis02))+
  geom_point(aes(color=organism, shape = as.factor(timepoint)), size=3)+
  geom_line(aes(group = pair), linewidth = 0.1) +
  #stat_ellipse(aes(color=organism)) +
  scale_color_see()+
  theme_classic()+
  theme(legend.position="bottom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/amazake/amazake_beta_genus.pdf", device = "pdf", width = 5, height = 5 , units = "in")


adonis2(counts.bray ~ metadata$organism * timepoint, data=metadata, permutations = 999, method="bray")



# 4. Plot the taxonomic composition
# fill the emtpy cells with NA

##################################################################################### AAO
# put NA where empty spaces
tax_count$A15 <- gsub(tax_count$A15, pattern = "^$", replacement = "NA")

tax_count_genus <- tax_count %>%
  # fill column 16 (the species column) with 15 if empty
  mutate(A16 = ifelse(A16 %in% "", A15, A16)) %>%
  dplyr::select(-c("A1", "A2", "A3","A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",  "A13", "A14", "A15")) %>%
  na.omit() %>%
  filter(!A16 == "NA") %>% # needed for previously named NA
  group_by(A16) %>%
  summarize_all(sum) %>%
  column_to_rownames(var = "A16")  %>%
  rownames_to_column(var = "Genus") 
  


amazake_tax_plot <- tax_count_genus  %>%
  #rownames_to_column(var = "phylum")%>%
  pivot_longer(!Genus, names_to= "Sample", values_to = "count") %>%
  mutate(
    fermentation = case_when(
      str_detect(Sample, "^AAO") ~ "AAO",
      TRUE ~ "BA")) %>%
  mutate(
    timepoint = case_when(
      str_detect(Sample, "0H") ~ "0", 
      TRUE ~ "16")) %>%
  dplyr::select(-c("Sample")) %>%
  group_by(timepoint, fermentation, Genus) %>%
  summarize(count = sum(count))


# divide in 4, one per fermentation and timepoint and calculate percentages
# AAO
AAO_phyla <- amazake_tax_plot %>%
  filter(fermentation== "AAO") 


########## AAO-0
AAO_0_phyla <- AAO_phyla %>%
  filter(timepoint== "0") %>%
  arrange(-count)

# separate into top 10 and other
AAO_0_phyla_top19<- AAO_0_phyla %>%
  slice(1:9)
AAO_0_phyla_other<- AAO_0_phyla%>%
  slice(10:669) 

# collapse other
AAO_0_phyla_other_total <- as.data.frame(t(c("0", "AAO", "Other", sum(AAO_0_phyla_other$count))))
colnames(AAO_0_phyla_other_total) <- colnames(AAO_0_phyla_top19)

# make "count" into same class so we can merge them
AAO_0_phyla_other_total$count <- as.integer(AAO_0_phyla_other_total$count)
AAO_0_phyla_top19$count <- as.integer(AAO_0_phyla_top19$count)

# merge top 10 and collapse other to make the table to plot
AAO_0_phyla <- AAO_0_phyla_top19 %>%
  bind_rows(AAO_0_phyla_other_total) %>%
  mutate(pct = count/sum(count)) 

# plot
# AAO-0
AAO_0_phyla$ymax <- cumsum(AAO_0_phyla$pct)
# compute bottom of each rectangle
AAO_0_phyla$ymin <- c(0, head(AAO_0_phyla$ymax, n=-1))
# compute label position
AAO_0_phyla$label_position <- (AAO_0_phyla$ymax + AAO_0_phyla$ymin) /2

ggplot(AAO_0_phyla, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Genus)) +
  geom_rect()+
  #coord_polar(theta = "y")+
  xlim(c(2,4))+
  scale_fill_see()+
  theme_void() 
ggsave("/Volumes/BunnyBike/amazake/AAO_0_tax.pdf", device = "pdf", width = 6, height = 8 , units = "in")




####### AAO-16
AAO_16_phyla <- AAO_phyla %>%
  filter(timepoint== "16") %>%
  arrange(-count)

# separate into top 10 and other
AAO_16_phyla_top19<- AAO_16_phyla %>%
  slice(1:9)
AAO_16_phyla_other<- AAO_16_phyla%>%
  slice(10:669) 

# collapse other
AAO_16_phyla_other_total <- as.data.frame(t(c("16", "AAO", "Other", sum(AAO_16_phyla_other$count))))
colnames(AAO_16_phyla_other_total) <- colnames(AAO_16_phyla_top19)

# make "count" into same class so we can merge them
AAO_16_phyla_other_total$count <- as.integer(AAO_16_phyla_other_total$count)
AAO_16_phyla_top19$count <- as.integer(AAO_16_phyla_top19$count)

# merge top 10 and collapse other to make the table to plot
AAO_16_phyla <- AAO_16_phyla_top19 %>%
  bind_rows(AAO_16_phyla_other_total) %>%
  mutate(pct = count/sum(count)) 



# AAO-16
AAO_16_phyla$ymax <- cumsum(AAO_16_phyla$pct)
# compute bottom of each rectangle
AAO_16_phyla$ymin <- c(0, head(AAO_16_phyla$ymax, n=-1))
# compute label position
AAO_16_phyla$label_position <- (AAO_16_phyla$ymax + AAO_16_phyla$ymin) /2

ggplot(AAO_16_phyla, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Genus)) +
  geom_rect()+
  #coord_polar(theta = "y")+
  xlim(c(2,4))+
  scale_fill_see()+
  theme_void() 
ggsave("/Volumes/BunnyBike/amazake/AAO_16_tax.pdf", device = "pdf", width = 6, height = 8 , units = "in")




###






#################################################################################### BA
#tax_count$A11 <- gsub(tax_count$A11, pattern = "^$", replacement = "NA")
tax_count_prefam <- tax_count %>%
  mutate(A11 = ifelse(A11 %in% "", A10, A11)) %>%
  mutate(A11 = ifelse(A11 %in% "^$", A10, A11)) %>%
  mutate(A11 = ifelse(A11 %in% " ", A10, A11)) %>%
  mutate(A11 = ifelse(is.na(A11), A10, A11)) %>%
  mutate(A11 = ifelse(A11 %in% "", A9, A11)) %>%
  mutate(A11 = ifelse(A11 %in% "^$", A9, A11)) %>%
  mutate(A11 = ifelse(A11 %in% " ", A9, A11)) %>%
  mutate(A11 = ifelse(is.na(A11), A9, A11)) %>%
  mutate(A11 = ifelse(A11 %in% "", A8, A11)) %>%
  mutate(A11 = ifelse(A11 %in% "^$", A8, A11)) %>%
  mutate(A11 = ifelse(A11 %in% " ", A8, A11)) %>%
  mutate(A11 = ifelse(is.na(A11), A8, A11)) %>%
  dplyr::select(-c("A1", "A2", "A3","A4", "A5", "A6", "A7", "A8", "A9", "A10", "A12",  "A13","A14", "A15", "A16")) %>%
  na.omit() %>%
  filter(!A11 == "NA") %>% # needed for previously named NA
  group_by(A11) %>%
  summarize_all(sum) %>%
  column_to_rownames(var = "A11")  %>%
  rownames_to_column(var = "Tax") 


amazake_tax_plot <- tax_count_prefam  %>%
  pivot_longer(!Tax, names_to= "Sample", values_to = "count") %>%
  mutate(
    fermentation = case_when(
      str_detect(Sample, "^AAO") ~ "AAO",
      TRUE ~ "BA")) %>%
  mutate(
    timepoint = case_when(
      str_detect(Sample, "0H") ~ "0", 
      TRUE ~ "16")) %>%
  dplyr::select(-c("Sample")) %>%
  group_by(timepoint, fermentation, Tax) %>%
  summarize(count = sum(count))




# BA
BA_phyla <- amazake_tax_plot %>%
  filter(fermentation== "BA") %>%
  arrange(-count)
  

# BA-0
BA_0_phyla <- BA_phyla %>%
  filter(timepoint== "0") %>%
  arrange(-count)

# separate into top 10 and other
BA_0_phyla_top19<- BA_0_phyla %>%
  slice(1:9)
BA_0_phyla_other<- BA_0_phyla%>%
  slice(10:816) 

# collapse other
BA_0_phyla_other_total <- as.data.frame(t(c("0", "BA", "Other", sum(BA_0_phyla_other$count))))
colnames(BA_0_phyla_other_total) <- colnames(BA_0_phyla_top19)

# make "count" into same class so we can merge them
BA_0_phyla_other_total$count <- as.integer(BA_0_phyla_other_total$count)
BA_0_phyla_top19$count <- as.integer(BA_0_phyla_top19$count)

# merge top 10 and collapse other to make the table to plot
BA_0_phyla <- BA_0_phyla_top19 %>%
  bind_rows(BA_0_phyla_other_total) %>%
  mutate(pct = count/sum(count)) 

# plot
BA_0_phyla$ymax <- cumsum(BA_0_phyla$pct)
# compute bottom of each rectangle
BA_0_phyla$ymin <- c(0, head(BA_0_phyla$ymax, n=-1))
# compute label position
BA_0_phyla$label_position <- (BA_0_phyla$ymax + BA_0_phyla$ymin) /2

ggplot(BA_0_phyla, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Tax)) +
  geom_rect()+
  #coord_polar(theta = "y")+
  xlim(c(2,4))+
  scale_fill_see()+
  theme_void() 
ggsave("/Volumes/BunnyBike/amazake/BA_0_tax.pdf", device = "pdf", width = 6, height = 8 , units = "in")




# BA-16
BA_16_phyla <- BA_phyla %>%
  filter(timepoint== "16") %>%
  arrange(-count)

# separate into top 10 and other
BA_16_phyla_top19<- BA_16_phyla %>%
  slice(1:9)
BA_16_phyla_other<- BA_16_phyla%>%
  slice(10:816) 

# collapse other
BA_16_phyla_other_total <- as.data.frame(t(c("16", "BA", "Other", sum(BA_16_phyla_other$count))))
colnames(BA_16_phyla_other_total) <- colnames(BA_16_phyla_top19)

# make "count" into same class so we can merge them
BA_16_phyla_other_total$count <- as.integer(BA_16_phyla_other_total$count)
BA_16_phyla_top19$count <- as.integer(BA_16_phyla_top19$count)

# merge top 10 and collapse other to make the table to plot
BA_16_phyla <- BA_16_phyla_top19 %>%
  bind_rows(BA_16_phyla_other_total) %>%
  mutate(pct = count/sum(count)) 



# plot
BA_16_phyla$ymax <- cumsum(BA_16_phyla$pct)
# compute bottom of each rectangle
BA_16_phyla$ymin <- c(0, head(BA_16_phyla$ymax, n=-1))
# compute label position
BA_16_phyla$label_position <- (BA_16_phyla$ymax + BA_16_phyla$ymin) /2

ggplot(BA_16_phyla, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Tax)) +
  geom_rect()+
  #coord_polar(theta = "y")+
  xlim(c(2,4))+
  scale_fill_see()+
  theme_void() 
ggsave("/Volumes/BunnyBike/amazake/BA_16_tax.pdf", device = "pdf", width = 6, height = 8 , units = "in")





### Alluvial plot
tax_alluvial <- tax_count %>%
  select(c("A2", "A3","A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",  "A13","A14")) %>%
  unite(col = "classification", A2:A14, sep = ";") %>%
  group_by(classification) %>%
  summarize(freq = n()) %>%
  separate(classification, sep = ";", into= c("A2", "A3","A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",  "A13","A14")) %>%
  select(-c("A2", "A3", "A5", "A7", "A9")) %>%
  arrange(desc(freq)) %>%
  slice_head(n = 100)


ggplot(as.data.frame(tax_alluvial),
       aes(y = freq,
           axis1 = A4, axis2 = A6, axis3 = A8, axis4 = A10,
           axis5 = A11, axis6 = A12, axis7 = A13,
           axis8= A14
           )) +
  geom_alluvium(aes(fill = A14),
                width = 1/8, knot.pos = 0, reverse = FALSE) +
  scale_fill_see()+
  guides(fill = "none") +
  geom_stratum(alpha = .25, width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:8, labels = 
                       c("A4", "A6", "A8","A10", "A11", "A12",  "A13","A14")) 
ggsave("/Volumes/BunnyBike/amazake/amazake_tax_alluvial_color_genus.pdf", device = "pdf", width = 15, height = 10 , units = "in")

# improve by: grouping the "other"
# improve by: divide by sample

