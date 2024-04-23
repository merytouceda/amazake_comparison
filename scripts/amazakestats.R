# Statistical analysis of amazake project


# ------------------------------------------------------------------------------------
# Set up
# ------------------------------------------------------------------------------------

# load required packages
library(ggplot2)
library(dplyr)

setwd("/Volumes/MeriTSHD/amazake")
# load data
datos <- read.csv("datos.csv" , header = T )

# create dataset with only low temperature
datosL <- datos[(datos$temperature == '55'),]


datosL_sumarized <- datosL  %>% group_by(organism, time) %>% 
  summarize_at(.vars = vars("gcon", "acon", "prot"), 
               .funs = c(mean="mean")) %>%
  as.data.frame()

# ------------------------------------------------------------------------------------
# Glucose
# ------------------------------------------------------------------------------------

# ----------------------------
# Visualization
# ----------------------------

# boxplot gcon by organism 
ggplot(datosL, aes(x = organism, y = gcon))+
  geom_boxplot(aes(color = organism))+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  geom_point(size = 4)+ 
  ylab("Reducing sugars concentration (g/L)")+
  xlab("")+
  theme_minimal() 
ggsave("boxplot_gcon.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# boxplot gcon by organism color based on  time
ggplot(datosL, aes(x = organism, y = gcon))+
  geom_boxplot()+
  geom_point(aes(color = as.factor(time)), size = 4)+ 
  scale_color_manual(name= "Timepoint", values = c("orange", "orange3", "orangered2", "orangered4"))+
  ylab("Reducing sugars concentration (g/L)")+
  theme_minimal() 
ggsave("boxplot_colorstime_gcon.pdf", device = "pdf", width = 9, height = 7 , units = "in")

#line graph of glucose concentration change with time
ggplot(datosL_sumarized, aes(x = as.factor(time), y = gcon_mean, group = organism)) +
  geom_line(aes(color = organism)) +
  geom_point(aes(color = organism))+ 
  xlab("Timepoint (h)")+
  ylab("Reducing sugars concentration (g/L)")+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  theme_minimal()
ggsave("line_gcon.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# ----------------------------
# Stats
# ----------------------------

# wilcoxon test between organisms
wilcox.test(gcon_mean~organism, data = datosL_sumarized)
# paired (es decir, cada timepoint comparado con el mismo timepoint del otro bicho)
wilcox.test(gcon_mean~organism, paired = TRUE, alternative = "two.sided", data = datosL_sumarized)


# manova between time points




# ------------------------------------------------------------------------------------
# Amylase
# ------------------------------------------------------------------------------------
# ----------------------------
# Visualization
# ----------------------------

# boxplot gcon by organism 
ggplot(datosL, aes(x = organism, y = acon))+
  geom_boxplot(aes(color = organism))+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  geom_point(size = 4)+ 
  ylab("Alpha-amylase activity (units/mL)")+
  xlab("")+
  theme_minimal() 
ggsave("boxplot_acon.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# boxplot gcon by organism color based on  time
ggplot(datosL, aes(x = organism, y = acon))+
  geom_boxplot()+
  geom_point(aes(color = as.factor(time)), size = 4)+ 
  scale_color_manual(name= "Timepoint", values = c("orange", "orange3", "orangered2", "orangered4"))+
  ylab("Alpha-amylase activity (units/mL)")+
  theme_minimal() 
ggsave("boxplot_colorstime_acon.pdf", device = "pdf", width = 9, height = 7 , units = "in")


#line graph of glucose concentration change with time
ggplot(datosL_sumarized, aes(x = as.factor(time), y = acon_mean, group = organism)) +
  geom_line(aes(color = organism)) +
  geom_point(aes(color = organism))+ 
  xlab("Timepoint (h)")+
  ylab("Alpha-amylase activity (units/mL)")+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  theme_minimal()
ggsave("line_acon.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# ----------------------------
# Stats
# ----------------------------

# wilcoxon test between organisms
wilcox.test(acon_mean~organism, data = datosL_sumarized)
wilcox.test(acon_mean~organism, paired = TRUE, alternative = "two.sided", data = datosL_sumarized)

# 


# ------------------------------------------------------------------------------------
# Protein
# ------------------------------------------------------------------------------------
# --
# Visualization
# --

# boxplot gcon by organism 
ggplot(datosL, aes(x = organism, y = prot))+
  geom_boxplot(aes(color = organism))+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  geom_point(size = 4)+ 
  ylab("Protease activity (units/mL)")+
  xlab("")+
  theme_minimal() 
ggsave("boxplot_prot.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# boxplot gcon by organism color based on  time
ggplot(datosL, aes(x = organism, y = prot))+
  geom_boxplot()+
  geom_point(aes(color = as.factor(time)), size = 4)+ 
  scale_color_manual(name= "Timepoint", values = c("orange", "orange3", "orangered2", "orangered4"))+
  ylab("Protease activity (units/mL)")+
  theme_minimal() 
ggsave("boxplot_colorstime_prot.pdf", device = "pdf", width = 9, height = 7 , units = "in")


#line graph of glucose concentration change with time
ggplot(datosL_sumarized, aes(x = as.factor(time), y = prot_mean, group = organism)) +
  geom_line(aes(color = organism)) +
  geom_point(aes(color = organism))+ 
  xlab("Timepoint (h)")+
  ylab("Protease activity (units/mL)")+
  scale_color_manual(name = "Organism", values = c("deepskyblue2", "seagreen2"))+
  theme_minimal()
ggsave("line_prot.pdf", device = "pdf", width = 9, height = 7 , units = "in")


# --
# Stats
# --

# wilcoxon test between organisms
wilcox.test(prot_mean~organism, data = datosL_sumarized)
wilcox.test(prot_mean~organism, paired = TRUE, alternative = "two.sided", data = datosL_sumarized)




# ------------------------------------------------------------------------------------
# Esperimento 2, Bacillus 2 temperaturas
# ------------------------------------------------------------------------------------
# create dataset with only bacillus
datosB <- datos[(datos$organism == 'Bacillus amyloliquefaciens'),]

## ESTO NO SALE! SOLO COGE LOS DE 55
datosB_sumarized <- datosL  %>% group_by(as.factor(temperature), time) %>% 
  summarize_at(.vars = vars("gcon", "acon", "prot"), 
               .funs = c(mean="mean")) %>%
  as.data.frame()

# boxplot gcon by temperature
ggplot(datosB, aes(x = as.factor(temperature), y = prot))+
  geom_boxplot(aes(color = as.factor(temperature)))+
  scale_color_manual(name = "Temperature (F)", values = c("dodgerblue2", "tomato2"))+
  geom_point(size = 4)+ 
  xlab("")+
  ylab("Protease activity (units/mL)")+
  theme_minimal() 
ggsave("boxplot_bacillus_temp.pdf", device = "pdf", width = 9, height = 7 , units = "in")

# boxplot gcon by temperature and timepoint
ggplot(datosB, aes(x = as.factor(temperature), y = prot))+
  geom_boxplot()+
  geom_point(aes(color = as.factor(time)), size = 4)+ 
  scale_color_manual(name= "Timepoint", values = c("orange", "orange3", "orangered2", "orangered4"))+
  ylab("Protease activity (units/mL)")+
  xlab("")+
  theme_minimal() 
ggsave("boxplot_bacillus_temp_time.pdf", device = "pdf", width = 9, height = 7 , units = "in")



#line graph of glucose concentration change with time inbacillus two temps
ggplot(datosB_sumarized, aes(x = as.factor(time), y = prot_mean, group = as.factor(temperature))) +
  geom_line(aes(color = as.factor(temperature))) +
  geom_point(aes(color = as.factor(temperature)))+ 
  xlab("Timepoint (h)")+
  ylab("Protease activity (units/mL)")+
  scale_color_manual(name = "Temperature (F)", values = c("dodgerblue2", "tomato2"))+
  theme_minimal()
ggsave("line_bacilus_temp.pdf", device = "pdf", width = 9, height = 7 , units = "in")






