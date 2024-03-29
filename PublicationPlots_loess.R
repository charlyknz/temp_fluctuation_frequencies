## Script to create publication ready plots ##
rm(list=ls()) #Empty the environment

#Load packages
library(readxl)
library(tidyverse)
library(car)
library(ggplot2)
library(vegan)
library(cowplot)
library(here)
library(ggpubr)

setwd("~/Desktop/MA/MA_Rcode")

#### import Mastertable  ####
Mastertable_fluctron <- read_delim("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.csv", 
                                   ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                                   trim_ws = TRUE)
#View(Mastertable_fluctron)
str(Mastertable_fluctron)
names(Mastertable_fluctron)



### subdataframe with nutrient informations only ####

nutrients_master <- Mastertable_fluctron %>% 
  dplyr::select( fluctuation, sampling, planktotron, carbon_umol_l, C_Zoo_µmol_l, Chl.a_invivo) %>%
  drop_na(C_Zoo_µmol_l) %>%
  mutate(day = sampling * 2,
         carbo = carbon_umol_l-C_Zoo_µmol_l)
names(nutrients_master)

#calculate how much percent zoo carbon amounts of total c
sum(nutrients_master$C_Zoo_µmol_l)/sum(nutrients_master$carbon_umol_l)

#change fluctuation to factor to adjust coloring
nutrients_master$fluctuation <- factor(as.factor(nutrients_master$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

rate = (mean(nutrients_master$carbo[nutrients_master$day == 4]) - mean(nutrients_master$carbo[nutrients_master$day == 0]))/4

#### Phyto- & Zooplankton biomass ####
# total C plot
total <- ggplot(nutrients_master, aes(x = day, y = carbon_umol_l, color = as.factor(fluctuation), group = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(pch = 21, size = 3)+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(Total~C ~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
total

#phytoplankton plot
phyto <- ggplot(nutrients_master, aes(x = day, y = carbo, color = as.factor(fluctuation), group = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(pch = 21, size = 3)+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(Phytoplankton~C~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
phyto
#ggsave(plot= phyto, width = 6, height = 4, file = here('images',  'growthrates_POC_loess.png'))

# see if levels are right for plotting
levels(as.factor(nutrients_master$fluctuation))

#zooplankton plot
zoo<-ggplot(nutrients_master, aes(x = day, y = log(C_Zoo_µmol_l), color = as.factor(fluctuation), group = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ x, size = 1)+  
  geom_point(pch = 21, size = 3)+
  scale_y_continuous(limits = c(0,7), breaks = c(1,3,5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(Ln~Zooplankton~C~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
zoo
#ggsave(plot= zoo, file = 'Zooplankton_growthrates_POC_lm.png', width = 6, height = 5)


# correlation plot zoo and phyto c
zoophytocorr<-nutrients_master %>%
  drop_na(C_Zoo_µmol_l) %>%
ggscatter(., y = "C_Zoo_µmol_l", x = "carbo", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman", xlab = expression(Phytoplankton~C~'['~mu*mol*~L^-1~']'),
          ylab = expression(Zooplankton~C~'['~mu*mol*~L^-1~']'))


#ggsave(zoophytocorr, file = 'AppendixFigure_zoophyto.png', width = 4, height = 3)
#########################################################################

#### Stoichiometry ####

#subdataframe 
nutrients <- dplyr::select(Mastertable_fluctron, fluctuation, sampling, planktotron, C_Zoo_µmol_l, carbon_umol_l,  nitrate_umol_l, 
                                  'diss_Nitrat+Nitrit_umol_l', 'diss_Nitrat+Nitrit_ug_l', diss_Silikat_umol_l, diss_Phosphat_umol_l,diss_Phosphat__ug_l,
                                  POP_micromol_l,  SiP_micromol_l, srp_micromol_l) %>%
  mutate(carbon_phyto = carbon_umol_l-C_Zoo_µmol_l)
  
names(nutrients)

#### calculate ratios ####

ratio <- nutrients %>%
  mutate(CN = carbon_umol_l/nitrate_umol_l,
         NP = nitrate_umol_l/POP_micromol_l,
         CP = carbon_umol_l/POP_micromol_l,
         CSi = carbon_umol_l/SiP_micromol_l,
         NSi = nitrate_umol_l/SiP_micromol_l,
         SiP = SiP_micromol_l/POP_micromol_l) %>%
  dplyr::select(fluctuation, sampling, planktotron,CN, NP, CP, CSi, NSi, NSi, SiP)%>%
  gather(key = 'ratio', value = 'value', -fluctuation, -sampling, -planktotron) %>%
  mutate(day = sampling *2) 
#### calculate RUE and Molar ratios ####
# RUE = Biomass (micromol) / TP bzw. TN 

# 1. calculate Total Phosphoros/ Nitrate
# diss Nut + Part Nut
# RUE = biomass / TP/N
# mean of RUE per treatment

## use molar values
RUE12 <- nutrients %>%
  dplyr::select(fluctuation, planktotron, sampling, srp_micromol_l, carbon_umol_l, nitrate_umol_l, POP_micromol_l, 'diss_Nitrat+Nitrit_umol_l', 'diss_Phosphat_umol_l') %>%
  dplyr::rename(diss_P = 'diss_Phosphat_umol_l',
                diss_N = 'diss_Nitrat+Nitrit_umol_l') %>%
  mutate(TP = diss_P + POP_micromol_l,
         TN = diss_N + nitrate_umol_l) %>%
  drop_na(carbon_umol_l) %>%
  mutate(P_RUE = carbon_umol_l/ TP,
         N_RUE = carbon_umol_l/ TN) %>%
  mutate(dissNP = diss_N/diss_P) %>%
  gather(key = 'nutrient', value = 'value',-dissNP, -fluctuation, -sampling, -planktotron) %>%
  mutate(day = sampling *2)

RUE12$fluctuation <- factor(as.factor(RUE12$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

dissNP <- ggplot(RUE12, aes(x = day, y = dissNP, color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch=21,size = 3,position = position_dodge2(width = .5))+
  # scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'N:P (dissolved nutrients)', color = 'Fluctuation frequency (in h)' )+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
dissNP



#### RUE plots ####

RUE_N <- ggplot(subset(RUE12, nutrient == 'N_RUE'), aes(x = day, y = value, color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch=21,size = 3,position = position_dodge2(width = .5))+
 # scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(RUE[N]), color = 'Fluctuation frequency (in h)' )+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
RUE_N
plot_grid(total, phyto, RUE_N,zoo,  label_size = 20, labels = c('A','B','C','D'), hjust= c(-3),ncol = 2, rel_heights = c(1,1,1,1), rel_widths = c(1,1,1,1))
#ggsave(plot = last_plot(), file = 'Biomass.jpeg', width = 13, height = 11)


RUE_P <- ggplot(subset(RUE12, nutrient == 'P_RUE'), aes(x = day, y = value, color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(RUE[P]), color = 'Frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=26))
RUE_P
#ggsave(plot = last_plot(), file = 'legendi.jpeg', width = 9, height = 5)

plot_grid(dissNP,RUE_P,  label_size = 20, labels = c('A','B'),  hjust = -3, vjust = 1,ncol = 2, rel_heights = c(1,1,1,1), rel_widths = c(1,1,1,1))
#ggsave(plot = last_plot(), file = 'Fig3S.jpeg', width = 10, height = 5)

#plot_grid(RUE_N, RUE_P,labels=c("A","B", 'C', 'D'),ncol = 2, label_size = 17.5, hjust = 0, vjust = 1.2)


#### Molar ratio plots####
ratio$fluctuation <- factor(as.factor(ratio$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))

CN <- ggplot(subset(ratio, ratio == 'CN') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = ' ', y = 'C:N ratio')+
  theme( panel.background = element_rect(fill = NA),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
CN

CP <- ggplot(subset(ratio, ratio == 'CP') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = '', y = 'C:P ratio')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
CP

CSi <- ggplot(subset(ratio, ratio == 'CSi') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ x, size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'C:Si ratio', color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
CSi

NP <- ggplot(subset(ratio, ratio == 'NP') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'N:P ratio', color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
NP


plot_grid(CN, CP,CSi,NP, labels=c("A","B", 'C', 'D'),ncol = 2, label_size = 20, hjust = 0, vjust = 1)
#ggsave(plot = last_plot(), file = 'MolarRatio.png', width = 11, height = 9)

## dissolved nutrients ####
diss <- Mastertable_fluctron %>%
  dplyr::select(fluctuation, planktotron, sampling, srp_micromol_l, 'diss_Nitrat+Nitrit_umol_l', 'diss_Phosphat_umol_l', "diss_Silikat_umol_l" ) %>%
  dplyr::rename(diss_P = 'diss_Phosphat_umol_l',
                diss_N = 'diss_Nitrat+Nitrit_umol_l',
                diss_Si = "diss_Silikat_umol_l" ) %>%
    gather(key = 'nutrient', value = 'value', -fluctuation, -sampling, -planktotron) %>%
  mutate(day = sampling *2) 

#check levels 
diss$fluctuation <- factor(as.factor(diss$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))


dissN <- ggplot(subset(diss, nutrient == 'diss_N') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(dissolved~N~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
dissN

dissP <- ggplot(subset(diss, nutrient == 'diss_P') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(dissolved~P~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
dissP

dissSi <- ggplot(subset(diss, nutrient == 'diss_Si') , aes( x = day, y = value,color = fluctuation))+
  geom_smooth(method = lm, se = F,formula =  y ~ splines::bs(x, 3), size = 1)+  
  geom_point(aes(color = fluctuation), pch = 21,  size = 3,position = position_dodge2(width = .5))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = expression(dissolved~Si~'['~mu*mol*~L^-1~']'), color = 'Fluctuation frequency (in h)',fill = 'Fluctuation frequency (in h)')+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
dissSi

plot_grid(dissN, dissP,dissSi, labels=c("A","B", 'C', 'D'),ncol = 3, label_size = 18, hjust = 0, vjust = 1)
#ggsave(plot = last_plot(), file = 'DissNutrients.png', width = 13, height = 5)

#### Phytoplankton Diversity ####

## import data ###
counts <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                     ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                     locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE) %>%
  drop_na(MC)
str(counts) #check import
names(counts)

#metadata
df <- data.frame( MC = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

#change variable format to numeric 
counts$cells_ml <- as.numeric(gsub(",", ".", counts$cells_ml))
counts$MC = as.character(counts$MC)
counts$cells_ml[is.na(counts$cells_ml)] <-0


#### Diversity indices ####
all_data <- left_join(counts, df, by = c('MC')) 

shannon_BV <- all_data %>%
  select(-phylum, -grid_length, -magnification, -'stripes/grids',-date, -volume, -counts) %>%
  mutate(cells_ml = gsub(',', '.', cells_ml),
         cells_ml = as.numeric(cells_ml)) %>%
  spread(key = species, value = cells_ml) %>%
  group_by(treatment_id, sampling) %>%
  separate(treatment_id, into = c('treatment', 'fluctuation'), '_')%>%  
  select(treatment, sampling, fluctuation, everything()) 

#remove NAs and exchange with 0
shannon_BV$fluctuation[is.na(shannon_BV$fluctuation)] <- 0
shannon_BV[is.na(shannon_BV)] <- 0

##calculate shannon and inverse simpson diversity index
shannon_BV$shan = diversity(shannon_BV[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan,everything() )
shannon_BV$simpson = diversity(shannon_BV[, -c(1:5)], MARGIN = 1, index='invsimpson') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan, simpson,everything() )

## calculate species richness and evenness
absence_presence <- decostand(shannon_BV[, -c(1:6)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon_BV$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)
shannon_BV$evenness = shannon_BV$shan/log(shannon_BV$no)


# Data Wrangling - reduce dataset in complexity #
diversity <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  mutate(day = 2*sampling,
         dayname = as.factor(day)) %>%
  mutate(log = log(no))

#change levels
diversity$fluctuation <- factor(as.factor(diversity$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))


#plot for species richness or shannon
specRich<-ggplot(diversity, aes(x = sampling, y = log, group = fluctuation)) +
  geom_point(aes(color = fluctuation), position = position_dodge(width = .5),pch =21, size=3)+
  geom_smooth(aes(color = fluctuation),method = lm, se = F,formula =  y ~ x, size = 1)+  
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'Ln species richness', color = 'Fluctuation frequency (in h)')+
  scale_y_continuous(limits = c(2,3), breaks = c(2,2.5,3))+
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
specRich
#ggsave(plot = specRich, width = 6, height = 4,file =  'species_richness_counting.png')

specSimp<-ggplot(diversity, aes(x = sampling, y = simpson, group = fluctuation)) +
  geom_point(aes(color = fluctuation), pch =21, size=3)+
  geom_smooth(aes(color = fluctuation),method = lm, se = F,formula =  y ~ x, size = 1)+  
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y = 'Inv Simpson index', color = 'Fluctuation frequency (in h)')+
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  theme( panel.background = element_rect(fill = NA),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
specSimp
#ggsave(plot = specSimp, width = 6, height = 4,file = 'SimpsonIndex_counting.png')


#### compositional turnover - bray curtis ####
## import data ###
counts <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                     ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                     locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE) %>%
  drop_na(MC)
str(counts) #check import
names(counts)

#metadata
df <- data.frame( MC = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

#change variable format to numeric 
counts$cells_ml <- as.numeric(gsub(",", ".", counts$cells_ml))
counts$MC = as.character(counts$MC)
counts$cells_ml[is.na(counts$cells_ml)] <-0

#merge data
all_data <- left_join(counts, df, by = c('MC')) 

#copy sampling 0 data to replicates
data0 <- filter(all_data,sampling == 0)
data0$MC[data0$MC == 4] <-1
data0$MC[data0$MC == 11] <-12
data0$MC[data0$MC == 8] <-9
data0$MC[data0$MC == 3] <-5
data0$MC[data0$MC == 6] <-7
data0$MC[data0$MC == 10] <-2

#new df to calculate mean abundance and data wrangling
pca_B <- all_data%>%
  bind_rows(., data0) %>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(cells_ml, na.rm  =T )) %>%
  ungroup()%>%
  mutate(rel = cells_ml/sum)%>%
  group_by(treatment_id, MC, sampling, species) %>%
  summarise(mean_V = mean(rel, na.rm = T),
            sd_V = sd(rel, na.rm = T),
            se_V = sd_V/sqrt(n())) %>% #calculate mean cells with biovolume
  dplyr::select(-sd_V) %>%
  drop_na(mean_V) %>%
  ungroup()%>%
  mutate(dummy = paste(treatment_id, MC, sampling, sep = ' ')) %>%
  dplyr::select(-treatment_id, -MC, -sampling, -se_V)%>%
  spread(key = species, value = mean_V) %>% #bring data in a wide format
  column_to_rownames('dummy')

pca_B[is.na(pca_B)] <- 0

#calculate bray distance using vegdist 
data.dist <- vegdist(pca_B, "bray") %>%
  broom::tidy() %>% #change to df
  separate(item1, into = c('treatment', 'MC', 'sampling'), ' ') %>%
  separate(item2, into = c('treatment2','MC2', 'sampling2'), ' ') %>%
  filter(treatment == treatment2 & MC == MC2)#%>% #extracts only treatment observations
data.dist1 <- data.dist %>% 
  filter( sampling == 0 & sampling2 != 0)  %>%
  group_by(treatment2, sampling2, MC) %>%
  summarise(mean.dist = mean(distance, na.rm = T))%>%
  separate(treatment2, into = c('treatment', 'Fluctuation'), '_') %>%
  mutate(sampling2 = as.numeric(sampling2), 
         day = 2* sampling2,
         log = log(mean.dist)) 

data.dist1$Fluctuation[is.na(data.dist1$Fluctuation)] <-0
data.dist1$Fluctuation=as.numeric(data.dist1$Fluctuation)
data.dist1$interval = 48/data.dist1$Fluctuation
data.dist1$interval[is.infinite(data.dist1$interval)] <-0
str(data.dist1)

data.plot <- data.dist1%>%
  group_by(MC,Fluctuation, sampling2, day) %>%
  summarise(mean = mean(mean.dist), sd = sd(mean.dist), se = sd/sqrt(n())) %>%
  mutate(lower.ci = mean - 1.96*se/sqrt(n()),
         upper.ci = mean + 1.96*se/sqrt(n()))
data.plot$Fluctuation <- factor(data.plot$Fluctuation, levels= c('0', '48', '36', '24', '12', '6'))
data.dist1$Fluctuation <- factor(data.dist1$Fluctuation, levels= c('0', '48', '36', '24', '12', '6'))
names(data.plot)
distance <- ggplot(data.dist1, aes(x = sampling2, y = mean.dist, group = Fluctuation))+
  #geom_errorbar(aes(ymin = mean-lower.ci, ymax = mean +upper.ci,color = Fluctuation), width = .5, position = position_dodge2(width = .5))+
  geom_point(aes(color = Fluctuation), pch =21, size=3)+
  geom_smooth(aes(color = Fluctuation),method = lm, se = F,formula =  y ~ x, size = 1)+  
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+
  # geom_line(aes(color = Fluctuation),linetype = 'dashed')+
  labs(col =  'Fluctuation frequency (in h)',x = 'Time [days]', y = 'Compositional turnover')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=26))
#save plots in grid
plot_grid(specRich, specSimp, distance, labels=c("A","B", 'C'),ncol = 3, label_size = 18, hjust = -0.5, vjust = 1.1)
#ggsave(plot = last_plot(), file = 'specDiv.jpeg', width = 12, height = 5)

#### composition: stacked barplot ####
#import  data 
counts <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                     ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                     locale = locale(decimal_mark = ","), 
                     trim_ws = TRUE) %>%
  drop_na(MC)
str(counts) #check import

#metadata
df <- data.frame( MC = c('1', '4', '2', '10', '3', '5', '6', '7', '8', '9', '11', '12'),
                  treatment_id = c('control', 'control', 'Fluctuating_36','Fluctuating_36',
                                   'Fluctuating_6', 'Fluctuating_6', 'Fluctuating_48', 'Fluctuating_48',
                                   'Fluctuating_24', 'Fluctuating_24', 'Fluctuating_12','Fluctuating_12') )

#change variable format to numeric 
counts$cells_ml <- as.numeric(gsub(",", ".", counts$cells_ml))
counts$MC = as.character(counts$MC)
counts$cells_ml[is.na(counts$cells_ml)] <-0


all_data <- left_join(counts, df, by = c('MC')) 
names(all_data)
#new df to calculate mean biovolume
rel_BV <- all_data%>%
  filter(species != 'Ciliate indet.')  %>%
  group_by(treatment_id, sampling)%>%
  mutate(sum = sum(cells_ml, na.rm = T)) %>%
  ungroup()%>%
  mutate(day = sampling*2,
         treatment_ID = treatment_id) %>%
  mutate(rel_V = cells_ml/sum *100) %>%
  separate(treatment_ID, into = c('treatment', 'fluctuation'), sep = '_')
rel_BV$fluctuation[is.na(rel_BV$fluctuation)] <-0


#adjust fluctuation type as factor for coloring
rel_BV$fluctuation <- factor(as.factor(rel_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
label = c('0' = "constant",'48' = "Fluctuating 48", '36'= "Fluctuating 36", '24'='Fluctuating 24', '12'='Fluctuating 12', '6'='Fluctuating 6')


#barplot
species <- ggplot(rel_BV, aes( x = day, y = rel_V))+
  geom_col(aes(fill = species))+
  scale_color_brewer()+
  scale_x_continuous(limits = c(-5,40), breaks = c(0,12, 20, 28, 36))+
  facet_wrap(~fluctuation, labeller= labeller(fluctuation = label))+
  labs(x = 'Time [days]', y = 'Rel. species abundance [%]', fill = 'Species')+
  theme(  panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border= element_rect(colour = "black", fill=NA, size=0.5),
          strip.text = element_text(face = 'bold'),
          legend.background = element_blank(),
          legend.position  ='bottom',
          legend.key = element_blank(),
          text = element_text(size=20))
species
#ggsave(plot=species, file = 'rel_ab_perspecies.png', width = 15, height = 9)
 
#new df for dominant groups
groups <- all_data%>%
  filter(species != 'Ciliate indet.')  %>%
  group_by(treatment_id, sampling)%>%
  mutate(sum = sum(cells_ml, na.rm = T)) %>%
  ungroup()%>%
  mutate(day = sampling*2,
         treatment_ID = treatment_id) %>%
  group_by(phylum) %>%
  mutate(rel_V = cells_ml/sum *100) %>%
  separate(treatment_ID, into = c('treatment', 'fluctuation'), sep = '_')
groups$fluctuation[is.na(groups$fluctuation)] <-0

#adjust fluctuation type as factor for coloring
groups$fluctuation <- factor(as.factor(groups$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
label = c('0' = "constant",'48' = "Fluctuating 48", '36'= "Fluctuating 36", '24'='Fluctuating 24', '12'='Fluctuating 12', '6'='Fluctuating 6')


#barplot
group <- ggplot(groups, aes( x = day, y = rel_V))+
  geom_col(aes(fill = phylum))+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_continuous(limits = c(-5,40), breaks = c(0,12, 20, 28, 36))+
  facet_wrap(~fluctuation, labeller= labeller(fluctuation = label))+
  labs(x = 'Time [days]', y = 'Rel. abundance [%]', fill = 'Taxonomic group')+
  theme(  panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border= element_rect(colour = "black", fill=NA, size=0.5),
          strip.text = element_text(face = 'bold'),
          legend.background = element_blank(),
          legend.position  ='bottom',
          legend.key = element_blank(),
          text = element_text(size=20))
#ggsave(plot=group, file = 'rel_group.png', width = 15, height = 9)


##############################################################################
#### Zooplankton species ####
Zoo_dominant_Groups_Jacuzzitron <- read_excel("project_data/Zoo_dominant_Groups_Jacuzzitron.xlsx", 
                                              na = "NA")
names(Zoo_dominant_Groups_Jacuzzitron)
### subdataframe with zoo informations only ####

Zoo_dominant_Groups <- Zoo_dominant_Groups_Jacuzzitron %>% 
  drop_na(copepods_l) %>%
  rename(sample = "...1" ) %>%
  rename(euplotes_l = ciliateA_l,
         strombidinopsis_l = ciliateB_l) %>%
  dplyr::select(-eggA_l, -eggs_l, -eggC_l, -eggB_l, -fuzzyball_l, -IDK_A_l, -fuzzything_l, -wormlikething_l, ) %>%
  tidyr::gather(key = 'group', value = 'abundance', -n_sampling,-fluct, -sample, -sampling,  -planktotron) %>%
  group_by(fluct, n_sampling) %>%
  mutate(sum = sum(abundance, na.rm = T)) %>%
  ungroup()%>%
  mutate(day = n_sampling*2) %>%
  group_by(group) %>%
  mutate(rel_ab = abundance/sum *100) %>%
    mutate(fluctuation = paste(ifelse(fluct == 'Con', 'constant', ifelse(fluct == 'F48', 'Fluctuating 48', ifelse(fluct == 'F36', 'Fluctuating 36', ifelse(fluct == 'F24', 'Fluctuating 24', ifelse(fluct == 'F12', 'Fluctuating 12', 'Fluctuating 6')))))))

str(Zoo_dominant_Groups)
levels(as.factor(Zoo_dominant_Groups$fluctuation))
Zoo_dominant_Groups$fluctuation <- factor(as.factor(Zoo_dominant_Groups$fluctuation),levels=c("constant", "Fluctuating 48", "Fluctuating 36", "Fluctuating 24", "Fluctuating 12" ,"Fluctuating 6" ))


ggplot(Zoo_dominant_Groups, aes( x = day, y = rel_ab))+
  geom_col(aes(fill = group))+
 # scale_x_continuous(limits = c(-5,40), breaks = c(0,12, 20, 28, 36))+
  facet_wrap(~fluctuation)+
  labs(x = 'Time [days]', y = 'Rel. abundance [%]', fill = 'Zooplankton group')+
  theme(  panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border= element_rect(colour = "black", fill=NA, size=0.5),
          strip.text = element_text(face = 'bold'),
          legend.background = element_blank(),
          legend.position  ='bottom',
          legend.key = element_blank())#
#ggsave(plot = last_plot(),'zooplanktongroups-relab.png', width = 8, height = 4)          
