#Clear all
rm(list=ls(all=TRUE))
dev.off()
#https://jkzorz.github.io/2019/06/06/NMDS.html
##### Script for ndms plots 
# Charlotte Kunze
##-----------------------------------------------------------------------------------------------##
#check working directory
getwd()
##-----------------------------------------------------------------------------------------------##
# Load required packages
library(tidyverse)
library(lubridate)
library(scales)
library(vegan)

####import data####

####PCA on pigment data only ####
pigments <- read.csv2("~/Desktop/MA/MA_Rcode/project_data/master_pigments.csv", sep = ";", dec = ',')
str(pigments)

#bring data in form
pigment_datas <- pigments %>%
  gather(key = 'pigments',value = 'value',-X, -no, -date,-treatment, -sampling, -planktotron )%>%
  group_by(treatment,sampling, planktotron)%>%
  mutate(sum = sum(value)) %>%
  ungroup()%>%
  spread(pigments, value ) %>%
  mutate(rel_allo = Allo/sum)%>%
  mutate(rel_anth = Anth/sum)%>% ###berechnung der relativen werte in neuer spalte nach Schlueter et al.
  mutate(rel_bb.Car = bb.Car/sum)%>%
  mutate(rel_cNeo = c.Neo/sum) %>%
  mutate(rel_cantha = Cantha/sum)%>%
  mutate(rel_chl.a = Chl.a/sum)%>%
  mutate(rel_chl.b = Chl.b/sum)%>%
  mutate(rel_chl.c1 = Chl.c1/sum)%>%
  mutate(rel_chl.c2 = Chl.c2/sum)%>%
  mutate(rel_diad = Diadino/sum)%>%
  mutate(rel_diato = Dino/sum)%>%
  mutate(rel_dino = Diato/sum)%>%
  mutate(rel_echin = Echin/sum)%>%
  mutate(rel_fuco = Fuco/sum)%>%
  mutate(rel_lut = Lut/sum)%>%
  mutate(rel_myxo = Myxo/sum)%>%
  mutate(rel_peri = Peri/sum)%>%
  mutate(rel_phe.a = Phe.a/sum)%>%
  mutate(rel_viola = Viola/sum) %>%
  mutate(rel_zea = Zea/sum)%>%
  dplyr::select(-c("X","Allo", "Anth","bb.Car", "c.Neo", "Cantha","Chl.a","Chl.c1",
                   "Chl.c2","Cryp","Diadino", "Diato","Dino","Echin","Lut","Myxo","Peri",       
                   "Phe.a","Phe.b", "Viola", "Zea" , 'Fuco', 'Chl.b')) %>%
  gather(key = 'pigment', value = 'value', -planktotron, - sampling, -treatment,-sum, -no, -date) %>%
  group_by(treatment, sampling, pigment) %>% 
  summarise(mean = mean(value, na.rm = T))%>%
  mutate(fluctuation = paste(ifelse(treatment == 'control', 0, ifelse(treatment == 'Fluctuating_48', 48, ifelse(treatment == 'Fluctuating_36', 36,
                                                                                                                ifelse(treatment == 'Fluctuating_24', 24, ifelse(treatment == 'Fluctuating_12', 12, 6)))))))
pigment_datas$fluctuation <- factor(as.factor(pigment_datas$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
pigment_datas$sampling = as.numeric(pigment_datas$sampling)
pigment_datas[is.na(pigment_datas)] <-0 #make sure we don't have NAs

names(pigment_datas)
# 1. decide which variables are active, remove character variables and turn them into row names 
pigments_dat <- pigment_datas %>% 
  ungroup()%>%
  dplyr::select(-treatment) %>% #remove columns which are not needed anymore
  spread(key = pigment, value = mean) %>%
  dplyr::select(-rel_bb.Car) #remove variable with only 0 entries

#fill all N.A.s with 0
pigments_dat[3:ncol(pigments_dat)][is.na(pigments_dat[3:ncol(pigments_dat)])] <-0          
str(pigments_dat)


com_p = pigments_dat[, 3:ncol(pigments_dat)]   #create a data frame from column 4 till 18 (all the species)
m_com_p = as.matrix(com_p)             #the data frame we just created should be seen as matrix

#nmds
nmds_p = metaMDS(m_com_p, k=2)         #k are the dimensions for the plots
plot(nmds_p)                         #not a nice plot, so we continue in order to plot it with ggplot
#within the plot red crosses = species, open circles = communities

#extract nmds scores
data.scores_p = as.data.frame(scores(nmds_p)) #obtain the coordinates for the axes and put them in a new data frame

#add columns with further information
data.scores_p$sampling = pigments_dat$sampling
data.scores_p$day = pigments_dat$sampling *2
#data.scores$date = data_nd$date
data.scores_p$treat = pigments_dat$fluctuation


head(data.scores_p)  #see our results (but only the first rows)

#plot with ggplot
nmds_plot <-ggplot(subset(data.scores_p, sampling %in% c(0)), 
                   aes(x = NMDS1, y = NMDS2, alpha = day, size = day))+           #treatments in same colour
  geom_point(aes(fill = treat), size =4, pch = 21, col = 'black')+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_x_continuous(limits = c(-0.35, 0.35), breaks = seq(-0.3,0.3,0.15))+
  scale_y_continuous(limits = c(-0.35, 0.35), breaks = seq(-0.3,0.3,0.15))+
  annotate('text', x = 0.0, y = 0.38, label = paste('Stress:', round(nmds_p$stress, 3)))+   #3 gives number of digits
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=18))
nmds_plot
ggsave(plot = nmds_plot, file = 'NMDS_pigmente_red0.tiff', width = 6, height = 5)  
## stress
# stress < 0.05 excellent representation
# stress < 0.1 great
# stress < 0.2 okay
# stress < 0.3 poor represented
# stress > 0.3 choose another method


################################################################################
####################################################################################

#### composition data ####
CountsBV <-read_delim("~/Desktop/MA/MA_Rcode/project_data/CountsBV.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)
names(CountsBV) 
str(CountsBV) 

#new df to calculate mean biovolume
nmds_BV <- CountsBV%>%
  group_by(sampling, treatment_id, MC) %>%
  mutate(sum = sum(volume, na.rm  =T ),
         rel = volume/sum)%>%
  group_by(treatment_id, sampling, species, MC) %>%
  summarise(mean_V = mean(rel, na.rm = T)) %>% #calculate mean cells with biovolume
  drop_na(mean_V) %>%
  separate(treatment_id, into = c('treat', 'fluctuation'), '_') %>%
  dplyr::select(-treat) %>%
  spread(key = species, value = mean_V) 
nmds_BV$fluctuation[is.na(nmds_BV$fluctuation)] <-0
#fill all N.A.s with 0
nmds_BV[4:ncol(nmds_BV)][is.na(nmds_BV[4:ncol(nmds_BV)])] <-0          
str(nmds_BV)


com1 = nmds_BV[, 4:ncol(nmds_BV)]   #create a data frame from column 4 till 18 (all the species)
m_com1 = as.matrix(com1)             #the data frame we just created should be seen as matrix

#nmds
nmds1 = metaMDS(m_com1, k=2)         #k are the dimensions for the plots
plot(nmds1)                         #not a nice plot, so we continue in order to plot it with ggplot
#within the plot red crosses = species, open circles = communities

#extract nmds scores
data.scores1 = as.data.frame(scores(nmds1)) #obtain the coordinates for the axes and put them in a new data frame

#add columns with further information
data.scores1$MC = nmds_BV$MC
data.scores1$sampling = nmds_BV$sampling
data.scores1$day = nmds_BV$sampling*2
#data.scores$date = data_nd$date
data.scores1$treatment = nmds_BV$fluctuation


head(data.scores1)  #see our results (but only the first rows)

#plot with ggplot
data.scores1$treatment <- factor(data.scores1$treatment, levels= c('0', '48', '36', '24', '12', '6'))

ggplot(subset(data.scores1, sampling == 18), aes(x = NMDS1, y = NMDS2, colour = treatment))+           #treatments in same colour
  geom_point(size = 4, aes(color = treatment), size = 4)+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  #scale_shape_manual(values = c(15,10,17))+
  labs(alpha = 'Time [days]', color = 'Fluctuation \nfrequency [h]')+  
  scale_x_continuous(limits = c(-0.5, 0.75), breaks = seq(-0.5, 0.5, 0.25))+
  scale_y_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.5, 0.5, 0.25))+
  annotate('text', x = 0.2, y = 0.48, label = paste('Stress:', round(nmds1$stress, 3)))+   #3 gives number of digits
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=1),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='right',
         legend.key = element_blank(),
         text = element_text(size=18))        
#ggsave(plot = last_plot(), file = 'NMDS_BV1.tiff', width = 8, height = 6)  
## stress
# stress < 0.05 excellent representation
# stress < 0.1 great
# stress < 0.2 okay
# stress < 0.3 poor represented
# stress > 0.3 choose another method

