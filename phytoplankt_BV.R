### Phytoplankton counts analysis ###
# Charlotte Kunze
# June 29, 2020
library(tidyverse)
library(vegan)
########################################################################################
CountsBV <-read_delim("~/Desktop/MA/MA_Rcode/project_data/CountsBV.csv", 
                      ";", escape_double = FALSE, locale = locale(decimal_mark = ","), 
                      trim_ws = TRUE)
names(CountsBV) 
str(CountsBV) 


#new df to calculate mean biovolume
rel_BV <- CountsBV%>%
  filter(species != 'Ciliate indet.') %>%
  mutate(fluctuation = paste(ifelse(treatment_id == 'control', 0, ifelse(treatment_id == 'Fluctuating_48', 48, ifelse(treatment_id == 'Fluctuating_36', 36,
                                                        ifelse(treatment_id == 'Fluctuating_24', 24, ifelse(treatment_id == 'Fluctuating_12', 12, 6))))))) %>%
  group_by(treatment_id, fluctuation, sampling)%>%
  mutate(sum = sum(volume, na.rm = T)) %>%
  ungroup()%>%
  mutate(day = sampling*2) %>%
  mutate(rel_V = volume/sum *100) %>%
  filter(rel_V >0.1)


rel_BV$fluctuation <- factor(as.factor(rel_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
label = c('0' = "constant",'48' = "Fluctuating 48", '36'= "Fluctuating 36", '24'='Fluctuating 24', '12'='Fluctuating 12', '6'='Fluctuating 6')


#plot
ggplot(rel_BV, aes( x = day, y = rel_V))+
  geom_col(aes(fill = species))+
  scale_color_brewer()+
  scale_x_continuous(limits = c(-5,40), breaks = c(0,12, 20, 28, 36))+
  facet_wrap(~fluctuation, labeller= labeller(fluctuation = label))+
  labs(x = 'Time [days]', y = 'mean contribution of species to BioVolume', fill = 'Species')+
  theme(  panel.background = element_rect(fill = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border= element_rect(colour = "black", fill=NA, size=0.5),
          #strip.background = element_rect(color ='black', fill = 'white'),
          strip.text = element_text(face = 'bold'),
          legend.background = element_blank(),
          legend.position  ='bottom',
          legend.key = element_blank(),
          text = element_text(size=12))
#ggsave(plot=last_plot(), file = 'rel_V_perspecies.tiff', width = 11, height = 8)


########################################################################################
#### indices based on BioV ####

shannon_BV <- all_data %>%
  mutate(species_id = paste(genus, species, sep = '_'))%>%
  select(-genus, -species,-id,-mean_BioV, -cells_ml) %>%
  #group_by(treatment, sampling, species_id) %>%
  # summarise(mean_V = mean(volume, na.rm = T)) %>%
  #  ungroup()%>%
  spread(key = species_id, value = volume) %>%
  group_by(treatment, sampling) %>%
  mutate(fluctuation = as.numeric(paste(ifelse(treatment == 'con', 0, 
                                               ifelse(treatment == 'F48', 48, ifelse(treatment == 'F36', 36, ifelse(treatment == 'F24', 24,
                                                                                                                    ifelse(treatment == 'F12',12,6)))))))) %>%
  select(treatment, sampling, fluctuation, everything()) 

#remove NAs and exchange with 0
shannon_BV[is.na(shannon_BV)] <- 0

##calculate shannon diversity index
shannon_BV$shan = diversity(shannon_BV[, -c(1:4)], MARGIN = 1, index='shannon') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan,everything() )
shannon_BV$simpson = diversity(shannon_BV[, -c(1:5)], MARGIN = 1, index='invsimpson') #new column containing the calculated shannon index
shannon_BV <- select(shannon_BV, treatment, sampling, fluctuation, MC, shan, simpson,everything() )
## calculate species richness
absence_presence <- decostand(shannon_BV[, -c(1:6)], method= 'pa', na.rm=T) #df giving absence/presence data using decostand function
shannon_BV$no = apply(absence_presence, MARGIN = 1, FUN = sum) #new column containing the sum of species present per side (by row = MARGIN = 1)

shannon_BV$evenness = shannon_BV$shan/log(shannon_BV$no)

diversity_BV <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  gather(key = 'index', value = 'value', -sampling, -fluctuation, -MC) %>%
  group_by(fluctuation, sampling, index) %>%
  summarise(mean_index = mean(value, na.rm = T),
            sd = sd(value, na.rm =T),
            se = sd/ sqrt(n()),
            n = n()) %>%
  mutate(lower.ci.mpg = mean_index - 1.96*se/sqrt(n),
         upper.ci.mpg = mean_index + 1.96*se/sqrt(n),
         day = sampling *2)

diversity_BV$fluctuation <- factor(as.factor(diversity_BV$fluctuation),levels=c("0", "48", "36", '24', '12', '6'))
even_p <- ggplot(subset(diversity_BV, index == 'evenness'), aes(x = day, y = mean_index, group = fluctuation)) +
  geom_line(linetype = 'dashed', size = 0.5, aes(col = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci.mpg , ymax = upper.ci.mpg, col = fluctuation), width = .5, position = position_dodge(width = .9))+
  geom_point(aes(fill = fluctuation), pch = 21,size= 3,col = '#000000', position = position_dodge(width = .9))+
  # scale_x_continuous(limits = c(-1, 20), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = '', y= 'Evenness species composition', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))
even_p

even_no <- ggplot(subset(diversity_BV, index == 'no'), aes(x = day, y = mean_index, group = fluctuation)) +
  geom_line(linetype = 'dashed', size = 0.5, aes(col = fluctuation))+
  geom_errorbar(aes(ymin = lower.ci.mpg , ymax = upper.ci.mpg, col = fluctuation), width = .5, position = position_dodge(width = .9))+
  geom_point(aes(fill = fluctuation), pch = 21,size= 3,col = '#000000', position = position_dodge(width = .9))+
  # scale_x_continuous(limits = c(-1, 20), breaks = seq(0,18,2))+
  scale_fill_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  scale_color_manual(values = c( '#000000','#0868ac','#41b6c4','#31a354','#addd8e','#fed976'))+
  labs(x = 'Time [days]', y= 'Species richness', fill = 'Fluctuation  \nfrequency [h]', col = 'Fluctuation  \nfrequency [h]')+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),legend.position = 'none', 
        panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
        #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
        panel.border= element_rect(colour = "black", fill=NA, size=1),
        strip.text = element_text(face = 'bold'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        text = element_text(size=17))
even_no

library(cowplot)
plot_grid(even_p, even_no,  labels=c("(a)","(b)"),ncol = 2, label_size = 18, hjust = 0, vjust = 0.95)
#ggsave(plot = last_plot(), file = 'even_no_specBV.tiff', width = 9, height = 5)



#ggsave(plot = last_plot(), file = 'shannon_div_BioV_time.png')  
#####################################################################################

#### STATS ####

diversity <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  mutate(day = 2*sampling,
         dayname = as.factor(day)) %>%
  mutate(fluctuation = as.numeric(fluctuation),
         interval = 48/fluctuation)
diversity$interval[!is.finite(diversity$interval)] <- 0 
even1 <- lmer(no ~ interval*day + (1|MC) + (1|dayname), data=diversity)
summary(even1)
anova(even1)  