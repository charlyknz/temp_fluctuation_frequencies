## Model formulation for thesis analysis
rm(list=ls()) #Empty the environment

# scatter plot with bray-curits diss and LRR of function (carbon)
# By Charlotte Kunze 06.08.2020

#### Load packages & import data ####
library(vegan)
library(tidyverse)
library(lme4)
library(nlme)
library(tidyverse)
library(readxl)
library(rr2)

#### Zooplankton data ####
Mastertable <- read_excel("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.xlsx")
str(Mastertable)

## merge with treatment information

master <- Mastertable %>%
  rename(C_Zoo_umol_l=C_Zoo_µmol_l,
         MC = planktotron)%>%
  dplyr::select(MC, sampling, fluctuation, C_Zoo_umol_l ) %>%
  drop_na(C_Zoo_umol_l) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day),
         MC = as.character(MC)) %>%
  mutate(interval = 48/fluctuation) %>%
  mutate(log = log(C_Zoo_umol_l))

master$interval[!is.finite(master$interval)] <- 0
str(master)

ggplot(master,aes(day, log,colour = interval))+geom_point(size=1.5)
str(master)

##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
M1 = lme(log ~ interval*day, random = ~1|MC, method = 'REML', data = master)
M2 = lme(log ~ interval*day, random = ~0+day|MC, method = 'REML', data = master)
M3 = lme(log ~ interval*day, random = ~day|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = master)


##model using dayname and MC as random effects

#compare models
anova(M1,M2, M3)
 
#M1 wins AIC 197.82

#save fitted values
master$fit_InterceptOnly2 <- predict(M1)


# mixed model vs. model without random component: gls 
M0=gls(log~interval*day, method="REML",data =master, na.action=na.omit)
anova(M1, M0) #gls is better AIC 195.86

#Autocorrelation test (data are not independent) - only if do not have nas
plot(ACF(M0), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(M0, type = "normalized"), ylab="residuales")
hist(resid(M0, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(M0),resid(M0, type = "normalized"),ylab="residuales")
qqnorm(resid(M0, type = "normalized"), main=""); qqline(resid(M0, type = "normalized"))


M5=lme(log ~ interval*day, random=~0+day|MC, method="REML",
         corr=corAR1(0.8,form=~day|MC),na.action=na.omit,data=master)
anova(M5, M0) #model with random effect is better


# add weights
ctrl <- lmeControl(opt='optim');

vf1 <-varIdent(form = ~ 1|MC)
M6=lme(log ~ interval*day, random=~0+day|MC, method="REML", 
       corr=corAR1(0.8,form=~day|MC),weights = vf1, na.action=na.omit,data=master)
anova(M5, M6)

vf3 <- varIdent(form =~ 1|day)
M7=lme(log ~ interval*day, random=~0+day|MC, method="REML", control = ctrl,
       corr=corAR1(0.8,form=~day|MC),weights = vf3, na.action=na.omit,data=master)

vf4 <- varExp(form =~ interval)
M8=lme(log ~ interval*day, random=~0+day|MC, method="REML", control = ctrl,
       corr=corAR1(0.8,form=~day|MC),weights = vf4, na.action=na.omit,data=master)
str(master)
anova(M6, M8, M7, M5)
summary(M7)
anova(M7)

R2.lik(M0)

#Autocorrelation test (data are not independent) - only if do not have nas
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(M7, type = "normalized"), ylab="residuales")
plot(fitted(M7),resid(M6, type = "normalized"),ylab="residuales")
hist(resid(M7, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
qqnorm(resid(M7, type = "normalized"), main=""); qqline(resid(M7, type = "normalized"))
summary(M7)
anova(M7)


#### Phytoplankton Diversity data ####

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


#### Compositional turnover ####
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
  group_by(Fluctuation, sampling2, day) %>%
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
  labs(col =  'Fluctuation frequency (in h)',x = 'Time [days]', y = 'Bray-Curtis distance')+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'italic'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=12))
ggsave(plot = last_plot(), file = 'BrayDistance_overTimePhytoCoutns.png', width = 7, height = 4)
hist((data.dist1$log))
qqnorm(data.dist1$log)
qqline(data.dist1$log)


##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(log ~ interval*day, random = ~1|MC, method = 'REML', data = data.dist1)
C_m2 = lme(log ~ interval*day, random = ~0+day|MC, method = 'REML', data = data.dist1)
C_m3 = lme(log ~ interval*day, random = ~day|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = data.dist1)


##model using dayname and MC as random effects

#compare models
anova(C_m1,C_m2, C_m3)

#save fitted values
data.dist1$fit_InterceptOnly2 <- predict(C_m2)


## fit model output ##
data.dist1$Fluctuation <- factor(data.dist1$Fluctuation, levels = c('0', '48','36','24','12','6'))
ggplot(data.dist1, aes(x = day, color = Fluctuation, y=log ))+
  geom_point()+
  geom_line(aes(y = fit_InterceptOnly2), size = 1) +
  #scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2))+
  theme_classic()


# mixed model vs. model without random component: gls 
C_m0=gls(log~interval*day, method="REML",data =data.dist1, na.action=na.omit)
anova(C_m2, C_m0) #gls is better

anova(C_m2)
#Autocorrelation test (data are not independent) - only if do not have nas
plot(ACF(C_m0), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(C_m0, type = "normalized"), ylab="residuales")
hist(resid(C_m0, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(C_m0),resid(C_m2, type = "normalized"),ylab="residuales")
qqnorm(resid(C_m0, type = "normalized"), main=""); qqline(resid(C_m0, type = "normalized"))

# gls with weighted variance 
M1<-gls(log ~ interval*day, data = data.dist1)

#Fixed variance
vf1Fixed <- varFixed(~interval)
M1.1<-gls(log ~ interval*day, data = data.dist1, weights = vf1Fixed)

#note: propably the VarIdent would be a good choice, however here we see the lm is the better choice
#Power of the variance covariate
vf3 <- varPower(form = ~interval)
M3<-gls(log ~ interval*day, data = data.dist1, weights = vf3)

vf4 <- varPower(form = ~ interval | day)
M4<-gls(log ~ interval*day, data = data.dist1, weights = vf4)

#If the variance covariate can take the value of zero
## the exponential variance structure:
vf5 <- varExp(form = ~interval)
M5<-gls(log ~ interval*day, data = data.dist1, weights = vf5)

#The varConstPower Variance Structure
#The varPower should not be used if the variance covariate takes the value of zero.
vf6 <- varConstPower(form = ~interval)
M6<-gls(log ~ interval*day,  weights = vf6,data = data.dist1)

vf7 <- varConstPower(form = ~interval| day)
M7<-gls(log ~ interval*day,  weights = vf7,data = data.dist1)

#Different variances per stratum
vf2 <- varIdent(form = ~ 1|day)
M2<-gls(log ~ interval*day, data = data.dist1, weights = vf2)
anova(M1, M2, M5) #compare models

#final model:
M.dist.lm<-lm(log ~ interval*day, data = data.dist1)
anova(M.dist.lm) # indicated time * richness effect 

R2.lik(M.dist.lm)
# plot
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(M.dist.lm, which = c(1), col = 1, add.smooth = FALSE,caption = "")
plot(diversity$day, resid(M.dist.lm), xlab = "Month",ylab = "Residuals")
plot(diversity$interval, resid(M.dist.lm), xlab = "DML",ylab = "Residuals")
par(op)

#find R2
library(rr2)
R2.lik(M.lm)

#### #######################################################################

#### Diversity indices calculation ####

shannon_BV <- all_data %>%
select(-phylum, -grid_length, -magnification, -'stripes/grids',-date, -volume, -counts) %>%
  spread(key = species, value = cells_ml) %>%
  group_by(treatment_id, sampling) %>%
  separate(treatment_id, into = c('treatment', 'fluctuation'), '_')%>%  
  select(treatment, sampling, fluctuation, everything()) 

#remove NAs and exchange with 0
shannon_BV$fluctuation[is.na(shannon_BV$fluctuation)] <- 0
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
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) 


# Data Wrangling for LMM #

diversity <- shannon_BV %>%
  ungroup() %>%
  select(MC, fluctuation, sampling, evenness, no, shan,simpson) %>%
  mutate(day = 2*sampling,
         dayname = as.factor(day)) %>%
  mutate(fluctuation = as.numeric(fluctuation),
         interval = 48/fluctuation)%>%
  mutate(log = log(no))
diversity$interval[!is.finite(diversity$interval)] <- 0 


#### Richness ####
hist(log(diversity$no))
qqnorm(diversity$no)
qqline(diversity$no)


##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(log ~ interval*day, random = ~1|MC, method = 'REML', data = diversity)
C_m2 = lme(log ~ interval*day, random = ~0+day|MC, method = 'REML', data = diversity)
C_m3 = lme(log ~ interval*day, random = ~day|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = diversity)


##model using dayname and MC as random effects

#compare models
anova(C_m1,C_m2, C_m3)


# mixed model vs. model without random component: gls 
C_m0=gls(log~interval*day, method="REML",data =diversity, na.action=na.omit)
anova(C_m2, C_m0) #model with random effect is better

C_m5=lme(log ~ interval*day, random=~0+day|MC, method="REML",
         corr=corAR1(0.6,form=~day|MC),na.action=na.omit,data=diversity)
anova(C_m5, C_m0) #model with random effect is better

#Autocorrelation test (data are not independent) - only if do not have nas
plot(ACF(C_m0), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(C_m5, type = "normalized"), ylab="residuales")
hist(resid(C_m5, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(C_m5),resid(C_m2, type = "normalized"),ylab="residuales")
qqnorm(resid(C_m5, type = "normalized"), main=""); qqline(resid(C_m5, type = "normalized"))
anova(C_m0)
anova(C_m5)

# gls with weighted variance 
#Fixed variance
vf1Fixed <- varFixed(~interval)
M1.1<-gls(no ~ interval*day, data = diversity, weights = vf1Fixed)

#note: propably the VarIdent would be a good choice, however here we see the lm is the better choice
#Power of the variance covariate
vf3 <- varPower(form = ~interval)
M3<-gls(no ~ interval*day, data = diversity, weights = vf3)

vf4 <- varPower(form = ~ interval | day)
M4<-gls(no ~ interval*day, data = diversity, weights = vf4)

#If the variance covariate can take the value of zero
## the exponential variance structure:
vf5 <- varExp(form = ~interval)
M5<-gls(no ~ interval*day, data = diversity, weights = vf5)

#The varConstPower Variance Structure
#The varPower should not be used if the variance covariate takes the value of zero.
vf6 <- varConstPower(form = ~interval)
M6<-gls(no ~ interval*day,  weights = vf6,data = diversity)

vf7 <- varConstPower(form = ~interval| day)
M7<-gls(no ~ interval*day,  weights = vf7,data = diversity)

#Different variances per stratum
vf2 <- varIdent(form = ~ 1|interval)
M2<-gls(no ~ interval*day, data = diversity, weights = vf2)
anova(C_m0, M2) #compare models

#final model:
M.lm<-lm(log ~ interval*day, data = diversity)
anova(M.lm) # significant time effect of species richness

#find R2
library(rr2)
R2.lik(M.lm)

op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(M.lm, which = c(1), col = 1, add.smooth = FALSE,caption = "")
plot(diversity$day, resid(M.lm), xlab = "Month",
     ylab = "Residuals")
plot(diversity$interval, resid(M.lm), xlab = "DML",
     ylab = "Residuals")
par(op)
#### inverse Simpson ####
diversity$logSimpson = log(diversity$simpson)
hist(log(diversity$logSimpson))
qqnorm(diversity$logSimpson)
qqline(diversity$logSimpson)


##### fixed and random effects ######
# Response variable: Carbon mg (measured over time)
# Explanatory variables: treatments (fluctuation hours) + time
# Random: MC number

#Step 1. Test possible cases:Fit models using REML (use ML in simplifications)
C_m1 = lme(logSimpson ~ interval*day, random = ~1|MC, method = 'REML', data = diversity)
C_m2 = lme(logSimpson ~ interval*day, random = ~0+day|MC, method = 'REML', data = diversity)
C_m3 = lme(logSimpson ~ interval*day, random = ~day|MC, method = 'REML',control= lmeControl(niterEM =5000, msMaxIter =5000, msMaxEval =5000), data = diversity)


##model using dayname and MC as random effects

#compare models
anova(C_m1,C_m2, C_m3)

# mixed model vs. model without random component: gls 
C_m0=gls(logSimpson~interval*day, method="REML",data =diversity, na.action=na.omit)
anova(C_m2, C_m0) #model with random effect is better


#Autocorrelation test (data are not independent) - only if do not have nas
plot(ACF(C_m2), alpha=0.05)

#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(C_m0, type = "normalized"), ylab="residuales")
hist(resid(C_m0, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(C_m0),resid(C_m2, type = "normalized"),ylab="residuales")
qqnorm(resid(C_m0, type = "normalized"), main=""); qqline(resid(C_m0, type = "normalized"))

#compare wiath model with autocorrelation
C_m5=lme(logSimpson ~ interval*day, random=~0+day|MC, method="REML",
         corr=corAR1(0.6,form=~day|MC),na.action=na.omit,data=diversity)
anova(C_m5, C_m0) #model without random effect is better


# linear regression
M1<-gls(logSimpson~interval*day, data =diversity)
anova(M1)

op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
plot(M1, which = c(1), col = 1, add.smooth = FALSE,caption = "")
plot(diversity$day, resid(M1), xlab = "Month",
       ylab = "Residuals")
plot(diversity$interval, resid(M1), xlab = "DML",
       ylab = "Residuals")
par(op)

# gls with weighted variance 
#Fixed variance
vf1Fixed <- varFixed(~interval)
M1.1<-gls(logSimpson ~ interval*day, data = diversity, weights = vf1Fixed)

#note: propably the VarIdent would be a good choice, however here we see the lm is the better choice
#Power of the variance covariate
vf3 <- varPower(form = ~interval)
M3<-gls(logSimpson ~ interval*day, data = diversity, weights = vf3)

vf4 <- varPower(form = ~ interval | day)
M4<-gls(logSimpson ~ interval*day, data = diversity, weights = vf4)

#If the variance covariate can take the value of zero
## the exponential variance structure:
vf5 <- varExp(form = ~interval)
M5<-gls(logSimpson ~ interval*day, data = diversity, weights = vf5)

#The varConstPower Variance Structure
#The varPower should not be used if the variance covariate takes the value of zero.
vf6 <- varConstPower(form = ~interval)
M6<-gls(logSimpson ~ interval*day,  weights = vf6,data = diversity)

vf7 <- varConstPower(form = ~interval| day)
M7<-gls(logSimpson ~ interval*day,  weights = vf7,data = diversity)

#Different variances per stratum
vf2 <- varIdent(form = ~ 1|interval)
M2<-gls(logSimpson ~ interval*day, data = diversity, weights = vf2)
anova(M1, M2, M5) #compare models


#check residuals for other model
#Residuals
par(mfrow=c(2,2),cex.axis=1.2, cex.lab=1.5)
plot(resid(M2, type = "normalized"), ylab="residuales")
hist(resid(M2, type = "normalized"), ylab="frecuencia",xlab="residuales", main="")
plot(fitted(M2),resid(C_m2, type = "normalized"),ylab="residuales")
qqnorm(resid(M2, type = "normalized"), main=""); qqline(resid(M2, type = "normalized"))

#final model after AIC:
M1.lm<-lm(logSimpson~interval*day, data =diversity)
anova(M1.lm) #no sign. effects
R2.lik(M1.lm)

#### ANOSIM AND SIMPER TEST ####
#https://jkzorz.github.io/2019/06/11/ANOSIM-test.html

# install packages
library(tidyverse)
library(vegan)

### Start ####
# first import your data and bring them in the correct format
# to implement them in your anosim, we need a matrix with our values (species/pigment names) as headers
# second we need a dataframe, which still contains all our explanatory variables such as sampling, treatment usw


#### Step 1: import datasets####
data <- read_delim("~/Desktop/MA/MA_Rcode/project_data/phyto_new.csv", 
                   ";", escape_double = FALSE, col_types = cols(date = col_character()), 
                   locale = locale(decimal_mark = ","), 
                   trim_ws = TRUE) %>%
  drop_na(MC)

#check your importt:
View(data)
str(data)

#change variable format to numeric 
data$cells_ml <- as.numeric(gsub(",", ".", data$cells_ml))
## merge with treatment information
treatments <- read.csv2('~/Desktop/MA/MA_Rcode/project_data/treatments_units.csv')
str(treatments)
names(treatments) = c('MC', 'fluctuation', 'treatment')
rich <- left_join(data, treatments, by = c('MC')) %>%
  separate(species, into = c('genus', 'species'), ' ')


##### Step 2: get data in the right format ####
# create matrix and dataframe with necessary informations

rel_BV <- rich%>%
  #filter(sampling ==6)%>%
  mutate(species_id = paste(genus, species, sep = '_'))%>% #create a species id column
  group_by(treatment, fluctuation, sampling, MC)%>% #groups after MC, treatment, sampling etc
  mutate(sum = sum(cells_ml, na.rm = T), #calculate sum and relative contribution of each species/pigment
         rel_V = cells_ml/sum) %>%
  ungroup()%>%
  dplyr::select(sampling, fluctuation, MC, species_id,rel_V) %>% #select only important columns
  drop_na(rel_V) %>% #remove NAs from the dataset
  mutate(interval = 48/fluctuation, #interval and data columns as explanatory variables
         day = sampling *2) %>% 
  spread(key = species_id, value = rel_V) #wide format

rel_BV[is.na(rel_BV)] <-0 #exchange NA with 0
rel_BV$interval[is.infinite(rel_BV$interval)] <-0 #infinite values shall be 0 for interval

#create my data Matrix
Data <- dplyr::select(rel_BV, -sampling,-day, -fluctuation, -MC, -interval)#remove grouping variables
Data <- as.matrix(Data) #create matrix

####Step 3: ANOSIM####
anosim(Data, rel_BV$day, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

### interpretation: The divisor is chosen so that R will be in the interval -1 … +1, value 0 indicating completely random grouping.
#An R value close to “1.0” suggests dissimilarity between groups while an R value close to “0” suggests an even distribution of high and low ranks within and between groups”
## R signif greater than 0.05, means that there is no statistical difference between the microbial communities in your groups.

# treatment as grouping
anosim(Data, rel_BV$fluctuation, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

# time as grouping
anosim(Data, rel_BV$day, permutations = 999, distance = "bray", strata = NULL,
       parallel = getOption("mc.cores"))

#An R value close to  “0” suggests an even distribution of high and low ranks within and between groups”
#My significance value is much lower than 0.05
#Therefore, there is a statistically significant difference in my microbial communities based on the grouping “Time”.

simper(Data, rel_BV$day, permutations = 0, trace = FALSE,
       parallel = getOption("mc.cores"))
# greater than 70% is required to say that groups are different from each other



#### Phytoplankton POC GAMM ####
## Model formulation using GAMM

# 1. load required packages
library(lme4)
library(nlme)
library(lme)
library(mgcv)
library(tidyverse)
library(devtools)
library(itsadug)
library(rr2)
library(readxl)
library(tidymv)

## import Mastertable
Mastertable <- read_excel("~/Desktop/MA/MA_Rcode/project_data/Mastertable_fluctron.xlsx")
str(Mastertable)
names(Mastertable)

## merge with treatment information
data1 <- Mastertable%>%
  rename(C_Zoo_umol_l=C_Zoo_µmol_l,
         MC = planktotron)%>%
drop_na(carbon_umol_l) %>%
 drop_na(C_Zoo_umol_l) %>%
  dplyr::select(MC, sampling, fluctuation, C_Zoo_umol_l, carbon_umol_l) %>%
  mutate(carbo = carbon_umol_l-C_Zoo_umol_l) %>%
  mutate(day = 2* sampling, 
         dayname = as.factor(day),
         MC = as.character(MC)) %>%
  mutate(interval = 48/fluctuation)

data1$interval[!is.finite(data1$interval)] <- 0 
str(data1)


### Model 1 total carbon 

#Additive mixed model 
M_total <- gamm(carbon_umol_l ~ interval + s(day), random = list(MC =~ 1), na.action=na.omit, method = "REML", data = data1)

summary(M_total $gam) #This gives detailed output on the smoother and parametric
#terms in the models.

anova(M_total $gam) #This command gives a more compact presentation of the
#results. The anova table is not doing sequential testing!
plot(M_total $gam) #This command plots the smoothers.
plot(M_total $lme) #This command plots the normalised residuals versus fitted
#values and can be used to assess homogeneity.
summary(M_total $lme) #Detailed output on the estimated variances. 


#plot residuals
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
plot(M_total $lme)

#show residuals
data1$res.M_total  <- residuals(M_total $lme, type = "pearson")
data1$fit.M_total  <- predict(M_total $gam, type = "link")
data1$fit2.M_total <- predict(M_total $gam, type = "response")
plot(data1$fit.M_total , data1$res.M_total )

#autocorrelation
plot(Variogram(M_total$lme, robust = TRUE, data = data1, form = ~day | MC))

#plot model and data
ggplot(data1, aes(x = day, y = carbon_umol_l, color = interval, group = interval))+
  geom_point()+
  geom_line(aes(x = day, y = fit.M_total, color = interval), linetype = 'dashed')+
  theme_bw()



### Model 2 corrected phytoplankton carbon 
#Additive mixed model 
M_phyto <- gamm(carbo ~ interval + s(day), random = list(MC =~ 1), na.action=na.omit, method = "REML", data = data1)

summary(M_phyto$gam) #This gives detailed output on the smoother and parametric
#terms in the models.

anova(M_phyto$gam) #This command gives a more compact presentation of the
#results. The anova table is not doing sequential testing!
plot(M_phyto$gam) #This command plots the smoothers.
plot(M_phyto$lme) #This command plots the normalised residuals versus fitted
#values and can be used to assess homogeneity.
summary(M_phyto$lme) #Detailed output on the estimated variances. 
#Not everything is relevant.

#plots
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
plot(M_phyto$lme)

# check residuals
data1$res.M_phyto <- residuals(M_phyto$lme, type = "pearson")
data1$fit.M_phyto <- predict(M_phyto$gam, type = "link")
data1$fit2.M_phyto <- predict(M_phyto$gam, type = "response")
plot(data1$fit.M_phyto, data1$res.M_phyto)

#autocorrelation
plot(Variogram(M_phyto$lme, robust = TRUE, data = data1, form = ~day | MC))

#plot model
ggplot(data1, aes(x = day, y = carbo, color = interval, group = interval))+
  geom_point()+
  geom_line(aes(x = day, y = fit.M_phyto, color = interval), linetype = 'dashed')+
  theme_bw()


