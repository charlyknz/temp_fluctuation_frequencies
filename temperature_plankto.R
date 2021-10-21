## R Script to gather temperature data with actual chlorophyll a data


#### required packages ####
library(tidyverse)
library(readxl)
library(lubridate)
library(hms)

#----------------------------------------------------------------------------------------------------------------#

#sets working directory where your R script is saved in
setwd("~/Desktop/MA/MA_data/abiotics/temperature_pH/Data")

## import excel files using a Loop

file.list <- list.files(pattern='*.xls')
df.list <- lapply(file.list, read_excel)
df <- bind_rows(df.list, .id = "id")
str(df)
#write.csv(x = df, file = 'all_temperatures.csv')
#df$actual_temptop <- as.numeric(gsub(",", ".", df$actual_temptop))
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##
#### START ####
all_temperatures <- read_csv("~/Desktop/MA/MA_Rcode/project_data/all_temperatures.csv", 
                             locale = locale())
View(all_temperatures)
## data wrangling

temp <- all_temperatures %>%
  filter(timestamp >= '2019-08-27 08:00:00') %>%
  select(unit, timestamp, contains('actual')) %>%
  mutate(datetime = timestamp) %>%
  separate(timestamp, into = c('date', 'time'), sep = ' ') %>%
  mutate(time1 = as.character(time)) %>%
  filter(stringr::str_detect(time1, c(':00:1', ':30:1', ':15:1'))) %>% #keepsonly entries with 00:1 or 30:1 minute/seconds
  filter(datetime > "2019-08-27 15:30:15", datetime < "2019-10-02 12:00:15") #%>% #datetime between two dates
#filter(datetime > "2019-08-27 15:30:15", datetime < "2019-08-30 18:00:15") #%>% #datetime between two dates
#filter(unit %in% c(1,2,5,6,8,11)) #keep only one replicate
  
dev.off()
# order MC after treatments
temp$unit = factor(temp$unit, levels = c(1,4,6,7,2,10,9,8, 11,12,5))

#for only one replicate: 
#temp$unit = factor(temp$unit, levels = c(1,6,2,8, 11,5))

label_temp <- c('1' = 'constant','4' = 'constant', '2' = 'Fluctuating 36 h','10' = 'Fluctuating 36 h', '5' = 'Fluctuating 6h', '3' = 'Fluctuating 6h',
                '6' = 'Fluctuating 48 h','7' = 'Fluctuating 48 h','8' = 'Fluctuating 24 h','9' = 'Fluctuating 24 h', '11' = 'Fluctuating 12 h','12' = 'Fluctuating 12 h')
plot <- ggplot(temp, aes(x = datetime, y = actual_tempmiddle))+
  #geom_hline(aes(yintercept = 18), col = 'darkgrey', linetype = 'dashed', size = 0.5)+
  geom_line(size = 1.0)+
  facet_wrap(~unit, labeller = labeller(unit = label_temp), ncol = 1)+
  labs( x = 'Date', y = 'Temperature (in °C)')+
  scale_y_continuous(limits = c(14, 22), breaks = c(15,18,21))+
  scale_x_datetime(breaks = '5 days')+
  theme_classic()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background = element_blank(),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         #strip.background = element_rect(color ='black', fill = 'white'),
         strip.text = element_text(face = 'bold'),
         legend.background = element_blank(),
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size=27))
plot
ggsave(plot = plot, file = 'Fig5STemp_curves.png', width = 20, height = 30)
## ------------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------------ ##

#new data frame with selected columns
temp_plankto <- tbl %>%
  select(timestamp, unit, actual_temptop) %>%
  mutate(date = timestamp) %>%
  separate(col = date, into = c('date', 'time'), ' ', fill = 'right') %>%
  mutate(time = format(time, format="%H:%M:%S")) %>%
           filter(time >= "10:00:00" & time < "10:59:00") %>%
  group_by(unit, date) %>%
  summarise(mean_temp = mean(actual_temptop),
         time = paste("10:00:00"))

str(temp_plankto)
temp_plankto$date = as.Date(temp_plankto$date)
temp_plankto$time = as.hms(temp_plankto$time)
names(temp_plankto) <- c('plankton', 'date','mean_temp',  'time')
  
# check strucutre
str(temp_plankto)


## ------------------------------------------------------------------------------ ##
master <- read_delim("~/Desktop/MA/Data/Jacuzzitron 2019 Charly/master.csv", ";", escape_double = FALSE, col_types = cols(date = col_date(format = "%d.%m.%y")), 
                               trim_ws = TRUE) %>%
  select(-no)
#change value type
master$date = as.Date(master$date)
names(master) <-c('sampling', 'plankton', 'date', 'time') #change names
str(master) #check structure

#merge df
all_temp <- full_join(temp_plankto, master, by = c('plankton', 'date', 'time')) %>%
  drop_na(sampling) #can be removed again
write.csv2(x = all_temp, file = 'temperatures.csv')
## ------------------------------------------------------------------------------ ##
# Improt chlorophyll data
chlorophylla_data <- read_csv("~/Desktop/MA/Chlorophyll_a/chlorophylla_sampling11.csv") %>%
  select(-X1)
names(chlorophylla_data) = c('plankton', 'treatment', 'frequency', 'sampling', 'wavelength', 'chl', 'treatment_id')
chl_master <- full_join(all_temp, chlorophylla_data, by = c('plankton', 'sampling'))

#new df 
chl <- chl_master %>% 
  group_by(plankton, treatment, treatment_id, frequency, sampling, date, mean_temp,wavelength) %>%
  summarise( mean_chl = mean(chl, na.rm = T)) %>%
  drop_na(mean_chl)

scl = 10
ggplot(chl, aes(x = as.Date(date)) )+
  geom_line(aes(y = mean_temp))+
  geom_point(aes(y = mean_chl*scl, size = 3, col = mean_temp))+
  scale_color_gradient(low = 'blue', high = 'yellow')+
  scale_y_continuous(sec.axis = sec_axis(~./scl, name = "chlorophyll"))+ 
  facet_wrap(~plankton)+
  theme_bw()

ggsave(plot=last_plot(), file = 'mean_chl_temp.png', width = 9, height = 8)






