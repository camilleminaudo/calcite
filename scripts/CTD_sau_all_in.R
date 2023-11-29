# Assemble CTD profiles into one dataframe
# Camille Minaudo
# June 2022
# ------------------- Load packages -------------------------


cat("\014")
rm(list = ls())


library(ggplot2)
library(gridExtra)
library(grid)
library(plot3D)
library(lubridate)
library(dplyr)
library(tidyr)
library(zoo)
library(egg)
library(openair)
library(rLakeAnalyzer)
library(cmocean)
library(stringi)




data_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/processed/"

lof = list.files(path = data_path, pattern = ".csv", full.names = T)

isF = T
for (f in lof){
  data.temp <- read.table(f, header = T, sep = ";")
  if (isF){
    isF = F
    data_CTD <- data.temp
  } else {
    data_CTD <- rbind(data_CTD, data.temp)
  }
}

data_CTD$profile_date <- as.factor(as.Date(data_CTD$date))

p_Temp <- ggplot(data_CTD, aes(Temp_degC, depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Water temperature [ºC]")+ylab("Depth [m]")
p_DO <- ggplot(data_CTD, aes(DO_mg_L, depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Dissolved Oxygen [mg/L]")+ylab("Depth [m]")
p_DOsat <- ggplot(data_CTD, aes(DOsat_., depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Dissolved Oxygen [%sat]")+ylab("Depth [m]")
p_Chla <- ggplot(data_CTD, aes(Chla_ug_L, depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Chlorophyll a [ug/L]")+ylab("Depth [m]")
p_turb <- ggplot(data_CTD, aes(turb_FTU, depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Turbidity [FTU]")+ylab("Depth [m]")
p_pH <- ggplot(data_CTD, aes(pH, depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("pH")+ylab("Depth [m]")



data.temp$profile_date <- as.factor(as.Date(data.temp$date))
ggplot(data.temp, aes(DO_mg_L, depth_m))+
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 40), fill = "grey75", alpha= 0.2)+
  geom_path()+scale_y_reverse()+theme_article()+
  ggtitle(paste("Panta de Sau", unique(data.temp$profile_date)))+
  xlab("Dissolved Oxygen [mg/L]")+ylab("Depth [m]")

ggplot(data.temp, aes(DOsat_., depth_m, colour = profile_date))+geom_path()+scale_y_reverse()+theme_article()+xlab("Dissolved Oxygen [%sat]")+ylab("Depth [m]")


p_Temp+facet_wrap(profile_date~.)



