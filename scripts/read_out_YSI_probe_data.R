
# Read YSI_sonde
# Camille Minaudo
# April 2022
# ------------------- Packages -------------------------


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

source("~/GitRepos/my_r_scripts/Projects/Thetis/get_dxdy.R")


# ------------------- extract and plot -------------------------


df_YSI <- read.table("C:/Users/camil/Documents/MariaZambrano/data/Sau/YSI_probe/SAU_2206.txt", header = F, sep = ",", skip = 5)
names(df_YSI) <- c("Date","Time","Temp","SpCond","Cond","pH","pH_mV","DOsat","DO","Battery")
df_YSI$date <- as.POSIXct(paste(df_YSI$Date, " ",df_YSI$Time, sep=""), format="%d/%m/%y %H:%M:%S", tz = "UTC")
df_YSI <- df_YSI[,-c(1,2)]

df_YSI_gath <- gather(df_YSI[-seq(1,20),], variable, value, -date)
ggplot(df_YSI_gath, aes(date, value))+geom_path()+
  theme_article()+
  scale_x_datetime(date_labels = "%d")+xlab("day in the month")+
  facet_wrap(variable~., scales = "free_y")





