

# Read HOBOs_Thermistor after bluetooth download
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

# ------------------- function -------------------------

# function to extract data
extractData <- function(data_path){
  setwd(data_path)
  lof <- list.files(path = data_path, pattern = ".txt")
  isF_file <- T
  for (f in lof){
    df_f <- read.table(file = f, header = F, skip = 3, sep = "\t", dec = ",")
    df_f <- df_f[,c(1,2)]
    names(df_f) <- c("date", "Temp")
    df_f$date <- as.POSIXct(df_f$date, format="%Y-%m-%d %H:%M:%S", tz = "CET")
    df_f$name = gsub(" ", "_",sub('\\ UnknownReadoutTime.txt$', '', f) )
    
    
    if(isF_file){
      isF_file=F
      df_out <- df_f
    } else {
      df_out <- rbind(df_out, df_f)
    }
  }
  
  return(df_out)
}



# ------------------- extract and plot -------------------------


data_path = "C:/Users/camil/Documents/MariaZambrano/data/Sau/HOBOs_T_chain/20220616_retrieved/"


df_out <- extractData(data_path)

ggplot(df_out, aes(date, Temp, colour = name))+geom_path()+theme_article()


ggplot(df_out, aes(date, Temp, colour = name))+geom_path()+theme_article()



df_out_sprd <- spread(df_out, name, Temp)



# identify and trim in-air period
df_out_sprd_10min <- as.data.frame(timeAverage(df_out_sprd, "10 min"))
df_out_sprd_10min$T_avg <- apply(as.matrix(df_out_sprd_10min[,-1]), 1, mean, na.rm = T)
df_out_sprd_10min$T_sd <- apply(as.matrix(df_out_sprd_10min[,-1]), 1, sd, na.rm = T)
df_out_sprd_10min$dT_sd_dt <- get_dxdy(as.numeric(df_out_sprd_10min$date), df_out_sprd_10min$T_sd)

df_out_sprd_10min$inAir <- F
df_out_sprd_10min$inAir[df_out_sprd_10min$T_avg > 30] <- T
df_out_sprd_10min$inAir[df_out_sprd_10min$T_sd<3] <- T
df_out_sprd_10min$inAir[abs(df_out_sprd_10min$dT_sd_dt)>1e-3] <- T

ggplot(df_out_sprd_10min, aes(date,T_avg, colour = inAir))+geom_point()+theme_article() 

df_out_sprd_10min_inWater <- df_out_sprd_10min[!df_out_sprd_10min$inAir,]
mystart <- min(df_out_sprd_10min_inWater$date)
mystop <- max(df_out_sprd_10min_inWater$date)

df_out_sprd_gath <- gather(df_out_sprd_10min_inWater, sensor, Temp, -date, -T_avg, -T_sd, -dT_sd_dt, -inAir)

ggplot(df_out_sprd_gath)+
  geom_point(aes(T_avg, Temp, colour = sensor ))+
  geom_abline(slope = 1, intercept = 0)+
  theme_article()+facet_wrap(sensor~.)



z_depth <- data.frame(name = c("Z1-20397895",  "Z10-20397907", "Z12-20397909", "Z13-20397910", "Z14-20397911", "Z15-20397912", 
                               "Z2-20397897", "Z3-20397898",  "Z4-20397899",  "Z5-20397900",  "Z6-20397901",  "Z7-20397902",  "Z8-20397903"),
                      depth = c(0.5,10,15,20,25,30,
                                1,2,3,4,5,6,7))

z_depth_sort <- sort(z_depth$depth)
z_depth_high_res <- seq(0,30,0.1)


df_out_sprd_1hour <- as.data.frame(timeAverage(df_out_sprd_10min_inWater, "1 hour"))

T_mat <- as.matrix(df_out_sprd_1hour[,seq(which(names(df_out_sprd_1hour)=="Z1-20397895"),which(names(df_out_sprd_1hour)=="Z8-20397903"))])
T_mat <- T_mat[,order(z_depth$depth)]
T_mat_apprx <- matrix(data = NA, nrow=dim(df_out_sprd_1hour)[1], ncol = length(z_depth_high_res))

for(ti in seq(1, dim(df_out_sprd_1hour)[1])){
  if(sum(!is.na(T_mat[ti,])) > 5 ){
    T_mat_apprx[ti,] <- approx(x = z_depth_sort, y = T_mat[ti,], 
                               xout = z_depth_high_res, method = "linear", rule = 2)$y
  }
  
}

date_stamp = paste(as.Date(min(df_out_sprd$date)), as.Date(max(df_out_sprd$date)), sep = " to ")

# do plot and save it

setwd("~/MariaZambrano/data/Sau/HOBOs_T_chain/")
png(paste("Sau_thermistor_chain_", 
           paste(as.Date(min(df_out_sprd$date)), as.Date(max(df_out_sprd$date)), sep = "_to_"),
           ".png",sep=""),
     width = 6, height = 4, units = 'in', res = 300)

image2D((T_mat_apprx), 
        y= z_depth_high_res, 
        x = (df_out_sprd_1hour$date), 
        ylim=c(30,0),
        col=cmocean("balance")(56), contour = T,
        xlab="",
        ylab="Depth [m]",
        clab="T [°C]",
        clim = c(8,28),
        ylab("depth [m]"),
        main = date_stamp, xaxt = "n")

# Add dates to x-axis
axis(1, as.POSIXct(unique(as.Date(df_out_sprd_1hour$date))), format( unique(as.Date(df_out_sprd_1hour$date)), "%d-%m"), )
dev.off()




# ------------------- COMPARE WITH YSI probe -------------------------


df_YSI <- read.table("C:/Users/camil/Documents/MariaZambrano/data/Sau/YSI_probe/test_YSI_20220531.csv", header = T, sep = ";")
df_YSI$date <- as.POSIXct(paste(df_YSI$ï..date, df_YSI$time, " "), format="%d.%m.%Y %H:%M:%S", tz = "UTC")
df_YSI <- df_YSI[,-c(1,2)]

df_YSI_gath <- gather(df_YSI[-seq(1,20),], variable, value, -date)
ggplot(df_YSI_gath, aes(date, value))+geom_path()+
  theme_article()+
  scale_x_datetime(date_labels = "%u")+xlab("day in the week")+
  facet_wrap(variable~., scales = "free_y")


# ggplot(df_YSI, aes(Cond_uS_cm, SpCond_uS_cm))+geom_point()+theme_article()


ggplot(df_out_sprd_gath, aes(date, Temp, colour = sensor))+
  geom_path(data = df_YSI, aes(date, Temp_degC, colour = "YSI"), color = "black")+
  geom_path()+
  theme_article()+facet_wrap(sensor~.)


df_out_sprd$YSI <- approx(df_YSI$date, df_YSI$Temp_degC, xout = df_out_sprd$date, method = "linear")$y

df_out_sprd_gath <- gather(df_out_sprd, sensor, Temp, -date, -T_avg, -YSI)

ggplot(df_out_sprd_gath)+
  geom_point(aes(YSI, Temp, colour = sensor ))+
  geom_abline(slope = 1, intercept = 0)+
  theme_article()+facet_wrap(sensor~.)

ggplot(df_out_sprd_gath)+
  geom_path(aes(date, Temp-T_avg, colour = sensor ))+
  geom_abline(slope = 0, intercept = 0)+
  theme_article()+facet_wrap(sensor~.)+
  theme(legend.position = 'none')+
  scale_x_datetime(date_labels = "%u")+xlab("day in the week")











