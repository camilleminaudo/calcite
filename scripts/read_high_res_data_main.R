# Read all high res data from Sau, merge the different variables into a single database.
# 
# Camille Minaudo
# July 2022
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
library(stringr)
library(EMD)


source("C:/Projects/myGit/my_r_scripts/Projects/Thetis/get_dxdy.R")
source("C:/Projects/myGit/my_r_scripts/Projects/Thetis/my_integrate.R")


plotpath = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/plots"


# ------------------- Functions -------------------------

extractData_YSI <- function(YSI_file){
  
  df_YSI <- read.table(YSI_file, header = F, sep = ",", skip = 5)
  names(df_YSI) <- c("Date","Time","Temp","SpCond","Cond","pH","pH_mV","DOsat","DO","Battery")
  df_YSI$date <- as.POSIXct(paste(df_YSI$Date, " ",df_YSI$Time, sep=""), format="%d/%m/%y %H:%M:%S", tz = "UTC")
  df_YSI <- df_YSI[,-c(1,2)]
  df_YSI <- df_YSI[df_YSI$Cond>300,]
  
  df_YSI <- df_YSI[order(as.numeric(df_YSI$date)),]
  
  return(df_YSI)
}

# function to extract data from HOBO thermistors
extractData_thermistors <- function(data_path){
  
  sensor_map <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_T_chain/sensor_map.csv", sep = ",", 
                           header = T)
  sensor_map$date_start <- as.POSIXct(sensor_map$date_start, format = "%d/%m/%Y %H:%M", tz = 'CET')
  sensor_map$date_stop <- as.POSIXct(sensor_map$date_stop, format = "%d/%m/%Y %H:%M", tz = 'CET')
  
  sensor_map$date_stop[is.na(sensor_map$date_stop)] <- as.POSIXct(Sys.Date())
  
  # z_depth <- data.frame(name = c("Z1-20397895","Z2-20397897", "Z3-20397898",  "Z4-20397899",  
  #                                "Z5-20397900",  "Z6-20397901",  "Z7-20397902",  "Z8-20397903", 
  #                                "Z9-20397905", "Z10-20397907", "Z12-20397909", "Z13-20397910", 
  #                                "Z14-20397911", "Z15-20397912"),
  #                       depth = c(0.5,1,2,3,
  #                                 4,5,6,7,
  #                                 8,9,10,15,
  #                                 20,30))
  

  setwd(data_path)
  lof <- list.files(path = data_path, pattern = ".txt", recursive = T)
  lof <- lof[!str_detect(lof, "test")]
  lof <- lof[!str_detect(lof, "~BROMIUM")]
  isF_file <- T
  for (f in lof){
    print(f)
    df_f <- read.table(file = f, header = F, skip = 3, sep = "\t", dec = ",")
    df_f <- df_f[,c(1,2)]
    names(df_f) <- c("date", "Temp")
    df_f$date <- as.POSIXct(df_f$date, format="%Y-%m-%d %H:%M:%S", tz = "CET")
    df_f$name = gsub(" ", "_",sub('\\ UnknownReadoutTime.txt$', '', basename(f)) )
    
    
    df_f$depth <- NA
    ind_corresp <- which(sensor_map$sensor == unique(df_f$name))
    # sensor_map[ind_corresp,]
    for (i in seq(1,length(ind_corresp))){
      df_f$depth[df_f$date>=sensor_map$date_start[ind_corresp[i]] & 
                   df_f$date < sensor_map$date_stop[ind_corresp[i]]] <- sensor_map$depth[ind_corresp[i]]
    }
    
    if(isF_file){
      isF_file=F
      df_out <- df_f
    } else {
      df_out <- rbind(df_out, df_f)
    }
  }
  # ind_depth_corresp_to_data <- match(df_out$name, sensor_map$sensor)
  # df_out$depth <- z_depth$depth[ind_depth_corresp_to_data]
  
  return(df_out)
}





extractData_DO <- function(myfile){
  
  mydepth_m <- unlist(strsplit(basename(myfile), "_"))[3]
  mydepth <- as.numeric(sub('m', '', mydepth_m))+0.5
  
  df_f <- read.table(file = myfile, header = F, skip = 2, sep = ",")
  df_f <- df_f[,c(2,3,4)]
  names(df_f) <- c("date", "O2", "Temp")
  df_f$date <- as.POSIXct(df_f$date, format="%d/%m/%y %H:%M:%S", tz = "UTC")
  df_f <- df_f[df_f$date > (min(df_f$date)+2*60*60) & df_f$date < (max(df_f$date)-2*60*60),]
  
  # identify abnormal in-air temperature on day of deployment
  df_f_depl <- df_f[as.Date(df_f$date) == min(as.Date(df_f$date)),]
  df_f_depl$dTdt <- get_dxdy(as.numeric(df_f_depl$date)/3600, df_f_depl$Temp)
  if (max(df_f_depl$Temp)>30){ # we consider 30?C as an impossible value in-water
    ind_min_dTdt <- which.min(df_f_depl$dTdt)
    date_start_in_water <- df_f_depl$date[ind_min_dTdt]
    df_f <- df_f[df_f$date > date_start_in_water,]
  }
  
  # identify abnormal in-air temperature on day of retrieval
  df_f_retrv <- df_f[as.Date(df_f$date) == max(as.Date(df_f$date)),]
  df_f_retrv$dTdt <- get_dxdy(as.numeric(df_f_retrv$date)/3600, df_f_retrv$Temp)
  if (max(df_f_retrv$Temp)>30){ # we consider 30?C as an impossible value in-water
    ind_min_dTdt <- which.min(df_f_retrv$dTdt)
    date_stop_in_water <- df_f_retrv$date[ind_min_dTdt]
    df_f <- df_f[df_f$date < date_start_in_water,]
  }
  
  df_f$depth <- mydepth
  
  
  return(df_f)
}


extractData_HOBO_light <- function(data_path){
  
  z_depth <- data.frame(name = c("INCOMING","UNDERW1", "UNDERW2"),
                        depth = c(0,6.5,9.5))
  
  
  setwd(data_path)
  lof <- list.files(path = data_path, pattern = ".csv", recursive = T)
  isF_file <- T
  for (f in lof){
    print(f)
    df_f <- read.table(file = f, header = F, skip = 2, sep = ",")
    df_f <- df_f[,c(2,3,4)]
    names(df_f) <- c("date", "Temp", "Lux")
    df_f$date <- as.POSIXct(df_f$date, format="%d/%m/%y %H:%M:%S", tz = "CET")
    df_f$name = gsub(".csv", "",sub('Sau_light_', '', basename(f)))
    
    df_f$sensor = unlist(strsplit(gsub(".csv", "",sub('Sau_light_', '', basename(f))), "_"))[1]
    
    if(isF_file){
      isF_file=F
      df_out <- df_f
    } else {
      df_out <- rbind(df_out, df_f)
    }
  }
  
  ind_depth_corresp_to_data <- match(df_out$sensor, z_depth$name)
  
  df_out$depth <- z_depth$depth[ind_depth_corresp_to_data]
  
  
  uniq_ID_ind <- !duplicated(paste(df_out$date, df_out$depth, "-"))
  df_out <- df_out[uniq_ID_ind, ]
  
  
  return(df_out)
}



get_EMD <- function(df, varName, my_max.imf){
  
  y <- df[[varName]][!is.na(df[[varName]])]
  t  <- as.numeric(df$date[!is.na(df[[varName]])])/3600/24 - as.numeric(df$date[!is.na(df[[varName]])])[1]/3600/24
  
  myEMD <- emd(y,t, boundary = 'periodic', max.imf = my_max.imf)
  # myEMD$nimf
  trend_emd = myEMD$residue
  detrended_emd = y-myEMD$residue
  
  par(mfrow=c(3,1))
  plot(t, y, type = 'l', main = paste("EMD of ", varName,sep=""))
  plot(t, trend_emd, type = 'l')
  plot(t, detrended_emd, type = 'l')
  abline(h=0)
  par(mfrow=c(1,1))
  
  df_out <- data.frame(time = t,
                       date = df$date[!is.na(df[[varName]])],
                       varName = varName,
                       y = y,
                       trend = trend_emd,
                       detrended = detrended_emd)
  return(df_out)
}


#---------------------  load YSI probe data --------------------- 


data_path_YSI <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/YSI_probe/"
setwd(data_path_YSI)
lof <- list.files(data_path_YSI, pattern = ".txt")
isF <- T
for(f in lof){
  print(f)
  df_YSI.temp <- extractData_YSI(YSI_file = f)
  # df_YSI.temp <- as.data.frame(timeAverage(df_YSI.temp, "1 hour"))
  df_YSI.temp$filename = f
  if (isF){
    isF = F
    df_YSI <- df_YSI.temp
  } else {
    df_YSI <- rbind(df_YSI, df_YSI.temp)
  }
}
df_YSI <- df_YSI[order(as.numeric(df_YSI$date)),]

uniq_ID_ind <- !duplicated(df_YSI$date)
df_YSI <- df_YSI[uniq_ID_ind, ]

df_YSI_gath <- gather(df_YSI[,-which(names(df_YSI) =="filename")], variable, value, -date)
ggplot(df_YSI_gath, aes(date, value))+geom_path()+
  theme_article()+
  # scale_x_datetime(date_labels = "%d")+xlab("day in the month")+
  facet_wrap(variable~., scales = "free_y")


df_YSI_hour <- as.data.frame(timeAverage(df_YSI, "1 hour"))
ggplot(df_YSI_hour, aes(date, DO))+geom_path()+theme_article()+ggtitle("DO at 7.5 m  - YSI probe")



#---------------------  load HOBO DO data --------------------- 

data_path_O2 <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_DO/"
setwd(data_path_O2)
lof <- list.files(data_path_O2, pattern = ".csv")
isF <- T
for(f in lof){
  print(f)
  df_DO.temp <- extractData_DO(f)
  df_DO.temp <- as.data.frame(timeAverage(df_DO.temp, "1 hour"))
  df_DO.temp$filename = f
  if (isF){
    isF = F
    df_DO <- df_DO.temp
  } else {
    df_DO <- rbind(df_DO, df_DO.temp)
  }
}

uniq_ID_ind <- !duplicated(paste(df_DO$date, df_DO$depth, "-"))
df_DO <- df_DO[uniq_ID_ind, ]
df_DO <- df_DO[df_DO$O2>=0,]

df_DO$sensor_depth <- paste("HOBO",df_DO$depth,sep="_")


df_YSI_formated <- data.frame(date = df_YSI$date,
                              O2 = df_YSI$DO,
                              Temp = df_YSI$Temp,
                              
                              depth = 7.5,
                              filename = "YSI_probe",
                              sensor_depth = "YSI_7.5")

df_DO <- rbind(df_DO, df_YSI_formated)

df_DO <- df_DO[df_DO$date > as.POSIXct("2022-06-01"),]

plt_DO <- ggplot(df_DO, aes(date, O2, colour = as.factor(sensor_depth)))+
  geom_path()+theme_article()+ggtitle(paste("DO in Sau", sep = ""))+
  theme(legend.position = 'none', legend.title = element_blank())+
  ylab("Dissolved oxygen [mg/L]")+xlab("")+scale_x_datetime(date_breaks = "1 month")+
  facet_wrap(as.factor(sensor_depth)~., nrow = 3)
plt_DO
ggsave(plot = plt_DO,  filename = "DO_timeseries_Sau.png", path = plotpath,
       width = 7, height = 6, units = "in", dpi = 300, bg = 'white')






df_DO_sprd <- tidyr::spread(df_DO[,c("depth","date", "O2")], key = depth, value = O2)

df_DO_sprd_30min <- as.data.frame(timeAverage(df_DO_sprd, "1 hour"))
df_DO_sprd_30min <- df_DO_sprd_30min[df_DO_sprd_30min$date>as.POSIXct("2022-06-02"),]




DO_mat <- as.matrix(df_DO_sprd_30min[,-1])
z_DO <- as.numeric(names(df_DO_sprd_30min[,-1]))


z_high_res_DO <- seq(0,30,0.5)

DO_mat_appx <- matrix(data = NA, nrow=dim(DO_mat)[1], ncol = length(z_high_res_DO))

for(ti in seq(1, dim(DO_mat)[1])){
  if(sum(!is.na(DO_mat[ti,])) > 1 ){
    DO_mat_appx[ti,] <- approx(x = z_DO, y = as.numeric(DO_mat[ti,]), 
                               xout = z_high_res_DO, method = "linear", rule = 1:1)$y
  }
  
}

date_stamp = paste(as.Date(min(df_DO_sprd_30min$date)), as.Date(max(df_DO_sprd_30min$date)), sep = " to ")


image2D((DO_mat_appx), 
        y= z_high_res_DO, 
        x = (df_DO_sprd_30min$date), 
        ylim=c(10,0),
        col=cmocean("delta")(56), contour = F,
        xlab="",
        ylab="Depth [m]",
        clab="DO [mg O2/L]",
        clim = c(floor(min(DO_mat_appx, na.rm=T))-1,round(max(DO_mat_appx, na.rm=T))+1),
        ylab("depth [m]"),
        main = date_stamp, xaxt = "n")
# Add dates to x-axis
axis(1, as.POSIXct(unique(as.Date(df_DO_sprd_30min$date))), format( unique(as.Date(df_DO_sprd_30min$date)), "%d-%m"), )



# require(LakeMetabolizer)
# 
# metab_test <- LakeMetabolizer::metab.bookkeep(do.obs = df_YSI$DO, 
#                                               do.sat = o2.at.sat.base(df_YSI$Temp, 
#                                                                       baro = rep(0, length(df_YSI$temp)),
#                                                                       altitude = 500, 
#                                                                       salinity = rep(0, length(df_YSI$Temp)), 
#                                                                       model = "garcia-benson"),
#                                               k.gas, z.mix, irr, ...)





#---------------------  load thermistor chain data --------------------- 


data_path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_T_chain/"
df_out <- extractData_thermistors(data_path)


# adding TºC info from DO HOBOs and YSI
df_DO_Tformatted <- data.frame(date = df_DO$date,
                               Temp =df_DO$Temp,
                               name = paste("DO_",df_DO$depth,"m",sep = ""),
                               depth = df_DO$depth)

df_out <- rbind(df_out, df_DO_Tformatted)

uniq_ID_ind <- !duplicated(paste(df_out$date, df_out$name, "-"))
df_out <- df_out[uniq_ID_ind, ]

df_out$ID_z <- paste(df_out$name,df_out$depth,sep = "_")


# ggplot(df_out, aes(date, Temp, colour = name))+geom_path()+theme_article()
df_out_sprd <- spread(df_out[,c("date", "ID_z", "Temp")], ID_z, Temp)

# identify and trim in-air period
df_out_sprd_10min <- as.data.frame(timeAverage(df_out_sprd, "30 min"))
df_out_sprd_10min$T_avg <- apply(as.matrix(df_out_sprd_10min[,-1]), 1, mean, na.rm = T)
df_out_sprd_10min$T_sd <- apply(as.matrix(df_out_sprd_10min[,-1]), 1, sd, na.rm = T)
df_out_sprd_10min$dT_sd_dt <- get_dxdy(as.numeric(df_out_sprd_10min$date), df_out_sprd_10min$T_sd)


df_out_sprd_10min_sd <- as.data.frame(timeAverage(df_out_sprd, "30 min", statistic = "sd"))
df_out_sprd_10min_sd$T_avg <- apply(as.matrix(df_out_sprd_10min_sd[,-1]), 1, mean, na.rm = T)
df_out_sprd_10min_sd$inAir <- F
df_out_sprd_10min_sd$inAir[df_out_sprd_10min_sd$T_avg > 0.15] <- T

ggplot(df_out_sprd_10min_sd, aes(date,T_avg))+geom_point()+
  geom_abline(slope = 0, intercept = 0.15)+
  theme_article()


df_out_sprd_10min$inAir <- F
df_out_sprd_10min$inAir[df_out_sprd_10min$T_avg > 28] <- T
# df_out_sprd_10min$inAir[df_out_sprd_10min$T_sd<0.3] <- T
df_out_sprd_10min$inAir[abs(df_out_sprd_10min$dT_sd_dt)>1.5e-4] <- T
df_out_sprd_10min$inAir[df_out_sprd_10min_sd$T_avg > 0.15] <- T

ggplot(df_out_sprd_10min, aes(date,abs(T_avg), colour = inAir))+geom_point()+
  geom_abline(slope = 0, intercept = 1.5e-4)+
  theme_article()



df_out_sprd_10min_inWater <- df_out_sprd_10min[!df_out_sprd_10min$inAir,]
mystart <- min(df_out_sprd_10min_inWater$date)
mystop <- max(df_out_sprd_10min_inWater$date)

df_out_sprd_gath <- gather(df_out_sprd_10min_inWater, sensor, Temp, -date, -T_avg, -T_sd, -dT_sd_dt, -inAir)

# ggplot(df_out_sprd_gath)+
#   geom_point(aes(T_avg, Temp, colour = sensor ))+
#   geom_abline(slope = 1, intercept = 0)+
#   theme_article()+facet_wrap(sensor~.)

df_out_sprd_10min_inAir <- df_out_sprd_10min[df_out_sprd_10min$inAir,]
df_out_sprd_gath_inAir <- gather(df_out_sprd_10min_inAir, sensor, Temp, -date, -T_avg, -T_sd, -dT_sd_dt, -inAir)

ggplot(df_out_sprd_gath_inAir)+
  geom_point(aes(T_avg, Temp, colour = sensor ))+
  geom_abline(slope = 1, intercept = 0)+
  theme_article()+facet_wrap(sensor~.)



df_out_sprd_10min_inWater <- df_out_sprd_10min[!df_out_sprd_10min$inAir,]
df_out_sprd_gath_inWater <- gather(df_out_sprd_10min_inWater, sensor, Temp, -date, -T_avg, -T_sd, -dT_sd_dt, -inAir)

ggplot(df_out_sprd_gath_inWater)+
  geom_point(aes(T_avg, Temp, colour = sensor ))+
  geom_abline(slope = 1, intercept = 0)+
  theme_article()+facet_wrap(sensor~.)

# grid to 1 hour matrix
z_depth_high_res <- seq(0,30,0.5)
df_out_sprd_1hour <- as.data.frame(timeAverage(df_out_sprd_10min_inWater, "1 hour"))

df_sel <- df_out_sprd_1hour[,match(unique(df_out$ID_z), names(df_out_sprd_1hour))]
df_depths <-  df_out$depth[match(names(df_sel), df_out$ID_z)]
df_ordered <- df_sel[,order(df_depths)]
df_depths_ordered <- df_depths[order(df_depths)]

T_mat <- as.matrix(df_ordered)
T_mat_apprx <- matrix(data = NA, nrow=dim(df_ordered)[1], ncol = length(z_depth_high_res))

for(ti in seq(1, dim(df_out_sprd_1hour)[1])){
  if(sum(!is.na(T_mat[ti,])) > 2 ){
    my_T_vect <- as.numeric(T_mat[ti,])
    my_depths <- df_depths_ordered[!is.na(my_T_vect)]
    my_T_vect <- my_T_vect[!is.na(my_T_vect)]
    T_mat_apprx[ti,] <- approx(x = my_depths, y = my_T_vect, 
                               xout = z_depth_high_res, method = "linear", rule = 2:1)$y   # rule = 2:1
    
    # plot(z_depth_high_res, T_mat_apprx[ti,] )
  }
  
}

date_stamp = paste(as.Date(min(df_out_sprd_1hour$date)), as.Date(max(df_out_sprd_1hour$date)), sep = " to ")


# do plot and save it
setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_T_chain/plots/")
png(paste("Sau_thermistor_chain_", 
          paste(as.Date(min(df_out_sprd_1hour$date)), as.Date(max(df_out_sprd_1hour$date)), sep = "_to_"),
          ".png",sep=""),
    width = 6, height = 4, units = 'in', res = 300)

image2D((T_mat_apprx), 
        y= z_depth_high_res, 
        x = (df_out_sprd_1hour$date), 
        ylim=c(30,0),
        col=cmocean("thermal")(56), contour = list(col = "black", alpha = 0.5, NAcol = "white", drawlabels = F, nlevels = 5, lwd = 0.5),
        xlab="",
        ylab="Depth [m]",
        clab="T [ºC]",
        clim = c(floor(min(T_mat_apprx, na.rm=T))-1,round(max(T_mat_apprx, na.rm=T))+1),
        ylab("depth [m]"),
        main = date_stamp, xaxt = "n")
# Add dates to x-axis
myXticks <- as.POSIXct(unique(as.Date(df_out_sprd_1hour$date)))
myXticks <- myXticks[seq(1,length(myXticks),5)]
axis(1, myXticks, format(myXticks, "%d-%m"), )
abline(v = as.POSIXct(c("2022-06-16, 12:00", "2022-07-12, 12:00", "2022-08-03, 12:00", "2022-08-29, 12:00", "2022-10-13, 12:00")))
dev.off()



#---------------------  load CTD data --------------------- 
data_CTD_path <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/processed/"

lof <- list.files(data_CTD_path, pattern = ".csv")
setwd(data_CTD_path)
isF = T
for(f in lof){
  print(f)
  df_CTD <- read.table(f, header = T, sep = ";")
  df_CTD <- df_CTD[order(df_CTD$depth_m),]
  if(isF){
    isF = F
    df_CTDs <- df_CTD
  }
  df_CTDs <- rbind(df_CTDs, df_CTD)
}
df_CTDs$date_profile <- as.factor(as.Date(df_CTDs$date))

p_CTD_T <- ggplot(df_CTDs, aes(depth_m, Temp_degC, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Water temperature [ºC]")+
  xlab("Depth [m]")

p_CTD_DO <- ggplot(df_CTDs, aes(depth_m, DO_mg_L, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Dissolved O2 [mg/L]")+
  xlab("Depth [m]")

p_CTD_DOsat <- ggplot(df_CTDs, aes(depth_m, DOsat_., colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Dissolved O2 [%sat]")+
  xlab("Depth [m]")

p_CTD_CHLa <- ggplot(df_CTDs, aes(depth_m, Chla_ug_L, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Chlorophyll a [ug/L]")+
  xlab("Depth [m]")

p_CTD_pH <- ggplot(df_CTDs, aes(depth_m, pH, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("pH [-]")+
  xlab("Depth [m]")

p_CTD_turb <- ggplot(df_CTDs, aes(depth_m, turb_FTU, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Turbidity [FTU]")+
  xlab("Depth [m]")


p_CTD_cond <- ggplot(df_CTDs, aes(depth_m, cond_uS_cm, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Conductivity [uS/cm]")+
  xlab("Depth [m]")


p_CTD_SpeCond <- ggplot(df_CTDs, aes(depth_m, spec_Cond_uS_cm, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+facet_wrap(date_profile~.)+theme(legend.position = "none")+
  ylab("Specific conductivity [uS/cm]")+
  xlab("Depth [m]")


p_CTD_PAR <- ggplot(df_CTDs, aes(depth_m, PAR, colour = date_profile))+geom_line()+theme_article()+coord_flip()+
  scale_x_reverse()+theme(legend.position = "none")+scale_y_log10()+
  ylab("PAR [*]")+
  xlab("Depth [m]")+facet_wrap(date_profile~.)


ggsave(plot = p_CTD_T, filename = "CTD_T_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_DO, filename = "CTD_DO_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_DOsat, filename = "CTD_DOsat_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_CHLa, filename = "CTD_CHLa_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_turb, filename = "CTD_turbidity_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_cond, filename = "CTD_conductivity_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_SpeCond, filename = "CTD_SpecificCond_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_pH, filename = "CTD_pH_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = p_CTD_PAR, filename = "CTD_PAR_Sau.png", path = plotpath, width = 5, height = 4, units = "in", dpi = 300, bg = 'white')


#---------------------  Blend CTD TºC data with high-res TºC data --------------------- 

T_mat_incl_CTD <- T_mat_apprx
for(ki in seq(1,length(df_CTDs$date))){
  zi <- which.min(abs(df_CTDs$depth_m[ki] - z_depth_high_res))
  ti <- which.min(abs(as.POSIXct(df_CTDs$date [ki], tz = 'UTC') - df_out_sprd_1hour$date))
  
  doIT <- T
  if(abs(z_depth_high_res[zi] - df_CTDs$depth_m[ki])>0.5){doIT = F}
  if(abs(as.numeric(difftime(df_out_sprd_1hour$date[ti], df_CTDs$date [ki], units = "hours")))>2){doIT = F}
  if(!is.na(T_mat_incl_CTD[ti, zi])){doIT = F}
  if (doIT){T_mat_incl_CTD[ti, zi] = df_CTDs$Temp_degC[ki]}
}

# fill the gaps temporally
T_mat_incl_CTD_apprx <- T_mat_incl_CTD
for(zi in seq(1, length(z_depth_high_res))){
  if(sum(!is.na(T_mat_incl_CTD[,zi])) > 5 ){
    ind_no_nans <- !is.na(T_mat_incl_CTD[,zi])
    T_mat_incl_CTD_apprx[,zi] <- approx(x = df_out_sprd_1hour$date, y = as.numeric(T_mat_incl_CTD[,zi]), 
                               xout = df_out_sprd_1hour$date, method = "linear", rule = 2:1)$y   # rule = 2:1
    ind_nans <- is.na(T_mat_incl_CTD_apprx[,zi])
    T_mat_incl_CTD_apprx[,zi][!ind_nans] <- smooth.spline(x = df_out_sprd_1hour$date[!ind_nans], y = T_mat_incl_CTD_apprx[,zi][!ind_nans], all.knots = F, spar = 0.9)$y
    T_mat_incl_CTD_apprx[,zi][ind_no_nans] <- T_mat_incl_CTD[,zi][ind_no_nans] 
    # plot(z_depth_high_res, T_mat_apprx[ti,] )
  }
}
# Smoothing vertically
T_mat_incl_CTD_apprx_smooth <- T_mat_incl_CTD_apprx
for(ti in seq(1, dim(T_mat_incl_CTD_apprx)[1])){
  ind_nans <- is.na(T_mat_incl_CTD_apprx_smooth[ti,])
  T_mat_incl_CTD_apprx_smooth[ti,][!ind_nans] <- smooth.spline(x = z_depth_high_res[!ind_nans], y = T_mat_incl_CTD_apprx[ti,][!ind_nans], all.knots = F, spar = 0.5)$y
}

# do plot and save it
setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_T_chain/plots/")
png(paste("Sau_thermistor_chain_", 
          paste(as.Date(min(df_out_sprd_1hour$date)), as.Date(max(df_out_sprd_1hour$date)), sep = "_to_"),
          "_interpolated.png",sep=""),
    width = 6, height = 4, units = 'in', res = 300)

image2D((T_mat_incl_CTD_apprx_smooth), 
        y= z_depth_high_res, 
        x = df_out_sprd_1hour$date, 
        ylim=c(30,0),
        col=cmocean("thermal")(56), contour = list(col = "black", alpha = 0.5, NAcol = "white", drawlabels = F, nlevels = 5, lwd = 0.5),
        xlab="",
        ylab="Depth [m]",
        clab="T [ºC]",
        clim = c(floor(min(T_mat_incl_CTD_apprx_smooth, na.rm=T))-1,round(max(T_mat_incl_CTD_apprx_smooth, na.rm=T))+1),
        ylab("depth [m]"),
        main = date_stamp, xaxt = "n")
# Add dates to x-axis
myXticks <- as.POSIXct(unique(as.Date(df_out_sprd_1hour$date)))
myXticks <- myXticks[seq(1,length(myXticks),5)]
axis(1, myXticks, format(myXticks, "%d-%m"), )
abline(v = as.POSIXct(unique(as.Date(df_CTDs$date))))
dev.off()


# Heat content
get_rho_w <- function(Temp){
  rho_water = 4e-5*Temp^3 - 0.0078*Temp^2 + 0.0646*Temp + 999.86  # kg.m-3
  return(rho_water)
}

get_Cp <- function(Temp){
  Cp = 1e+3*(4.2174 - 3.6608e-3*Temp + 1.3129e-4*Temp^2 - 2.210e-6*Temp^3 + 1.508e-8*Temp^4) # J kg-1 deg-1
  return(Cp)
}

Heat_content <- T_avg <- NA*seq(1, dim(T_mat_incl_CTD_apprx_smooth)[1])
for(ti in seq(1, dim(T_mat_incl_CTD_apprx_smooth)[1])){
  Tin <- T_mat_incl_CTD_apprx_smooth[ti,]
  Heat_profile <- (Tin)*get_Cp(Temp = Tin)*get_rho_w(Temp = Tin) # J/m3
  
  Heat_content[ti] <- my_integrate(x = z_depth_high_res, y = Heat_profile, from = 0, to = 30)/1e+6 # MJ/m2
  T_avg[ti] <- mean(Tin, na.rm = T)
}


p_T_selected_z <- ggplot()+
  geom_line(data = data.frame(date = df_out_sprd_1hour$date, Temp = T_mat_incl_CTD_apprx_smooth[,1]), aes(date, Temp, colour= "0 m"))+
  geom_line(data = data.frame(date = df_out_sprd_1hour$date, Temp = T_mat_incl_CTD_apprx_smooth[,21]), aes(date, Temp, colour= "10 m"))+
  geom_line(data = data.frame(date = df_out_sprd_1hour$date, Temp = T_mat_incl_CTD_apprx_smooth[,41]), aes(date, Temp, colour= "20 m"))+
  geom_line(data = data.frame(date = df_out_sprd_1hour$date, Temp = T_mat_incl_CTD_apprx_smooth[,61]), aes(date, Temp, colour= "30 m"))+
  theme_article()+theme(legend.title = element_blank())+ #, legend.position = c(0.9,0.9)
  ylab("Water temperature [ºC]")+
  xlab("")


# compute thermocline and mixed layer depth
require(rLakeAnalyzer)
isF_d <- T
for(di in unique(as.Date(df_out_sprd_1hour$date))){
  mat_T_that_day <- T_mat_incl_CTD_apprx_smooth[which(as.Date(df_out_sprd_1hour$date) == di),]
  T_avg_that_day <- apply(mat_T_that_day, MARGIN = 2, mean, na.rm=T)
  # ind_no_NA <- !is.na(T_avg_that_day)
  
  df_thermo.depth_temp <- data.frame(date = di,
                                z_th = thermo.depth(wtr = T_avg_that_day[seq(1,30)], depths = z_depth_high_res[seq(1,30)], Smin = 0.05, seasonal = F),
                                z_ML = z_depth_high_res[which.min(abs(T_avg_that_day - (T_avg_that_day[1] - 0.4)))])
  if(isF_d){
    isF_d = F
    df_thermo.depth <- df_thermo.depth_temp
  } else {
    df_thermo.depth <- rbind(df_thermo.depth, df_thermo.depth_temp)
  }
}
df_thermo.depth$date <- as.Date(df_thermo.depth$date)

p_z_th <- ggplot(df_thermo.depth)+
  geom_line(aes(date, z_th, colour = "thermocline"))+
  geom_line(aes(date, z_ML, colour = "mixed layer"))+
  ylim(c(40,0))+theme_article()+theme(legend.title = element_blank())+
  ylab("Depth [m]")+
  xlab("")


ggsave(plot = ggarrange(p_T_selected_z, p_z_th, nrow = 2), filename = "T_selected_z_Sau.png", path = plotpath, width = 9, height = 6, units = "in", dpi = 300, bg = 'white')




#---------------------  Blend CTD DO data with high-res DO data --------------------- 

DO_mat_incl_CTD <- DO_mat_appx
for(ki in seq(1,length(df_CTDs$date))){
  zi <- which.min(abs(df_CTDs$depth_m[ki] - z_high_res_DO))
  ti <- which.min(abs(as.POSIXct(df_CTDs$date [ki], tz = 'UTC') - df_DO_sprd_30min$date))
  
  doIT <- T
  if(abs(z_high_res_DO[zi] - df_CTDs$depth_m[ki])>0.5){doIT = F}
  if(abs(as.numeric(difftime(df_DO_sprd_30min$date[ti], df_CTDs$date [ki], units = "hours")))>2){doIT = F}
  if(!is.na(DO_mat_incl_CTD[ti, zi])){doIT = F}
  if (doIT){DO_mat_incl_CTD[ti, zi] = df_CTDs$DO_mg_L[ki]}
}

# fill the gaps temporally
DO_mat_incl_CTD_apprx <- DO_mat_incl_CTD
for(zi in seq(1, length(z_high_res_DO))){
  if(sum(!is.na(DO_mat_incl_CTD[,zi])) > 5 ){
    ind_no_nans <- !is.na(DO_mat_incl_CTD[,zi])
    DO_mat_incl_CTD_apprx[,zi] <- approx(x = df_DO_sprd_30min$date, y = as.numeric(DO_mat_incl_CTD[,zi]), 
                                        xout = df_DO_sprd_30min$date, method = "linear", rule = 2:1)$y   # rule = 2:1
    ind_nans <- is.na(DO_mat_incl_CTD_apprx[,zi])
    DO_mat_incl_CTD_apprx[,zi][!ind_nans] <- smooth.spline(x = df_DO_sprd_30min$date[!ind_nans], y = DO_mat_incl_CTD_apprx[,zi][!ind_nans], all.knots = F, spar = 0.9)$y
    DO_mat_incl_CTD_apprx[,zi][ind_no_nans] <- DO_mat_incl_CTD[,zi][ind_no_nans] 
    # plot(z_high_res_DO, T_mat_apprx[ti,] )
  }
}
# Smoothing vertically
DO_mat_incl_CTD_apprx_smooth <- DO_mat_incl_CTD_apprx
for(ti in seq(1, dim(DO_mat_incl_CTD_apprx)[1])){
  ind_nans <- is.na(DO_mat_incl_CTD_apprx_smooth[ti,])
  if(sum(ind_nans) < length(ind_nans)){
    DO_mat_incl_CTD_apprx_smooth[ti,][!ind_nans] <- smooth.spline(x = z_high_res_DO[!ind_nans], 
                                                                  y = DO_mat_incl_CTD_apprx[ti,][!ind_nans], 
                                                                  all.knots = F, spar = 0.5)$y
  }
}


DO_mat_incl_CTD_apprx_smooth[DO_mat_incl_CTD_apprx_smooth<0] <- 0

date_stamp = paste(as.Date(min(df_DO_sprd_30min$date)), as.Date(max(df_DO_sprd_30min$date)), sep = " to ")

# do plot and save it
setwd(plotpath)
png(paste("Sau_DO_", date_stamp,
          "_interpolated.png",sep=""),
    width = 6, height = 4, units = 'in', res = 300)

image2D((DO_mat_incl_CTD_apprx_smooth), 
        y= z_high_res_DO, 
        x = (df_DO_sprd_30min$date), 
        ylim=c(30,0),
        col=cmocean("delta")(56), contour = F,
        xlab="",
        ylab="Depth [m]",
        clab="DO [mg O2/L]",
        clim = c(0,round(max(DO_mat_incl_CTD_apprx_smooth, na.rm=T))+1),
        ylab("depth [m]"),
        main = date_stamp, xaxt = "n")
# Add dates to x-axis
myXticks <- as.POSIXct(unique(as.Date(df_DO_sprd_30min$date)))
myXticks <- myXticks[seq(1,length(myXticks),5)]
axis(1, myXticks, format(myXticks, "%d-%m"), )
abline(v = as.POSIXct(unique(as.Date(df_CTDs$date))))
dev.off()


#---------------------  Load light data --------------------- 


data_path = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/HOBOs_light/"
df_Lux <- extractData_HOBO_light(data_path)

df_Lux$depth_char <- paste(df_Lux$depth,"m",sep = "")

df_Lux_sprd <- spread(df_Lux[,c("depth_char","date", "Lux")], depth_char, Lux)
df_Lux_sprd_h <- as.data.frame(timeAverage(df_Lux_sprd, "1 hour", statistic = "max"))

ggplot(df_Lux_sprd_h)+
  geom_path(aes(date, `0m`, colour = "0+"))+
  geom_path(aes(date, `6.5m`, colour = "6.5"))+
  geom_path(aes(date, `9.5m`, colour = "9.5"))+
  theme_article()+scale_y_log10()


df_DO_sprd <- spread(df_DO[,c("depth","date", "O2")], depth, O2)

df_DO_sprd_30min <- as.data.frame(timeAverage(df_DO_sprd, "1 hour"))
df_DO_sprd_30min <- df_DO_sprd_30min[df_DO_sprd_30min$date>as.POSIXct("2022-06-02"),]


#---------------------  EMD YSI probe analysis --------------------- 

EMD_Temp <- get_EMD(df = df_YSI_hour, varName = "Temp", my_max.imf = 6)
EMD_pH <- get_EMD(df = df_YSI_hour, varName = "pH", my_max.imf = 6)
EMD_DO <- get_EMD(df = df_YSI_hour, varName = "DO", my_max.imf = 5)
EMD_DOsat <- get_EMD(df = df_YSI_hour, varName = "DOsat", my_max.imf = 5)
EMD_SpCond <- get_EMD(df = df_YSI_hour, varName = "SpCond", my_max.imf = 5)
EMD_Cond <- get_EMD(df = df_YSI_hour, varName = "Cond", my_max.imf = 6)


p1.T <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, Temp))+
  geom_path(data = EMD_Temp, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("Water temperature")+xlab("")+ylab("[ºC]")
p2.T <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_Temp, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[ºC]")

pT <- ggarrange(p1.T, p2.T)



p1.DO <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, DO))+
  geom_path(data = EMD_DO, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("Dissolved O2")+xlab("")+ylab("[mg/L]")
p2.DO <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_DO, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[mg/L]")

pDO <- ggarrange(p1.DO, p2.DO)


p1.DOsat <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, DOsat))+
  geom_path(data = EMD_DOsat, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("Dissolved O2")+xlab("")+ylab("[%sat]")
p2.DOsat <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_DOsat, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[%sat]")

pDOsat <- ggarrange(p1.DOsat, p2.DOsat)



p1.pH <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, pH))+
  geom_path(data = EMD_pH, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("pH")+xlab("")+ylab("[-]")
p2.pH <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_pH, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[-]")

ppH <- ggarrange(p1.pH, p2.pH)



p1.Cond <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, Cond))+
  geom_path(data = EMD_Cond, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("Conductivity")+xlab("")+ylab("[uS/cm]")
p2.Cond <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_Cond, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[uS/cm]")

pCond <- ggarrange(p1.Cond, p2.Cond)



p1.SpCond <- ggplot()+
  geom_path(data = df_YSI_hour, aes(date, SpCond))+
  geom_path(data = EMD_SpCond, aes(date, trend), size=2, alpha = 0.3)+
  theme_article()+ggtitle("Specific conductivity")+xlab("")+ylab("[uS/cm]")
p2.SpCond <- ggplot()+
  geom_abline(slope = 0, intercept = 0)+
  geom_path(data = EMD_SpCond, aes(date, detrended))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[uS/cm]")

pSpCond <- ggarrange(p1.SpCond, p2.SpCond)


ggsave(plot = pT, filename = "T_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = pDO, filename = "DO_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = pDOsat, filename = "DOsat_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = ppH, filename = "pH_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = pCond, filename = "Cond_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')
ggsave(plot = pSpCond, filename = "SpCond_EMD_YSI.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')


# plot(EMD_DOsat$detrended, EMD_SpCond$detrended)
# plot(EMD_Temp$detrended, EMD_SpCond$detrended)


df_YSI_sd_daily <- as.data.frame(timeAverage(df_YSI, "1 day", statistic = "sd"))
# ggplot(df_YSI_sd_daily, aes(date, SpCond))+geom_path()+theme_article()


df_YSI_sd_daily_gath <- gather(df_YSI_sd_daily[,-c(which(names(df_YSI_sd_daily)=="Battery"), 
                                                   which(names(df_YSI_sd_daily)=="pH_mV"))], 
                               variable, value, -date)
ggplot(df_YSI_sd_daily_gath, aes(date, value))+geom_path()+
  theme_article()+
  facet_wrap(variable~., scales = "free_y")+ylab("Daily amplitude")+xlab("")


ggplot(df_YSI_sd_daily, aes(pH, DOsat))+geom_point()+theme_article()



# ggplot(EMD_DO, aes(hour(date), detrended))+geom_line()+theme_article()+facet_wrap(as.Date(date)~.)+geom_abline(slope = 0, intercept = 0)






require(biwavelet)
SpCond_7.5 = biwavelet::wt(data.frame(seq(1,length(EMD_SpCond$date)), EMD_SpCond$detrended),max.scale = 240)
plot(SpCond_7.5, zlim = c(0,12), main="SpCond at 7.5 m", ncol = 32, plot.cb = T)

SpCond_7.5_near_diel_power <- SpCond_7.5$power[SpCond_7.5$period>20 & SpCond_7.5$period<27,]
mean_SpCond_7.5_near_diel_power <- apply(SpCond_7.5_near_diel_power, 2, mean)
EMD_SpCond$isDiel <- F
EMD_SpCond$isDiel[mean_SpCond_7.5_near_diel_power>100] <- T


p_SpCond_emd <- ggplot()+
  geom_vline(xintercept = EMD_SpCond$date[EMD_SpCond$isDiel], size=1.2, colour = "lightblue", alpha=0.25)+
  # geom_point(data = EMD_SpCond[EMD_SpCond$isDiel,], aes(date, y, colour = "diel"), size=2, alpha=0.5)+
  geom_path(data = df_YSI_hour, aes(date, SpCond))+
  geom_path(data = EMD_SpCond, aes(date, trend), size=2, alpha = 0.2)+
  theme_article()+ylab("Specific cond. [uS/cm]")+xlab("")+
  scale_x_datetime(breaks = "1 month", minor_breaks = "1 week", limits = as.POSIXct(c("2022-06-01", "2022-12-31")))+
  scale_colour_viridis_d(option = "C", begin = 0.2, end = 0.8, direction = -1)+theme(legend.position = c(0.1,0.9))
p_SpCond_emd


# ggsave(plot = p_SpCond_emd, filename = "SpCond_EMD_YSI_dielID.png", path = plotpath, width = 6, height = 4, units = "in", dpi = 300, bg = 'white')



DOsat_7.5_wl = biwavelet::wt(data.frame(seq(1,length(EMD_DOsat$date)), EMD_DOsat$detrended),max.scale = 240)
plot(DOsat_7.5_wl, zlim = c(0,12), main="DOsat at 7.5 m", ncol = 32, plot.cb = T)

DOsat_7.5_near_diel_power <- DOsat_7.5_wl$power[DOsat_7.5_wl$period>20 & DOsat_7.5_wl$period<27,]
mean_DOsat_7.5_near_diel_power <- apply(DOsat_7.5_near_diel_power, 2, mean)
EMD_DOsat$isDiel <- F
EMD_DOsat$isDiel[mean_DOsat_7.5_near_diel_power>700] <- T


p_DOsat_emd <- ggplot()+
  geom_vline(xintercept = EMD_DOsat$date[EMD_DOsat$isDiel], size=1.2, colour = "lightblue", alpha=0.25)+
  # geom_point(data = EMD_DOsat[EMD_DOsat$isDiel,], aes(date, y, colour = "diel"), size=2, alpha=0.5)+
  geom_path(data = df_YSI_hour, aes(date, DOsat))+
  geom_path(data = EMD_DOsat, aes(date, trend), size=2, alpha = 0.2)+
  theme_article()+ylab("DOsat [%sat]")+xlab("")+
  scale_x_datetime(breaks = "1 month", minor_breaks = "1 week", limits = as.POSIXct(c("2022-06-01", "2022-12-31")))+
  scale_colour_viridis_d(option = "C", begin = 0.2, end = 0.8, direction = -1)+theme(legend.position = c(0.1,0.9))
p_DOsat_emd


p_SpCond_DO_dielID <- ggarrange(p_DOsat_emd, p_SpCond_emd, nrow = 2)

ggsave(plot = p_SpCond_DO_dielID, filename = "SpCond_DO_EMD_YSI_dielID.png", 
       path = plotpath, scale = 0.5, width = 15, height = 8, units = "in", dpi = 300, bg = 'white')








Temp_7.5_wl = biwavelet::wt(data.frame(seq(1,length(EMD_Temp$date)), EMD_Temp$detrended),max.scale = 240)
plot(Temp_7.5_wl, zlim = c(0,12), main="Temp at 7.5 m", ncol = 32, plot.cb = T)

Temp_7.5_near_diel_power <- Temp_7.5_wl$power[Temp_7.5_wl$period>20 & Temp_7.5_wl$period<27,]
mean_Temp_7.5_near_diel_power <- apply(Temp_7.5_near_diel_power, 2, mean)
plot(mean_Temp_7.5_near_diel_power)
EMD_Temp$isDiel <- F
EMD_Temp$isDiel[mean_Temp_7.5_near_diel_power>700] <- T


p_Temp_emd <- ggplot()+
  geom_vline(xintercept = EMD_Temp$date[EMD_Temp$isDiel], size=1.2, colour = "lightblue", alpha=0.25)+
  # geom_point(data = EMD_Temp[EMD_Temp$isDiel,], aes(date, y, colour = "diel"), size=2, alpha=0.5)+
  geom_path(data = df_YSI_hour, aes(date, Temp))+
  geom_path(data = EMD_Temp, aes(date, trend), size=2, alpha = 0.2)+
  theme_article()+ylab("Temp [%sat]")+xlab("")+
  scale_x_datetime(breaks = "1 month", minor_breaks = "1 week")+
  scale_colour_viridis_d(option = "C", begin = 0.2, end = 0.8, direction = -1)+theme(legend.position = c(0.1,0.9))
p_Temp_emd




#---------------------  EMD DO data analysis --------------------- 

# 


EMD_DO.2.5 <- get_EMD(df = df_DO_sprd_30min, varName = "2.5", my_max.imf = 4)
EMD_DO.5.5 <- get_EMD(df = df_DO_sprd_30min, varName = "5.5", my_max.imf = 5)
EMD_DO.7.5 <- get_EMD(df = df_DO_sprd_30min, varName = "7.5", my_max.imf = 5)


EMD_DOs <- rbind(EMD_DO.2.5, EMD_DO.5.5, EMD_DO.7.5)

plt_EMD_DOs_detrended <- ggplot(EMD_DOs)+
  geom_abline(slope = 0, intercept = 0)+
  geom_line(aes(date, detrended, colour = varName))+
  theme_article()+ggtitle("Detrended")+xlab("")+ylab("[mg/L]")+
  theme(legend.position = 'none', legend.title = element_blank())+
  facet_wrap(varName~., nrow = 3)


plt_EMD_DOs_trend <- ggplot(EMD_DOs)+
  geom_abline(slope = 0, intercept = 0)+
  geom_line(aes(date, trend, colour = varName), size=1.2)+
  geom_line(aes(date, y, colour = varName), alpha=0.4)+
  theme_article()+ggtitle("Trend")+xlab("")+ylab("[mg/L]")+
  theme(legend.position = 'none', legend.title = element_blank())+
  facet_wrap(varName~., nrow = 3)

plt_EMD_DOs <- ggarrange(plt_EMD_DOs_trend, plt_EMD_DOs_detrended, nrow = 1)

ggsave(plot = plt_EMD_DOs, filename = "EMD_DOs_Sau.png", path = plotpath,
       width = 10, height = 6, units = "in", dpi = 300, bg = 'white')


ggplot(EMD_DOs)+
  # geom_abline(slope = 0, intercept = 0)+
  geom_line(aes(date, trend, colour = varName))+
  theme_article()+ggtitle("Trend")+xlab("")+ylab("[mg/L]")+
  theme(legend.position = 'none', legend.title = element_blank())+
  facet_wrap(varName~., nrow = 3)


# ggplot(EMD_DO.7.5, aes(hour(date), detrended))+geom_path()+facet_wrap(as.Date(date)~.)+theme_article()


require(biwavelet)
wt_2.5 = biwavelet::wt(data.frame(seq(1,length(EMD_DO.2.5$date)), EMD_DO.2.5$detrended),max.scale = 240)
wt_5.5 = biwavelet::wt(data.frame(seq(1,length(EMD_DO.5.5$date)), EMD_DO.5.5$detrended),max.scale = 240)
wt_7.5 = biwavelet::wt(data.frame(seq(1,length(EMD_DO.7.5$date)), EMD_DO.7.5$detrended),max.scale = 240)

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/plots/")
png("wavelet_DO_2.5.png", width = 8, height = 6, units = 'in', bg = 'white', res = 300)
plot(wt_2.5, zlim = c(0,10), main="DO at 2.5 m")
dev.off()

png("wavelet_DO_5.5.png", width = 8, height = 6, units = 'in', bg = 'white', res = 300)
plot(wt_5.5, zlim = c(0,10), main="DO at 5.5 m")
dev.off()

png("wavelet_DO_7.5.png", width = 8, height = 6, units = 'in', bg = 'white', res = 300)
plot(wt_7.5, zlim = c(0,10), main="DO at 7.5 m")
dev.off()


# ggplot(df_DO[df_DO$date > as.POSIXct("2022-06-02"),], aes(Temp, O2, colour = date))+geom_point()+facet_wrap(.~sensor_depth)+theme_article()

# ggplot(df_YSI_formated, aes(Temp, O2))+geom_point()+geom_smooth(method = "loess", colour = 'red')+
#   facet_wrap(.~as.Date(date))+theme_article()


# df_DO_sd_daily <- as.data.frame(timeAverage(df_DO_sprd_30min, "1 day", statistic = "sd"))
# 
# df_DO_sd_daily_gath <- gather(df_DO_sd_daily, variable, value, -date)
# ggplot(df_DO_sd_daily_gath, aes(date, value))+geom_path()+
#   theme_article()+
#   facet_wrap(variable~., scales = "free_y")+ylab("Daily amplitude")+xlab("")


#---------------------   --------------------- 

