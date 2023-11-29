# Extract CTD data from SBE19plus v2 ICRA (Rafa)
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

source("C:/Projects/my_r_scripts/Projects/Thetis/get_dxdy.R")

# -------------------  function to extract data -------------------------

extractCTD <- function(myfile){
  data_raw <- readLines(myfile)
  rows_first_6 <- substr(data_raw, start = 1, stop = 6)
  cast_ID_vect <- as.vector(unlist(strsplit(data_raw[which(rows_first_6=="* cast")], " ")))
  Sys.setlocale("LC_TIME", "English")
  my_date <- as.POSIXct(paste(cast_ID_vect[seq(5,8)], sep = " ", collapse=" "), format = "%d %b %Y %H:%M:%S", tz = "UTC")
  
  data = read.table(myfile, header = F, skip = 563)
  # name 0 = chloroflTC0: Chlorophyll, Turner Cyclops [ug/l]
  # name 1 = c0uS/cm: Conductivity [uS/cm]
  # name 2 = depFM: Depth [fresh water, m]
  # name 3 = sbeox0Mg/L: Oxygen, SBE 43 [mg/l]
  # name 4 = sbeox0PS: Oxygen, SBE 43 [% saturation]
  # name 5 = par: PAR/Irradiance, Biospherical/Licor
  # name 6 = ph: pH
  # name 7 = prdE: Pressure, Strain Gauge [psi]
  # name 8 = specc: Specific Conductance [uS/cm]
  # name 9 = tv290C: Temperature [ITS-90, deg C]
  # name 10 = timeS: Time, Elapsed [seconds]
  # name 11 = seaTurbMtr: Turbidity, Seapoint [FTU]
  # name 12 = flag:  0.000e+00
  
  data <- data[,-13]
  names(data) <- c("chla", "cond","depth","DO","DOsat", "PAR", "pH", "press", "spec_Cond", "temp", "time", "turb")
  
  # select only the downward cast
  i <- seq(1,length(data$depth))
  z_smooth <- smooth.spline(i,data$depth, spar = 0.3)$y
  dzdt <- get_dxdy(i,z_smooth)
  ind_downcast <- which(dzdt > 0.02)
  data_keep <- data[ind_downcast,]
  
  ind_z_max = which.max(data_keep$depth)
  data_keep <- data_keep[seq(1, ind_z_max),]
  
  # additionnal selection based on conductivity
  i <- seq(1,length(data_keep$depth))
  # cond_smooth <- smooth.spline(i,data_keep$cond, spar = 0.05)$y
  dconddt <- get_dxdy(i,data_keep$cond/mean(data_keep$cond))
  ind_outlier <- which(abs(dconddt)>0.005)
  if(length(ind_outlier)>0){
    data_keep <- data_keep[-ind_outlier,]
  }
  
  
  plot(data$time, data$depth, type = 'l')
  lines(data_keep$time, data_keep$depth, col='red', lwd=3)
  
  # 
  # # additionnal selection
  # i <- seq(1,length(data_keep$depth))
  # z_smooth <- smooth.spline(i,data_keep$depth, spar = 0.3)$y
  # dzdt <- get_dxdy(i,z_smooth)
  # # d2zdt2 <- get_dxdy(i,dzdt)
  # ind_downcast <- which(dzdt > 0.02)
  # data_keep <- data_keep[seq(ind_downcast[1], max(ind_downcast)),]
  # data_keep <- data_keep[data_keep$depth > 0,]
  
  data_keep$date <- my_date+data_keep$time
  return(data_keep)
}


# psa_file <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/config_files/DatCnv_Sau.psa"
# config_file <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/config_files/SBE19plus-6858-April2021_rafa.xmlcon"
# 
# sbe_cmd <- function(cmd, in_file, out_dir, xmlcon, psa){
#   
#   exec_str <- paste(cmd, in_file, out_dir, xmlcon, psa, sep = " ")
#   shell(exec_str)
# }

# ------------------- Loading data -------------------------

datapath = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/raw_temp/"
dataplot = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/plots/"
processpath = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/processed/"

files_hex = list.files(datapath, pattern = ".hex", full.names = T)
# myfile_hex = files_hex[2]

# path_to_sbe <- "C:\\Program Files (x86)\\Sea-Bird\\SBEDataProcessing-Win32"
cmd = paste("DatCnvW.exe")
system(cmd)
# sbe_cmd(cmd = cmd, in_file = myfile_hex, out_dir = datapath, xmlcon = config_file, psa = psa_file)


files = list.files(datapath, pattern = ".cnv", full.names = T)
myfile = files[1]

data = extractCTD(myfile)

names(data) <-  c("Chla_ug_L", "cond_uS_cm","depth_m","DO_mg_L","DOsat_%", "PAR", "pH", "press_PSI", "spec_Cond_uS_cm", 
                  "Temp_degC", "time", "turb_FTU", "date")


# reorganize table
mydata <- data[,c("depth_m",
                "date",
                "press_PSI", "Temp_degC",
                "cond_uS_cm", "spec_Cond_uS_cm",
                "DO_mg_L", "DOsat_%",
                "PAR", "pH",
                "Chla_ug_L", "turb_FTU"),]

# write table
setwd(processpath)
write.table(mydata, file = paste("CTD_Sau_",as.Date(min(mydata$date)),".csv", sep = ""), sep = ";", row.names = F, col.names = T)


data.gath <- gather(mydata[,-2], var, value, -depth_m)
data.gath$var <- factor(data.gath$var, levels = unique(data.gath$var))
# make plot
plt <- ggplot(data.gath, aes(value, depth_m))+geom_path()+facet_wrap(.~var, scales = "free_x", nrow = 2)+
  ggtitle(paste("Sau Reservoir, ",min(data$date), " UTC",sep=""))+
  scale_y_reverse()+theme_article()
ggsave(plt, filename = paste("CTD_Sau_",as.Date(min(data$date)),".png",sep=""), 
       path = dataplot, width = 12, height = 8, dpi = 300, units = 'in', bg = "white")
plt

# copy raw files into raw_backup
filestocopy <- list.files(datapath, full.names = T, recursive = T)
targetdir <- "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/CTD/raw_backup/"
file.copy(from=filestocopy, to=targetdir, 
          overwrite = T, recursive = T, 
          copy.mode = TRUE)

# remove raw files from raw_temp
file.remove(filestocopy)


