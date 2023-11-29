# Exploring CaCO3 precip in Sau
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


source("C:/Projects/myGit/my_r_scripts/Projects/Thetis/get_dxdy.R")


# function to extract data
extractData <- function(data, my_depth.vect, variable, depth, showPlot, showPoints){
  if (missing(showPlot)){showPlot = F}
  if (missing(showPoints)){showPoints = F}
  data.subs <- data[!is.na(data[[variable]]),]
  
  date.vect <- unique(data.subs$date)
  date.vect <- date.vect[order(date.vect)]
  
  depth.vect <- unique(data.subs[[depth]])
  depth.vect <- depth.vect[order(depth.vect)]
  
  matOut <- matrix(data = NA, ncol = length(date.vect), nrow = length(my_depth.vect))
  for (d in seq(1, length(date.vect))){
    data.thatdate <- data.subs[data.subs$date == date.vect[d],]
    if (dim(data.thatdate)[1]>3 & length(unique(data.thatdate[[depth]]))>3){
      matOut[, d] <- approx(x = data.thatdate[[depth]], y = data.thatdate[[variable]], xout = my_depth.vect)$y
    }
  }
  
  matOut <- rbind(decimal_date(date.vect), matOut)
  
  if (showPlot){
    image2D(t(matOut[-1,]),
            x = matOut[1,],
            y = my_depth.vect, ylim = c(max(my_depth.vect),0),
            xlab = "", ylab = "depth (m)",
            col = jet.col(n = 100, alpha = 1),
            main = variable)
    if (showPoints){
      for (d in seq(1, length(date.vect))){
        data.thatdate <- data.subs[data.subs$date == date.vect[d],]
        matPins.temp <- cbind(decimal_date(date.vect[d]), data.thatdate[[depth]])
        
        if (d == 1){
          matPins = matPins.temp
        } else {
          matPins = rbind(matPins, matPins.temp)
        }
      }
      points(matPins[,1], matPins[,2], col='grey', pch = 20, cex = 0.1)
    }
  }
  return(matOut)
}



# ------------------- Loading data -------------------------

setwd("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/Sau/Historical_data/")



data_alk <- read.table(file = "ALK.csv", header = T, sep = ",")
data_ca <- read.table(file = "Ca.csv", header = T, sep = ",")
data_chla <- read.table(file = "Chl-a.csv", header = T, sep = ",")
data_cond <- read.table(file = "COND.csv", header = T, sep = ",")
data_ph <- read.table(file = "pH.csv", header = T, sep = ",")
# data_phyto <- read.table(file = "Phyto_reshape.csv", header = T, sep = ",", quote="\"")


data_alk$date = as.Date(data_alk$Fecha)
data_ca$date = as.Date(data_ca$Fecha)
data_chla$date = as.Date(data_chla$Fecha)
data_cond$date = as.Date(data_cond$Fecha)
data_ph$date = as.Date(data_ph$Fecha)
# data_phyto$date = as.Date(data_phyto$Fecha)




# ggplot(data_alk[!is.na(data_alk$Valor),])+
#   geom_tile(aes(date, Depth, fill = Valor))+
#   scale_fill_viridis_c(option = "A", end = 1, direction = -1)+
#   theme_article()

depth.max <- 60
depth.step <- 0.5
my_depth.vect <- seq(0,depth.max, by = depth.step)


alk_mat <- extractData(data = data_alk, my_depth.vect = my_depth.vect, 
                        variable = "Valor", depth = "Depth", 
                        showPlot = T)
ca_mat <- extractData(data = data_ca, my_depth.vect = my_depth.vect, 
                       variable = "Valor", depth = "Depth", 
                       showPlot = T)
chla_mat <- extractData(data = data_chla, my_depth.vect = my_depth.vect, 
                      variable = "Valor", depth = "Depth", 
                      showPlot = T)
cond_mat <- extractData(data = data_cond, my_depth.vect = my_depth.vect, 
                        variable = "Valor", depth = "Depth", 
                        showPlot = T)
ph_mat <- extractData(data = data_ph, my_depth.vect = my_depth.vect, 
                        variable = "Valor", depth = "Depth", 
                        showPlot = T)



image2D(t(alk_mat[-1,]),
        x = alk_mat[1,],
        y = my_depth.vect, ylim = c(max(my_depth.vect),0),
        xlab = "", ylab = "depth (m)",
        col = cmocean("balance")(56),
        main = "alk", clim = c(0,6))






get_top_Xm <- function(mat, Xm){
  ind_select_top_Xm <- which(my_depth.vect<Xm)+1
  
  mat_avg <- apply(as.matrix(mat[ind_select_top_Xm,]), 2, mean, na.rm = T)
  df_Xm <- data.frame(date = as.Date(date_decimal(mat[1,])),
                          doy = as.numeric(strftime(as.Date(date_decimal(mat[1,])), format = "%j")),
                          values = mat_avg)
  
  return(df_Xm)
}

df_alk_10m <- get_top_Xm(mat = alk_mat, Xm = 10)
ggplot(df_alk_10m, aes(date, values))+geom_path(size=0.6)+geom_point(size=2)+theme_article()+ggtitle(paste("Alk"))

df_alk_10m_noNA <- df_alk_10m[!is.na(df_alk_10m$values),]
df_alk_10m_noNA$smoothed <- smooth.spline(x = df_alk_10m_noNA$date, df_alk_10m_noNA$values, spar = 0.4, all.knots = T)$y
df_alk_10m_noNA$dydt <- get_dxdy(decimal_date(df_alk_10m_noNA$date), df_alk_10m_noNA$smoothed)

plot(df_alk_10m_noNA$date, df_alk_10m_noNA$dydt, type = 'l')


df_Ca_10m <- get_top_Xm(mat = ca_mat, Xm = 10)
df_chla_10m <- get_top_Xm(mat = chla_mat, Xm = 10)





