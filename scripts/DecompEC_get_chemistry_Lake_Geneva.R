
# Camille, June 2022

# -------------------  LONG TERM LEMAN PROJECT ------------------- #
#
# This program
# ------------------------------------------------------------------- #

cat("\014")
rm(list = ls())

library(ggplot2)
library(gridExtra)
library(grid)
library(plot3D)
library(hexView)
library(lubridate)
library(dplyr)
library(tidyr)
library(zoo)
library(egg)

# source('C:/Users/camil/Documents/Programs/Thetis/find_ML_depth.R')

# -------------- FUNCTIONS ------------------

# function to extract data
extractData <- function(data, my_depth.vect, variable, showPlot, showPoints){
  if (missing(showPlot)){showPlot = F}
  if (missing(showPoints)){showPoints = F}
  data.subs <- data[!is.na(data[[variable]]),]
  
  date.vect <- unique(data.subs$sampling.date)
  date.vect <- date.vect[order(date.vect)]
  
  depth.vect <- unique(data.subs$depth)
  depth.vect <- depth.vect[order(depth.vect)]
  
  matOut <- matrix(data = NA, ncol = length(date.vect), nrow = length(my_depth.vect))
  for (d in seq(1, length(date.vect))){
    data.thatdate <- data.subs[data.subs$sampling.date == date.vect[d],]
    if (dim(data.thatdate)[1]>3 & length(unique(data.thatdate$depth))>3){
      matOut[, d] <- approx(x = data.thatdate$depth, y = data.thatdate[[variable]], xout = my_depth.vect)$y
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
        data.thatdate <- data.subs[data.subs$sampling.date == date.vect[d],]
        matPins.temp <- cbind(decimal_date(date.vect[d]), data.thatdate$depth)
        
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


compare2vars_annual <- function(mat1, mat2, label1, label2){
  year_v1 <- year(date_decimal(mat1[1,]))
  year_v2 <- year(date_decimal(mat2[1,]))
  
  mean_v1 <- apply(mat1[-1,], 2, mean, na.rm = T)
  mean_v2 <- apply(mat2[-1,], 2, mean, na.rm = T)
  
  isFirst = T
  for (my_y in seq(max(min(year_v1), min(year_v2)), min(max(year_v1), max(year_v2)))){
    df.temp <- data.frame(year = my_y,
                          v1_mean = mean(mean_v1[year_v1 == my_y], na.rm = T),
                          v2_mean = mean(mean_v2[year_v2 == my_y], na.rm = T),
                          v1_min = min(mean_v1[year_v1 == my_y], na.rm = T),
                          v2_min = min(mean_v2[year_v2 == my_y], na.rm = T),
                          v1_max = max(mean_v1[year_v1 == my_y], na.rm = T),
                          v2_max = max(mean_v2[year_v2 == my_y], na.rm = T))
    if (isFirst){
      df <- df.temp
      isFirst = F
    } else {
      df <- rbind(df, df.temp)
    }
  }
  p1 <- ggplot(df)+
    geom_line(aes(year, v1_min), color = 'grey')+
    geom_line(aes(year, v1_max), color = 'grey')+
    geom_line(aes(year, v1_mean), color = 'red')+
    ylab(label1)+
    theme_bw()
  p2 <- ggplot(df)+
    geom_line(aes(year, v2_min), color = 'grey')+
    geom_line(aes(year, v2_max), color = 'grey')+
    geom_line(aes(year, v2_mean), color = 'red')+
    ylab(label2)+
    theme_bw()
  
  p3 <- ggplot(df)+
    geom_point(aes(v1_mean, v2_mean, colour = year))+
    xlab(label1)+
    ylab(label2)+
    theme_bw()
  
  p.all <- grid.arrange(p1, p2, p3, ncol = 2)
}


# -------------------- SETTINGS ---------------------------

depth.max <- 300
depth.step <- 0.5
my_depth.vect <- seq(0,depth.max, by = depth.step)



# -------------------- Loading long term PHYSIC-CHEMICAL DATA ---------------------------


dataChem <- read.table("~/Data/CIPEL/extractions/extraction_physico_chimie_03-09-2019-14-43-08-031/leman__1567514588378-648095.csv",
                       header = T, sep = ';')
dataChem$sampling.date <- as.Date(dataChem$sampling.date, format="%d/%m/%Y", tz = 'GMT')

names(dataChem) <- c("project.name","site.name", "platform.name", "sampling.date","sampling.tool" ,"measuring.tool" ,
                     "min.depth", "max.depth" , "depth","Temp",
                     "TN","PON","NO3" ,"NH4", "NO2", "pH", "SiO2",
                     "Conductivity", "Total.Alkalinity", "ionic.balance", "TOC","DOC",  "POC",
                     "Ca", "Mg", "Na", "K" , "Cl",
                     "SO4", "O2" , "TSM", "TP", "TDP", "PP", "SRP")
units.vect <- c("-","-", "-", "date","-" ,"-" ,
                "m", "m" , "m", "°C",
                "mgN.L-1","mgN.L-1","mgN.L-1" ,"mgN.L-1",
                "mgN.L-1", "-","mg.L-1",
                "uS.cm-1", "meq.L-1", "%", "mgC.L-1","mgC.L-1",  "mgC.L-1"  ,
                "mg.L-1" , "mg.L-1" , "mg.L-1", "mg.L-1", "mg.L-1",
                "mg.L-1", "mgO2.L-1" ,"mg.L-1","mgP.L-1", "mgP.L-1", "mgP.L-1", "mgP.L-1")
dataChem$depth2 <- (dataChem$min.depth + dataChem$max.depth)/2
dataChem$depth[is.na(dataChem$depth)] <- dataChem$depth2[is.na(dataChem$depth)] 
dataChem <- dataChem[!is.na(dataChem$depth),-which(names(dataChem) == "depth2")]


dataChem_EC <- dataChem[!is.na(dataChem$Conductivity),]
dataChem_temp <- dataChem[!is.na(dataChem$Temp),]
dataChem_EC$Temp <- approx(dataChem_temp$sampling.date, dataChem_temp$Temp, dataChem_EC$sampling.date)$y # °C

# select some majore elements
data_majorE <- dataChem_EC[dataChem_EC$depth<30,
                           c("Temp","pH","Total.Alkalinity","NO3", "SO4","Cl", "Na", "Mg", "K" , "Ca")]

IDs <- paste("A",seq(1,length(data_majorE$pH)), sep = "")
data_majorE <- cbind(IDs, data_majorE)

write.table(data_majorE, file = "~/MariaZambrano/phreeqc/major_Ions_Lake_Geneva.txt",
            sep = "\t", row.names = F, col.names = T, quote = F, na = "\t")


dataChem_EC$EC25 <- dataChem_EC$Conductivity/ (1.8831 -0.054506*dataChem_EC$Temp + 
                                                 0.00095782*(dataChem_EC$Temp)^2 - 7.6262e-6 * (dataChem_EC$Temp)^3)



# -------------------- read out PHREEQC results ---------------------------
df_EC <- read.table("~/MariaZambrano/phreeqc/output_EC.dat", header = F, sep = "\t", skip = 1)
names(df_EC) <- c("Sample_ID", "EC_tot", "EC_majel", "EC_error", "EC_H", "EC_Ca",
                   "EC_Cl", "EC_K", "EC_NO3", "EC_Na", "EC_Mg", "EC_SO4", "EC_HCO3", "charge_error", "void")
df_EC <- df_EC[,-which(names(df_EC) == "void")]
l_df <- length(df_EC$Sample_ID)
df_EC$date <- dataChem_EC$sampling.date[seq(1,l_df)]
df_EC$depth <- dataChem_EC$depth[seq(1,l_df)]
df_EC$EC25 <- dataChem_EC$EC25[seq(1,l_df)]
df_EC$doy <- as.numeric(strftime(df_EC$date, format = "%j"))

ggplot(df_EC, aes(date, EC_majel/EC_tot*100, colour = depth))+
  scale_colour_viridis_c(begin = 0.1, end = 0.9, option = "B")+
  geom_point()+
  ylab("contribution of major ions to total EC")+
  theme_article()



# ggplot(df_EC, aes(EC_majel/EC_tot*100, depth))+
#   geom_point()+
#   ylab("contribution of major ions to total EC")+
#   theme_article()+scale_y_reverse()+facet_wrap(date~.)


df_EC_gath <- gather(df_EC, species, EC, -Sample_ID, -EC_tot, -EC_error, -charge_error, -date, -doy, -depth, -EC25)

ggplot(df_EC_gath, aes(date, EC/EC_tot*100, colour = species))+
  # scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "B")+
  geom_point()+theme_article()

ggplot(df_EC_gath, aes(date, EC/EC25*100, colour = species))+
  # scale_colour_viridis_d(begin = 0.1, end = 0.9, option = "B")+
  geom_point()+theme_article()


ggplot(df_EC_gath, aes(doy, EC/EC_tot*100, colour = year(date)))+
  geom_point()+theme_article()+facet_wrap(species~.)


ggplot(df_EC_gath, aes(EC, depth, colour = species))+
  geom_path()+
  theme_article()+scale_y_reverse()+facet_wrap(date~.)






# a coefficients
df_a <- read.table("~/MariaZambrano/phreeqc/output_a.dat", header = F, sep = "\t", skip = 1)
names(df_a) <- c("Sample_ID", "a_H", "a_Ca","a_Cl", "a_K", "a_NO3", "a_Na", "a_Mg", "a_SO4", "a_HCO3", "void")
df_a <- df_a[,-which(names(df_a) == "void")]
l_df <- length(df_a$Sample_ID)
df_a$date <- dataChem_EC$sampling.date[seq(1,l_df)]
df_a$depth <- dataChem_EC$depth[seq(1,l_df)]
df_a$EC25 <- dataChem_EC$Conductivity[seq(1,l_df)]
df_a$doy <- as.numeric(strftime(df_a$date, format = "%j"))

ggplot(df_a, aes(date, a_H))+
  geom_point()+
  geom_path()+
  theme_article()


# ggplot(df_a, aes(a_Ca, depth))+
#   geom_point()+
#   theme_article()+scale_y_reverse()+facet_wrap(date~.)


tab_a_mean <- apply(as.matrix(df_a[,seq(2,10)]), 2, mean)
matrix_a_mean <- matrix(data = t(as.vector(tab_a_mean)), nrow = dim(df_a)[1], ncol = length(tab_a_mean), byrow = T)

df_a_rel2mean <- df_a
df_a_rel2mean[,seq(2,10)] <- (df_a_rel2mean[,seq(2,10)]-matrix_a_mean)/(matrix_a_mean)*100

ggplot(df_a_rel2mean, aes(date, a_Ca))+
  geom_point()+
  theme_article()+scale_y_reverse()


df_a_gath <- gather(df_a_rel2mean, species, a_coeff, -Sample_ID, -date, -doy, -depth, -EC25)

ggplot(df_a_gath, aes(doy, a_coeff, colour = species))+
  geom_point()+
  geom_path()+
  theme_article()+facet_wrap(species~.)





