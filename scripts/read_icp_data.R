

library(ggplot2)
library(stringi)
library(stringr)
library(egg)
library(lubridate)
library(openair)
library(tidyverse)


# ------------------- Functions -------------------------

extractData_YSI <- function(YSI_file){

  df_YSI <- read.table(YSI_file, header = F, sep = ",", skip = 5)
  names(df_YSI) <- c("Date","Time","Temp","SpCond","Cond","pH","pH_mV","DOsat","DO","Battery")
  df_YSI$date <- as.POSIXct(paste(df_YSI$Date, " ",df_YSI$Time, sep=""), format="%d/%m/%y %H:%M:%S", tz = "UTC")
  df_YSI <- df_YSI[,-c(1,2)]
  df_YSI <- df_YSI[df_YSI$Cond>300,]

  return(df_YSI)
}

# ------------------- Loading data -------------------------

#
icp <- read.table("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/Sed_water_sampling/analyses_ICP_sau_gergal.csv",
                  skip = 8, header = T, sep = ",")
icp$date_sampling <- as.Date(icp$date_sampling, format = "%d/%m/%Y")

icp_avg <- data.frame(site = icp$site,
                      ID = icp$ID_real,
                      ID_lab = icp$ID,
                      depth = as.factor(icp$depth),
                      date = icp$date_sampling,
                      sample_type = icp$sample_type,
                      Ca = (icp$Ca+icp$Ca.1)/2, #ppm
                      Na = (icp$Na+icp$Na.1)/2,
                      Mg = (icp$Mg+icp$Mg.1)/2,
                      K = icp$K,
                      Sr = (icp$Sr+icp$Sr.1)/2,
                      Al = (icp$Al+icp$Al.1+icp$Al.2)/3,
                      Fe = (icp$Fe+icp$Fe.1+icp$Fe.2)/3)

icp_avg$Ca_mmol <- icp_avg$Ca/40.0780      # mmol Ca
icp_avg$Na_mmol <- icp_avg$Na/22.989769280 # mmol Na
icp_avg$Mg_mmol <- icp_avg$Mg/24.305       # mmol Mg
icp_avg$K_mmol <- icp_avg$K/39.09830       # mmol K
icp_avg$Sr_mmol <- icp_avg$Sr/87.6200      # mmol Sr
icp_avg$Al_mmol <- icp_avg$Al/26.98153860  # mmol Al
icp_avg$Fe_mmol <- icp_avg$Fe/55.845       # mmol Fe


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

uniq_ID_ind <- !duplicated(df_YSI$date)
df_YSI <- df_YSI[uniq_ID_ind, ]


# ------------------- Do plots -------------------------


icp_gath <- gather(icp_avg, values, variable, -ID, -ID_lab, -depth, -date, -sample_type)

# ggplot(icp_gath, aes(date, values, colour = depth))+
#   # geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
#   geom_point(size=2, alpha = 0.5)+
#   xlab("")+
#   theme_article()+facet_wrap(variable~sample_type, nrow = 2, scales = "free_y")



p_Ca <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, Ca, colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  ylab("Ca [mg/L]")+
  theme_article()+theme(legend.position = c(0.1,0.7))+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)

p_Mg <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, Mg, colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  ylab("Mg [mg/L]")+
  theme_article()+theme(legend.position = 'none')+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)

p_Sr <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, Sr, colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  ylab("Sr [mg/L]")+
  theme_article()+theme(legend.position = 'none')+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)

p_Sr_Ca <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, (Sr_mmol)/(Ca_mmol)*1000, colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  ylab("Sr/Ca [mmol/mol]")+
  theme_article()+theme(legend.position = 'none')+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)

p_Mg_Ca <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, (Mg_mmol)/(Ca_mmol), colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  ylab("Mg/Ca [mol/mol]")+
  theme_article()+theme(legend.position = 'none')+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)

p_Mg_Sr_Ca <- ggplot(icp_avg[icp_avg$sample_type=="water",], aes((Sr_mmol)/(Ca_mmol), (Mg_mmol)/(Ca_mmol), colour = depth))+
  geom_smooth(method = "lm", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("Sr/Ca [mol/mol]")+
  ylab("Mg/Ca [mol/mol]")+
  theme_article()+theme(legend.position = 'none')+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)



plt_all_icp <- ggarrange(p_Ca, p_Mg, p_Sr, p_Sr_Ca, p_Mg_Ca, p_Mg_Sr_Ca)


ggsave(filename = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/plots/Ca_sau_gergal.png",
       plot = p_Ca, scale = 1, width = 5, height = 4, units = 'in', dpi = 300)

ggsave(filename = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/plots/icp_sau_gergal.png",
       plot = plt_all_icp, width = 10, height = 8, units = 'in', dpi = 300)



ggplot(icp_avg[icp_avg$sample_type=="water",], aes(date, Na, colour = depth))+
  geom_smooth(method = "loess", aes(fill = depth), alpha = 0.2)+
  geom_point(size=2, alpha = 0.5)+
  xlab("")+
  # ylab("Ca [ppm]")+
  theme_article()+theme(legend.position = c(0.1,0.7))+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8)+facet_wrap(.~site, nrow = 1)




df_YSI_daily <- as.data.frame(timeAverage(df_YSI, "1 day"))
p_cond <- ggplot(df_YSI_daily, aes(date, SpCond))+geom_path()+theme_article()

ggarrange(p_Ca , p_cond, nrow = 2)


ggplot()+
  geom_smooth(data = icp_avg[icp_avg$sample_type=="water",], aes(date, Ca, colour = depth), method = "loess", alpha = 0.2)+
  geom_point(data = icp_avg[icp_avg$sample_type=="water",], aes(date, Ca, colour = depth), size=2, alpha = 0.5)+
  geom_line(data = df_YSI_daily, aes(as.Date(date), SpCond/10))+
  theme_article()+facet_wrap(.~sample_type, nrow = 2, scales = "free_y")















# ------------------- Link with DOC and alk data -------------------------

require(readxl)
alk_doc <- read_xlsx("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/Sed_water_sampling/DOC_Alk.xlsx",
                     sheet = "RESUMEN DATOS")

alk_doc$uniqID <- paste(as.Date(alk_doc$`Date of sampling`),"-",alk_doc$`Depth (m)`,"m",sep = "")
icp_avg$uniqID <- paste(as.Date(alk_doc$`Date of sampling`),"-",alk_doc$`Depth (m)`,"m",sep = "")


plot(alk_doc$`Date of sampling`, alk_doc$`Alcalinity (meq/L)`)






