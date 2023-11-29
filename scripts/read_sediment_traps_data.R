
library(ggplot2)
library(stringi)
library(stringr)
library(egg)
library(lubridate)


SedTraps <- readxl::read_xlsx("C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/Sed_water_sampling/Sed_traps_Water_sampling_sau.xlsx", sheet = "SedTraps")
SedTraps$depth <- as.factor(SedTraps$depth)

SedTraps$middle_date <- SedTraps$deploy_date+SedTraps$n_days/2*3600*24

plt_rate <- ggplot(SedTraps[SedTraps$middle_date<as.Date("2022-11-01"),], aes(middle_date, settling_rate_g_m2_d, colour = depth))+
  geom_smooth(method = "loess", alpha = 0.2, aes(fill = depth))+
  ylab("Sedimentation rate [g/m2/day]")+
  xlab("")+
  geom_point()+theme_article()+
  scale_colour_viridis_d(option = "A", begin = 0.2, end = 0.8)+
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.8, direction = -1)+theme(legend.position = 'none')

plt_rate


ggsave(filename = "C:/Users/Camille Minaudo/OneDrive - Universitat de Barcelona/Documentos/data/sau/plots/settling_rate_sau.png", plot = plt_rate,
       bg = 'white', scale = 1, width = 5, height = 4, units = 'in', dpi = 300)


