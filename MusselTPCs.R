#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(patchwork)

library(ggplot2)

#load TPCs
#https://link.springer.com/article/10.1007/s00442-012-2486-6
#use M. trossulus cleanance rate?

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PLoSextremes/MusselData/")
cr= read.csv("FlyHilbishTPC.csv")

#omit edulis for simplicity
mussel.seas.fig= ggplot(data=cr[cr$species %in% c("M. galloprovincialis", "M. trossulus"),], 
       aes(x=temp_C, y =clearance_ml_min, color=season))+
  facet_wrap(~species)+
  geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  scale_color_viridis_d()+
  theme(legend.position = c(0.8, 0.85))+
  xlab("temperature (°C)")+
  ylab("clearance rate (ml/min)")

#fit TPCs in MusselTPCs_fits.R

#estimate performance
perf= beta_2012(1:40, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])
plot(1:40, perf, type="l")

#--------
#Look for data from M. californicus, if not: trosollus or gallo?

#Use TPC data from Monaco et al. 2016? But only SMR data
#https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1230

#rough Gilman et al. data here: https://doi.org/10.1093/icb/icx031
#assimilation

#--------------------
#load variable TPCs
#https://doi.org/10.1111/1365-2435.13889
#feeding rate, Fig 2d
#shell length growth, Fig 2a
#Mytilus edulis is the dominant species, small fractions of M. galloprovincialis and M. trossulus

#data
#https://doi.pangaea.de/10.1594/PANGAEA.933828
#https://doi.pangaea.de/10.1594/PANGAEA.897938

#load TPC data
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PLoSextremes/MusselData/")

st= read.csv("Shortterm_experiment.csv")
lt= read.csv("Longterm_growth.csv")

#short term
st.fig<- ggplot(data=st, aes(x=Temp_C_x, y =WS_feed_J_per_h_S, color=replicate, lty=trial_name))+
  geom_point()+ theme_classic()

#long term
lt.mean= lt %>%
  group_by(Mean.temperature...C.,Fluctuation.scenario) %>%
  summarise(mean=mean(Length.growth..mm.day.), sd=sd(Length.growth..mm.day.), n=n(), se=sd/sqrt(n) )

lt.fig<- ggplot(data=lt.mean, aes(x=Mean.temperature...C., y =mean, color=Fluctuation.scenario))+
  geom_point()+ geom_line(alpha=0.8, lwd=1) +theme_classic(base_size = 20)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  ylab("shell length growth (mm/day)")+
  xlab("thermal average (°C)")+
  scale_color_viridis_d("fluctuation")+
  theme(legend.position = c(0.3, 0.3))
#change to GAM

#--------
#combine
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PLoSextremes/figures/")
pdf("Fig1_musselTPCs.pdf",height = 6, width = 11)
mussel.seas.fig+ lt.fig + plot_layout(ncol = 2, widths=c(2,1.5))+ 
  plot_annotation(tag_levels = 'A')
dev.off()

#------------------
#Apply TPC to predict performance

#Bernhardt: https://doi.org/10.1098/rspb.2018.1076
#https://github.com/JoeyBernhardt/thermal-variability
#Text S1

#?? Calculate ith derivative of polynomial TPC, Calculate first 10 moments of temperature time series
#just calcualte performance using TPC?



