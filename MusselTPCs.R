#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

library(ggplot2)

#load TPCs
#https://link.springer.com/article/10.1007/s00442-012-2486-6
#use M. trossulus cleanance rate?

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PLoSextremes/MusselData/")
cr= read.csv("FlyHilbishTPC.csv")

ggplot(data=cr, aes(x=temp_C, y =clearance_ml_min, color=season))+
  facet_wrap(~species)+
  geom_point()+ geom_line(alpha=0.8) +theme_classic()

#fit TPCs in MusselTPCs_fits.R

#estimate performance
perf= beta_2012(1:40, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])
plot(1:40, perf, type="l")

#--------
#Look for data from M. californicus, if not: trosollus or gallo?

#Use TPC data from Monaco et al. 2016? But only SMR data
#https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1230

#rough data here: https://doi.org/10.1093/icb/icx031


#--------------------
#load variable TPCs
#https://doi.org/10.1111/1365-2435.13889
#feeding rate, Fig 2d
#shell length growth, Fig 2a

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
  geom_point()+ geom_line(alpha=0.8) +theme_classic()+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  ylab("Shelllength growth (mm/day)")+
  xlab("Thermal average (°C)")
#change to GAM

#------------------
#Apply TPC to predict performance

#Bernhardt: https://doi.org/10.1098/rspb.2018.1076
#https://github.com/JoeyBernhardt/thermal-variability
#Text S1

#?? Calculate ith derivative of polynomial TPC, Calculate first 10 moments of temperature time series
#just calcualte performance using TPC?



