#load libraries
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)

library(ggplot2)

#load TPCs
#https://link.springer.com/article/10.1007/s00442-012-2486-6


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

#short term

